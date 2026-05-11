"""Emit a runnable ``scripts/run_gtdbtk.sh`` for MAGs that failed 16S extraction.

The Phase 1.5 SSU extraction step labels every MAG with a routing decision:

* :data:`src.ssu_extractor.SSU_DECISION_SILVA` - run SILVA SINA on the
  per-MAG 16S FASTA bundle.
* :data:`src.ssu_extractor.SSU_DECISION_GTDBTK` - run GTDB-Tk's genome-wide
  marker classification because no 16S record was recoverable.

This module emits the shell script for the GTDB-Tk path. Mirrors the
contract of :mod:`src.checkm2_shellscript_generator`: we never invoke the
external tool ourselves, only generate a ``set -euo pipefail`` script the
user runs on whatever machine has GTDB-Tk + the GTDB reference data
installed.

The generator is a pure transformation of:

* the Phase 1.5 routing log (``data/ssu/routing_log.tsv``), filtered to
  ``decision == "gtdbtk"``;
* the assembly FASTA paths recorded in that log.

A separate ``gtdbtk_batchfile.tsv`` is also written; ``gtdbtk classify_wf``
is invoked with ``--batchfile`` so adding or removing MAGs is a one-line
edit, not a re-run of this generator.
"""

from __future__ import annotations

import csv
import gzip
import shutil
import stat
from dataclasses import dataclass, field
from pathlib import Path

from .ssu_extractor import SSU_DECISION_GTDBTK


GTDBTK_THREADS_DEFAULT: int = 8
"""Thread count baked into the generated ``gtdbtk classify_wf`` invocation.
Conservative default that fits a typical 16-core HPC node."""

GTDBTK_RESULTS_DIRNAME: str = "gtdbtk_results"
"""Directory under ``data/`` that GTDB-Tk writes into."""

GTDBTK_BATCHFILE_NAME: str = "gtdbtk_batchfile.tsv"
"""Filename for the GTDB-Tk batchfile; lives alongside the routing log."""

GTDBTK_SCRIPT_RELPATH: str = "scripts/run_gtdbtk.sh"
"""Conventional location of the generated script, relative to the repo root."""

GTDBTK_FASTA_EXTENSION: str = "fna"
"""Default ``--extension`` argument for GTDB-Tk. We materialize every input
FASTA to a plain ``.fna`` (un-gzipping ENA ``.fasta.gz`` outputs) so this
extension works for the whole batch."""


@dataclass
class GtdbtkScriptResult:
    """Outcome of one :func:`generate_gtdbtk_script` invocation.

    Attributes:
        script_path: Path the script was written to.
        batchfile_path: Path the GTDB-Tk batchfile was written to.
        n_assemblies: Number of MAGs the script will run GTDB-Tk on.
        skipped_accessions: GTDB-Tk-routed MAGs that could not be included
            (FASTA missing on disk, etc.). Surfaced for the CLI to log.
    """

    script_path: Path
    batchfile_path: Path
    n_assemblies: int
    skipped_accessions: list[str] = field(default_factory=list)


def _read_routing_log(path: Path) -> list[dict[str, str]]:
    """Read the Phase 1.5 routing log into header-keyed dicts."""
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _materialize_plain_fasta(
    src_path: Path,
    accession: str,
    staged_dir: Path,
) -> Path:
    """Return a path to a plain (non-gzipped) FASTA usable by GTDB-Tk.

    GTDB-Tk's ``classify_wf`` reads plain FASTA files; ENA ``.fasta.gz``
    inputs are decompressed once into ``staged_dir``. Plain FASTAs are
    returned unchanged.

    Args:
        src_path: Input assembly FASTA (plain or gzipped).
        accession: MAG accession (used as the staged filename stem).
        staged_dir: Directory holding the staged plain FASTAs.

    Returns:
        Absolute path to the plain FASTA. The staged file is created on
        demand and reused on subsequent runs (idempotent).
    """
    src_path = src_path.resolve()
    if src_path.suffix != ".gz":
        return src_path
    staged_dir.mkdir(parents=True, exist_ok=True)
    staged = (staged_dir / f"{accession}.{GTDBTK_FASTA_EXTENSION}").resolve()
    if staged.exists() and staged.stat().st_size > 0:
        return staged
    with gzip.open(src_path, "rb") as src_fh, staged.open("wb") as dst_fh:
        shutil.copyfileobj(src_fh, dst_fh, length=1024 * 1024)
    return staged


def _write_batchfile(
    rows: list[tuple[Path, str]],
    batchfile_path: Path,
) -> None:
    """Write the GTDB-Tk batchfile (one ``<fasta>\\t<accession>`` per row)."""
    batchfile_path.parent.mkdir(parents=True, exist_ok=True)
    with batchfile_path.open("w", encoding="utf-8") as fh:
        for fasta_path, accession in rows:
            fh.write(f"{fasta_path}\t{accession}\n")


def _make_executable(path: Path) -> None:
    """Set the user/group/other execute bits on a file."""
    current = path.stat().st_mode
    path.chmod(current | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def generate_gtdbtk_script(
    routing_log_path: Path,
    data_root: Path,
    script_path: Path,
    *,
    threads: int = GTDBTK_THREADS_DEFAULT,
    extension: str = GTDBTK_FASTA_EXTENSION,
) -> GtdbtkScriptResult:
    """Render a runnable ``run_gtdbtk.sh`` from the Phase 1.5 routing log.

    Args:
        routing_log_path: Path to ``data/ssu/routing_log.tsv``.
        data_root: Root data directory (used to root the staged FASTA cache,
            the batchfile path, and the GTDB-Tk output directory).
        script_path: Output path for the shell script.
        threads: Threads passed to ``gtdbtk classify_wf --cpus``.
        extension: ``--extension`` argument; defaults to ``fna`` because
            every staged FASTA lands as ``.fna``.

    Returns:
        A populated :class:`GtdbtkScriptResult`.
    """
    rows = _read_routing_log(routing_log_path)
    gtdbtk_rows = [
        r for r in rows
        if (r.get("decision") or "") == SSU_DECISION_GTDBTK
        and (r.get("accession") or "")
    ]

    staged_dir = data_root / "ssu" / "gtdbtk_input"
    batchfile_path = data_root / "ssu" / GTDBTK_BATCHFILE_NAME
    results_dir = data_root / GTDBTK_RESULTS_DIRNAME

    batch_rows: list[tuple[Path, str]] = []
    skipped: list[str] = []
    for row in gtdbtk_rows:
        accession = row["accession"]
        fasta_str = row.get("fasta_path", "")
        if not fasta_str:
            skipped.append(accession)
            continue
        fasta_path = Path(fasta_str)
        if not fasta_path.exists():
            skipped.append(accession)
            continue
        plain_fasta = _materialize_plain_fasta(fasta_path, accession, staged_dir)
        batch_rows.append((plain_fasta, accession))

    _write_batchfile(batch_rows, batchfile_path)

    body_lines = [
        "#!/usr/bin/env bash",
        "# Generated by src/gtdbtk_shellscript_generator.py.",
        "# Re-run `python -m src.pipeline_phase1_5 gen-gtdbtk-script` to refresh.",
        "#",
        "# Genome-wide marker classification for MAGs that failed 16S extraction in",
        "# Phase 1.5 (no barrnap-detected SSU). Requires GTDB-Tk + the GTDB reference",
        "# data (~141 GB at R226). Export GTDBTK_DATA_PATH to that reference data",
        "# directory before running this script.",
        "set -euo pipefail",
        "",
        ': "${GTDBTK_DATA_PATH:?required - export GTDBTK_DATA_PATH pointing to '
        'the GTDB-Tk reference data before running}"',
        "",
        f'mkdir -p "{results_dir}"',
        "",
        f'echo "Running GTDB-Tk classify_wf for {len(batch_rows)} MAG(s) ..."',
        "gtdbtk classify_wf \\",
        f"    --batchfile \"{batchfile_path}\" \\",
        f"    --out_dir \"{results_dir}\" \\",
        f"    --extension {extension} \\",
        f"    --cpus {threads} \\",
        "    --skip_ani_screen",
        "",
        'echo "GTDB-Tk classification complete. Results: "',
        f'echo "  {results_dir}/gtdbtk.bac120.summary.tsv"',
        "",
    ]

    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text("\n".join(body_lines), encoding="utf-8")
    _make_executable(script_path)

    return GtdbtkScriptResult(
        script_path=script_path,
        batchfile_path=batchfile_path,
        n_assemblies=len(batch_rows),
        skipped_accessions=skipped,
    )


__all__ = [
    "GTDBTK_BATCHFILE_NAME",
    "GTDBTK_FASTA_EXTENSION",
    "GTDBTK_RESULTS_DIRNAME",
    "GTDBTK_SCRIPT_RELPATH",
    "GTDBTK_THREADS_DEFAULT",
    "GtdbtkScriptResult",
    "generate_gtdbtk_script",
]
