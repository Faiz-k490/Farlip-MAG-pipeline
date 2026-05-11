"""Emit a runnable ``scripts/run_checkm2.sh`` for opt-in CheckM2 cross-validation.

CheckM2 is the gold-standard tool for MAG completeness / contamination but
ships a ~3 GB DIAMOND reference database and depends on a conda environment we
deliberately do not pin to the dev box (Phase 3 plan, "Hybrid" tool strategy).
This module emits a self-contained shell script the user can run on whatever
machine has CheckM2 installed; results land in ``data/checkm2_results/`` and
are ingested by :mod:`src.checkm2_ingestor`.

Generated script properties:

* ``set -euo pipefail`` so a single failure aborts the run rather than masking
  errors with an exit-zero outer script.
* ``CHECKM2DB`` is sourced from the shell environment via
  ``${CHECKM2DB:?required}``; we never embed an absolute path.
* Each MAG is invoked with ``checkm2 predict --genes`` against the per-MAG
  ``data/orfs/<accession>.faa`` we already produced in Phase 2 - CheckM2 skips
  its own Prodigal gene-prediction step, which is otherwise the slowest part
  of the run on a small machine.
* The output directory is per-MAG (``data/checkm2_results/<accession>/``) so
  re-runs of one MAG don't clobber the rest.

The generator is a pure transformation of ``data/manifest.tsv`` plus the
existing ``data/orfs/<accession>.faa`` files - no execution, no network, no
CheckM2 import. Tests verify the script syntax with ``bash -n`` only.
"""

from __future__ import annotations

import stat
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .manifest import read_manifest


CHECKM2_THREADS_DEFAULT: int = 8
"""Thread count baked into every ``checkm2 predict`` line. Conservative default
that fits a typical 16-core HPC node."""

CHECKM2_RESULTS_DIRNAME: str = "checkm2_results"
"""Subdirectory of ``data/`` that CheckM2 writes into. Mirrored by
:mod:`src.checkm2_ingestor`."""

NON_ASSEMBLY_SOURCES: frozenset[str] = frozenset({"SILVA", "GTDB", "QC_MARKERS"})
"""Manifest sources that are not genome assemblies and must be skipped when
listing MAGs for CheckM2."""

CHECKM2_SCRIPT_RELPATH: str = "scripts/run_checkm2.sh"
"""Conventional location of the generated script, relative to the repo root."""


@dataclass(frozen=True)
class CheckM2ScriptResult:
    """Outcome of one :func:`generate_checkm2_script` invocation.

    Attributes:
        script_path: Path the script was written to.
        n_assemblies: Number of CheckM2 invocations baked into the script.
        skipped_accessions: MAGs that were eligible by manifest source but had
            no per-MAG ``.faa`` on disk yet (so CheckM2 cannot run on them
            until Phase 2 ``orf-predict`` is run).
    """

    script_path: Path
    n_assemblies: int
    skipped_accessions: list[str]


def _eligible_rows(manifest_rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    """Filter Phase 1 manifest rows to genome-assembly entries only.

    Args:
        manifest_rows: Rows as produced by :func:`src.manifest.read_manifest`.

    Returns:
        List of rows whose ``source`` is not in :data:`NON_ASSEMBLY_SOURCES`
        and whose status indicates the file is on disk.
    """
    out: list[dict[str, str]] = []
    for row in manifest_rows:
        if row.get("status") not in {"ok", "skipped"}:
            continue
        if row.get("source", "") in NON_ASSEMBLY_SOURCES:
            continue
        if not row.get("accession"):
            continue
        out.append(row)
    return out


def _format_invocation(
    accession: str,
    faa_path: Path,
    output_dir: Path,
    *,
    threads: int,
) -> str:
    """Render one ``checkm2 predict`` shell command line.

    Args:
        accession: MAG accession (used only for the progress ``echo``).
        faa_path: Per-MAG protein FASTA path.
        output_dir: Per-MAG output directory.
        threads: Thread count.

    Returns:
        A multi-line shell snippet (echo + invocation).
    """
    return (
        f'echo "[{accession}] running checkm2 predict --genes ..."\n'
        f"checkm2 predict \\\n"
        f"    --threads {threads} \\\n"
        f"    --genes \\\n"
        f'    --input "{faa_path}" \\\n'
        f'    --output-directory "{output_dir}" \\\n'
        f'    --database_path "$CHECKM2DB"'
    )


def generate_checkm2_script(
    manifest_path: Path,
    data_root: Path,
    script_path: Path,
    *,
    threads: int = CHECKM2_THREADS_DEFAULT,
) -> CheckM2ScriptResult:
    """Render a runnable ``run_checkm2.sh`` from the Phase 1 manifest.

    Args:
        manifest_path: Path to ``data/manifest.tsv``.
        data_root: Root data directory (used to locate
            ``data/orfs/<acc>.faa`` and to root the per-MAG output dirs).
        script_path: Output path for the shell script.
        threads: Thread count baked into each ``checkm2 predict`` line.

    Returns:
        A :class:`CheckM2ScriptResult` with the script path and counts.
    """
    rows = _eligible_rows(read_manifest(manifest_path))
    orfs_dir = data_root / "orfs"
    results_dir = data_root / CHECKM2_RESULTS_DIRNAME

    invocations: list[str] = []
    skipped: list[str] = []
    for row in rows:
        accession = row["accession"]
        faa_path = orfs_dir / f"{accession}.faa"
        if not faa_path.exists():
            skipped.append(accession)
            continue
        per_mag_out = results_dir / accession
        invocations.append(
            _format_invocation(
                accession,
                faa_path,
                per_mag_out,
                threads=threads,
            )
        )

    body_lines = [
        "#!/usr/bin/env bash",
        "# Generated by src/checkm2_shellscript_generator.py.",
        "# Re-run `python -m src.pipeline_phase3 gen-checkm2-script` to refresh.",
        "#",
        "# This script invokes CheckM2 (Chklovski et al. 2023, Nat Methods 20:1203)",
        "# in --genes mode against the Phase 2 protein FASTAs. CheckM2 must already",
        "# be installed and the DIAMOND reference database extracted; export CHECKM2DB",
        "# pointing to that database before invoking this script.",
        "set -euo pipefail",
        "",
        ': "${CHECKM2DB:?required - export CHECKM2DB pointing to the CheckM2 reference database before running}"',
        "",
        f'mkdir -p "{results_dir}"',
        "",
    ]
    body_lines.append(f'echo "Running CheckM2 for {len(invocations)} assemblies ..."')
    body_lines.append("")
    body_lines.extend(invocations)
    body_lines.append("")
    body_lines.append(
        'echo "All assemblies complete. Run: python -m src.pipeline_phase3 ingest-checkm2"'
    )
    body_lines.append("")

    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text("\n".join(body_lines), encoding="utf-8")
    _make_executable(script_path)

    return CheckM2ScriptResult(
        script_path=script_path,
        n_assemblies=len(invocations),
        skipped_accessions=skipped,
    )


def _make_executable(path: Path) -> None:
    """Set the user/group/other execute bit on a file (idempotent)."""
    current = path.stat().st_mode
    path.chmod(current | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


__all__ = [
    "CHECKM2_RESULTS_DIRNAME",
    "CHECKM2_SCRIPT_RELPATH",
    "CHECKM2_THREADS_DEFAULT",
    "CheckM2ScriptResult",
    "NON_ASSEMBLY_SOURCES",
    "generate_checkm2_script",
]
