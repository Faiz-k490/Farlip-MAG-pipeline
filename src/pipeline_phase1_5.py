"""CLI orchestrator for the Phase 1.5 SSU + taxonomic-routing pipeline.

Reads the Phase 1 manifest (``data/manifest.tsv``) and routes every assembly
row through 16S extraction:

1. **extract-16S** - run :func:`src.ssu_extractor.assess_mag_for_ssu` per
   assembly. Each call writes ``data/ssu/<acc>.barrnap.gff`` plus
   ``data/ssu/<acc>.16S.fasta`` (empty file when no 16S is found) and an
   idempotency stamp; the routing log is upserted in place.
2. **bundle-silva** - concatenate every ``decision == silva`` MAG's 16S
   FASTA into ``data/ssu/silva_sina_input.fasta`` for SILVA SINA submission.
3. **gen-gtdbtk-script** - render ``scripts/run_gtdbtk.sh`` plus
   ``data/ssu/gtdbtk_batchfile.tsv`` for every ``decision == gtdbtk`` MAG.
4. **all** - run 1 -> 2 -> 3 in order.
5. **status** - one-screen summary of routing-log decisions and asset
   presence.

Every subcommand is idempotent. Phase 1.5 progress is recorded in
``data/phase1_5_manifest.tsv`` (one row per ``(accession, stage)``) so a
partial ``all`` run can resume from where it stopped.

The pipeline never invokes SILVA SINA or GTDB-Tk itself - those tools require
large reference data and/or web submission. Phase 1.5 only **prepares** their
inputs.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from .gtdbtk_shellscript_generator import (
    GTDBTK_SCRIPT_RELPATH,
    GTDBTK_THREADS_DEFAULT,
    generate_gtdbtk_script,
)
from .manifest import read_manifest
from .ssu_bundler import (
    ROUTING_LOG_FILENAME,
    SILVA_SINA_BUNDLE_FILENAME,
    bundle_silva_sina,
    read_routing_log,
    routing_decision_from_assessment,
    upsert_routing_log,
)
from .ssu_extractor import (
    BARRNAP_DEFAULT_EVALUE,
    BARRNAP_DEFAULT_REJECT,
    BARRNAP_DEFAULT_THREADS,
    BARRNAP_KINGDOM_BACTERIA,
    BarrnapError,
    SSU_DECISION_GTDBTK,
    SSU_DECISION_SILVA,
    assess_mag_for_ssu,
    is_barrnap_available,
)


REPO_ROOT: Path = Path(__file__).resolve().parent.parent
"""Repository root, derived from this module's location."""

DEFAULT_DATA_ROOT: Path = REPO_ROOT / "data"
"""Root of all pipeline outputs."""

DEFAULT_MANIFEST: Path = DEFAULT_DATA_ROOT / "manifest.tsv"
"""Phase 1 manifest, read-only here."""

SSU_DIRNAME: str = "ssu"
"""Subdirectory under ``data/`` that Phase 1.5 writes into."""

PHASE1_5_MANIFEST_COLUMNS: tuple[str, ...] = (
    "accession",
    "source",
    "stage",
    "status",
    "output_path",
    "n_units",
    "elapsed_s",
    "completed_at",
    "note",
)
"""Column order for ``data/phase1_5_manifest.tsv``."""

NON_ASSEMBLY_SOURCES: frozenset[str] = frozenset({"SILVA", "GTDB", "QC_MARKERS"})
"""Manifest sources whose rows are skipped during 16S extraction."""


@dataclass
class _PhaseRow:
    """One row written to ``data/phase1_5_manifest.tsv``."""

    accession: str
    source: str
    stage: str
    status: str
    output_path: str
    n_units: int
    elapsed_s: float
    completed_at: str
    note: str = ""


def _utc_now() -> str:
    """Return a UTC ISO-8601 timestamp matching the Phase 1/2/3 manifests."""
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _read_phase_manifest(path: Path) -> list[dict[str, str]]:
    """Load ``data/phase1_5_manifest.tsv`` if it exists."""
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_phase_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    """Write the Phase 1.5 manifest TSV using the canonical column order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=PHASE1_5_MANIFEST_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in PHASE1_5_MANIFEST_COLUMNS})


def _upsert_phase_rows(path: Path, new_rows: list[_PhaseRow]) -> None:
    """Merge new rows into the Phase 1.5 manifest, replacing on (accession, stage)."""
    existing = _read_phase_manifest(path)
    keyed: dict[tuple[str, str], dict[str, str]] = {
        (r.get("accession", ""), r.get("stage", "")): r for r in existing
    }
    for r in new_rows:
        keyed[(r.accession, r.stage)] = {
            "accession": r.accession,
            "source": r.source,
            "stage": r.stage,
            "status": r.status,
            "output_path": r.output_path,
            "n_units": str(r.n_units),
            "elapsed_s": f"{r.elapsed_s:.2f}",
            "completed_at": r.completed_at,
            "note": r.note,
        }
    _write_phase_manifest(path, list(keyed.values()))


def _iter_assembly_rows(
    manifest_path: Path,
    accession_filter: list[str] | None = None,
) -> list[dict[str, str]]:
    """Yield every Phase 1 manifest row eligible for 16S extraction."""
    rows = read_manifest(manifest_path)
    out: list[dict[str, str]] = []
    filter_set = set(accession_filter or [])
    seen: set[str] = set()
    for row in rows:
        accession = row.get("accession", "")
        source = row.get("source", "")
        status = row.get("status", "")
        if not accession or accession in seen:
            continue
        if source in NON_ASSEMBLY_SOURCES:
            continue
        if status not in {"ok", "skipped"}:
            continue
        if filter_set and accession not in filter_set:
            continue
        seen.add(accession)
        out.append(row)
    return out


def _resolve_assembly_fasta(row: dict[str, str], data_root: Path) -> Path | None:
    """Resolve the FASTA path for a Phase 1 manifest row.

    NCBI assemblies arrive as ZIP bundles; the unpacked genomic FASTA lives
    under ``data/.../<accession>/ncbi_dataset/data/<accession>/*_genomic.fna``
    or similar. ENA assemblies arrive as ``<accession>.fasta.gz`` written
    directly to ``local_path``.

    Args:
        row: Manifest row.
        data_root: Root data directory (used to resolve relative
            ``extracted_dir`` columns).

    Returns:
        Absolute path to a FASTA file, or ``None`` if nothing usable was
        found on disk.
    """
    source = (row.get("source") or "").upper()
    local_path = Path(row.get("local_path") or "")
    if source == "ENA":
        if local_path.exists():
            return local_path.resolve()
        return None
    if source == "NCBI":
        extracted_dir_str = row.get("extracted_dir") or ""
        if extracted_dir_str:
            extracted_dir = Path(extracted_dir_str)
            if not extracted_dir.is_absolute():
                extracted_dir = (data_root / extracted_dir).resolve()
            candidates = list(extracted_dir.rglob("*_genomic.fna")) + list(
                extracted_dir.rglob("*_genomic.fna.gz")
            )
            if candidates:
                # Prefer the largest file (the assembly FASTA is the only
                # large *_genomic.fna under the bundle).
                candidates.sort(key=lambda p: p.stat().st_size, reverse=True)
                return candidates[0].resolve()
        if local_path.exists() and local_path.suffix in {".fna", ".fasta", ".gz"}:
            return local_path.resolve()
        return None
    # Unknown source - fall back to whatever local_path was recorded.
    if local_path.exists():
        return local_path.resolve()
    return None


def cmd_extract_16S(
    manifest_path: Path,
    data_root: Path,
    *,
    overwrite: bool,
    accession_filter: list[str] | None,
    threads: int,
    kingdom: str,
    evalue: float,
    reject: float,
) -> int:
    """Run barrnap-based 16S extraction for every eligible MAG.

    Returns:
        Exit code: 0 on success, 1 if any MAG failed, 2 if barrnap is
        missing from PATH (the user is shown a clear error).
    """
    if not is_barrnap_available():
        sys.stderr.write(
            "[extract-16S] barrnap not found on PATH. "
            "Install via conda-forge / bioconda before running Phase 1.5.\n"
        )
        return 2

    rows = _iter_assembly_rows(manifest_path, accession_filter)
    ssu_dir = data_root / SSU_DIRNAME
    routing_log_path = ssu_dir / ROUTING_LOG_FILENAME
    phase_manifest = data_root / "phase1_5_manifest.tsv"

    sys.stderr.write(
        f"[extract-16S] {len(rows)} assembly row(s) eligible "
        f"(barrnap kingdom={kingdom}, threads={threads})\n"
    )

    new_phase_rows: list[_PhaseRow] = []
    routing_decisions = []
    failures = 0
    for row in rows:
        accession = row["accession"]
        source = row.get("source", "")
        fasta_path = _resolve_assembly_fasta(row, data_root)
        if fasta_path is None:
            sys.stderr.write(
                f"  ! {accession}: could not resolve a FASTA on disk; skipping\n"
            )
            failures += 1
            new_phase_rows.append(
                _PhaseRow(
                    accession=accession,
                    source=source,
                    stage="extract-16S",
                    status="failed",
                    output_path="",
                    n_units=0,
                    elapsed_s=0.0,
                    completed_at=_utc_now(),
                    note="FASTA not resolved on disk",
                )
            )
            continue
        t0 = time.perf_counter()
        try:
            assessment = assess_mag_for_ssu(
                accession=accession,
                fasta_path=fasta_path,
                output_dir=ssu_dir,
                kingdom=kingdom,
                threads=threads,
                evalue=evalue,
                reject=reject,
                overwrite=overwrite,
            )
        except (BarrnapError, OSError) as exc:
            failures += 1
            elapsed = time.perf_counter() - t0
            sys.stderr.write(f"  ! {accession}: {exc}\n")
            new_phase_rows.append(
                _PhaseRow(
                    accession=accession,
                    source=source,
                    stage="extract-16S",
                    status="failed",
                    output_path="",
                    n_units=0,
                    elapsed_s=elapsed,
                    completed_at=_utc_now(),
                    note=str(exc),
                )
            )
            continue
        elapsed = time.perf_counter() - t0
        marker = "=" if assessment.skipped else "+"
        sys.stderr.write(
            f"  {marker} {accession}: {assessment.n_16S_predicted} 16S record(s); "
            f"decision={assessment.decision} ({elapsed:.1f}s)\n"
        )
        routing_decisions.append(
            routing_decision_from_assessment(
                assessment,
                source=source,
                notes="skipped (stamp matched)" if assessment.skipped else "",
            )
        )
        new_phase_rows.append(
            _PhaseRow(
                accession=accession,
                source=source,
                stage="extract-16S",
                status="skipped" if assessment.skipped else "ok",
                output_path=str(assessment.ssu_fasta_path),
                n_units=assessment.n_16S_predicted,
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note=f"decision={assessment.decision}",
            )
        )

    if routing_decisions:
        upsert_routing_log(routing_log_path, routing_decisions)
    _upsert_phase_rows(phase_manifest, new_phase_rows)
    return 1 if failures else 0


def cmd_bundle_silva(data_root: Path) -> int:
    """Concatenate per-MAG 16S FASTAs into the SILVA SINA input bundle."""
    ssu_dir = data_root / SSU_DIRNAME
    routing_log_path = ssu_dir / ROUTING_LOG_FILENAME
    bundle_path = ssu_dir / SILVA_SINA_BUNDLE_FILENAME
    phase_manifest = data_root / "phase1_5_manifest.tsv"

    if not routing_log_path.exists():
        sys.stderr.write(
            f"[bundle-silva] routing log missing at {routing_log_path}; "
            "run extract-16S first.\n"
        )
        return 2

    t0 = time.perf_counter()
    result = bundle_silva_sina(routing_log_path, ssu_dir, bundle_path)
    elapsed = time.perf_counter() - t0
    sys.stderr.write(
        f"[bundle-silva] wrote {bundle_path} "
        f"(MAGs included={result.n_mags_included}, "
        f"records={result.n_records_total}, "
        f"skipped={len(result.skipped_accessions)}) "
        f"({elapsed:.1f}s)\n"
    )
    _upsert_phase_rows(
        phase_manifest,
        [
            _PhaseRow(
                accession="*",
                source="SILVA",
                stage="bundle-silva",
                status="ok",
                output_path=str(bundle_path),
                n_units=result.n_records_total,
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note=f"mags={result.n_mags_included}; "
                f"skipped={len(result.skipped_accessions)}",
            )
        ],
    )
    return 0


def cmd_gen_gtdbtk_script(data_root: Path, *, threads: int) -> int:
    """Render ``scripts/run_gtdbtk.sh`` for every GTDB-Tk-routed MAG."""
    ssu_dir = data_root / SSU_DIRNAME
    routing_log_path = ssu_dir / ROUTING_LOG_FILENAME
    script_path = REPO_ROOT / GTDBTK_SCRIPT_RELPATH
    phase_manifest = data_root / "phase1_5_manifest.tsv"

    if not routing_log_path.exists():
        sys.stderr.write(
            f"[gen-gtdbtk-script] routing log missing at {routing_log_path}; "
            "run extract-16S first.\n"
        )
        return 2

    t0 = time.perf_counter()
    result = generate_gtdbtk_script(
        routing_log_path=routing_log_path,
        data_root=data_root,
        script_path=script_path,
        threads=threads,
    )
    elapsed = time.perf_counter() - t0
    sys.stderr.write(
        f"[gen-gtdbtk-script] wrote {result.script_path} "
        f"(MAGs={result.n_assemblies}, skipped={len(result.skipped_accessions)}) "
        f"({elapsed:.1f}s)\n"
    )
    if result.skipped_accessions:
        sys.stderr.write(
            "  ! skipped (FASTA missing on disk): "
            + ", ".join(result.skipped_accessions[:5])
            + ("..." if len(result.skipped_accessions) > 5 else "")
            + "\n"
        )
    _upsert_phase_rows(
        phase_manifest,
        [
            _PhaseRow(
                accession="*",
                source="GTDB",
                stage="gen-gtdbtk-script",
                status="ok",
                output_path=str(result.script_path),
                n_units=result.n_assemblies,
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note=f"batchfile={result.batchfile_path}; "
                f"skipped={len(result.skipped_accessions)}",
            )
        ],
    )
    return 0


def cmd_status(data_root: Path) -> int:
    """Print a one-screen summary of Phase 1.5 state."""
    ssu_dir = data_root / SSU_DIRNAME
    routing_log_path = ssu_dir / ROUTING_LOG_FILENAME
    bundle_path = ssu_dir / SILVA_SINA_BUNDLE_FILENAME
    script_path = REPO_ROOT / GTDBTK_SCRIPT_RELPATH

    rows = read_routing_log(routing_log_path)
    silva_count = sum(1 for r in rows if r.get("decision") == SSU_DECISION_SILVA)
    gtdbtk_count = sum(1 for r in rows if r.get("decision") == SSU_DECISION_GTDBTK)

    sys.stdout.write(
        "FaRLiP Phase 1.5 status\n"
        f"  routing log               : "
        f"{'present' if routing_log_path.exists() else 'absent'} "
        f"({routing_log_path})\n"
        f"  routed -> SILVA           : {silva_count}\n"
        f"  routed -> GTDB-Tk         : {gtdbtk_count}\n"
        f"  silva_sina_input.fasta    : "
        f"{'present' if bundle_path.exists() else 'absent'}\n"
        f"  scripts/run_gtdbtk.sh     : "
        f"{'present' if script_path.exists() else 'absent'}\n"
        f"  barrnap on PATH           : "
        f"{'yes' if is_barrnap_available() else 'no'}\n"
    )
    return 0


def _build_parser() -> argparse.ArgumentParser:
    """Construct the argparse tree for ``python -m src.pipeline_phase1_5``."""
    p = argparse.ArgumentParser(
        prog="src.pipeline_phase1_5",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help=f"Phase 1 manifest path (default: {DEFAULT_MANIFEST})",
    )
    p.add_argument(
        "--data-root",
        type=Path,
        default=DEFAULT_DATA_ROOT,
        help=f"Data root directory (default: {DEFAULT_DATA_ROOT})",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    common_filter = argparse.ArgumentParser(add_help=False)
    common_filter.add_argument(
        "--accession",
        action="append",
        default=None,
        help="Restrict to this accession (repeatable).",
    )
    common_overwrite = argparse.ArgumentParser(add_help=False)
    common_overwrite.add_argument(
        "--overwrite",
        action="store_true",
        help="Force re-run even when the stamp matches.",
    )

    barrnap_args = argparse.ArgumentParser(add_help=False)
    barrnap_args.add_argument(
        "--threads",
        type=int,
        default=BARRNAP_DEFAULT_THREADS,
        help=f"Threads for barrnap (default: {BARRNAP_DEFAULT_THREADS}).",
    )
    barrnap_args.add_argument(
        "--kingdom",
        choices=("bac", "arc", "euk", "mito"),
        default=BARRNAP_KINGDOM_BACTERIA,
        help="Kingdom passed to barrnap (default: bac).",
    )
    barrnap_args.add_argument(
        "--evalue",
        type=float,
        default=BARRNAP_DEFAULT_EVALUE,
        help=f"barrnap --evalue (default: {BARRNAP_DEFAULT_EVALUE:g}).",
    )
    barrnap_args.add_argument(
        "--reject",
        type=float,
        default=BARRNAP_DEFAULT_REJECT,
        help=f"barrnap --reject (default: {BARRNAP_DEFAULT_REJECT:g}).",
    )

    sub.add_parser(
        "extract-16S",
        parents=[common_overwrite, common_filter, barrnap_args],
        help="Run barrnap on every eligible assembly + write the routing log.",
    )
    sub.add_parser(
        "bundle-silva",
        help="Concatenate per-MAG 16S FASTAs into silva_sina_input.fasta.",
    )
    gtdb_p = sub.add_parser(
        "gen-gtdbtk-script",
        help="Render scripts/run_gtdbtk.sh for every GTDB-Tk-routed MAG.",
    )
    gtdb_p.add_argument(
        "--threads",
        type=int,
        default=GTDBTK_THREADS_DEFAULT,
        help=f"Threads for `gtdbtk classify_wf` (default: {GTDBTK_THREADS_DEFAULT}).",
    )
    sub.add_parser(
        "all",
        parents=[common_overwrite, barrnap_args],
        help="Run extract-16S -> bundle-silva -> gen-gtdbtk-script.",
    )
    sub.add_parser(
        "status",
        help="One-screen summary of Phase 1.5 state.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    """Entry point for ``python -m src.pipeline_phase1_5``."""
    args = _build_parser().parse_args(argv)
    manifest_path: Path = args.manifest
    data_root: Path = args.data_root

    cmd = args.cmd
    accession_filter = getattr(args, "accession", None) or None
    overwrite = bool(getattr(args, "overwrite", False))
    threads = int(getattr(args, "threads", BARRNAP_DEFAULT_THREADS))
    kingdom = getattr(args, "kingdom", BARRNAP_KINGDOM_BACTERIA)
    evalue = float(getattr(args, "evalue", BARRNAP_DEFAULT_EVALUE))
    reject = float(getattr(args, "reject", BARRNAP_DEFAULT_REJECT))

    if cmd == "extract-16S":
        return cmd_extract_16S(
            manifest_path,
            data_root,
            overwrite=overwrite,
            accession_filter=accession_filter,
            threads=threads,
            kingdom=kingdom,
            evalue=evalue,
            reject=reject,
        )
    if cmd == "bundle-silva":
        return cmd_bundle_silva(data_root)
    if cmd == "gen-gtdbtk-script":
        return cmd_gen_gtdbtk_script(data_root, threads=threads)
    if cmd == "all":
        rc = cmd_extract_16S(
            manifest_path,
            data_root,
            overwrite=overwrite,
            accession_filter=None,
            threads=threads,
            kingdom=kingdom,
            evalue=evalue,
            reject=reject,
        )
        if rc != 0:
            return rc
        rc = cmd_bundle_silva(data_root)
        if rc != 0:
            return rc
        return cmd_gen_gtdbtk_script(data_root, threads=GTDBTK_THREADS_DEFAULT)
    if cmd == "status":
        return cmd_status(data_root)
    sys.stderr.write(f"unknown command: {cmd}\n")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
