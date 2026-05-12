"""CLI orchestrator for the FaRLiP detection pipeline (Phase 3).

Reads the Phase 1 manifest (``data/manifest.tsv``) and runs - in order, or as
individually-invoked subcommands - the four-stage detection pipeline:

1. **orf-predict** - :func:`src.orf_predictor.predict_orfs` for every assembly row
   in the manifest (NCBI + ENA sources; SILVA and GTDB rows are skipped).
2. **build-profiles** - :func:`src.profile_builder.build_profiles` against the JSC-1
   reference assembly. Produces ``data/hmm_profiles/farlip_21.hmm``.
3. **search** - :func:`src.hmm_searcher.run_search` for every per-MAG ``.faa``.
4. **parse-synteny** - :func:`src.synteny_parser.classify_mag` and reports.
5. **export-vault** - :func:`src.obsidian_exporter.export_vault`.

Every subcommand is independently re-runnable and idempotent. Phase 3 progress is
tracked in ``data/phase3_manifest.tsv`` (one row per accession per stage), so a
partial ``all`` run can resume from where it stopped.

Subcommands (see ``python -m src.pipeline_phase3 --help``):

* ``orf-predict   [--overwrite] [--accession ...]``
* ``build-profiles [--overwrite]``
* ``search         [--overwrite] [--accession ...]``
* ``parse-synteny  [--overwrite]``
* ``export-vault   [--overwrite]``
* ``all            [--overwrite]``
* ``status``
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from .manifest import read_manifest
from .orf_predictor import predict_orfs
from .profile_builder import build_profiles
from .hmm_searcher import run_search
from .synteny_parser import Classification, classify_mag, write_reports


_CLASSIFICATION_DISPLAY_PRIORITY: dict[str, int] = {
    Classification.FULL_CLUSTER.value: 3,
    Classification.EDGE_TRUNCATED.value: 2,
    Classification.PARTIAL_CLUSTER.value: 1,
    Classification.NONE.value: 0,
}
"""Priority used to pick a per-MAG headline classification for stderr logging."""
from .obsidian_exporter import export_vault


REPO_ROOT: Path = Path(__file__).resolve().parent.parent
"""Repository root, derived from this module's location."""

DEFAULT_DATA_ROOT: Path = REPO_ROOT / "data"
"""Root of all pipeline outputs. Phase 1 already wrote into this directory."""

DEFAULT_MANIFEST: Path = DEFAULT_DATA_ROOT / "manifest.tsv"
"""Phase 1 manifest, consumed read-only by Phase 3."""

DEFAULT_CONFIG_DIR: Path = REPO_ROOT / "config"
"""Directory containing ``farlip_genes.yaml``."""

REFERENCE_ACCESSION: str = "GCF_000733415.1"
"""Primary FaRLiP reference (JSC-1). Profiles are built from this assembly."""

PHASE3_MANIFEST_COLUMNS: tuple[str, ...] = (
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
"""Column order for ``data/phase3_manifest.tsv``. Stable; extended only by appending."""

NON_ASSEMBLY_SOURCES: frozenset[str] = frozenset({"SILVA", "GTDB"})
"""Manifest sources whose rows are skipped during ORF prediction (these are reference
databases, not genome assemblies)."""


@dataclass
class _Phase3Row:
    """One row written to ``data/phase3_manifest.tsv``."""

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
    """Return a UTC ISO-8601 timestamp matching the Phase 1 manifest convention."""
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _read_phase3_manifest(path: Path) -> list[dict[str, str]]:
    """Load ``data/phase3_manifest.tsv`` if it exists.

    Args:
        path: Manifest path.

    Returns:
        Rows as dicts (empty list if file absent).
    """
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_phase3_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    """Write a Phase 3 manifest TSV with the canonical column order.

    Args:
        path: Output path. Parent created if missing.
        rows: Iterable of row dicts. Missing columns are filled with empty strings.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=PHASE3_MANIFEST_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in PHASE3_MANIFEST_COLUMNS})


def _upsert_phase3_rows(path: Path, new_rows: list[_Phase3Row]) -> None:
    """Merge new rows into the Phase 3 manifest, replacing on ``(accession, stage)``.

    Args:
        path: Manifest path.
        new_rows: Rows to merge in.
    """
    existing = _read_phase3_manifest(path)
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
    _write_phase3_manifest(path, list(keyed.values()))


def _resolve_genomic_fasta(row: dict[str, str]) -> Path | None:
    """Locate the genomic FASTA for one Phase 1 manifest row.

    NCBI rows use the unpacked ``ncbi_dataset/data/<accession>/*_genomic.fna``
    inside ``extracted_dir``; ENA rows use ``local_path`` directly.

    Args:
        row: One Phase 1 manifest row.

    Returns:
        Absolute path to the genomic FASTA, or ``None`` if the file cannot be
        located. Returning ``None`` lets the caller skip cleanly rather than
        crashing the whole sweep on one bad row.
    """
    source = row.get("source", "")
    accession = row.get("accession", "")
    if source == "NCBI":
        extracted_dir = row.get("extracted_dir", "")
        if not extracted_dir:
            return None
        base = Path(extracted_dir) / "ncbi_dataset" / "data" / accession
        if not base.exists():
            return None
        candidates = sorted(base.glob("*_genomic.fna"))
        if candidates:
            return candidates[0]
        candidates = sorted(base.glob("*.fna"))
        return candidates[0] if candidates else None
    if source == "ENA":
        local_path = row.get("local_path", "")
        if not local_path:
            return None
        p = Path(local_path)
        return p if p.exists() else None
    return None


def _iter_assembly_rows(
    manifest_path: Path,
    accession_filter: list[str] | None = None,
) -> list[dict[str, str]]:
    """Yield every manifest row eligible for ORF prediction / search.

    Args:
        manifest_path: Phase 1 manifest path.
        accession_filter: If non-empty, restrict to these accessions.

    Returns:
        List of manifest rows in file order.
    """
    rows = read_manifest(manifest_path)
    out: list[dict[str, str]] = []
    filter_set = set(accession_filter or [])
    for row in rows:
        # ``ok`` = freshly downloaded this run; ``skipped`` = already on disk and
        # left in place by an idempotent re-fetch. Both states mean the artifact
        # exists locally and is eligible for Phase 3 processing.
        if row.get("status") not in {"ok", "skipped"}:
            continue
        if row.get("source", "") in NON_ASSEMBLY_SOURCES:
            continue
        if filter_set and row.get("accession", "") not in filter_set:
            continue
        out.append(row)
    return out


def cmd_orf_predict(
    manifest_path: Path,
    data_root: Path,
    *,
    overwrite: bool,
    accession_filter: list[str] | None,
) -> int:
    """Run ORF prediction for every eligible manifest row.

    Args:
        manifest_path: Phase 1 manifest path.
        data_root: Root data directory (outputs land under ``data/orfs``).
        overwrite: If True, re-run even when stamp matches.
        accession_filter: If non-empty, restrict to these accessions.

    Returns:
        Process exit code (0 on success, non-zero if any row errored fatally).
    """
    rows = _iter_assembly_rows(manifest_path, accession_filter)
    output_dir = data_root / "orfs"
    phase3_manifest = data_root / "phase3_manifest.tsv"
    sys.stderr.write(f"[orf-predict] {len(rows)} assembly row(s) eligible\n")
    new_rows: list[_Phase3Row] = []
    failures = 0
    for row in rows:
        accession = row["accession"]
        fasta_path = _resolve_genomic_fasta(row)
        if fasta_path is None:
            sys.stderr.write(f"  ! {accession}: genomic FASTA not found, skipping\n")
            new_rows.append(
                _Phase3Row(
                    accession=accession,
                    source=row.get("source", ""),
                    stage="orf-predict",
                    status="failed",
                    output_path="",
                    n_units=0,
                    elapsed_s=0.0,
                    completed_at=_utc_now(),
                    note="genomic FASTA not located",
                )
            )
            failures += 1
            continue
        t0 = time.perf_counter()
        try:
            result = predict_orfs(
                accession=accession,
                fasta_path=fasta_path,
                output_dir=output_dir,
                overwrite=overwrite,
            )
        except Exception as exc:
            failures += 1
            elapsed = time.perf_counter() - t0
            sys.stderr.write(f"  ! {accession}: {exc}\n")
            new_rows.append(
                _Phase3Row(
                    accession=accession,
                    source=row.get("source", ""),
                    stage="orf-predict",
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
        marker = "=" if result.skipped else "+"
        sys.stderr.write(
            f"  {marker} {accession}: n_contigs={result.n_contigs} n_orfs={result.n_orfs} "
            f"({elapsed:.1f}s)\n"
        )
        new_rows.append(
            _Phase3Row(
                accession=accession,
                source=row.get("source", ""),
                stage="orf-predict",
                status="skipped" if result.skipped else "ok",
                output_path=str(result.faa_path),
                n_units=result.n_orfs,
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note=f"n_contigs={result.n_contigs}",
            )
        )
    _upsert_phase3_rows(phase3_manifest, new_rows)
    return 1 if failures else 0


def cmd_build_profiles(
    data_root: Path,
    config_dir: Path,
    *,
    overwrite: bool,
) -> int:
    """Build the 21 FaRLiP HMM profiles from JSC-1.

    Args:
        data_root: Root data directory.
        config_dir: Directory containing ``farlip_genes.yaml``.
        overwrite: If True, rebuild even when the HMM file exists.

    Returns:
        Process exit code (0 on success).
    """
    output_dir = data_root / "hmm_profiles"
    phase3_manifest = data_root / "phase3_manifest.tsv"
    jsc1_dir = (
        data_root
        / "references"
        / REFERENCE_ACCESSION
        / "ncbi_dataset"
        / "data"
        / REFERENCE_ACCESSION
    )
    gff_path = jsc1_dir / "genomic.gff"
    faa_path = jsc1_dir / "protein.faa"
    catalog_path = config_dir / "farlip_genes.yaml"

    if not gff_path.exists() or not faa_path.exists():
        sys.stderr.write(
            f"[build-profiles] missing JSC-1 GFF or protein FASTA under {jsc1_dir}\n"
        )
        return 2

    t0 = time.perf_counter()
    resolved = build_profiles(
        catalog_path=catalog_path,
        gff_path=gff_path,
        protein_faa_path=faa_path,
        output_dir=output_dir,
        overwrite=overwrite,
    )
    elapsed = time.perf_counter() - t0
    sys.stderr.write(
        f"[build-profiles] resolved {len(resolved)} catalog entries -> "
        f"{output_dir / 'farlip_21.hmm'} ({elapsed:.1f}s)\n"
    )

    _upsert_phase3_rows(
        phase3_manifest,
        [
            _Phase3Row(
                accession=REFERENCE_ACCESSION,
                source="JSC-1",
                stage="build-profiles",
                status="ok",
                output_path=str(output_dir / "farlip_21.hmm"),
                n_units=len(resolved),
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note=f"{len([g for g in resolved if g.in_cluster])} of {len(resolved)} in cluster",
            )
        ],
    )
    return 0


def cmd_search(
    manifest_path: Path,
    data_root: Path,
    *,
    overwrite: bool,
    accession_filter: list[str] | None,
) -> int:
    """Run hmmsearch for every accession with an existing ``.faa``.

    Args:
        manifest_path: Phase 1 manifest path.
        data_root: Root data directory.
        overwrite: If True, re-run even when ``.domtblout`` exists.
        accession_filter: If non-empty, restrict to these accessions.

    Returns:
        Process exit code (0 on success).
    """
    rows = _iter_assembly_rows(manifest_path, accession_filter)
    orfs_dir = data_root / "orfs"
    output_dir = data_root / "hmm_results"
    hmm_path = data_root / "hmm_profiles" / "farlip_21.hmm"
    phase3_manifest = data_root / "phase3_manifest.tsv"

    if not hmm_path.exists():
        sys.stderr.write(
            f"[search] HMM file missing - run build-profiles first ({hmm_path})\n"
        )
        return 2

    sys.stderr.write(f"[search] {len(rows)} assembly row(s) eligible\n")
    new_rows: list[_Phase3Row] = []
    failures = 0
    for row in rows:
        accession = row["accession"]
        faa_path = orfs_dir / f"{accession}.faa"
        if not faa_path.exists():
            sys.stderr.write(
                f"  ! {accession}: missing {faa_path.name} - run orf-predict first\n"
            )
            failures += 1
            continue
        t0 = time.perf_counter()
        try:
            result = run_search(
                accession=accession,
                faa_path=faa_path,
                hmm_path=hmm_path,
                output_dir=output_dir,
                overwrite=overwrite,
            )
        except Exception as exc:
            failures += 1
            elapsed = time.perf_counter() - t0
            sys.stderr.write(f"  ! {accession}: {exc}\n")
            new_rows.append(
                _Phase3Row(
                    accession=accession,
                    source=row.get("source", ""),
                    stage="search",
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
        marker = "=" if result.skipped else "+"
        sys.stderr.write(
            f"  {marker} {accession}: hits={result.n_hits_reported} "
            f"hmms_with_hits={result.n_hmms_with_hits} ({elapsed:.1f}s)\n"
        )
        new_rows.append(
            _Phase3Row(
                accession=accession,
                source=row.get("source", ""),
                stage="search",
                status="skipped" if result.skipped else "ok",
                output_path=str(result.domtblout_path),
                n_units=result.n_hits_reported,
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note=f"hmms_with_hits={result.n_hmms_with_hits}",
            )
        )
    _upsert_phase3_rows(phase3_manifest, new_rows)
    return 1 if failures else 0


def cmd_parse_synteny(
    manifest_path: Path,
    data_root: Path,
    *,
    overwrite: bool,
) -> int:
    """Classify every MAG's HMM hits and write the two synteny reports.

    Args:
        manifest_path: Phase 1 manifest path.
        data_root: Root data directory.
        overwrite: Currently unused (the reports are always rewritten end-to-end -
            the flag is accepted for CLI symmetry with the other subcommands).

    Returns:
        Process exit code (0 on success).
    """
    del overwrite
    rows = _iter_assembly_rows(manifest_path, None)
    orfs_dir = data_root / "orfs"
    hmm_results_dir = data_root / "hmm_results"
    report_path = data_root / "farlip_hits_report.tsv"
    per_mag_path = data_root / "farlip_hits_per_mag.tsv"
    phase3_manifest = data_root / "phase3_manifest.tsv"

    sys.stderr.write(f"[parse-synteny] {len(rows)} assembly row(s) eligible\n")
    all_calls = {}
    new_rows: list[_Phase3Row] = []
    for row in rows:
        accession = row["accession"]
        coords_path = orfs_dir / f"{accession}.coords.tsv"
        domtblout_path = hmm_results_dir / f"{accession}.domtblout"
        if not coords_path.exists() or not domtblout_path.exists():
            sys.stderr.write(
                f"  ! {accession}: missing coords or domtblout - skipping\n"
            )
            continue
        t0 = time.perf_counter()
        calls = classify_mag(accession, coords_path, domtblout_path)
        elapsed = time.perf_counter() - t0
        all_calls[accession] = calls
        top = (
            max(
                (c.classification.value for c in calls),
                key=lambda v: _CLASSIFICATION_DISPLAY_PRIORITY.get(v, -1),
            )
            if calls
            else Classification.NONE.value
        )
        sys.stderr.write(
            f"  + {accession}: candidate_contigs={len(calls)} top={top} ({elapsed:.2f}s)\n"
        )
        new_rows.append(
            _Phase3Row(
                accession=accession,
                source=row.get("source", ""),
                stage="parse-synteny",
                status="ok",
                output_path=str(report_path),
                n_units=len(calls),
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note=f"top_classification={top}",
            )
        )

    write_reports(all_calls, report_path, per_mag_path)
    sys.stderr.write(
        f"[parse-synteny] wrote {report_path} and {per_mag_path}\n"
    )
    _upsert_phase3_rows(phase3_manifest, new_rows)
    return 0


def cmd_export_vault(
    manifest_path: Path,
    data_root: Path,
    *,
    overwrite: bool,
) -> int:
    """Render the Obsidian vault from the synteny reports.

    Args:
        manifest_path: Phase 1 manifest (read for scientific_name lookup).
        data_root: Root data directory.
        overwrite: Currently unused (the exporter always overwrites; flag kept for
            CLI symmetry).

    Returns:
        Process exit code (0 on success).
    """
    del overwrite
    report_path = data_root / "farlip_hits_report.tsv"
    per_mag_path = data_root / "farlip_hits_per_mag.tsv"
    vault_dir = data_root / "obsidian_export"
    phase3_manifest = data_root / "phase3_manifest.tsv"

    if not report_path.exists() or not per_mag_path.exists():
        sys.stderr.write(
            "[export-vault] synteny reports missing - run parse-synteny first\n"
        )
        return 2

    t0 = time.perf_counter()
    n_written = export_vault(report_path, per_mag_path, manifest_path, vault_dir)
    elapsed = time.perf_counter() - t0
    sys.stderr.write(
        f"[export-vault] wrote {n_written} per-MAG markdown files + master index "
        f"to {vault_dir} ({elapsed:.1f}s)\n"
    )
    _upsert_phase3_rows(
        phase3_manifest,
        [
            _Phase3Row(
                accession="*",
                source="export",
                stage="export-vault",
                status="ok",
                output_path=str(vault_dir),
                n_units=n_written,
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note="",
            )
        ],
    )
    return 0


def cmd_status(manifest_path: Path, data_root: Path) -> int:
    """Print a one-screen summary of pipeline state.

    Args:
        manifest_path: Phase 1 manifest.
        data_root: Root data directory.

    Returns:
        Always 0.
    """
    rows = _iter_assembly_rows(manifest_path, None)
    orfs_dir = data_root / "orfs"
    hmm_results_dir = data_root / "hmm_results"
    hmm_path = data_root / "hmm_profiles" / "farlip_21.hmm"
    report_path = data_root / "farlip_hits_report.tsv"
    vault_dir = data_root / "obsidian_export"

    n_with_orfs = sum(
        1 for r in rows if (orfs_dir / f"{r['accession']}.faa").exists()
    )
    n_with_hits = sum(
        1
        for r in rows
        if (hmm_results_dir / f"{r['accession']}.domtblout").exists()
    )

    sys.stdout.write(
        "FaRLiP Phase 3 status\n"
        f"  assembly rows in manifest : {len(rows)}\n"
        f"  ORFs predicted (.faa)      : {n_with_orfs}\n"
        f"  hmmsearch outputs          : {n_with_hits}\n"
        f"  farlip_21.hmm present      : {'yes' if hmm_path.exists() else 'no'}\n"
        f"  synteny report present     : {'yes' if report_path.exists() else 'no'}\n"
        f"  Obsidian vault populated   : "
        f"{'yes' if vault_dir.exists() and any(vault_dir.iterdir()) else 'no'}\n"
    )
    return 0


def _build_parser() -> argparse.ArgumentParser:
    """Construct the argparse tree for ``python -m src.pipeline_phase3``."""
    p = argparse.ArgumentParser(
        prog="src.pipeline_phase3",
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
    p.add_argument(
        "--config-dir",
        type=Path,
        default=DEFAULT_CONFIG_DIR,
        help=f"Config directory (default: {DEFAULT_CONFIG_DIR})",
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
        help="Force re-run even when outputs already exist.",
    )

    sub.add_parser(
        "orf-predict",
        parents=[common_overwrite, common_filter],
        help="Run pyrodigal meta-mode ORF prediction for every assembly row.",
    )
    sub.add_parser(
        "build-profiles",
        parents=[common_overwrite],
        help="Build the 21 FaRLiP HMM profiles from JSC-1.",
    )
    sub.add_parser(
        "search",
        parents=[common_overwrite, common_filter],
        help="Run hmmsearch for every per-MAG protein set.",
    )
    sub.add_parser(
        "parse-synteny",
        parents=[common_overwrite],
        help="Classify hits with the FaRLiP synteny rule and write reports.",
    )
    sub.add_parser(
        "export-vault",
        parents=[common_overwrite],
        help="Render the Obsidian vault from the synteny reports.",
    )
    sub.add_parser(
        "all",
        parents=[common_overwrite],
        help="Run orf-predict -> build-profiles -> search -> parse-synteny -> export-vault.",
    )
    sub.add_parser(
        "status",
        help="Print a one-screen summary of Phase 3 progress.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    """Entry point for ``python -m src.pipeline_phase3``.

    Args:
        argv: Override for sys.argv[1:] (used by tests).

    Returns:
        Process exit code.
    """
    args = _build_parser().parse_args(argv)
    manifest_path: Path = args.manifest
    data_root: Path = args.data_root
    config_dir: Path = args.config_dir

    cmd = args.cmd
    accession_filter = getattr(args, "accession", None) or None
    overwrite = bool(getattr(args, "overwrite", False))

    if cmd == "orf-predict":
        return cmd_orf_predict(
            manifest_path, data_root, overwrite=overwrite, accession_filter=accession_filter
        )
    if cmd == "build-profiles":
        return cmd_build_profiles(data_root, config_dir, overwrite=overwrite)
    if cmd == "search":
        return cmd_search(
            manifest_path, data_root, overwrite=overwrite, accession_filter=accession_filter
        )
    if cmd == "parse-synteny":
        return cmd_parse_synteny(manifest_path, data_root, overwrite=overwrite)
    if cmd == "export-vault":
        return cmd_export_vault(manifest_path, data_root, overwrite=overwrite)
    if cmd == "all":
        for step in (
            lambda: cmd_orf_predict(manifest_path, data_root, overwrite=overwrite, accession_filter=None),
            lambda: cmd_build_profiles(data_root, config_dir, overwrite=overwrite),
            lambda: cmd_search(manifest_path, data_root, overwrite=overwrite, accession_filter=None),
            lambda: cmd_parse_synteny(manifest_path, data_root, overwrite=overwrite),
            lambda: cmd_export_vault(manifest_path, data_root, overwrite=overwrite),
        ):
            rc = step()
            if rc != 0:
                return rc
        return 0
    if cmd == "status":
        return cmd_status(manifest_path, data_root)
    sys.stderr.write(f"unknown command: {cmd}\n")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
