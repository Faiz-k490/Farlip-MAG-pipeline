"""CLI orchestrator for the genome-quality assessment pipeline (Phase 3).

Reads the Phase 1 manifest (``data/manifest.tsv``) and the Phase 2 FaRLiP
roll-up (``data/farlip_hits_per_mag.tsv``) and runs - in order, or as
individually-invoked subcommands - the QC pipeline:

1. **fetch-markers** - :func:`src.qc_markers_fetcher.fetch_default_marker_set`
   downloads the configured single-copy marker HMM set into
   ``data/qc_markers/<label>.hmm``.
2. **assess** - :func:`src.qc_assessor.assess_assembly` for every assembly row,
   writing one ``data/qc_assessor/<accession>.qc.tsv`` per MAG.
3. **gen-checkm2-script** - :func:`src.checkm2_shellscript_generator.generate_checkm2_script`
   emits ``scripts/run_checkm2.sh`` for opt-in cross-validation.
4. **ingest-checkm2** - :func:`src.checkm2_ingestor.ingest_all_checkm2`
   normalizes any ``data/checkm2_results/<accession>/quality_report.tsv``
   files into ``data/checkm2_ingest/<accession>.qc.tsv``.
5. **report** - :func:`src.qc_reporter.write_qc_reports` joins the QC sources
   with the manifest and FaRLiP roll-up, writing ``data/qc_report.tsv`` and
   ``data/qc_filtered_mags.tsv``.
6. **update-vault** - re-runs :func:`src.obsidian_exporter.export_vault` with
   the QC report path so per-MAG markdown gets the QC frontmatter and the
   master index gets the Quality subsection.

Every subcommand is idempotent. Phase 3 progress is tracked in
``data/phase3_manifest.tsv`` (one row per ``(accession, stage)``) so a partial
``all`` run can resume from where it stopped.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from .checkm2_ingestor import (
    CHECKM2_INGEST_DIRNAME,
    CHECKM2_RESULTS_DIRNAME,
    ingest_all_checkm2,
)
from .checkm2_shellscript_generator import (
    CHECKM2_SCRIPT_RELPATH,
    CHECKM2_THREADS_DEFAULT,
    generate_checkm2_script,
)
from .extremophile_vault_exporter import (
    VAULT_DIRNAME as EXTREMOPHILE_VAULT_DIRNAME,
    export_extremophile_vault,
)
from .manifest import read_manifest
from .obsidian_exporter import export_vault
from .qc_assessor import assess_assembly
from .qc_markers_fetcher import (
    QC_MARKERS_DIR_NAME,
    fetch_default_marker_set,
    resolve_default_spec,
)
from .qc_reporter import write_qc_reports
from .ssu_bundler import ROUTING_LOG_FILENAME


REPO_ROOT: Path = Path(__file__).resolve().parent.parent
"""Repository root, derived from this module's location."""

DEFAULT_DATA_ROOT: Path = REPO_ROOT / "data"
"""Root of all pipeline outputs. Phases 1 and 2 already wrote into this directory."""

DEFAULT_MANIFEST: Path = DEFAULT_DATA_ROOT / "manifest.tsv"
"""Phase 1 manifest, consumed read-only by Phase 3."""

DEFAULT_CONFIG_DIR: Path = REPO_ROOT / "config"
"""Directory containing ``qc_markers.yaml``."""

DEFAULT_QC_MARKERS_CONFIG: Path = DEFAULT_CONFIG_DIR / "qc_markers.yaml"
"""Default location of the QC marker-set config."""

DEFAULT_SCRIPTS_DIR: Path = REPO_ROOT / "scripts"
"""Where the generated CheckM2 shell script lands by default."""

QC_ASSESSOR_DIRNAME: str = "qc_assessor"
"""Subdirectory under ``data/`` for per-MAG pyhmmer QC TSVs."""

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

NON_ASSEMBLY_SOURCES: frozenset[str] = frozenset({"SILVA", "GTDB", "QC_MARKERS"})
"""Manifest sources whose rows are skipped during QC assessment."""


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
    """Return a UTC ISO-8601 timestamp matching the Phase 1/2 manifests."""
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _read_phase3_manifest(path: Path) -> list[dict[str, str]]:
    """Load ``data/phase3_manifest.tsv`` if it exists."""
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_phase3_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    """Write a Phase 3 manifest TSV with the canonical column order."""
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
    """Merge new rows into the Phase 3 manifest, replacing on (accession, stage)."""
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


def _iter_assembly_rows(
    manifest_path: Path,
    accession_filter: list[str] | None = None,
) -> list[dict[str, str]]:
    """Yield every Phase 1 manifest row eligible for QC.

    Args:
        manifest_path: Phase 1 manifest path.
        accession_filter: If non-empty, restrict to these accessions.

    Returns:
        Eligible rows in manifest order.
    """
    rows = read_manifest(manifest_path)
    out: list[dict[str, str]] = []
    filter_set = set(accession_filter or [])
    for row in rows:
        if row.get("status") not in {"ok", "skipped"}:
            continue
        if row.get("source", "") in NON_ASSEMBLY_SOURCES:
            continue
        if filter_set and row.get("accession", "") not in filter_set:
            continue
        out.append(row)
    return out


def cmd_fetch_markers(
    config_dir: Path,
    data_root: Path,
    *,
    overwrite: bool,
    manifest_path: Path,
) -> int:
    """Download the configured default marker set and concatenate it.

    Args:
        config_dir: Config directory containing ``qc_markers.yaml``.
        data_root: Root data directory (output lands under ``data/qc_markers/``).
        overwrite: If True, re-download even when the stamp matches.
        manifest_path: Phase 1 manifest path; the marker artifact gets a row.

    Returns:
        Exit code (0 on success, 2 on configuration / download failure).
    """
    config_path = config_dir / "qc_markers.yaml"
    phase3_manifest = data_root / "phase3_manifest.tsv"
    if not config_path.exists():
        sys.stderr.write(f"[fetch-markers] missing {config_path}\n")
        return 2
    spec = resolve_default_spec(config_path)
    sys.stderr.write(
        f"[fetch-markers] default marker set = {spec.label} ({spec.archive_format})\n"
    )
    t0 = time.perf_counter()
    try:
        result = fetch_default_marker_set(
            config_path,
            data_root,
            overwrite=overwrite,
            manifest_path=manifest_path,
        )
    except Exception as exc:
        sys.stderr.write(f"[fetch-markers] FAILED: {exc}\n")
        _upsert_phase3_rows(
            phase3_manifest,
            [
                _Phase3Row(
                    accession=spec.label,
                    source="QC_MARKERS",
                    stage="fetch-markers",
                    status="failed",
                    output_path="",
                    n_units=0,
                    elapsed_s=time.perf_counter() - t0,
                    completed_at=_utc_now(),
                    note=str(exc),
                )
            ],
        )
        return 2
    elapsed = time.perf_counter() - t0
    marker = "=" if result.skipped else "+"
    sys.stderr.write(
        f"  {marker} {result.label}: n_markers={result.n_markers} "
        f"-> {result.hmm_db_path} ({elapsed:.1f}s)\n"
    )
    _upsert_phase3_rows(
        phase3_manifest,
        [
            _Phase3Row(
                accession=result.label,
                source="QC_MARKERS",
                stage="fetch-markers",
                status="skipped" if result.skipped else "ok",
                output_path=str(result.hmm_db_path),
                n_units=result.n_markers,
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note=f"sha256={result.archive_sha256[:16]}" if result.archive_sha256 else "",
            )
        ],
    )
    return 0


def cmd_assess(
    manifest_path: Path,
    data_root: Path,
    config_dir: Path,
    *,
    overwrite: bool,
    accession_filter: list[str] | None,
) -> int:
    """Run the in-process QC assessor for every eligible assembly.

    Args:
        manifest_path: Phase 1 manifest path.
        data_root: Root data directory.
        config_dir: Config directory containing ``qc_markers.yaml``.
        overwrite: If True, re-run even when the stamp matches.
        accession_filter: If non-empty, restrict to these accessions.

    Returns:
        Exit code (0 on success, 1 if any assembly failed, 2 on missing
        marker DB).
    """
    config_path = config_dir / "qc_markers.yaml"
    spec = resolve_default_spec(config_path)
    hmm_db_path = data_root / QC_MARKERS_DIR_NAME / f"{spec.key}.hmm"
    if not hmm_db_path.exists():
        sys.stderr.write(
            f"[assess] marker HMM database missing - run fetch-markers first ({hmm_db_path})\n"
        )
        return 2

    scores_cutoff_path: Path | None = (
        data_root / QC_MARKERS_DIR_NAME / f"{spec.key}.scores_cutoff.tsv"
    )
    if not scores_cutoff_path.exists():
        scores_cutoff_path = None

    rows = _iter_assembly_rows(manifest_path, accession_filter)
    orfs_dir = data_root / "orfs"
    output_dir = data_root / QC_ASSESSOR_DIRNAME
    phase3_manifest = data_root / "phase3_manifest.tsv"
    cutoff_note = (
        "with per-HMM bit-score cutoffs"
        if scores_cutoff_path is not None
        else "E-value-only mode (no scores_cutoff sidecar)"
    )
    sys.stderr.write(
        f"[assess] {len(rows)} assembly row(s) eligible "
        f"(markers={spec.label}, {cutoff_note})\n"
    )
    new_rows: list[_Phase3Row] = []
    failures = 0
    for row in rows:
        accession = row["accession"]
        faa_path = orfs_dir / f"{accession}.faa"
        coords_path = orfs_dir / f"{accession}.coords.tsv"
        if not faa_path.exists() or not coords_path.exists():
            sys.stderr.write(
                f"  ! {accession}: missing .faa or .coords.tsv - run "
                f"`python -m src.pipeline_phase2 orf-predict` first\n"
            )
            failures += 1
            new_rows.append(
                _Phase3Row(
                    accession=accession,
                    source=row.get("source", ""),
                    stage="assess",
                    status="failed",
                    output_path="",
                    n_units=0,
                    elapsed_s=0.0,
                    completed_at=_utc_now(),
                    note="missing .faa or .coords.tsv",
                )
            )
            continue
        t0 = time.perf_counter()
        try:
            result = assess_assembly(
                accession=accession,
                faa_path=faa_path,
                coords_path=coords_path,
                hmm_db_path=hmm_db_path,
                output_dir=output_dir,
                marker_set_label=spec.label,
                overwrite=overwrite,
                scores_cutoff_path=scores_cutoff_path,
            )
        except Exception as exc:
            failures += 1
            elapsed = time.perf_counter() - t0
            sys.stderr.write(f"  ! {accession}: {exc}\n")
            new_rows.append(
                _Phase3Row(
                    accession=accession,
                    source=row.get("source", ""),
                    stage="assess",
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
            f"  {marker} {accession}: completeness={result.completeness:.1f} "
            f"contamination={result.contamination:.1f} N50={result.n50_contig} "
            f"quality={result.qc_quality} ({elapsed:.1f}s)\n"
        )
        new_rows.append(
            _Phase3Row(
                accession=accession,
                source=row.get("source", ""),
                stage="assess",
                status="skipped" if result.skipped else "ok",
                output_path=str(result.qc_tsv_path),
                n_units=result.n_markers_found,
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note=(
                    f"completeness={result.completeness:.2f}; "
                    f"contamination={result.contamination:.2f}; "
                    f"qc_quality={result.qc_quality}"
                ),
            )
        )
    _upsert_phase3_rows(phase3_manifest, new_rows)
    return 1 if failures else 0


def cmd_gen_checkm2_script(
    manifest_path: Path,
    data_root: Path,
    *,
    threads: int,
) -> int:
    """Render ``scripts/run_checkm2.sh`` for opt-in CheckM2 cross-validation."""
    script_path = REPO_ROOT / CHECKM2_SCRIPT_RELPATH
    phase3_manifest = data_root / "phase3_manifest.tsv"
    t0 = time.perf_counter()
    result = generate_checkm2_script(
        manifest_path=manifest_path,
        data_root=data_root,
        script_path=script_path,
        threads=threads,
    )
    elapsed = time.perf_counter() - t0
    sys.stderr.write(
        f"[gen-checkm2-script] wrote {result.script_path} with "
        f"{result.n_assemblies} assemblies ({elapsed:.1f}s)\n"
    )
    if result.skipped_accessions:
        sys.stderr.write(
            f"  ! {len(result.skipped_accessions)} assemblies were skipped "
            f"(no .faa yet): {', '.join(result.skipped_accessions[:5])}"
            f"{'...' if len(result.skipped_accessions) > 5 else ''}\n"
        )
    _upsert_phase3_rows(
        phase3_manifest,
        [
            _Phase3Row(
                accession="*",
                source="checkm2",
                stage="gen-checkm2-script",
                status="ok",
                output_path=str(result.script_path),
                n_units=result.n_assemblies,
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note=f"skipped={len(result.skipped_accessions)}",
            )
        ],
    )
    return 0


def cmd_ingest_checkm2(
    data_root: Path,
    *,
    accession_filter: list[str] | None,
) -> int:
    """Normalize CheckM2 ``quality_report.tsv`` files into the Phase 3 schema."""
    phase3_manifest = data_root / "phase3_manifest.tsv"
    results_root = data_root / CHECKM2_RESULTS_DIRNAME
    output_dir = data_root / CHECKM2_INGEST_DIRNAME

    t0 = time.perf_counter()
    ingested = ingest_all_checkm2(
        results_root=results_root,
        output_dir=output_dir,
        accession_filter=accession_filter,
    )
    elapsed = time.perf_counter() - t0
    sys.stderr.write(
        f"[ingest-checkm2] ingested {len(ingested)} report(s) from {results_root} "
        f"({elapsed:.1f}s)\n"
    )
    new_rows: list[_Phase3Row] = []
    for r in ingested:
        new_rows.append(
            _Phase3Row(
                accession=r.accession,
                source="checkm2",
                stage="ingest-checkm2",
                status="ok",
                output_path=str(r.normalized_path),
                n_units=1,
                elapsed_s=0.0,
                completed_at=_utc_now(),
                note=(
                    f"completeness={r.completeness:.2f}; "
                    f"contamination={r.contamination:.2f}; "
                    f"qc_quality={r.qc_quality}"
                ),
            )
        )
    if new_rows:
        _upsert_phase3_rows(phase3_manifest, new_rows)
    return 0


def cmd_report(
    manifest_path: Path,
    data_root: Path,
) -> int:
    """Build ``qc_report.tsv`` and ``qc_filtered_mags.tsv``."""
    phase3_manifest = data_root / "phase3_manifest.tsv"
    pyhmmer_dir = data_root / QC_ASSESSOR_DIRNAME
    checkm2_dir = data_root / CHECKM2_INGEST_DIRNAME
    farlip_per_mag_path = data_root / "farlip_hits_per_mag.tsv"
    report_path = data_root / "qc_report.tsv"
    filtered_path = data_root / "qc_filtered_mags.tsv"

    t0 = time.perf_counter()
    n_total, n_filtered = write_qc_reports(
        manifest_path=manifest_path,
        pyhmmer_qc_dir=pyhmmer_dir,
        checkm2_ingest_dir=checkm2_dir,
        farlip_per_mag_path=farlip_per_mag_path,
        report_path=report_path,
        filtered_path=filtered_path,
    )
    elapsed = time.perf_counter() - t0
    sys.stderr.write(
        f"[report] wrote {report_path} ({n_total} rows) and {filtered_path} "
        f"({n_filtered} rows) ({elapsed:.1f}s)\n"
    )
    _upsert_phase3_rows(
        phase3_manifest,
        [
            _Phase3Row(
                accession="*",
                source="report",
                stage="report",
                status="ok",
                output_path=str(report_path),
                n_units=n_total,
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note=f"recommended_for_phylogeny={n_filtered}",
            )
        ],
    )
    return 0


def cmd_update_vault(
    manifest_path: Path,
    data_root: Path,
) -> int:
    """Re-run the Obsidian export with the QC report path supplied."""
    phase3_manifest = data_root / "phase3_manifest.tsv"
    report_path = data_root / "farlip_hits_report.tsv"
    per_mag_path = data_root / "farlip_hits_per_mag.tsv"
    qc_report_path = data_root / "qc_report.tsv"
    vault_dir = data_root / "obsidian_export"

    if not report_path.exists() or not per_mag_path.exists():
        sys.stderr.write(
            "[update-vault] FaRLiP synteny reports missing - run "
            "`python -m src.pipeline_phase2 parse-synteny` first\n"
        )
        return 2
    if not qc_report_path.exists():
        sys.stderr.write(
            "[update-vault] qc_report.tsv missing - run report first\n"
        )
        return 2

    t0 = time.perf_counter()
    n_written = export_vault(
        report_path=report_path,
        per_mag_path=per_mag_path,
        manifest_path=manifest_path,
        vault_dir=vault_dir,
        qc_report_path=qc_report_path,
    )
    elapsed = time.perf_counter() - t0
    sys.stderr.write(
        f"[update-vault] wrote {n_written} per-MAG markdown files (QC-aware) "
        f"to {vault_dir} ({elapsed:.1f}s)\n"
    )
    _upsert_phase3_rows(
        phase3_manifest,
        [
            _Phase3Row(
                accession="*",
                source="export",
                stage="update-vault",
                status="ok",
                output_path=str(vault_dir),
                n_units=n_written,
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note="qc-aware",
            )
        ],
    )
    return 0


def cmd_export_extremophile_vault(
    data_root: Path,
    *,
    routing_log_path: Path | None = None,
) -> int:
    """Render the high-quality extremophile Obsidian vault.

    Reads ``data/qc_filtered_mags.tsv`` (the Phase 3
    ``recommended_for_phylogeny`` view) plus, if present, the Phase 1.5
    routing log ``data/ssu/routing_log.tsv`` and writes a per-MAG markdown
    note under ``data/extremophile_vault/`` plus the
    ``High_Quality_Extremophile_Genomes.md`` master index.

    Idempotent (per-MAG files are overwritten).

    Args:
        data_root: Root data directory.
        routing_log_path: Override for the Phase 1.5 routing log path.
            Defaults to ``data_root/ssu/routing_log.tsv``; when the file is
            missing, ``Taxonomic_Placement`` slots are left empty.

    Returns:
        Exit code (0 on success, 2 when the filtered TSV is missing).
    """
    phase3_manifest = data_root / "phase3_manifest.tsv"
    filtered_path = data_root / "qc_filtered_mags.tsv"
    vault_dir = data_root / EXTREMOPHILE_VAULT_DIRNAME

    if not filtered_path.exists():
        sys.stderr.write(
            f"[export-extremophile-vault] {filtered_path} missing - run report first.\n"
        )
        return 2

    if routing_log_path is None:
        routing_log_path = data_root / "ssu" / ROUTING_LOG_FILENAME
    routing_arg = routing_log_path if routing_log_path.exists() else None
    if routing_arg is None:
        sys.stderr.write(
            f"[export-extremophile-vault] routing log {routing_log_path} not present; "
            "Taxonomic_Placement frontmatter will be left blank.\n"
        )

    t0 = time.perf_counter()
    n_written = export_extremophile_vault(
        filtered_csv_path=filtered_path,
        vault_dir=vault_dir,
        routing_log_path=routing_arg,
    )
    elapsed = time.perf_counter() - t0
    sys.stderr.write(
        f"[export-extremophile-vault] wrote {n_written} per-MAG markdown file(s) "
        f"to {vault_dir} ({elapsed:.1f}s)\n"
    )
    _upsert_phase3_rows(
        phase3_manifest,
        [
            _Phase3Row(
                accession="*",
                source="export",
                stage="export-extremophile-vault",
                status="ok",
                output_path=str(vault_dir),
                n_units=n_written,
                elapsed_s=elapsed,
                completed_at=_utc_now(),
                note=(
                    "routing-aware"
                    if routing_arg is not None
                    else "placement placeholders empty"
                ),
            )
        ],
    )
    return 0


def cmd_status(manifest_path: Path, data_root: Path) -> int:
    """Print a one-screen summary of Phase 3 state."""
    rows = _iter_assembly_rows(manifest_path, None)
    qc_dir = data_root / QC_ASSESSOR_DIRNAME
    ingest_dir = data_root / CHECKM2_INGEST_DIRNAME
    qc_markers_dir = data_root / QC_MARKERS_DIR_NAME

    n_pyhmmer = sum(1 for r in rows if (qc_dir / f"{r['accession']}.qc.tsv").exists())
    n_checkm2 = sum(
        1 for r in rows if (ingest_dir / f"{r['accession']}.qc.tsv").exists()
    )
    has_markers = any(qc_markers_dir.glob("*.hmm")) if qc_markers_dir.exists() else False
    report_path = data_root / "qc_report.tsv"
    filtered_path = data_root / "qc_filtered_mags.tsv"
    extremophile_dir = data_root / EXTREMOPHILE_VAULT_DIRNAME
    extremophile_index = (
        extremophile_dir / "High_Quality_Extremophile_Genomes.md"
    )

    sys.stdout.write(
        "FaRLiP Phase 3 status\n"
        f"  assembly rows in manifest    : {len(rows)}\n"
        f"  marker HMM database present  : {'yes' if has_markers else 'no'}\n"
        f"  pyhmmer QC TSVs              : {n_pyhmmer}\n"
        f"  CheckM2 QC TSVs (ingested)   : {n_checkm2}\n"
        f"  qc_report.tsv present        : {'yes' if report_path.exists() else 'no'}\n"
        f"  qc_filtered_mags.tsv present : {'yes' if filtered_path.exists() else 'no'}\n"
        f"  extremophile vault index     : "
        f"{'yes' if extremophile_index.exists() else 'no'}\n"
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
        "fetch-markers",
        parents=[common_overwrite],
        help="Download the configured single-copy marker HMM set.",
    )
    sub.add_parser(
        "assess",
        parents=[common_overwrite, common_filter],
        help="Run the in-process pyhmmer QC assessor for every assembly.",
    )
    gen_p = sub.add_parser(
        "gen-checkm2-script",
        help="Generate scripts/run_checkm2.sh for opt-in CheckM2 cross-validation.",
    )
    gen_p.add_argument(
        "--threads",
        type=int,
        default=CHECKM2_THREADS_DEFAULT,
        help=f"Threads for each `checkm2 predict` call (default: {CHECKM2_THREADS_DEFAULT}).",
    )
    sub.add_parser(
        "ingest-checkm2",
        parents=[common_filter],
        help="Normalize CheckM2 quality_report.tsv outputs into the QC schema.",
    )
    sub.add_parser(
        "report",
        help="Join QC sources with the manifest + FaRLiP roll-up into qc_report.tsv.",
    )
    sub.add_parser(
        "update-vault",
        help="Re-render the Obsidian vault with QC frontmatter merged in.",
    )
    extremophile_p = sub.add_parser(
        "export-extremophile-vault",
        help=(
            "Render data/extremophile_vault/ "
            "(High_Quality_Extremophile_Genomes.md + per-MAG notes) from "
            "data/qc_filtered_mags.tsv."
        ),
    )
    extremophile_p.add_argument(
        "--routing-log",
        type=Path,
        default=None,
        help=(
            "Override path to the Phase 1.5 routing log (default: "
            "<data-root>/ssu/routing_log.tsv). When missing, "
            "Taxonomic_Placement is left blank."
        ),
    )
    sub.add_parser(
        "all",
        parents=[common_overwrite],
        help=(
            "Run fetch-markers -> assess -> gen-checkm2-script -> report -> "
            "update-vault -> export-extremophile-vault."
        ),
    )
    sub.add_parser(
        "status",
        help="Print a one-screen summary of Phase 3 progress.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    """Entry point for ``python -m src.pipeline_phase3``.

    Args:
        argv: Override for ``sys.argv[1:]`` (used by tests).

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
    threads = int(getattr(args, "threads", CHECKM2_THREADS_DEFAULT))

    if cmd == "fetch-markers":
        return cmd_fetch_markers(
            config_dir,
            data_root,
            overwrite=overwrite,
            manifest_path=manifest_path,
        )
    if cmd == "assess":
        return cmd_assess(
            manifest_path,
            data_root,
            config_dir,
            overwrite=overwrite,
            accession_filter=accession_filter,
        )
    if cmd == "gen-checkm2-script":
        return cmd_gen_checkm2_script(
            manifest_path,
            data_root,
            threads=threads,
        )
    if cmd == "ingest-checkm2":
        return cmd_ingest_checkm2(
            data_root,
            accession_filter=accession_filter,
        )
    if cmd == "report":
        return cmd_report(manifest_path, data_root)
    if cmd == "update-vault":
        return cmd_update_vault(manifest_path, data_root)
    if cmd == "export-extremophile-vault":
        return cmd_export_extremophile_vault(
            data_root,
            routing_log_path=getattr(args, "routing_log", None),
        )
    if cmd == "all":
        for step in (
            lambda: cmd_fetch_markers(
                config_dir, data_root, overwrite=overwrite, manifest_path=manifest_path
            ),
            lambda: cmd_assess(
                manifest_path,
                data_root,
                config_dir,
                overwrite=overwrite,
                accession_filter=None,
            ),
            lambda: cmd_gen_checkm2_script(
                manifest_path, data_root, threads=CHECKM2_THREADS_DEFAULT
            ),
            lambda: cmd_report(manifest_path, data_root),
            lambda: cmd_update_vault(manifest_path, data_root),
            lambda: cmd_export_extremophile_vault(data_root),
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
