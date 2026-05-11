"""Join per-MAG QC TSVs with the Phase 1 manifest and Phase 2 FaRLiP rollup.

Produces two files:

* ``data/qc_report.tsv`` - one row per MAG in the Phase 1 manifest, with the
  preferred QC source merged in (CheckM2 wins over pyhmmer when both are
  present for the same accession), the FaRLiP top classification copied over
  from ``data/farlip_hits_per_mag.tsv``, and a ``recommended_for_phylogeny``
  boolean computed as ``qc_quality in {HIGH, MEDIUM} AND
  farlip_top_classification != NONE``.
* ``data/qc_filtered_mags.tsv`` - subset of the above where
  ``recommended_for_phylogeny = true``. The filter is a sidecar view; no row is
  ever deleted from ``qc_report.tsv``.

The reporter does no biology of its own and never re-runs hmmsearch / CheckM2;
it is a pure TSV-join step.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .qc_assessor import (
    QC_QUALITY_HIGH,
    QC_QUALITY_LOW,
    QC_QUALITY_MEDIUM,
    QC_TSV_COLUMNS,
)


REPORT_COLUMNS: tuple[str, ...] = (
    "accession",
    "source",
    "scientific_name",
    "qc_source",
    "marker_set_label",
    "completeness",
    "contamination",
    "n50_contig",
    "genome_size_bp",
    "n_contigs",
    "longest_contig_bp",
    "qc_quality",
    "n_farlip_distinct_markers",
    "farlip_top_classification",
    "recommended_for_phylogeny",
)
"""Schema of ``data/qc_report.tsv`` and ``data/qc_filtered_mags.tsv``.

Stable column order. The Obsidian QC-aware exporter reads this file and assumes
these column names; downstream consumers must add new columns at the right.
"""

NON_ASSEMBLY_SOURCES: frozenset[str] = frozenset({"SILVA", "GTDB", "QC_MARKERS"})
"""Manifest sources that are not genome assemblies and must be excluded from
the per-MAG report. Mirrors the filter applied in
:mod:`src.checkm2_shellscript_generator`."""


_FARLIP_NONE_TAG: str = "NONE"
"""Top-classification value emitted by :mod:`src.synteny_parser` when no
candidate contig was retained."""


@dataclass(frozen=True)
class _QcRecord:
    """Per-accession QC data assembled before TSV serialization."""

    accession: str
    qc_source: str
    marker_set_label: str
    completeness: float
    contamination: float
    n50_contig: int
    genome_size_bp: int
    n_contigs: str
    longest_contig_bp: str
    qc_quality: str


_QC_SOURCE_PRIORITY: dict[str, int] = {
    "checkm2": 2,
    "pyhmmer": 1,
}
"""Higher = preferred when more than one QC source has a row for one MAG."""


def _read_qc_dir(qc_dir: Path) -> dict[str, _QcRecord]:
    """Load every ``<accession>.qc.tsv`` in a directory into a dict.

    Args:
        qc_dir: Directory containing per-MAG ``.qc.tsv`` files.

    Returns:
        Mapping of accession -> :class:`_QcRecord`. Empty if the directory
        does not exist or contains no TSVs.
    """
    out: dict[str, _QcRecord] = {}
    if not qc_dir.exists():
        return out
    for path in sorted(qc_dir.glob("*.qc.tsv")):
        with path.open("r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                accession = row.get("accession") or path.stem.replace(".qc", "")
                out[accession] = _QcRecord(
                    accession=accession,
                    qc_source=row.get("qc_source", ""),
                    marker_set_label=row.get("marker_set_label", ""),
                    completeness=_safe_float(row.get("completeness", "")),
                    contamination=_safe_float(row.get("contamination", "")),
                    n50_contig=_safe_int(row.get("n50_contig", "")),
                    genome_size_bp=_safe_int(row.get("genome_size_bp", "")),
                    n_contigs=row.get("n_contigs", ""),
                    longest_contig_bp=row.get("longest_contig_bp", ""),
                    qc_quality=row.get("qc_quality", QC_QUALITY_LOW),
                )
                break
    return out


def _safe_float(value: str) -> float:
    """Parse a numeric cell tolerating empty / missing values."""
    if not value:
        return 0.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def _safe_int(value: str) -> int:
    """Parse an integer cell tolerating empty / missing values."""
    if not value:
        return 0
    try:
        return int(round(float(value)))
    except ValueError:
        return 0


def _merge_sources(
    pyhmmer_rows: dict[str, _QcRecord],
    checkm2_rows: dict[str, _QcRecord],
) -> dict[str, _QcRecord]:
    """Merge the two QC sources, preferring the higher-priority source per MAG.

    Args:
        pyhmmer_rows: pyhmmer-derived QC rows by accession.
        checkm2_rows: CheckM2-derived QC rows by accession.

    Returns:
        One row per accession from whichever source has the higher priority in
        :data:`_QC_SOURCE_PRIORITY`. Identical priorities (shouldn't happen)
        favor pyhmmer for stability.
    """
    merged: dict[str, _QcRecord] = dict(pyhmmer_rows)
    for accession, record in checkm2_rows.items():
        existing = merged.get(accession)
        if existing is None:
            merged[accession] = record
            continue
        if (
            _QC_SOURCE_PRIORITY.get(record.qc_source, 0)
            > _QC_SOURCE_PRIORITY.get(existing.qc_source, 0)
        ):
            merged[accession] = record
    return merged


def _read_per_mag_farlip(per_mag_path: Path) -> dict[str, dict[str, str]]:
    """Load ``data/farlip_hits_per_mag.tsv`` indexed by accession."""
    out: dict[str, dict[str, str]] = {}
    if not per_mag_path.exists():
        return out
    with per_mag_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            accession = row.get("mag_accession", "")
            if accession:
                out[accession] = row
    return out


def _read_manifest_assemblies(manifest_path: Path) -> list[dict[str, str]]:
    """Read manifest rows that represent genome assemblies (filter non-assembly sources)."""
    if not manifest_path.exists():
        return []
    with manifest_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    out: list[dict[str, str]] = []
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
        seen.add(accession)
        out.append(row)
    return out


def _is_recommended(qc_quality: str, farlip_top: str) -> bool:
    """Encode the Phase 4 phylogeny gate.

    A MAG is recommended if its QC quality is at least MEDIUM **and** Phase 2
    flagged at least one FaRLiP candidate contig (i.e. ``farlip_top`` is not
    ``NONE``).
    """
    if qc_quality not in {QC_QUALITY_HIGH, QC_QUALITY_MEDIUM}:
        return False
    if not farlip_top or farlip_top == _FARLIP_NONE_TAG:
        return False
    return True


def write_qc_reports(
    manifest_path: Path,
    pyhmmer_qc_dir: Path,
    checkm2_ingest_dir: Path,
    farlip_per_mag_path: Path,
    report_path: Path,
    filtered_path: Path,
) -> tuple[int, int]:
    """Build ``qc_report.tsv`` and ``qc_filtered_mags.tsv``.

    Args:
        manifest_path: Path to ``data/manifest.tsv``.
        pyhmmer_qc_dir: Directory of ``<accession>.qc.tsv`` written by
            :mod:`src.qc_assessor` (``data/qc_assessor``).
        checkm2_ingest_dir: Directory of normalized ``<accession>.qc.tsv``
            written by :mod:`src.checkm2_ingestor` (``data/checkm2_ingest``).
        farlip_per_mag_path: Path to ``data/farlip_hits_per_mag.tsv``.
        report_path: Output path for the joined per-MAG report.
        filtered_path: Output path for the recommended-only subset.

    Returns:
        Tuple of ``(rows_in_report, rows_in_filtered)``.
    """
    manifest_rows = _read_manifest_assemblies(manifest_path)
    pyhmmer_rows = _read_qc_dir(pyhmmer_qc_dir)
    checkm2_rows = _read_qc_dir(checkm2_ingest_dir)
    qc_by_acc = _merge_sources(pyhmmer_rows, checkm2_rows)
    farlip_by_acc = _read_per_mag_farlip(farlip_per_mag_path)

    report_rows: list[dict[str, str]] = []
    filtered_rows: list[dict[str, str]] = []

    for row in manifest_rows:
        accession = row["accession"]
        qc = qc_by_acc.get(accession)
        farlip = farlip_by_acc.get(accession, {})
        farlip_top = farlip.get("top_classification", _FARLIP_NONE_TAG) or _FARLIP_NONE_TAG
        farlip_distinct = farlip.get("n_distinct_markers", "0") or "0"

        if qc is None:
            qc_record = _QcRecord(
                accession=accession,
                qc_source="",
                marker_set_label="",
                completeness=0.0,
                contamination=0.0,
                n50_contig=0,
                genome_size_bp=0,
                n_contigs="",
                longest_contig_bp="",
                qc_quality=QC_QUALITY_LOW,
            )
        else:
            qc_record = qc

        recommended = _is_recommended(qc_record.qc_quality, farlip_top)
        out_row = {
            "accession": accession,
            "source": row.get("source", ""),
            "scientific_name": row.get("scientific_name", ""),
            "qc_source": qc_record.qc_source,
            "marker_set_label": qc_record.marker_set_label,
            "completeness": f"{qc_record.completeness:.2f}",
            "contamination": f"{qc_record.contamination:.2f}",
            "n50_contig": str(qc_record.n50_contig),
            "genome_size_bp": str(qc_record.genome_size_bp),
            "n_contigs": qc_record.n_contigs,
            "longest_contig_bp": qc_record.longest_contig_bp,
            "qc_quality": qc_record.qc_quality,
            "n_farlip_distinct_markers": farlip_distinct,
            "farlip_top_classification": farlip_top,
            "recommended_for_phylogeny": "true" if recommended else "false",
        }
        report_rows.append(out_row)
        if recommended:
            filtered_rows.append(out_row)

    _write_report(report_path, report_rows)
    _write_report(filtered_path, filtered_rows)
    return len(report_rows), len(filtered_rows)


def _write_report(path: Path, rows: Iterable[dict[str, str]]) -> None:
    """Write a TSV using the canonical :data:`REPORT_COLUMNS` order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=REPORT_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in REPORT_COLUMNS})


__all__ = [
    "NON_ASSEMBLY_SOURCES",
    "REPORT_COLUMNS",
    "write_qc_reports",
]


# Re-export so external users do not need a second import to read the schema
# of the per-MAG TSVs that feed this module.
__qc_tsv_columns_passthrough__ = QC_TSV_COLUMNS
