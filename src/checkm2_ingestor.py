"""Normalize CheckM2 ``quality_report.tsv`` outputs into the Phase 4 QC schema.

CheckM2 (Chklovski et al. 2023, *Nat Methods* 20:1203) is invoked via the
shell script emitted by :mod:`src.checkm2_shellscript_generator`. Each
``checkm2 predict`` call writes per-MAG output under
``data/checkm2_results/<accession>/``, with the headline metrics in a
``quality_report.tsv`` whose canonical CheckM2 schema is::

    Name<TAB>Completeness<TAB>Contamination<TAB>Completeness_Model_Used<TAB>
    Translation_Table_Used<TAB>Coding_Density<TAB>Contig_N50<TAB>
    Average_Gene_Length<TAB>Genome_Size<TAB>GC_Content<TAB>...

This module reads those rows and re-emits one
``data/checkm2_ingest/<accession>.qc.tsv`` per assembly, using the same
:data:`src.qc_assessor.QC_TSV_COLUMNS` schema as the in-process pyhmmer
assessor so the downstream :mod:`src.qc_reporter` can join either source
without branching.

We do not run CheckM2 ourselves; we just translate its outputs. The module is
purely I/O-bound and trivially testable on a synthetic ``quality_report.tsv``.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .qc_assessor import (
    MIMAG_HIGH_QUALITY_COMPLETENESS,
    MIMAG_HIGH_QUALITY_CONTAMINATION,
    MIMAG_MEDIUM_QUALITY_COMPLETENESS,
    MIMAG_MEDIUM_QUALITY_CONTAMINATION,
    QC_QUALITY_HIGH,
    QC_QUALITY_LOW,
    QC_QUALITY_MEDIUM,
    QC_TSV_COLUMNS,
)


CHECKM2_QUALITY_REPORT_FILENAME: str = "quality_report.tsv"
"""CheckM2's headline output filename inside each per-MAG output dir."""

CHECKM2_RESULTS_DIRNAME: str = "checkm2_results"
"""Subdirectory of ``data/`` written by ``scripts/run_checkm2.sh``. Matches
:data:`src.checkm2_shellscript_generator.CHECKM2_RESULTS_DIRNAME`."""

CHECKM2_INGEST_DIRNAME: str = "checkm2_ingest"
"""Subdirectory of ``data/`` where this module writes the normalized
``.qc.tsv`` files."""

QC_SOURCE_CHECKM2: str = "checkm2"
"""``qc_source`` tag emitted by this module. The reporter prefers ``checkm2``
rows over ``pyhmmer`` rows when both are present for the same accession."""

CHECKM2_MARKER_SET_LABEL: str = "CheckM2_uniref100_diamond"
"""Marker-set label written into the normalized TSV. Differentiates CheckM2
output from BUSCO/GTDB pyhmmer outputs in the joined report."""


_CHECKM2_NAME_COL: str = "Name"
_CHECKM2_COMPLETENESS_COL: str = "Completeness"
_CHECKM2_CONTAMINATION_COL: str = "Contamination"
_CHECKM2_N50_COL: str = "Contig_N50"
_CHECKM2_GENOME_SIZE_COL: str = "Genome_Size"


@dataclass(frozen=True)
class CheckM2Ingested:
    """Outcome of one :func:`ingest_checkm2_report` invocation.

    Attributes:
        accession: MAG accession (taken from the directory name).
        source_path: ``quality_report.tsv`` that was parsed.
        normalized_path: Path to the written normalized ``.qc.tsv``.
        completeness: Completeness percent.
        contamination: Contamination percent.
        n50_contig: Contig N50 in bp.
        genome_size_bp: Genome size in bp.
        qc_quality: MIMAG bin (HIGH / MEDIUM / LOW).
    """

    accession: str
    source_path: Path
    normalized_path: Path
    completeness: float
    contamination: float
    n50_contig: int
    genome_size_bp: int
    qc_quality: str


def _classify_quality(completeness: float, contamination: float) -> str:
    """MIMAG ladder identical to :func:`src.qc_assessor._classify_quality`.

    Reproduced here (rather than imported) so this module's contract is
    self-contained for unit testing - the constants come from the assessor so
    a single source of truth is preserved.

    Args:
        completeness: 0-100.
        contamination: 0-100.

    Returns:
        :data:`src.qc_assessor.QC_QUALITY_HIGH`,
        :data:`src.qc_assessor.QC_QUALITY_MEDIUM`, or
        :data:`src.qc_assessor.QC_QUALITY_LOW`.
    """
    if (
        completeness >= MIMAG_HIGH_QUALITY_COMPLETENESS
        and contamination <= MIMAG_HIGH_QUALITY_CONTAMINATION
    ):
        return QC_QUALITY_HIGH
    if (
        completeness >= MIMAG_MEDIUM_QUALITY_COMPLETENESS
        and contamination < MIMAG_MEDIUM_QUALITY_CONTAMINATION
    ):
        return QC_QUALITY_MEDIUM
    return QC_QUALITY_LOW


def _read_checkm2_row(path: Path) -> dict[str, str] | None:
    """Read the (single) data row out of a CheckM2 ``quality_report.tsv``.

    CheckM2 emits one row per input genome; we always pass exactly one
    ``--input`` per call, so a single row is expected.

    Args:
        path: Path to the ``quality_report.tsv``.

    Returns:
        First data row as a dict, or ``None`` if the file is missing / empty.
    """
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            return dict(row)
    return None


def _to_float(value: str) -> float:
    """Parse a CheckM2 numeric cell, tolerating empty / ``NA`` values."""
    if value is None:
        return 0.0
    s = value.strip()
    if not s or s.upper() in {"NA", "N/A"}:
        return 0.0
    try:
        return float(s)
    except ValueError:
        return 0.0


def _to_int(value: str) -> int:
    """Parse a CheckM2 numeric cell as int (rounded), tolerating empty values."""
    f = _to_float(value)
    return int(round(f))


def ingest_checkm2_report(
    accession: str,
    quality_report_path: Path,
    output_dir: Path,
) -> CheckM2Ingested:
    """Translate one CheckM2 ``quality_report.tsv`` into the Phase 4 schema.

    Args:
        accession: MAG accession (drives the output filename).
        quality_report_path: Path to CheckM2's ``quality_report.tsv``.
        output_dir: Directory under which ``<accession>.qc.tsv`` is written.

    Returns:
        A :class:`CheckM2Ingested` describing what was written and computed.

    Raises:
        FileNotFoundError: If ``quality_report_path`` does not exist.
        RuntimeError: If the report exists but has no data rows.
    """
    if not quality_report_path.exists():
        raise FileNotFoundError(
            f"CheckM2 quality_report.tsv not found: {quality_report_path}"
        )
    row = _read_checkm2_row(quality_report_path)
    if row is None:
        raise RuntimeError(
            f"CheckM2 quality_report.tsv at {quality_report_path} has no data rows"
        )

    completeness = _to_float(row.get(_CHECKM2_COMPLETENESS_COL, ""))
    contamination = _to_float(row.get(_CHECKM2_CONTAMINATION_COL, ""))
    n50 = _to_int(row.get(_CHECKM2_N50_COL, ""))
    genome_size = _to_int(row.get(_CHECKM2_GENOME_SIZE_COL, ""))
    quality = _classify_quality(completeness, contamination)

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{accession}.qc.tsv"

    output_row = {
        "accession": accession,
        "qc_source": QC_SOURCE_CHECKM2,
        "marker_set_label": CHECKM2_MARKER_SET_LABEL,
        "n_markers_total": "",
        "n_markers_found": "",
        "n_markers_duplicated": "",
        "completeness": f"{completeness:.2f}",
        "contamination": f"{contamination:.2f}",
        "n50_contig": str(n50),
        "genome_size_bp": str(genome_size),
        "n_contigs": "",
        "longest_contig_bp": "",
        "qc_quality": quality,
    }
    with out_path.open("w", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=QC_TSV_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerow({col: output_row.get(col, "") for col in QC_TSV_COLUMNS})

    return CheckM2Ingested(
        accession=accession,
        source_path=quality_report_path,
        normalized_path=out_path,
        completeness=completeness,
        contamination=contamination,
        n50_contig=n50,
        genome_size_bp=genome_size,
        qc_quality=quality,
    )


def ingest_all_checkm2(
    results_root: Path,
    output_dir: Path,
    *,
    accession_filter: Iterable[str] | None = None,
) -> list[CheckM2Ingested]:
    """Walk ``data/checkm2_results/`` and ingest every per-MAG report present.

    Args:
        results_root: Directory containing ``<accession>/quality_report.tsv``
            subdirectories (``data/checkm2_results``).
        output_dir: Directory under which normalized ``.qc.tsv`` files are
            written (``data/checkm2_ingest``).
        accession_filter: Optional iterable of accessions to restrict to.

    Returns:
        One :class:`CheckM2Ingested` per ingested directory, in sorted
        accession order. Sub-directories that lack a ``quality_report.tsv``
        are silently skipped (the user may still be running CheckM2 on them).
    """
    if not results_root.exists():
        return []
    keep: set[str] | None = set(accession_filter) if accession_filter else None
    out: list[CheckM2Ingested] = []
    for child in sorted(results_root.iterdir()):
        if not child.is_dir():
            continue
        if keep is not None and child.name not in keep:
            continue
        report = child / CHECKM2_QUALITY_REPORT_FILENAME
        if not report.exists():
            continue
        out.append(
            ingest_checkm2_report(
                accession=child.name,
                quality_report_path=report,
                output_dir=output_dir,
            )
        )
    return out


__all__ = [
    "CHECKM2_INGEST_DIRNAME",
    "CHECKM2_MARKER_SET_LABEL",
    "CHECKM2_QUALITY_REPORT_FILENAME",
    "CHECKM2_RESULTS_DIRNAME",
    "CheckM2Ingested",
    "QC_SOURCE_CHECKM2",
    "ingest_all_checkm2",
    "ingest_checkm2_report",
]
