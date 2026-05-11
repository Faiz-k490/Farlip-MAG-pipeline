"""Tests for ``src.checkm2_ingestor``.

We materialize a synthetic CheckM2 ``quality_report.tsv`` with the canonical
column names and confirm the normalized output matches the
:data:`src.qc_assessor.QC_TSV_COLUMNS` schema exactly.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from src.checkm2_ingestor import (
    CHECKM2_MARKER_SET_LABEL,
    CHECKM2_QUALITY_REPORT_FILENAME,
    QC_SOURCE_CHECKM2,
    ingest_all_checkm2,
    ingest_checkm2_report,
)
from src.qc_assessor import QC_QUALITY_HIGH, QC_QUALITY_LOW, QC_QUALITY_MEDIUM, QC_TSV_COLUMNS


_CHECKM2_HEADER = (
    "Name\tCompleteness\tContamination\tCompleteness_Model_Used\t"
    "Translation_Table_Used\tCoding_Density\tContig_N50\tAverage_Gene_Length\t"
    "Genome_Size\tGC_Content\tTotal_Coding_Sequences\n"
)


def _write_quality_report(
    path: Path,
    name: str,
    completeness: str,
    contamination: str,
    n50: str,
    genome_size: str,
) -> None:
    """Write a one-row CheckM2 quality_report.tsv at the given path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    body = _CHECKM2_HEADER + (
        f"{name}\t{completeness}\t{contamination}\tNeural Network\t11\t0.92\t"
        f"{n50}\t302\t{genome_size}\t40.5\t3850\n"
    )
    path.write_text(body, encoding="utf-8")


def test_ingest_high_quality(tmp_path: Path):
    """A 95/2 row maps to MIMAG HIGH."""
    src = tmp_path / "GCA_TEST" / CHECKM2_QUALITY_REPORT_FILENAME
    _write_quality_report(
        src,
        name="GCA_TEST.faa",
        completeness="95.4",
        contamination="1.2",
        n50="184050",
        genome_size="4923118",
    )
    out_dir = tmp_path / "ingest"
    result = ingest_checkm2_report(
        accession="GCA_TEST",
        quality_report_path=src,
        output_dir=out_dir,
    )
    assert result.qc_quality == QC_QUALITY_HIGH
    assert result.completeness == pytest.approx(95.4)
    assert result.contamination == pytest.approx(1.2)
    assert result.n50_contig == 184050
    assert result.genome_size_bp == 4923118

    with result.normalized_path.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    for col in QC_TSV_COLUMNS:
        assert col in row
    assert row["accession"] == "GCA_TEST"
    assert row["qc_source"] == QC_SOURCE_CHECKM2
    assert row["marker_set_label"] == CHECKM2_MARKER_SET_LABEL


def test_ingest_medium_quality(tmp_path: Path):
    """60/8 maps to MEDIUM."""
    src = tmp_path / "GCA_MED" / CHECKM2_QUALITY_REPORT_FILENAME
    _write_quality_report(
        src,
        name="GCA_MED.faa",
        completeness="60.0",
        contamination="8.0",
        n50="50000",
        genome_size="3000000",
    )
    result = ingest_checkm2_report("GCA_MED", src, tmp_path / "out")
    assert result.qc_quality == QC_QUALITY_MEDIUM


def test_ingest_low_quality(tmp_path: Path):
    """30/3 maps to LOW because completeness fails MEDIUM."""
    src = tmp_path / "GCA_LOW" / CHECKM2_QUALITY_REPORT_FILENAME
    _write_quality_report(
        src,
        name="GCA_LOW.faa",
        completeness="30.0",
        contamination="3.0",
        n50="10000",
        genome_size="2000000",
    )
    result = ingest_checkm2_report("GCA_LOW", src, tmp_path / "out")
    assert result.qc_quality == QC_QUALITY_LOW


def test_ingest_tolerates_na_cells(tmp_path: Path):
    """``NA`` and empty cells parse to zero rather than exploding."""
    src = tmp_path / "GCA_NA" / CHECKM2_QUALITY_REPORT_FILENAME
    _write_quality_report(
        src,
        name="GCA_NA.faa",
        completeness="NA",
        contamination="",
        n50="NA",
        genome_size="",
    )
    result = ingest_checkm2_report("GCA_NA", src, tmp_path / "out")
    assert result.completeness == 0.0
    assert result.contamination == 0.0
    assert result.n50_contig == 0
    assert result.genome_size_bp == 0
    assert result.qc_quality == QC_QUALITY_LOW


def test_ingest_missing_file_raises(tmp_path: Path):
    """An absent file raises ``FileNotFoundError`` up-front."""
    with pytest.raises(FileNotFoundError):
        ingest_checkm2_report(
            "X",
            tmp_path / "missing.tsv",
            tmp_path / "out",
        )


def test_ingest_all_skips_dirs_without_report(tmp_path: Path):
    """Sub-directories that lack ``quality_report.tsv`` are silently skipped."""
    results_root = tmp_path / "checkm2_results"
    # GCA_OK has a report; GCA_PENDING does not.
    _write_quality_report(
        results_root / "GCA_OK" / CHECKM2_QUALITY_REPORT_FILENAME,
        name="GCA_OK.faa",
        completeness="92",
        contamination="1",
        n50="100000",
        genome_size="3000000",
    )
    (results_root / "GCA_PENDING").mkdir(parents=True)

    out_dir = tmp_path / "ingest"
    ingested = ingest_all_checkm2(results_root=results_root, output_dir=out_dir)
    accessions = {i.accession for i in ingested}
    assert accessions == {"GCA_OK"}


def test_ingest_all_filter(tmp_path: Path):
    """The ``accession_filter`` argument is honored."""
    results_root = tmp_path / "checkm2_results"
    for acc in ("GCA_A", "GCA_B"):
        _write_quality_report(
            results_root / acc / CHECKM2_QUALITY_REPORT_FILENAME,
            name=f"{acc}.faa",
            completeness="80",
            contamination="2",
            n50="100000",
            genome_size="2500000",
        )
    out_dir = tmp_path / "ingest"
    ingested = ingest_all_checkm2(
        results_root=results_root,
        output_dir=out_dir,
        accession_filter=["GCA_A"],
    )
    accessions = {i.accession for i in ingested}
    assert accessions == {"GCA_A"}
