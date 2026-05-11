"""Tests for ``src.qc_reporter``.

The reporter is a pure TSV-join step. We materialize tiny manifest, per-MAG
QC, and FaRLiP roll-up TSVs and verify that the merged ``qc_report.tsv`` and
``qc_filtered_mags.tsv`` honour the source-priority and recommendation rules.
"""

from __future__ import annotations

import csv
from pathlib import Path

from src.qc_assessor import (
    QC_QUALITY_HIGH,
    QC_QUALITY_LOW,
    QC_QUALITY_MEDIUM,
    QC_TSV_COLUMNS,
)
from src.qc_reporter import REPORT_COLUMNS, write_qc_reports


_MANIFEST_HEADER = (
    "accession\tsource\turl\tlocal_path\tsha256\tbytes_downloaded\tfetched_at\t"
    "status\tnote\tlabel\tscientific_name\tsample_accession\tassembly_level\t"
    "format\tmd5\textracted_dir\n"
)
_PER_MAG_HEADER = (
    "mag_accession\tn_candidate_contigs\tn_total_markers\tn_distinct_markers\t"
    "has_edge_truncated\ttop_classification\n"
)


def _write_manifest(path: Path, rows: list[str]) -> None:
    path.write_text(_MANIFEST_HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def _write_per_mag(path: Path, rows: list[str]) -> None:
    path.write_text(_PER_MAG_HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def _write_qc_tsv(
    path: Path,
    *,
    accession: str,
    qc_source: str,
    completeness: float,
    contamination: float,
    quality: str,
    n50: int = 100000,
    genome_size: int = 3000000,
) -> None:
    """Write a single-row per-MAG QC TSV in the canonical schema."""
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "accession": accession,
        "qc_source": qc_source,
        "marker_set_label": qc_source.upper(),
        "n_markers_total": "120",
        "n_markers_found": str(int(120 * completeness / 100)),
        "n_markers_duplicated": "0",
        "completeness": f"{completeness:.2f}",
        "contamination": f"{contamination:.2f}",
        "n50_contig": str(n50),
        "genome_size_bp": str(genome_size),
        "n_contigs": "10",
        "longest_contig_bp": str(n50 * 2),
        "qc_quality": quality,
    }
    with path.open("w", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=QC_TSV_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerow({c: row.get(c, "") for c in QC_TSV_COLUMNS})


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_report_columns_stable():
    """Column order is part of the public contract."""
    assert REPORT_COLUMNS[0] == "accession"
    assert "completeness" in REPORT_COLUMNS
    assert "qc_quality" in REPORT_COLUMNS
    assert REPORT_COLUMNS[-1] == "recommended_for_phylogeny"


def test_report_joins_pyhmmer_and_farlip(tmp_path: Path):
    """A high-quality MAG with a FULL_CLUSTER call is recommended."""
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(
        manifest,
        [
            "GCA_HIT\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\tCandidatus Test\tSAM_X\t-\t-\t-\t-",
        ],
    )
    pyhmmer_dir = tmp_path / "qc_assessor"
    _write_qc_tsv(
        pyhmmer_dir / "GCA_HIT.qc.tsv",
        accession="GCA_HIT",
        qc_source="pyhmmer",
        completeness=95.0,
        contamination=1.0,
        quality=QC_QUALITY_HIGH,
    )
    per_mag = tmp_path / "farlip_hits_per_mag.tsv"
    _write_per_mag(per_mag, ["GCA_HIT\t1\t20\t20\tfalse\tFULL_CLUSTER"])

    report_path = tmp_path / "qc_report.tsv"
    filtered_path = tmp_path / "qc_filtered_mags.tsv"
    n_total, n_filtered = write_qc_reports(
        manifest_path=manifest,
        pyhmmer_qc_dir=pyhmmer_dir,
        checkm2_ingest_dir=tmp_path / "checkm2_ingest",
        farlip_per_mag_path=per_mag,
        report_path=report_path,
        filtered_path=filtered_path,
    )
    assert n_total == 1
    assert n_filtered == 1
    rows = _read_tsv_rows(report_path)
    assert rows[0]["accession"] == "GCA_HIT"
    assert rows[0]["qc_quality"] == QC_QUALITY_HIGH
    assert rows[0]["recommended_for_phylogeny"] == "true"
    assert rows[0]["farlip_top_classification"] == "FULL_CLUSTER"


def test_report_excludes_low_quality_from_filtered(tmp_path: Path):
    """LOW-quality MAGs land in the report but not the recommended subset."""
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(
        manifest,
        [
            "GCA_LOW\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\t-\t-\t-\t-\t-\t-",
        ],
    )
    pyhmmer_dir = tmp_path / "qc_assessor"
    _write_qc_tsv(
        pyhmmer_dir / "GCA_LOW.qc.tsv",
        accession="GCA_LOW",
        qc_source="pyhmmer",
        completeness=30.0,
        contamination=2.0,
        quality=QC_QUALITY_LOW,
    )
    per_mag = tmp_path / "farlip_hits_per_mag.tsv"
    _write_per_mag(per_mag, ["GCA_LOW\t1\t10\t10\tfalse\tPARTIAL_CLUSTER"])

    report_path = tmp_path / "qc_report.tsv"
    filtered_path = tmp_path / "qc_filtered_mags.tsv"
    n_total, n_filtered = write_qc_reports(
        manifest_path=manifest,
        pyhmmer_qc_dir=pyhmmer_dir,
        checkm2_ingest_dir=tmp_path / "checkm2_ingest",
        farlip_per_mag_path=per_mag,
        report_path=report_path,
        filtered_path=filtered_path,
    )
    assert n_total == 1
    assert n_filtered == 0
    rows = _read_tsv_rows(report_path)
    assert rows[0]["recommended_for_phylogeny"] == "false"


def test_report_excludes_no_farlip_from_filtered(tmp_path: Path):
    """A high-quality MAG without a FaRLiP candidate is not recommended."""
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(
        manifest,
        [
            "GCA_NOFAR\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\t-\t-\t-\t-\t-\t-",
        ],
    )
    pyhmmer_dir = tmp_path / "qc_assessor"
    _write_qc_tsv(
        pyhmmer_dir / "GCA_NOFAR.qc.tsv",
        accession="GCA_NOFAR",
        qc_source="pyhmmer",
        completeness=95.0,
        contamination=1.0,
        quality=QC_QUALITY_HIGH,
    )
    per_mag = tmp_path / "farlip_hits_per_mag.tsv"
    _write_per_mag(per_mag, ["GCA_NOFAR\t0\t0\t0\tfalse\tNONE"])

    report_path = tmp_path / "qc_report.tsv"
    filtered_path = tmp_path / "qc_filtered_mags.tsv"
    write_qc_reports(
        manifest_path=manifest,
        pyhmmer_qc_dir=pyhmmer_dir,
        checkm2_ingest_dir=tmp_path / "checkm2_ingest",
        farlip_per_mag_path=per_mag,
        report_path=report_path,
        filtered_path=filtered_path,
    )
    rows = _read_tsv_rows(report_path)
    assert rows[0]["recommended_for_phylogeny"] == "false"
    filtered_rows = _read_tsv_rows(filtered_path)
    assert filtered_rows == []


def test_report_prefers_checkm2_over_pyhmmer(tmp_path: Path):
    """When both sources report on a MAG, CheckM2 wins."""
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(
        manifest,
        [
            "GCA_BOTH\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\t-\t-\t-\t-\t-\t-",
        ],
    )
    pyhmmer_dir = tmp_path / "qc_assessor"
    _write_qc_tsv(
        pyhmmer_dir / "GCA_BOTH.qc.tsv",
        accession="GCA_BOTH",
        qc_source="pyhmmer",
        completeness=70.0,
        contamination=2.0,
        quality=QC_QUALITY_MEDIUM,
    )
    checkm2_dir = tmp_path / "checkm2_ingest"
    _write_qc_tsv(
        checkm2_dir / "GCA_BOTH.qc.tsv",
        accession="GCA_BOTH",
        qc_source="checkm2",
        completeness=95.0,
        contamination=1.0,
        quality=QC_QUALITY_HIGH,
    )
    per_mag = tmp_path / "farlip_hits_per_mag.tsv"
    _write_per_mag(per_mag, ["GCA_BOTH\t1\t10\t10\tfalse\tFULL_CLUSTER"])

    report_path = tmp_path / "qc_report.tsv"
    filtered_path = tmp_path / "qc_filtered_mags.tsv"
    write_qc_reports(
        manifest_path=manifest,
        pyhmmer_qc_dir=pyhmmer_dir,
        checkm2_ingest_dir=checkm2_dir,
        farlip_per_mag_path=per_mag,
        report_path=report_path,
        filtered_path=filtered_path,
    )
    rows = _read_tsv_rows(report_path)
    assert rows[0]["qc_source"] == "checkm2"
    assert rows[0]["completeness"] == "95.00"
    assert rows[0]["qc_quality"] == QC_QUALITY_HIGH


def test_report_handles_missing_qc(tmp_path: Path):
    """A MAG with no QC source still appears in the report (LOW quality)."""
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(
        manifest,
        [
            "GCA_NOQC\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\t-\t-\t-\t-\t-\t-",
        ],
    )
    per_mag = tmp_path / "farlip_hits_per_mag.tsv"
    _write_per_mag(per_mag, ["GCA_NOQC\t1\t10\t10\tfalse\tPARTIAL_CLUSTER"])
    report_path = tmp_path / "qc_report.tsv"
    filtered_path = tmp_path / "qc_filtered_mags.tsv"
    write_qc_reports(
        manifest_path=manifest,
        pyhmmer_qc_dir=tmp_path / "qc_assessor",
        checkm2_ingest_dir=tmp_path / "checkm2_ingest",
        farlip_per_mag_path=per_mag,
        report_path=report_path,
        filtered_path=filtered_path,
    )
    rows = _read_tsv_rows(report_path)
    assert len(rows) == 1
    assert rows[0]["qc_source"] == ""
    assert rows[0]["qc_quality"] == QC_QUALITY_LOW
    assert rows[0]["recommended_for_phylogeny"] == "false"


def test_report_excludes_non_assembly_sources(tmp_path: Path):
    """SILVA / GTDB / QC_MARKERS rows do not contribute report rows."""
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(
        manifest,
        [
            "GCA_OK\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\t-\t-\t-\t-\t-\t-",
            "SILVA_NR99\tSILVA\thttps://s\t/tmp/s\t-\t-\t-\tok\t-\tlbl\t-\t-\t-\t-\t-\t-",
            "BUSCO\tQC_MARKERS\thttps://m\t/tmp/m\t-\t-\t-\tok\t-\tlbl\t-\t-\t-\t-\t-\t-",
        ],
    )
    per_mag = tmp_path / "farlip_hits_per_mag.tsv"
    _write_per_mag(per_mag, ["GCA_OK\t0\t0\t0\tfalse\tNONE"])
    report_path = tmp_path / "qc_report.tsv"
    filtered_path = tmp_path / "qc_filtered_mags.tsv"
    write_qc_reports(
        manifest_path=manifest,
        pyhmmer_qc_dir=tmp_path / "qc_assessor",
        checkm2_ingest_dir=tmp_path / "checkm2_ingest",
        farlip_per_mag_path=per_mag,
        report_path=report_path,
        filtered_path=filtered_path,
    )
    rows = _read_tsv_rows(report_path)
    assert {r["accession"] for r in rows} == {"GCA_OK"}
