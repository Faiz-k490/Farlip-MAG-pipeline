"""Tests for ``src.obsidian_exporter``.

The exporter is a pure TSV-to-markdown transformation; the tests fabricate the
two synteny TSVs plus a tiny manifest and verify that:

* One per-MAG markdown file is written per candidate with parseable YAML
  frontmatter.
* The master index lists candidates grouped by classification.
* MAGs with zero candidate contigs are excluded from the per-MAG file count.
* YAML escaping handles values containing commas and colons.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest
import yaml

from src.obsidian_exporter import (
    CLASSIFICATION_ORDER,
    HITS_TABLE_HEADER,
    MASTER_INDEX_FILENAME,
    _format_frontmatter,
    _format_hits_table,
    _yaml_escape,
    export_vault,
)


_REPORT_HEADER = (
    "mag_accession\tcontig_id\tcontig_length\tn_markers\tn_total_hits_on_contig\t"
    "marker_names\tmarker_positions\tmin_evalue\tedge_truncated\tclassification\n"
)
_PER_MAG_HEADER = (
    "mag_accession\tn_candidate_contigs\tn_total_markers\tn_distinct_markers\t"
    "has_edge_truncated\ttop_classification\n"
)
_MANIFEST_HEADER = (
    "accession\tsource\turl\tlocal_path\tsha256\tbytes_downloaded\tfetched_at\t"
    "status\tnote\tlabel\tscientific_name\tsample_accession\tassembly_level\t"
    "format\tmd5\textracted_dir\n"
)


def _write_minimal_reports(
    tmp_path: Path,
    *,
    contig_rows: list[str],
    per_mag_rows: list[str],
    manifest_rows: list[str],
) -> tuple[Path, Path, Path, Path]:
    """Materialize the three input TSVs the exporter expects."""
    report = tmp_path / "report.tsv"
    per_mag = tmp_path / "per_mag.tsv"
    manifest = tmp_path / "manifest.tsv"
    vault = tmp_path / "vault"
    report.write_text(_REPORT_HEADER + "\n".join(contig_rows) + "\n", encoding="utf-8")
    per_mag.write_text(_PER_MAG_HEADER + "\n".join(per_mag_rows) + "\n", encoding="utf-8")
    manifest.write_text(_MANIFEST_HEADER + "\n".join(manifest_rows) + "\n", encoding="utf-8")
    return report, per_mag, manifest, vault


def test_constants_match_classification_order():
    """The display order must contain every Classification value."""
    assert "FULL_CLUSTER" in CLASSIFICATION_ORDER
    assert "EDGE_TRUNCATED" in CLASSIFICATION_ORDER
    assert "PARTIAL_CLUSTER" in CLASSIFICATION_ORDER
    assert MASTER_INDEX_FILENAME.endswith(".md")
    assert HITS_TABLE_HEADER[0] == "contig_id"


def test_yaml_escape_quotes_dangerous_chars():
    """Colons and commas in scalar values trigger single-quoting."""
    assert _yaml_escape("simple") == "simple"
    assert _yaml_escape("with:colon") == "'with:colon'"
    assert _yaml_escape("a, b") == "'a, b'"
    assert _yaml_escape("") == '""'


def test_format_frontmatter_yields_parseable_yaml():
    """The frontmatter block must round-trip through ``yaml.safe_load``."""
    block = _format_frontmatter(
        {
            "accession": "GCA_TEST",
            "n_farlip_markers": 21,
            "edge_truncated": False,
            "contig_ids": ["ctg1", "ctg2:weird"],
            "scientific_name": "Test sp. PCC9999",
        }
    )
    payload = block.strip("-\n")
    parsed = yaml.safe_load(payload)
    assert parsed["accession"] == "GCA_TEST"
    assert parsed["n_farlip_markers"] == 21
    assert parsed["edge_truncated"] is False
    assert parsed["contig_ids"] == ["ctg1", "ctg2:weird"]


def test_format_hits_table_renders_with_header():
    """The hits table starts with the canonical column header."""
    from src.obsidian_exporter import _ContigRow

    rows = [
        _ContigRow(
            mag_accession="GCA_TEST",
            contig_id="ctg1",
            contig_length=12345,
            n_markers=4,
            marker_names="rfpA,psbA2,apcE2,psaA2",
            marker_positions="1-100,200-300,400-500,600-700",
            min_evalue="1.00e-30",
            edge_truncated=False,
            classification="PARTIAL_CLUSTER",
        )
    ]
    table = _format_hits_table(rows)
    assert table.splitlines()[0].startswith("| contig_id |")
    assert "PARTIAL_CLUSTER" in table


def test_export_vault_skips_non_candidate_mags(tmp_path: Path):
    """MAGs with zero candidate contigs do not produce a markdown file."""
    contig_rows = [
        "GCA_HIT\tctg1\t12345\t5\t10\tA,B,C,D,E\t1-100,200-300,400-500,600-700,800-900"
        "\t1.00e-30\tfalse\tPARTIAL_CLUSTER",
    ]
    per_mag_rows = [
        "GCA_HIT\t1\t5\t5\tfalse\tPARTIAL_CLUSTER",
        "GCA_MISS\t0\t0\t0\tfalse\tNONE",
    ]
    manifest_rows = [
        "GCA_HIT\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\tTest sp. nov.\tSRS_X\t-\tfasta.gz\t-\t-",
    ]
    report, per_mag, manifest, vault = _write_minimal_reports(
        tmp_path,
        contig_rows=contig_rows,
        per_mag_rows=per_mag_rows,
        manifest_rows=manifest_rows,
    )
    n = export_vault(report, per_mag, manifest, vault)
    assert n == 1
    assert (vault / "GCA_HIT.md").exists()
    assert not (vault / "GCA_MISS.md").exists()
    assert (vault / MASTER_INDEX_FILENAME).exists()


def test_export_vault_frontmatter_is_parseable(tmp_path: Path):
    """The generated frontmatter must parse cleanly through ``yaml.safe_load``."""
    contig_rows = [
        "GCA_HIT\tctg1\t12345\t5\t10\tA,B,C,D,E\t1-100,200-300,400-500,600-700,800-900"
        "\t1.00e-30\tfalse\tPARTIAL_CLUSTER",
    ]
    per_mag_rows = ["GCA_HIT\t1\t5\t5\tfalse\tPARTIAL_CLUSTER"]
    manifest_rows = [
        "GCA_HIT\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\tCandidatus, Test sp.\tSRS_X\t-\tfasta.gz\t-\t-",
    ]
    report, per_mag, manifest, vault = _write_minimal_reports(
        tmp_path,
        contig_rows=contig_rows,
        per_mag_rows=per_mag_rows,
        manifest_rows=manifest_rows,
    )
    export_vault(report, per_mag, manifest, vault)
    md = (vault / "GCA_HIT.md").read_text(encoding="utf-8")
    front, _, _ = md.partition("\n---\n")
    payload = front.strip("-\n")
    parsed = yaml.safe_load(payload)
    assert parsed["accession"] == "GCA_HIT"
    assert parsed["scientific_name"] == "Candidatus, Test sp."
    assert parsed["classification"] == "PARTIAL_CLUSTER"
    assert parsed["contig_ids"] == ["ctg1"]


_QC_REPORT_HEADER = (
    "accession\tsource\tscientific_name\tqc_source\tmarker_set_label\t"
    "completeness\tcontamination\tn50_contig\tgenome_size_bp\t"
    "n_contigs\tlongest_contig_bp\tqc_quality\tn_farlip_distinct_markers\t"
    "farlip_top_classification\trecommended_for_phylogeny\n"
)


def _write_qc_report(tmp_path: Path, rows: list[str]) -> Path:
    """Materialize a Phase 3 ``qc_report.tsv`` fixture and return its path."""
    p = tmp_path / "qc_report.tsv"
    p.write_text(_QC_REPORT_HEADER + "\n".join(rows) + "\n", encoding="utf-8")
    return p


def test_export_vault_qc_aware_frontmatter(tmp_path: Path):
    """Passing ``qc_report_path`` adds the QC fields to the per-MAG frontmatter."""
    contig_rows = [
        "GCA_HIT\tctg1\t12345\t5\t10\tA,B,C,D,E\t1-100,200-300,400-500,600-700,800-900"
        "\t1.00e-30\tfalse\tPARTIAL_CLUSTER",
    ]
    per_mag_rows = ["GCA_HIT\t1\t5\t5\tfalse\tPARTIAL_CLUSTER"]
    manifest_rows = [
        "GCA_HIT\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\tTest sp.\tSAM\t-\t-\t-\t-",
    ]
    report, per_mag, manifest, vault = _write_minimal_reports(
        tmp_path,
        contig_rows=contig_rows,
        per_mag_rows=per_mag_rows,
        manifest_rows=manifest_rows,
    )
    qc_report = _write_qc_report(
        tmp_path,
        [
            "GCA_HIT\tENA\tTest sp.\tpyhmmer\tBUSCO_bacteria_odb10\t"
            "92.40\t1.20\t184050\t4923118\t10\t300000\tHIGH\t5\tPARTIAL_CLUSTER\ttrue",
        ],
    )
    n = export_vault(report, per_mag, manifest, vault, qc_report_path=qc_report)
    assert n == 1
    md = (vault / "GCA_HIT.md").read_text(encoding="utf-8")
    front, _, _ = md.partition("\n---\n")
    payload = front.strip("-\n")
    parsed = yaml.safe_load(payload)
    assert parsed["qc_quality"] == "HIGH"
    assert parsed["qc_source"] == "pyhmmer"
    assert parsed["completeness"] == 92.4
    assert parsed["contamination"] == 1.2
    assert parsed["n50_contig"] == 184050
    assert parsed["genome_size_bp"] == 4923118
    assert parsed["recommended_for_phylogeny"] is True
    # Body retains the FaRLiP table and gains a Quality assessment section.
    assert "## Quality assessment" in md


def test_export_vault_master_index_includes_quality_section(tmp_path: Path):
    """A QC-aware export adds a ``## Quality`` master-index section."""
    contig_rows = [
        "GCA_HIT\tctg1\t12345\t5\t10\tA,B,C,D,E\t1-100,200-300,400-500,600-700,800-900"
        "\t1.00e-30\tfalse\tPARTIAL_CLUSTER",
    ]
    per_mag_rows = ["GCA_HIT\t1\t5\t5\tfalse\tPARTIAL_CLUSTER"]
    manifest_rows = [
        "GCA_HIT\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\tTest\tSAM\t-\t-\t-\t-",
    ]
    report, per_mag, manifest, vault = _write_minimal_reports(
        tmp_path,
        contig_rows=contig_rows,
        per_mag_rows=per_mag_rows,
        manifest_rows=manifest_rows,
    )
    qc_report = _write_qc_report(
        tmp_path,
        [
            "GCA_HIT\tENA\tTest\tpyhmmer\tBUSCO\t92.40\t1.20\t100000\t3000000\t"
            "10\t100000\tHIGH\t5\tPARTIAL_CLUSTER\ttrue",
            "GCA_DROP\tENA\tFiller\tpyhmmer\tBUSCO\t30.00\t1.00\t10000\t1500000\t"
            "100\t10000\tLOW\t0\tNONE\tfalse",
        ],
    )
    export_vault(report, per_mag, manifest, vault, qc_report_path=qc_report)
    index = (vault / MASTER_INDEX_FILENAME).read_text(encoding="utf-8")
    assert "## Quality" in index
    assert "### HIGH" in index
    assert "[[GCA_HIT]]" in index
    # MAG with no FaRLiP candidate but a QC row still appears under Quality.
    assert "[[GCA_DROP]]" in index


def test_export_vault_omits_qc_when_path_none(tmp_path: Path):
    """Without a QC report, the Phase 2 layout is unchanged."""
    contig_rows = [
        "GCA_HIT\tctg1\t12345\t5\t10\tA,B,C,D,E\t1-100,200-300,400-500,600-700,800-900"
        "\t1.00e-30\tfalse\tPARTIAL_CLUSTER",
    ]
    per_mag_rows = ["GCA_HIT\t1\t5\t5\tfalse\tPARTIAL_CLUSTER"]
    manifest_rows = [
        "GCA_HIT\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\t-\t-\t-\t-\t-\t-",
    ]
    report, per_mag, manifest, vault = _write_minimal_reports(
        tmp_path,
        contig_rows=contig_rows,
        per_mag_rows=per_mag_rows,
        manifest_rows=manifest_rows,
    )
    export_vault(report, per_mag, manifest, vault)
    md = (vault / "GCA_HIT.md").read_text(encoding="utf-8")
    assert "qc_quality" not in md
    assert "## Quality assessment" not in md
    index = (vault / MASTER_INDEX_FILENAME).read_text(encoding="utf-8")
    assert "## Quality" not in index


def test_export_vault_master_index_groups_by_classification(tmp_path: Path):
    """The master index lists candidates under classification headings."""
    contig_rows = [
        "GCA_FULL\tctg1\t30000\t20\t30\t" + ",".join(f"H{i}" for i in range(20)) + "\t"
        + ",".join("1-100" for _ in range(20))
        + "\t1.00e-50\tfalse\tFULL_CLUSTER",
        "GCA_PART\tctg1\t10000\t4\t8\tH1,H2,H3,H4\t1-100,200-300,400-500,600-700"
        + "\t1.00e-20\tfalse\tPARTIAL_CLUSTER",
    ]
    per_mag_rows = [
        "GCA_FULL\t1\t20\t20\tfalse\tFULL_CLUSTER",
        "GCA_PART\t1\t4\t4\tfalse\tPARTIAL_CLUSTER",
    ]
    manifest_rows = [
        "GCA_FULL\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\t-\t-\t-\t-\t-\t-",
        "GCA_PART\tENA\thttps://y\t/tmp/y\t-\t-\t-\tok\t-\tlbl\t-\t-\t-\t-\t-\t-",
    ]
    report, per_mag, manifest, vault = _write_minimal_reports(
        tmp_path,
        contig_rows=contig_rows,
        per_mag_rows=per_mag_rows,
        manifest_rows=manifest_rows,
    )
    export_vault(report, per_mag, manifest, vault)
    index = (vault / MASTER_INDEX_FILENAME).read_text(encoding="utf-8")
    full_pos = index.find("## FULL_CLUSTER")
    part_pos = index.find("## PARTIAL_CLUSTER")
    assert 0 <= full_pos < part_pos
    assert "[[GCA_FULL]]" in index
    assert "[[GCA_PART]]" in index
