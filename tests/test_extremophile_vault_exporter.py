"""Tests for ``src.extremophile_vault_exporter``.

Pure transformation tests - the exporter is a TSV -> markdown step, so we
hand it a fake ``qc_filtered_mags.tsv`` and (optionally) a fake
``routing_log.tsv`` and assert on the rendered output.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from src.extremophile_vault_exporter import (
    FRONTMATTER_KEYS,
    MASTER_INDEX_FILENAME,
    PLACEMENT_PENDING_GTDBTK,
    PLACEMENT_PENDING_SILVA,
    PLACEMENT_UNKNOWN,
    export_extremophile_vault,
    list_filtered_accessions,
)
from src.ssu_bundler import RoutingDecision, write_routing_log
from src.ssu_extractor import SSU_DECISION_GTDBTK, SSU_DECISION_SILVA


_FILTERED_HEADER = (
    "accession\tsource\tscientific_name\tqc_source\tmarker_set_label\t"
    "completeness\tcontamination\tn50_contig\tgenome_size_bp\t"
    "n_contigs\tlongest_contig_bp\tqc_quality\tn_farlip_distinct_markers\t"
    "farlip_top_classification\trecommended_for_phylogeny\n"
)


def _write_filtered(path: Path, rows: list[str]) -> None:
    """Write a Phase 3 ``qc_filtered_mags.tsv`` fixture."""
    path.write_text(_FILTERED_HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def _decision(accession: str, decision: str, fasta: Path) -> RoutingDecision:
    """Build a fake routing decision."""
    return RoutingDecision(
        accession=accession,
        source="ENA",
        fasta_path=fasta,
        n_contigs_scanned=10,
        n_16S_predicted=1 if decision == SSU_DECISION_SILVA else 0,
        decision=decision,
        output_path=fasta,
        barrnap_version="barrnap 0.9",
        notes="",
    )


def _parse_frontmatter(md: str) -> dict:
    """Pull the YAML frontmatter out of a markdown file."""
    front, _, _ = md.partition("\n---\n")
    payload = front.strip("-\n")
    return yaml.safe_load(payload)


def test_constants_match_spec():
    """The five frontmatter keys land in the exact mandated order/casing."""
    assert FRONTMATTER_KEYS == (
        "Accession",
        "Completeness",
        "Contamination",
        "Taxonomic_Placement",
        "FaRLiP_Cluster_Present",
    )
    assert MASTER_INDEX_FILENAME == "High_Quality_Extremophile_Genomes.md"


def test_export_writes_per_mag_files_with_exact_frontmatter(tmp_path: Path):
    """Each per-MAG file carries exactly the five expected frontmatter keys."""
    filtered = tmp_path / "qc_filtered_mags.tsv"
    _write_filtered(
        filtered,
        [
            "GCA_AAA\tENA\tCandidatus, sp.\tpyhmmer\tBUSCO\t92.40\t1.20\t"
            "184050\t4923118\t10\t300000\tHIGH\t5\tPARTIAL_CLUSTER\ttrue",
        ],
    )

    vault = tmp_path / "extremophile_vault"
    n = export_extremophile_vault(filtered, vault)
    assert n == 1

    md = (vault / "GCA_AAA.md").read_text(encoding="utf-8")
    parsed = _parse_frontmatter(md)
    # The five keys must be present and in the canonical casing.
    assert set(parsed.keys()) == set(FRONTMATTER_KEYS)
    assert parsed["Accession"] == "GCA_AAA"
    assert parsed["Completeness"] == 92.40
    assert parsed["Contamination"] == 1.20
    # FaRLiP_Cluster_Present is intentionally blank.
    assert parsed["FaRLiP_Cluster_Present"] is None
    # No routing log => Taxonomic_Placement is empty.
    assert parsed["Taxonomic_Placement"] is None or parsed["Taxonomic_Placement"] == PLACEMENT_UNKNOWN


def test_master_index_lists_every_filtered_mag(tmp_path: Path):
    """The master index has a wikilink for every per-MAG note."""
    filtered = tmp_path / "qc_filtered_mags.tsv"
    _write_filtered(
        filtered,
        [
            "GCA_AAA\tENA\tFoo\tpyhmmer\tBUSCO\t92.40\t1.20\t100000\t3000000\t"
            "10\t100000\tHIGH\t5\tPARTIAL_CLUSTER\ttrue",
            "GCA_BBB\tENA\tBar\tpyhmmer\tBUSCO\t72.10\t3.30\t50000\t2500000\t"
            "20\t80000\tMEDIUM\t3\tEDGE_TRUNCATED\ttrue",
        ],
    )
    vault = tmp_path / "extremophile_vault"
    export_extremophile_vault(filtered, vault)
    index = (vault / MASTER_INDEX_FILENAME).read_text(encoding="utf-8")
    assert "[[GCA_AAA]]" in index
    assert "[[GCA_BBB]]" in index
    assert "## HIGH" in index
    assert "## MEDIUM" in index


def test_taxonomic_placement_uses_routing_decision(tmp_path: Path):
    """When a routing log is supplied, the placement field tracks the decision."""
    filtered = tmp_path / "qc_filtered_mags.tsv"
    _write_filtered(
        filtered,
        [
            "GCA_AAA\tENA\tFoo\tpyhmmer\tBUSCO\t92.40\t1.20\t100000\t3000000\t"
            "10\t100000\tHIGH\t5\tPARTIAL_CLUSTER\ttrue",
            "GCA_BBB\tENA\tBar\tpyhmmer\tBUSCO\t72.10\t3.30\t50000\t2500000\t"
            "20\t80000\tMEDIUM\t3\tEDGE_TRUNCATED\ttrue",
        ],
    )
    log = tmp_path / "routing_log.tsv"
    write_routing_log(
        [
            _decision("GCA_AAA", SSU_DECISION_SILVA, tmp_path / "GCA_AAA.fna"),
            _decision("GCA_BBB", SSU_DECISION_GTDBTK, tmp_path / "GCA_BBB.fna"),
        ],
        log,
    )

    vault = tmp_path / "extremophile_vault"
    export_extremophile_vault(filtered, vault, routing_log_path=log)

    md_a = _parse_frontmatter((vault / "GCA_AAA.md").read_text(encoding="utf-8"))
    md_b = _parse_frontmatter((vault / "GCA_BBB.md").read_text(encoding="utf-8"))
    assert md_a["Taxonomic_Placement"] == PLACEMENT_PENDING_SILVA
    assert md_b["Taxonomic_Placement"] == PLACEMENT_PENDING_GTDBTK


def test_export_handles_empty_filtered_tsv(tmp_path: Path):
    """An empty filter still creates the master index (with zero entries)."""
    filtered = tmp_path / "qc_filtered_mags.tsv"
    _write_filtered(filtered, [])
    vault = tmp_path / "extremophile_vault"
    n = export_extremophile_vault(filtered, vault)
    assert n == 0
    index = vault / MASTER_INDEX_FILENAME
    assert index.exists()
    body = index.read_text(encoding="utf-8")
    assert "Total candidates:** 0" in body


def test_list_filtered_accessions_matches_input(tmp_path: Path):
    """The convenience helper enumerates every filtered MAG."""
    filtered = tmp_path / "qc_filtered_mags.tsv"
    _write_filtered(
        filtered,
        [
            "GCA_BBB\tENA\tFoo\tpyhmmer\tBUSCO\t90\t2\t1\t1\t1\t1\tHIGH\t1\tFOO\ttrue",
            "GCA_AAA\tENA\tBar\tpyhmmer\tBUSCO\t91\t3\t1\t1\t1\t1\tHIGH\t1\tFOO\ttrue",
        ],
    )
    assert list_filtered_accessions(filtered) == ["GCA_AAA", "GCA_BBB"]


def test_pipe_in_scientific_name_does_not_break_yaml(tmp_path: Path):
    """A scientific name with awkward characters still yields parseable YAML."""
    filtered = tmp_path / "qc_filtered_mags.tsv"
    _write_filtered(
        filtered,
        [
            "GCA_AAA\tENA\tFunky: name, with chars\tpyhmmer\tBUSCO\t90\t2\t1\t1\t1\t1\tHIGH\t1\tFOO\ttrue",
        ],
    )
    vault = tmp_path / "extremophile_vault"
    export_extremophile_vault(filtered, vault)
    parsed = _parse_frontmatter((vault / "GCA_AAA.md").read_text(encoding="utf-8"))
    assert parsed["Accession"] == "GCA_AAA"
