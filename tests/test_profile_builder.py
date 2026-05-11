"""Tests for ``src.profile_builder``.

Offline tests with in-memory mini-GFFs and protein FASTAs. They cover:

* GFF attribute parsing including URL-encoded values.
* Catalog YAML loading and missing-key error paths.
* Resolution priority (locus_tag > protein_id > fuzzy with cluster tiebreak).
* HMM emission (one HMM per catalog entry, sidecar TSV written).
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from src.profile_builder import (
    CATALOG_TSV_HEADER,
    JSC1_FARLIP_CLUSTER_CONTIG,
    JSC1_FARLIP_CLUSTER_END,
    JSC1_FARLIP_CLUSTER_START,
    CatalogEntry,
    GffCDS,
    ProfileBuildError,
    _hmm_name,
    _is_in_cluster,
    _parse_gff_attributes,
    _resolve_entry,
    build_profiles,
    load_catalog,
    parse_jsc1_gff,
    parse_jsc1_proteins,
)


_FAKE_PROTEIN_SEQ_A = "MAKAAVLLAAVAATALVAAAPALAQAEPAQAVQYPLKQGNTIRVAVLPF" + "A" * 60
_FAKE_PROTEIN_SEQ_B = "MGEELELELEELELE" + "L" * 80
_FAKE_PROTEIN_SEQ_C = "MTRADVLAIQEEVRARLAEAGYDPAALKLDAS" + "T" * 60


def _write_demo_gff(path: Path) -> None:
    """Write a tiny GFF with three CDS: two inside the cluster window, one outside."""
    body = textwrap.dedent(
        f"""\
        ##gff-version 3
        {JSC1_FARLIP_CLUSTER_CONTIG}\tRefSeq\tCDS\t1095910\t1096998\t.\t+\t0\tID=cds1;locus_tag=LOC_A;gene=psbA;product=photosystem II q(b) protein;protein_id=PROT_A
        {JSC1_FARLIP_CLUSTER_CONTIG}\tRefSeq\tCDS\t1097106\t1097582\t.\t+\t0\tID=cds2;locus_tag=LOC_B;product=allophycocyanin;protein_id=PROT_B
        OTHER_CONTIG\tRefSeq\tCDS\t100\t1099\t.\t+\t0\tID=cds3;locus_tag=LOC_C;gene=psbA;product=photosystem II q(b) protein;protein_id=PROT_C
        {JSC1_FARLIP_CLUSTER_CONTIG}\tRefSeq\tCDS\t1098000\t1098500\t.\t+\t0\tID=cds_pseudo;locus_tag=LOC_PSEUDO;product=pseudogene;pseudo=true
        """
    )
    path.write_text(body, encoding="utf-8")


def _write_demo_faa(path: Path) -> None:
    path.write_text(
        f">PROT_A description\n{_FAKE_PROTEIN_SEQ_A}\n"
        f">PROT_B\n{_FAKE_PROTEIN_SEQ_B}\n"
        f">PROT_C\n{_FAKE_PROTEIN_SEQ_C}\n",
        encoding="utf-8",
    )


def _write_demo_catalog(path: Path) -> None:
    path.write_text(
        textwrap.dedent(
            """\
            schema_version: "1.0"
            source: "test"
            genes:
              - name: psbA2
                role: D1 paralog
                match:
                  locus_tag: LOC_A
              - name: apcD2
                role: allophycocyanin paralog
                match:
                  protein_id: PROT_B
            """
        ),
        encoding="utf-8",
    )


def test_constants_define_real_cluster_window():
    """The constants must describe a valid, non-empty cluster window."""
    assert JSC1_FARLIP_CLUSTER_CONTIG.startswith("NZ_")
    assert JSC1_FARLIP_CLUSTER_END > JSC1_FARLIP_CLUSTER_START > 0
    assert JSC1_FARLIP_CLUSTER_END - JSC1_FARLIP_CLUSTER_START < 50_000
    assert CATALOG_TSV_HEADER[0] == "hmm_name"


def test_parse_gff_attributes_decodes_url_escapes():
    """URL-escaped commas / semicolons in GFF values must round-trip."""
    attrs = _parse_gff_attributes("ID=cds1;product=foo%2C bar%3B baz;locus_tag=L1")
    assert attrs["product"] == "foo, bar; baz"
    assert attrs["locus_tag"] == "L1"


def test_parse_jsc1_gff_skips_pseudogenes(tmp_path: Path):
    """Rows with ``pseudo=true`` are not returned (no matching protein.faa entry)."""
    gff_path = tmp_path / "demo.gff"
    _write_demo_gff(gff_path)
    rows = parse_jsc1_gff(gff_path)
    locus_tags = [r.locus_tag for r in rows]
    assert "LOC_A" in locus_tags
    assert "LOC_B" in locus_tags
    assert "LOC_C" in locus_tags
    assert "LOC_PSEUDO" not in locus_tags


def test_parse_jsc1_proteins_yields_first_token_id(tmp_path: Path):
    """Protein FASTA parsing keeps only the first whitespace-delimited header token."""
    faa = tmp_path / "demo.faa"
    _write_demo_faa(faa)
    proteins = parse_jsc1_proteins(faa)
    assert set(proteins.keys()) == {"PROT_A", "PROT_B", "PROT_C"}
    assert proteins["PROT_A"] == _FAKE_PROTEIN_SEQ_A


def test_is_in_cluster_window():
    """The cluster-window predicate is correct for both contig and coordinates."""
    inside = GffCDS(
        contig=JSC1_FARLIP_CLUSTER_CONTIG,
        start=JSC1_FARLIP_CLUSTER_START + 100,
        end=JSC1_FARLIP_CLUSTER_START + 200,
        strand="+",
    )
    outside_contig = GffCDS(
        contig="OTHER",
        start=JSC1_FARLIP_CLUSTER_START + 100,
        end=JSC1_FARLIP_CLUSTER_START + 200,
        strand="+",
    )
    outside_coords = GffCDS(
        contig=JSC1_FARLIP_CLUSTER_CONTIG,
        start=10,
        end=200,
        strand="+",
    )
    assert _is_in_cluster(inside) is True
    assert _is_in_cluster(outside_contig) is False
    assert _is_in_cluster(outside_coords) is False


def test_resolve_entry_prefers_locus_tag(tmp_path: Path):
    """When ``locus_tag`` matches, it wins regardless of other hints."""
    gff_path = tmp_path / "demo.gff"
    faa_path = tmp_path / "demo.faa"
    _write_demo_gff(gff_path)
    _write_demo_faa(faa_path)
    cds = parse_jsc1_gff(gff_path)
    proteins = parse_jsc1_proteins(faa_path)
    entry = CatalogEntry(
        name="psbA2",
        role="",
        locus_tag="LOC_A",
        protein_id=None,
        gene_tag=None,
        product_substring=None,
    )
    resolved = _resolve_entry(entry, cds, proteins)
    assert resolved.protein_id == "PROT_A"
    assert resolved.resolution == "locus_tag"
    assert resolved.in_cluster is True


def test_resolve_entry_falls_back_to_protein_id(tmp_path: Path):
    """When ``locus_tag`` does not match, ``protein_id`` is used next."""
    gff_path = tmp_path / "demo.gff"
    faa_path = tmp_path / "demo.faa"
    _write_demo_gff(gff_path)
    _write_demo_faa(faa_path)
    cds = parse_jsc1_gff(gff_path)
    proteins = parse_jsc1_proteins(faa_path)
    entry = CatalogEntry(
        name="apcD2",
        role="",
        locus_tag="DOES_NOT_EXIST",
        protein_id="PROT_B",
        gene_tag=None,
        product_substring=None,
    )
    resolved = _resolve_entry(entry, cds, proteins)
    assert resolved.protein_id == "PROT_B"
    assert resolved.resolution == "protein_id"


def test_resolve_entry_fuzzy_with_cluster_tiebreak(tmp_path: Path):
    """When multiple candidates qualify, the one inside the cluster wins."""
    gff_path = tmp_path / "demo.gff"
    faa_path = tmp_path / "demo.faa"
    _write_demo_gff(gff_path)
    _write_demo_faa(faa_path)
    cds = parse_jsc1_gff(gff_path)
    proteins = parse_jsc1_proteins(faa_path)
    entry = CatalogEntry(
        name="psbA2_fuzzy",
        role="",
        locus_tag=None,
        protein_id=None,
        gene_tag="psbA",
        product_substring="photosystem II q(b)",
    )
    resolved = _resolve_entry(entry, cds, proteins)
    assert resolved.protein_id == "PROT_A"
    assert resolved.in_cluster is True


def test_resolve_entry_raises_when_unresolvable(tmp_path: Path):
    """A catalog entry with no matching CDS triggers ProfileBuildError."""
    gff_path = tmp_path / "demo.gff"
    faa_path = tmp_path / "demo.faa"
    _write_demo_gff(gff_path)
    _write_demo_faa(faa_path)
    cds = parse_jsc1_gff(gff_path)
    proteins = parse_jsc1_proteins(faa_path)
    entry = CatalogEntry(
        name="ghost",
        role="",
        locus_tag="NO_SUCH_LOCUS",
        protein_id="NO_SUCH_PROTEIN",
        gene_tag=None,
        product_substring=None,
    )
    with pytest.raises(ProfileBuildError):
        _resolve_entry(entry, cds, proteins)


def test_load_catalog_missing_genes_key_raises(tmp_path: Path):
    """A YAML without a top-level ``genes`` key must fail loudly."""
    bad = tmp_path / "bad.yaml"
    bad.write_text("schema_version: '1.0'\n", encoding="utf-8")
    with pytest.raises(ProfileBuildError):
        load_catalog(bad)


def test_hmm_name_sanitizer():
    """The HMM name sanitizer must collapse disallowed characters."""
    assert _hmm_name("rfpA") == "rfpA"
    assert _hmm_name("psbA2") == "psbA2"
    assert _hmm_name("cluster:special name") == "cluster_special_name"


def test_build_profiles_writes_hmm_and_sidecar(tmp_path: Path):
    """End-to-end with toy inputs writes the HMM file and the catalog TSV."""
    gff_path = tmp_path / "demo.gff"
    faa_path = tmp_path / "demo.faa"
    catalog_path = tmp_path / "catalog.yaml"
    _write_demo_gff(gff_path)
    _write_demo_faa(faa_path)
    _write_demo_catalog(catalog_path)
    output_dir = tmp_path / "hmm_profiles"
    resolved = build_profiles(catalog_path, gff_path, faa_path, output_dir)
    assert len(resolved) == 2
    assert (output_dir / "farlip_21.hmm").exists()
    sidecar = output_dir / "farlip_21.catalog.tsv"
    assert sidecar.exists()
    header = sidecar.read_text(encoding="utf-8").splitlines()[0]
    assert header == "\t".join(CATALOG_TSV_HEADER)
