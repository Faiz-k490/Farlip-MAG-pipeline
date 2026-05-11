"""Tests for ``src/config_loader.py``.

These tests exercise the schema validation paths against the real, committed config
files plus several in-memory bad-shape fixtures. They are intentionally pure-Python
and require no network access.
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest

from src.config_loader import (
    ConfigError,
    DEFAULT_CONFIG_DIR,
    load_all,
    load_databases,
    load_extremophiles,
    load_references,
)


REPO_ROOT = Path(__file__).resolve().parent.parent


def test_load_references_real_file():
    cfg = load_references(REPO_ROOT / "config" / "references.json")
    assert cfg.primary.refseq_accession == "GCF_000733415.1"
    assert cfg.primary.label == "Leptolyngbya_sp_JSC-1"
    assert cfg.primary.source == "NCBI"
    assert cfg.secondary == []


def test_load_extremophiles_real_file():
    cfg = load_extremophiles(REPO_ROOT / "config" / "extremophiles.json")
    assert len(cfg.bioprojects) == 1
    assert cfg.bioprojects[0].accession == "PRJEB59431"
    assert cfg.bioprojects[0].source == "ENA"
    assert cfg.bioprojects[0].expected_assembly_count == 37
    # Template entry under cave_arctic_pending is silently dropped.
    assert cfg.cave_arctic_pending == []
    assert cfg.individual_accessions == []


def test_load_databases_real_file():
    cfg = load_databases(REPO_ROOT / "config" / "databases.yaml")
    assert cfg.silva_version == "138.2"
    assert cfg.gtdb_version == "R226"
    assert len(cfg.silva_assets) >= 1
    assert len(cfg.gtdb_assets) >= 2
    silva_urls = {a.url for a in cfg.silva_assets}
    assert any("SILVA_138.2_SSURef_NR99" in u for u in silva_urls)
    gtdb_urls = {a.url for a in cfg.gtdb_assets}
    assert any("bac120_taxonomy_r226" in u for u in gtdb_urls)


def test_load_all_real_files():
    bundle = load_all()
    assert bundle.references.primary.refseq_accession == "GCF_000733415.1"
    assert bundle.databases.silva_version == "138.2"


def test_references_missing_primary(tmp_path: Path):
    p = tmp_path / "references.json"
    p.write_text(json.dumps({"secondary": []}))
    with pytest.raises(ConfigError, match="Missing required field 'primary'"):
        load_references(p)


def test_references_missing_required_subfield(tmp_path: Path):
    p = tmp_path / "references.json"
    p.write_text(json.dumps({"primary": {"label": "x", "source": "NCBI"}}))
    with pytest.raises(ConfigError, match="refseq_accession"):
        load_references(p)


def test_references_bad_json(tmp_path: Path):
    p = tmp_path / "references.json"
    p.write_text("{ this is not json")
    with pytest.raises(ConfigError, match="Invalid JSON"):
        load_references(p)


def test_references_missing_file(tmp_path: Path):
    with pytest.raises(ConfigError, match="Config file not found"):
        load_references(tmp_path / "nope.json")


def test_extremophiles_empty_lists_ok(tmp_path: Path):
    p = tmp_path / "extremophiles.json"
    p.write_text(json.dumps({"bioprojects": [], "individual_accessions": [], "cave_arctic_pending": []}))
    cfg = load_extremophiles(p)
    assert cfg.bioprojects == []
    assert cfg.individual_accessions == []
    assert cfg.cave_arctic_pending == []


def test_extremophiles_pending_with_real_entry(tmp_path: Path):
    """When Dr. Barton drops a real cave accession in, it must parse cleanly."""
    p = tmp_path / "extremophiles.json"
    p.write_text(json.dumps({
        "bioprojects": [],
        "individual_accessions": [],
        "cave_arctic_pending": [
            {"accession": "GCA_111111111.1", "source": "NCBI", "label": "cave_x", "biome": "cave"}
        ],
    }))
    cfg = load_extremophiles(p)
    assert len(cfg.cave_arctic_pending) == 1
    assert cfg.cave_arctic_pending[0].biome == "cave"


def test_databases_missing_silva(tmp_path: Path):
    p = tmp_path / "databases.yaml"
    p.write_text(textwrap.dedent("""
        gtdb:
          version: "R226"
          bac120:
            label: "x"
            url: "https://example.com/x.tsv.gz"
            format: "tsv.gz"
    """))
    with pytest.raises(ConfigError, match="Missing required field 'silva'"):
        load_databases(p)


def test_databases_silva_with_no_assets(tmp_path: Path):
    p = tmp_path / "databases.yaml"
    p.write_text(textwrap.dedent("""
        silva:
          version: "138.2"
        gtdb:
          version: "R226"
          bac120:
            label: "x"
            url: "https://example.com/x.tsv.gz"
            format: "tsv.gz"
    """))
    with pytest.raises(ConfigError, match="No downloadable SILVA assets"):
        load_databases(p)


def test_default_config_dir_resolves_to_repo():
    assert (DEFAULT_CONFIG_DIR / "references.json").exists()
    assert (DEFAULT_CONFIG_DIR / "extremophiles.json").exists()
    assert (DEFAULT_CONFIG_DIR / "databases.yaml").exists()
