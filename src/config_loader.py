"""Configuration loader for the FaRLiP ingestion pipeline.

Reads the three configuration files that drive every fetcher:

* ``config/references.json``    - FaRLiP reference anchors (JSC-1 + secondary).
* ``config/extremophiles.json`` - MAG BioProjects, individual accessions, and the
  ``cave_arctic_pending`` placeholder list.
* ``config/databases.yaml``     - SILVA + GTDB release versions and download URLs.

Configs are validated against minimal schemas in this module - we do not pull in
a heavyweight dependency like ``jsonschema`` because the schemas are small and the
error messages are clearer when raised by hand.

No environment variables are read here. Network-facing modules consume the loaded
config dataclasses directly.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG_DIR = REPO_ROOT / "config"


class ConfigError(ValueError):
    """Raised when a configuration file is missing required fields or has a bad shape."""


@dataclass(frozen=True)
class ReferenceEntry:
    """A single FaRLiP reference organism."""

    label: str
    refseq_accession: str
    source: str
    bioproject: str | None = None
    genbank_accession: str | None = None
    taxonomy_id: int | None = None
    rationale: str | None = None


@dataclass(frozen=True)
class ReferencesConfig:
    """Parsed contents of ``config/references.json``."""

    primary: ReferenceEntry
    secondary: list[ReferenceEntry] = field(default_factory=list)

    def all_entries(self) -> list[ReferenceEntry]:
        """Return primary followed by every secondary reference in declaration order."""
        return [self.primary, *self.secondary]


@dataclass(frozen=True)
class BioProjectEntry:
    """A BioProject whose every member assembly should be fetched."""

    accession: str
    source: str
    label: str
    notes: str | None = None
    expected_assembly_count: int | None = None


@dataclass(frozen=True)
class IndividualAccession:
    """A single one-off assembly accession (not bound to a BioProject)."""

    accession: str
    source: str
    label: str
    biome: str | None = None
    notes: str | None = None


@dataclass(frozen=True)
class ExtremophilesConfig:
    """Parsed contents of ``config/extremophiles.json``."""

    bioprojects: list[BioProjectEntry]
    individual_accessions: list[IndividualAccession]
    cave_arctic_pending: list[IndividualAccession]


@dataclass(frozen=True)
class DatabaseAsset:
    """A single downloadable database file (SILVA SSU, GTDB taxonomy, etc.)."""

    label: str
    url: str
    fmt: str
    md5_url: str | None = None
    notes: str | None = None


@dataclass(frozen=True)
class DatabasesConfig:
    """Parsed contents of ``config/databases.yaml``."""

    silva_version: str
    silva_assets: list[DatabaseAsset]
    gtdb_version: str
    gtdb_assets: list[DatabaseAsset]

    def all_assets(self) -> list[DatabaseAsset]:
        """Return SILVA assets followed by GTDB assets."""
        return [*self.silva_assets, *self.gtdb_assets]


def _read_json(path: Path) -> dict[str, Any]:
    """Read and decode a JSON file with a clear error if it is missing or malformed."""
    if not path.exists():
        raise ConfigError(f"Config file not found: {path}")
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ConfigError(f"Invalid JSON in {path}: {exc}") from exc


def _read_yaml(path: Path) -> dict[str, Any]:
    """Read and decode a YAML file with a clear error if it is missing or malformed."""
    if not path.exists():
        raise ConfigError(f"Config file not found: {path}")
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigError(f"Invalid YAML in {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigError(f"Top-level YAML in {path} must be a mapping, got {type(data).__name__}")
    return data


def _require(d: dict[str, Any], key: str, where: str) -> Any:
    """Fetch ``d[key]`` or raise ConfigError naming the missing field and its location."""
    if key not in d:
        raise ConfigError(f"Missing required field '{key}' in {where}")
    return d[key]


def _parse_reference_entry(raw: dict[str, Any], where: str) -> ReferenceEntry:
    """Validate and parse one reference organism entry."""
    return ReferenceEntry(
        label=_require(raw, "label", where),
        refseq_accession=_require(raw, "refseq_accession", where),
        source=_require(raw, "source", where),
        bioproject=raw.get("bioproject"),
        genbank_accession=raw.get("genbank_accession"),
        taxonomy_id=raw.get("taxonomy_id"),
        rationale=raw.get("rationale"),
    )


def load_references(path: Path | None = None) -> ReferencesConfig:
    """Load and validate ``config/references.json``.

    Args:
        path: Optional explicit path. Defaults to ``<repo>/config/references.json``.

    Returns:
        A ``ReferencesConfig`` with at least one primary anchor and zero or more secondary.

    Raises:
        ConfigError: If the file is missing, malformed, or omits required fields.
    """
    path = path or (DEFAULT_CONFIG_DIR / "references.json")
    raw = _read_json(path)
    primary_raw = _require(raw, "primary", str(path))
    if not isinstance(primary_raw, dict):
        raise ConfigError(f"'primary' must be a mapping in {path}")
    primary = _parse_reference_entry(primary_raw, f"{path}::primary")

    secondary_raw = raw.get("secondary", [])
    if not isinstance(secondary_raw, list):
        raise ConfigError(f"'secondary' must be a list in {path}")
    secondary = [
        _parse_reference_entry(item, f"{path}::secondary[{i}]")
        for i, item in enumerate(secondary_raw)
    ]
    return ReferencesConfig(primary=primary, secondary=secondary)


def _parse_bioproject(raw: dict[str, Any], where: str) -> BioProjectEntry:
    return BioProjectEntry(
        accession=_require(raw, "accession", where),
        source=_require(raw, "source", where),
        label=_require(raw, "label", where),
        notes=raw.get("notes"),
        expected_assembly_count=raw.get("expected_assembly_count"),
    )


def _parse_individual(raw: dict[str, Any], where: str) -> IndividualAccession | None:
    """Parse an individual-accession entry, skipping placeholder/template entries.

    Template entries (used in extremophiles.json to document the expected shape) carry
    ``$comment`` or ``$example`` keys instead of real fields. We silently drop them so the
    config file can self-document without polluting the parsed result.
    """
    if "$comment" in raw or "$example" in raw:
        return None
    return IndividualAccession(
        accession=_require(raw, "accession", where),
        source=_require(raw, "source", where),
        label=_require(raw, "label", where),
        biome=raw.get("biome"),
        notes=raw.get("notes"),
    )


def load_extremophiles(path: Path | None = None) -> ExtremophilesConfig:
    """Load and validate ``config/extremophiles.json``.

    Args:
        path: Optional explicit path. Defaults to ``<repo>/config/extremophiles.json``.

    Returns:
        An ``ExtremophilesConfig``. Empty lists are allowed everywhere - the pipeline
        is designed to run with zero MAGs (skips that stage) so the config file can be
        committed with placeholders until accessions arrive.

    Raises:
        ConfigError: If the file is missing, malformed, or omits required fields.
    """
    path = path or (DEFAULT_CONFIG_DIR / "extremophiles.json")
    raw = _read_json(path)

    bp_raw = raw.get("bioprojects", [])
    ind_raw = raw.get("individual_accessions", [])
    pending_raw = raw.get("cave_arctic_pending", [])

    for key, value in [("bioprojects", bp_raw), ("individual_accessions", ind_raw), ("cave_arctic_pending", pending_raw)]:
        if not isinstance(value, list):
            raise ConfigError(f"'{key}' must be a list in {path}")

    bioprojects = [_parse_bioproject(item, f"{path}::bioprojects[{i}]") for i, item in enumerate(bp_raw)]
    individual = [
        parsed
        for i, item in enumerate(ind_raw)
        if (parsed := _parse_individual(item, f"{path}::individual_accessions[{i}]")) is not None
    ]
    pending = [
        parsed
        for i, item in enumerate(pending_raw)
        if (parsed := _parse_individual(item, f"{path}::cave_arctic_pending[{i}]")) is not None
    ]
    return ExtremophilesConfig(
        bioprojects=bioprojects,
        individual_accessions=individual,
        cave_arctic_pending=pending,
    )


def _parse_asset(raw: dict[str, Any], where: str) -> DatabaseAsset:
    return DatabaseAsset(
        label=_require(raw, "label", where),
        url=_require(raw, "url", where),
        fmt=_require(raw, "format", where),
        md5_url=raw.get("md5_url"),
        notes=raw.get("notes"),
    )


def load_databases(path: Path | None = None) -> DatabasesConfig:
    """Load and validate ``config/databases.yaml``.

    Args:
        path: Optional explicit path. Defaults to ``<repo>/config/databases.yaml``.

    Returns:
        A ``DatabasesConfig`` with SILVA and GTDB version strings and one or more
        downloadable assets per database.

    Raises:
        ConfigError: If the file is missing, malformed, or omits required sections.
    """
    path = path or (DEFAULT_CONFIG_DIR / "databases.yaml")
    raw = _read_yaml(path)

    silva = _require(raw, "silva", str(path))
    gtdb = _require(raw, "gtdb", str(path))

    # Every key under silva/gtdb that maps to a dict with a 'url' field is treated as
    # a downloadable asset. This lets us add new SILVA exports (LSU, taxonomy, etc.)
    # without touching loader code.
    silva_assets = [
        _parse_asset(v, f"{path}::silva.{k}")
        for k, v in silva.items()
        if isinstance(v, dict) and "url" in v
    ]
    gtdb_assets = [
        _parse_asset(v, f"{path}::gtdb.{k}")
        for k, v in gtdb.items()
        if isinstance(v, dict) and "url" in v
    ]

    if not silva_assets:
        raise ConfigError(f"No downloadable SILVA assets found in {path}")
    if not gtdb_assets:
        raise ConfigError(f"No downloadable GTDB assets found in {path}")

    return DatabasesConfig(
        silva_version=str(_require(silva, "version", f"{path}::silva")),
        silva_assets=silva_assets,
        gtdb_version=str(_require(gtdb, "version", f"{path}::gtdb")),
        gtdb_assets=gtdb_assets,
    )


@dataclass(frozen=True)
class PipelineConfig:
    """Aggregate config bundle used by ``src/pipeline.py``."""

    references: ReferencesConfig
    extremophiles: ExtremophilesConfig
    databases: DatabasesConfig


def load_all(config_dir: Path | None = None) -> PipelineConfig:
    """Load every config file and return a single bundle.

    Args:
        config_dir: Optional directory containing the three config files. Defaults to
            ``<repo>/config``.

    Returns:
        A ``PipelineConfig`` aggregating references, extremophiles, and databases.

    Raises:
        ConfigError: If any of the three config files fails to load or validate.
    """
    cfg_dir = config_dir or DEFAULT_CONFIG_DIR
    return PipelineConfig(
        references=load_references(cfg_dir / "references.json"),
        extremophiles=load_extremophiles(cfg_dir / "extremophiles.json"),
        databases=load_databases(cfg_dir / "databases.yaml"),
    )
