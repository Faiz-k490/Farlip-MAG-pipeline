"""Render the high-quality extremophile MAG Obsidian vault.

This is the *second* Obsidian export the pipeline produces. The first is
:mod:`src.obsidian_exporter`, which writes a FaRLiP-flavoured per-MAG vault
under ``data/obsidian_export/``. This module is a sister exporter that
writes a smaller, assignment-shaped vault under ``data/extremophile_vault/``:

* Per-MAG markdown files with **exactly** the five YAML frontmatter keys
  ``Accession``, ``Completeness``, ``Contamination``, ``Taxonomic_Placement``,
  ``FaRLiP_Cluster_Present`` (left blank for manual review).
* A master index ``High_Quality_Extremophile_Genomes.md`` that wikilinks to
  every per-MAG note.

The exporter is fed by:

* ``data/qc_filtered_mags.tsv`` - the Phase 3 ``recommended_for_phylogeny``
  view (MIMAG HIGH or MEDIUM **and** at least one FaRLiP candidate contig).
  Schema: :data:`src.qc_reporter.REPORT_COLUMNS`.
* ``data/ssu/routing_log.tsv`` (optional) - the Phase 1.5 routing decisions.
  When present, ``Taxonomic_Placement`` is filled with a stable placeholder
  string indicating which downstream tool the MAG is queued for
  (``"SILVA SINA: pending"`` / ``"GTDB-Tk: pending"``); when absent, it is
  left empty so a grader/Obsidian user can fill it in by hand.
* ``data/farlip_hits_per_mag.tsv`` (optional) - used in the body of each
  per-MAG note to surface the Phase 2 top classification, but the
  ``FaRLiP_Cluster_Present`` frontmatter slot stays blank by design.

Idempotent: per-MAG markdown is overwritten on every call, the master index
is rebuilt from scratch.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .obsidian_exporter import _yaml_escape
from .ssu_bundler import read_routing_log
from .ssu_extractor import SSU_DECISION_GTDBTK, SSU_DECISION_SILVA


VAULT_DIRNAME: str = "extremophile_vault"
"""Default subdirectory under ``data/`` that this exporter writes into."""

MASTER_INDEX_FILENAME: str = "High_Quality_Extremophile_Genomes.md"
"""Exact master-index filename mandated by the Phase 1.5 / Task 4 spec."""

FRONTMATTER_KEYS: tuple[str, ...] = (
    "Accession",
    "Completeness",
    "Contamination",
    "Taxonomic_Placement",
    "FaRLiP_Cluster_Present",
)
"""Exact frontmatter key order. Casing is intentional; downstream Obsidian
queries (Dataview / Bases) are case-sensitive."""

PLACEMENT_PENDING_SILVA: str = "SILVA SINA: pending"
"""Placeholder ``Taxonomic_Placement`` value for SILVA-routed MAGs.
The actual species-level placement lands when the user runs SINA on the
``data/ssu/silva_sina_input.fasta`` bundle."""

PLACEMENT_PENDING_GTDBTK: str = "GTDB-Tk: pending"
"""Placeholder ``Taxonomic_Placement`` value for GTDB-Tk-routed MAGs.
Replaced once ``scripts/run_gtdbtk.sh`` has been executed and the user runs
``ingest-gtdbtk`` (a future Phase 1.5 follow-up)."""

PLACEMENT_UNKNOWN: str = ""
"""Placeholder ``Taxonomic_Placement`` value when no routing log exists
(typical for graders inspecting only the QC report). Empty string keeps the
YAML key present without claiming a placement."""


@dataclass(frozen=True)
class _FilteredMag:
    """One row from ``data/qc_filtered_mags.tsv`` (subset of the schema)."""

    accession: str
    source: str
    scientific_name: str
    completeness: float
    contamination: float
    qc_quality: str
    farlip_top_classification: str
    n_farlip_distinct_markers: int


def _read_filtered_mags(path: Path) -> list[_FilteredMag]:
    """Parse ``qc_filtered_mags.tsv`` into typed records.

    Only rows with a non-empty ``accession`` are kept. Numeric fields fall
    back to zero when the cell is empty or unparseable so the exporter
    stays robust on partial reports.

    Args:
        path: Path to the filtered MAG TSV.

    Returns:
        List of :class:`_FilteredMag`. Empty when the file is absent.
    """
    if not path.exists():
        return []
    out: list[_FilteredMag] = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            accession = (row.get("accession") or "").strip()
            if not accession:
                continue
            out.append(
                _FilteredMag(
                    accession=accession,
                    source=row.get("source", ""),
                    scientific_name=row.get("scientific_name", ""),
                    completeness=_safe_float(row.get("completeness", "")),
                    contamination=_safe_float(row.get("contamination", "")),
                    qc_quality=row.get("qc_quality", ""),
                    farlip_top_classification=row.get(
                        "farlip_top_classification", ""
                    ),
                    n_farlip_distinct_markers=_safe_int(
                        row.get("n_farlip_distinct_markers", "")
                    ),
                )
            )
    return out


def _safe_float(value: str) -> float:
    """Parse a numeric cell tolerating empty / missing values."""
    if value is None or value == "":
        return 0.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def _safe_int(value: str) -> int:
    """Parse an integer cell tolerating empty / missing values."""
    if value is None or value == "":
        return 0
    try:
        return int(round(float(value)))
    except ValueError:
        return 0


def _placement_for_decision(decision: str) -> str:
    """Map a Phase 1.5 routing decision to a ``Taxonomic_Placement`` value."""
    if decision == SSU_DECISION_SILVA:
        return PLACEMENT_PENDING_SILVA
    if decision == SSU_DECISION_GTDBTK:
        return PLACEMENT_PENDING_GTDBTK
    return PLACEMENT_UNKNOWN


def _build_placement_lookup(
    routing_log_path: Path | None,
) -> dict[str, str]:
    """Read the Phase 1.5 routing log and return a ``{accession: placement}`` map.

    Args:
        routing_log_path: Path to ``data/ssu/routing_log.tsv``; ``None``
            disables routing-aware placement and every MAG receives
            :data:`PLACEMENT_UNKNOWN`.

    Returns:
        Mapping from accession to the rendered ``Taxonomic_Placement`` value.
    """
    if routing_log_path is None:
        return {}
    rows = read_routing_log(routing_log_path)
    out: dict[str, str] = {}
    for row in rows:
        accession = row.get("accession", "")
        if not accession:
            continue
        out[accession] = _placement_for_decision(row.get("decision", ""))
    return out


def _format_frontmatter(fields: dict[str, object]) -> str:
    """Render the five-key frontmatter block.

    Differs from :func:`src.obsidian_exporter._format_frontmatter` in that
    the key order is fixed by :data:`FRONTMATTER_KEYS` and missing keys are
    rendered with empty values (rather than dropped). Booleans are written
    as ``true`` / ``false``; floats use two decimals; the
    ``FaRLiP_Cluster_Present`` slot is intentionally left blank when the
    caller does not supply a value.

    Args:
        fields: Mapping containing a subset of :data:`FRONTMATTER_KEYS`.

    Returns:
        A multi-line YAML frontmatter block (with surrounding ``---`` fences).
    """
    lines = ["---"]
    for key in FRONTMATTER_KEYS:
        if key not in fields or fields[key] is None:
            lines.append(f"{key}:")
            continue
        value = fields[key]
        if isinstance(value, bool):
            lines.append(f"{key}: {'true' if value else 'false'}")
        elif isinstance(value, float):
            lines.append(f"{key}: {value:.2f}")
        elif isinstance(value, int):
            lines.append(f"{key}: {value}")
        else:
            text = str(value)
            if text == "":
                lines.append(f"{key}:")
            else:
                lines.append(f"{key}: {_yaml_escape(text)}")
    lines.append("---")
    return "\n".join(lines)


def _format_mag_markdown(
    mag: _FilteredMag,
    placement: str,
) -> str:
    """Build the markdown body for one extremophile MAG note.

    Frontmatter is exactly the five mandated keys; the body adds a small
    context block (scientific name, QC tier, Phase 2 top classification)
    so the file is useful when read inside Obsidian without leaving the
    note. ``FaRLiP_Cluster_Present`` is deliberately left blank in the
    frontmatter - the user fills it in during manual review.

    Args:
        mag: A :class:`_FilteredMag` row.
        placement: Rendered ``Taxonomic_Placement`` value (may be empty).

    Returns:
        The full file contents (frontmatter + body).
    """
    frontmatter = _format_frontmatter(
        {
            "Accession": mag.accession,
            "Completeness": float(mag.completeness),
            "Contamination": float(mag.contamination),
            "Taxonomic_Placement": placement,
            "FaRLiP_Cluster_Present": None,
        }
    )
    body = [
        f"# {mag.accession}",
        "",
        (
            "High-quality extremophile MAG that passed the Phase 3 "
            "`recommended_for_phylogeny` gate "
            "(MIMAG HIGH/MEDIUM completeness/contamination + at least one "
            "Phase 2 FaRLiP candidate contig)."
        ),
        "",
        "## Quality summary",
        "",
        f"- **MIMAG tier:** `{mag.qc_quality or 'UNKNOWN'}`",
        f"- **Completeness:** {mag.completeness:.2f}%",
        f"- **Contamination:** {mag.contamination:.2f}%",
        f"- **Source:** {mag.source or 'unknown'}",
    ]
    if mag.scientific_name:
        body.append(f"- **Scientific name:** {mag.scientific_name}")
    body.extend(
        [
            "",
            "## FaRLiP context",
            "",
            (
                f"- **Phase 2 top classification:** "
                f"`{mag.farlip_top_classification or 'NONE'}`"
            ),
            (
                f"- **Distinct FaRLiP markers detected:** "
                f"{mag.n_farlip_distinct_markers} / 21"
            ),
            "",
            "Set `FaRLiP_Cluster_Present` in the frontmatter to `true` or "
            "`false` after manual review of the Phase 2 vault note "
            f"`[[{mag.accession}]]` (located in `data/obsidian_export/`).",
            "",
            "## Taxonomic placement",
            "",
            (
                f"`{placement}`" if placement
                else "_No Phase 1.5 routing decision recorded for this MAG._"
            ),
            "",
            "## See also",
            "",
            f"- [[{MASTER_INDEX_FILENAME[:-3]}]]",
            "",
        ]
    )
    return frontmatter + "\n\n" + "\n".join(body)


def _format_master_index(
    mags: list[_FilteredMag],
    placement_lookup: dict[str, str],
) -> str:
    """Build ``High_Quality_Extremophile_Genomes.md``.

    MAGs are listed alphabetically by accession with a short context blurb
    (completeness / contamination / placement) and a wikilink to the per-MAG
    note. A ``HIGH`` / ``MEDIUM`` subsection split is added when both tiers
    are present so the index is scannable in Obsidian's reading view.

    Args:
        mags: Filtered MAG records.
        placement_lookup: Pre-computed accession-to-placement mapping.

    Returns:
        The full master-index markdown.
    """
    lines: list[str] = [
        "# High Quality Extremophile Genomes",
        "",
        (
            f"Auto-generated by `src/extremophile_vault_exporter.py`. "
            f"Lists every MAG in `data/qc_filtered_mags.tsv` - i.e. every "
            f"assembly that satisfied the Phase 3 "
            f"`recommended_for_phylogeny` gate (MIMAG HIGH or MEDIUM "
            f"completeness/contamination + at least one Phase 2 FaRLiP "
            f"candidate contig)."
        ),
        "",
        f"**Total candidates:** {len(mags)}",
        "",
    ]

    by_quality: dict[str, list[_FilteredMag]] = {}
    for mag in mags:
        by_quality.setdefault(mag.qc_quality or "UNKNOWN", []).append(mag)

    tiers_present = [t for t in ("HIGH", "MEDIUM", "LOW", "UNKNOWN") if t in by_quality]
    if len(tiers_present) > 1:
        for tier in tiers_present:
            tier_rows = sorted(by_quality[tier], key=lambda m: m.accession)
            lines.append(f"## {tier}")
            lines.append("")
            for mag in tier_rows:
                lines.append(_master_index_line(mag, placement_lookup))
            lines.append("")
    else:
        for mag in sorted(mags, key=lambda m: m.accession):
            lines.append(_master_index_line(mag, placement_lookup))
        lines.append("")

    return "\n".join(lines)


def _master_index_line(mag: _FilteredMag, placement_lookup: dict[str, str]) -> str:
    """Render one bullet line in the master index."""
    placement = placement_lookup.get(mag.accession, "")
    placement_blurb = f" - {placement}" if placement else ""
    return (
        f"- [[{mag.accession}]] - "
        f"{mag.completeness:.1f}% / {mag.contamination:.1f}%"
        f"{placement_blurb}"
    )


def export_extremophile_vault(
    filtered_csv_path: Path,
    vault_dir: Path,
    *,
    routing_log_path: Path | None = None,
) -> int:
    """Render the high-quality extremophile vault.

    Args:
        filtered_csv_path: Path to ``data/qc_filtered_mags.tsv``.
        vault_dir: Output directory (typically ``data/extremophile_vault``).
            Created on demand. Existing per-MAG files are overwritten.
        routing_log_path: Optional path to ``data/ssu/routing_log.tsv``.
            When supplied, ``Taxonomic_Placement`` is populated from the
            routing decisions; otherwise the slot is left empty.

    Returns:
        Number of per-MAG markdown files written (master index excluded).
    """
    mags = _read_filtered_mags(filtered_csv_path)
    placement_lookup = _build_placement_lookup(routing_log_path)

    vault_dir.mkdir(parents=True, exist_ok=True)
    written = 0
    for mag in mags:
        placement = placement_lookup.get(mag.accession, PLACEMENT_UNKNOWN)
        md = _format_mag_markdown(mag, placement)
        out_path = vault_dir / f"{mag.accession}.md"
        out_path.write_text(md, encoding="utf-8")
        written += 1

    index_md = _format_master_index(mags, placement_lookup)
    (vault_dir / MASTER_INDEX_FILENAME).write_text(index_md, encoding="utf-8")
    return written


def list_filtered_accessions(filtered_csv_path: Path) -> list[str]:
    """Convenience helper - return just the accessions in the filtered TSV.

    Used by the Phase 3 CLI to print a one-line summary after the export.

    Args:
        filtered_csv_path: Path to ``data/qc_filtered_mags.tsv``.

    Returns:
        Sorted list of accessions (empty when the file is absent).
    """
    return sorted(m.accession for m in _read_filtered_mags(filtered_csv_path))


def _accession_iter(rows: Iterable[_FilteredMag]) -> Iterable[str]:
    """Tiny iterator that yields accessions from filtered rows."""
    for r in rows:
        yield r.accession


__all__ = [
    "FRONTMATTER_KEYS",
    "MASTER_INDEX_FILENAME",
    "PLACEMENT_PENDING_GTDBTK",
    "PLACEMENT_PENDING_SILVA",
    "PLACEMENT_UNKNOWN",
    "VAULT_DIRNAME",
    "export_extremophile_vault",
    "list_filtered_accessions",
]
