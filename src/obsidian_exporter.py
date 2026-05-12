"""Export the per-MAG FaRLiP candidate report into an Obsidian-friendly markdown vault.

Each candidate MAG (``n_candidate_contigs >= 1`` in ``farlip_hits_per_mag.tsv``) gets
its own ``<accession>.md`` with a YAML frontmatter block (parseable by ``yaml.safe_load``)
plus a markdown table of every per-contig hit. A ``Phase_3_Master_Index.md`` at the
vault root groups MAGs by classification and links to each candidate file via
wikilink syntax (``[[<accession>]]``) so the vault is browsable directly in
Obsidian's default reading view.

The exporter is purely a transformation step on the TSV reports written by
:mod:`src.synteny_parser` plus the MAG metadata in ``data/manifest.tsv``. It does
no biology of its own.

Phase 4 extension:
    Passing ``qc_report_path`` to :func:`export_vault` joins the QC report
    written by :mod:`src.qc_reporter` so the per-MAG frontmatter gains
    ``completeness``, ``contamination``, ``n50_contig``, ``genome_size_bp``,
    ``qc_quality``, ``qc_source`` and ``recommended_for_phylogeny``, and the
    master index gains a ``## Quality`` subsection grouping MAGs by tier.
    Backward-compatible: omitting ``qc_report_path`` keeps the Phase 3 layout
    untouched.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


MASTER_INDEX_FILENAME: str = "Phase_3_Master_Index.md"
"""Filename for the top-level vault index. Convention: ``Phase_3_*`` so future
phases get their own indices without renaming this one."""

YAML_FENCE: str = "---"
"""Markdown YAML frontmatter fence. Obsidian recognizes this for inline metadata."""

CLASSIFICATION_ORDER: tuple[str, ...] = (
    "FULL_CLUSTER",
    "EDGE_TRUNCATED",
    "PARTIAL_CLUSTER",
    "NONE",
)
"""Display order for the master index sections (most-promising-first)."""

QC_QUALITY_ORDER: tuple[str, ...] = (
    "HIGH",
    "MEDIUM",
    "LOW",
)
"""Display order for the master-index Quality subsection (best-first)."""

HITS_TABLE_HEADER: tuple[str, ...] = (
    "contig_id",
    "contig_length",
    "n_markers",
    "marker_names",
    "min_evalue",
    "edge_truncated",
    "classification",
)
"""Columns surfaced in the per-MAG markdown hits table. A subset of the per-contig
TSV - we drop ``marker_positions`` because the long position list is not useful in
a human-facing markdown table (the full TSV remains the source of truth for it)."""


@dataclass(frozen=True)
class _ContigRow:
    """One row parsed out of ``farlip_hits_report.tsv``."""

    mag_accession: str
    contig_id: str
    contig_length: int
    n_markers: int
    marker_names: str
    marker_positions: str
    min_evalue: str
    edge_truncated: bool
    classification: str


@dataclass(frozen=True)
class _PerMagRow:
    """One row parsed out of ``farlip_hits_per_mag.tsv``."""

    mag_accession: str
    n_candidate_contigs: int
    n_total_markers: int
    n_distinct_markers: int
    has_edge_truncated: bool
    top_classification: str


@dataclass(frozen=True)
class _QcRow:
    """One row parsed out of ``data/qc_report.tsv`` (Phase 4 join).

    Attributes:
        accession: MAG accession.
        qc_source: Source of the QC numbers (``pyhmmer`` or ``checkm2``).
        marker_set_label: Human-readable marker-set name.
        completeness: 0-100.
        contamination: 0-100.
        n50_contig: Contig N50 in bp.
        genome_size_bp: Total assembly size.
        qc_quality: ``HIGH`` / ``MEDIUM`` / ``LOW``.
        recommended_for_phylogeny: Pre-computed gate for Phase 5.
    """

    accession: str
    qc_source: str
    marker_set_label: str
    completeness: float
    contamination: float
    n50_contig: int
    genome_size_bp: int
    qc_quality: str
    recommended_for_phylogeny: bool


def _read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into a list of header-keyed row dicts.

    Args:
        path: TSV file to read. Returns an empty list if missing.

    Returns:
        Rows as dicts (empty list if file absent).
    """
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def _parse_contig_rows(rows: list[dict[str, str]]) -> dict[str, list[_ContigRow]]:
    """Group ``farlip_hits_report.tsv`` rows by MAG accession.

    Args:
        rows: Output of ``_read_tsv`` on the per-contig report.

    Returns:
        Mapping of MAG accession -> list of contig rows.
    """
    by_mag: dict[str, list[_ContigRow]] = {}
    for row in rows:
        cr = _ContigRow(
            mag_accession=row["mag_accession"],
            contig_id=row["contig_id"],
            contig_length=int(row["contig_length"]),
            n_markers=int(row["n_markers"]),
            marker_names=row["marker_names"],
            marker_positions=row["marker_positions"],
            min_evalue=row["min_evalue"],
            edge_truncated=row["edge_truncated"].lower() == "true",
            classification=row["classification"],
        )
        by_mag.setdefault(cr.mag_accession, []).append(cr)
    return by_mag


def _parse_per_mag_rows(rows: list[dict[str, str]]) -> dict[str, _PerMagRow]:
    """Index ``farlip_hits_per_mag.tsv`` rows by MAG accession.

    Args:
        rows: Output of ``_read_tsv`` on the per-MAG report.

    Returns:
        Mapping of MAG accession -> :class:`_PerMagRow`.
    """
    out: dict[str, _PerMagRow] = {}
    for row in rows:
        pmr = _PerMagRow(
            mag_accession=row["mag_accession"],
            n_candidate_contigs=int(row["n_candidate_contigs"]),
            n_total_markers=int(row["n_total_markers"]),
            n_distinct_markers=int(row["n_distinct_markers"]),
            has_edge_truncated=row["has_edge_truncated"].lower() == "true",
            top_classification=row["top_classification"],
        )
        out[pmr.mag_accession] = pmr
    return out


def _parse_manifest_lookup(manifest_path: Path) -> dict[str, dict[str, str]]:
    """Index ``data/manifest.tsv`` by accession.

    Args:
        manifest_path: Path to the Phase 1 manifest TSV.

    Returns:
        Mapping of accession -> full row dict (or empty dict if manifest absent).
    """
    rows = _read_tsv(manifest_path)
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        acc = row.get("accession")
        if not acc:
            continue
        out[acc] = row
    return out


def _parse_qc_lookup(qc_report_path: Path | None) -> dict[str, _QcRow]:
    """Index ``data/qc_report.tsv`` by accession.

    Args:
        qc_report_path: Path to the Phase 4 QC report (or ``None`` to opt out
            of QC-aware rendering).

    Returns:
        Mapping of accession -> :class:`_QcRow`. Empty dict when the path is
        ``None`` or the file does not exist.
    """
    if qc_report_path is None:
        return {}
    rows = _read_tsv(qc_report_path)
    out: dict[str, _QcRow] = {}
    for row in rows:
        accession = row.get("accession", "")
        if not accession:
            continue
        try:
            completeness = float(row.get("completeness", "0") or 0.0)
        except ValueError:
            completeness = 0.0
        try:
            contamination = float(row.get("contamination", "0") or 0.0)
        except ValueError:
            contamination = 0.0
        try:
            n50 = int(row.get("n50_contig", "0") or 0)
        except ValueError:
            n50 = 0
        try:
            genome_size = int(row.get("genome_size_bp", "0") or 0)
        except ValueError:
            genome_size = 0
        out[accession] = _QcRow(
            accession=accession,
            qc_source=row.get("qc_source", ""),
            marker_set_label=row.get("marker_set_label", ""),
            completeness=completeness,
            contamination=contamination,
            n50_contig=n50,
            genome_size_bp=genome_size,
            qc_quality=row.get("qc_quality", ""),
            recommended_for_phylogeny=(
                row.get("recommended_for_phylogeny", "false").lower() == "true"
            ),
        )
    return out


def _yaml_escape(value: str) -> str:
    """Quote a YAML scalar if it contains characters that would be ambiguous unquoted.

    Args:
        value: Scalar string to emit.

    Returns:
        Either the bare value or a single-quoted version with internal quotes
        doubled per YAML 1.2.
    """
    if value == "":
        return '""'
    if any(c in value for c in ":#{}[],&*!|>%@`?\n\t"):
        return "'" + value.replace("'", "''") + "'"
    return value


def _format_frontmatter(fields: dict[str, object]) -> str:
    """Render a flat ``{key: value}`` mapping as a YAML frontmatter block.

    Args:
        fields: Mapping. Lists are rendered as flow-style sequences; booleans as
            ``true``/``false``; everything else is coerced to string and escaped.

    Returns:
        Multi-line string starting and ending with ``---`` fences.
    """
    out = [YAML_FENCE]
    for key, val in fields.items():
        if isinstance(val, bool):
            out.append(f"{key}: {'true' if val else 'false'}")
        elif isinstance(val, list):
            inner = ", ".join(_yaml_escape(str(v)) for v in val)
            out.append(f"{key}: [{inner}]")
        elif val is None:
            out.append(f"{key}: null")
        else:
            out.append(f"{key}: {_yaml_escape(str(val))}")
    out.append(YAML_FENCE)
    return "\n".join(out)


def _escape_md_cell(value: str) -> str:
    """Escape characters that would otherwise break a GFM table cell.

    Pipe characters appear in ENA-style contig IDs (``ENA|GCA_XXXX|contig.1``) and
    would prematurely close the table cell. GFM specifies ``\\|`` as the escape.

    Args:
        value: Cell text.

    Returns:
        The same text with pipes escaped.
    """
    return value.replace("|", "\\|")


def _format_hits_table(contig_rows: Iterable[_ContigRow]) -> str:
    """Render the per-MAG contig table in GitHub-flavored markdown.

    Args:
        contig_rows: Iterable of contig rows to include in the table.

    Returns:
        Multi-line markdown string. Returns an empty-table sentinel if no rows.
    """
    rows = list(contig_rows)
    if not rows:
        return "_No candidate contigs were retained for this MAG._"
    lines = ["| " + " | ".join(HITS_TABLE_HEADER) + " |"]
    lines.append("|" + "|".join(["---"] * len(HITS_TABLE_HEADER)) + "|")
    for r in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    _escape_md_cell(r.contig_id),
                    str(r.contig_length),
                    str(r.n_markers),
                    _escape_md_cell(r.marker_names),
                    r.min_evalue,
                    "true" if r.edge_truncated else "false",
                    r.classification,
                ]
            )
            + " |"
        )
    return "\n".join(lines)


def _format_mag_markdown(
    mag_accession: str,
    per_mag: _PerMagRow,
    contig_rows: list[_ContigRow],
    manifest_row: dict[str, str],
    qc_row: _QcRow | None = None,
) -> str:
    """Build the markdown body for one candidate MAG.

    Args:
        mag_accession: MAG accession key.
        per_mag: The roll-up row for this MAG.
        contig_rows: All per-contig rows for this MAG (already filtered to
            candidates).
        manifest_row: The Phase 1 manifest row for this MAG (may be empty if the
            MAG was scanned without a matching manifest entry, e.g. JSC-1).
        qc_row: Optional Phase 4 QC row. When supplied, completeness /
            contamination / N50 / genome size / quality fields are added to
            the YAML frontmatter and a "Quality assessment" section is added
            to the body.

    Returns:
        The full file contents (frontmatter + body) ready to write.
    """
    contig_ids = [r.contig_id for r in contig_rows]
    scientific_name = manifest_row.get("scientific_name") or manifest_row.get("label") or ""
    source = manifest_row.get("source") or ""
    sample_accession = manifest_row.get("sample_accession") or ""
    frontmatter_fields: dict[str, object] = {
        "accession": mag_accession,
        "scientific_name": scientific_name,
        "source": source,
        "sample_accession": sample_accession,
        "n_farlip_markers": per_mag.n_distinct_markers,
        "n_candidate_contigs": per_mag.n_candidate_contigs,
        "contig_ids": contig_ids,
        "classification": per_mag.top_classification,
        "edge_truncated": per_mag.has_edge_truncated,
        "manual_review": per_mag.top_classification != "FULL_CLUSTER",
    }
    if qc_row is not None:
        frontmatter_fields.update(
            {
                "completeness": round(qc_row.completeness, 2),
                "contamination": round(qc_row.contamination, 2),
                "n50_contig": qc_row.n50_contig,
                "genome_size_bp": qc_row.genome_size_bp,
                "qc_quality": qc_row.qc_quality,
                "qc_source": qc_row.qc_source,
                "recommended_for_phylogeny": qc_row.recommended_for_phylogeny,
            }
        )

    body = [
        f"# {mag_accession} - FaRLiP candidate",
        "",
        (
            f"**Top classification:** `{per_mag.top_classification}`  "
            f"**Distinct markers:** {per_mag.n_distinct_markers} / 21  "
            f"**Candidate contigs:** {per_mag.n_candidate_contigs}  "
            f"**Edge-truncated:** {'yes' if per_mag.has_edge_truncated else 'no'}"
        ),
        "",
    ]

    if qc_row is not None:
        body.extend(
            [
                "## Quality assessment",
                "",
                (
                    f"**Quality tier:** `{qc_row.qc_quality}`  "
                    f"**Source:** {qc_row.qc_source or 'unknown'}  "
                    f"**Completeness:** {qc_row.completeness:.2f}%  "
                    f"**Contamination:** {qc_row.contamination:.2f}%  "
                    f"**N50:** {qc_row.n50_contig:,} bp  "
                    f"**Genome size:** {qc_row.genome_size_bp:,} bp  "
                    f"**Recommended for phylogeny:** "
                    f"{'yes' if qc_row.recommended_for_phylogeny else 'no'}"
                ),
                "",
            ]
        )

    body.extend(
        [
            "## Per-contig FaRLiP hits",
            "",
            _format_hits_table(contig_rows),
            "",
            "## See also",
            "",
            f"- [[{MASTER_INDEX_FILENAME[:-3]}]]",
            "",
        ]
    )
    return _format_frontmatter(frontmatter_fields) + "\n\n" + "\n".join(body)


def _format_master_index(
    per_mag_rows: dict[str, _PerMagRow],
    qc_lookup: dict[str, _QcRow] | None = None,
) -> str:
    """Build the contents of ``Phase_3_Master_Index.md``.

    MAGs are grouped under headings by classification, in
    :data:`CLASSIFICATION_ORDER` (most promising first). Within each section
    MAGs are sorted by distinct marker count descending, then accession
    ascending, so the most-evidenced candidates appear at the top.

    When ``qc_lookup`` is supplied, an additional ``## Quality`` section is
    appended grouping every MAG with a QC row by tier (HIGH / MEDIUM / LOW),
    so the same index doubles as a Phase 4 quality dashboard.

    Args:
        per_mag_rows: Per-MAG roll-up indexed by accession.
        qc_lookup: Optional Phase 4 QC rows by accession.

    Returns:
        Markdown string.
    """
    by_class: dict[str, list[_PerMagRow]] = {c: [] for c in CLASSIFICATION_ORDER}
    for row in per_mag_rows.values():
        by_class.setdefault(row.top_classification, []).append(row)
    lines = [
        "# Phase 3 Master Index - FaRLiP candidate MAGs",
        "",
        (
            "Auto-generated by `src/obsidian_exporter.py`. Re-run "
            "`python -m src.pipeline_phase3 export-vault` (or "
            "`python -m src.pipeline_phase4 update-vault` for the QC-aware "
            "render) to refresh."
        ),
        "",
    ]
    for classification in CLASSIFICATION_ORDER:
        section_rows = by_class.get(classification) or []
        if not section_rows:
            continue
        section_rows.sort(key=lambda r: (-r.n_distinct_markers, r.mag_accession))
        lines.append(f"## {classification}")
        lines.append("")
        for r in section_rows:
            edge = " (edge-truncated)" if r.has_edge_truncated else ""
            qc = (qc_lookup or {}).get(r.mag_accession)
            qc_blurb = ""
            if qc is not None and qc.qc_quality:
                qc_blurb = (
                    f" - QC `{qc.qc_quality}` "
                    f"({qc.completeness:.1f}% / {qc.contamination:.1f}%)"
                )
            lines.append(
                f"- [[{r.mag_accession}]] - {r.n_distinct_markers}/21 markers"
                f" across {r.n_candidate_contigs} contig(s){edge}{qc_blurb}"
            )
        lines.append("")

    if qc_lookup:
        lines.append("## Quality")
        lines.append("")
        by_quality: dict[str, list[_QcRow]] = {q: [] for q in QC_QUALITY_ORDER}
        for qc in qc_lookup.values():
            by_quality.setdefault(qc.qc_quality or "LOW", []).append(qc)
        for tier in QC_QUALITY_ORDER:
            tier_rows = by_quality.get(tier) or []
            if not tier_rows:
                continue
            tier_rows.sort(key=lambda r: (-r.completeness, r.accession))
            lines.append(f"### {tier}")
            lines.append("")
            for qc in tier_rows:
                rec = " (recommended for phylogeny)" if qc.recommended_for_phylogeny else ""
                lines.append(
                    f"- [[{qc.accession}]] - {qc.completeness:.1f}% complete, "
                    f"{qc.contamination:.1f}% contamination, "
                    f"N50 {qc.n50_contig:,} bp{rec}"
                )
            lines.append("")
    return "\n".join(lines)


def export_vault(
    report_path: Path,
    per_mag_path: Path,
    manifest_path: Path,
    vault_dir: Path,
    qc_report_path: Path | None = None,
) -> int:
    """Generate the full Obsidian vault from the Phase 3 reports.

    The vault directory is created if missing. Existing files are overwritten
    so the export is idempotent.

    Args:
        report_path: ``data/farlip_hits_report.tsv``.
        per_mag_path: ``data/farlip_hits_per_mag.tsv``.
        manifest_path: ``data/manifest.tsv``.
        vault_dir: Output directory (typically ``data/obsidian_export``).
        qc_report_path: Optional ``data/qc_report.tsv`` produced by
            :func:`src.qc_reporter.write_qc_reports`. When supplied, every
            per-MAG markdown gains QC frontmatter fields and the master index
            gains a ``## Quality`` section. When ``None``, behaviour is
            identical to the Phase 3 baseline.

    Returns:
        Number of per-MAG markdown files written (not counting the master index).
    """
    vault_dir.mkdir(parents=True, exist_ok=True)

    contig_rows_by_mag = _parse_contig_rows(_read_tsv(report_path))
    per_mag_rows = _parse_per_mag_rows(_read_tsv(per_mag_path))
    manifest_lookup = _parse_manifest_lookup(manifest_path)
    qc_lookup = _parse_qc_lookup(qc_report_path)

    written = 0
    for mag_accession, per_mag in per_mag_rows.items():
        if per_mag.n_candidate_contigs == 0:
            continue
        contig_rows = contig_rows_by_mag.get(mag_accession, [])
        contig_rows.sort(key=lambda r: -r.n_markers)
        manifest_row = manifest_lookup.get(mag_accession, {})
        md = _format_mag_markdown(
            mag_accession,
            per_mag,
            contig_rows,
            manifest_row,
            qc_row=qc_lookup.get(mag_accession),
        )
        out_path = vault_dir / f"{mag_accession}.md"
        out_path.write_text(md, encoding="utf-8")
        written += 1

    index_md = _format_master_index(per_mag_rows, qc_lookup=qc_lookup or None)
    (vault_dir / MASTER_INDEX_FILENAME).write_text(index_md, encoding="utf-8")

    return written
