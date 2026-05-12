"""Phase 2 SILVA SINA bundle + routing log writer.

Two responsibilities:

1. **Concatenate** every per-MAG ``data/ssu/<accession>.16S.fasta`` whose
   :func:`src.ssu_extractor.assess_mag_for_ssu` decision was
   :data:`src.ssu_extractor.SSU_DECISION_SILVA` into a single multi-FASTA
   ``data/ssu/silva_sina_input.fasta``. That file is the input to a SILVA
   SINA submission (``https://www.arb-silva.de/aligner/``) - SINA is the
   only step that actually contacts SILVA, and the user runs it themselves
   (web upload or the SINA CLI). This module never touches the network.

2. **Write the routing log** ``data/ssu/routing_log.tsv``. The schema is
   defined by :data:`src.ssu_extractor.ROUTING_LOG_COLUMNS`. The log is the
   single source of truth for which MAG goes to which downstream tool, and
   is consumed verbatim by :mod:`src.gtdbtk_shellscript_generator` and the
   :mod:`src.extremophile_vault_exporter`.

The bundler is a pure transformation - it never re-runs barrnap, never
re-extracts a 16S record, and never deletes anything. Re-running it is
idempotent: existing routing-log rows for a given accession are replaced in
place, the SILVA bundle is rebuilt from whatever SSU FASTAs are currently on
disk.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from .ssu_extractor import (
    ROUTING_LOG_COLUMNS,
    SSU_DECISION_GTDBTK,
    SSU_DECISION_SILVA,
    SsuAssessment,
)


SILVA_SINA_BUNDLE_FILENAME: str = "silva_sina_input.fasta"
"""Filename for the concatenated multi-FASTA. Lives at
``<data_root>/ssu/silva_sina_input.fasta``."""

ROUTING_LOG_FILENAME: str = "routing_log.tsv"
"""Filename for the routing log. Lives at ``<data_root>/ssu/routing_log.tsv``."""


@dataclass(frozen=True)
class RoutingDecision:
    """One row destined for ``data/ssu/routing_log.tsv``.

    Attributes:
        accession: MAG accession.
        source: ``NCBI`` / ``ENA`` / ... (copied from the Phase 1 manifest).
        fasta_path: Path to the assembly FASTA the SSU extractor scanned.
        n_contigs_scanned: From :class:`SsuAssessment`.
        n_16S_predicted: From :class:`SsuAssessment`.
        decision: :data:`SSU_DECISION_SILVA` or :data:`SSU_DECISION_GTDBTK`.
        output_path: For SILVA-routed MAGs, the per-MAG 16S FASTA. For
            GTDB-Tk-routed MAGs, the assembly FASTA (since that is what
            gets handed to ``gtdbtk classify_wf``).
        barrnap_version: Captured from ``barrnap --version`` (best effort).
        notes: Free-text. The bundler appends short status notes here
            (e.g. ``"barrnap binary missing"`` for failures the CLI surfaced).
    """

    accession: str
    source: str
    fasta_path: Path
    n_contigs_scanned: int
    n_16S_predicted: int
    decision: str
    output_path: Path
    barrnap_version: str = ""
    notes: str = ""


@dataclass
class BundleResult:
    """Outcome of one :func:`bundle_silva_sina` call.

    Attributes:
        bundle_path: Path the multi-FASTA was written to.
        n_mags_included: Number of MAGs whose 16S FASTA contributed at
            least one record.
        n_records_total: Total number of FASTA records emitted.
        skipped_accessions: MAGs whose decision was ``silva`` but whose
            per-MAG 16S FASTA was missing or empty (typical when the
            extraction stamp is stale).
    """

    bundle_path: Path
    n_mags_included: int = 0
    n_records_total: int = 0
    skipped_accessions: list[str] = field(default_factory=list)


def routing_decision_from_assessment(
    assessment: SsuAssessment,
    *,
    source: str = "",
    notes: str = "",
) -> RoutingDecision:
    """Build a :class:`RoutingDecision` from a Phase 2 SSU assessment.

    The ``output_path`` field is set to the per-MAG 16S FASTA for SILVA-routed
    MAGs and to the input assembly FASTA otherwise; the GTDB-Tk script
    generator joins on this column to assemble its batchfile.

    Args:
        assessment: An :class:`SsuAssessment` returned by
            :func:`src.ssu_extractor.assess_mag_for_ssu`.
        source: Manifest ``source`` column for the routing log row.
        notes: Optional free-text annotation.

    Returns:
        A populated :class:`RoutingDecision`.
    """
    if assessment.decision == SSU_DECISION_SILVA:
        out_path = assessment.ssu_fasta_path
    else:
        out_path = assessment.fasta_path
    return RoutingDecision(
        accession=assessment.accession,
        source=source,
        fasta_path=assessment.fasta_path,
        n_contigs_scanned=assessment.n_contigs_scanned,
        n_16S_predicted=assessment.n_16S_predicted,
        decision=assessment.decision,
        output_path=out_path,
        barrnap_version=assessment.barrnap_version,
        notes=notes,
    )


def write_routing_log(rows: Iterable[RoutingDecision], path: Path) -> int:
    """Write or rewrite the Phase 2 routing log.

    The log is overwritten each call (not appended) so the contents always
    reflect the *current* state of ``data/ssu/``. Callers who need to merge
    new rows with old ones should call :func:`upsert_routing_log` instead.

    Args:
        rows: Iterable of :class:`RoutingDecision` to write.
        path: Output TSV path. Parent directories are created on demand.

    Returns:
        The number of rows written.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    rows_list = list(rows)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=ROUTING_LOG_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for r in rows_list:
            writer.writerow(_row_to_dict(r))
    return len(rows_list)


def read_routing_log(path: Path) -> list[dict[str, str]]:
    """Read the routing log into a list of header-keyed dicts.

    Args:
        path: TSV path; an empty list is returned when the file is missing.

    Returns:
        Raw rows; numeric fields stay as strings (parse on demand).
    """
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def upsert_routing_log(path: Path, new_rows: Iterable[RoutingDecision]) -> int:
    """Merge new rows into an existing routing log.

    Existing rows whose ``accession`` matches an incoming row are replaced;
    other rows are preserved. New rows are appended in iteration order.

    Args:
        path: Routing log path. Created if missing.
        new_rows: Iterable of :class:`RoutingDecision`.

    Returns:
        Total number of rows in the merged log.
    """
    incoming = list(new_rows)
    keyed: dict[str, dict[str, str]] = {}
    for existing in read_routing_log(path):
        accession = existing.get("accession", "")
        if accession:
            keyed[accession] = existing
    for r in incoming:
        keyed[r.accession] = _row_to_dict(r)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=ROUTING_LOG_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in keyed.values():
            writer.writerow({col: row.get(col, "") for col in ROUTING_LOG_COLUMNS})
    return len(keyed)


def _row_to_dict(r: RoutingDecision) -> dict[str, str]:
    """Flatten a :class:`RoutingDecision` into a TSV-ready dict."""
    return {
        "accession": r.accession,
        "source": r.source,
        "fasta_path": str(r.fasta_path),
        "n_contigs_scanned": str(r.n_contigs_scanned),
        "n_16S_predicted": str(r.n_16S_predicted),
        "decision": r.decision,
        "output_path": str(r.output_path),
        "barrnap_version": r.barrnap_version,
        "notes": r.notes,
    }


def bundle_silva_sina(
    routing_log_path: Path,
    ssu_dir: Path,
    output_path: Path,
) -> BundleResult:
    """Concatenate every SILVA-routed MAG's 16S FASTA into a single multi-FASTA.

    The file is rebuilt from scratch each call. Headers are preserved verbatim
    from the per-MAG FASTAs (which already prefix the MAG accession), so the
    bundle is suitable for direct upload to SILVA SINA.

    Args:
        routing_log_path: Path to ``data/ssu/routing_log.tsv``.
        ssu_dir: Directory containing the per-MAG ``.16S.fasta`` files.
        output_path: Where to write the bundle (typically
            ``ssu_dir / "silva_sina_input.fasta"``).

    Returns:
        A :class:`BundleResult` describing what landed in the bundle.
    """
    rows = read_routing_log(routing_log_path)
    silva_rows = [r for r in rows if r.get("decision") == SSU_DECISION_SILVA]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_mags = 0
    n_records = 0
    skipped: list[str] = []
    with output_path.open("w", encoding="utf-8") as out_fh:
        for row in silva_rows:
            accession = row.get("accession", "")
            per_mag_path = Path(row.get("output_path", ""))
            if not per_mag_path.is_absolute() and not per_mag_path.exists():
                # routing log may carry relative paths in tests
                per_mag_path = (ssu_dir / per_mag_path.name).resolve()
            if not per_mag_path.exists() or per_mag_path.stat().st_size == 0:
                skipped.append(accession)
                continue
            n_in_file = 0
            with per_mag_path.open("r", encoding="utf-8") as in_fh:
                for line in in_fh:
                    out_fh.write(line)
                    if line.startswith(">"):
                        n_in_file += 1
            if n_in_file == 0:
                skipped.append(accession)
                continue
            n_mags += 1
            n_records += n_in_file
    return BundleResult(
        bundle_path=output_path,
        n_mags_included=n_mags,
        n_records_total=n_records,
        skipped_accessions=skipped,
    )


__all__ = [
    "BundleResult",
    "ROUTING_LOG_FILENAME",
    "RoutingDecision",
    "SILVA_SINA_BUNDLE_FILENAME",
    "bundle_silva_sina",
    "read_routing_log",
    "routing_decision_from_assessment",
    "upsert_routing_log",
    "write_routing_log",
]
