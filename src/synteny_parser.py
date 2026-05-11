"""Synteny-aware classification of FaRLiP HMM hits.

Joins each MAG's ``<accession>.domtblout`` (from :mod:`src.hmm_searcher`) with the
corresponding ``<accession>.coords.tsv`` (from :mod:`src.orf_predictor`) and applies
two rules:

* **Candidate rule.** Any window of :data:`SYNTENY_WINDOW_BP` that contains
  :data:`MIN_MARKERS_FOR_CANDIDATE` or more distinct FaRLiP marker HMMs marks the
  hosting contig as a candidate. Each ORF is allowed to count for at most one
  marker - the HMM with the lowest i-Evalue wins ties.
* **Edge-truncation rule.** A candidate marker lying within
  :data:`EDGE_PROXIMITY_BP` of either end of its host contig flags the contig as
  edge-truncated. This is the key affordance for fragmented MAGs (a partial
  cluster split across multiple contigs is a high-priority manual-review hit, not
  noise).

Outputs:

* ``data/farlip_hits_report.tsv``  - one row per **candidate contig**.
* ``data/farlip_hits_per_mag.tsv`` - one row per MAG (rolled up).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path


SYNTENY_WINDOW_BP: int = 10_000
"""Sliding-window size for the **candidate trigger** (Phase 2 spec). A contig is
a FaRLiP candidate only if some 10 kb window contains
:data:`MIN_MARKERS_FOR_CANDIDATE`+ distinct HMM names. Tight by design - this
filters out chromosome-wide scattered paralog hits."""

MIN_MARKERS_FOR_CANDIDATE: int = 3
"""Minimum distinct FaRLiP marker HMMs that must co-localize inside one
:data:`SYNTENY_WINDOW_BP` window to flag a contig as a candidate."""

CLUSTER_SPAN_BP: int = 30_000
"""Window size used for **classification** of an already-flagged candidate.

The canonical FaRLiP cluster spans ~26 kb (Gan et al. 2014; JSC-1 cluster is
1,095,910 - 1,121,582 on NZ_KL662192.1 = 25.7 kb). A 30 kb window comfortably
captures the full cluster on a clean contig while still excluding scattered
paralog hits elsewhere on the chromosome. This is the window inside which we
count distinct HMM names for the ``n_markers`` field."""

EDGE_PROXIMITY_BP: int = 2_000
"""Distance from either contig end inside which a marker triggers the
``EDGE_TRUNCATED`` flag (Phase 2 spec)."""

FULL_CLUSTER_MARKER_THRESHOLD: int = 15
"""Markers (distinct HMMs) required inside the densest :data:`CLUSTER_SPAN_BP`
window to classify a contig as ``FULL_CLUSTER``. Derived: 15 / 21 markers is
the typical threshold quoted in the FaRLiP literature for "intact" cluster
calls."""

CONTIG_REPORT_COLUMNS: tuple[str, ...] = (
    "mag_accession",
    "contig_id",
    "contig_length",
    "n_markers",
    "n_total_hits_on_contig",
    "marker_names",
    "marker_positions",
    "min_evalue",
    "edge_truncated",
    "classification",
)
"""Column order for ``data/farlip_hits_report.tsv``.

* ``n_markers`` is the **distinct-HMM** count inside the densest synteny window
  (the value used for classification).
* ``n_total_hits_on_contig`` is the total ORF-level hit count anywhere on the
  contig - useful for diagnosing scattered-paralog false positives."""

PER_MAG_REPORT_COLUMNS: tuple[str, ...] = (
    "mag_accession",
    "n_candidate_contigs",
    "n_total_markers",
    "n_distinct_markers",
    "has_edge_truncated",
    "top_classification",
)
"""Column order for ``data/farlip_hits_per_mag.tsv``."""


class Classification(str, Enum):
    """Per-contig FaRLiP classification tag."""

    FULL_CLUSTER = "FULL_CLUSTER"
    PARTIAL_CLUSTER = "PARTIAL_CLUSTER"
    EDGE_TRUNCATED = "EDGE_TRUNCATED"
    NONE = "NONE"


_CLASSIFICATION_PRIORITY: dict[Classification, int] = {
    Classification.FULL_CLUSTER: 3,
    Classification.EDGE_TRUNCATED: 2,
    Classification.PARTIAL_CLUSTER: 1,
    Classification.NONE: 0,
}
"""Used to pick a MAG-level summary classification from a list of contig calls."""


@dataclass
class _OrfCoord:
    """One row from ``<accession>.coords.tsv``, lightly typed."""

    contig_id: str
    start: int
    end: int
    strand: str
    contig_length: int


@dataclass
class _BestHit:
    """Best (lowest i-Evalue) HMM hit recorded for one ORF."""

    orf_id: str
    hmm_name: str
    i_evalue: float
    full_evalue: float
    score: float


@dataclass
class ContigCall:
    """Per-contig synteny result before TSV serialization.

    ``markers`` is the subset of hits that fall inside the **densest** synteny
    window on the contig - that is, the 10 kb window with the largest number of
    distinct FaRLiP HMM names. Counting markers this way enforces the
    "co-localized cluster" intent of the FULL_CLUSTER threshold; scattered hits
    from promiscuous HMMs (e.g. the rfpA phytochrome cross-matching unrelated
    GAF-domain regulators across the chromosome) are deliberately excluded.
    """

    mag_accession: str
    contig_id: str
    contig_length: int
    markers: list[tuple[str, int, int, float]]
    """Each tuple is (hmm_name, start, end, i_evalue), restricted to the densest
    synteny window and sorted by start."""
    n_total_hits_on_contig: int
    """Number of ORF hits anywhere on the contig (informational - includes scattered
    paralog hits that are not counted toward ``n_markers``)."""
    edge_truncated: bool
    classification: Classification

    @property
    def n_markers(self) -> int:
        """Number of distinct FaRLiP marker HMMs inside the densest window."""
        return len({m[0] for m in self.markers})

    @property
    def min_evalue(self) -> float:
        """Smallest (best) i-Evalue across markers inside the densest window."""
        return min(m[3] for m in self.markers) if self.markers else float("inf")


def parse_coords_tsv(path: Path) -> dict[str, _OrfCoord]:
    """Load ``<accession>.coords.tsv`` into an ``{orf_id: _OrfCoord}`` map.

    Args:
        path: Path to a coords TSV emitted by :func:`src.orf_predictor.predict_orfs`.

    Returns:
        Mapping keyed by the full ORF id (``accession|contig|orf_N``).
    """
    out: dict[str, _OrfCoord] = {}
    with path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out[row["orf_id"]] = _OrfCoord(
                contig_id=row["contig_id"],
                start=int(row["start"]),
                end=int(row["end"]),
                strand=row["strand"],
                contig_length=int(row["contig_length"]),
            )
    return out


def parse_domtblout(path: Path) -> dict[str, _BestHit]:
    """Load a HMMER ``.domtblout`` into an ``{orf_id: _BestHit}`` map.

    For each target ORF only the hit with the lowest i-Evalue is retained, so
    one ORF cannot count for two different markers in the downstream synteny rule.

    Args:
        path: Path to a ``.domtblout`` produced by :func:`src.hmm_searcher.run_search`.

    Returns:
        Mapping of ORF id to its best (lowest i-Evalue) HMM hit.
    """
    best: dict[str, _BestHit] = {}
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 23:
                continue
            orf_id = parts[0]
            hmm_name = parts[3]
            try:
                full_evalue = float(parts[6])
                full_score = float(parts[7])
                i_evalue = float(parts[12])
            except ValueError:
                continue
            prior = best.get(orf_id)
            if prior is None or i_evalue < prior.i_evalue:
                best[orf_id] = _BestHit(
                    orf_id=orf_id,
                    hmm_name=hmm_name,
                    i_evalue=i_evalue,
                    full_evalue=full_evalue,
                    score=full_score,
                )
    return best


def _best_window_indices(
    positions: list[tuple[str, int]], window_bp: int
) -> tuple[int, int, set[str]]:
    """Identify the densest fixed-width window on a contig.

    Used twice: once with :data:`SYNTENY_WINDOW_BP` for the candidate trigger,
    and once with :data:`CLUSTER_SPAN_BP` for the classification count.

    Args:
        positions: ``(hmm_name, midpoint_bp)`` tuples, sorted ascending by midpoint.
        window_bp: Window width in base pairs.

    Returns:
        ``(left, right, distinct_names)`` where ``left``/``right`` are inclusive
        indices into ``positions`` bounding the densest window and
        ``distinct_names`` is the set of HMM names found in it. If
        ``positions`` is empty, the indices are ``(0, -1)`` and the set is empty.
    """
    if not positions:
        return 0, -1, set()
    best_left = 0
    best_right = 0
    best_names: set[str] = {positions[0][0]}
    left = 0
    for right in range(len(positions)):
        while positions[right][1] - positions[left][1] > window_bp:
            left += 1
        names_in_window = {positions[i][0] for i in range(left, right + 1)}
        if len(names_in_window) > len(best_names):
            best_left = left
            best_right = right
            best_names = names_in_window
    return best_left, best_right, best_names


def _classify(n_distinct: int, edge_truncated: bool) -> Classification:
    """Apply the classification ladder defined in the Phase 2 plan.

    Args:
        n_distinct: Number of distinct FaRLiP marker HMMs on this contig.
        edge_truncated: True if any marker is within
            :data:`EDGE_PROXIMITY_BP` of either contig end.

    Returns:
        The resulting :class:`Classification` tag.
    """
    if n_distinct >= FULL_CLUSTER_MARKER_THRESHOLD:
        return Classification.FULL_CLUSTER
    if n_distinct >= MIN_MARKERS_FOR_CANDIDATE and edge_truncated:
        return Classification.EDGE_TRUNCATED
    if n_distinct >= MIN_MARKERS_FOR_CANDIDATE:
        return Classification.PARTIAL_CLUSTER
    return Classification.NONE


def classify_mag(
    mag_accession: str,
    coords_path: Path,
    domtblout_path: Path,
) -> list[ContigCall]:
    """Classify every contig in one MAG against the FaRLiP synteny rule.

    Args:
        mag_accession: MAG / assembly accession key.
        coords_path: Path to ``<accession>.coords.tsv``.
        domtblout_path: Path to ``<accession>.domtblout``.

    Returns:
        :class:`ContigCall` objects, one per **candidate or non-NONE** contig
        (contigs with no FaRLiP hits at all are omitted from the per-contig output;
        their absence is the reason their MAG ends up with ``n_candidate_contigs=0``
        in the per-MAG roll-up).
    """
    coords = parse_coords_tsv(coords_path)
    best_hits = parse_domtblout(domtblout_path)

    by_contig: dict[str, list[tuple[str, int, int, float]]] = {}
    contig_lengths: dict[str, int] = {}

    for orf_id, hit in best_hits.items():
        coord = coords.get(orf_id)
        if coord is None:
            continue
        by_contig.setdefault(coord.contig_id, []).append(
            (hit.hmm_name, coord.start, coord.end, hit.i_evalue)
        )
        contig_lengths[coord.contig_id] = coord.contig_length

    calls: list[ContigCall] = []
    for contig_id, all_hits in by_contig.items():
        all_hits.sort(key=lambda m: m[1])
        positions = [(m[0], (m[1] + m[2]) // 2) for m in all_hits]

        _, _, trigger_names = _best_window_indices(positions, SYNTENY_WINDOW_BP)
        if len(trigger_names) < MIN_MARKERS_FOR_CANDIDATE:
            continue

        left, right, cluster_names = _best_window_indices(positions, CLUSTER_SPAN_BP)
        window_markers = all_hits[left : right + 1]
        contig_length = contig_lengths[contig_id]
        edge_truncated = any(
            (m[1] <= EDGE_PROXIMITY_BP)
            or (contig_length - m[2] <= EDGE_PROXIMITY_BP)
            for m in window_markers
        )
        classification = _classify(len(cluster_names), edge_truncated)
        if classification == Classification.NONE:
            continue
        calls.append(
            ContigCall(
                mag_accession=mag_accession,
                contig_id=contig_id,
                contig_length=contig_length,
                markers=window_markers,
                n_total_hits_on_contig=len(all_hits),
                edge_truncated=edge_truncated,
                classification=classification,
            )
        )

    calls.sort(key=lambda c: (-c.n_markers, c.contig_id))
    return calls


def _pick_top_classification(calls: list[ContigCall]) -> Classification:
    """Return the highest-priority classification across a MAG's contig calls."""
    if not calls:
        return Classification.NONE
    return max(calls, key=lambda c: _CLASSIFICATION_PRIORITY[c.classification]).classification


def write_reports(
    all_calls: dict[str, list[ContigCall]],
    report_path: Path,
    per_mag_path: Path,
) -> None:
    """Serialize per-contig and per-MAG synteny reports.

    Args:
        all_calls: Mapping of MAG accession -> list of :class:`ContigCall`.
            MAGs with an empty list are still recorded in the per-MAG report (with
            ``top_classification = NONE``) so a downstream consumer can see that
            every MAG in the manifest was scanned.
        report_path: Output path for the per-contig TSV.
        per_mag_path: Output path for the per-MAG TSV.
    """
    report_path.parent.mkdir(parents=True, exist_ok=True)
    per_mag_path.parent.mkdir(parents=True, exist_ok=True)

    with report_path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(CONTIG_REPORT_COLUMNS) + "\n")
        for mag, calls in sorted(all_calls.items()):
            for call in calls:
                marker_names = ",".join(m[0] for m in call.markers)
                marker_positions = ",".join(f"{m[1]}-{m[2]}" for m in call.markers)
                fh.write(
                    "\t".join(
                        [
                            mag,
                            call.contig_id,
                            str(call.contig_length),
                            str(call.n_markers),
                            str(call.n_total_hits_on_contig),
                            marker_names,
                            marker_positions,
                            f"{call.min_evalue:.2e}",
                            "true" if call.edge_truncated else "false",
                            call.classification.value,
                        ]
                    )
                    + "\n"
                )

    with per_mag_path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(PER_MAG_REPORT_COLUMNS) + "\n")
        for mag, calls in sorted(all_calls.items()):
            n_candidates = len(calls)
            total_markers = sum(c.n_markers for c in calls)
            distinct_markers = len({m[0] for c in calls for m in c.markers})
            has_edge = any(c.edge_truncated for c in calls)
            top = _pick_top_classification(calls)
            fh.write(
                "\t".join(
                    [
                        mag,
                        str(n_candidates),
                        str(total_markers),
                        str(distinct_markers),
                        "true" if has_edge else "false",
                        top.value,
                    ]
                )
                + "\n"
            )
