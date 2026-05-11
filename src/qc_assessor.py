"""In-process genome-quality assessor for Phase 3.

Scans every per-MAG protein FASTA (``data/orfs/<accession>.faa`` written by
:mod:`src.orf_predictor`) against the bacterial single-copy ortholog (SCO)
marker HMM database (``data/qc_markers/<label>.hmm`` written by
:mod:`src.qc_markers_fetcher`) using ``pyhmmer.hmmsearch``, and computes:

* **completeness**  = ``100 * (# markers with >= 1 reported hit) / total_markers``
* **contamination** = ``100 * (# markers with >= 2 reported hits) / total_markers``
* **n50_contig**    = standard N50 over the distinct contig lengths recorded in
                       ``data/orfs/<accession>.coords.tsv``.
* **genome_size_bp** = sum of the distinct contig lengths.
* **qc_quality**    = MIMAG ladder (HIGH / MEDIUM / LOW) per Bowers et al. 2017
                       (Nat Biotechnol 35:725) using the constants below.

Idempotency mirrors :mod:`src.orf_predictor`: each per-MAG ``.qc.tsv`` carries
a sidecar stamp file recording the SHA-256 of both inputs (``.faa`` and the
marker DB). A re-run with both hashes unchanged short-circuits.

Constants are module-level by convention - no biological threshold is buried
inside the function body.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from pathlib import Path

import pyhmmer

from .orf_predictor import _sha256_file


QC_HMM_EVALUE_THRESHOLD: float = 1e-10
"""Sequence-level E-value floor used for marker presence/duplication calls.

This is a backstop. The preferred discriminator is the per-HMM bit-score
cutoff sidecar shipped by BUSCO (``scores_cutoff``) - flat E-value alone
cannot separate the true single-copy ortholog from divergent paralogs that
share a conserved domain. When a cutoff sidecar is available, hits below the
HMM-specific bit score are dropped before the completeness / contamination
tally; the E-value still applies as a coarse pre-filter."""

QC_HMM_DOMAIN_EVALUE_THRESHOLD: float = 1e-10
"""Per-domain E-value cutoff. Same rationale as :data:`QC_HMM_EVALUE_THRESHOLD`."""

QC_HMM_CPU_COUNT: int = max(1, (os.cpu_count() or 2) - 1)
"""Number of worker threads handed to ``pyhmmer.hmmsearch`` per MAG. Matches the
Phase 2 default so the QC sweep does not contend with whatever else is running
on the dev box."""

MIMAG_HIGH_QUALITY_COMPLETENESS: float = 90.0
"""MIMAG HIGH-quality threshold (Bowers et al. 2017): completeness >= 90%."""

MIMAG_HIGH_QUALITY_CONTAMINATION: float = 5.0
"""MIMAG HIGH-quality threshold (Bowers et al. 2017): contamination <= 5%."""

MIMAG_MEDIUM_QUALITY_COMPLETENESS: float = 50.0
"""MIMAG MEDIUM-quality threshold (Bowers et al. 2017): completeness >= 50%."""

MIMAG_MEDIUM_QUALITY_CONTAMINATION: float = 10.0
"""MIMAG MEDIUM-quality threshold (Bowers et al. 2017): contamination < 10%."""

QC_TSV_COLUMNS: tuple[str, ...] = (
    "accession",
    "qc_source",
    "marker_set_label",
    "n_markers_total",
    "n_markers_found",
    "n_markers_duplicated",
    "completeness",
    "contamination",
    "n50_contig",
    "genome_size_bp",
    "n_contigs",
    "longest_contig_bp",
    "qc_quality",
)
"""Schema of every per-MAG ``data/qc_assessor/<accession>.qc.tsv`` row.

Stable - the reporter joins on these column names. Adding new columns must
happen at the right of the tuple."""


QC_QUALITY_HIGH: str = "HIGH"
"""MIMAG high-quality bin (>= :data:`MIMAG_HIGH_QUALITY_COMPLETENESS` and
<= :data:`MIMAG_HIGH_QUALITY_CONTAMINATION`)."""

QC_QUALITY_MEDIUM: str = "MEDIUM"
"""MIMAG medium-quality bin."""

QC_QUALITY_LOW: str = "LOW"
"""Anything below medium quality (still recorded so the row count matches the
manifest)."""

QC_SOURCE_PYHMMER: str = "pyhmmer"
"""Default ``qc_source`` tag emitted by this module. ``checkm2`` rows come from
:mod:`src.checkm2_ingestor` and are preferred by the reporter when both exist
for the same accession."""


@dataclass(frozen=True)
class QcAssessment:
    """Outcome of one :func:`assess_assembly` invocation.

    Attributes:
        accession: MAG / assembly accession key.
        qc_tsv_path: Path to the written ``.qc.tsv``.
        stamp_path: Path to the sidecar stamp file.
        n_markers_total: Total markers in the loaded HMM database.
        n_markers_found: Markers with at least one reported hit.
        n_markers_duplicated: Markers with two or more reported hits.
        completeness: Percent of markers found.
        contamination: Percent of markers duplicated.
        n50_contig: N50 over distinct contig lengths.
        genome_size_bp: Sum of distinct contig lengths.
        n_contigs: Number of distinct contigs.
        longest_contig_bp: Length of the largest contig.
        qc_quality: MIMAG bin (HIGH / MEDIUM / LOW).
        skipped: True when the run short-circuited via the stamp.
    """

    accession: str
    qc_tsv_path: Path
    stamp_path: Path
    n_markers_total: int
    n_markers_found: int
    n_markers_duplicated: int
    completeness: float
    contamination: float
    n50_contig: int
    genome_size_bp: int
    n_contigs: int
    longest_contig_bp: int
    qc_quality: str
    skipped: bool


def _read_stamp(stamp_path: Path) -> dict[str, str]:
    """Parse a stamp file into a ``{key: value}`` dict (or ``{}`` if missing)."""
    if not stamp_path.exists():
        return {}
    out: dict[str, str] = {}
    for line in stamp_path.read_text(encoding="utf-8").splitlines():
        if "\t" not in line:
            continue
        k, v = line.split("\t", 1)
        out[k] = v
    return out


def _write_stamp(
    stamp_path: Path,
    *,
    accession: str,
    faa_sha256: str,
    hmm_db_sha256: str,
    marker_set_label: str,
    n_markers_total: int,
    n_markers_found: int,
    scores_cutoff_sha256: str = "",
) -> None:
    """Write a stamp file recording what inputs produced this ``.qc.tsv``.

    Args:
        stamp_path: Output path.
        accession: MAG accession.
        faa_sha256: SHA-256 of the input ``.faa`` at run time.
        hmm_db_sha256: SHA-256 of the marker HMM database at run time.
        marker_set_label: Human-readable label of the marker set used.
        n_markers_total: Total markers in the database.
        n_markers_found: Markers detected in this MAG.
        scores_cutoff_sha256: SHA-256 of the per-HMM bit-score cutoff sidecar
            (empty string if no cutoff sidecar was applied).
    """
    stamp_path.parent.mkdir(parents=True, exist_ok=True)
    payload = [
        f"accession\t{accession}",
        f"faa_sha256\t{faa_sha256}",
        f"hmm_db_sha256\t{hmm_db_sha256}",
        f"scores_cutoff_sha256\t{scores_cutoff_sha256}",
        f"marker_set_label\t{marker_set_label}",
        f"n_markers_total\t{n_markers_total}",
        f"n_markers_found\t{n_markers_found}",
    ]
    stamp_path.write_text("\n".join(payload) + "\n", encoding="utf-8")


def _read_distinct_contig_lengths(coords_path: Path) -> list[int]:
    """Extract one length per distinct contig from a coords TSV.

    Args:
        coords_path: Path to ``<accession>.coords.tsv``.

    Returns:
        List of contig lengths in coords-file order. An empty list is returned
        if the coords file is missing or empty so the assessor can still emit
        a ``n50_contig=0`` row instead of crashing the whole sweep.
    """
    if not coords_path.exists():
        return []
    seen: dict[str, int] = {}
    with coords_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cid = row.get("contig_id", "")
            if not cid or cid in seen:
                continue
            try:
                seen[cid] = int(row.get("contig_length", "0") or 0)
            except ValueError:
                continue
    return [length for length in seen.values() if length > 0]


def _compute_n50(lengths: list[int]) -> int:
    """Compute the N50 statistic on a list of contig lengths.

    N50 is the contig length such that contigs of at least that length together
    account for >= 50% of the assembly size.

    Args:
        lengths: Contig lengths (any order).

    Returns:
        The N50, or 0 if ``lengths`` is empty.
    """
    if not lengths:
        return 0
    total = sum(lengths)
    half = total / 2.0
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= half:
            return length
    return lengths[-1]


def _classify_quality(completeness: float, contamination: float) -> str:
    """Apply the MIMAG ladder.

    Args:
        completeness: Percent (0-100).
        contamination: Percent (0-100).

    Returns:
        :data:`QC_QUALITY_HIGH`, :data:`QC_QUALITY_MEDIUM`, or
        :data:`QC_QUALITY_LOW`.
    """
    if (
        completeness >= MIMAG_HIGH_QUALITY_COMPLETENESS
        and contamination <= MIMAG_HIGH_QUALITY_CONTAMINATION
    ):
        return QC_QUALITY_HIGH
    if (
        completeness >= MIMAG_MEDIUM_QUALITY_COMPLETENESS
        and contamination < MIMAG_MEDIUM_QUALITY_CONTAMINATION
    ):
        return QC_QUALITY_MEDIUM
    return QC_QUALITY_LOW


def load_score_cutoffs(path: Path | None) -> dict[str, float]:
    """Parse a BUSCO-style ``scores_cutoff`` file into ``{hmm_name: bit_score}``.

    Args:
        path: Path to a TSV with one ``<hmm_name>\\t<bit_score>`` per line, or
            ``None`` to opt out.

    Returns:
        Mapping. Empty dict when ``path`` is ``None`` or absent. Lines that
        don't parse cleanly are silently skipped.
    """
    if path is None or not path.exists():
        return {}
    out: dict[str, float] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            out[parts[0]] = float(parts[1])
        except ValueError:
            continue
    return out


def _count_marker_hits(
    hmms: list[pyhmmer.plan7.HMM],
    sequences: list,
    *,
    cpus: int,
    evalue: float,
    domain_evalue: float,
    score_cutoffs: dict[str, float] | None,
) -> tuple[int, int]:
    """Run hmmsearch and tally how many markers are present / duplicated.

    When ``score_cutoffs`` is non-empty, hits whose full-sequence bit score
    falls below the per-HMM cutoff are dropped from the tally - the canonical
    BUSCO discrimination rule.

    Args:
        hmms: The marker HMM list (one entry per single-copy ortholog).
        sequences: Digital-mode sequences loaded from the per-MAG ``.faa``.
        cpus: Worker thread count.
        evalue: Sequence-level E-value cutoff.
        domain_evalue: Per-domain E-value cutoff.
        score_cutoffs: Optional ``{hmm_name: min_bit_score}`` mapping. Pass
            ``None`` or an empty dict to disable per-HMM gating (E-value-only
            mode).

    Returns:
        Tuple of ``(n_markers_found, n_markers_duplicated)``. ``found`` counts
        markers with any retained hit; ``duplicated`` counts markers with two
        or more retained hits.
    """
    cutoffs = score_cutoffs or {}
    n_found = 0
    n_dup = 0
    for top_hits in pyhmmer.hmmsearch(
        hmms,
        sequences,
        cpus=cpus,
        E=evalue,
        domE=domain_evalue,
    ):
        query_name = top_hits.query.name
        hmm_name = (
            query_name.decode("ascii", errors="replace")
            if isinstance(query_name, (bytes, bytearray))
            else str(query_name)
        )
        cutoff = cutoffs.get(hmm_name)
        if cutoff is None:
            retained = sum(1 for h in top_hits if h.reported)
        else:
            retained = sum(
                1 for h in top_hits if h.reported and h.score >= cutoff
            )
        if retained >= 1:
            n_found += 1
        if retained >= 2:
            n_dup += 1
    return n_found, n_dup


def assess_assembly(
    accession: str,
    faa_path: Path,
    coords_path: Path,
    hmm_db_path: Path,
    output_dir: Path,
    *,
    marker_set_label: str,
    overwrite: bool = False,
    cpus: int = QC_HMM_CPU_COUNT,
    evalue: float = QC_HMM_EVALUE_THRESHOLD,
    domain_evalue: float = QC_HMM_DOMAIN_EVALUE_THRESHOLD,
    scores_cutoff_path: Path | None = None,
) -> QcAssessment:
    """Score one MAG for completeness + contamination + N50 + genome size.

    Idempotency:
        Re-running with the same ``faa_path``, ``hmm_db_path`` and (if
        provided) ``scores_cutoff_path`` short-circuits when the stamp
        records all hashes unchanged. Set ``overwrite=True`` to force
        re-assessment.

    Args:
        accession: MAG / assembly accession key.
        faa_path: Path to ``data/orfs/<accession>.faa``.
        coords_path: Path to ``data/orfs/<accession>.coords.tsv``.
        hmm_db_path: Path to the concatenated marker HMM database.
        output_dir: Directory under which ``<accession>.qc.tsv`` and
            ``.qc.stamp`` are written.
        marker_set_label: Human-readable label written into the TSV (e.g.
            ``"BUSCO_bacteria_odb10"``).
        overwrite: If True, ignore any existing stamp.
        cpus: Worker threads handed to ``pyhmmer.hmmsearch``.
        evalue: Sequence-level E-value cutoff (coarse pre-filter when
            per-HMM bit-score cutoffs are also supplied).
        domain_evalue: Per-domain E-value cutoff.
        scores_cutoff_path: Optional path to a BUSCO-style
            ``scores_cutoff`` TSV. When supplied, hits below the per-HMM
            bit score are dropped before the tally (recommended).

    Returns:
        A :class:`QcAssessment` describing what was written and computed.

    Raises:
        FileNotFoundError: If any required input file is missing.
        RuntimeError: If the marker HMM database is empty.
    """
    if not faa_path.exists():
        raise FileNotFoundError(f"protein FASTA not found: {faa_path}")
    if not hmm_db_path.exists():
        raise FileNotFoundError(f"marker HMM database not found: {hmm_db_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    qc_tsv_path = output_dir / f"{accession}.qc.tsv"
    stamp_path = output_dir / f"{accession}.qc.stamp"

    faa_sha256 = _sha256_file(faa_path)
    hmm_db_sha256 = _sha256_file(hmm_db_path)
    scores_cutoff_sha256 = (
        _sha256_file(scores_cutoff_path)
        if scores_cutoff_path is not None and scores_cutoff_path.exists()
        else ""
    )

    if not overwrite and qc_tsv_path.exists():
        prior = _read_stamp(stamp_path)
        if (
            prior.get("faa_sha256") == faa_sha256
            and prior.get("hmm_db_sha256") == hmm_db_sha256
            and prior.get("scores_cutoff_sha256", "") == scores_cutoff_sha256
        ):
            row = _read_qc_row(qc_tsv_path)
            if row is not None:
                return QcAssessment(
                    accession=accession,
                    qc_tsv_path=qc_tsv_path,
                    stamp_path=stamp_path,
                    n_markers_total=int(row.get("n_markers_total", 0) or 0),
                    n_markers_found=int(row.get("n_markers_found", 0) or 0),
                    n_markers_duplicated=int(row.get("n_markers_duplicated", 0) or 0),
                    completeness=float(row.get("completeness", 0.0) or 0.0),
                    contamination=float(row.get("contamination", 0.0) or 0.0),
                    n50_contig=int(row.get("n50_contig", 0) or 0),
                    genome_size_bp=int(row.get("genome_size_bp", 0) or 0),
                    n_contigs=int(row.get("n_contigs", 0) or 0),
                    longest_contig_bp=int(row.get("longest_contig_bp", 0) or 0),
                    qc_quality=str(row.get("qc_quality", QC_QUALITY_LOW) or QC_QUALITY_LOW),
                    skipped=True,
                )

    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.plan7.HMMFile(hmm_db_path) as hmm_file:
        hmms = list(hmm_file)
    if not hmms:
        raise RuntimeError(f"no HMMs loaded from marker database {hmm_db_path}")
    n_markers_total = len(hmms)

    with pyhmmer.easel.SequenceFile(faa_path, digital=True, alphabet=alphabet) as seq_file:
        sequences = list(seq_file)

    score_cutoffs = load_score_cutoffs(scores_cutoff_path)
    n_found, n_dup = _count_marker_hits(
        hmms,
        sequences,
        cpus=cpus,
        evalue=evalue,
        domain_evalue=domain_evalue,
        score_cutoffs=score_cutoffs,
    )
    completeness = 100.0 * n_found / n_markers_total if n_markers_total else 0.0
    contamination = 100.0 * n_dup / n_markers_total if n_markers_total else 0.0

    contig_lengths = _read_distinct_contig_lengths(coords_path)
    n50 = _compute_n50(contig_lengths)
    genome_size = sum(contig_lengths)
    longest = max(contig_lengths) if contig_lengths else 0
    quality = _classify_quality(completeness, contamination)

    _write_qc_row(
        qc_tsv_path,
        {
            "accession": accession,
            "qc_source": QC_SOURCE_PYHMMER,
            "marker_set_label": marker_set_label,
            "n_markers_total": str(n_markers_total),
            "n_markers_found": str(n_found),
            "n_markers_duplicated": str(n_dup),
            "completeness": f"{completeness:.2f}",
            "contamination": f"{contamination:.2f}",
            "n50_contig": str(n50),
            "genome_size_bp": str(genome_size),
            "n_contigs": str(len(contig_lengths)),
            "longest_contig_bp": str(longest),
            "qc_quality": quality,
        },
    )

    _write_stamp(
        stamp_path,
        accession=accession,
        faa_sha256=faa_sha256,
        hmm_db_sha256=hmm_db_sha256,
        marker_set_label=marker_set_label,
        n_markers_total=n_markers_total,
        n_markers_found=n_found,
        scores_cutoff_sha256=scores_cutoff_sha256,
    )

    return QcAssessment(
        accession=accession,
        qc_tsv_path=qc_tsv_path,
        stamp_path=stamp_path,
        n_markers_total=n_markers_total,
        n_markers_found=n_found,
        n_markers_duplicated=n_dup,
        completeness=completeness,
        contamination=contamination,
        n50_contig=n50,
        genome_size_bp=genome_size,
        n_contigs=len(contig_lengths),
        longest_contig_bp=longest,
        qc_quality=quality,
        skipped=False,
    )


def _read_qc_row(path: Path) -> dict[str, str] | None:
    """Read the single data row out of a per-MAG ``.qc.tsv``.

    Args:
        path: Path to a ``<accession>.qc.tsv`` file.

    Returns:
        The first data row as a dict, or ``None`` if the file is missing or
        contains no data row.
    """
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            return dict(row)
    return None


def _write_qc_row(path: Path, row: dict[str, str]) -> None:
    """Write one ``.qc.tsv`` row using the canonical :data:`QC_TSV_COLUMNS` order.

    Args:
        path: Output path.
        row: Mapping containing every key in :data:`QC_TSV_COLUMNS`.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=QC_TSV_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerow({col: row.get(col, "") for col in QC_TSV_COLUMNS})


__all__ = [
    "MIMAG_HIGH_QUALITY_COMPLETENESS",
    "MIMAG_HIGH_QUALITY_CONTAMINATION",
    "MIMAG_MEDIUM_QUALITY_COMPLETENESS",
    "MIMAG_MEDIUM_QUALITY_CONTAMINATION",
    "QC_HMM_CPU_COUNT",
    "QC_HMM_DOMAIN_EVALUE_THRESHOLD",
    "QC_HMM_EVALUE_THRESHOLD",
    "QC_QUALITY_HIGH",
    "QC_QUALITY_LOW",
    "QC_QUALITY_MEDIUM",
    "QC_SOURCE_PYHMMER",
    "QC_TSV_COLUMNS",
    "QcAssessment",
    "assess_assembly",
]
