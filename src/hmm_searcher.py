"""Run pyhmmer's ``hmmsearch`` across every MAG protein set.

This module wraps :func:`pyhmmer.hmmsearch` and writes a standard HMMER domain
table (``.domtblout``) per MAG so downstream tools can consume the output without
any pipeline-specific format. The synteny parser then joins these tables with the
per-MAG ``coords.tsv`` files emitted by :mod:`src.orf_predictor`.

E-value thresholds are declared as module-level constants so they show up in code
review whenever the team revisits noise tolerance. Compute parallelism is also a
module-level constant, defaulted to ``cpu_count() - 1`` so the host stays
responsive during long searches.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import pyhmmer


HMM_EVALUE_THRESHOLD: float = 1e-5
"""Sequence-level E-value cutoff applied at search time.

Hits above this threshold are not even reported. Synteny parsing applies the same
cutoff defensively, so changing this value here will not silently relax the
downstream filter."""

HMM_DOMAIN_EVALUE_THRESHOLD: float = 1e-5
"""Per-domain E-value cutoff applied at search time."""

HMM_BIT_SCORE_FLOOR: float | None = None
"""Optional bit-score floor. ``None`` means rely on the E-value cutoffs alone.

Set this to a positive float (e.g. ``25.0``) only if specific divergent paralogs
slip past the E-value cutoff on the strength of long-but-noisy alignments."""

HMM_CPU_COUNT: int = max(1, (os.cpu_count() or 2) - 1)
"""Number of worker threads pyhmmer is allowed to use per search call.

``cpu_count() - 1`` keeps the machine usable during a 38-MAG sweep. A value of 0
would let pyhmmer pick its own default (all cores); we make the choice explicit."""


@dataclass(frozen=True)
class HmmSearchResult:
    """Summary of one ``run_search`` invocation."""

    accession: str
    faa_path: Path
    domtblout_path: Path
    n_hits_reported: int
    n_hmms_with_hits: int
    skipped: bool


def run_search(
    accession: str,
    faa_path: Path,
    hmm_path: Path,
    output_dir: Path,
    *,
    overwrite: bool = False,
    cpus: int = HMM_CPU_COUNT,
    evalue_threshold: float = HMM_EVALUE_THRESHOLD,
    domain_evalue_threshold: float = HMM_DOMAIN_EVALUE_THRESHOLD,
) -> HmmSearchResult:
    """Search every HMM in ``hmm_path`` against the proteins in ``faa_path``.

    Args:
        accession: MAG / assembly accession key. Drives the output filename.
        faa_path: Path to ``<accession>.faa`` (output of :func:`predict_orfs`).
        hmm_path: Path to the concatenated HMM file (output of
            :func:`build_profiles`).
        output_dir: Directory under which ``<accession>.domtblout`` is written.
        overwrite: If True, re-run even if ``<accession>.domtblout`` exists.
        cpus: Worker thread count handed to pyhmmer. Defaults to
            :data:`HMM_CPU_COUNT`.
        evalue_threshold: Sequence-level E-value cutoff. Defaults to
            :data:`HMM_EVALUE_THRESHOLD`.
        domain_evalue_threshold: Per-domain E-value cutoff. Defaults to
            :data:`HMM_DOMAIN_EVALUE_THRESHOLD`.

    Returns:
        An :class:`HmmSearchResult` describing the run.

    Raises:
        FileNotFoundError: If either input file does not exist.
    """
    if not faa_path.exists():
        raise FileNotFoundError(f"protein FASTA not found: {faa_path}")
    if not hmm_path.exists():
        raise FileNotFoundError(f"HMM file not found: {hmm_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    domtblout_path = output_dir / f"{accession}.domtblout"

    if domtblout_path.exists() and not overwrite:
        n_hits, n_hmms = _count_domtblout_hits(domtblout_path)
        return HmmSearchResult(
            accession=accession,
            faa_path=faa_path,
            domtblout_path=domtblout_path,
            n_hits_reported=n_hits,
            n_hmms_with_hits=n_hmms,
            skipped=True,
        )

    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        hmms = list(hmm_file)

    if not hmms:
        raise RuntimeError(f"no HMMs loaded from {hmm_path}")

    with pyhmmer.easel.SequenceFile(faa_path, digital=True, alphabet=alphabet) as seq_file:
        sequences = list(seq_file)

    n_hits = 0
    n_hmms_with_hits = 0
    with domtblout_path.open("wb") as fh:
        for top_hits in pyhmmer.hmmsearch(
            hmms,
            sequences,
            cpus=cpus,
            E=evalue_threshold,
            domE=domain_evalue_threshold,
        ):
            top_hits.write(fh, format="domains", header=False)
            reported = sum(1 for h in top_hits if h.reported)
            n_hits += reported
            if reported:
                n_hmms_with_hits += 1

    return HmmSearchResult(
        accession=accession,
        faa_path=faa_path,
        domtblout_path=domtblout_path,
        n_hits_reported=n_hits,
        n_hmms_with_hits=n_hmms_with_hits,
        skipped=False,
    )


def _count_domtblout_hits(path: Path) -> tuple[int, int]:
    """Count non-comment lines and unique HMM names in an existing ``.domtblout``.

    Used to populate :class:`HmmSearchResult` when a re-run short-circuits.

    Args:
        path: Path to a HMMER domain-table file.

    Returns:
        Tuple of ``(total_hit_rows, unique_hmm_names_with_hits)``.
    """
    total = 0
    hmms: set[str] = set()
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            total += 1
            parts = line.split()
            if len(parts) >= 4:
                hmms.add(parts[3])
    return total, len(hmms)
