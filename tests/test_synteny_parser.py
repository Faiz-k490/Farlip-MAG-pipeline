"""Tests for ``src.synteny_parser``.

The synteny rules and edge-truncation logic are pure-Python and ideal for unit
testing. We build small synthetic domtblout + coords.tsv files in tmp_path and
exercise every classification branch:

* No markers -> no candidate.
* 2 markers in 10 kb -> below the trigger -> NONE.
* 3+ markers tightly co-located -> PARTIAL_CLUSTER.
* 3+ markers near contig end -> EDGE_TRUNCATED.
* 15+ markers in cluster span -> FULL_CLUSTER.
* Best-hit-per-ORF tiebreak: same ORF appearing in two HMMs counts once.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.synteny_parser import (
    CLUSTER_SPAN_BP,
    CONTIG_REPORT_COLUMNS,
    EDGE_PROXIMITY_BP,
    FULL_CLUSTER_MARKER_THRESHOLD,
    MIN_MARKERS_FOR_CANDIDATE,
    PER_MAG_REPORT_COLUMNS,
    SYNTENY_WINDOW_BP,
    Classification,
    classify_mag,
    parse_coords_tsv,
    parse_domtblout,
    write_reports,
)


def _write_coords(path: Path, rows: list[tuple[str, str, int, int, int]]) -> None:
    """Write a minimal coords.tsv fixture.

    Args:
        path: Output path.
        rows: ``(contig_id, orf_id, start, end, contig_length)`` tuples. Strand
            is fixed to ``+`` in this fixture.
    """
    with path.open("w", encoding="utf-8") as fh:
        fh.write("contig_id\torf_id\tstart\tend\tstrand\tpartial_left\tpartial_right\tcontig_length\n")
        for contig_id, orf_id, start, end, contig_length in rows:
            fh.write(
                f"{contig_id}\t{orf_id}\t{start}\t{end}\t+\t0\t0\t{contig_length}\n"
            )


def _write_domtblout(path: Path, rows: list[tuple[str, str, float]]) -> None:
    """Write a minimal HMMER-format domtblout fixture.

    Args:
        path: Output path.
        rows: ``(orf_id, hmm_name, i_evalue)`` tuples. All other fields are
            placeholders sufficient to satisfy the parser.
    """
    with path.open("w", encoding="utf-8") as fh:
        for orf_id, hmm_name, ievalue in rows:
            # 23+ whitespace-delimited columns; only columns 1 (target), 4 (query)
            # 7 (full-Evalue), 8 (full-score), 13 (i-Evalue) are read by the parser.
            fh.write(
                f"{orf_id} - 100 {hmm_name} - 100 {ievalue:g} 100.0 0.0 "
                f"1 1 {ievalue:g} {ievalue:g} 100.0 0.0 1 100 1 100 1 100 0.95 -\n"
            )


def test_constants_are_consistent():
    """Constants must form a sensible ladder."""
    assert MIN_MARKERS_FOR_CANDIDATE >= 2
    assert FULL_CLUSTER_MARKER_THRESHOLD > MIN_MARKERS_FOR_CANDIDATE
    assert CLUSTER_SPAN_BP >= SYNTENY_WINDOW_BP
    assert EDGE_PROXIMITY_BP > 0


def test_no_markers_returns_no_calls(tmp_path: Path):
    """An empty domtblout yields zero candidate contigs."""
    coords = tmp_path / "coords.tsv"
    dom = tmp_path / "dom.tsv"
    _write_coords(coords, [("ctg1", "ORF1", 100, 1000, 100_000)])
    _write_domtblout(dom, [])
    calls = classify_mag("DEMO", coords, dom)
    assert calls == []


def test_two_markers_below_trigger(tmp_path: Path):
    """Two markers in 10 kb is below MIN_MARKERS_FOR_CANDIDATE -> no candidate."""
    coords = tmp_path / "coords.tsv"
    dom = tmp_path / "dom.tsv"
    _write_coords(
        coords,
        [
            ("ctg1", "ORF1", 1000, 2000, 100_000),
            ("ctg1", "ORF2", 3000, 4000, 100_000),
        ],
    )
    _write_domtblout(
        dom,
        [
            ("ORF1", "rfpA", 1e-30),
            ("ORF2", "psbA2", 1e-30),
        ],
    )
    calls = classify_mag("DEMO", coords, dom)
    assert calls == []


def test_three_markers_in_window_yields_partial(tmp_path: Path):
    """Three distinct markers inside 10 kb -> PARTIAL_CLUSTER."""
    coords = tmp_path / "coords.tsv"
    dom = tmp_path / "dom.tsv"
    _write_coords(
        coords,
        [
            ("ctg1", "ORF1", 50_000, 51_000, 200_000),
            ("ctg1", "ORF2", 52_000, 53_000, 200_000),
            ("ctg1", "ORF3", 54_000, 55_000, 200_000),
        ],
    )
    _write_domtblout(
        dom,
        [
            ("ORF1", "rfpA", 1e-30),
            ("ORF2", "psbA2", 1e-30),
            ("ORF3", "apcE2", 1e-30),
        ],
    )
    calls = classify_mag("DEMO", coords, dom)
    assert len(calls) == 1
    assert calls[0].classification == Classification.PARTIAL_CLUSTER
    assert calls[0].n_markers == 3
    assert calls[0].edge_truncated is False


def test_edge_proximity_flips_to_edge_truncated(tmp_path: Path):
    """Three markers within EDGE_PROXIMITY_BP of contig end -> EDGE_TRUNCATED."""
    coords = tmp_path / "coords.tsv"
    dom = tmp_path / "dom.tsv"
    contig_length = 20_000
    _write_coords(
        coords,
        [
            ("ctg1", "ORF1", 500, 1500, contig_length),
            ("ctg1", "ORF2", 2000, 3000, contig_length),
            ("ctg1", "ORF3", 4000, 5000, contig_length),
        ],
    )
    _write_domtblout(
        dom,
        [
            ("ORF1", "rfpA", 1e-30),
            ("ORF2", "psbA2", 1e-30),
            ("ORF3", "apcE2", 1e-30),
        ],
    )
    calls = classify_mag("DEMO", coords, dom)
    assert len(calls) == 1
    assert calls[0].edge_truncated is True
    assert calls[0].classification == Classification.EDGE_TRUNCATED


def test_full_cluster_threshold(tmp_path: Path):
    """FULL_CLUSTER_MARKER_THRESHOLD distinct markers in CLUSTER_SPAN_BP -> FULL_CLUSTER."""
    coords_rows: list[tuple[str, str, int, int, int]] = []
    dom_rows: list[tuple[str, str, float]] = []
    for i in range(FULL_CLUSTER_MARKER_THRESHOLD):
        orf_id = f"ORF{i:03d}"
        start = 50_000 + i * 1500
        end = start + 1000
        coords_rows.append(("ctg1", orf_id, start, end, 500_000))
        dom_rows.append((orf_id, f"marker_{i:02d}", 1e-30))
    coords = tmp_path / "coords.tsv"
    dom = tmp_path / "dom.tsv"
    _write_coords(coords, coords_rows)
    _write_domtblout(dom, dom_rows)
    calls = classify_mag("DEMO", coords, dom)
    assert len(calls) == 1
    assert calls[0].classification == Classification.FULL_CLUSTER
    assert calls[0].n_markers >= FULL_CLUSTER_MARKER_THRESHOLD


def test_best_hit_per_orf_dedup(tmp_path: Path):
    """One ORF appearing under two HMMs must count for only the better hit."""
    coords = tmp_path / "coords.tsv"
    dom = tmp_path / "dom.tsv"
    _write_coords(
        coords,
        [
            ("ctg1", "ORF1", 50_000, 51_000, 200_000),
            ("ctg1", "ORF2", 52_000, 53_000, 200_000),
        ],
    )
    _write_domtblout(
        dom,
        [
            ("ORF1", "rfpA", 1e-10),
            ("ORF1", "rfpB", 1e-30),
            ("ORF2", "psbA2", 1e-30),
        ],
    )
    calls = classify_mag("DEMO", coords, dom)
    assert calls == []


def test_parse_coords_tsv_round_trip(tmp_path: Path):
    """``parse_coords_tsv`` returns one entry per row, keyed by ORF id."""
    coords = tmp_path / "coords.tsv"
    _write_coords(coords, [("ctg1", "ORF1", 100, 200, 1000)])
    out = parse_coords_tsv(coords)
    assert "ORF1" in out
    assert out["ORF1"].contig_length == 1000


def test_parse_domtblout_keeps_best_hit(tmp_path: Path):
    """``parse_domtblout`` retains only the lowest-i-Evalue hit per ORF."""
    dom = tmp_path / "dom.tsv"
    _write_domtblout(
        dom,
        [
            ("ORF1", "rfpA", 1e-5),
            ("ORF1", "rfpB", 1e-30),
        ],
    )
    out = parse_domtblout(dom)
    assert out["ORF1"].hmm_name == "rfpB"


def test_write_reports_headers(tmp_path: Path):
    """The TSV headers must match the declared column order."""
    report_path = tmp_path / "report.tsv"
    per_mag_path = tmp_path / "per_mag.tsv"
    write_reports({}, report_path, per_mag_path)
    assert report_path.read_text(encoding="utf-8").splitlines()[0] == "\t".join(
        CONTIG_REPORT_COLUMNS
    )
    assert per_mag_path.read_text(encoding="utf-8").splitlines()[0] == "\t".join(
        PER_MAG_REPORT_COLUMNS
    )
