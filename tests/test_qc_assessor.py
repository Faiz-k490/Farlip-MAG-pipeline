"""Tests for ``src.qc_assessor``.

We build tiny single-sequence HMMs in-memory with pyhmmer and search them
against the same sequences (and against a contaminated set with extra copies)
to verify the completeness / contamination calculation. The tests also cover
the N50 / genome-size join, MIMAG ladder thresholds, and idempotency.
"""

from __future__ import annotations

from pathlib import Path

import pyhmmer
import pytest

from src.qc_assessor import (
    MIMAG_HIGH_QUALITY_COMPLETENESS,
    MIMAG_HIGH_QUALITY_CONTAMINATION,
    MIMAG_MEDIUM_QUALITY_COMPLETENESS,
    QC_HMM_CPU_COUNT,
    QC_HMM_DOMAIN_EVALUE_THRESHOLD,
    QC_HMM_EVALUE_THRESHOLD,
    QC_QUALITY_HIGH,
    QC_QUALITY_LOW,
    QC_QUALITY_MEDIUM,
    QC_TSV_COLUMNS,
    _classify_quality,
    _compute_n50,
    _read_distinct_contig_lengths,
    assess_assembly,
    load_score_cutoffs,
)


# Three distinct ~80-residue protein sequences. They share no significant
# similarity, so each HMM will only match its own protein under HMMER's
# default thresholds.
_MARKER_SEQUENCES: dict[str, str] = {
    "M1": (
        "MAKQLTPETLNRINAALLDPNVKAVEYDGDQVQVRYDLQGKEELAGRGAAGTVMSEHESVDLRTKLP"
        "VQAGAEELASRWPADLAGTGRLELQGKE"
    ),
    "M2": (
        "MFRTILVPRTGNTAEDLILRDGGSAVLQATGEVLAYRSDDFGRTGLISTNAEHGAVDPSGLATKVAK"
        "GTRDTYLNEVDSAGTRPSDFGSEAGFNV"
    ),
    "M3": (
        "MGSIDETKAKLRSAIDEQGYDIIDLEKAGGAVAVKINGREYRAEKGSDDLLPYQGVRVNTPALEAEI"
        "QSRARQAVAFAEGTAYGTGFRSDLPGLA"
    ),
}


def _build_marker_hmm_db(path: Path) -> int:
    """Build a tiny single-sequence HMM database for the test markers.

    Args:
        path: Output path for the concatenated HMM file.

    Returns:
        Number of HMMs written.
    """
    alphabet = pyhmmer.easel.Alphabet.amino()
    background = pyhmmer.plan7.Background(alphabet)
    builder = pyhmmer.plan7.Builder(alphabet)
    hmms = []
    for name, seq in _MARKER_SEQUENCES.items():
        digital = pyhmmer.easel.TextSequence(
            name=name.encode("ascii"),
            sequence=seq,
        ).digitize(alphabet)
        hmm, _profile, _optprofile = builder.build(digital, background)
        hmm.name = name.encode("ascii")
        hmms.append(hmm)
    with path.open("wb") as fh:
        for h in hmms:
            h.write(fh)
    return len(hmms)


def _write_faa(path: Path, records: list[tuple[str, str]]) -> None:
    """Write a list of ``(orf_id, sequence)`` pairs as a FASTA file."""
    with path.open("w", encoding="utf-8") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


def _write_coords(path: Path, rows: list[dict[str, str]]) -> None:
    """Write a minimal coords TSV (only the columns the assessor reads)."""
    headers = (
        "contig_id",
        "orf_id",
        "start",
        "end",
        "strand",
        "partial_left",
        "partial_right",
        "contig_length",
    )
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(headers) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(h, "") for h in headers) + "\n")


def test_constants_match_mimag():
    """Constants encode the Bowers 2017 MIMAG thresholds verbatim."""
    assert MIMAG_HIGH_QUALITY_COMPLETENESS == 90.0
    assert MIMAG_HIGH_QUALITY_CONTAMINATION == 5.0
    assert MIMAG_MEDIUM_QUALITY_COMPLETENESS == 50.0
    assert QC_HMM_EVALUE_THRESHOLD == 1e-10
    assert QC_HMM_DOMAIN_EVALUE_THRESHOLD == 1e-10
    assert QC_HMM_CPU_COUNT >= 1
    assert QC_TSV_COLUMNS[0] == "accession"
    assert "completeness" in QC_TSV_COLUMNS
    assert "qc_quality" in QC_TSV_COLUMNS


def test_classify_quality_ladder():
    """The MIMAG classifier assigns the right tier across the boundary."""
    assert _classify_quality(95.0, 1.0) == QC_QUALITY_HIGH
    assert _classify_quality(90.0, 5.0) == QC_QUALITY_HIGH
    assert _classify_quality(89.99, 1.0) == QC_QUALITY_MEDIUM
    assert _classify_quality(60.0, 9.0) == QC_QUALITY_MEDIUM
    assert _classify_quality(60.0, 10.0) == QC_QUALITY_LOW
    assert _classify_quality(49.0, 1.0) == QC_QUALITY_LOW


def test_compute_n50_basic():
    """N50 returns the mid-cumulative contig length."""
    assert _compute_n50([100, 50, 50, 30, 20]) == 50
    assert _compute_n50([1000]) == 1000
    assert _compute_n50([]) == 0


def test_assess_full_completeness(tmp_path: Path):
    """A FAA containing every marker scores 100% completeness, 0% contamination."""
    hmm_path = tmp_path / "markers.hmm"
    n = _build_marker_hmm_db(hmm_path)
    assert n == 3
    faa = tmp_path / "GCA_TEST.faa"
    _write_faa(
        faa,
        [
            ("GCA_TEST|ctg1|orf_000001", _MARKER_SEQUENCES["M1"]),
            ("GCA_TEST|ctg1|orf_000002", _MARKER_SEQUENCES["M2"]),
            ("GCA_TEST|ctg2|orf_000003", _MARKER_SEQUENCES["M3"]),
        ],
    )
    coords = tmp_path / "GCA_TEST.coords.tsv"
    _write_coords(
        coords,
        [
            {
                "contig_id": "ctg1",
                "orf_id": "GCA_TEST|ctg1|orf_000001",
                "start": "1",
                "end": "240",
                "strand": "+",
                "partial_left": "0",
                "partial_right": "0",
                "contig_length": "1500",
            },
            {
                "contig_id": "ctg1",
                "orf_id": "GCA_TEST|ctg1|orf_000002",
                "start": "300",
                "end": "540",
                "strand": "+",
                "partial_left": "0",
                "partial_right": "0",
                "contig_length": "1500",
            },
            {
                "contig_id": "ctg2",
                "orf_id": "GCA_TEST|ctg2|orf_000003",
                "start": "1",
                "end": "240",
                "strand": "+",
                "partial_left": "0",
                "partial_right": "0",
                "contig_length": "500",
            },
        ],
    )

    out_dir = tmp_path / "qc"
    result = assess_assembly(
        accession="GCA_TEST",
        faa_path=faa,
        coords_path=coords,
        hmm_db_path=hmm_path,
        output_dir=out_dir,
        marker_set_label="TestMarkers",
    )
    assert result.skipped is False
    assert result.n_markers_total == 3
    assert result.n_markers_found == 3
    assert result.n_markers_duplicated == 0
    assert result.completeness == pytest.approx(100.0)
    assert result.contamination == pytest.approx(0.0)
    assert result.qc_quality == QC_QUALITY_HIGH
    # N50 over [1500, 500] = 1500 (>= 50% cumulative on first contig).
    assert result.n50_contig == 1500
    assert result.genome_size_bp == 2000
    assert result.n_contigs == 2
    assert result.qc_tsv_path.exists()
    assert result.stamp_path.exists()


def test_assess_partial_and_duplicated(tmp_path: Path):
    """Missing markers reduce completeness; duplicated markers raise contamination."""
    hmm_path = tmp_path / "markers.hmm"
    _build_marker_hmm_db(hmm_path)
    faa = tmp_path / "GCA_PART.faa"
    # Two markers found, one missing; M1 is present twice (a duplicate).
    _write_faa(
        faa,
        [
            ("GCA_PART|ctg1|orf_000001", _MARKER_SEQUENCES["M1"]),
            ("GCA_PART|ctg1|orf_000002", _MARKER_SEQUENCES["M1"]),
            ("GCA_PART|ctg1|orf_000003", _MARKER_SEQUENCES["M2"]),
        ],
    )
    coords = tmp_path / "GCA_PART.coords.tsv"
    _write_coords(
        coords,
        [
            {
                "contig_id": "ctg1",
                "orf_id": f"GCA_PART|ctg1|orf_{i:06d}",
                "start": "1",
                "end": "240",
                "strand": "+",
                "partial_left": "0",
                "partial_right": "0",
                "contig_length": "800",
            }
            for i in range(1, 4)
        ],
    )
    out_dir = tmp_path / "qc"
    result = assess_assembly(
        accession="GCA_PART",
        faa_path=faa,
        coords_path=coords,
        hmm_db_path=hmm_path,
        output_dir=out_dir,
        marker_set_label="TestMarkers",
    )
    assert result.n_markers_total == 3
    assert result.n_markers_found == 2  # M1 + M2 found, M3 missing
    assert result.n_markers_duplicated == 1  # M1 duplicated
    assert result.completeness == pytest.approx(2 / 3 * 100.0, rel=1e-3)
    assert result.contamination == pytest.approx(1 / 3 * 100.0, rel=1e-3)


def test_assess_idempotent(tmp_path: Path):
    """A second invocation with unchanged inputs short-circuits via the stamp."""
    hmm_path = tmp_path / "markers.hmm"
    _build_marker_hmm_db(hmm_path)
    faa = tmp_path / "GCA_X.faa"
    _write_faa(faa, [("orf_x", _MARKER_SEQUENCES["M1"])])
    coords = tmp_path / "GCA_X.coords.tsv"
    _write_coords(
        coords,
        [
            {
                "contig_id": "ctg",
                "orf_id": "orf_x",
                "start": "1",
                "end": "240",
                "strand": "+",
                "partial_left": "0",
                "partial_right": "0",
                "contig_length": "1000",
            }
        ],
    )
    out_dir = tmp_path / "qc"
    first = assess_assembly(
        accession="GCA_X",
        faa_path=faa,
        coords_path=coords,
        hmm_db_path=hmm_path,
        output_dir=out_dir,
        marker_set_label="TestMarkers",
    )
    assert first.skipped is False
    second = assess_assembly(
        accession="GCA_X",
        faa_path=faa,
        coords_path=coords,
        hmm_db_path=hmm_path,
        output_dir=out_dir,
        marker_set_label="TestMarkers",
    )
    assert second.skipped is True
    # Re-loaded values come from the rounded TSV cell; tolerate the {:.2f} round-trip.
    assert second.completeness == pytest.approx(first.completeness, abs=0.01)


def test_assess_overwrite(tmp_path: Path):
    """``overwrite=True`` re-runs even when the stamp matches."""
    hmm_path = tmp_path / "markers.hmm"
    _build_marker_hmm_db(hmm_path)
    faa = tmp_path / "GCA_X.faa"
    _write_faa(faa, [("orf_x", _MARKER_SEQUENCES["M1"])])
    coords = tmp_path / "GCA_X.coords.tsv"
    _write_coords(
        coords,
        [
            {
                "contig_id": "c",
                "orf_id": "orf_x",
                "start": "1",
                "end": "240",
                "strand": "+",
                "partial_left": "0",
                "partial_right": "0",
                "contig_length": "100",
            }
        ],
    )
    out_dir = tmp_path / "qc"
    assess_assembly(
        accession="GCA_X",
        faa_path=faa,
        coords_path=coords,
        hmm_db_path=hmm_path,
        output_dir=out_dir,
        marker_set_label="TestMarkers",
    )
    forced = assess_assembly(
        accession="GCA_X",
        faa_path=faa,
        coords_path=coords,
        hmm_db_path=hmm_path,
        output_dir=out_dir,
        marker_set_label="TestMarkers",
        overwrite=True,
    )
    assert forced.skipped is False


def test_read_distinct_contig_lengths_dedupes(tmp_path: Path):
    """Repeated contig_id rows contribute exactly one length."""
    coords = tmp_path / "c.tsv"
    _write_coords(
        coords,
        [
            {
                "contig_id": "ctg1",
                "orf_id": "o1",
                "start": "1",
                "end": "100",
                "strand": "+",
                "partial_left": "0",
                "partial_right": "0",
                "contig_length": "1000",
            },
            {
                "contig_id": "ctg1",
                "orf_id": "o2",
                "start": "200",
                "end": "300",
                "strand": "+",
                "partial_left": "0",
                "partial_right": "0",
                "contig_length": "1000",
            },
            {
                "contig_id": "ctg2",
                "orf_id": "o3",
                "start": "1",
                "end": "100",
                "strand": "+",
                "partial_left": "0",
                "partial_right": "0",
                "contig_length": "500",
            },
        ],
    )
    lengths = _read_distinct_contig_lengths(coords)
    assert sorted(lengths) == [500, 1000]


def test_load_score_cutoffs_parses_tsv(tmp_path: Path):
    """The cutoff loader parses BUSCO-style ``<name>\\t<bit>`` lines."""
    p = tmp_path / "scores_cutoff"
    p.write_text("M0\t10.5\nM1\t100.0\n# comment\nbad_line_no_score\n", encoding="utf-8")
    cutoffs = load_score_cutoffs(p)
    assert cutoffs == {"M0": 10.5, "M1": 100.0}


def test_load_score_cutoffs_handles_none(tmp_path: Path):
    """``None`` and missing files yield an empty mapping."""
    assert load_score_cutoffs(None) == {}
    assert load_score_cutoffs(tmp_path / "missing") == {}


def test_assess_with_score_cutoffs_filters_weak_hits(tmp_path: Path):
    """A high per-HMM cutoff drops every hit and zeros out completeness."""
    hmm_path = tmp_path / "markers.hmm"
    _build_marker_hmm_db(hmm_path)
    faa = tmp_path / "GCA_T.faa"
    _write_faa(
        faa,
        [
            ("orf1", _MARKER_SEQUENCES["M1"]),
            ("orf2", _MARKER_SEQUENCES["M2"]),
            ("orf3", _MARKER_SEQUENCES["M3"]),
        ],
    )
    coords = tmp_path / "GCA_T.coords.tsv"
    _write_coords(
        coords,
        [
            {
                "contig_id": "ctg",
                "orf_id": f"orf{i}",
                "start": "1",
                "end": "240",
                "strand": "+",
                "partial_left": "0",
                "partial_right": "0",
                "contig_length": "1000",
            }
            for i in range(1, 4)
        ],
    )
    cutoff_file = tmp_path / "scores_cutoff"
    # Cutoffs higher than the bit-scores any of these single-sequence HMMs can
    # achieve when searched against themselves; the gate suppresses every hit.
    cutoff_file.write_text("M1\t9999\nM2\t9999\nM3\t9999\n", encoding="utf-8")
    out_dir = tmp_path / "qc"
    result = assess_assembly(
        accession="GCA_T",
        faa_path=faa,
        coords_path=coords,
        hmm_db_path=hmm_path,
        output_dir=out_dir,
        marker_set_label="TestMarkers",
        scores_cutoff_path=cutoff_file,
    )
    assert result.n_markers_found == 0
    assert result.completeness == 0.0
    assert result.qc_quality == QC_QUALITY_LOW


def test_assess_with_lenient_cutoffs_allows_hits(tmp_path: Path):
    """A low per-HMM cutoff lets every reported hit through, matching no-cutoff behavior."""
    hmm_path = tmp_path / "markers.hmm"
    _build_marker_hmm_db(hmm_path)
    faa = tmp_path / "GCA_T.faa"
    _write_faa(
        faa,
        [
            ("orf1", _MARKER_SEQUENCES["M1"]),
            ("orf2", _MARKER_SEQUENCES["M2"]),
            ("orf3", _MARKER_SEQUENCES["M3"]),
        ],
    )
    coords = tmp_path / "GCA_T.coords.tsv"
    _write_coords(
        coords,
        [
            {
                "contig_id": "ctg",
                "orf_id": f"orf{i}",
                "start": "1",
                "end": "240",
                "strand": "+",
                "partial_left": "0",
                "partial_right": "0",
                "contig_length": "1000",
            }
            for i in range(1, 4)
        ],
    )
    cutoff_file = tmp_path / "scores_cutoff"
    cutoff_file.write_text("M1\t0.0\nM2\t0.0\nM3\t0.0\n", encoding="utf-8")
    out_dir = tmp_path / "qc"
    result = assess_assembly(
        accession="GCA_T",
        faa_path=faa,
        coords_path=coords,
        hmm_db_path=hmm_path,
        output_dir=out_dir,
        marker_set_label="TestMarkers",
        scores_cutoff_path=cutoff_file,
    )
    assert result.n_markers_found == 3
    assert result.completeness == pytest.approx(100.0)


def test_assess_missing_inputs(tmp_path: Path):
    """Missing inputs raise FileNotFoundError up-front."""
    hmm_path = tmp_path / "markers.hmm"
    _build_marker_hmm_db(hmm_path)
    out_dir = tmp_path / "qc"
    with pytest.raises(FileNotFoundError):
        assess_assembly(
            accession="X",
            faa_path=tmp_path / "missing.faa",
            coords_path=tmp_path / "missing.coords.tsv",
            hmm_db_path=hmm_path,
            output_dir=out_dir,
            marker_set_label="X",
        )
