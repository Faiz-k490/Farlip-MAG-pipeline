"""Tests for ``src.orf_predictor``.

These tests are offline and use small in-memory FASTAs so they finish in well
under a second. The goal is to verify the I/O contract (output paths, stamp
file, idempotent re-run, gzip transparency, empty-input edge case) and the
fundamental pyrodigal integration - not to re-test pyrodigal itself.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from src.orf_predictor import (
    CHUNK_BYTES,
    COORDS_HEADER,
    FAA_LINE_WIDTH,
    PYRODIGAL_META_MODE,
    OrfPredictionResult,
    predict_orfs,
)


_SAMPLE_CONTIG_SEQ = (
    "ATGAAAGCAGGTATTGCTGCTGCAGCAGCTGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCT"
    "GCAGCAGCAGCAGCAGCAGCAGCTAAGCGTGTTGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCT"
    "GCAGCTGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATAAGGCCGGCTGAATAAA"
    "ATGGGAGAACGTGCTGAATTAGAAGAATTAATTAATGAATTAGAATTAATCGAACTAGAACTAGAA"
    "GAAGAACTAGAAGAACTAGAAGAACTAGAACTAGAAGAACTAGAACTAGAAGAACTAGAACTATAA"
) * 50


def _write_sample_fasta(path: Path, contigs: dict[str, str]) -> None:
    """Write a small multi-contig FASTA fixture."""
    with path.open("w", encoding="utf-8") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i : i + 70] + "\n")


def test_constants_unchanged():
    """Constants serve as the only knobs - guard their values so changes are explicit."""
    assert PYRODIGAL_META_MODE is True
    assert FAA_LINE_WIDTH == 60
    assert CHUNK_BYTES == 1024 * 1024
    assert COORDS_HEADER[0] == "contig_id"
    assert COORDS_HEADER[1] == "orf_id"
    assert COORDS_HEADER[-1] == "contig_length"


def test_predict_orfs_writes_three_outputs(tmp_path: Path):
    """A normal run writes .faa, .coords.tsv, and .stamp under the requested dir."""
    fasta = tmp_path / "demo.fna"
    _write_sample_fasta(fasta, {"ctgA": _SAMPLE_CONTIG_SEQ, "ctgB": _SAMPLE_CONTIG_SEQ})
    out_dir = tmp_path / "orfs"
    result = predict_orfs("DEMO", fasta, out_dir)
    assert isinstance(result, OrfPredictionResult)
    assert result.skipped is False
    assert result.faa_path.exists()
    assert result.coords_path.exists()
    assert result.stamp_path.exists()
    assert result.n_contigs == 2
    assert result.n_orfs > 0
    header_line = result.coords_path.read_text(encoding="utf-8").splitlines()[0]
    assert header_line == "\t".join(COORDS_HEADER)


def test_predict_orfs_is_idempotent(tmp_path: Path):
    """A second invocation with the same input short-circuits via the stamp."""
    fasta = tmp_path / "demo.fna"
    _write_sample_fasta(fasta, {"ctgA": _SAMPLE_CONTIG_SEQ})
    out_dir = tmp_path / "orfs"
    first = predict_orfs("DEMO", fasta, out_dir)
    assert first.skipped is False
    second = predict_orfs("DEMO", fasta, out_dir)
    assert second.skipped is True
    assert second.n_orfs == first.n_orfs
    assert second.input_sha256 == first.input_sha256


def test_predict_orfs_reruns_after_input_change(tmp_path: Path):
    """Mutating the input FASTA must force a re-run even with the stamp present."""
    fasta = tmp_path / "demo.fna"
    _write_sample_fasta(fasta, {"ctgA": _SAMPLE_CONTIG_SEQ})
    out_dir = tmp_path / "orfs"
    first = predict_orfs("DEMO", fasta, out_dir)
    _write_sample_fasta(
        fasta, {"ctgA": _SAMPLE_CONTIG_SEQ, "ctgB": _SAMPLE_CONTIG_SEQ}
    )
    second = predict_orfs("DEMO", fasta, out_dir)
    assert second.skipped is False
    assert second.n_contigs == 2
    assert second.input_sha256 != first.input_sha256


def test_predict_orfs_overwrite_flag(tmp_path: Path):
    """``overwrite=True`` re-runs even when the stamp matches."""
    fasta = tmp_path / "demo.fna"
    _write_sample_fasta(fasta, {"ctgA": _SAMPLE_CONTIG_SEQ})
    out_dir = tmp_path / "orfs"
    predict_orfs("DEMO", fasta, out_dir)
    forced = predict_orfs("DEMO", fasta, out_dir, overwrite=True)
    assert forced.skipped is False


def test_predict_orfs_handles_gzipped_input(tmp_path: Path):
    """``.fna.gz`` inputs are opened transparently and produce identical ORF counts."""
    plain = tmp_path / "demo.fna"
    gzipped = tmp_path / "demo.fna.gz"
    _write_sample_fasta(plain, {"ctgA": _SAMPLE_CONTIG_SEQ})
    with gzip.open(gzipped, "wt") as fh:
        fh.write(plain.read_text(encoding="utf-8"))
    out_plain = predict_orfs("PLAIN", plain, tmp_path / "orfs_plain")
    out_gz = predict_orfs("GZ", gzipped, tmp_path / "orfs_gz")
    assert out_plain.n_orfs == out_gz.n_orfs


def test_predict_orfs_missing_input_raises(tmp_path: Path):
    """A missing input path raises FileNotFoundError up-front."""
    with pytest.raises(FileNotFoundError):
        predict_orfs("DEMO", tmp_path / "nonexistent.fna", tmp_path / "orfs")


def test_predict_orfs_empty_contig_skipped(tmp_path: Path):
    """A FASTA whose contig has zero length contributes zero ORFs and zero contigs."""
    fasta = tmp_path / "empty.fna"
    fasta.write_text(">ctg_empty\n\n", encoding="utf-8")
    result = predict_orfs("EMPTY", fasta, tmp_path / "orfs")
    assert result.n_contigs == 0
    assert result.n_orfs == 0
