"""Tests for ``src.ssu_extractor``.

The barrnap binary is mocked via ``monkeypatch.setattr(subprocess, "run", ...)``
so the suite stays offline and barrnap-free. A small synthetic FASTA + GFF3
fixture exercises every code path: 16S found, 16S not found, gzipped input
support, reverse-complement slicing, and idempotent re-runs via the stamp.
"""

from __future__ import annotations

import gzip
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest

from src import ssu_extractor as ssu
from src.ssu_extractor import (
    BARRNAP_BINARY,
    BarrnapError,
    SSU_DECISION_GTDBTK,
    SSU_DECISION_SILVA,
    Ssu16SRecord,
    assess_mag_for_ssu,
    extract_16S,
    run_barrnap,
    write_ssu_fasta,
)


_FAKE_BARRNAP_BIN = "fake-barrnap"
"""Stub binary name used in tests; never invoked because we monkeypatch
``shutil.which`` and ``subprocess.run``."""


def _seq() -> str:
    """Return a tiny synthetic contig with a fake 16S region embedded.

    The 16S coordinates the GFF will reference are 11..30 (1-based inclusive).
    """
    prefix = "ACGTACGTAC"  # 10 bp upstream
    ssu_region = "GATCGATCGATCGATCGATC"  # 20 bp "16S"
    suffix = "TTTTTTTTTT"  # 10 bp downstream
    return prefix + ssu_region + suffix  # total 40 bp


def _write_test_fasta(path: Path, *, gzipped: bool = False) -> None:
    """Write a one-contig test FASTA (plain or gzipped)."""
    body = f">contig_1 description\n{_seq()}\n"
    if gzipped:
        with gzip.open(path, "wt") as fh:
            fh.write(body)
    else:
        path.write_text(body, encoding="utf-8")


def _make_gff_text(*, with_16S: bool = True) -> str:
    """Render the barrnap GFF3 stdout fixture (with or without a 16S row)."""
    header = "##gff-version 3\n"
    if not with_16S:
        return header + (
            "contig_1\tbarrnap:0.9\trRNA\t1\t10\t1e-30\t+\t.\t"
            "Name=5S_rRNA;product=5S ribosomal RNA\n"
        )
    return header + (
        "contig_1\tbarrnap:0.9\trRNA\t11\t30\t1e-50\t+\t.\t"
        "Name=16S_rRNA;product=16S ribosomal RNA\n"
    )


@dataclass
class _FakeCompleted:
    """Stand-in for ``subprocess.CompletedProcess``."""

    returncode: int
    stdout: str
    stderr: str = ""


def _patch_barrnap(
    monkeypatch: pytest.MonkeyPatch,
    *,
    gff_text: str,
    returncode: int = 0,
    stderr: str = "",
    available: bool = True,
) -> list[list[str]]:
    """Monkeypatch ``shutil.which`` and ``subprocess.run`` so barrnap is a stub.

    Args:
        monkeypatch: pytest fixture.
        gff_text: GFF3 stdout the fake barrnap will return.
        returncode: Exit code of the fake invocation.
        stderr: stderr the fake invocation will return.
        available: Whether the fake binary is reported as on PATH.

    Returns:
        A captured list of every argv handed to the fake ``subprocess.run``.
        Useful for asserting that barrnap got the right arguments.
    """
    captured: list[list[str]] = []

    def fake_which(name: str) -> str | None:
        if name == _FAKE_BARRNAP_BIN and available:
            return f"/usr/local/bin/{name}"
        if name == BARRNAP_BINARY and available:
            return f"/usr/local/bin/{name}"
        return None

    def fake_run(argv, **kwargs):
        captured.append(list(argv))
        is_version_call = len(argv) >= 2 and argv[1] == "--version"
        if is_version_call:
            return _FakeCompleted(returncode=0, stdout="barrnap 0.9", stderr="")
        return _FakeCompleted(returncode=returncode, stdout=gff_text, stderr=stderr)

    monkeypatch.setattr(ssu.shutil, "which", fake_which)
    monkeypatch.setattr(ssu.subprocess, "run", fake_run)
    return captured


def test_run_barrnap_writes_gff_and_counts_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Happy path: GFF stdout lands on disk and 16S rows are counted."""
    fasta = tmp_path / "mag.fasta"
    _write_test_fasta(fasta)
    captured = _patch_barrnap(monkeypatch, gff_text=_make_gff_text())
    out_gff = tmp_path / "out.gff"

    result = run_barrnap(fasta, out_gff, barrnap_binary=_FAKE_BARRNAP_BIN)

    assert out_gff.exists()
    assert "16S_rRNA" in out_gff.read_text()
    assert result.n_rows_total == 1
    assert result.n_rows_16S == 1
    invocation = next(c for c in captured if c[0] == _FAKE_BARRNAP_BIN and "--version" not in c)
    assert invocation[0] == _FAKE_BARRNAP_BIN
    assert "--kingdom" in invocation
    assert "--threads" in invocation
    assert "--quiet" in invocation


def test_run_barrnap_handles_gzipped_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Gzipped FASTAs are decompressed to a tempfile before the call."""
    fasta = tmp_path / "mag.fasta.gz"
    _write_test_fasta(fasta, gzipped=True)
    captured = _patch_barrnap(monkeypatch, gff_text=_make_gff_text())
    out_gff = tmp_path / "out.gff"

    run_barrnap(fasta, out_gff, barrnap_binary=_FAKE_BARRNAP_BIN)

    invocation = next(c for c in captured if c[0] == _FAKE_BARRNAP_BIN and "--version" not in c)
    fasta_arg = invocation[-1]
    assert fasta_arg.endswith(".fasta")
    # The temp file should not still exist after run_barrnap returns.
    assert not Path(fasta_arg).exists()


def test_run_barrnap_raises_when_binary_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """A missing barrnap binary raises a clear BarrnapError."""
    fasta = tmp_path / "mag.fasta"
    _write_test_fasta(fasta)
    _patch_barrnap(
        monkeypatch, gff_text=_make_gff_text(), available=False
    )
    with pytest.raises(BarrnapError):
        run_barrnap(fasta, tmp_path / "out.gff", barrnap_binary=_FAKE_BARRNAP_BIN)


def test_run_barrnap_raises_on_nonzero_exit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """A non-zero barrnap exit code becomes a BarrnapError."""
    fasta = tmp_path / "mag.fasta"
    _write_test_fasta(fasta)
    _patch_barrnap(
        monkeypatch,
        gff_text="",
        returncode=2,
        stderr="boom",
    )
    with pytest.raises(BarrnapError):
        run_barrnap(fasta, tmp_path / "out.gff", barrnap_binary=_FAKE_BARRNAP_BIN)


def test_extract_16S_returns_correct_slice(tmp_path: Path):
    """The slice and metadata match the GFF coordinates."""
    fasta = tmp_path / "mag.fasta"
    _write_test_fasta(fasta)
    gff = tmp_path / "barrnap.gff"
    gff.write_text(_make_gff_text(), encoding="utf-8")

    records = extract_16S(fasta, gff)

    assert len(records) == 1
    rec = records[0]
    assert isinstance(rec, Ssu16SRecord)
    assert rec.contig_id == "contig_1"
    assert rec.start == 11
    assert rec.end == 30
    assert rec.length_bp == 20
    assert rec.strand == "+"
    assert rec.sequence == _seq()[10:30].upper()


def test_extract_16S_reverse_complements_minus_strand(tmp_path: Path):
    """A ``-`` strand row yields the reverse complement of the slice."""
    fasta = tmp_path / "mag.fasta"
    _write_test_fasta(fasta)
    gff = tmp_path / "barrnap.gff"
    gff.write_text(
        "##gff-version 3\n"
        "contig_1\tbarrnap:0.9\trRNA\t11\t30\t1e-50\t-\t.\t"
        "Name=16S_rRNA;product=16S ribosomal RNA\n",
        encoding="utf-8",
    )
    records = extract_16S(fasta, gff)
    forward = _seq()[10:30].upper()
    expected = forward.translate(
        str.maketrans(
            "ACGTUNRYSWKMBDHVacgtunryswkmbdhv",
            "TGCAANYRSWMKVHDBtgcaanyrswmkvhdb",
        )
    )[::-1]
    assert records[0].sequence == expected


def test_extract_16S_empty_when_no_16S_rows(tmp_path: Path):
    """A barrnap GFF with no 16S rows produces an empty record list."""
    fasta = tmp_path / "mag.fasta"
    _write_test_fasta(fasta)
    gff = tmp_path / "barrnap.gff"
    gff.write_text(_make_gff_text(with_16S=False), encoding="utf-8")
    assert extract_16S(fasta, gff) == []


def test_write_ssu_fasta_emits_one_record_per_input(tmp_path: Path):
    """The FASTA writer emits one header per record and is always created."""
    out = tmp_path / "MAG.16S.fasta"
    rec = Ssu16SRecord(
        contig_id="ctg1",
        start=1,
        end=20,
        strand="+",
        length_bp=20,
        sequence="ACGT" * 5,
        attributes="Name=16S_rRNA",
    )
    n = write_ssu_fasta([rec], out, mag_accession="GCA_TEST")
    assert n == 1
    body = out.read_text(encoding="utf-8")
    assert body.startswith(">GCA_TEST|ctg1|1-20|+|len=20")
    assert "ACGT" in body


def test_write_ssu_fasta_creates_empty_file_for_empty_records(tmp_path: Path):
    """A zero-record list still produces an (empty) FASTA file."""
    out = tmp_path / "MAG.16S.fasta"
    n = write_ssu_fasta([], out, mag_accession="GCA_TEST")
    assert n == 0
    assert out.exists()
    assert out.read_text(encoding="utf-8") == ""


def test_assess_mag_for_ssu_routes_to_silva_when_16S_found(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """At least one 16S record yields a SILVA decision."""
    fasta = tmp_path / "GCA_TEST.fasta"
    _write_test_fasta(fasta)
    _patch_barrnap(monkeypatch, gff_text=_make_gff_text())

    out_dir = tmp_path / "ssu"
    assessment = assess_mag_for_ssu(
        accession="GCA_TEST",
        fasta_path=fasta,
        output_dir=out_dir,
        barrnap_binary=_FAKE_BARRNAP_BIN,
    )

    assert assessment.decision == SSU_DECISION_SILVA
    assert assessment.n_16S_predicted == 1
    assert assessment.skipped is False
    assert assessment.ssu_fasta_path.exists()
    assert assessment.gff_path.exists()
    assert assessment.stamp_path.exists()


def test_assess_mag_for_ssu_routes_to_gtdbtk_when_no_16S(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Zero 16S records yield a GTDB-Tk decision and an empty 16S FASTA."""
    fasta = tmp_path / "GCA_TEST.fasta"
    _write_test_fasta(fasta)
    _patch_barrnap(monkeypatch, gff_text=_make_gff_text(with_16S=False))

    assessment = assess_mag_for_ssu(
        accession="GCA_TEST",
        fasta_path=fasta,
        output_dir=tmp_path / "ssu",
        barrnap_binary=_FAKE_BARRNAP_BIN,
    )

    assert assessment.decision == SSU_DECISION_GTDBTK
    assert assessment.n_16S_predicted == 0
    assert assessment.ssu_fasta_path.exists()
    assert assessment.ssu_fasta_path.read_text(encoding="utf-8") == ""


def test_assess_mag_for_ssu_idempotent_re_run_short_circuits(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """A second run with the same FASTA + parameters skips barrnap entirely."""
    fasta = tmp_path / "GCA_TEST.fasta"
    _write_test_fasta(fasta)
    captured = _patch_barrnap(monkeypatch, gff_text=_make_gff_text())

    out_dir = tmp_path / "ssu"
    first = assess_mag_for_ssu(
        accession="GCA_TEST",
        fasta_path=fasta,
        output_dir=out_dir,
        barrnap_binary=_FAKE_BARRNAP_BIN,
    )
    assert first.skipped is False
    n_calls_after_first = sum(1 for c in captured if "--version" not in c)

    second = assess_mag_for_ssu(
        accession="GCA_TEST",
        fasta_path=fasta,
        output_dir=out_dir,
        barrnap_binary=_FAKE_BARRNAP_BIN,
    )
    assert second.skipped is True
    assert second.decision == SSU_DECISION_SILVA
    n_calls_after_second = sum(1 for c in captured if "--version" not in c)
    assert n_calls_after_second == n_calls_after_first  # no extra barrnap call
