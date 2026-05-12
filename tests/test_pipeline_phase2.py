"""Tests for the Phase 2 CLI orchestrator (``src.pipeline_phase2``).

Targets argparse wiring + exit-code semantics. Subprocess (barrnap) is mocked
the same way ``test_ssu_extractor`` does it, so the suite stays offline and
barrnap-free.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import pytest

from src import pipeline_phase2 as p15
from src import ssu_extractor as ssu


_MANIFEST_HEADER = (
    "accession\tsource\turl\tlocal_path\tsha256\tbytes_downloaded\tfetched_at\t"
    "status\tnote\tlabel\tscientific_name\tsample_accession\tassembly_level\t"
    "format\tmd5\textracted_dir\n"
)


def _write_manifest(path: Path, rows: list[str]) -> None:
    """Write a Phase 1 manifest fixture."""
    path.write_text(_MANIFEST_HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def _write_fasta(path: Path, *, gzipped: bool = False) -> None:
    """Write a tiny one-contig FASTA."""
    body = ">contig_1\nACGTACGTACGATCGATCGATCGATCGATCTTTTTTTTTT\n"
    if gzipped:
        with gzip.open(path, "wt") as fh:
            fh.write(body)
    else:
        path.write_text(body, encoding="utf-8")


@dataclass
class _FakeCompleted:
    returncode: int
    stdout: str
    stderr: str = ""


def _patch_barrnap(
    monkeypatch: pytest.MonkeyPatch,
    *,
    gff_text: str = "",
    available: bool = True,
) -> None:
    """Same shape as the helper in ``test_ssu_extractor`` but local."""
    def fake_which(name: str) -> str | None:
        return "/usr/local/bin/barrnap" if available else None

    def fake_run(argv, **kwargs):
        if len(argv) >= 2 and argv[1] == "--version":
            return _FakeCompleted(returncode=0, stdout="barrnap 0.9", stderr="")
        return _FakeCompleted(returncode=0, stdout=gff_text, stderr="")

    monkeypatch.setattr(ssu.shutil, "which", fake_which)
    monkeypatch.setattr(ssu.subprocess, "run", fake_run)


def test_extract_16S_returns_2_when_barrnap_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """The CLI surfaces a clear exit code 2 when barrnap is absent."""
    manifest = tmp_path / "manifest.tsv"
    fasta = tmp_path / "GCA_TEST.fasta"
    _write_fasta(fasta)
    _write_manifest(
        manifest,
        [
            f"GCA_TEST\tENA\thttps://x\t{fasta}\t-\t-\t-\tok\t-\tlbl\tsci\t-\t-\t-\t-\t-",
        ],
    )
    _patch_barrnap(monkeypatch, available=False)

    rc = p15.cmd_extract_16S(
        manifest_path=manifest,
        data_root=tmp_path / "data",
        overwrite=False,
        accession_filter=None,
        threads=1,
        kingdom="bac",
        evalue=1e-6,
        reject=0.25,
    )
    assert rc == 2


def test_extract_16S_writes_routing_log_and_per_mag_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """A successful run writes the routing log and a per-MAG 16S FASTA."""
    manifest = tmp_path / "manifest.tsv"
    fasta = tmp_path / "GCA_TEST.fasta"
    _write_fasta(fasta)
    _write_manifest(
        manifest,
        [
            f"GCA_TEST\tENA\thttps://x\t{fasta}\t-\t-\t-\tok\t-\tlbl\tsci\t-\t-\t-\t-\t-",
        ],
    )
    gff_text = (
        "##gff-version 3\n"
        "contig_1\tbarrnap:0.9\trRNA\t11\t30\t1e-50\t+\t.\t"
        "Name=16S_rRNA;product=16S ribosomal RNA\n"
    )
    _patch_barrnap(monkeypatch, gff_text=gff_text)

    data_root = tmp_path / "data"
    rc = p15.cmd_extract_16S(
        manifest_path=manifest,
        data_root=data_root,
        overwrite=False,
        accession_filter=None,
        threads=1,
        kingdom="bac",
        evalue=1e-6,
        reject=0.25,
    )
    assert rc == 0
    assert (data_root / "ssu" / "GCA_TEST.16S.fasta").exists()
    assert (data_root / "ssu" / "routing_log.tsv").exists()
    log_body = (data_root / "ssu" / "routing_log.tsv").read_text(encoding="utf-8")
    assert "GCA_TEST" in log_body
    assert ssu.SSU_DECISION_SILVA in log_body
    # Phase 2 manifest gets a row, too.
    phase15 = data_root / "phase2_manifest.tsv"
    assert phase15.exists()
    assert "extract-16S" in phase15.read_text(encoding="utf-8")


def test_status_smoke_runs_with_no_outputs(tmp_path: Path, capsys):
    """``status`` runs cleanly even when nothing has been generated."""
    rc = p15.cmd_status(tmp_path / "data")
    assert rc == 0
    out = capsys.readouterr().out
    assert "Phase 2 status" in out


def test_bundle_silva_returns_2_without_routing_log(tmp_path: Path):
    """``bundle-silva`` exits 2 if the routing log is missing."""
    rc = p15.cmd_bundle_silva(tmp_path / "data")
    assert rc == 2


def test_gen_gtdbtk_script_returns_2_without_routing_log(tmp_path: Path):
    """``gen-gtdbtk-script`` exits 2 if the routing log is missing."""
    rc = p15.cmd_gen_gtdbtk_script(tmp_path / "data", threads=4)
    assert rc == 2


def test_main_argparse_smoke():
    """The argparse builder parses every advertised subcommand."""
    parser = p15._build_parser()
    for argv in (
        ["extract-16S"],
        ["bundle-silva"],
        ["gen-gtdbtk-script"],
        ["all"],
        ["status"],
    ):
        ns = parser.parse_args(argv)
        assert ns.cmd == argv[0]
