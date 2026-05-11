"""Tests for ``src.gtdbtk_shellscript_generator``.

We never invoke GTDB-Tk; the generated script is syntax-checked with
``bash -n``. Routing-log fixtures are fabricated in ``tmp_path``.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from src.gtdbtk_shellscript_generator import (
    GTDBTK_BATCHFILE_NAME,
    GTDBTK_FASTA_EXTENSION,
    GTDBTK_RESULTS_DIRNAME,
    GTDBTK_SCRIPT_RELPATH,
    GTDBTK_THREADS_DEFAULT,
    generate_gtdbtk_script,
)
from src.ssu_bundler import RoutingDecision, write_routing_log
from src.ssu_extractor import SSU_DECISION_GTDBTK, SSU_DECISION_SILVA


def _decision(
    accession: str,
    fasta_path: Path,
    decision: str,
) -> RoutingDecision:
    """Construct a fake routing decision for fixtures."""
    return RoutingDecision(
        accession=accession,
        source="ENA",
        fasta_path=fasta_path,
        n_contigs_scanned=10,
        n_16S_predicted=0 if decision == SSU_DECISION_GTDBTK else 1,
        decision=decision,
        output_path=fasta_path,
        barrnap_version="barrnap 0.9",
        notes="",
    )


def _write_fixture_fasta(path: Path, *, gzipped: bool = False) -> None:
    """Write a one-byte-of-sequence FASTA fixture (plain or gzipped)."""
    body = ">contig\nACGT\n"
    if gzipped:
        with gzip.open(path, "wt") as fh:
            fh.write(body)
    else:
        path.write_text(body, encoding="utf-8")


def test_constants_consistent():
    """Module-level constants are sane."""
    assert GTDBTK_THREADS_DEFAULT >= 1
    assert GTDBTK_RESULTS_DIRNAME == "gtdbtk_results"
    assert GTDBTK_BATCHFILE_NAME.endswith(".tsv")
    assert GTDBTK_FASTA_EXTENSION == "fna"
    assert GTDBTK_SCRIPT_RELPATH.endswith("run_gtdbtk.sh")


def test_generate_script_filters_to_gtdbtk_routed_only(tmp_path: Path):
    """SILVA-routed MAGs are excluded; GTDB-Tk-routed MAGs land in the batchfile."""
    data_root = tmp_path / "data"
    fasta_a = tmp_path / "MAG_A.fna"
    fasta_b = tmp_path / "MAG_B.fna"
    _write_fixture_fasta(fasta_a)
    _write_fixture_fasta(fasta_b)

    log = data_root / "ssu" / "routing_log.tsv"
    write_routing_log(
        [
            _decision("MAG_A", fasta_a, SSU_DECISION_GTDBTK),
            _decision("MAG_B", fasta_b, SSU_DECISION_SILVA),
        ],
        log,
    )

    script_path = tmp_path / "scripts" / "run_gtdbtk.sh"
    result = generate_gtdbtk_script(
        routing_log_path=log,
        data_root=data_root,
        script_path=script_path,
    )

    assert result.n_assemblies == 1
    body = script_path.read_text(encoding="utf-8")
    batchfile = result.batchfile_path.read_text(encoding="utf-8")
    assert "MAG_A" in batchfile
    assert "MAG_B" not in batchfile
    assert "gtdbtk classify_wf" in body
    assert "GTDBTK_DATA_PATH" in body
    assert "set -euo pipefail" in body


def test_generate_script_decompresses_gzipped_inputs(tmp_path: Path):
    """A ``.fasta.gz`` row is materialised to ``data/ssu/gtdbtk_input/<acc>.fna``."""
    data_root = tmp_path / "data"
    fasta = tmp_path / "MAG_A.fasta.gz"
    _write_fixture_fasta(fasta, gzipped=True)
    log = data_root / "ssu" / "routing_log.tsv"
    write_routing_log(
        [_decision("MAG_A", fasta, SSU_DECISION_GTDBTK)],
        log,
    )

    script_path = tmp_path / "scripts" / "run_gtdbtk.sh"
    result = generate_gtdbtk_script(
        routing_log_path=log,
        data_root=data_root,
        script_path=script_path,
    )

    assert result.n_assemblies == 1
    staged = data_root / "ssu" / "gtdbtk_input" / "MAG_A.fna"
    assert staged.exists()
    assert staged.read_text(encoding="utf-8").startswith(">contig")
    batchfile = result.batchfile_path.read_text(encoding="utf-8")
    assert str(staged.resolve()) in batchfile


def test_generate_script_skips_missing_fastas(tmp_path: Path):
    """Routing rows with a missing FASTA are reported but do not crash."""
    data_root = tmp_path / "data"
    log = data_root / "ssu" / "routing_log.tsv"
    write_routing_log(
        [_decision("MAG_GHOST", tmp_path / "missing.fna", SSU_DECISION_GTDBTK)],
        log,
    )
    script_path = tmp_path / "scripts" / "run_gtdbtk.sh"
    result = generate_gtdbtk_script(
        routing_log_path=log,
        data_root=data_root,
        script_path=script_path,
    )
    assert result.n_assemblies == 0
    assert result.skipped_accessions == ["MAG_GHOST"]


def test_generated_script_passes_bash_syntax_check(tmp_path: Path):
    """``bash -n`` accepts the generated script."""
    if not shutil.which("bash"):
        pytest.skip("bash not on PATH")
    data_root = tmp_path / "data"
    fasta = tmp_path / "MAG_A.fna"
    _write_fixture_fasta(fasta)
    log = data_root / "ssu" / "routing_log.tsv"
    write_routing_log(
        [_decision("MAG_A", fasta, SSU_DECISION_GTDBTK)],
        log,
    )
    script_path = tmp_path / "scripts" / "run_gtdbtk.sh"
    generate_gtdbtk_script(
        routing_log_path=log,
        data_root=data_root,
        script_path=script_path,
    )
    completed = subprocess.run(
        ["bash", "-n", str(script_path)],
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr


def test_generated_script_is_executable(tmp_path: Path):
    """The script gets the user/group/other execute bit."""
    data_root = tmp_path / "data"
    fasta = tmp_path / "MAG_A.fna"
    _write_fixture_fasta(fasta)
    log = data_root / "ssu" / "routing_log.tsv"
    write_routing_log(
        [_decision("MAG_A", fasta, SSU_DECISION_GTDBTK)],
        log,
    )
    script_path = tmp_path / "scripts" / "run_gtdbtk.sh"
    generate_gtdbtk_script(
        routing_log_path=log,
        data_root=data_root,
        script_path=script_path,
    )
    mode = os.stat(script_path).st_mode
    assert mode & 0o111
