"""Tests for ``src.checkm2_shellscript_generator``.

Pure transformation tests - we never invoke CheckM2. The generated script is
syntax-checked with ``bash -n`` only.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from src.checkm2_shellscript_generator import (
    CHECKM2_RESULTS_DIRNAME,
    CHECKM2_THREADS_DEFAULT,
    NON_ASSEMBLY_SOURCES,
    generate_checkm2_script,
)


_MANIFEST_HEADER = (
    "accession\tsource\turl\tlocal_path\tsha256\tbytes_downloaded\tfetched_at\t"
    "status\tnote\tlabel\tscientific_name\tsample_accession\tassembly_level\t"
    "format\tmd5\textracted_dir\n"
)


def _write_manifest(path: Path, rows: list[str]) -> None:
    """Write a small Phase 1 manifest for the generator."""
    path.write_text(_MANIFEST_HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def _touch_faa(orfs_dir: Path, accession: str) -> None:
    """Create a stub ``<accession>.faa`` so the generator does not skip the row."""
    orfs_dir.mkdir(parents=True, exist_ok=True)
    (orfs_dir / f"{accession}.faa").write_text(">orf\nM\n", encoding="utf-8")


def test_constants_consistent():
    """Module-level constants are sane."""
    assert CHECKM2_THREADS_DEFAULT >= 1
    assert CHECKM2_RESULTS_DIRNAME == "checkm2_results"
    assert "SILVA" in NON_ASSEMBLY_SOURCES
    assert "GTDB" in NON_ASSEMBLY_SOURCES


def test_generate_script_lists_only_assemblies(tmp_path: Path):
    """Manifest rows from non-assembly sources are excluded."""
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(
        manifest,
        [
            "GCA_AAA\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\tsci\tsample\t-\t-\t-\t-",
            "GCA_BBB\tENA\thttps://y\t/tmp/y\t-\t-\t-\tok\t-\tlbl\tsci\tsample\t-\t-\t-\t-",
            "SILVA_NR99\tSILVA\thttps://s\t/tmp/s\t-\t-\t-\tok\t-\tlbl\t-\t-\t-\t-\t-\t-",
            "GCA_FAILED\tENA\thttps://z\t/tmp/z\t-\t-\t-\tfailed\terr\tlbl\tsci\tsample\t-\t-\t-\t-",
        ],
    )
    data_root = tmp_path / "data"
    _touch_faa(data_root / "orfs", "GCA_AAA")
    _touch_faa(data_root / "orfs", "GCA_BBB")
    script_path = tmp_path / "scripts" / "run_checkm2.sh"
    result = generate_checkm2_script(
        manifest_path=manifest,
        data_root=data_root,
        script_path=script_path,
    )
    assert result.n_assemblies == 2
    body = script_path.read_text(encoding="utf-8")
    assert "GCA_AAA" in body
    assert "GCA_BBB" in body
    assert "SILVA_NR99" not in body
    assert "GCA_FAILED" not in body
    assert "checkm2 predict" in body


def test_generated_script_passes_bash_syntax_check(tmp_path: Path):
    """The script must parse cleanly with ``bash -n`` (no execution)."""
    if not shutil.which("bash"):
        pytest.skip("bash not on PATH")
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(
        manifest,
        [
            "GCA_TEST\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\tsci\tsample\t-\t-\t-\t-",
        ],
    )
    data_root = tmp_path / "data"
    _touch_faa(data_root / "orfs", "GCA_TEST")
    script_path = tmp_path / "scripts" / "run_checkm2.sh"
    generate_checkm2_script(manifest_path=manifest, data_root=data_root, script_path=script_path)

    completed = subprocess.run(
        ["bash", "-n", str(script_path)],
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr


def test_skips_assemblies_without_faa(tmp_path: Path):
    """MAGs without a ``.faa`` are reported in ``skipped_accessions``."""
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(
        manifest,
        [
            "GCA_HAVE_FAA\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\tsci\tsample\t-\t-\t-\t-",
            "GCA_NO_FAA\tENA\thttps://y\t/tmp/y\t-\t-\t-\tok\t-\tlbl\tsci\tsample\t-\t-\t-\t-",
        ],
    )
    data_root = tmp_path / "data"
    _touch_faa(data_root / "orfs", "GCA_HAVE_FAA")
    script_path = tmp_path / "scripts" / "run_checkm2.sh"
    result = generate_checkm2_script(
        manifest_path=manifest,
        data_root=data_root,
        script_path=script_path,
    )
    assert result.n_assemblies == 1
    assert result.skipped_accessions == ["GCA_NO_FAA"]


def test_script_is_executable(tmp_path: Path):
    """The generated script gets the user/group/other execute bit."""
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(
        manifest,
        [
            "GCA_TEST\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\tsci\tsample\t-\t-\t-\t-",
        ],
    )
    data_root = tmp_path / "data"
    _touch_faa(data_root / "orfs", "GCA_TEST")
    script_path = tmp_path / "scripts" / "run_checkm2.sh"
    generate_checkm2_script(manifest_path=manifest, data_root=data_root, script_path=script_path)
    mode = os.stat(script_path).st_mode
    assert mode & 0o111  # at least one of u/g/o execute


def test_script_demands_checkm2db_env_var(tmp_path: Path):
    """The script body references ``CHECKM2DB`` so the run aborts without it."""
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(
        manifest,
        [
            "GCA_TEST\tENA\thttps://x\t/tmp/x\t-\t-\t-\tok\t-\tlbl\tsci\tsample\t-\t-\t-\t-",
        ],
    )
    data_root = tmp_path / "data"
    _touch_faa(data_root / "orfs", "GCA_TEST")
    script_path = tmp_path / "scripts" / "run_checkm2.sh"
    generate_checkm2_script(manifest_path=manifest, data_root=data_root, script_path=script_path)
    body = script_path.read_text(encoding="utf-8")
    assert "CHECKM2DB" in body
    assert "set -euo pipefail" in body
