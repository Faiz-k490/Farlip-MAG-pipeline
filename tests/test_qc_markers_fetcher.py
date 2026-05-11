"""Tests for ``src.qc_markers_fetcher``.

Offline-only. We avoid hitting the BUSCO mirror by:

* For the ``local_directory`` archive format: pointing the spec at a tmpdir of
  fake ``.hmm`` files - this exercises the ``_gather_hmm_files`` and
  ``_concatenate_hmms`` paths without any network.
* For the ``tar.gz`` format: building a small in-memory tar fixture, then
  monkey-patching :func:`src.qc_markers_fetcher.stream_download` so it copies
  that fixture into place rather than calling out to the internet.
"""

from __future__ import annotations

import hashlib
import io
import tarfile
from pathlib import Path

import pytest

from src import qc_markers_fetcher
from src.qc_markers_fetcher import (
    ARCHIVE_FORMAT_LOCAL_DIR,
    ARCHIVE_FORMAT_TAR_GZ,
    MarkerFetcherError,
    MarkerSetSpec,
    fetch_marker_set,
    load_marker_config,
    resolve_default_spec,
)


_CONFIG_YAML = """
schema_version: "1.0"
default_marker_set: "test_busco"
marker_sets:
  test_busco:
    label: "Test_BUSCO"
    version: "1.0"
    url: "https://example.invalid/test.tar.gz"
    sha256: ""
    expected_marker_count: 3
    archive_format: "tar.gz"
    hmm_subdir: "lineage/hmms"
    description: "Test marker set"
  local_only:
    label: "Local_Only"
    archive_format: "local_directory"
    expected_marker_count: 2
    description: "User-supplied local marker dir"
"""


def _write_fake_hmm(path: Path, name: str) -> None:
    """Write a tiny but well-formed HMMER3 stub file.

    The fetcher only concatenates and stamps; pyhmmer is not invoked here so
    we don't need a real-format HMM, just a non-empty file that ends with a
    newline.
    """
    path.write_text(
        f"HMMER3/f [3.4 | placeholder]\nNAME  {name}\n//\n",
        encoding="utf-8",
    )


def _make_local_marker_dir(tmp_path: Path, n: int = 3) -> Path:
    """Materialize ``n`` fake .hmm files in a fresh subdirectory."""
    d = tmp_path / "local_markers"
    d.mkdir()
    for i in range(n):
        _write_fake_hmm(d / f"marker_{i}.hmm", f"M{i}")
    return d


def _write_config(tmp_path: Path) -> Path:
    """Write the test ``qc_markers.yaml`` and return its path."""
    p = tmp_path / "qc_markers.yaml"
    p.write_text(_CONFIG_YAML, encoding="utf-8")
    return p


def _make_test_tarball(tmp_path: Path, n: int = 3) -> tuple[Path, str]:
    """Create a tar.gz containing ``lineage/hmms/marker_*.hmm``.

    Returns:
        (archive path, sha256 hex digest of the archive bytes).
    """
    raw_buf = io.BytesIO()
    with tarfile.open(fileobj=raw_buf, mode="w:gz") as tf:
        for i in range(n):
            content = f"HMMER3/f [3.4 | placeholder]\nNAME  M{i}\n//\n".encode("utf-8")
            info = tarfile.TarInfo(name=f"lineage/hmms/marker_{i}.hmm")
            info.size = len(content)
            tf.addfile(info, io.BytesIO(content))
    archive = tmp_path / "fake.tar.gz"
    data = raw_buf.getvalue()
    archive.write_bytes(data)
    return archive, hashlib.sha256(data).hexdigest()


def test_load_marker_config_round_trips(tmp_path: Path):
    """The YAML file parses into one :class:`MarkerSetSpec` per entry."""
    config = _write_config(tmp_path)
    specs = load_marker_config(config)
    assert set(specs) == {"test_busco", "local_only"}
    assert specs["test_busco"].url.endswith("test.tar.gz")
    assert specs["local_only"].archive_format == ARCHIVE_FORMAT_LOCAL_DIR


def test_resolve_default_spec(tmp_path: Path):
    """``default_marker_set`` selects the correct entry."""
    config = _write_config(tmp_path)
    spec = resolve_default_spec(config)
    assert spec.key == "test_busco"


def test_resolve_default_missing_default(tmp_path: Path):
    """A YAML with no ``default_marker_set`` raises a clear error."""
    bad = tmp_path / "qc.yaml"
    bad.write_text("schema_version: '1.0'\nmarker_sets: {}\n", encoding="utf-8")
    with pytest.raises(MarkerFetcherError):
        resolve_default_spec(bad)


def _local_spec(local_dir: Path | str = "", n: int = 4) -> MarkerSetSpec:
    """Build a ``local_directory`` :class:`MarkerSetSpec` for tests."""
    return MarkerSetSpec(
        key="local_only",
        label="Local_Only",
        version="dev",
        url="",
        sha256="",
        expected_marker_count=n,
        archive_format=ARCHIVE_FORMAT_LOCAL_DIR,
        hmm_subdir="",
        local_marker_dir=str(local_dir),
        scores_cutoff_relpath="",
        description="test",
    )


def _tarball_spec(*, sha256: str, n: int = 5) -> MarkerSetSpec:
    """Build a ``tar.gz`` :class:`MarkerSetSpec` for tests."""
    return MarkerSetSpec(
        key="test_busco",
        label="Test_BUSCO",
        version="1.0",
        url="https://example.invalid/test.tar.gz",
        sha256=sha256,
        expected_marker_count=n,
        archive_format=ARCHIVE_FORMAT_TAR_GZ,
        hmm_subdir="lineage/hmms",
        local_marker_dir="",
        scores_cutoff_relpath="",
        description="test",
    )


def test_fetch_local_directory_format(tmp_path: Path):
    """``local_directory`` archives are concatenated without network access."""
    local_dir = _make_local_marker_dir(tmp_path, n=4)
    spec = _local_spec(local_dir=local_dir, n=4)
    data_root = tmp_path / "data"
    result = fetch_marker_set(spec, data_root)
    assert result.skipped is False
    assert result.n_markers == 4
    assert result.hmm_db_path.exists()
    assert result.stamp_path.exists()
    assert result.archive_sha256 == ""
    assert result.scores_cutoff_path is None
    again = fetch_marker_set(spec, data_root)
    assert again.skipped is True


def test_fetch_local_directory_requires_dir(tmp_path: Path):
    """An unset ``local_marker_dir`` raises a clear error."""
    spec = _local_spec(local_dir="", n=4)
    with pytest.raises(MarkerFetcherError):
        fetch_marker_set(spec, tmp_path / "data")


def test_fetch_local_directory_picks_up_scores_cutoff(tmp_path: Path):
    """A ``scores_cutoff`` file in the local directory is copied beside the DB."""
    local_dir = _make_local_marker_dir(tmp_path, n=3)
    (local_dir / "scores_cutoff").write_text("M0\t100.0\nM1\t50.5\nM2\t75.0\n", encoding="utf-8")
    spec = _local_spec(local_dir=local_dir, n=3)
    data_root = tmp_path / "data"
    result = fetch_marker_set(spec, data_root)
    assert result.scores_cutoff_path is not None
    assert result.scores_cutoff_path.exists()
    assert "M0\t100.0" in result.scores_cutoff_path.read_text(encoding="utf-8")


def test_fetch_tarball_format(tmp_path: Path, monkeypatch):
    """``tar.gz`` archives extract individual ``.hmm`` files into the database."""
    archive_src, archive_sha = _make_test_tarball(tmp_path, n=5)

    def _fake_stream_download(url, dest, **kwargs):
        """Pretend we downloaded ``url`` by copying our local fixture into ``dest``."""
        dest.parent.mkdir(parents=True, exist_ok=True)
        data = archive_src.read_bytes()
        dest.write_bytes(data)
        return hashlib.sha256(data).hexdigest(), len(data)

    monkeypatch.setattr(qc_markers_fetcher, "stream_download", _fake_stream_download)
    spec = _tarball_spec(sha256=archive_sha, n=5)
    data_root = tmp_path / "data"
    result = fetch_marker_set(spec, data_root)
    assert result.skipped is False
    assert result.n_markers == 5
    assert result.hmm_db_path.exists()
    assert result.archive_sha256 == archive_sha

    again = fetch_marker_set(spec, data_root)
    assert again.skipped is True


def test_fetch_tarball_picks_up_scores_cutoff(tmp_path: Path, monkeypatch):
    """A ``scores_cutoff`` member inside the tarball is auto-detected and copied."""
    raw_buf = io.BytesIO()
    cutoff_text = "M0\t100.0\nM1\t50.5\nM2\t75.0\n"
    with tarfile.open(fileobj=raw_buf, mode="w:gz") as tf:
        for i in range(3):
            content = f"HMMER3/f [3.4]\nNAME  M{i}\n//\n".encode("utf-8")
            info = tarfile.TarInfo(name=f"lineage/hmms/marker_{i}.hmm")
            info.size = len(content)
            tf.addfile(info, io.BytesIO(content))
        cutoff_bytes = cutoff_text.encode("utf-8")
        info = tarfile.TarInfo(name="lineage/scores_cutoff")
        info.size = len(cutoff_bytes)
        tf.addfile(info, io.BytesIO(cutoff_bytes))
    archive = tmp_path / "fixture.tar.gz"
    archive.write_bytes(raw_buf.getvalue())
    archive_sha = hashlib.sha256(raw_buf.getvalue()).hexdigest()

    def _fake_stream_download(url, dest, **kwargs):
        dest.parent.mkdir(parents=True, exist_ok=True)
        data = archive.read_bytes()
        dest.write_bytes(data)
        return hashlib.sha256(data).hexdigest(), len(data)

    monkeypatch.setattr(qc_markers_fetcher, "stream_download", _fake_stream_download)
    spec = MarkerSetSpec(
        key="test_busco",
        label="Test_BUSCO",
        version="1.0",
        url="https://example.invalid/test.tar.gz",
        sha256=archive_sha,
        expected_marker_count=3,
        archive_format=ARCHIVE_FORMAT_TAR_GZ,
        hmm_subdir="lineage/hmms",
        local_marker_dir="",
        scores_cutoff_relpath="lineage/scores_cutoff",
        description="test",
    )
    data_root = tmp_path / "data"
    result = fetch_marker_set(spec, data_root)
    assert result.scores_cutoff_path is not None
    assert result.scores_cutoff_path.exists()
    assert result.scores_cutoff_path.read_text(encoding="utf-8") == cutoff_text


def test_fetch_tarball_sha_mismatch_raises(tmp_path: Path, monkeypatch):
    """A wrong configured SHA-256 aborts the fetch with a clear error."""
    archive_src, _archive_sha = _make_test_tarball(tmp_path, n=2)

    def _fake_stream_download(url, dest, **kwargs):
        dest.parent.mkdir(parents=True, exist_ok=True)
        data = archive_src.read_bytes()
        dest.write_bytes(data)
        return hashlib.sha256(data).hexdigest(), len(data)

    monkeypatch.setattr(qc_markers_fetcher, "stream_download", _fake_stream_download)
    spec = _tarball_spec(sha256="0" * 64, n=2)
    with pytest.raises(MarkerFetcherError):
        fetch_marker_set(spec, tmp_path / "data")


def test_fetch_unsupported_archive_format(tmp_path: Path):
    """An unknown ``archive_format`` raises a clear error."""
    spec = MarkerSetSpec(
        key="weird",
        label="Weird",
        version="dev",
        url="",
        sha256="",
        expected_marker_count=0,
        archive_format="zip",
        hmm_subdir="",
        local_marker_dir="",
        scores_cutoff_relpath="",
        description="",
    )
    with pytest.raises(MarkerFetcherError):
        fetch_marker_set(spec, tmp_path / "data")


def test_fetch_with_manifest_upserts_row(tmp_path: Path):
    """Passing ``manifest_path`` writes one row into the unified manifest."""
    local_dir = _make_local_marker_dir(tmp_path, n=2)
    spec = _local_spec(local_dir=local_dir, n=2)
    manifest = tmp_path / "manifest.tsv"
    fetch_marker_set(spec, tmp_path / "data", manifest_path=manifest)
    assert manifest.exists()
    body = manifest.read_text(encoding="utf-8")
    assert "Local_Only" in body
    assert "QC_MARKERS" in body
