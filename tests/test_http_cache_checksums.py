"""A cached download must satisfy the same upstream checksum as new bytes."""
from __future__ import annotations

import hashlib
from pathlib import Path
from unittest.mock import patch

import pytest

from src._http import FetchError, stream_download
from src.config_loader import DatabaseAsset
from src.silva_fetcher import fetch_silva_asset


def test_matching_cached_download_is_reused_without_network(tmp_path: Path):
    data = b">reference\nACGT\n"
    dest = tmp_path / "reference.fasta"
    dest.write_bytes(data)
    with patch("src._http.requests.get", side_effect=AssertionError("unexpected download")):
        digest, size = stream_download("https://example.test/reference.fasta", dest,
                                       expected_md5=hashlib.md5(data).hexdigest().upper())
    assert (digest, size) == (hashlib.sha256(data).hexdigest(), len(data))


def test_mismatched_cached_download_is_rejected_and_preserved(tmp_path: Path):
    dest = tmp_path / "reference.fasta"
    dest.write_bytes(b"truncated reference")
    with patch("src._http.requests.get", side_effect=AssertionError("unexpected download")):
        with pytest.raises(FetchError, match="MD5 mismatch"):
            stream_download("https://example.test/reference.fasta", dest,
                            expected_md5=hashlib.md5(b"complete reference").hexdigest())
    assert dest.read_bytes() == b"truncated reference"


def test_cache_without_an_expected_checksum_keeps_existing_behavior(tmp_path: Path):
    data = b"cached bytes"
    dest = tmp_path / "reference.fasta"
    dest.write_bytes(data)
    with patch("src._http.requests.get", side_effect=AssertionError("unexpected download")):
        assert stream_download("https://example.test/reference.fasta", dest) == (
            hashlib.sha256(data).hexdigest(), len(data))


def test_silva_does_not_report_corrupt_cached_bytes_as_success(tmp_path: Path):
    data = b"complete reference"
    dest = tmp_path / "reference.fasta"
    dest.write_bytes(b"corrupt reference")
    asset = DatabaseAsset(label="reference", url="https://example.test/reference.fasta",
                          md5_url="https://example.test/reference.fasta.md5", fmt="fasta")
    with patch("src.silva_fetcher.fetch_md5_checksum", return_value=hashlib.md5(data).hexdigest()), \
         patch("src._http.requests.get", side_effect=AssertionError("unexpected download")):
        with pytest.raises(FetchError, match="MD5 mismatch"):
            fetch_silva_asset(asset, tmp_path)
