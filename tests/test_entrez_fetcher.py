"""Tests for ``src/entrez_fetcher.py``.

We mock both the streaming HTTP layer and the Biopython Entrez handle so the file
behavior is fully exercised without contacting NCBI.
"""

from __future__ import annotations

import io
import zipfile
from pathlib import Path
from unittest.mock import MagicMock, patch

from src import entrez_fetcher
from src._http import FetchResult


def _make_zip_bytes(contents: dict[str, str]) -> bytes:
    """Build an in-memory ZIP archive matching the NCBI Datasets layout."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for name, body in contents.items():
            zf.writestr(name, body)
    return buf.getvalue()


def test_fetch_assembly_writes_and_extracts(tmp_path: Path):
    zip_bytes = _make_zip_bytes(
        {
            "ncbi_dataset/data/GCF_000733415.1/GCF_000733415.1_genomic.fna": ">contig1\nACGT\n",
            "ncbi_dataset/data/assembly_data_report.jsonl": '{"accession": "GCF_000733415.1"}\n',
        }
    )

    def _fake_stream_download(url, dest, *, expected_md5=None, overwrite=False, **kwargs):
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(zip_bytes)
        import hashlib

        return hashlib.sha256(zip_bytes).hexdigest(), len(zip_bytes)

    with patch.object(entrez_fetcher, "stream_download", side_effect=_fake_stream_download):
        result = entrez_fetcher.fetch_assembly("GCF_000733415.1", tmp_path, label="JSC-1")

    assert result.status == "ok"
    assert result.source == "NCBI"
    assert result.accession == "GCF_000733415.1"
    assert result.local_path.exists()
    assert result.bytes_downloaded == len(zip_bytes)
    extracted = tmp_path / "GCF_000733415.1"
    assert extracted.exists()
    assert (extracted / "ncbi_dataset" / "data" / "GCF_000733415.1" / "GCF_000733415.1_genomic.fna").exists()


def test_fetch_assembly_skips_when_present(tmp_path: Path):
    extracted = tmp_path / "GCF_000733415.1"
    extracted.mkdir()
    (extracted / "placeholder.txt").write_text("present")
    zip_path = tmp_path / "GCF_000733415.1.zip"
    zip_path.write_bytes(_make_zip_bytes({"x": "y"}))

    with patch.object(entrez_fetcher, "stream_download") as mock_dl:
        result = entrez_fetcher.fetch_assembly("GCF_000733415.1", tmp_path)
        mock_dl.assert_not_called()

    assert result.status == "skipped"
    assert "already present" in (result.note or "")


def test_fetch_genbank_record_writes_text(tmp_path: Path):
    fake_handle = MagicMock()
    fake_handle.read.return_value = "LOCUS       FAKE       100 bp    DNA\n//\n"
    fake_handle.close.return_value = None

    with patch.object(entrez_fetcher, "Entrez") as mock_entrez:
        mock_entrez.efetch.return_value = fake_handle
        result = entrez_fetcher.fetch_genbank_record("CP003610", tmp_path)

    assert result.status == "ok"
    assert isinstance(result, FetchResult)
    assert result.local_path.read_text().startswith("LOCUS")
    assert result.bytes_downloaded > 0
