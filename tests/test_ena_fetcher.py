"""Tests for ``src/ena_fetcher.py``.

These tests mock the ENA Portal API so the parser is exercised without network access.
The fixture mirrors the verified PRJEB59431 response (37 cyanobacterial MAGs).
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from src import ena_fetcher
from src._http import FetchError
from src.ena_fetcher import ENAAssembly, list_bioproject_assemblies


PRJEB59431_FIXTURE = """accession\tsample_accession\tscientific_name\tassembly_level
GCA_949127635\tSAMEA112467283\tNostocales cyanobacterium PMM_0069\tcontig
GCA_949127645\tSAMEA112467274\tSynechococcales cyanobacterium PMM_0036\tcontig
GCA_949127655\tSAMEA112467292\tChroococcales cyanobacterium PMM_0086\tcontig
GCA_949127665\tSAMEA112467289\tSynechococcales cyanobacterium PMM_0083\tcontig
GCA_949127685\tSAMEA112467282\tCandidatus Sivonenia alaskensis PMM_0068\tcontig
"""


class _MockResponse:
    def __init__(self, text: str, status_code: int = 200):
        self.text = text
        self.status_code = status_code

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            import requests
            raise requests.HTTPError(f"status {self.status_code}")


def test_list_assemblies_parses_tsv():
    with patch.object(ena_fetcher.requests, "get", return_value=_MockResponse(PRJEB59431_FIXTURE)):
        out = list_bioproject_assemblies("PRJEB59431")
    assert len(out) == 5
    assert out[0] == ENAAssembly(
        accession="GCA_949127635",
        sample_accession="SAMEA112467283",
        scientific_name="Nostocales cyanobacterium PMM_0069",
        assembly_level="contig",
    )
    assert out[-1].scientific_name.startswith("Candidatus Sivonenia")


def test_list_assemblies_skips_blank_rows():
    body = "accession\tsample_accession\tscientific_name\tassembly_level\n\nGCA_111111111\tS1\tFoo bar\tcontig\n"
    with patch.object(ena_fetcher.requests, "get", return_value=_MockResponse(body)):
        out = list_bioproject_assemblies("PRJEB99999")
    assert len(out) == 1
    assert out[0].accession == "GCA_111111111"


def test_list_assemblies_http_error():
    with patch.object(ena_fetcher.requests, "get", return_value=_MockResponse("", status_code=503)):
        with pytest.raises(FetchError, match="filereport request failed"):
            list_bioproject_assemblies("PRJEB59431")


def test_fetch_bioproject_records_failures(tmp_path: Path):
    """A failing per-assembly download must NOT abort the batch; it produces a
    failed FetchResult row instead."""
    with patch.object(ena_fetcher.requests, "get", return_value=_MockResponse(PRJEB59431_FIXTURE)):
        def _boom(*_args, **_kwargs):
            raise FetchError("simulated network failure")

        with patch.object(ena_fetcher, "fetch_assembly_fasta", side_effect=_boom):
            results = ena_fetcher.fetch_bioproject("PRJEB59431", tmp_path, limit=3)

    assert len(results) == 3
    assert all(r.status == "failed" for r in results)
    assert all("simulated network failure" in (r.note or "") for r in results)
