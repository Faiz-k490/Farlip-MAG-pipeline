"""NCBI assembly fetcher.

Talks to two NCBI services, in priority order:

1. **NCBI Datasets v2** (preferred): one HTTPS call returns a reproducible ZIP bundle
   containing the genomic FASTA, GFF3 annotation, and a JSON assembly report.
   Endpoint: ``https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download``
   Documentation: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/

2. **Biopython Entrez** (fallback): when only raw GenBank flat-files are needed and
   the Datasets bundle is unavailable. Used as a fallback so the pipeline keeps
   working if the Datasets service has an outage.
   Documentation: https://biopython.org/docs/latest/api/Bio.Entrez.html

Environment variables read here:

* ``NCBI_EMAIL`` *(recommended)* - NCBI's E-utility usage policy requires an email so
  they can contact the user before hitting any rate limit. Defaults to a placeholder
  with a stderr warning.
* ``NCBI_API_KEY`` *(optional)* - raises the per-second request budget from 3 to 10.
"""

from __future__ import annotations

import os
import sys
import zipfile
from pathlib import Path

from Bio import Entrez

from ._http import (
    FetchError,
    FetchResult,
    compute_sha256,
    stream_download,
    utcnow_iso,
)


NCBI_DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
"""Base URL for NCBI Datasets v2 REST API. See:
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/rest-api/"""

DATASETS_DOWNLOAD_TEMPLATE = (
    NCBI_DATASETS_BASE
    + "/genome/accession/{accession}/download"
    + "?include_annotation_type=GENOME_FASTA"
    + "&include_annotation_type=GENOME_GFF"
    + "&include_annotation_type=PROT_FASTA"
    + "&hydrated=FULLY_HYDRATED"
)
"""Datasets v2 download URL. Requests genome FASTA + GFF3 + protein FASTA so iterations
2-3 (HMMER search, MSA) have everything they need without re-fetching."""

_ENV_EMAIL = "NCBI_EMAIL"
_ENV_API_KEY = "NCBI_API_KEY"


def _configure_entrez() -> None:
    """Apply NCBI_EMAIL and NCBI_API_KEY to the Biopython Entrez module-level state."""
    email = os.environ.get(_ENV_EMAIL)
    if not email:
        sys.stderr.write(
            "[entrez_fetcher] WARNING: NCBI_EMAIL is not set. NCBI requires a contact "
            "email; set the env var before running large batches.\n"
        )
        email = "anonymous@example.invalid"
    Entrez.email = email
    api_key = os.environ.get(_ENV_API_KEY)
    if api_key:
        Entrez.api_key = api_key


def fetch_assembly(
    accession: str,
    dest_dir: Path,
    *,
    label: str | None = None,
    overwrite: bool = False,
) -> FetchResult:
    """Download an NCBI assembly bundle via the Datasets v2 service.

    The returned bundle is a ZIP archive. We download it to ``<dest_dir>/<accession>.zip``,
    SHA-256 hash the archive bytes, then extract into ``<dest_dir>/<accession>/`` so
    downstream tools can read the unpacked FASTA/GFF directly.

    Args:
        accession: RefSeq or GenBank assembly accession (e.g. ``"GCF_000733415.1"``).
        dest_dir: Directory under which to place the ZIP and the extracted folder.
        label: Optional human-readable label written into the FetchResult.
        overwrite: If True, re-download even if the ZIP and extracted dir already exist.

    Returns:
        FetchResult describing the ZIP archive (the sha256 of the ZIP, not the FASTA).

    Raises:
        FetchError: On HTTP failure, malformed ZIP, or empty download.
    """
    _configure_entrez()
    dest_dir.mkdir(parents=True, exist_ok=True)
    zip_path = dest_dir / f"{accession}.zip"
    extracted_dir = dest_dir / accession
    url = DATASETS_DOWNLOAD_TEMPLATE.format(accession=accession)

    if not overwrite and zip_path.exists() and extracted_dir.exists() and any(extracted_dir.iterdir()):
        sha = compute_sha256(zip_path)
        return FetchResult(
            accession=accession,
            source="NCBI",
            url=url,
            local_path=zip_path,
            sha256=sha,
            bytes_downloaded=zip_path.stat().st_size,
            fetched_at=utcnow_iso(),
            status="skipped",
            note="already present; sha256 recomputed",
            extra={"label": label or "", "extracted_dir": str(extracted_dir)},
        )

    sha, n_bytes = stream_download(url, zip_path, overwrite=overwrite)

    try:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(extracted_dir)
    except zipfile.BadZipFile as exc:
        raise FetchError(f"Downloaded archive at {zip_path} is not a valid ZIP: {exc}") from exc

    return FetchResult(
        accession=accession,
        source="NCBI",
        url=url,
        local_path=zip_path,
        sha256=sha,
        bytes_downloaded=n_bytes,
        fetched_at=utcnow_iso(),
        status="ok",
        extra={"label": label or "", "extracted_dir": str(extracted_dir)},
    )


def fetch_genbank_record(
    accession: str,
    dest_dir: Path,
    *,
    rettype: str = "gbwithparts",
    retmode: str = "text",
    overwrite: bool = False,
) -> FetchResult:
    """Fetch a raw GenBank flat-file via Biopython Entrez (fallback path).

    Used when only the annotated GenBank record is required (e.g. for a quick gene
    presence check) or when the Datasets v2 service is unreachable.

    Args:
        accession: Any nucleotide accession addressable by ``efetch(db='nuccore')``
            (e.g. ``"CP003610"``).
        dest_dir: Directory to write the ``.gb`` file into.
        rettype: Entrez rettype; defaults to ``gbwithparts`` for complete annotation.
        retmode: Entrez retmode; defaults to ``text``.
        overwrite: If True, re-fetch even if the file already exists.

    Returns:
        FetchResult pointing at the local GenBank file.

    Raises:
        FetchError: On Entrez/HTTP failure or empty response.
    """
    _configure_entrez()
    dest_dir.mkdir(parents=True, exist_ok=True)
    out_path = dest_dir / f"{accession}.gb"
    if not overwrite and out_path.exists() and out_path.stat().st_size > 0:
        sha = compute_sha256(out_path)
        return FetchResult(
            accession=accession,
            source="NCBI",
            url=f"entrez://nuccore/{accession}?rettype={rettype}",
            local_path=out_path,
            sha256=sha,
            bytes_downloaded=out_path.stat().st_size,
            fetched_at=utcnow_iso(),
            status="skipped",
            note="already present",
        )

    try:
        handle = Entrez.efetch(db="nuccore", id=accession, rettype=rettype, retmode=retmode)
        body = handle.read()
        handle.close()
    except Exception as exc:
        raise FetchError(f"Entrez efetch failed for {accession}: {exc}") from exc

    if isinstance(body, bytes):
        out_path.write_bytes(body)
    else:
        out_path.write_text(body, encoding="utf-8")

    if out_path.stat().st_size == 0:
        raise FetchError(f"Entrez efetch returned an empty body for {accession}")

    return FetchResult(
        accession=accession,
        source="NCBI",
        url=f"entrez://nuccore/{accession}?rettype={rettype}",
        local_path=out_path,
        sha256=compute_sha256(out_path),
        bytes_downloaded=out_path.stat().st_size,
        fetched_at=utcnow_iso(),
        status="ok",
    )
