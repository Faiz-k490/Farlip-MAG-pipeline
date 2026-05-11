"""ENA (European Nucleotide Archive) BioProject fetcher.

Two responsibilities:

1. **Enumerate** every assembly under a BioProject (e.g. PRJEB59431) using the ENA
   Portal API. One HTTPS GET returns a TSV with one row per assembly.
   Endpoint: ``https://www.ebi.ac.uk/ena/portal/api/filereport``
   Documentation: https://www.ebi.ac.uk/ena/portal/api/

2. **Download** the assembled-contigs FASTA for each ``GCA_*`` assembly via the ENA
   Browser API.
   Endpoint: ``https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}``
   Documentation: https://www.ebi.ac.uk/ena/browser/api/

No authentication is needed for either endpoint. No environment variables are read.

For PRJEB59431 (the CCAMBIO polar cyanobacteria BioProject locked in for iteration 1)
the enumeration call has been verified to return 37 assembly rows.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from pathlib import Path

import requests

from ._http import (
    DEFAULT_TIMEOUT,
    USER_AGENT,
    FetchError,
    FetchResult,
    compute_sha256,
    stream_download,
    utcnow_iso,
)


ENA_PORTAL_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
"""ENA Portal API filereport endpoint. The same URL serves reads, assemblies, sequences,
analyses - selectivity is controlled by the ``result=`` query parameter."""

ENA_BROWSER_FASTA_TEMPLATE = "https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}?download=true&gzip=true"
"""ENA Browser API FASTA endpoint. For a GCA_* accession this returns every contig in the
assembly concatenated into a single gzipped FASTA stream."""

ENA_ASSEMBLY_FIELDS = ("accession", "sample_accession", "scientific_name", "assembly_level")
"""Default fields requested from the assembly file report. ``assembly_level`` is included
so the caller can distinguish chromosome / scaffold / contig MAGs."""


@dataclass(frozen=True)
class ENAAssembly:
    """One row from the ENA assembly file report."""

    accession: str
    sample_accession: str
    scientific_name: str
    assembly_level: str


def list_bioproject_assemblies(
    bioproject: str,
    fields: tuple[str, ...] = ENA_ASSEMBLY_FIELDS,
) -> list[ENAAssembly]:
    """Return every assembly registered under an ENA BioProject.

    The Portal API's ``result=assembly`` view is used (not ``result=read_run``) because
    PRJEB59431 has both raw reads and assembled MAGs - we want only the MAGs.

    Args:
        bioproject: ENA / INSDC BioProject accession (e.g. ``"PRJEB59431"``).
        fields: Comma-joined TSV column list. Defaults cover label, sample, name, level.

    Returns:
        List of :class:`ENAAssembly`. Empty list if the project has no assemblies
        (uncommon, but a possibility for read-only submissions).

    Raises:
        FetchError: On HTTP error or unparseable TSV.
    """
    params = {
        "accession": bioproject,
        "result": "assembly",
        "fields": ",".join(fields),
        "format": "tsv",
    }
    try:
        r = requests.get(
            ENA_PORTAL_FILEREPORT,
            params=params,
            timeout=DEFAULT_TIMEOUT,
            headers={"User-Agent": USER_AGENT, "Accept": "text/tab-separated-values"},
        )
        r.raise_for_status()
    except requests.RequestException as exc:
        raise FetchError(f"ENA Portal filereport request failed for {bioproject}: {exc}") from exc

    reader = csv.DictReader(io.StringIO(r.text), delimiter="\t")
    out: list[ENAAssembly] = []
    for row in reader:
        if not row.get("accession"):
            continue
        out.append(
            ENAAssembly(
                accession=row["accession"],
                sample_accession=row.get("sample_accession", ""),
                scientific_name=row.get("scientific_name", ""),
                assembly_level=row.get("assembly_level", ""),
            )
        )
    return out


def fetch_assembly_fasta(
    assembly: ENAAssembly | str,
    dest_dir: Path,
    *,
    overwrite: bool = False,
) -> FetchResult:
    """Download one assembly's contigs FASTA from the ENA Browser API.

    Args:
        assembly: An :class:`ENAAssembly` or a raw accession string. Strings are wrapped
            into a stub assembly so the function works in both enumerate-then-download
            and one-off-by-accession workflows.
        dest_dir: Directory to write ``<accession>.fasta.gz`` into.
        overwrite: If True, re-download even if the file already exists.

    Returns:
        FetchResult describing the downloaded gzipped FASTA.

    Raises:
        FetchError: On HTTP failure or empty body.
    """
    if isinstance(assembly, str):
        assembly = ENAAssembly(
            accession=assembly, sample_accession="", scientific_name="", assembly_level=""
        )

    dest_dir.mkdir(parents=True, exist_ok=True)
    out_path = dest_dir / f"{assembly.accession}.fasta.gz"
    url = ENA_BROWSER_FASTA_TEMPLATE.format(accession=assembly.accession)

    if not overwrite and out_path.exists() and out_path.stat().st_size > 0:
        sha = compute_sha256(out_path)
        return FetchResult(
            accession=assembly.accession,
            source="ENA",
            url=url,
            local_path=out_path,
            sha256=sha,
            bytes_downloaded=out_path.stat().st_size,
            fetched_at=utcnow_iso(),
            status="skipped",
            note="already present",
            extra={
                "scientific_name": assembly.scientific_name,
                "sample_accession": assembly.sample_accession,
                "assembly_level": assembly.assembly_level,
            },
        )

    sha, n_bytes = stream_download(url, out_path, overwrite=overwrite)
    if n_bytes == 0:
        raise FetchError(f"ENA returned an empty body for {assembly.accession}")

    return FetchResult(
        accession=assembly.accession,
        source="ENA",
        url=url,
        local_path=out_path,
        sha256=sha,
        bytes_downloaded=n_bytes,
        fetched_at=utcnow_iso(),
        status="ok",
        extra={
            "scientific_name": assembly.scientific_name,
            "sample_accession": assembly.sample_accession,
            "assembly_level": assembly.assembly_level,
        },
    )


def fetch_bioproject(
    bioproject: str,
    dest_dir: Path,
    *,
    overwrite: bool = False,
    limit: int | None = None,
) -> list[FetchResult]:
    """Enumerate a BioProject and download every assembly.

    Args:
        bioproject: ENA / INSDC BioProject accession (e.g. ``"PRJEB59431"``).
        dest_dir: Directory under which assemblies are written.
            One file per assembly: ``<dest_dir>/<accession>.fasta.gz``.
        overwrite: Forwarded to :func:`fetch_assembly_fasta`.
        limit: If set, only fetch the first ``limit`` assemblies. Useful for smoke tests
            and the iteration-1 partial-validation gate.

    Returns:
        One FetchResult per assembly attempted, in enumeration order.

    Raises:
        FetchError: If the enumeration call itself fails. Per-assembly failures are
            captured as FetchResult rows with ``status="failed"`` so a single bad MAG
            does not abort the whole batch.
    """
    assemblies = list_bioproject_assemblies(bioproject)
    if limit is not None:
        assemblies = assemblies[:limit]

    results: list[FetchResult] = []
    for asm in assemblies:
        try:
            results.append(fetch_assembly_fasta(asm, dest_dir, overwrite=overwrite))
        except FetchError as exc:
            results.append(
                FetchResult(
                    accession=asm.accession,
                    source="ENA",
                    url=ENA_BROWSER_FASTA_TEMPLATE.format(accession=asm.accession),
                    local_path=dest_dir / f"{asm.accession}.fasta.gz",
                    fetched_at=utcnow_iso(),
                    status="failed",
                    note=str(exc),
                    extra={
                        "scientific_name": asm.scientific_name,
                        "sample_accession": asm.sample_accession,
                        "assembly_level": asm.assembly_level,
                    },
                )
            )
    return results
