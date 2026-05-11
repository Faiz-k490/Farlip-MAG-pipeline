"""SILVA ribosomal RNA database fetcher.

Downloads the SILVA SSU Ref NR99 FASTA (and optionally any other SILVA assets listed in
``config/databases.yaml``). Each asset has a sibling ``.md5`` URL which is fetched and
verified before the download is considered successful.

Source: https://www.arb-silva.de/ - the FASTA "Exports" directory under each release.

No environment variables are read. No authentication is required.
"""

from __future__ import annotations

from pathlib import Path

from ._http import FetchError, FetchResult, fetch_md5_checksum, stream_download, utcnow_iso
from .config_loader import DatabaseAsset, DatabasesConfig


def fetch_silva_asset(
    asset: DatabaseAsset,
    dest_dir: Path,
    *,
    overwrite: bool = False,
) -> FetchResult:
    """Download a single SILVA asset and verify its MD5 sidecar.

    Args:
        asset: A SILVA :class:`DatabaseAsset` from the parsed databases config.
        dest_dir: Directory to write the file into.
            Filename is derived from the asset URL.
        overwrite: If True, re-download even when the file exists and hashes match.

    Returns:
        FetchResult capturing the SHA-256 of the downloaded file and the upstream MD5
        (recorded in ``extra['md5']`` for cross-reference).

    Raises:
        FetchError: On HTTP failure or MD5 mismatch.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    filename = asset.url.rsplit("/", 1)[-1]
    out_path = dest_dir / filename

    expected_md5: str | None = None
    if asset.md5_url:
        try:
            expected_md5 = fetch_md5_checksum(asset.md5_url)
        except FetchError as exc:
            raise FetchError(
                f"Could not retrieve SILVA MD5 sidecar for {asset.label} ({asset.md5_url}): {exc}"
            ) from exc

    sha, n_bytes = stream_download(asset.url, out_path, expected_md5=expected_md5, overwrite=overwrite)

    return FetchResult(
        accession=asset.label,
        source="SILVA",
        url=asset.url,
        local_path=out_path,
        sha256=sha,
        bytes_downloaded=n_bytes,
        fetched_at=utcnow_iso(),
        status="ok",
        extra={"md5": expected_md5 or "", "format": asset.fmt},
    )


def fetch_silva(databases: DatabasesConfig, dest_root: Path, *, overwrite: bool = False) -> list[FetchResult]:
    """Download every SILVA asset declared in the databases config.

    Args:
        databases: Parsed :class:`DatabasesConfig`.
        dest_root: Root output directory. Files are written under
            ``<dest_root>/silva/<version>/``.
        overwrite: Forwarded to :func:`fetch_silva_asset`.

    Returns:
        One FetchResult per asset (in declaration order).

    Raises:
        FetchError: Propagated from any failing asset download.
    """
    dest = dest_root / "silva" / databases.silva_version
    return [fetch_silva_asset(a, dest, overwrite=overwrite) for a in databases.silva_assets]
