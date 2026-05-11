"""GTDB (Genome Taxonomy Database) fetcher.

Downloads the bacterial taxonomy and metadata files for the locked-in GTDB release.
GTDB does not publish per-file MD5 sidecars on the main download path, so we record only
the locally-computed SHA-256 in the manifest. The release number itself (e.g. ``R226``)
is the version anchor.

Source: https://gtdb.ecogenomic.org/downloads
No authentication is required.
"""

from __future__ import annotations

from pathlib import Path

from ._http import FetchResult, stream_download, utcnow_iso
from .config_loader import DatabaseAsset, DatabasesConfig


def fetch_gtdb_asset(
    asset: DatabaseAsset,
    dest_dir: Path,
    *,
    overwrite: bool = False,
) -> FetchResult:
    """Download a single GTDB asset and record its SHA-256.

    Args:
        asset: A GTDB :class:`DatabaseAsset` from the parsed databases config.
        dest_dir: Directory to write the file into.
        overwrite: If True, re-download even when the file exists.

    Returns:
        FetchResult with the local SHA-256 recorded. No upstream MD5 is verified
        because GTDB does not publish one alongside these files.

    Raises:
        FetchError: On HTTP failure.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    filename = asset.url.rsplit("/", 1)[-1]
    out_path = dest_dir / filename

    sha, n_bytes = stream_download(asset.url, out_path, overwrite=overwrite)
    return FetchResult(
        accession=asset.label,
        source="GTDB",
        url=asset.url,
        local_path=out_path,
        sha256=sha,
        bytes_downloaded=n_bytes,
        fetched_at=utcnow_iso(),
        status="ok",
        extra={"format": asset.fmt},
    )


def fetch_gtdb(databases: DatabasesConfig, dest_root: Path, *, overwrite: bool = False) -> list[FetchResult]:
    """Download every GTDB asset declared in the databases config.

    Args:
        databases: Parsed :class:`DatabasesConfig`.
        dest_root: Root output directory. Files land under ``<dest_root>/gtdb/<version>/``.
        overwrite: Forwarded to :func:`fetch_gtdb_asset`.

    Returns:
        One FetchResult per asset (in declaration order).
    """
    dest = dest_root / "gtdb" / databases.gtdb_version
    return [fetch_gtdb_asset(a, dest, overwrite=overwrite) for a in databases.gtdb_assets]
