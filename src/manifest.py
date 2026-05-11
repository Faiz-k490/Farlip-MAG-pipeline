"""Single-source-of-truth manifest writer for the ingestion pipeline.

Every successful (or failed) fetch produces one row in ``data/manifest.tsv``. Iterations
2+ (HMM search, phylogeny, etc.) read this file to know what is on disk and where.

The manifest is a TSV (not CSV) because some upstream fields - notably
``scientific_name`` - frequently contain commas and quotation marks. Tabs are far less
likely to collide with biological text.

Columns are stable and additive. New optional fields must be appended to the right;
never reordered, never removed without a major version bump.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

from ._http import FetchResult


MANIFEST_COLUMNS: tuple[str, ...] = (
    "accession",
    "source",
    "url",
    "local_path",
    "sha256",
    "bytes_downloaded",
    "fetched_at",
    "status",
    "note",
    "label",
    "scientific_name",
    "sample_accession",
    "assembly_level",
    "format",
    "md5",
    "extracted_dir",
)
"""Canonical manifest column order. Additive only - never remove or reorder."""


def _row_from_result(result: FetchResult) -> dict[str, str]:
    """Flatten a FetchResult (including its ``extra`` dict) into a TSV-friendly mapping."""
    base = {
        "accession": result.accession,
        "source": result.source,
        "url": result.url,
        "local_path": str(result.local_path),
        "sha256": result.sha256,
        "bytes_downloaded": str(result.bytes_downloaded),
        "fetched_at": result.fetched_at,
        "status": result.status,
        "note": result.note or "",
    }
    # ``extra`` carries source-specific metadata (scientific_name from ENA, etc.). We
    # merge it after the base fields so it cannot accidentally overwrite required columns.
    for k, v in result.extra.items():
        if k in base:
            continue
        base[k] = str(v)
    return base


def read_manifest(path: Path) -> list[dict[str, str]]:
    """Read an existing manifest TSV into a list of row dicts.

    Args:
        path: Manifest path. If absent, returns an empty list.

    Returns:
        List of dicts keyed by the headers actually present in the file.
    """
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def write_manifest(path: Path, rows: Iterable[dict[str, str]]) -> None:
    """Write a manifest TSV using the canonical column order.

    Args:
        path: Output path. Parent directories are created if missing.
        rows: Iterable of row dicts. Missing columns are filled with empty strings;
            unknown columns are silently dropped (so future schema extensions in
            FetchResult ``extra`` don't blow up older readers).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            full_row = {col: row.get(col, "") for col in MANIFEST_COLUMNS}
            writer.writerow(full_row)


def upsert_results(path: Path, results: Iterable[FetchResult]) -> int:
    """Merge new FetchResults into the manifest, replacing rows with the same
    (accession, url) key.

    The (accession, url) compound key makes the pipeline idempotent: re-running the
    fetch updates the existing row in place, so the manifest never grows duplicate
    entries for the same artifact.

    Args:
        path: Manifest path.
        results: Iterable of new FetchResults to merge in.

    Returns:
        Number of rows written (existing + new, after dedupe).
    """
    new_rows = [_row_from_result(r) for r in results]
    existing = read_manifest(path)

    keyed: dict[tuple[str, str], dict[str, str]] = {
        (row.get("accession", ""), row.get("url", "")): row for row in existing
    }
    for row in new_rows:
        keyed[(row["accession"], row["url"])] = row

    merged = list(keyed.values())
    write_manifest(path, merged)
    return len(merged)
