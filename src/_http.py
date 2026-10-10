"""Shared HTTP, checksum, and FetchResult primitives used by every fetcher.

This module is intentionally private (leading underscore). It exposes one dataclass
(:class:`FetchResult`) and a handful of small utilities that all four fetcher modules
(NCBI Entrez, ENA Portal, SILVA, GTDB) rely on:

* :func:`stream_download` - resumable streaming HTTPS GET with SHA-256 hashing.
* :func:`compute_sha256` - hash an already-downloaded file in chunks.
* :func:`fetch_md5_checksum` - parse the single hex digest from a ``.md5`` URL.
* :func:`utcnow_iso` - timestamp string used in every manifest row.

No upstream provider is contacted unless one of the functions in this module is called.
"""

from __future__ import annotations

import hashlib
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import requests
from tqdm import tqdm


DEFAULT_TIMEOUT = (30, 120)
"""(connect, read) timeout pair in seconds. Read timeout is generous because SILVA and
GTDB downloads can stall briefly mid-stream on the mirror."""

DEFAULT_CHUNK_SIZE = 1024 * 1024  # 1 MiB
"""Streaming download chunk size. Balances progress-bar granularity against syscall cost."""

USER_AGENT = "farlip-ingestion/0.1 (BARTON LAB CLASS; +https://github.com/)"
"""Sent with every outbound request so upstream maintainers can identify the client."""


class FetchError(RuntimeError):
    """Raised when a fetch step fails (network, HTTP, checksum, or parsing)."""


@dataclass
class FetchResult:
    """Outcome of a single fetch operation, written verbatim into ``data/manifest.tsv``.

    Attributes:
        accession: The accession or label identifying the fetched object.
        source: Upstream provider name (``"NCBI"``, ``"ENA"``, ``"SILVA"``, ``"GTDB"``).
        url: The exact URL that was hit.
        local_path: Filesystem location of the downloaded artifact.
        sha256: Hex-encoded SHA-256 of the downloaded bytes.
        bytes_downloaded: Size of the artifact on disk.
        fetched_at: ISO-8601 UTC timestamp recorded after a successful write.
        status: One of ``"ok"``, ``"skipped"``, ``"failed"``.
        note: Optional free-text note (e.g. "skipped: sha256 match" or error detail).
    """

    accession: str
    source: str
    url: str
    local_path: Path
    sha256: str = ""
    bytes_downloaded: int = 0
    fetched_at: str = ""
    status: str = "ok"
    note: str | None = None
    extra: dict[str, str] = field(default_factory=dict)


def utcnow_iso() -> str:
    """Return the current UTC time as an ISO-8601 string with second precision."""
    return datetime.now(tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def compute_sha256(path: Path, chunk_size: int = DEFAULT_CHUNK_SIZE) -> str:
    """Compute the SHA-256 of an existing file by streaming it in chunks."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def stream_download(
    url: str,
    dest: Path,
    *,
    expected_md5: str | None = None,
    overwrite: bool = False,
    show_progress: bool = True,
    headers: dict[str, str] | None = None,
    max_retries: int = 3,
    retry_backoff: float = 2.0,
) -> tuple[str, int]:
    """Stream a URL to ``dest`` with retries, progress, and SHA-256 hashing.

    The download is atomic in the success path: bytes go to ``dest.with_suffix(dest.suffix + '.part')``
    and are renamed to ``dest`` only after a clean close. If ``expected_md5`` is provided,
    it is verified before the rename and a mismatch raises :class:`FetchError`.

    Args:
        url: HTTPS URL to fetch.
        dest: Final on-disk path. Parent directories are created if missing.
        expected_md5: Optional hex MD5 to verify against. Mismatch raises FetchError.
        overwrite: If False and ``dest`` already exists with a non-zero size, hash
            and verify the existing file. A checksum mismatch raises FetchError;
            use overwrite=True to replace the cached file.
        show_progress: Display a tqdm progress bar.
        headers: Extra HTTP headers to merge into the request.
        max_retries: How many times to retry transient errors before giving up.
        retry_backoff: Sleep this many seconds * attempt number between retries.

    Returns:
        Tuple of ``(sha256_hex, bytes_downloaded)``.

    Raises:
        FetchError: On HTTP errors after retries, or on MD5 mismatch.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not overwrite and dest.exists() and dest.stat().st_size > 0:
        if expected_md5 is None:
            return compute_sha256(dest), dest.stat().st_size
        sha = hashlib.sha256()
        md5 = hashlib.md5()
        bytes_read = 0
        with dest.open("rb") as fh:
            for chunk in iter(lambda: fh.read(DEFAULT_CHUNK_SIZE), b""):
                sha.update(chunk)
                md5.update(chunk)
                bytes_read += len(chunk)
        got, want = md5.hexdigest().lower(), expected_md5.lower()
        if got != want:
            raise FetchError(
                f"MD5 mismatch for cached file {dest}: expected {want}, got {got}; "
                "use overwrite=True to download it again"
            )
        return sha.hexdigest(), bytes_read

    part = dest.with_name(dest.name + ".part")
    request_headers = {"User-Agent": USER_AGENT}
    if headers:
        request_headers.update(headers)

    last_error: Exception | None = None
    for attempt in range(1, max_retries + 1):
        try:
            with requests.get(url, stream=True, headers=request_headers, timeout=DEFAULT_TIMEOUT) as r:
                r.raise_for_status()
                total = int(r.headers.get("Content-Length", "0")) or None
                sha = hashlib.sha256()
                md5 = hashlib.md5() if expected_md5 else None
                bytes_written = 0
                # tqdm context guards against ZeroDivisionError when total is unknown.
                bar = tqdm(
                    total=total,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    desc=dest.name,
                    disable=not show_progress,
                    leave=False,
                )
                with bar, part.open("wb") as fh:
                    for chunk in r.iter_content(chunk_size=DEFAULT_CHUNK_SIZE):
                        if not chunk:
                            continue
                        fh.write(chunk)
                        sha.update(chunk)
                        if md5 is not None:
                            md5.update(chunk)
                        bytes_written += len(chunk)
                        bar.update(len(chunk))

            if md5 is not None and expected_md5 is not None:
                got = md5.hexdigest().lower()
                want = expected_md5.lower()
                if got != want:
                    part.unlink(missing_ok=True)
                    raise FetchError(
                        f"MD5 mismatch for {url}: expected {want}, got {got}"
                    )

            part.replace(dest)
            return sha.hexdigest(), bytes_written

        except (requests.RequestException, FetchError) as exc:
            last_error = exc
            part.unlink(missing_ok=True)
            if attempt < max_retries:
                time.sleep(retry_backoff * attempt)
                continue
            raise FetchError(f"Failed to download {url} after {max_retries} attempts: {exc}") from exc

    raise FetchError(f"Unreachable retry exit for {url}: {last_error}")


def fetch_md5_checksum(md5_url: str) -> str:
    """Fetch a SILVA/EBI-style ``.md5`` file and return only the hex digest.

    These files typically contain ``<hexdigest>  <filename>`` on one line; we ignore
    the filename. If the body shape is unexpected we raise rather than guess.

    Args:
        md5_url: URL of the .md5 sidecar file.

    Returns:
        Lowercased 32-character MD5 hex digest.

    Raises:
        FetchError: On HTTP failure or if the file does not contain a hex digest.
    """
    try:
        r = requests.get(md5_url, timeout=DEFAULT_TIMEOUT, headers={"User-Agent": USER_AGENT})
        r.raise_for_status()
    except requests.RequestException as exc:
        raise FetchError(f"Could not fetch MD5 sidecar {md5_url}: {exc}") from exc

    token = r.text.strip().split()[0] if r.text.strip() else ""
    if len(token) != 32 or not all(c in "0123456789abcdefABCDEF" for c in token):
        raise FetchError(f"MD5 sidecar at {md5_url} did not yield a hex digest: {r.text!r}")
    return token.lower()


def iter_with_user_agent(urls: Iterable[str]) -> Iterable[tuple[str, requests.Response]]:
    """Yield ``(url, response)`` pairs for a sequence of URLs, raising on HTTP errors."""
    for url in urls:
        r = requests.get(url, timeout=DEFAULT_TIMEOUT, headers={"User-Agent": USER_AGENT})
        r.raise_for_status()
        yield url, r
