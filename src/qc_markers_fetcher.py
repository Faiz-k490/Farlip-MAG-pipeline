"""Marker-HMM fetcher for the Phase 4 quality assessor.

Fetches the bacterial single-copy ortholog (SCO) marker HMM set declared in
[`config/qc_markers.yaml`](../config/qc_markers.yaml), verifies the upstream
checksum, extracts the individual ``.hmm`` files, and concatenates them into
one combined database (``data/qc_markers/<label>.hmm``) that
:func:`pyhmmer.plan7.HMMFile` reads in a single open call.

Two marker-set archive formats are supported:

* ``tar.gz`` (default - BUSCO ``bacteria_odb10``). Streamed download with the
  shared :func:`src._http.stream_download` helper, SHA-256 verified against the
  config, then extracted into a per-marker-set scratch directory under
  ``data/qc_markers/<label>/``.
* ``local_directory`` (manual escape hatch - GTDB R226 bac120 inside a
  pre-extracted GTDB-Tk reference package). The user supplies the path via the
  ``local_marker_dir`` field; we copy nothing - we just glob the directory and
  concatenate.

The module deliberately does **not** attempt to download multi-hundred-GB
archives such as the full GTDB-Tk reference package; entries flagged
``archive_format: local_directory`` raise :class:`MarkerFetcherError` if no
``local_marker_dir`` is configured.

Outputs:

* ``data/qc_markers/<label>.hmm``  - concatenated marker HMM database used by
  :mod:`src.qc_assessor`.
* ``data/qc_markers/<label>.stamp`` - tiny key/value file recording the source
  URL, the SHA-256 of the concatenated HMM file, and the marker count, used by
  the assessor and the Phase 4 manifest writer.
* One row in ``data/manifest.tsv`` keyed on ``(accession=label, url)``.
"""

from __future__ import annotations

import hashlib
import shutil
import tarfile
from dataclasses import dataclass
from pathlib import Path

import yaml

from ._http import (
    FetchError,
    FetchResult,
    compute_sha256,
    stream_download,
    utcnow_iso,
)
from .manifest import upsert_results


QC_MARKERS_DIR_NAME: str = "qc_markers"
"""Subdirectory under ``data/`` that holds extracted marker sets and the
concatenated HMM database the assessor reads."""

MAX_AUTO_DOWNLOAD_BYTES: int = 200 * 1024 * 1024
"""Hard cap on the size of any tarball we will auto-download. The BUSCO
bacteria_odb10 archive is ~6 MB, so 200 MB leaves comfortable headroom for
larger BUSCO lineages while refusing to silently pull a multi-GB GTDB
reference package."""

ARCHIVE_FORMAT_TAR_GZ: str = "tar.gz"
"""Marker-set archive format that is streamed and extracted automatically."""

ARCHIVE_FORMAT_LOCAL_DIR: str = "local_directory"
"""Marker-set archive format that requires a user-supplied local directory."""


class MarkerFetcherError(RuntimeError):
    """Raised when the marker fetcher cannot produce a usable HMM database."""


@dataclass(frozen=True)
class MarkerSetSpec:
    """Validated representation of one entry under ``marker_sets:`` in the YAML.

    Attributes:
        key: The YAML key (e.g. ``"busco_bacteria_odb10"``).
        label: Human-readable label written into the manifest.
        version: Upstream version / release tag.
        url: Download URL (only required for ``tar.gz`` archives).
        sha256: Expected SHA-256 of the downloaded archive (optional but
            strongly recommended).
        expected_marker_count: How many ``.hmm`` files we expect to find. Used
            as the denominator in the completeness calculation.
        archive_format: One of :data:`ARCHIVE_FORMAT_TAR_GZ` or
            :data:`ARCHIVE_FORMAT_LOCAL_DIR`.
        hmm_subdir: Path inside the extracted archive that contains the
            ``.hmm`` files (e.g. ``"bacteria_odb10/hmms"``). Ignored for
            ``local_directory``.
        local_marker_dir: Filesystem path to a directory of ``.hmm`` files;
            required for the ``local_directory`` archive format and ignored
            otherwise.
        scores_cutoff_relpath: Optional path inside the extracted archive
            pointing at a BUSCO-style ``scores_cutoff`` file (one
            ``<hmm_name>\\t<bit_score>`` per line). When present, the assessor
            applies these per-HMM bit-score cutoffs to discriminate true
            single-copy markers from divergent paralogs that an E-value alone
            would let through.
        description: Free-text note copied into the manifest.
    """

    key: str
    label: str
    version: str
    url: str
    sha256: str
    expected_marker_count: int
    archive_format: str
    hmm_subdir: str
    local_marker_dir: str
    scores_cutoff_relpath: str
    description: str


@dataclass(frozen=True)
class MarkerFetchResult:
    """Outcome of one :func:`fetch_marker_set` invocation.

    Attributes:
        label: Marker-set label from the YAML.
        hmm_db_path: Path to the concatenated ``<label>.hmm`` database.
        stamp_path: Path to the sidecar stamp file.
        n_markers: Number of individual ``.hmm`` profiles concatenated into the
            database.
        archive_path: Path to the downloaded archive (or to the user-supplied
            local directory for the ``local_directory`` format).
        archive_sha256: SHA-256 of the downloaded archive, or empty string for
            the ``local_directory`` format.
        scores_cutoff_path: Path to the per-HMM bit-score cutoff TSV (BUSCO
            ``scores_cutoff`` format). ``None`` if the marker set does not
            ship one - the assessor then falls back to the global E-value
            threshold alone.
        skipped: True if the call short-circuited because the database already
            existed and matched the expected source SHA.
    """

    label: str
    hmm_db_path: Path
    stamp_path: Path
    n_markers: int
    archive_path: Path
    archive_sha256: str
    scores_cutoff_path: Path | None
    skipped: bool


def load_marker_config(config_path: Path) -> dict[str, MarkerSetSpec]:
    """Parse ``config/qc_markers.yaml`` into a ``{key: MarkerSetSpec}`` map.

    Args:
        config_path: Path to ``qc_markers.yaml``.

    Returns:
        Mapping of YAML keys to validated :class:`MarkerSetSpec` instances.

    Raises:
        MarkerFetcherError: If the file is missing or malformed.
    """
    if not config_path.exists():
        raise MarkerFetcherError(f"qc_markers config not found: {config_path}")
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    sets = raw.get("marker_sets") or {}
    if not isinstance(sets, dict) or not sets:
        raise MarkerFetcherError(
            f"qc_markers config has no `marker_sets:` block: {config_path}"
        )
    out: dict[str, MarkerSetSpec] = {}
    for key, entry in sets.items():
        if not isinstance(entry, dict):
            raise MarkerFetcherError(f"marker set {key!r} is not a mapping")
        out[key] = MarkerSetSpec(
            key=str(key),
            label=str(entry.get("label") or key),
            version=str(entry.get("version") or ""),
            url=str(entry.get("url") or ""),
            sha256=str(entry.get("sha256") or ""),
            expected_marker_count=int(entry.get("expected_marker_count") or 0),
            archive_format=str(entry.get("archive_format") or ARCHIVE_FORMAT_TAR_GZ),
            hmm_subdir=str(entry.get("hmm_subdir") or ""),
            local_marker_dir=str(entry.get("local_marker_dir") or ""),
            scores_cutoff_relpath=str(entry.get("scores_cutoff_relpath") or ""),
            description=str(entry.get("description") or ""),
        )
    return out


def resolve_default_spec(config_path: Path) -> MarkerSetSpec:
    """Pick the marker-set spec named by ``default_marker_set:`` in the YAML.

    Args:
        config_path: Path to ``qc_markers.yaml``.

    Returns:
        The :class:`MarkerSetSpec` for the configured default marker set.

    Raises:
        MarkerFetcherError: If the YAML does not name a default or names a
            marker set that is not declared.
    """
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    default_key = raw.get("default_marker_set")
    if not default_key:
        raise MarkerFetcherError(
            f"qc_markers config does not set `default_marker_set`: {config_path}"
        )
    specs = load_marker_config(config_path)
    if default_key not in specs:
        raise MarkerFetcherError(
            f"default_marker_set={default_key!r} not declared in marker_sets:"
        )
    return specs[default_key]


def _read_stamp(stamp_path: Path) -> dict[str, str]:
    """Parse a stamp file into a ``{key: value}`` dict (or ``{}`` if missing)."""
    if not stamp_path.exists():
        return {}
    out: dict[str, str] = {}
    for line in stamp_path.read_text(encoding="utf-8").splitlines():
        if "\t" not in line:
            continue
        k, v = line.split("\t", 1)
        out[k] = v
    return out


def _write_stamp(
    stamp_path: Path,
    *,
    label: str,
    archive_sha256: str,
    hmm_db_sha256: str,
    n_markers: int,
    source_url: str,
    scores_cutoff_path: Path | None,
) -> None:
    """Write a small key/value stamp recording the marker DB provenance.

    Args:
        stamp_path: Output path.
        label: Marker-set label.
        archive_sha256: SHA-256 of the downloaded source archive (empty for
            ``local_directory`` marker sets).
        hmm_db_sha256: SHA-256 of the concatenated HMM file we just wrote.
        n_markers: Number of HMM profiles concatenated.
        source_url: URL or local path the markers were sourced from.
        scores_cutoff_path: Optional path to the per-HMM bit-score cutoff
            sidecar; written into the stamp so the assessor can locate it
            without re-parsing the YAML.
    """
    stamp_path.parent.mkdir(parents=True, exist_ok=True)
    payload = [
        f"label\t{label}",
        f"archive_sha256\t{archive_sha256}",
        f"hmm_db_sha256\t{hmm_db_sha256}",
        f"n_markers\t{n_markers}",
        f"source\t{source_url}",
        f"scores_cutoff_path\t{scores_cutoff_path or ''}",
        f"generated_at\t{utcnow_iso()}",
    ]
    stamp_path.write_text("\n".join(payload) + "\n", encoding="utf-8")


def _gather_hmm_files(directory: Path) -> list[Path]:
    """Return a deterministic, sorted list of every ``.hmm`` file under ``directory``."""
    return sorted(p for p in directory.rglob("*.hmm") if p.is_file())


def _concatenate_hmms(hmm_files: list[Path], dest: Path) -> str:
    """Concatenate every HMM file into ``dest`` and return the SHA-256 of the result.

    Each pyhmmer-readable ``.hmm`` file is itself the textual HMMER3 plan7
    format; concatenating them is a valid multi-HMM database. We defensively
    insert a newline between profiles when the previous file does not already
    end in one, so a single missing trailing newline upstream cannot fuse two
    profile records.

    Args:
        hmm_files: Per-profile ``.hmm`` paths in the order to write.
        dest: Output path.

    Returns:
        Hex SHA-256 of the concatenated database.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    hasher = hashlib.sha256()
    with dest.open("wb") as out_fh:
        for f in hmm_files:
            last_byte: bytes = b""
            with f.open("rb") as in_fh:
                for chunk in iter(lambda: in_fh.read(1024 * 1024), b""):
                    out_fh.write(chunk)
                    hasher.update(chunk)
                    last_byte = chunk[-1:] or last_byte
            if last_byte and last_byte != b"\n":
                out_fh.write(b"\n")
                hasher.update(b"\n")
    return hasher.hexdigest()


def _autodetect_scores_cutoff(extract_dir: Path) -> Path | None:
    """Look for a BUSCO-style ``scores_cutoff`` file inside an extracted archive.

    BUSCO lineage tarballs always ship one at ``<lineage>/scores_cutoff``.
    Other marker-set archives may or may not include one; we just walk the
    extraction tree and pick the first match, which is unambiguous in
    practice (BUSCO directories contain exactly one).

    Args:
        extract_dir: Root of the extracted archive.

    Returns:
        Path to the first ``scores_cutoff`` (or ``*.scores_cutoff``) file
        found, or ``None`` if none exists.
    """
    candidates: list[Path] = []
    for name in ("scores_cutoff",):
        candidates.extend(p for p in extract_dir.rglob(name) if p.is_file())
    if candidates:
        return sorted(candidates)[0]
    return None


def _copy_scores_cutoff(src: Path, dest: Path) -> Path:
    """Copy a per-HMM bit-score cutoff file into the qc_markers directory.

    Args:
        src: Source path inside the extracted archive.
        dest: Destination path (typically
            ``data/qc_markers/<key>.scores_cutoff.tsv``).

    Returns:
        The destination path.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(src.read_bytes())
    return dest


def _safe_extract_tarball(archive_path: Path, dest_dir: Path) -> None:
    """Extract a ``.tar.gz`` archive guarding against path-traversal members.

    Args:
        archive_path: Path to the tarball.
        dest_dir: Directory to extract into. Created if missing.

    Raises:
        MarkerFetcherError: If a member path would escape ``dest_dir``.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, mode="r:gz") as tf:
        for member in tf.getmembers():
            target = (dest_dir / member.name).resolve()
            if not str(target).startswith(str(dest_dir.resolve())):
                raise MarkerFetcherError(
                    f"refusing to extract path-traversal entry {member.name!r}"
                )
        # Python 3.12+ allows declaring the extraction filter explicitly.
        # ``data`` rejects absolute paths, parent traversal, and special
        # files; the loop above is kept as a belt-and-braces defence.
        tf.extractall(dest_dir, filter="data")


def fetch_marker_set(
    spec: MarkerSetSpec,
    data_root: Path,
    *,
    overwrite: bool = False,
    manifest_path: Path | None = None,
) -> MarkerFetchResult:
    """Download (if needed), extract, and concatenate one marker set.

    Args:
        spec: Validated marker set spec from :func:`load_marker_config`.
        data_root: Root data directory. Outputs land under
            ``<data_root>/qc_markers/``.
        overwrite: If True, re-download / re-extract even when the stamp already
            records the same source SHA.
        manifest_path: Optional path to ``data/manifest.tsv``. When supplied the
            fetcher upserts a row keyed on ``(label, url-or-localdir)`` so the
            marker set is discoverable from the unified ledger.

    Returns:
        A :class:`MarkerFetchResult` describing the artifacts on disk.

    Raises:
        MarkerFetcherError: If the archive format is unsupported, the upstream
            archive exceeds :data:`MAX_AUTO_DOWNLOAD_BYTES`, the SHA-256 does
            not match, or no ``.hmm`` files are found after extraction.
        FetchError: On HTTP failure during download.
    """
    qc_root = data_root / QC_MARKERS_DIR_NAME
    qc_root.mkdir(parents=True, exist_ok=True)
    hmm_db_path = qc_root / f"{spec.key}.hmm"
    stamp_path = qc_root / f"{spec.key}.stamp"
    scores_cutoff_dest = qc_root / f"{spec.key}.scores_cutoff.tsv"
    extract_dir = qc_root / spec.key

    if spec.archive_format == ARCHIVE_FORMAT_LOCAL_DIR:
        if not spec.local_marker_dir:
            raise MarkerFetcherError(
                f"marker set {spec.key!r} declares archive_format={ARCHIVE_FORMAT_LOCAL_DIR!r} "
                f"but no `local_marker_dir` is set in the config"
            )
        local_dir = Path(spec.local_marker_dir)
        if not local_dir.is_dir():
            raise MarkerFetcherError(
                f"local_marker_dir for {spec.key!r} is not a directory: {local_dir}"
            )
        hmm_files = _gather_hmm_files(local_dir)
        if not hmm_files:
            raise MarkerFetcherError(
                f"no .hmm files found under local_marker_dir={local_dir}"
            )
        scores_cutoff_path = _resolve_local_scores_cutoff(local_dir, spec)
        prior = _read_stamp(stamp_path) if not overwrite else {}
        if (
            prior.get("source") == str(local_dir)
            and prior.get("n_markers") == str(len(hmm_files))
            and hmm_db_path.exists()
        ):
            return MarkerFetchResult(
                label=spec.label,
                hmm_db_path=hmm_db_path,
                stamp_path=stamp_path,
                n_markers=len(hmm_files),
                archive_path=local_dir,
                archive_sha256="",
                scores_cutoff_path=Path(prior["scores_cutoff_path"])
                if prior.get("scores_cutoff_path")
                else None,
                skipped=True,
            )
        hmm_db_sha = _concatenate_hmms(hmm_files, hmm_db_path)
        scores_cutoff_final = (
            _copy_scores_cutoff(scores_cutoff_path, scores_cutoff_dest)
            if scores_cutoff_path is not None
            else None
        )
        _write_stamp(
            stamp_path,
            label=spec.label,
            archive_sha256="",
            hmm_db_sha256=hmm_db_sha,
            n_markers=len(hmm_files),
            source_url=str(local_dir),
            scores_cutoff_path=scores_cutoff_final,
        )
        if manifest_path is not None:
            _upsert_manifest_row(
                manifest_path,
                spec=spec,
                source_path=str(local_dir),
                hmm_db_path=hmm_db_path,
                archive_sha256="",
                hmm_db_sha256=hmm_db_sha,
                n_markers=len(hmm_files),
            )
        return MarkerFetchResult(
            label=spec.label,
            hmm_db_path=hmm_db_path,
            stamp_path=stamp_path,
            n_markers=len(hmm_files),
            archive_path=local_dir,
            archive_sha256="",
            scores_cutoff_path=scores_cutoff_final,
            skipped=False,
        )

    if spec.archive_format != ARCHIVE_FORMAT_TAR_GZ:
        raise MarkerFetcherError(
            f"unsupported archive_format for {spec.key!r}: {spec.archive_format!r}"
        )
    if not spec.url:
        raise MarkerFetcherError(
            f"marker set {spec.key!r} declares archive_format={ARCHIVE_FORMAT_TAR_GZ!r} "
            f"but no `url` is set"
        )

    archive_path = qc_root / f"{spec.key}.tar.gz"

    prior = _read_stamp(stamp_path) if not overwrite else {}
    if (
        not overwrite
        and hmm_db_path.exists()
        and prior.get("archive_sha256")
        and (not spec.sha256 or prior["archive_sha256"].lower() == spec.sha256.lower())
    ):
        n_markers = int(prior.get("n_markers", "0") or 0)
        return MarkerFetchResult(
            label=spec.label,
            hmm_db_path=hmm_db_path,
            stamp_path=stamp_path,
            n_markers=n_markers,
            archive_path=archive_path,
            archive_sha256=prior["archive_sha256"],
            scores_cutoff_path=Path(prior["scores_cutoff_path"])
            if prior.get("scores_cutoff_path")
            else None,
            skipped=True,
        )

    archive_sha, n_bytes = stream_download(
        spec.url,
        archive_path,
        overwrite=overwrite,
        show_progress=True,
    )
    if n_bytes > MAX_AUTO_DOWNLOAD_BYTES:
        archive_path.unlink(missing_ok=True)
        raise MarkerFetcherError(
            f"refusing to use {spec.key!r} archive: {n_bytes} bytes exceeds "
            f"MAX_AUTO_DOWNLOAD_BYTES={MAX_AUTO_DOWNLOAD_BYTES}"
        )
    if spec.sha256 and archive_sha.lower() != spec.sha256.lower():
        raise MarkerFetcherError(
            f"SHA-256 mismatch for {spec.key!r}: expected {spec.sha256}, got {archive_sha}"
        )

    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    _safe_extract_tarball(archive_path, extract_dir)

    search_root = extract_dir / spec.hmm_subdir if spec.hmm_subdir else extract_dir
    if not search_root.exists():
        # Some BUSCO mirrors sometimes ship without a top-level dir; fall back
        # to recursive search across the whole extraction.
        search_root = extract_dir
    hmm_files = _gather_hmm_files(search_root)
    if not hmm_files:
        raise MarkerFetcherError(
            f"no .hmm files found in extracted archive: {search_root}"
        )

    hmm_db_sha = _concatenate_hmms(hmm_files, hmm_db_path)

    scores_cutoff_src: Path | None = None
    if spec.scores_cutoff_relpath:
        candidate = extract_dir / spec.scores_cutoff_relpath
        if candidate.is_file():
            scores_cutoff_src = candidate
    if scores_cutoff_src is None:
        scores_cutoff_src = _autodetect_scores_cutoff(extract_dir)
    scores_cutoff_final = (
        _copy_scores_cutoff(scores_cutoff_src, scores_cutoff_dest)
        if scores_cutoff_src is not None
        else None
    )

    _write_stamp(
        stamp_path,
        label=spec.label,
        archive_sha256=archive_sha,
        hmm_db_sha256=hmm_db_sha,
        n_markers=len(hmm_files),
        source_url=spec.url,
        scores_cutoff_path=scores_cutoff_final,
    )

    if manifest_path is not None:
        _upsert_manifest_row(
            manifest_path,
            spec=spec,
            source_path=spec.url,
            hmm_db_path=hmm_db_path,
            archive_sha256=archive_sha,
            hmm_db_sha256=hmm_db_sha,
            n_markers=len(hmm_files),
        )

    return MarkerFetchResult(
        label=spec.label,
        hmm_db_path=hmm_db_path,
        stamp_path=stamp_path,
        n_markers=len(hmm_files),
        archive_path=archive_path,
        archive_sha256=archive_sha,
        scores_cutoff_path=scores_cutoff_final,
        skipped=False,
    )


def _resolve_local_scores_cutoff(local_dir: Path, spec: MarkerSetSpec) -> Path | None:
    """Find a ``scores_cutoff`` file when using a ``local_directory`` marker set.

    Args:
        local_dir: The user-supplied marker directory.
        spec: The marker-set spec (used for an optional explicit
            ``scores_cutoff_relpath`` override).

    Returns:
        Path to the cutoff file, or ``None`` if no usable file is found.
    """
    if spec.scores_cutoff_relpath:
        explicit = local_dir / spec.scores_cutoff_relpath
        if explicit.is_file():
            return explicit
    direct = local_dir / "scores_cutoff"
    if direct.is_file():
        return direct
    return _autodetect_scores_cutoff(local_dir)


def _upsert_manifest_row(
    manifest_path: Path,
    *,
    spec: MarkerSetSpec,
    source_path: str,
    hmm_db_path: Path,
    archive_sha256: str,
    hmm_db_sha256: str,
    n_markers: int,
) -> None:
    """Record the marker-set artifact in the unified Phase 1 manifest.

    The accession is the marker-set label so the row is easy to find with
    ``grep`` or pandas; the URL field doubles as the local-path fallback for
    ``local_directory`` marker sets.

    Args:
        manifest_path: Path to ``data/manifest.tsv``.
        spec: The validated marker set spec.
        source_path: Either the upstream URL or the local directory path.
        hmm_db_path: Path to the concatenated HMM database.
        archive_sha256: SHA-256 of the downloaded archive (empty for local
            directories).
        hmm_db_sha256: SHA-256 of the concatenated HMM file.
        n_markers: Number of HMM profiles concatenated.
    """
    bytes_on_disk = hmm_db_path.stat().st_size if hmm_db_path.exists() else 0
    note_parts = [
        f"n_markers={n_markers}",
        f"hmm_db_sha256={hmm_db_sha256}",
    ]
    if archive_sha256:
        note_parts.append(f"archive_sha256={archive_sha256}")
    result = FetchResult(
        accession=spec.label,
        source="QC_MARKERS",
        url=source_path,
        local_path=hmm_db_path,
        sha256=hmm_db_sha256,
        bytes_downloaded=bytes_on_disk,
        fetched_at=utcnow_iso(),
        status="ok",
        note="; ".join(note_parts),
        extra={
            "label": spec.label,
            "format": spec.archive_format,
            "extracted_dir": str(hmm_db_path.parent / spec.key),
        },
    )
    upsert_results(manifest_path, [result])


def fetch_default_marker_set(
    config_path: Path,
    data_root: Path,
    *,
    overwrite: bool = False,
    manifest_path: Path | None = None,
) -> MarkerFetchResult:
    """Convenience wrapper for the most common Phase 3 entry point.

    Resolves the YAML's ``default_marker_set`` and fetches it.

    Args:
        config_path: Path to ``config/qc_markers.yaml``.
        data_root: Root data directory.
        overwrite: Forwarded to :func:`fetch_marker_set`.
        manifest_path: Forwarded to :func:`fetch_marker_set`.

    Returns:
        The :class:`MarkerFetchResult` from the underlying fetch.
    """
    spec = resolve_default_spec(config_path)
    return fetch_marker_set(
        spec,
        data_root,
        overwrite=overwrite,
        manifest_path=manifest_path,
    )


def hmm_db_sha256(hmm_db_path: Path) -> str:
    """Compute the SHA-256 of an existing concatenated HMM database.

    Convenience helper used by :mod:`src.qc_assessor` to construct stamp keys.

    Args:
        hmm_db_path: Path to the concatenated HMM file.

    Returns:
        Hex digest, or an empty string if the file is missing.
    """
    if not hmm_db_path.exists():
        return ""
    return compute_sha256(hmm_db_path)


__all__ = [
    "ARCHIVE_FORMAT_LOCAL_DIR",
    "ARCHIVE_FORMAT_TAR_GZ",
    "MAX_AUTO_DOWNLOAD_BYTES",
    "MarkerFetchResult",
    "MarkerFetcherError",
    "MarkerSetSpec",
    "QC_MARKERS_DIR_NAME",
    "FetchError",
    "fetch_default_marker_set",
    "fetch_marker_set",
    "hmm_db_sha256",
    "load_marker_config",
    "resolve_default_spec",
]
