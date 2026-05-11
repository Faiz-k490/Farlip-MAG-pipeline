"""CLI orchestrator for the FaRLiP ingestion pipeline.

Subcommands (see ``python -m src.pipeline --help``):

* ``fetch references`` - pull every reference anchor declared in
  ``config/references.json`` via NCBI Datasets v2.
* ``fetch mags`` - pull every BioProject in ``config/extremophiles.json``
  (enumerate then download) plus any ``individual_accessions`` or already-populated
  ``cave_arctic_pending`` rows. Both NCBI and ENA accessions are dispatched to the
  correct fetcher based on the ``source`` field.
* ``fetch databases`` - pull SILVA + GTDB assets declared in ``config/databases.yaml``.
* ``fetch all`` - run all three of the above in sequence.
* ``verify`` - recompute SHA-256 for every manifest row and flag changes.
* ``status`` - print a one-screen summary of what is on disk.

Every successful or failed fetch results in a row in ``data/manifest.tsv`` so iterations
2+ have an authoritative listing of available artifacts. Re-running any subcommand is
idempotent: existing artifacts are skipped (unless ``--overwrite`` is passed).
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

from ._http import FetchError, FetchResult, compute_sha256, utcnow_iso
from .config_loader import (
    DEFAULT_CONFIG_DIR,
    BioProjectEntry,
    IndividualAccession,
    PipelineConfig,
    load_all,
)
from .ena_fetcher import fetch_assembly_fasta, fetch_bioproject
from .entrez_fetcher import fetch_assembly as fetch_ncbi_assembly
from .gtdb_fetcher import fetch_gtdb
from .manifest import MANIFEST_COLUMNS, read_manifest, upsert_results, write_manifest
from .silva_fetcher import fetch_silva


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATA_ROOT = REPO_ROOT / "data"
DEFAULT_MANIFEST = DEFAULT_DATA_ROOT / "manifest.tsv"


def _print_results(results: list[FetchResult], header: str) -> None:
    """Pretty-print a list of FetchResults grouped by status to stderr."""
    by_status: Counter[str] = Counter(r.status for r in results)
    sys.stderr.write(f"\n[{header}] {len(results)} result(s): {dict(by_status)}\n")
    for r in results:
        marker = {"ok": "+", "skipped": "=", "failed": "!"}.get(r.status, "?")
        sys.stderr.write(f"  {marker} {r.source:5s} {r.accession:24s} {r.status:8s} {r.local_path}\n")


def fetch_references(cfg: PipelineConfig, data_root: Path, *, overwrite: bool) -> list[FetchResult]:
    """Fetch the primary and any secondary reference anchors via NCBI Datasets v2."""
    dest = data_root / "references"
    results: list[FetchResult] = []
    for ref in cfg.references.all_entries():
        try:
            results.append(
                fetch_ncbi_assembly(ref.refseq_accession, dest, label=ref.label, overwrite=overwrite)
            )
        except FetchError as exc:
            results.append(
                FetchResult(
                    accession=ref.refseq_accession,
                    source="NCBI",
                    url="(failed before request)",
                    local_path=dest / f"{ref.refseq_accession}.zip",
                    fetched_at=utcnow_iso(),
                    status="failed",
                    note=str(exc),
                    extra={"label": ref.label},
                )
            )
    return results


def _fetch_individual(entry: IndividualAccession, data_root: Path, *, overwrite: bool) -> FetchResult:
    """Dispatch one individual accession to the right fetcher based on its ``source``."""
    if entry.source.upper() == "NCBI":
        dest = data_root / "mags" / "individual"
        return fetch_ncbi_assembly(entry.accession, dest, label=entry.label, overwrite=overwrite)
    if entry.source.upper() == "ENA":
        dest = data_root / "mags" / "individual"
        return fetch_assembly_fasta(entry.accession, dest, overwrite=overwrite)
    raise FetchError(f"Unknown source '{entry.source}' for accession {entry.accession}")


def fetch_mags(
    cfg: PipelineConfig,
    data_root: Path,
    *,
    overwrite: bool,
    limit: int | None = None,
) -> list[FetchResult]:
    """Fetch every BioProject + individual + (already-populated) cave_arctic_pending entry.

    The ``limit`` argument restricts per-BioProject fetches; individual accessions are
    always downloaded in full (they're already cherry-picked).
    """
    results: list[FetchResult] = []
    for bp in cfg.extremophiles.bioprojects:
        dest = data_root / "mags" / bp.accession
        if bp.source.upper() != "ENA":
            results.append(
                FetchResult(
                    accession=bp.accession,
                    source=bp.source,
                    url="(not implemented)",
                    local_path=dest,
                    fetched_at=utcnow_iso(),
                    status="failed",
                    note=f"BioProject source '{bp.source}' is not yet supported; only ENA is implemented in iteration 1.",
                    extra={"label": bp.label},
                )
            )
            continue
        results.extend(fetch_bioproject(bp.accession, dest, overwrite=overwrite, limit=limit))

    for entry in cfg.extremophiles.individual_accessions:
        try:
            results.append(_fetch_individual(entry, data_root, overwrite=overwrite))
        except FetchError as exc:
            results.append(
                FetchResult(
                    accession=entry.accession,
                    source=entry.source,
                    url="(failed before request)",
                    local_path=data_root / "mags" / "individual" / entry.accession,
                    fetched_at=utcnow_iso(),
                    status="failed",
                    note=str(exc),
                    extra={"label": entry.label, "biome": entry.biome or ""},
                )
            )

    for entry in cfg.extremophiles.cave_arctic_pending:
        try:
            results.append(_fetch_individual(entry, data_root, overwrite=overwrite))
        except FetchError as exc:
            results.append(
                FetchResult(
                    accession=entry.accession,
                    source=entry.source,
                    url="(failed before request)",
                    local_path=data_root / "mags" / "individual" / entry.accession,
                    fetched_at=utcnow_iso(),
                    status="failed",
                    note=str(exc),
                    extra={"label": entry.label, "biome": entry.biome or ""},
                )
            )
    return results


def fetch_databases(cfg: PipelineConfig, data_root: Path, *, overwrite: bool) -> list[FetchResult]:
    """Fetch SILVA + GTDB assets in that order."""
    results: list[FetchResult] = []
    results.extend(fetch_silva(cfg.databases, data_root, overwrite=overwrite))
    results.extend(fetch_gtdb(cfg.databases, data_root, overwrite=overwrite))
    return results


def cmd_fetch(args: argparse.Namespace) -> int:
    """Handler for ``fetch references|mags|databases|all``."""
    cfg = load_all(Path(args.config_dir) if args.config_dir else None)
    data_root = Path(args.data_root)

    results: list[FetchResult] = []
    target = args.target

    if target in {"references", "all"}:
        out = fetch_references(cfg, data_root, overwrite=args.overwrite)
        _print_results(out, "references")
        results.extend(out)
    if target in {"mags", "all"}:
        out = fetch_mags(cfg, data_root, overwrite=args.overwrite, limit=args.limit)
        _print_results(out, "mags")
        results.extend(out)
    if target in {"databases", "all"}:
        out = fetch_databases(cfg, data_root, overwrite=args.overwrite)
        _print_results(out, "databases")
        results.extend(out)

    n = upsert_results(Path(args.manifest), results)
    sys.stderr.write(f"\nManifest now has {n} row(s) at {args.manifest}\n")
    return 0 if all(r.status != "failed" for r in results) else 2


def cmd_verify(args: argparse.Namespace) -> int:
    """Recompute SHA-256 for every manifest row that has a local file; flag mismatches."""
    manifest_path = Path(args.manifest)
    rows = read_manifest(manifest_path)
    if not rows:
        sys.stderr.write(f"Manifest is empty or missing: {manifest_path}\n")
        return 0

    n_ok = n_changed = n_missing = 0
    for row in rows:
        path = Path(row.get("local_path", ""))
        recorded = row.get("sha256", "")
        if not path.exists():
            row["status"] = "missing"
            row["note"] = f"local_path absent at verify time ({utcnow_iso()})"
            n_missing += 1
            continue
        actual = compute_sha256(path)
        if actual != recorded:
            row["status"] = "changed"
            row["note"] = f"sha256 drift: was {recorded[:12]}..., now {actual[:12]}..."
            row["sha256"] = actual
            n_changed += 1
        else:
            n_ok += 1

    write_manifest(manifest_path, rows)
    sys.stderr.write(f"verify: ok={n_ok} changed={n_changed} missing={n_missing}\n")
    return 0 if n_changed == 0 and n_missing == 0 else 2


def cmd_status(args: argparse.Namespace) -> int:
    """Print a compact summary of the manifest grouped by source and status."""
    manifest_path = Path(args.manifest)
    rows = read_manifest(manifest_path)
    if not rows:
        sys.stdout.write(f"Manifest empty or missing: {manifest_path}\n")
        return 0

    by_source: Counter[str] = Counter(r.get("source", "?") for r in rows)
    by_status: Counter[str] = Counter(r.get("status", "?") for r in rows)
    sys.stdout.write(f"Manifest: {manifest_path}\n")
    sys.stdout.write(f"Total rows: {len(rows)}\n")
    sys.stdout.write(f"By source : {dict(by_source)}\n")
    sys.stdout.write(f"By status : {dict(by_status)}\n")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level argparse parser with all subcommands wired up."""
    p = argparse.ArgumentParser(
        prog="python -m src.pipeline",
        description="FaRLiP ingestion pipeline (iteration 1: data acquisition).",
    )
    p.add_argument("--config-dir", default=str(DEFAULT_CONFIG_DIR), help="Path to config/ directory.")
    p.add_argument("--data-root", default=str(DEFAULT_DATA_ROOT), help="Root output directory for downloads.")
    p.add_argument("--manifest", default=str(DEFAULT_MANIFEST), help="Path to manifest TSV.")

    sub = p.add_subparsers(dest="cmd", required=True)

    fetch_p = sub.add_parser("fetch", help="Download data.")
    fetch_p.add_argument(
        "target",
        choices=["references", "mags", "databases", "all"],
        help="Which slice of the pipeline to fetch.",
    )
    fetch_p.add_argument("--overwrite", action="store_true", help="Re-download even if already present.")
    fetch_p.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Cap the number of assemblies fetched per BioProject. Useful for smoke tests.",
    )
    fetch_p.set_defaults(func=cmd_fetch)

    verify_p = sub.add_parser("verify", help="Recompute SHA-256 for every manifest row.")
    verify_p.set_defaults(func=cmd_verify)

    status_p = sub.add_parser("status", help="Print a compact manifest summary.")
    status_p.set_defaults(func=cmd_status)

    return p


def main(argv: list[str] | None = None) -> int:
    """Entrypoint for ``python -m src.pipeline``."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
