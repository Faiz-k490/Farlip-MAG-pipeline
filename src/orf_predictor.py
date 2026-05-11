"""ORF prediction for every assembly in the Phase 1 manifest.

Wraps `pyrodigal <https://pyrodigal.readthedocs.io/>`_ (a Cython binding for Prodigal)
in meta mode (the correct mode for unbinned and lightly-binned MAGs - it uses Prodigal's
pre-trained universal model rather than per-genome training, which would mis-train on
short contigs).

For each genomic FASTA the module emits three artifacts under ``data/orfs/``:

* ``<accession>.faa``        - amino-acid FASTA of every predicted ORF. Headers carry
                                ``accession|contig_id|orf_N`` so the downstream
                                synteny parser can join hits back to contig coordinates
                                in a single pass.
* ``<accession>.coords.tsv`` - one row per ORF with ``contig_id, orf_id, start, end,
                                strand, partial_left, partial_right, contig_length``.
                                A native TSV is faster to join than the GFF that
                                pyrodigal can otherwise emit, and is the only
                                format the synteny parser needs.
* ``<accession>.stamp``      - tiny key/value file recording the SHA-256 of the input
                                FASTA at run time. Re-runs short-circuit if the input
                                hash matches what is recorded here.

The module is **purely I/O-bound at the boundaries** and **purely in-memory between
them** - no intermediate FASTAs are written just to be re-read. Pyrodigal returns
amino-acid sequences directly from the C extension; we stream them straight to the
output FAA.
"""

from __future__ import annotations

import gzip
import hashlib
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

import pyrodigal


PYRODIGAL_META_MODE: bool = True
"""Pyrodigal is invoked in meta mode for every MAG. Single-genome training is unsafe
on fragmented metagenome assemblies (Prodigal's training step needs ~100 kb of
contiguous sequence per genome to converge; many polar MAG contigs are shorter)."""

FAA_LINE_WIDTH: int = 60
"""Standard NCBI FASTA wrap width. Keeps the FAA portable for downstream tools."""

CHUNK_BYTES: int = 1024 * 1024
"""Read chunk size for SHA-256 hashing of input FASTAs (1 MiB)."""

COORDS_HEADER: tuple[str, ...] = (
    "contig_id",
    "orf_id",
    "start",
    "end",
    "strand",
    "partial_left",
    "partial_right",
    "contig_length",
)
"""Column order for ``<accession>.coords.tsv``. Stable - synteny parser joins on it."""


@dataclass(frozen=True)
class OrfPredictionResult:
    """Summary of one ``predict_orfs`` invocation.

    Attributes:
        accession: The MAG / assembly accession key.
        faa_path: Path to the written amino-acid FASTA.
        coords_path: Path to the written coords TSV.
        stamp_path: Path to the written stamp file.
        n_contigs: Number of contigs processed.
        n_orfs: Number of ORFs predicted across all contigs.
        skipped: True if the call short-circuited because the stamp matched.
        input_sha256: SHA-256 of the input genomic FASTA.
    """

    accession: str
    faa_path: Path
    coords_path: Path
    stamp_path: Path
    n_contigs: int
    n_orfs: int
    skipped: bool
    input_sha256: str


def _open_fasta(path: Path):
    """Open a genomic FASTA transparently for both ``.fna`` and ``.fna.gz`` inputs.

    Args:
        path: Path to either a plain or gzipped FASTA.

    Returns:
        A text-mode file-like object. Caller is responsible for closing it.
    """
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def _sha256_file(path: Path) -> str:
    """Stream-hash a file with SHA-256.

    Args:
        path: Any file on disk (compressed or otherwise - we hash the bytes as-is).

    Returns:
        Hex-encoded SHA-256 digest of the file's contents.
    """
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK_BYTES), b""):
            h.update(chunk)
    return h.hexdigest()


def _iter_fasta_records(handle) -> Iterator[tuple[str, str]]:
    """Yield ``(contig_id, sequence)`` tuples from a FASTA text stream.

    We hand-roll the parser to avoid a Biopython dependency on the hot path - we only
    need the first whitespace-delimited token of the header and the upper-case nucleotide
    sequence joined into one string.

    Args:
        handle: Text-mode file object positioned at the start of a FASTA file.

    Yields:
        Pairs of ``(contig_id, sequence_string)``.
    """
    contig_id: str | None = None
    seq_parts: list[str] = []
    for raw in handle:
        line = raw.rstrip()
        if not line:
            continue
        if line.startswith(">"):
            if contig_id is not None:
                yield contig_id, "".join(seq_parts)
            contig_id = line[1:].split()[0]
            seq_parts = []
        else:
            seq_parts.append(line)
    if contig_id is not None:
        yield contig_id, "".join(seq_parts)


def _wrap_fasta(seq: str, width: int = FAA_LINE_WIDTH) -> str:
    """Wrap a single-line protein string to the requested column width.

    Args:
        seq: Amino-acid sequence.
        width: Column width. Defaults to :data:`FAA_LINE_WIDTH`.

    Returns:
        The sequence with newlines inserted every ``width`` residues.
    """
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def _read_stamp(stamp_path: Path) -> dict[str, str]:
    """Parse a stamp file into a ``{key: value}`` dict.

    Args:
        stamp_path: Path to a stamp file (returns empty dict if missing).

    Returns:
        Key/value mapping. Lines that don't contain a tab are silently ignored.
    """
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
    input_sha256: str,
    input_path: Path,
    n_contigs: int,
    n_orfs: int,
) -> None:
    """Write a stamp file recording what was processed and when.

    Args:
        stamp_path: Output path.
        input_sha256: SHA-256 of the input FASTA at the time the stamp was written.
        input_path: Path to the input FASTA.
        n_contigs: Number of contigs processed.
        n_orfs: Number of ORFs predicted.
    """
    stamp_path.parent.mkdir(parents=True, exist_ok=True)
    payload = [
        f"input_sha256\t{input_sha256}",
        f"input_path\t{input_path}",
        f"n_contigs\t{n_contigs}",
        f"n_orfs\t{n_orfs}",
        f"pyrodigal_version\t{pyrodigal.__version__}",
        f"generated_at\t{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
    ]
    stamp_path.write_text("\n".join(payload) + "\n", encoding="utf-8")


def predict_orfs(
    accession: str,
    fasta_path: Path,
    output_dir: Path,
    *,
    overwrite: bool = False,
) -> OrfPredictionResult:
    """Predict ORFs for one assembly with pyrodigal and write FAA + coords.tsv + stamp.

    Idempotency:
        Re-running with the same input FASTA short-circuits when the stamp file
        already exists and records the same SHA-256. Set ``overwrite=True`` to force
        a fresh run.

    Args:
        accession: MAG / assembly accession key. Drives the output filenames.
        fasta_path: Path to the genomic FASTA (plain ``.fna`` or gzipped ``.fna.gz``).
        output_dir: Directory under which ``<accession>.faa``, ``.coords.tsv`` and
            ``.stamp`` are written. Created if it doesn't exist.
        overwrite: If True, ignore any existing stamp and re-run.

    Returns:
        An :class:`OrfPredictionResult` describing what was written.

    Raises:
        FileNotFoundError: If ``fasta_path`` does not exist.
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"genomic FASTA not found: {fasta_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    faa_path = output_dir / f"{accession}.faa"
    coords_path = output_dir / f"{accession}.coords.tsv"
    stamp_path = output_dir / f"{accession}.stamp"

    input_sha256 = _sha256_file(fasta_path)

    if not overwrite and faa_path.exists() and coords_path.exists():
        prior = _read_stamp(stamp_path)
        if prior.get("input_sha256") == input_sha256:
            return OrfPredictionResult(
                accession=accession,
                faa_path=faa_path,
                coords_path=coords_path,
                stamp_path=stamp_path,
                n_contigs=int(prior.get("n_contigs", 0)),
                n_orfs=int(prior.get("n_orfs", 0)),
                skipped=True,
                input_sha256=input_sha256,
            )

    gene_finder = pyrodigal.GeneFinder(meta=PYRODIGAL_META_MODE)

    n_contigs = 0
    n_orfs = 0
    orf_serial = 0

    with (
        _open_fasta(fasta_path) as fh,
        faa_path.open("w", encoding="utf-8") as faa_out,
        coords_path.open("w", encoding="utf-8") as coords_out,
    ):
        coords_out.write("\t".join(COORDS_HEADER) + "\n")

        for contig_id, seq in _iter_fasta_records(fh):
            if not seq:
                continue
            n_contigs += 1
            contig_length = len(seq)
            genes = gene_finder.find_genes(seq.encode("ascii"))
            for gene in genes:
                orf_serial += 1
                n_orfs += 1
                orf_id = f"{accession}|{contig_id}|orf_{orf_serial:06d}"
                protein = gene.translate()
                if protein.endswith("*"):
                    protein = protein[:-1]
                faa_out.write(f">{orf_id}\n{_wrap_fasta(protein)}\n")
                strand_char = "+" if gene.strand == 1 else "-"
                coords_out.write(
                    "\t".join(
                        [
                            contig_id,
                            orf_id,
                            str(gene.begin),
                            str(gene.end),
                            strand_char,
                            "1" if gene.partial_begin else "0",
                            "1" if gene.partial_end else "0",
                            str(contig_length),
                        ]
                    )
                    + "\n"
                )

    _write_stamp(
        stamp_path,
        input_sha256=input_sha256,
        input_path=fasta_path,
        n_contigs=n_contigs,
        n_orfs=n_orfs,
    )

    return OrfPredictionResult(
        accession=accession,
        faa_path=faa_path,
        coords_path=coords_path,
        stamp_path=stamp_path,
        n_contigs=n_contigs,
        n_orfs=n_orfs,
        skipped=False,
        input_sha256=input_sha256,
    )
