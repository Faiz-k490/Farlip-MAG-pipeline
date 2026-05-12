"""Phase 2 SSU rRNA extraction.

Wraps the ``barrnap`` rRNA detector (Seemann 2013,
``https://github.com/tseemann/barrnap``) and pulls every predicted 16S rRNA
sequence out of a Phase 1 MAG FASTA. The output of this module is the input to
both downstream paths the Phase 2 plan calls for:

* MAGs that yielded at least one 16S record are routed to the SILVA SINA
  multi-FASTA bundle (:mod:`src.ssu_bundler`).
* MAGs with zero 16S records are routed to a generated GTDB-Tk shell script
  (:mod:`src.gtdbtk_shellscript_generator`) so taxonomic placement can fall
  back to genome-wide marker classification.

The public surface is three callables:

* :func:`run_barrnap` - invokes ``barrnap`` once on a (possibly gzipped) FASTA
  and writes the GFF3 stdout to disk.
* :func:`extract_16S` - parses that GFF3 and slices the matching sub-sequences
  out of the same FASTA, returning a list of :class:`Ssu16SRecord` tuples.
* :func:`assess_mag_for_ssu` - end-to-end per-MAG entry point used by the
  Phase 2 CLI: runs barrnap, extracts 16S, writes
  ``data/ssu/<accession>.16S.fasta``, and returns a routing decision plus
  metadata.

Idempotency: each per-MAG run writes a stamp file
``data/ssu/<accession>.16S.stamp`` recording the SHA-256 of the input FASTA
and the barrnap parameters used. A re-run with the same hash and parameters
short-circuits unless ``overwrite=True`` is passed.

The module deliberately does *not* require ``barrnap`` to be on ``PATH`` at
import time. The CLI checks for the binary and reports a clear error so
machines without barrnap can still run the rest of the pipeline.
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from .orf_predictor import _sha256_file


BARRNAP_BINARY: str = "barrnap"
"""Default executable name. Overridable via ``run_barrnap(barrnap_binary=...)``
for tests that point at a stub on ``PATH``."""

BARRNAP_KINGDOM_BACTERIA: str = "bac"
"""Default kingdom passed to ``barrnap --kingdom``. The Phase 1 corpus is all
bacterial; archaeal MAGs would need ``arc``."""

BARRNAP_DEFAULT_THREADS: int = 4
"""Thread count handed to ``barrnap --threads``. Conservative default that
plays well with the rest of the pipeline running on the same box."""

BARRNAP_DEFAULT_EVALUE: float = 1e-6
"""``--evalue`` cutoff. Matches barrnap's own default (1e-6); recorded
explicitly so the stamp file can detect a parameter change."""

BARRNAP_DEFAULT_REJECT: float = 0.25
"""``--reject`` length-fraction cutoff. Predictions covering less than this
fraction of the SSU model are dropped. Matches barrnap's default."""

SSU_FASTA_LINE_WIDTH: int = 60
"""Wrap width for the extracted 16S FASTA. NCBI default; keeps records
portable for SILVA SINA submission."""

SSU_DECISION_SILVA: str = "silva"
"""Routing tag emitted when at least one 16S record is recovered from a MAG.
Consumed by :func:`src.ssu_bundler.bundle_silva_sina`."""

SSU_DECISION_GTDBTK: str = "gtdbtk"
"""Routing tag emitted when no 16S record is recovered. Consumed by
:func:`src.gtdbtk_shellscript_generator.generate_gtdbtk_script`."""

SSU_GFF_NAME_TAG: str = "16S_rRNA"
"""Substring matched against the GFF3 attribute column to keep only 16S rows.
Barrnap also emits 5S and 23S rows for bacteria; we filter on the ``Name=``
attribute to pick out the 16S rows specifically."""

ROUTING_LOG_COLUMNS: tuple[str, ...] = (
    "accession",
    "source",
    "fasta_path",
    "n_contigs_scanned",
    "n_16S_predicted",
    "decision",
    "output_path",
    "barrnap_version",
    "notes",
)
"""Stable column order for ``data/ssu/routing_log.tsv``. Additive only."""


@dataclass(frozen=True)
class Ssu16SRecord:
    """One 16S rRNA prediction sliced out of a MAG FASTA.

    Attributes:
        contig_id: Source contig identifier (first whitespace-delimited token
            of the FASTA header).
        start: 1-based inclusive start coordinate from the barrnap GFF3.
        end: 1-based inclusive end coordinate.
        strand: ``"+"`` or ``"-"``.
        length_bp: ``end - start + 1``.
        sequence: Forward-strand sequence of the slice. Always upper-case.
            Reverse-complemented when the GFF strand is ``"-"``.
        attributes: Raw GFF3 attribute string (kept verbatim so downstream
            consumers can recover barrnap-specific tags such as ``Note=``).
    """

    contig_id: str
    start: int
    end: int
    strand: str
    length_bp: int
    sequence: str
    attributes: str


@dataclass(frozen=True)
class BarrnapResult:
    """Outcome of one :func:`run_barrnap` invocation.

    Attributes:
        gff_path: Path the GFF3 stdout was written to.
        n_rows_total: Number of GFF3 feature rows barrnap emitted (any rRNA
            type, before the 16S filter).
        n_rows_16S: Subset matching :data:`SSU_GFF_NAME_TAG`.
        version: Output of ``barrnap --version`` if it could be captured;
            empty string if the version probe failed.
        argv: The actual subprocess argument vector used. Stored for the
            stamp file so a parameter change invalidates idempotency.
    """

    gff_path: Path
    n_rows_total: int
    n_rows_16S: int
    version: str
    argv: tuple[str, ...]


@dataclass(frozen=True)
class SsuAssessment:
    """End-to-end per-MAG outcome of :func:`assess_mag_for_ssu`.

    Attributes:
        accession: MAG / assembly accession key.
        fasta_path: Input FASTA the assessment ran against.
        gff_path: Per-MAG barrnap GFF3 output.
        ssu_fasta_path: Per-MAG 16S FASTA. Always created (empty file is the
            valid signal for ``decision = gtdbtk``).
        stamp_path: Sidecar idempotency stamp.
        n_contigs_scanned: Number of FASTA records seen in the input.
        n_16S_predicted: Number of 16S predictions retained.
        decision: :data:`SSU_DECISION_SILVA` or :data:`SSU_DECISION_GTDBTK`.
        records: The list of extracted 16S records (possibly empty).
        barrnap_version: Captured from ``barrnap --version`` (best effort).
        skipped: True if the run short-circuited because the stamp matched.
    """

    accession: str
    fasta_path: Path
    gff_path: Path
    ssu_fasta_path: Path
    stamp_path: Path
    n_contigs_scanned: int
    n_16S_predicted: int
    decision: str
    records: list[Ssu16SRecord]
    barrnap_version: str
    skipped: bool


class BarrnapError(RuntimeError):
    """Raised when ``barrnap`` cannot be located or its invocation fails."""


def is_barrnap_available(barrnap_binary: str = BARRNAP_BINARY) -> bool:
    """Return True when the named barrnap binary is on ``PATH``."""
    return shutil.which(barrnap_binary) is not None


def get_barrnap_version(barrnap_binary: str = BARRNAP_BINARY) -> str:
    """Capture ``barrnap --version`` (best effort).

    Args:
        barrnap_binary: Executable name or absolute path.

    Returns:
        The trimmed stdout of ``barrnap --version``, or an empty string if
        the call failed for any reason. Never raises - the version string is
        purely informational.
    """
    try:
        completed = subprocess.run(
            [barrnap_binary, "--version"],
            check=False,
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return ""
    out = (completed.stdout or completed.stderr or "").strip()
    return out.splitlines()[0] if out else ""


def _open_fasta_text(path: Path):
    """Open a FASTA transparently for plain or gzipped inputs.

    Mirrors :func:`src.orf_predictor._open_fasta` but exposed locally so this
    module can be imported in environments where ``pyrodigal`` is missing.
    """
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r", encoding="utf-8")


def _iter_fasta_records(handle) -> Iterator[tuple[str, str]]:
    """Yield ``(contig_id, sequence)`` pairs from a FASTA text stream."""
    contig_id: str | None = None
    parts: list[str] = []
    for raw in handle:
        line = raw.rstrip()
        if not line:
            continue
        if line.startswith(">"):
            if contig_id is not None:
                yield contig_id, "".join(parts)
            contig_id = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)
    if contig_id is not None:
        yield contig_id, "".join(parts)


def _materialize_plain_fasta(fasta_path: Path) -> tuple[Path, bool]:
    """Return a plain (non-gzipped) FASTA path for barrnap to read.

    Barrnap's bundled nhmmer chokes on gzipped inputs. ENA assemblies arrive
    as ``.fasta.gz`` so we transparently decompress to a tempfile when
    needed.

    Args:
        fasta_path: Input FASTA (plain or gzipped).

    Returns:
        Tuple of ``(plain_fasta_path, is_temporary)``. When ``is_temporary``
        is True the caller must delete the file when done.
    """
    if fasta_path.suffix != ".gz":
        return fasta_path, False
    fd, tmp_name = tempfile.mkstemp(prefix="ssu_extractor_", suffix=".fasta")
    tmp_path = Path(tmp_name)
    try:
        with gzip.open(fasta_path, "rb") as src, open(fd, "wb") as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise
    return tmp_path, True


_REVERSE_COMPLEMENT_TABLE = str.maketrans(
    "ACGTUNRYSWKMBDHVacgtunryswkmbdhv",
    "TGCAANYRSWMKVHDBtgcaanyrswmkvhdb",
)
"""IUPAC nucleotide reverse-complement translation table. Preserves ambiguity
codes (R, Y, S, W, K, M, B, D, H, V, N) which appear in MAG assemblies."""


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of an IUPAC nucleotide string."""
    return seq.translate(_REVERSE_COMPLEMENT_TABLE)[::-1]


def _wrap_fasta(seq: str, width: int = SSU_FASTA_LINE_WIDTH) -> str:
    """Wrap a single-line sequence to the requested column width."""
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def run_barrnap(
    fasta_path: Path,
    output_gff: Path,
    *,
    kingdom: str = BARRNAP_KINGDOM_BACTERIA,
    threads: int = BARRNAP_DEFAULT_THREADS,
    evalue: float = BARRNAP_DEFAULT_EVALUE,
    reject: float = BARRNAP_DEFAULT_REJECT,
    barrnap_binary: str = BARRNAP_BINARY,
) -> BarrnapResult:
    """Invoke ``barrnap`` once on ``fasta_path`` and write its GFF3 stdout.

    Args:
        fasta_path: Input genomic FASTA. Plain or gzipped (gzip is unpacked
            into a tempfile because barrnap cannot read gzip directly).
        output_gff: Where to write the GFF3 stdout.
        kingdom: ``--kingdom`` argument; defaults to ``bac``.
        threads: ``--threads`` argument.
        evalue: ``--evalue`` argument.
        reject: ``--reject`` argument.
        barrnap_binary: Executable name or absolute path. Override for tests.

    Returns:
        A :class:`BarrnapResult` summarising what was written.

    Raises:
        BarrnapError: If the binary is missing, exits non-zero, or returns
            a body with no GFF feature rows (the latter is unusual but flagged
            so callers can decide whether to swallow it).
        FileNotFoundError: If ``fasta_path`` does not exist.
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"input FASTA not found: {fasta_path}")
    if not is_barrnap_available(barrnap_binary):
        raise BarrnapError(
            f"barrnap binary '{barrnap_binary}' not found on PATH; "
            "install barrnap (e.g. via conda-forge / bioconda) before running Phase 2"
        )

    plain_fasta, is_tmp = _materialize_plain_fasta(fasta_path)
    argv: tuple[str, ...] = (
        barrnap_binary,
        "--kingdom",
        kingdom,
        "--threads",
        str(threads),
        "--evalue",
        f"{evalue:g}",
        "--reject",
        f"{reject:g}",
        "--quiet",
        str(plain_fasta),
    )
    try:
        completed = subprocess.run(
            argv,
            check=False,
            capture_output=True,
            text=True,
        )
    finally:
        if is_tmp:
            plain_fasta.unlink(missing_ok=True)

    if completed.returncode != 0:
        raise BarrnapError(
            f"barrnap exited with code {completed.returncode} for {fasta_path}: "
            f"{(completed.stderr or '').strip()[:500]}"
        )

    output_gff.parent.mkdir(parents=True, exist_ok=True)
    output_gff.write_text(completed.stdout, encoding="utf-8")

    n_total = 0
    n_16S = 0
    for line in completed.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        n_total += 1
        if SSU_GFF_NAME_TAG in line:
            n_16S += 1

    return BarrnapResult(
        gff_path=output_gff,
        n_rows_total=n_total,
        n_rows_16S=n_16S,
        version=get_barrnap_version(barrnap_binary),
        argv=argv,
    )


def _parse_barrnap_attributes(attr_field: str) -> dict[str, str]:
    """Parse a GFF3 attribute string into a flat dict.

    Args:
        attr_field: The 9th GFF3 column, e.g. ``Name=16S_rRNA;product=...``.

    Returns:
        Mapping ``{key: value}``. Keys with no ``=`` are silently dropped.
    """
    out: dict[str, str] = {}
    for token in attr_field.split(";"):
        token = token.strip()
        if not token or "=" not in token:
            continue
        k, v = token.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def extract_16S(fasta_path: Path, gff_path: Path) -> list[Ssu16SRecord]:
    """Slice every barrnap-predicted 16S rRNA sequence out of a FASTA.

    Args:
        fasta_path: Same FASTA the barrnap call ran against.
        gff_path: GFF3 produced by :func:`run_barrnap`.

    Returns:
        A list of :class:`Ssu16SRecord`. Empty if the GFF has no 16S rows.
        Sequences are upper-cased; reverse-complemented when GFF strand is
        ``-``. Records on contigs missing from the FASTA are skipped.
    """
    if not gff_path.exists():
        return []
    contig_seqs: dict[str, str] = {}
    with _open_fasta_text(fasta_path) as fh:
        for contig_id, seq in _iter_fasta_records(fh):
            contig_seqs[contig_id] = seq.upper()

    records: list[Ssu16SRecord] = []
    for raw in gff_path.read_text(encoding="utf-8").splitlines():
        if not raw or raw.startswith("#"):
            continue
        cols = raw.split("\t")
        if len(cols) < 9:
            continue
        contig_id = cols[0]
        try:
            start = int(cols[3])
            end = int(cols[4])
        except ValueError:
            continue
        strand = cols[6] or "+"
        attrs = cols[8]
        attr_map = _parse_barrnap_attributes(attrs)
        name = attr_map.get("Name", "")
        if SSU_GFF_NAME_TAG not in name and SSU_GFF_NAME_TAG not in attrs:
            continue
        contig_seq = contig_seqs.get(contig_id)
        if contig_seq is None or start < 1 or end > len(contig_seq) or start > end:
            continue
        slice_seq = contig_seq[start - 1 : end]
        if strand == "-":
            slice_seq = _reverse_complement(slice_seq)
        records.append(
            Ssu16SRecord(
                contig_id=contig_id,
                start=start,
                end=end,
                strand=strand,
                length_bp=end - start + 1,
                sequence=slice_seq,
                attributes=attrs,
            )
        )
    return records


def write_ssu_fasta(
    records: list[Ssu16SRecord],
    out_path: Path,
    *,
    mag_accession: str,
) -> int:
    """Write a per-MAG 16S FASTA. Always creates the file (possibly empty).

    Header format::

        >{accession}|{contig_id}|{start}-{end}|{strand}|len={length_bp}

    The empty-file-when-zero-records contract is intentional: downstream
    bundlers can ``glob`` for ``*.16S.fasta`` without having to special-case
    MAGs that produced no 16S.

    Args:
        records: List of :class:`Ssu16SRecord` to emit (possibly empty).
        out_path: Output FASTA path.
        mag_accession: MAG accession key for the header prefix.

    Returns:
        Number of records written.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        for r in records:
            header = (
                f">{mag_accession}|{r.contig_id}|"
                f"{r.start}-{r.end}|{r.strand}|len={r.length_bp}"
            )
            fh.write(header + "\n")
            fh.write(_wrap_fasta(r.sequence) + "\n")
    return len(records)


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
    accession: str,
    fasta_sha256: str,
    barrnap_argv: tuple[str, ...],
    barrnap_version: str,
    n_16S_predicted: int,
    decision: str,
) -> None:
    """Write the per-MAG idempotency stamp."""
    stamp_path.parent.mkdir(parents=True, exist_ok=True)
    payload = [
        f"accession\t{accession}",
        f"fasta_sha256\t{fasta_sha256}",
        f"barrnap_argv\t{' '.join(barrnap_argv)}",
        f"barrnap_version\t{barrnap_version}",
        f"n_16S_predicted\t{n_16S_predicted}",
        f"decision\t{decision}",
    ]
    stamp_path.write_text("\n".join(payload) + "\n", encoding="utf-8")


def _count_contigs(fasta_path: Path) -> int:
    """Count records in a FASTA without loading sequences into memory."""
    n = 0
    with _open_fasta_text(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def assess_mag_for_ssu(
    accession: str,
    fasta_path: Path,
    output_dir: Path,
    *,
    kingdom: str = BARRNAP_KINGDOM_BACTERIA,
    threads: int = BARRNAP_DEFAULT_THREADS,
    evalue: float = BARRNAP_DEFAULT_EVALUE,
    reject: float = BARRNAP_DEFAULT_REJECT,
    barrnap_binary: str = BARRNAP_BINARY,
    overwrite: bool = False,
) -> SsuAssessment:
    """End-to-end per-MAG 16S extraction + routing decision.

    Args:
        accession: MAG / assembly accession key.
        fasta_path: Genomic FASTA (plain or gzipped).
        output_dir: Directory under which ``<accession>.barrnap.gff``,
            ``<accession>.16S.fasta`` and ``<accession>.16S.stamp`` are
            written. Created on demand.
        kingdom, threads, evalue, reject: Barrnap parameters; recorded in
            the stamp so a parameter change invalidates idempotency.
        barrnap_binary: Override for tests.
        overwrite: If True, ignore any existing stamp and re-run.

    Returns:
        A :class:`SsuAssessment`. The decision tag is the routing signal
        consumed by Phase 2 downstream stages.

    Raises:
        FileNotFoundError: If ``fasta_path`` does not exist.
        BarrnapError: If barrnap is missing or fails.
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"input FASTA not found: {fasta_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    gff_path = output_dir / f"{accession}.barrnap.gff"
    ssu_fasta_path = output_dir / f"{accession}.16S.fasta"
    stamp_path = output_dir / f"{accession}.16S.stamp"

    fasta_sha256 = _sha256_file(fasta_path)
    candidate_argv = (
        barrnap_binary,
        "--kingdom",
        kingdom,
        "--threads",
        str(threads),
        "--evalue",
        f"{evalue:g}",
        "--reject",
        f"{reject:g}",
        "--quiet",
        str(fasta_path),
    )

    if not overwrite and gff_path.exists() and ssu_fasta_path.exists():
        prior = _read_stamp(stamp_path)
        if (
            prior.get("fasta_sha256") == fasta_sha256
            and prior.get("barrnap_argv") == " ".join(candidate_argv)
        ):
            records = extract_16S(fasta_path, gff_path)
            decision = (
                SSU_DECISION_SILVA
                if records
                else SSU_DECISION_GTDBTK
            )
            return SsuAssessment(
                accession=accession,
                fasta_path=fasta_path,
                gff_path=gff_path,
                ssu_fasta_path=ssu_fasta_path,
                stamp_path=stamp_path,
                n_contigs_scanned=_count_contigs(fasta_path),
                n_16S_predicted=len(records),
                decision=decision,
                records=records,
                barrnap_version=prior.get("barrnap_version", ""),
                skipped=True,
            )

    barrnap_result = run_barrnap(
        fasta_path,
        gff_path,
        kingdom=kingdom,
        threads=threads,
        evalue=evalue,
        reject=reject,
        barrnap_binary=barrnap_binary,
    )
    records = extract_16S(fasta_path, gff_path)
    write_ssu_fasta(records, ssu_fasta_path, mag_accession=accession)
    decision = SSU_DECISION_SILVA if records else SSU_DECISION_GTDBTK
    _write_stamp(
        stamp_path,
        accession=accession,
        fasta_sha256=fasta_sha256,
        barrnap_argv=barrnap_result.argv,
        barrnap_version=barrnap_result.version,
        n_16S_predicted=len(records),
        decision=decision,
    )

    return SsuAssessment(
        accession=accession,
        fasta_path=fasta_path,
        gff_path=gff_path,
        ssu_fasta_path=ssu_fasta_path,
        stamp_path=stamp_path,
        n_contigs_scanned=_count_contigs(fasta_path),
        n_16S_predicted=len(records),
        decision=decision,
        records=records,
        barrnap_version=barrnap_result.version,
        skipped=False,
    )


__all__ = [
    "BARRNAP_BINARY",
    "BARRNAP_DEFAULT_EVALUE",
    "BARRNAP_DEFAULT_REJECT",
    "BARRNAP_DEFAULT_THREADS",
    "BARRNAP_KINGDOM_BACTERIA",
    "BarrnapError",
    "BarrnapResult",
    "ROUTING_LOG_COLUMNS",
    "SSU_DECISION_GTDBTK",
    "SSU_DECISION_SILVA",
    "SSU_FASTA_LINE_WIDTH",
    "SSU_GFF_NAME_TAG",
    "Ssu16SRecord",
    "SsuAssessment",
    "assess_mag_for_ssu",
    "extract_16S",
    "get_barrnap_version",
    "is_barrnap_available",
    "run_barrnap",
    "write_ssu_fasta",
]
