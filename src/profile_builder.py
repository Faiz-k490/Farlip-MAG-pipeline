"""Build HMM profiles for the 21 canonical FaRLiP genes from JSC-1.

For Phase 3 iteration 1 we have a single FaRLiP anchor (Leptolyngbya sp. JSC-1,
RefSeq assembly GCF_000733415.1). With one sequence per gene we cannot run a
multiple-sequence alignment, so each HMM is built from a one-row "MSA" containing
just that sequence. This produces a reasonable starting profile thanks to pyhmmer's
default position-specific priors; sensitivity for highly divergent paralogs will
improve once the secondary anchor (Chlorogloeopsis fritschii PCC 9212) is integrated
in Phase 3.1.

Resolution priority for each catalog entry (highest first):

1. ``match.locus_tag`` - exact locus tag in the JSC-1 GFF.
2. ``match.protein_id`` - exact RefSeq protein accession.
3. ``match.gene_tag`` + ``match.product_substring`` - fuzzy match, with a
   cluster-locality tiebreak: candidates inside the canonical FaRLiP cluster
   window on contig ``NZ_KL662192.1`` win over candidates elsewhere in the genome.

Outputs:

* ``data/hmm_profiles/farlip_21.hmm``         - concatenated HMM file in HMMER3 format.
* ``data/hmm_profiles/farlip_21.catalog.tsv`` - sidecar mapping ``hmm_name`` to
  ``protein_id``, ``locus_tag``, ``gene`` (catalog name), ``role``, ``contig``,
  ``start``, ``end``, ``strand``, ``length_aa``, and ``in_cluster`` (boolean) so
  the synteny parser can label hits without re-parsing the GFF.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import unquote

import pyhmmer
import yaml


JSC1_FARLIP_CLUSTER_CONTIG: str = "NZ_KL662192.1"
"""JSC-1 RefSeq contig hosting the canonical 21-gene FaRLiP cluster.

The cluster spans approximately 1,095,910 - 1,121,582 bp on this contig (Gan et al.
2014; verified against the JSC-1 GFF, GCF_000733415.1, on 2026-05-10)."""

JSC1_FARLIP_CLUSTER_START: int = 1_095_910
"""Lower bound (1-based, inclusive) of the JSC-1 FaRLiP cluster window."""

JSC1_FARLIP_CLUSTER_END: int = 1_121_582
"""Upper bound (1-based, inclusive) of the JSC-1 FaRLiP cluster window."""

CATALOG_TSV_HEADER: tuple[str, ...] = (
    "hmm_name",
    "gene",
    "role",
    "protein_id",
    "locus_tag",
    "contig",
    "start",
    "end",
    "strand",
    "length_aa",
    "in_cluster",
    "resolution",
)
"""Column order for ``farlip_21.catalog.tsv``. Stable - synteny parser joins on this."""


class ProfileBuildError(RuntimeError):
    """Raised when the catalog cannot be resolved into a usable set of HMMs."""


@dataclass(frozen=True)
class CatalogEntry:
    """One gene entry as parsed from ``config/farlip_genes.yaml``."""

    name: str
    role: str
    locus_tag: str | None
    protein_id: str | None
    gene_tag: str | None
    product_substring: str | None


@dataclass(frozen=True)
class ResolvedGene:
    """A catalog entry resolved against the JSC-1 GFF + protein.faa."""

    name: str
    role: str
    protein_id: str
    locus_tag: str
    contig: str
    start: int
    end: int
    strand: str
    sequence: str
    in_cluster: bool
    resolution: str
    """How the gene was located: ``locus_tag``, ``protein_id``, or ``fuzzy``."""


@dataclass
class GffCDS:
    """One CDS row from a RefSeq GFF, normalized to the attributes we care about."""

    contig: str
    start: int
    end: int
    strand: str
    locus_tag: str | None = None
    protein_id: str | None = None
    gene_tag: str | None = None
    product: str | None = None


def load_catalog(path: Path) -> list[CatalogEntry]:
    """Load the FaRLiP gene catalog YAML into a list of :class:`CatalogEntry`.

    Args:
        path: Path to ``config/farlip_genes.yaml``.

    Returns:
        Catalog entries in the order they appear in the file.

    Raises:
        ProfileBuildError: If the YAML is missing the ``genes`` key or any entry
            is missing a ``name`` field.
    """
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or "genes" not in raw:
        raise ProfileBuildError(f"catalog YAML missing 'genes' key: {path}")
    out: list[CatalogEntry] = []
    for idx, entry in enumerate(raw["genes"]):
        if "name" not in entry:
            raise ProfileBuildError(f"catalog entry {idx} missing 'name'")
        match = entry.get("match", {}) or {}
        out.append(
            CatalogEntry(
                name=str(entry["name"]),
                role=str(entry.get("role", "")),
                locus_tag=match.get("locus_tag"),
                protein_id=match.get("protein_id"),
                gene_tag=match.get("gene_tag"),
                product_substring=match.get("product_substring"),
            )
        )
    return out


def _parse_gff_attributes(attr_field: str) -> dict[str, str]:
    """Parse a GFF3 attribute field into a key/value dict.

    Handles URL-escaped values (e.g. ``%2C`` for comma, ``%3B`` for semicolon).

    Args:
        attr_field: The 9th column of a GFF3 record.

    Returns:
        Mapping of attribute name to (URL-decoded) value.
    """
    out: dict[str, str] = {}
    for piece in attr_field.split(";"):
        if "=" not in piece:
            continue
        k, v = piece.split("=", 1)
        out[k.strip()] = unquote(v.strip())
    return out


def parse_jsc1_gff(gff_path: Path) -> list[GffCDS]:
    """Parse JSC-1's ``genomic.gff`` into a list of :class:`GffCDS`.

    Pseudogenes (rows with ``pseudo=true``) are skipped because they don't have a
    matching entry in ``protein.faa``.

    Args:
        gff_path: Path to ``genomic.gff`` inside the RefSeq Datasets bundle.

    Returns:
        One :class:`GffCDS` per non-pseudo CDS, in file order.
    """
    out: list[GffCDS] = []
    with gff_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            attrs = _parse_gff_attributes(cols[8])
            if attrs.get("pseudo") == "true":
                continue
            out.append(
                GffCDS(
                    contig=cols[0],
                    start=int(cols[3]),
                    end=int(cols[4]),
                    strand=cols[6],
                    locus_tag=attrs.get("locus_tag"),
                    protein_id=attrs.get("protein_id"),
                    gene_tag=attrs.get("gene"),
                    product=attrs.get("product"),
                )
            )
    return out


def parse_jsc1_proteins(faa_path: Path) -> dict[str, str]:
    """Parse JSC-1's ``protein.faa`` into a ``{protein_id: sequence}`` dict.

    Args:
        faa_path: Path to ``protein.faa`` inside the RefSeq Datasets bundle.

    Returns:
        Mapping of RefSeq protein accession (e.g. ``WP_036011370.1``) to its
        amino-acid sequence as a single uppercase string.
    """
    out: dict[str, str] = {}
    current_id: str | None = None
    buf: list[str] = []
    with faa_path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    out[current_id] = "".join(buf)
                current_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if current_id is not None:
            out[current_id] = "".join(buf)
    return out


def _is_in_cluster(cds: GffCDS) -> bool:
    """Return True iff a CDS sits inside the canonical FaRLiP cluster window."""
    return (
        cds.contig == JSC1_FARLIP_CLUSTER_CONTIG
        and cds.start >= JSC1_FARLIP_CLUSTER_START
        and cds.end <= JSC1_FARLIP_CLUSTER_END
    )


def _resolve_entry(
    entry: CatalogEntry,
    cds_list: list[GffCDS],
    proteins: dict[str, str],
) -> ResolvedGene:
    """Resolve one catalog entry to a concrete CDS + protein sequence.

    Resolution order:

    1. ``locus_tag`` exact match.
    2. ``protein_id`` exact match.
    3. ``gene_tag`` + ``product_substring`` fuzzy match, with cluster-locality
       tiebreak when more than one CDS qualifies.

    Args:
        entry: Catalog entry to resolve.
        cds_list: Parsed JSC-1 GFF CDS rows.
        proteins: ``{protein_id: sequence}`` dict from ``protein.faa``.

    Returns:
        A :class:`ResolvedGene` carrying everything the HMM build needs.

    Raises:
        ProfileBuildError: If no CDS can be located or the resolved CDS has no
            matching sequence in ``protein.faa``.
    """
    resolved: GffCDS | None = None
    resolution: str = ""

    if entry.locus_tag:
        for cds in cds_list:
            if cds.locus_tag == entry.locus_tag:
                resolved = cds
                resolution = "locus_tag"
                break

    if resolved is None and entry.protein_id:
        for cds in cds_list:
            if cds.protein_id == entry.protein_id:
                resolved = cds
                resolution = "protein_id"
                break

    if resolved is None and (entry.gene_tag or entry.product_substring):
        candidates: list[GffCDS] = []
        sub = (entry.product_substring or "").lower()
        for cds in cds_list:
            if entry.gene_tag and cds.gene_tag != entry.gene_tag:
                continue
            if sub and (cds.product is None or sub not in cds.product.lower()):
                continue
            candidates.append(cds)
        in_cluster = [c for c in candidates if _is_in_cluster(c)]
        if len(in_cluster) == 1:
            resolved = in_cluster[0]
            resolution = "fuzzy_cluster"
        elif len(in_cluster) > 1:
            resolved = in_cluster[0]
            resolution = "fuzzy_cluster_first"
        elif len(candidates) == 1:
            resolved = candidates[0]
            resolution = "fuzzy_genome"

    if resolved is None or resolved.protein_id is None:
        raise ProfileBuildError(
            f"could not resolve catalog entry {entry.name!r} (locus_tag="
            f"{entry.locus_tag!r}, protein_id={entry.protein_id!r}, "
            f"gene_tag={entry.gene_tag!r}, product_substring={entry.product_substring!r})"
        )

    sequence = proteins.get(resolved.protein_id)
    if sequence is None:
        raise ProfileBuildError(
            f"catalog entry {entry.name!r} resolved to {resolved.protein_id!r} "
            f"but no matching sequence in protein.faa"
        )

    return ResolvedGene(
        name=entry.name,
        role=entry.role,
        protein_id=resolved.protein_id,
        locus_tag=resolved.locus_tag or "",
        contig=resolved.contig,
        start=resolved.start,
        end=resolved.end,
        strand=resolved.strand,
        sequence=sequence,
        in_cluster=_is_in_cluster(resolved),
        resolution=resolution,
    )


_HMM_NAME_SANITIZER = re.compile(r"[^A-Za-z0-9_.-]")


def _hmm_name(gene_name: str) -> str:
    """Sanitize a gene name into a value HMMER will accept as an HMM identifier.

    HMMER allows letters, digits, underscore, hyphen, and period; everything else
    is collapsed to underscore.

    Args:
        gene_name: Gene name from the catalog (e.g. ``rfpA``, ``apcE2``).

    Returns:
        Sanitized name suitable for the HMMER ``NAME`` field.
    """
    return _HMM_NAME_SANITIZER.sub("_", gene_name)


def _build_single_sequence_hmm(
    name: str,
    sequence: str,
    builder: pyhmmer.plan7.Builder,
    background: pyhmmer.plan7.Background,
    alphabet: pyhmmer.easel.Alphabet,
) -> pyhmmer.plan7.HMM:
    """Build a single-sequence HMM via pyhmmer's plan7 builder.

    Args:
        name: HMM name (becomes the ``NAME`` field in the HMMER file).
        sequence: Amino-acid sequence (single string, no whitespace).
        builder: A configured :class:`pyhmmer.plan7.Builder`.
        background: A configured :class:`pyhmmer.plan7.Background`.
        alphabet: The amino-acid alphabet shared with ``builder`` and ``background``.

    Returns:
        The built :class:`pyhmmer.plan7.HMM` with its ``name`` set.
    """
    seq = pyhmmer.easel.TextSequence(
        name=name.encode("ascii"),
        sequence=sequence,
    ).digitize(alphabet)
    hmm, _profile, _optprofile = builder.build(seq, background)
    hmm.name = name.encode("ascii")
    return hmm


def build_profiles(
    catalog_path: Path,
    gff_path: Path,
    protein_faa_path: Path,
    output_dir: Path,
    *,
    overwrite: bool = False,
) -> list[ResolvedGene]:
    """Resolve every catalog entry against JSC-1 and write the combined HMM file.

    Args:
        catalog_path: Path to ``config/farlip_genes.yaml``.
        gff_path: Path to JSC-1's ``genomic.gff``.
        protein_faa_path: Path to JSC-1's ``protein.faa``.
        output_dir: Directory under which ``farlip_21.hmm`` and
            ``farlip_21.catalog.tsv`` will be written. Created if absent.
        overwrite: If True, rebuild even if the HMM file already exists.

    Returns:
        List of :class:`ResolvedGene` in catalog order.

    Raises:
        ProfileBuildError: If any catalog entry cannot be resolved.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    hmm_path = output_dir / "farlip_21.hmm"
    catalog_tsv_path = output_dir / "farlip_21.catalog.tsv"

    catalog = load_catalog(catalog_path)
    cds_list = parse_jsc1_gff(gff_path)
    proteins = parse_jsc1_proteins(protein_faa_path)

    resolved: list[ResolvedGene] = [
        _resolve_entry(entry, cds_list, proteins) for entry in catalog
    ]

    if hmm_path.exists() and not overwrite:
        if catalog_tsv_path.exists():
            return resolved

    alphabet = pyhmmer.easel.Alphabet.amino()
    background = pyhmmer.plan7.Background(alphabet)
    builder = pyhmmer.plan7.Builder(alphabet)

    hmms: list[pyhmmer.plan7.HMM] = []
    for gene in resolved:
        hmm_name = _hmm_name(gene.name)
        hmm = _build_single_sequence_hmm(
            hmm_name,
            gene.sequence,
            builder=builder,
            background=background,
            alphabet=alphabet,
        )
        hmms.append(hmm)

    with hmm_path.open("wb") as fh:
        for hmm in hmms:
            hmm.write(fh)

    with catalog_tsv_path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(CATALOG_TSV_HEADER) + "\n")
        for gene in resolved:
            fh.write(
                "\t".join(
                    [
                        _hmm_name(gene.name),
                        gene.name,
                        gene.role,
                        gene.protein_id,
                        gene.locus_tag,
                        gene.contig,
                        str(gene.start),
                        str(gene.end),
                        gene.strand,
                        str(len(gene.sequence)),
                        "true" if gene.in_cluster else "false",
                        gene.resolution,
                    ]
                )
                + "\n"
            )

    return resolved
