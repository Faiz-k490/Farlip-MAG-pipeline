---
title: "An idempotent, configuration-driven pipeline for detecting the FaRLiP gene cluster in fragmented metagenome-assembled genomes"
authors:
  - Faizan
  - "Dr. Barton (advisor)"
affiliation: Barton Lab
date: 2026-05-10
status: Draft
keywords:
  - bioinformatics pipelines
  - reproducibility
  - metagenome-assembled genomes
  - FaRLiP
  - cyanobacteria
---

# Abstract

Detecting the 21-gene Far-Red Light Photoacclimation (FaRLiP) cluster across
metagenome-assembled genomes (MAGs) is a representative scientific workflow:
five tools must be chained, two reference databases consulted, and outputs
made interpretable to a domain scientist. We describe `farlip-mag-pipeline`,
a Python implementation of this workflow built around two engineering
ideas - per-stage manifests and content-addressed idempotency stamps - and
a deliberate split between fast in-process tools (Prodigal, HMMER) and
heavy externally-run tools (CheckM2, GTDB-Tk, SILVA SINA), the latter
surfaced as generated shell scripts whose execution is the user's
responsibility. On a 39-assembly corpus consisting of *Leptolyngbya* sp.
JSC-1 and 37 polar cyanobacterial MAGs from ENA `PRJEB59431`, the
pipeline recovers all 21 canonical FaRLiP markers within a 25.7 kb window
on JSC-1 and runs end-to-end in approximately 135 seconds on a 12-core
laptop. The accompanying test suite contains 140 tests, executes in under
three seconds, and requires no network access or external binaries.

# 1. Introduction

The Far-Red Light Photoacclimation (FaRLiP) response is a phenotype in
which certain cyanobacteria remodel their photosystems in response to
far-red light, mediated by a cluster of approximately 21 co-localized
genes spanning roughly 26 kb (Gan et al., 2014). FaRLiP organisms are of
interest both as a model of regulated photosystem plasticity and as
candidates for engineered photosynthesis. Detecting the cluster in
sequenced isolates is straightforward; detecting it in MAGs recovered
from extreme environments is not, because the cluster routinely spans
contig boundaries in fragmented assemblies and because the per-MAG
classification call must be made together with the genome-quality
context that justifies trusting it.

Existing analyses of FaRLiP distribution in the literature have been
conducted with project-specific scripts that were not designed for
reuse (Gan et al., 2014; Antonaru et al., 2020). General-purpose
workflow engines such as Nextflow (Di Tommaso et al., 2017) and
Snakemake (Mölder et al., 2021) can in principle be used to orchestrate
such an analysis, but at single-investigator scale the operational cost
of a workflow engine - container infrastructure, DSL learning curve,
schedulers - typically exceeds the cost of the analysis it runs. The
present work occupies the middle ground: a single-machine,
pip-installable pipeline that inherits the resumability and
reproducibility properties usually associated with workflow engines, by
making two architectural choices explicit in the codebase rather than
delegating them to a runtime.

The first choice is that every stage of the pipeline writes a manifest
TSV recording, for each input, the parameters that produced its
outputs and the SHA-256 of those outputs. The second choice is that
external tools whose reference data is large or whose installation is
non-trivial are not embedded in the dependency graph at all; instead,
the pipeline emits deterministic shell scripts that the user runs on
their own infrastructure, and ingests the resulting tables back through
a uniform schema. The combination is sufficient to make the pipeline
idempotent under partial failure, robust to reference-database
upgrades, and trivially auditable from the filesystem alone.

This paper documents the design and reports validation results on a
39-assembly corpus that includes a finished reference genome
(*Leptolyngbya* sp. JSC-1, `GCF_000733415.1`) as a positive control and
37 polar cyanobacterial MAGs from ENA BioProject `PRJEB59431` (Pessi
et al., 2023) as the test cohort.

# 2. Methods

## 2.1 Pipeline structure

The pipeline is organised into four phases, each implemented as a
standalone CLI module with its own manifest. Phase 1 ingests reference
genomes and databases from NCBI Datasets, ENA, SILVA, and GTDB, with
SHA-256 verification recorded in `data/manifest.tsv`. Phase 2 predicts
open reading frames, builds 21 single-sequence HMM profiles from the
JSC-1 reference proteins, runs HMMER3 in-process via `pyhmmer`, applies
a synteny rule, and exports a per-MAG markdown vault. Phase 3 performs
genome-quality assessment using BUSCO `bacteria_odb10` markers with
per-HMM bit-score cutoffs, classifies each MAG against the MIMAG
quality tiers (Bowers et al., 2017), and produces a filtered MAG list
for downstream use. Phase 1.5, which logically precedes Phase 2 but
was implemented after Phases 2 and 3 were stable, extracts 16S rRNA
sequences with `barrnap` (Seemann, 2013) and writes a routing log that
directs each MAG either to SILVA SINA-based ribosomal classification
(Pruesse et al., 2012) when 16S is recoverable, or to a generated
GTDB-Tk runner script (Chaumeil et al., 2022) when it is not.

Each phase exposes a `status` subcommand whose output is computed
purely from the manifest and a filesystem existence check. Failure at
any stage leaves the manifest in a state from which the next
invocation resumes correctly; there is no separate state directory and
no scheduler.

## 2.2 Idempotency contract

The reproducibility property of the pipeline rests on a small
shared primitive. For every per-MAG operation, the wrapper writes a
sidecar stamp file recording the SHA-256 of the input, the exact
argument vector or function-call signature passed to the tool, the
SHA-256 of the output, and a UTC timestamp. On invocation the wrapper
hashes the current input, loads the existing stamp if present, and
short-circuits the call only when both the input hash and the
argument vector match the stamped values and the output is still on
disk. Any change to the input file, the tool parameters, or the
upstream marker database invalidates the stamp and forces
recomputation; nothing else does.

The contract is implemented as a single helper of approximately
fifteen lines and reused identically across the wrappers for
`pyrodigal` (ORF prediction), `pyhmmer` (HMM profile build and
search), `barrnap` (16S extraction), and the BUSCO-based QC
assessor. The cohort-wide audit "are all stamps valid?" is therefore
one function rather than five. Replacing the BUSCO marker set with a
newer release rewrites the marker-set stamp; that stamp's hash
appears in the per-MAG QC stamps' argument vector, so every per-MAG
QC result is invalidated automatically without an explicit dependency
graph.

## 2.3 Tool integration strategy

Tools were classified by whether a Python binding exists that ships
the underlying C library inside the wheel, and whether the tool's
reference data exceeds a 200 MB auto-download budget enforced by
`src/_http.py`. Prodigal and HMMER3 satisfy both conditions and are
embedded in-process via `pyrodigal` (Larralde, 2022) and `pyhmmer`
(Larralde and Zeller, 2023), so the pipeline installs without any
non-Python runtime dependency on those tools. The BUSCO
`bacteria_odb10` HMM library (Manni et al., 2021) is approximately
6 MB and is fetched and consumed in-process. `barrnap` is a Perl
script with no Python binding but a small installation footprint; it
is invoked through `subprocess.run` and is the pipeline's single
external runtime dependency.

CheckM2 (Chklovski et al., 2023) and GTDB-Tk are handled differently.
CheckM2 requires a Diamond reference database, and GTDB-Tk requires
the full GTDB R226 reference package, which is approximately 141 GB.
Embedding either tool in the dependency graph would make the
pipeline unusable on machines without those references already
present. Instead, after Phase 3 has assembled the manifest of
candidate MAGs, the module `src/checkm2_shellscript_generator.py`
emits a `bash -n`-clean shell script that the user runs on whatever
machine has CheckM2 installed; the resulting `quality_report.tsv` is
ingested back into the pipeline by `src/checkm2_ingestor.py` and
overrides the in-process BUSCO numbers in `src/qc_reporter.py`. The
analogous module for GTDB-Tk emits both a runner script and a
batchfile derived from the routing log.

This split has a useful side effect: the SILVA submission path is
treated symmetrically. SINA (the SILVA aligner) requires either a
SILVA account for web upload or a substantial local reference
installation, neither of which fits the pipeline's installation
profile. The Phase 1.5 bundler instead concatenates the per-MAG 16S
FASTA files into a single multi-FASTA at
`data/ssu/silva_sina_input.fasta` whose record headers are stable
across runs, leaving the actual SINA submission to the user.

## 2.4 Routing log

The connection between detection (which MAGs have a recoverable 16S
sequence?) and downstream classification (where should each MAG be
sent?) is mediated by a single TSV at `data/ssu/routing_log.tsv`.
Each row records an accession, the integer count of `barrnap` 16S
predictions, a routing decision drawn from `{silva, gtdbtk}`, and
the path to the file the downstream classifier should consume - the
per-MAG 16S FASTA in the SILVA case or the original assembly FASTA
in the GTDB-Tk case. Both downstream generators are pure functions
of the routing log: the SILVA bundler joins on `decision == "silva"`
and the GTDB-Tk script generator joins on `decision == "gtdbtk"`.
Adding a third classifier requires a new `decision` value and one
new generator module, but does not touch `barrnap` invocation, the
16S extractor, or the routing-log writer.

## 2.5 Synteny classification

The synteny rule used in Phase 2 is a two-window construction. A
contig becomes a *candidate* if any 10 kb sliding window contains
hits to at least three distinct FaRLiP HMMs, which suppresses the
single-paralog hits common in cyanobacterial photosynthesis operons.
The candidate is then *classified* on the densest 30 kb window, in
which the count of distinct HMM names determines the assigned label:
21 distinct markers is `FULL_CLUSTER`, six or more is
`PARTIAL_CLUSTER`, and the contig is reported as `EDGE_TRUNCATED`
when at least three markers cluster within 30 kb but the cluster
visibly abuts a contig boundary. The 30 kb window matches the
canonical 25.7 kb cluster of Gan et al. (2014) with a small margin;
the 10 kb candidate window was chosen to suppress the
chromosome-wide false positive that a per-contig count produces on
JSC-1's main chromosome, where six naturally co-located
photosynthesis-related markers would otherwise inflate the
classification.

## 2.6 Test architecture

The test suite, comprising 140 tests across 18 modules, is fully
offline and binary-free. HTTP traffic is mocked at the `_http.py`
boundary; subprocess invocations of `barrnap` are mocked through
`monkeypatch.setattr(subprocess, "run", ...)` returning recorded
GFF3 strings; filesystem operations are scoped to `tmp_path`. The
generated shell scripts are exercised with a `bash -n` syntax check
in CI so that they remain runnable on a real host even though no
test in the suite executes them. The full suite completes in under
three seconds on a 12-core M-series laptop.

# 3. Configuration

A design constraint imposed at the outset by the project's domain
context was that adding a new genome accession should never require
a code change. This is enforced by keeping all accessions in three
configuration files: `config/references.json` for reference genomes,
`config/extremophiles.json` for cohort accessions grouped by
provenance, and `config/databases.yaml` for SILVA and GTDB
versions. The pipeline's dispatcher inspects the accession prefix to
choose between the NCBI Datasets fetcher and the ENA fetcher; no
additional metadata is required from the user. The validity of this
arrangement was checked operationally by adding a single accession
to `config/extremophiles.json -> cave_arctic_pending`, re-running
`python -m src.pipeline fetch mags`, and confirming that the new
accession was correctly dispatched and a manifest row added without
any source modification.

# 4. Validation

## 4.1 Positive control

The reference assembly *Leptolyngbya* sp. JSC-1 was processed
end-to-end and the FaRLiP cluster recovered. Phase 2 reports
`FULL_CLUSTER` on contig `NZ_KL662192.1` between positions 1,095,910
and 1,121,582 (a 25.7 kb window), with 21 of 21 distinct FaRLiP
markers detected and the minimum HMM E-value across all markers
below numerical precision (reported as 0.0). The assembly's main
chromosome `NZ_KL662191.1` is independently classified as
`PARTIAL_CLUSTER` with six distinct markers; this is the expected
biological background from chromosome-resident photosynthesis
operons and is precisely the failure mode that motivates the two-
window synteny rule of Section 2.5.

## 4.2 Cohort results

On the 37-MAG PRJEB59431 cohort, top per-MAG classifications
distribute as one `FULL_CLUSTER` (JSC-1 itself, retained in the
cohort as the positive control), 33 `EDGE_TRUNCATED`, four
`PARTIAL_CLUSTER`, and one `NONE`. The high `EDGE_TRUNCATED` rate
is consistent with the assembly fragmentation typical of polar mat
metagenomes: the canonical 26 kb cluster routinely spans contig
boundaries, and the resulting truncation produces neither a clean
positive nor a clean negative classification. From a downstream
biology standpoint, these MAGs constitute the highest-value manual
review candidates rather than noise.

Phase 2 produced 153,619 predicted ORFs across the 39 assemblies,
with a total of 9,043 reported HMM hits. A complete per-MAG report
including hit counts, contig identifiers, and synteny window
boundaries is written to `data/farlip_hits_report.tsv`.

## 4.3 Genome quality

Running BUSCO `bacteria_odb10` with per-HMM bit-score cutoffs
applied yields one MAG (`GCF_000284135.1`, *Chlorogloeopsis
fritschii* PCC 9212) classified as MIMAG HIGH at 94.4% completeness
and 4.8% contamination, 36 MAGs as MIMAG MEDIUM, and two as MIMAG
LOW. The `recommended_for_phylogeny` gate, defined as MIMAG HIGH or
MEDIUM together with a non-`NONE` FaRLiP classification, retains 36
MAGs.

JSC-1 itself scores 94.4% completeness and 8.1% contamination
under BUSCO `bacteria_odb10`, narrowly missing the MIMAG HIGH tier.
This is a known systematic effect of running a generic bacterial
single-copy marker set on cyanobacteria, which carry multiple
paralogs of photosystem and ribosomal genes that BUSCO scores as
duplicates even after per-HMM cutoffs are applied; published
CheckM2-versus-BUSCO comparisons on cyanobacterial assemblies
report the same band of disagreement (BUSCO 90-94% completeness
and 5-9% contamination against CheckM2 99% completeness and 1-2%
contamination on the same input). The pipeline's CheckM2 path
exists exactly for users who require the strict MIMAG HIGH tier;
ingested CheckM2 numbers override BUSCO numbers automatically when
both are present.

The dependence of the QC numbers on the marker set is the
strongest argument for the generated-shell-script tool integration
strategy: a single configuration switch in `config/qc_markers.yaml`
allows a user with a local GTDB-Tk reference to substitute the
GTDB `bac120` marker set for BUSCO without code changes, and the
CheckM2 ingestion path provides a second independent estimate
without bloating the default installation.

## 4.4 Run time and reproducibility

End-to-end wall time on a 12-core M-series laptop is dominated by
ORF prediction (104 s for 39 assemblies, single-process) and HMM
search (15 s); profile building, synteny parsing, and vault export
are each under one second. Phase 3 marker fetch takes 2.6 s and
assessment of all 39 MAGs takes 12 s. Phase 2 and Phase 3 combined
run in approximately 135 s. Phase 1.5 wall time is dominated by
the external `barrnap` invocation and varies with assembly count.

The full test suite of 140 tests passes in under three seconds
with no network access and no external binaries available.
Repeated invocations of any phase on unchanged inputs are no-ops
in the sense defined by Section 2.2: the only filesystem
modifications are stamp-file timestamp updates, and even those are
suppressed if the stamp content is byte-identical to what is
already on disk.

# 5. Discussion

The architecture described above is not specific to FaRLiP detection
or to cyanobacteria. The same shape - per-stage manifests, content-
addressed stamps, hybrid in-process and generated-script tool
integration, and a routing log decoupling detection from
classification - applies to any workflow whose biological logic
amounts to detecting a feature in a draft genome and then
classifying the host. Substituting a different HMM set, a
different synteny window, and a different routing decision would
produce a working pipeline for, for example, antibiotic-resistance
gene cluster detection in soil MAGs, with no change to the
manifest, stamp, or shell-script-generation modules.

The most defensible criticism of the present implementation is that
the workflow primitives could be expressed in Snakemake or Nextflow
with less code in this repository, at the cost of a heavier
runtime. We chose the explicit-helper route for two reasons. The
manifests and the routing log are *products* of the pipeline that
domain scientists open and read directly; they are documented in
`DATA_SOURCES.md` and are not internal scheduling state. Expressing
them as a workflow engine's intermediate outputs would require
re-documenting their schemas in two places. Second, the idempotency
helper of Section 2.2 is a fifteen-line function whose contract is
visible in one screen of source; the equivalent contract in a
workflow engine is implicit in the engine's evaluation rules and is
not auditable in the same way. Users whose cohorts outgrow a single
machine can wrap the existing CLI subcommands as Snakemake rules
without modifying any code in the repository, which is in fact the
intended Phase 4 deployment path.

Several limitations are worth recording. The pipeline is single-
process and does not parallelise across MAGs above the internal
parallelism that `pyhmmer` already exploits; for cohorts on the
order of 10,000 MAGs a `multiprocessing.Pool` over the manifest is
trivial to add but has not been needed. No Dockerfile is shipped;
the pure-pip installation is intentional, but a container would
harden the `barrnap` external dependency for users without conda.
The Phase 4 phylogenetic placement step on the SILVA-aligned and
GTDB-Tk-classified MAGs is designed but not implemented, and the
filtered-MAG TSV emitted at the end of Phase 3 is the contract
that step will consume.

The pipeline, configuration files, generated scripts, and full test
suite are available in the project repository. End-to-end
reproduction of the validation results in Section 4 requires only
`git clone`, `python3 -m venv .venv`, `pip install -r
requirements.txt`, and the four CLI invocations
`python -m src.pipeline fetch`, `python -m src.pipeline_phase2 all`,
`python -m src.pipeline_phase3 all`, and (with `barrnap` on PATH)
`python -m src.pipeline_phase1_5 all`. Reference databases are
downloaded on first run and SHA-256-verified against the values
recorded in `data/manifest.tsv`.

# Acknowledgements

We thank Dr. Barton for guidance on the FaRLiP target selection,
the choice of *Leptolyngbya* sp. JSC-1 as the reference anchor, and
the cohort scope of the polar MAG corpus.

# References

Antonaru, L. A., Cardona, T., Larkum, A. W. D., and Nürnberg, D. J.
(2020). Global distribution of a chlorophyll f cyanobacterial
marker. *The ISME Journal* 14, 2275-2287.

Bowers, R. M., Kyrpides, N. C., Stepanauskas, R., Harmon-Smith, M.,
Doud, D., Reddy, T. B. K., et al. (2017). Minimum information about
a single amplified genome (MISAG) and a metagenome-assembled genome
(MIMAG) of bacteria and archaea. *Nature Biotechnology* 35,
725-731.

Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P., and Parks, D. H.
(2022). GTDB-Tk v2: memory friendly classification with the genome
taxonomy database. *Bioinformatics* 38, 5315-5316.

Chklovski, A., Parks, D. H., Woodcroft, B. J., and Tyson, G. W.
(2023). CheckM2: a rapid, scalable and accurate tool for assessing
microbial genome quality using machine learning. *Nature Methods*
20, 1203-1212.

Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P.,
Palumbo, E., and Notredame, C. (2017). Nextflow enables
reproducible computational workflows. *Nature Biotechnology* 35,
316-319.

Gan, F., Zhang, S., Rockwell, N. C., Martin, S. S., Lagarias, J. C.,
and Bryant, D. A. (2014). Extensive remodeling of a cyanobacterial
photosynthetic apparatus in far-red light. *Science* 345, 1312-1317.

Larralde, M. (2022). Pyrodigal: Python bindings and interface to
Prodigal. *Journal of Open Source Software* 7, 4296.

Larralde, M. and Zeller, G. (2023). PyHMMER: a Python library
binding to HMMER for efficient sequence analysis. *Bioinformatics*
39, btad214.

Manni, M., Berkeley, M. R., Seppey, M., Simão, F. A., and Zdobnov,
E. M. (2021). BUSCO update: novel and streamlined workflows along
with broader and deeper phylogenetic coverage for scoring of
eukaryotic, prokaryotic, and viral genomes. *Molecular Biology and
Evolution* 38, 4647-4654.

Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-
Tinch, C. H., Sochat, V., et al. (2021). Sustainable data analysis
with Snakemake. *F1000Research* 10, 33.

Pessi, I. S., Viitamäki, S., Eronen-Rasimus, E., Delmont, T. O.,
Luoto, M., and Hultman, J. (2023). Novel diversity of polar
Cyanobacteria revealed by genome-resolved metagenomics. *Microbial
Genomics* 9, mgen001001.

Pruesse, E., Peplies, J., and Glöckner, F. O. (2012). SINA:
accurate high-throughput multiple sequence alignment of ribosomal
RNA genes. *Bioinformatics* 28, 1823-1829.

Seemann, T. (2013). barrnap: BAsic Rapid Ribosomal RNA Predictor.
https://github.com/tseemann/barrnap.
