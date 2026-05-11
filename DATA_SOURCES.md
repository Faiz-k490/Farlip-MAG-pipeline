# Data Sources — FaRLiP Comparative Genomics Pipeline

This document is the canonical lock-in of every database, organism, and accession the
ingestion pipeline depends on. It is the human-readable counterpart to the machine-readable
configs in [`config/`](config/). When anything here changes, update the matching config file
*and* append a new entry to [`DEVLOG.md`](DEVLOG.md).

**Iteration:** 3 + Phase 1.5 (Phase 1 ingestion + Phase 1.5 SSU/taxonomic routing +
Phase 2 FaRLiP detection + Phase 3 QC + extremophile vault).
**Last updated:** 2026-05-10.

---

## 1. Primary FaRLiP Reference Anchor

| Field | Value |
|---|---|
| Organism | *Leptolyngbya* sp. strain JSC-1 |
| RefSeq assembly accession | `GCF_000733415.1` |
| GenBank assembly accession | `GCA_000733415.1` |
| NCBI BioProject | `PRJNA62993` |
| NCBI Taxonomy ID | `1547898` |
| Assembly size | ~7.87 Mb, 2 contigs, 6,056 genes |
| Source | NCBI RefSeq |
| Datasets v2 download URL | `https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000733415.1/download` |

**Why this organism.** JSC-1 is the original discovery organism for Far-Red Light
Photoacclimation. The canonical 21-gene FaRLiP cluster — including paralogs of *psbA, psbB,
psbC, psbD, psbH*, the red/far-red knotless phytochrome, and response regulators — was first
characterised here. It is the cleanest available evolutionary anchor for downstream cluster
alignment and gene-presence logic.

**Citations.**
- Gan, F., Zhang, S., Rockwell, N. C., Martin, S. S., Lagarias, J. C., & Bryant, D. A. (2014).
  Extensive remodeling of a cyanobacterial photosynthetic apparatus in far-red light.
  *Science*, 345(6202), 1312–1317. https://doi.org/10.1126/science.1256963
- Gan, F., Shen, G., & Bryant, D. A. (2015). Occurrence of Far-Red Light Photoacclimation
  (FaRLiP) in Diverse Cyanobacteria. *Life*, 5(1), 4–24. https://doi.org/10.3390/life5010004

**Secondary anchors (deferred to iteration 2).** *Chlorogloeopsis fritschii* PCC 9212
(`GCF_000317105.1`) introduces heterocyst-specific metabolic variables that would complicate
the iteration-1 structural mapping. The config schema accepts a list of secondary references
so this slots in later without a refactor.

---

## 2. Extremophile MAG Corpus (Test Dataset)

| Field | Value |
|---|---|
| Repository | ENA (European Nucleotide Archive) |
| BioProject | `PRJEB59431` |
| Secondary ID | `ERP144462` |
| Submitter | Finnish Environment Institute (CCAMBIO project) |
| First public | 2023-02-03 |
| Assembly count | 37 cyanobacterial / candidate-bacterial MAGs |
| Accession range | `GCA_949127635` … `GCA_949128045` |
| Enumeration call | `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB59431&result=assembly&fields=accession,sample_accession,scientific_name,assembly_level&format=tsv` |

**Why this BioProject.** It is the largest published genome-resolved metagenomics study of
polar cyanobacterial mats (Arctic, sub-Antarctic, Antarctic), and most MAGs are only
distantly related to previously sequenced cyanobacteria — exactly the evolutionary diversity
needed to test for FaRLiP cluster presence/absence under cold-stress conditions.

**Citation.**
- Pessi, I. S., Popin, R. V., Durieu, B., Lara, Y., Tytgat, B., Savaglia, V., … Wilmotte, A.
  (2023). Novel diversity of polar Cyanobacteria revealed by genome-resolved metagenomics.
  *Microbial Genomics*, 9(4), mgen001001. https://doi.org/10.1099/mgen.0.001001

**Notable lineages in the 37 MAGs.** *Pseudanabaena*, *Leptolyngbya*, *Microcoleus / Tychonema*,
*Phormidium*, *Crinalium*, *Chamaesiphon*, an enigmatic *Chroococcales* lineage distantly
related to *Microcystis*, and a cold-restricted *Gloeobacterales* lineage proposed as
*Candidatus* Sivonenia alaskensis.

---

## 3. Phylogenetic Database — SILVA

| Field | Value |
|---|---|
| Release | SILVA 138.2 |
| Release date | 2024-07-11 |
| Subset used | SSU Ref NR99 (16S, 99% redundancy filtered, ~510k sequences) |
| FASTA download URL | `https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz` |
| MD5 URL | `https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz.md5` |
| File size | ~150 MB (compressed FASTA) |

**Why SILVA 138.2.** This is the current production release as of 2026-05-10; no newer
release has been published. The NR99 subset is the standard for 16S phylogenetic placement
(eliminates near-duplicate sequences without losing diversity), and the FASTA export is
preferred over the binary ARB file because downstream tools in iterations 2–3 (HMMER,
mafft, IQ-TREE) consume FASTA directly.

**Citation.**
- Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., Peplies, J.,
  & Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: improved data
  processing and web-based tools. *Nucleic Acids Research*, 41(D1), D590–D596.
  https://doi.org/10.1093/nar/gks1219

---

## 4. Taxonomic Classification — GTDB

| Field | Value |
|---|---|
| Release | GTDB R226 |
| Bacterial taxonomy file | `bac120_taxonomy_r226.tsv.gz` |
| Bacterial metadata file | `bac120_metadata_r226.tsv.gz` |
| Taxonomy download URL | `https://data.gtdb.ecogenomic.org/releases/release226/226.0/bac120_taxonomy_r226.tsv.gz` |
| Metadata download URL | `https://data.gtdb.ecogenomic.org/releases/release226/226.0/bac120_metadata_r226.tsv.gz` |

**Why GTDB R226.** Explicitly requested by Dr. Barton. GTDB resolves long-standing
cyanobacterial classification problems by using genome-wide concatenated marker alignments
(`bac120`), giving more stable taxonomic placement than 16S alone for the polar MAGs.

**Citation.**
- Parks, D. H., Chuvochina, M., Rinke, C., Mussig, A. J., Chaumeil, P.-A., & Hugenholtz, P.
  (2022). GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically
  consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Research*,
  50(D1), D785–D794. https://doi.org/10.1093/nar/gkab776

---

## 5. Cave & Additional Arctic Genomes — TBD

Placeholder list in [`config/extremophiles.json`](config/extremophiles.json) under
`cave_arctic_pending: []`. Dr. Barton will supply the exact NCBI / ENA accessions; the
pipeline is config-driven, so populating this list is a one-line edit — no code changes.

**Action item.** Bring to next meeting with Dr. Barton: *"Please provide the exact cave +
additional Arctic genome accessions (NCBI BioProject, ENA BioProject, or individual
assembly accessions) so we can populate `cave_arctic_pending` and re-run the ingest."*

---

## 6. Verification Checklist

Before declaring iteration 1 complete, every entry above must satisfy:

- [ ] Listed in the matching `config/` file with the exact URL recorded here.
- [ ] Successfully downloaded by `python -m src.pipeline fetch all`.
- [ ] A row appears in `data/manifest.tsv` with a recorded SHA-256.
- [ ] The MD5 / checksum (where the upstream provides one) is verified.
- [ ] The relevant `DEVLOG.md` entry references this document.

---

## 7. FaRLiP Gene Catalog (Phase 2)

The 21-gene catalog driving HMM profile construction is committed at
[`config/farlip_genes.yaml`](config/farlip_genes.yaml). Each entry maps a FaRLiP
gene name to its exact JSC-1 locus tag, RefSeq protein accession (`WP_*`),
gene_tag (where annotated), and a product-substring fallback.

**Source.** Gan, F., Zhang, S., Rockwell, N. C., Martin, S. S., Lagarias, J. C.,
& Bryant, D. A. (2014). Extensive remodeling of a cyanobacterial photosynthetic
apparatus in far-red light. *Science*, 345(6202), 1312-1317.
https://doi.org/10.1126/science.1256963

**JSC-1 cluster window.** Verified against the JSC-1 RefSeq GFF
(`GCF_000733415.1`, downloaded 2026-05-10):

- Contig: `NZ_KL662192.1`
- Start (1-based, inclusive): `1,095,910`
- End   (1-based, inclusive): `1,121,582`
- Span: 25,672 bp (matches the ~25-30 kb cluster size reported by Gan 2014)
- 21 consecutive non-pseudo CDS rows, all resolved by exact locus_tag.

**Resolved 21-gene cluster (cluster order, 5' to 3' on `NZ_KL662192.1`):**

| # | Gene | Locus tag | Protein ID | Role |
|---|---|---|---|---|
| 1 | psbA2 | OSCCY1_RS30350 | WP_036011370.1 | PSII D1 paralog |
| 2 | apcD2 | OSCCY1_RS30355 | WP_036011372.1 | Allophycocyanin paralog |
| 3 | apcD3 | OSCCY1_RS30360 | WP_036011375.1 | Allophycocyanin paralog |
| 4 | apcB2 | OSCCY1_RS30365 | WP_036011378.1 | Allophycocyanin beta paralog |
| 5 | cluster_hypothetical_RS30370 | OSCCY1_RS30370 | WP_051926078.1 | Small cluster ORF (hypothetical) |
| 6 | apcE2 | OSCCY1_RS30375 | WP_051926079.1 | Phycobilisome core-membrane linker LCM |
| 7 | apcD4 | OSCCY1_RS30380 | WP_036011381.1 | Allophycocyanin paralog |
| 8 | psbD3 | OSCCY1_RS30385 | WP_036011384.1 | PSII D2 paralog |
| 9 | psbC2 | OSCCY1_RS30390 | WP_036011388.1 | PSII CP43 paralog |
| 10 | psbB2 | OSCCY1_RS30395 | WP_036011391.1 | PSII CP47 paralog |
| 11 | psbH2 | OSCCY1_RS33520 | WP_063834837.1 | PSII PsbH paralog |
| 12 | rfpC | OSCCY1_RS35750 | WP_036011392.1 | Response regulator (cluster regulator 1) |
| 13 | rfpA | OSCCY1_RS30410 | WP_051926080.1 | Knotless phytochrome (master regulator) |
| 14 | rfpB | OSCCY1_RS30415 | WP_036011395.1 | Response regulator (cluster regulator 2) |
| 15 | psbA3 | OSCCY1_RS30420 | WP_081915074.1 | PSII D1 paralog (second copy) |
| 16 | psaA2 | OSCCY1_RS30425 | WP_036011398.1 | PSI core protein PsaA paralog |
| 17 | psaB2 | OSCCY1_RS30430 | WP_036011400.1 | PSI core protein PsaB paralog |
| 18 | psaL2 | OSCCY1_RS30435 | WP_036011402.1 | PSI subunit XI paralog |
| 19 | psaI2 | OSCCY1_RS30440 | WP_036011405.1 | PSI subunit VIII paralog |
| 20 | psaF2 | OSCCY1_RS30445 | WP_036012501.1 | PSI subunit III paralog |
| 21 | psaJ2 | OSCCY1_RS30450 | WP_036011408.1 | PSI subunit IX paralog |

The catalog itself remains the authoritative source - this table is a
human-readable echo for reviewers; if the YAML and this table disagree, the YAML
wins and this table should be regenerated.

**Phase 2.1 plan.** Chlorogloeopsis fritschii PCC 9212 (`GCF_000019485.1`) will
be added as `anchors:` per gene in `config/farlip_genes.yaml`. Combined with
mafft (via subprocess) the per-gene HMMs become true MSA-based profiles, which
will substantially improve sensitivity on divergent FaRLiP paralogs in
fragmented MAGs.

---

## 8. Single-Copy Marker HMMs for Genome Quality Assessment (Phase 3)

The Phase 3 quality-assessment pipeline (`src.qc_assessor`) needs a bacterial
single-copy ortholog (SCO) marker HMM database to score MAG completeness and
contamination on the MIMAG quality tiers (Bowers et al. 2017,
*Nat Biotechnol* 35:725 - High-quality: completeness >= 90% and contamination
<= 5%; Medium-quality: completeness >= 50% and contamination < 10%).

The marker-set definitions are committed at
[`config/qc_markers.yaml`](config/qc_markers.yaml). The fetcher is
[`src/qc_markers_fetcher.py`](src/qc_markers_fetcher.py); when the assessor
runs, it concatenates every `.hmm` it found into a single
`data/qc_markers/<key>.hmm` file (the format `pyhmmer.plan7.HMMFile` reads in
one open call) and copies the BUSCO `scores_cutoff` sidecar to
`data/qc_markers/<key>.scores_cutoff.tsv` for per-HMM bit-score gating.

| Field | Value |
|---|---|
| Default marker set | `BUSCO_bacteria_odb10` |
| Source URL | <https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2024-01-08.tar.gz> |
| Tarball SHA-256 | `a8a625ee04510c33fb42d39fd730a82b7cf93df5894f4566e7fe05dff75b45cd` |
| Compressed size | ~6 MB |
| Number of HMMs | 124 single-copy bacterial markers (built against OrthoDB 10) |
| Per-HMM cutoffs | Yes (`bacteria_odb10/scores_cutoff` inside the archive) |
| Citation | Manni et al. 2021, *Mol Biol Evol* 38:4647 (BUSCO v5) |

**GTDB bac120 escape hatch (`gtdb_bac120_local`).** GTDB R226 ships its
bac120 markers only inside the full GTDB-Tk reference package (~141 GB at
R226). That is too large to require for a QC step, so the fetcher refuses to
auto-download anything past `MAX_AUTO_DOWNLOAD_BYTES = 200 MB`. Users who
already have the GTDB-Tk package extracted can flip the YAML to point
`local_marker_dir` at their `markers/pfam_and_tigrfam` directory; the
assessor then concatenates that directory exactly the same way.

**Validation against JSC-1.** BUSCO `bacteria_odb10` is a generic bacterial
marker set; cyanobacteria carry multiple paralogs of photosystem and
ribosomal genes that BUSCO scores as duplicates even after applying the
per-HMM `scores_cutoff` gate. JSC-1 (`GCF_000733415.1`) scores 94.4% / 8.1%
against this marker set (MIMAG MEDIUM), which sits inside the published
BUSCO-vs-CheckM2 envelope for cyanobacterial RefSeq assemblies. Users who
need the strict MIMAG HIGH tier should run the optional
`scripts/run_checkm2.sh` cross-validation path (see DEVLOG entry for
2026-05-10 Phase 3).

**MIMAG threshold citation.** Bowers, R., Kyrpides, N., Stepanauskas, R.,
Harmon-Smith, M., Doud, D., Reddy, T., et al. 2017. Minimum information about
a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG)
of bacteria and archaea. *Nat Biotechnol* 35:725-731.

---

## 9. Phase 1.5 - SSU Extraction and Taxonomic-Placement Routing

Phase 1.5 wires the Phase 1 manifest into species-level taxonomic placement.
Every assembly is scanned for 16S rRNA with `barrnap`; MAGs that yield at
least one 16S record are bundled for SILVA SINA submission, MAGs that do
not are routed to a generated GTDB-Tk shell script for genome-wide marker
classification. Phase 1.5 only prepares inputs - it never contacts SILVA
SINA or runs GTDB-Tk itself (the latter requires the ~141 GB GTDB R226
reference data and is impractical to bundle into the pipeline).

### 9.1 `barrnap` - rRNA detection

| Field | Value |
|---|---|
| Tool | `barrnap` (BAsic Rapid Ribosomal RNA Predictor) |
| Upstream | <https://github.com/tseemann/barrnap> |
| Default version | `0.9` (any 0.9+ release on PATH is accepted) |
| Default kingdom | `bac` (bacterial 16S/23S/5S HMM library bundled with barrnap) |
| Default threads | 4 |
| Default `--evalue` | `1e-6` (matches barrnap's own default) |
| Default `--reject` | `0.25` (matches barrnap's own default) |
| Distribution | conda-forge / bioconda; not auto-installed by this pipeline |

**Why barrnap.** It uses a curated nhmmer / pyhmmer HMM library for the
three rRNA classes and ships with sensible defaults for bacterial MAGs.
Compared to RNAmmer (license-restricted) and Infernal+Rfam (much heavier
configuration), barrnap is the lowest-friction option for the assignment.

**Citation.** Seemann, T. (2013). barrnap: BAsic Rapid Ribosomal RNA
Predictor. <https://github.com/tseemann/barrnap>

### 9.2 SILVA SINA submission (deferred to user)

| Field | Value |
|---|---|
| Tool | SILVA Incremental Aligner (SINA) |
| Web upload | <https://www.arb-silva.de/aligner/> |
| Reference DB | SILVA 138.2 NR99 (already on disk at `data/silva/138.2/`) |
| Input | `data/ssu/silva_sina_input.fasta` (multi-FASTA, built by `bundle-silva`) |

The pipeline's job is to produce the bundle. The user runs SINA (either via
the web upload or the SINA CLI against the on-disk SILVA 138.2 NR99 FASTA)
and feeds the species-level placements back into Phase 4.

**Citation.** Pruesse, E., Peplies, J., & Glöckner, F. O. (2012). SINA:
Accurate high-throughput multiple sequence alignment of ribosomal RNA
genes. *Bioinformatics* 28:1823-1829. <https://doi.org/10.1093/bioinformatics/bts252>

### 9.3 GTDB-Tk classify_wf (deferred to user, generator emits the script)

| Field | Value |
|---|---|
| Tool | `gtdbtk classify_wf` |
| Upstream | <https://github.com/Ecogenomics/GTDBTk> |
| Reference data | GTDB R226 GTDB-Tk reference package (~141 GB; user-provided via `GTDBTK_DATA_PATH`) |
| Input | `data/ssu/gtdbtk_batchfile.tsv` (assembly path + accession; built by `gen-gtdbtk-script`) |
| Output | `data/gtdbtk_results/gtdbtk.bac120.summary.tsv` |
| Default threads | 8 (overridable via `--threads`) |

The generator `src.gtdbtk_shellscript_generator` writes `scripts/run_gtdbtk.sh`
with `set -euo pipefail`, requires `GTDBTK_DATA_PATH` to be exported, and
materialises every gzipped ENA assembly into `data/ssu/gtdbtk_input/<acc>.fna`
so a single `--extension fna` argument works for the whole batch.

**Citation.** Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P., & Parks, D. H.
(2022). GTDB-Tk v2: memory-friendly classification with the Genome Taxonomy
Database. *Bioinformatics* 38:5315-5316.
<https://doi.org/10.1093/bioinformatics/btac672>

### 9.4 Routing log schema

`data/ssu/routing_log.tsv` is the single source of truth for Phase 1.5
decisions. Columns (`src.ssu_extractor.ROUTING_LOG_COLUMNS`):

```
accession  source  fasta_path  n_contigs_scanned  n_16S_predicted
decision   output_path  barrnap_version  notes
```

`decision` is one of `silva` / `gtdbtk`. The Phase 3
`extremophile_vault_exporter` joins on `accession` to set the
`Taxonomic_Placement` frontmatter slot in
`data/extremophile_vault/<accession>.md`.

---

## 10. Phase 3 Extremophile Vault Export

`data/extremophile_vault/` is rendered by
`src/extremophile_vault_exporter.py` from `data/qc_filtered_mags.tsv` (the
`recommended_for_phylogeny == True` view) plus, optionally, the Phase 1.5
routing log. Each per-MAG note carries exactly the five-key frontmatter
mandated by the assignment:

```yaml
---
Accession: GCA_949127XYZ
Completeness: 92.40
Contamination: 1.20
Taxonomic_Placement: SILVA SINA: pending
FaRLiP_Cluster_Present:
---
```

`FaRLiP_Cluster_Present` is intentionally **blank** in the frontmatter -
filled in by manual review after consulting the Phase 2 vault note for the
same accession (`data/obsidian_export/<accession>.md`). The master index
`High_Quality_Extremophile_Genomes.md` wikilinks to every per-MAG note and
groups by MIMAG tier when multiple tiers are present.
