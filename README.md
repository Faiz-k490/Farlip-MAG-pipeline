# FaRLiP Comparative Genomics Pipeline

A config-driven pipeline that locks in the canonical 21-gene Far-Red Light
Photoacclimation (FaRLiP) cluster organism (*Leptolyngbya* sp. JSC-1) plus a polar
cyanobacterial MAG corpus, predicts ORFs in every assembly, builds HMM profiles
for the 21 FaRLiP genes, runs synteny-aware HMM searches across the MAGs, and
exports a browsable Obsidian markdown vault of candidates - ready to ingest the
cave/Arctic genomes Dr. Barton will provide.

**Phase 1** (ingestion), **Phase 1.5** (16S extraction + SILVA/GTDB-Tk
taxonomic-placement routing + extremophile vault), **Phase 2** (FaRLiP
detection), and **Phase 3** (genome quality assessment) are all shipped. See
[`DEVLOG.md`](DEVLOG.md) for the change history and
[`DATA_SOURCES.md`](DATA_SOURCES.md) for the canonical lock-in spec.

| | |
|---|---|
| **Primary anchor** | *Leptolyngbya* sp. JSC-1 (`GCF_000733415.1`) |
| **MAG test corpus** | ENA `PRJEB59431` (37 polar cyanobacterial MAGs) |
| **Phylogenetic DB** | SILVA 138.2 SSU Ref NR99 |
| **Taxonomy DB** | GTDB R226 |
| **Cave / Arctic** | TBD - awaiting Dr. Barton |

Full lock-in spec lives in [`DATA_SOURCES.md`](DATA_SOURCES.md). Every change is
journalled in [`DEVLOG.md`](DEVLOG.md).

---

## Quickstart

```bash
# 1. Set up environment
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 2. Tell NCBI who you are (their usage policy)
export NCBI_EMAIL="you@your-institution.edu"
# Optional: 10 req/s instead of 3 req/s
export NCBI_API_KEY="..."

# 3. Phase 1: ingest references + MAGs + databases (~450 MB)
python -m src.pipeline fetch all
python -m src.pipeline status
python -m src.pipeline verify

# 4. Phase 2: run FaRLiP detection on everything that landed in the manifest
python -m src.pipeline_phase2 all
python -m src.pipeline_phase2 status

# 5. Phase 3: score every assembly for completeness + contamination + N50,
#    refresh the Obsidian vault with QC frontmatter, and emit the high-quality
#    extremophile vault (data/extremophile_vault/).
python -m src.pipeline_phase3 all
python -m src.pipeline_phase3 status

# 6. Phase 1.5: 16S extraction (requires barrnap on PATH) + SILVA SINA bundle
#    + scripts/run_gtdbtk.sh for MAGs without recoverable 16S. SILVA SINA and
#    GTDB-Tk themselves are run by the user (see scripts/run_gtdbtk.sh).
python -m src.pipeline_phase1_5 all
python -m src.pipeline_phase1_5 status

# 7. Re-render the extremophile vault so Taxonomic_Placement frontmatter is
#    routing-aware now that data/ssu/routing_log.tsv exists.
python -m src.pipeline_phase3 export-extremophile-vault
```

Or use the wrapper:

```bash
./scripts/fetch_all.sh                    # same as `fetch all`
./scripts/fetch_all.sh --limit 3          # smoke test: 3 MAGs per BioProject
```

Output lands under `data/` (gitignored). Every artifact is recorded in `data/manifest.tsv`.

---

## CLI Reference

### Phase 1 - ingestion (`src.pipeline`)

```
python -m src.pipeline fetch {references|mags|databases|all} [--overwrite] [--limit N]
python -m src.pipeline verify
python -m src.pipeline status
```

* `--overwrite` re-downloads even if the file already exists with a matching hash.
* `--limit N` caps the number of assemblies fetched per BioProject. Use it for smoke
  testing before committing to the full ~450 MB download.
* The pipeline is **idempotent**: artifacts already present (with non-zero size) are
  hashed in place and the manifest row is updated without redownloading.

### Phase 2 - FaRLiP detection (`src.pipeline_phase2`)

```
python -m src.pipeline_phase2 orf-predict   [--overwrite] [--accession ACC]...
python -m src.pipeline_phase2 build-profiles [--overwrite]
python -m src.pipeline_phase2 search         [--overwrite] [--accession ACC]...
python -m src.pipeline_phase2 parse-synteny  [--overwrite]
python -m src.pipeline_phase2 export-vault   [--overwrite]
python -m src.pipeline_phase2 all            [--overwrite]
python -m src.pipeline_phase2 status
```

* Subcommands are independently re-runnable. Stamp files under `data/orfs/` and
  the SHA-256 of every input FASTA make ORF prediction strictly idempotent;
  hmmsearch and the synteny / vault stages overwrite their outputs on each run.
* `data/phase2_manifest.tsv` records every stage's status, output path, and
  wall-time so `all` can resume from a partial run and so downstream auditors can
  reconstruct exactly what produced which file.

### Phase 1.5 - SSU extraction + taxonomic-placement routing (`src.pipeline_phase1_5`)

```
python -m src.pipeline_phase1_5 extract-16S       [--overwrite] [--accession ACC]...
                                                  [--threads N] [--kingdom {bac|arc|euk|mito}]
                                                  [--evalue E] [--reject R]
python -m src.pipeline_phase1_5 bundle-silva
python -m src.pipeline_phase1_5 gen-gtdbtk-script [--threads N]
python -m src.pipeline_phase1_5 all               [--overwrite]
python -m src.pipeline_phase1_5 status
```

* **Requires `barrnap` on PATH.** Install via `conda install -c bioconda
  barrnap` (or any equivalent). The pipeline never auto-installs `barrnap`;
  `extract-16S` exits 2 with a clear message when it cannot find the binary.
* **`extract-16S`** runs `barrnap` once per assembly listed in
  `data/manifest.tsv`, writes `data/ssu/<acc>.barrnap.gff` + `<acc>.16S.fasta`
  + an idempotency stamp, and upserts a row into `data/ssu/routing_log.tsv`.
  Each MAG receives a `decision` tag of `silva` (>=1 16S record) or `gtdbtk`
  (no 16S). Per-MAG stamps key on the input FASTA SHA-256 and the exact
  `barrnap` argv so a parameter change forces a re-run.
* **`bundle-silva`** concatenates every `decision == silva` MAG's
  `<acc>.16S.fasta` into a single `data/ssu/silva_sina_input.fasta` ready
  for SILVA SINA submission (web upload at <https://www.arb-silva.de/aligner/>
  or the SINA CLI against the on-disk SILVA 138.2 NR99 FASTA).
* **`gen-gtdbtk-script`** renders `scripts/run_gtdbtk.sh` (+
  `data/ssu/gtdbtk_batchfile.tsv`) for every `decision == gtdbtk` MAG.
  Gzipped ENA assemblies are decompressed once into
  `data/ssu/gtdbtk_input/<acc>.fna` so `--extension fna` works for the
  whole batch. The script begins with `set -euo pipefail`, requires
  `GTDBTK_DATA_PATH` to be exported, and is `bash -n`-clean.
* **`status`** prints whether the routing log, SILVA bundle, and GTDB-Tk
  script are present plus a count of MAGs routed each way.

Phase 1.5 progress is recorded in `data/phase1_5_manifest.tsv` (per-stage
status + wall time) for resumability mirroring Phase 2 / Phase 3.

### Phase 3 - genome quality assessment (`src.pipeline_phase3`)

```
python -m src.pipeline_phase3 fetch-markers              [--overwrite]
python -m src.pipeline_phase3 assess                     [--overwrite] [--accession ACC]...
python -m src.pipeline_phase3 gen-checkm2-script
python -m src.pipeline_phase3 ingest-checkm2             [--accession ACC]...
python -m src.pipeline_phase3 report                     [--require-farlip] [--include-low]
python -m src.pipeline_phase3 update-vault               [--overwrite]
python -m src.pipeline_phase3 export-extremophile-vault  [--routing-log PATH]
python -m src.pipeline_phase3 all                        [--overwrite]
python -m src.pipeline_phase3 status
```

* **Default (in-process) QC path.** `assess` runs `pyhmmer` against a bacterial
  single-copy ortholog (SCO) HMM database (BUSCO `bacteria_odb10`, 124 markers,
  ~6 MB) fetched once by `fetch-markers`. Each MAG receives a per-marker bit
  score gated by BUSCO's own `scores_cutoff` file (the canonical paralog
  filter) plus an E-value floor. Outputs land in `data/qc_assessor/<acc>.qc.tsv`.
* **Opt-in CheckM2 cross-check.** `gen-checkm2-script` emits
  `scripts/run_checkm2.sh`. You run the script yourself on any machine that has
  CheckM2 + the CheckM2 diamond DB installed (`export CHECKM2DB=...`); the
  resulting `data/checkm2_results/<acc>/quality_report.tsv` files are then
  picked up by `ingest-checkm2`. When both sources exist for the same MAG,
  `report` prefers CheckM2.
* **`report`** joins everything (pyhmmer / CheckM2 QC + Phase 2 FaRLiP
  per-MAG roll-up + Phase 1 manifest) into `data/qc_report.tsv` (every MAG)
  and `data/qc_filtered_mags.tsv` (only `recommended_for_phylogeny == True`,
  i.e. MIMAG HIGH/MEDIUM **and** a FaRLiP candidate). Pass `--require-farlip`
  to drop MAGs with no FaRLiP markers, `--include-low` to keep MIMAG LOW.
* **`update-vault`** re-renders the Phase 2 Obsidian vault, adding
  `completeness`, `contamination`, `n50_contig`, `genome_size_bp`,
  `qc_quality`, `qc_source`, and `recommended_for_phylogeny` to every
  per-MAG YAML frontmatter, plus a "## Quality" section grouped by tier
  in `Phase_2_Master_Index.md`.
* **`export-extremophile-vault`** renders `data/extremophile_vault/`:
  one per-MAG markdown file per `recommended_for_phylogeny == True` row in
  `qc_filtered_mags.tsv` plus the master index
  `High_Quality_Extremophile_Genomes.md`. Each per-MAG note carries the
  exact five-key frontmatter mandated by the assignment - `Accession`,
  `Completeness`, `Contamination`, `Taxonomic_Placement`,
  `FaRLiP_Cluster_Present` - with `FaRLiP_Cluster_Present` left **blank**
  for manual review. `Taxonomic_Placement` is filled with
  `"SILVA SINA: pending"` / `"GTDB-Tk: pending"` placeholders when the
  Phase 1.5 routing log exists at `data/ssu/routing_log.tsv`; otherwise the
  slot is left empty. The exporter is also part of `all`, so the standard
  Phase 3 workflow produces both vaults.
* `data/phase3_manifest.tsv` records every stage's status, output path, and
  wall time, mirroring Phase 2 for resumability.

---

## Configuration

Three config files drive every fetcher. Editing them changes what is fetched without
touching any Python.

| File | Purpose |
|---|---|
| [`config/references.json`](config/references.json) | FaRLiP reference anchors. `primary` is the canonical 21-gene cluster organism; `secondary` is a list (possibly empty) for cross-checks. |
| [`config/extremophiles.json`](config/extremophiles.json) | MAG sources. `bioprojects` = whole-project enumeration; `individual_accessions` = one-offs; `cave_arctic_pending` = the placeholder slot for Dr. Barton's accessions. |
| [`config/databases.yaml`](config/databases.yaml) | SILVA + GTDB versions and download URLs. |
| [`config/farlip_genes.yaml`](config/farlip_genes.yaml) | The 21-gene FaRLiP catalog used by `src/profile_builder.py`. Each entry pins a JSC-1 locus_tag + RefSeq protein_id with fuzzy fallbacks (gene_tag, product_substring) for future anchors. |
| [`config/qc_markers.yaml`](config/qc_markers.yaml) | Single-copy ortholog marker HMM sources for Phase 3 QC. Default = BUSCO `bacteria_odb10` (124 markers, ~6 MB, stable URL + SHA-256). `gtdb_bac120_local` provides an opt-in escape hatch for users who already have GTDB-Tk's bac120 markers extracted locally. |

### Adding a single accession (for example, when Dr. Barton sends a cave genome)

Open `config/extremophiles.json`, append to `cave_arctic_pending`:

```json
{
  "accession": "GCA_XXXXXXXXX.X",
  "source": "NCBI",
  "label": "cave_pool_descriptor",
  "biome": "cave",
  "notes": "Provided by Dr. Barton on YYYY-MM-DD"
}
```

Then re-run `python -m src.pipeline fetch mags`. The pipeline auto-dispatches NCBI
accessions to the Datasets v2 API and ENA accessions to the ENA Browser API; existing
artifacts are skipped.

---

## Repository Layout

```
.
+- README.md                    (this file)
+- DATA_SOURCES.md               (canonical lock-in: who, what, why, citations)
+- DEVLOG.md                     (reverse-chronological dev journal)
+- requirements.txt
+- config/
|  +- references.json
|  +- extremophiles.json
|  +- databases.yaml
+- src/
|  +- _http.py                   (shared HTTP + checksum primitives)
|  +- config_loader.py
|  +- entrez_fetcher.py          (NCBI Datasets v2 + Biopython Entrez fallback)
|  +- ena_fetcher.py             (ENA Portal API + Browser API)
|  +- silva_fetcher.py
|  +- gtdb_fetcher.py
|  +- manifest.py                (data/manifest.tsv writer)
|  +- pipeline.py                (Phase 1 CLI orchestrator)
|  +- orf_predictor.py           (Phase 2: pyrodigal meta-mode ORF prediction)
|  +- profile_builder.py         (Phase 2: parses JSC-1 GFF, builds 21 HMMs via pyhmmer)
|  +- hmm_searcher.py            (Phase 2: pyhmmer hmmsearch with module-level E-value thresholds)
|  +- synteny_parser.py          (Phase 2: 3-in-10kb candidate rule, 30kb classification window, EDGE_TRUNCATED tagging)
|  +- obsidian_exporter.py       (Phase 2: TSV reports -> Obsidian vault; Phase 3 QC-aware)
|  +- pipeline_phase2.py         (Phase 2 CLI orchestrator)
|  +- qc_markers_fetcher.py      (Phase 3: BUSCO/GTDB SCO marker HMM downloader)
|  +- qc_assessor.py             (Phase 3: pyhmmer in-process completeness/contamination)
|  +- checkm2_shellscript_generator.py  (Phase 3: emits scripts/run_checkm2.sh)
|  +- checkm2_ingestor.py        (Phase 3: parse CheckM2 quality_report.tsv -> common schema)
|  +- qc_reporter.py             (Phase 3: joins QC + FaRLiP + manifest into qc_report.tsv)
|  +- extremophile_vault_exporter.py    (Phase 3: data/extremophile_vault/ with 5-key frontmatter)
|  +- pipeline_phase3.py         (Phase 3 CLI orchestrator)
|  +- ssu_extractor.py           (Phase 1.5: barrnap subprocess wrapper + 16S extraction)
|  +- ssu_bundler.py             (Phase 1.5: SILVA SINA multi-FASTA + routing_log.tsv writer)
|  +- gtdbtk_shellscript_generator.py   (Phase 1.5: emits scripts/run_gtdbtk.sh)
|  +- pipeline_phase1_5.py       (Phase 1.5 CLI orchestrator)
+- tests/                        (pytest, no network)
+- scripts/fetch_all.sh
+- data/                         (gitignored output)
   +- references/<accession>/    (extracted NCBI Datasets bundle)
   +- mags/<bioproject>/         (per-BioProject ENA FASTAs)
   +- mags/individual/           (one-off accessions, NCBI or ENA)
   +- silva/<version>/
   +- gtdb/<version>/
   +- manifest.tsv               (Phase 1 source of truth)
   +- orfs/                      (Phase 2: <accession>.faa + <accession>.coords.tsv + <accession>.stamp)
   +- hmm_profiles/              (Phase 2: farlip_21.hmm + farlip_21.catalog.tsv)
   +- hmm_results/               (Phase 2: <accession>.domtblout per assembly)
   +- farlip_hits_report.tsv     (Phase 2: per-candidate-contig report)
   +- farlip_hits_per_mag.tsv    (Phase 2: per-MAG roll-up)
   +- obsidian_export/           (Phase 2: <accession>.md + Phase_2_Master_Index.md)
   +- phase2_manifest.tsv        (Phase 2: per-stage status, output paths, wall times)
   +- qc_markers/                (Phase 3: <key>.hmm + <key>.scores_cutoff.tsv)
   +- qc_assessor/               (Phase 3: <accession>.qc.tsv from pyhmmer SCO scan)
   +- checkm2_results/           (Phase 3: user-populated; per-accession CheckM2 outputs)
   +- checkm2_ingest/            (Phase 3: normalized <accession>.qc.tsv from CheckM2)
   +- qc_report.tsv              (Phase 3: joined QC + FaRLiP report)
   +- qc_filtered_mags.tsv       (Phase 3: recommended_for_phylogeny == True only)
   +- phase3_manifest.tsv        (Phase 3: per-stage status, output paths, wall times)
   +- extremophile_vault/        (Phase 3: <accession>.md + High_Quality_Extremophile_Genomes.md)
   +- ssu/                       (Phase 1.5: <accession>.barrnap.gff + <accession>.16S.fasta
   |                              + routing_log.tsv + silva_sina_input.fasta + gtdbtk_batchfile.tsv)
   +- ssu/gtdbtk_input/          (Phase 1.5: plain-FASTA staging for GTDB-Tk classify_wf)
   +- gtdbtk_results/            (Phase 1.5: user-populated; GTDB-Tk classify_wf outputs)
   +- phase1_5_manifest.tsv      (Phase 1.5: per-stage status, output paths, wall times)
+- scripts/
   +- fetch_all.sh
   +- run_checkm2.sh             (Phase 3: optional CheckM2 cross-check, generated)
   +- run_gtdbtk.sh              (Phase 1.5: GTDB-Tk classify_wf for non-16S MAGs, generated)
```

### `data/manifest.tsv` columns

`accession`, `source`, `url`, `local_path`, `sha256`, `bytes_downloaded`,
`fetched_at`, `status` (`ok` / `skipped` / `failed` / `changed` / `missing`),
`note`, plus source-specific fields (`label`, `scientific_name`, `sample_accession`,
`assembly_level`, `format`, `md5`, `extracted_dir`).

Iterations 2+ read this manifest to know what is on disk and where.

---

## Environment Variables

| Variable | Used by | Effect |
|---|---|---|
| `NCBI_EMAIL` | `entrez_fetcher` | Required by NCBI usage policy. Defaults to a placeholder with a stderr warning. |
| `NCBI_API_KEY` | `entrez_fetcher` | Optional. Raises rate limit from 3 to 10 requests / second. |

ENA, SILVA, and GTDB endpoints are unauthenticated.

---

## Tests

```bash
.venv/bin/python -m pytest tests/ -v
```

140 tests, all offline-only (20 Phase 1 + 37 Phase 2 + 47 Phase 3 + 36 Phase
1.5 / extremophile vault). Phase 1 mocks NCBI Datasets and the ENA Portal;
Phase 2 uses tiny in-memory FASTA / GFF fixtures plus hand-rolled HMMER
domtblout fixtures; Phase 3 mocks the BUSCO tarball, the per-MAG `pyhmmer`
scan, and CheckM2 `quality_report.tsv` parsing. Phase 1.5 mocks `subprocess.run`
so `barrnap` is never actually invoked - a stub GFF3 fixture drives the
extraction / routing / bundling code paths, and the generated GTDB-Tk shell
script is syntax-checked with `bash -n`. The suite works with no network and
no third-party binaries on PATH.

---

## Roadmap (more to come)

| Iteration | Status | Scope |
|---|---|---|
| **1** | shipped | Ingestion + manifest. Reference, MAGs, SILVA, GTDB on disk, SHA-256-verified. |
| **1.5** | shipped | `barrnap` 16S extraction per MAG + per-decision routing log (`silva` / `gtdbtk`) + SILVA SINA multi-FASTA bundle + generated `scripts/run_gtdbtk.sh` for genome-wide marker classification on MAGs without recoverable 16S. SILVA SINA and GTDB-Tk themselves remain user-run (network upload / 141 GB reference data). |
| **1.6** | planned | `ingest-silva` / `ingest-gtdbtk` that parse the SINA result file and the GTDB-Tk `bac120.summary.tsv` and overwrite `Taxonomic_Placement` in every per-MAG extremophile-vault note from `pending` to the resolved species / lineage. Symmetric to Phase 3's `ingest-checkm2`. |
| **2** | shipped | ORF prediction (pyrodigal meta) + 21 single-sequence HMMs from JSC-1 + pyhmmer hmmsearch + 3-in-10kb candidate trigger / 30kb classification window with EDGE_TRUNCATED tagging + Obsidian export. **Validation:** JSC-1 -> FULL_CLUSTER (21/21 markers in 30 kb on contig NZ_KL662192.1). |
| **2.1** | planned | Integrate Chlorogloeopsis fritschii PCC 9212 as secondary anchor. Switch `profile_builder` from single-sequence to MSA-based via mafft subprocess. |
| **3** | shipped | Per-MAG completeness / contamination / N50 / genome-size via in-process pyhmmer against BUSCO `bacteria_odb10` (124 SCOs, with per-HMM `scores_cutoff` paralog gate) + MIMAG HIGH/MEDIUM/LOW tier assignment + optional CheckM2 shell-script cross-check + QC-aware Obsidian re-render. **Validation:** JSC-1 -> 94.4% completeness / 8.1% contamination (MIMAG MEDIUM, within published BUSCO-vs-CheckM2 envelope for cyanobacterial RefSeq assemblies); 39-MAG corpus -> 1 HIGH / 36 MEDIUM / 2 LOW; 36 recommended for phylogeny. |
| **3.5** | planned | Per-domain marker sets (archaeal `ar122` etc.) for non-bacterial MAGs once they land in the corpus. |
| **4** | planned | Phylogenetic placement of `qc_filtered_mags.tsv` via SILVA 16S (Phase 1 already has SILVA 138.2 NR99) and/or GTDB-Tk placement (Phase 1 already has GTDB R226 metadata). Tree builder + biome cross-validation. |

The Phase 1 + Phase 2 + Phase 3 design choices specifically support the later stages:

* NCBI Datasets bundles include protein FASTA -> already consumed in Phase 2.
* ENA assemblies arrive as gzipped FASTA -> `orf_predictor` reads gzip natively.
* `data/manifest.tsv` is the unified Phase 1 input; `data/phase2_manifest.tsv`
  and `data/phase3_manifest.tsv` record every later stage's status and are the
  resumability anchors for their respective `all` subcommands.
* Phase 3 reuses the Phase 2 `.faa` outputs (no re-prediction) and writes both
  CheckM2 and pyhmmer rows into the *same* schema, so the reporter does not
  need to know which engine ran.

---

## Action Items

* **Bring to next meeting with Dr. Barton:**
  *"Please provide the exact cave + additional Arctic genome accessions (NCBI BioProject,
  ENA BioProject, or individual assembly accessions) to populate
  `config/extremophiles.json -> cave_arctic_pending`."*

---

## Citations

The primary references are tabulated in [`DATA_SOURCES.md`](DATA_SOURCES.md).
Key papers:

* Gan, F. *et al.* (2014). Extensive remodeling of a cyanobacterial photosynthetic
  apparatus in far-red light. *Science* 345:1312-1317.
* Pessi, I. S. *et al.* (2023). Novel diversity of polar Cyanobacteria revealed by
  genome-resolved metagenomics. *Microbial Genomics* 9:mgen001001.
* Quast, C. *et al.* (2013). The SILVA ribosomal RNA gene database project.
  *Nucleic Acids Research* 41:D590-D596.
* Parks, D. H. *et al.* (2022). GTDB: an ongoing census of bacterial and archaeal
  diversity. *Nucleic Acids Research* 50:D785-D794.
* Bowers, R. M. *et al.* (2017). Minimum information about a single amplified
  genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and
  archaea. *Nature Biotechnology* 35:725-731.
* Manni, M. *et al.* (2021). BUSCO update: novel and streamlined workflows along
  with broader and deeper phylogenetic coverage for scoring of eukaryotic,
  prokaryotic, and viral genomes. *Molecular Biology and Evolution* 38:4647-4654.
* Chklovski, A. *et al.* (2023). CheckM2: a rapid, scalable and accurate tool for
  assessing microbial genome quality using machine learning. *Nature Methods*
  20:1203-1212.
* Seemann, T. (2013). barrnap: BAsic Rapid Ribosomal RNA Predictor.
  GitHub: <https://github.com/tseemann/barrnap>
* Pruesse, E., Peplies, J., & Glöckner, F. O. (2012). SINA: Accurate
  high-throughput multiple sequence alignment of ribosomal RNA genes.
  *Bioinformatics* 28:1823-1829.
* Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P., & Parks, D. H. (2022).
  GTDB-Tk v2: memory-friendly classification with the Genome Taxonomy
  Database. *Bioinformatics* 38:5315-5316.
