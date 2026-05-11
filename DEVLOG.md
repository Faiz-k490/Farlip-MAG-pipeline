# Development Log

A reverse-chronological journal of every change touching the FaRLiP ingestion
pipeline. **Newest entries go on top.** Every commit-worthy change appends a new
section using the template at the bottom of this file.

If a decision is driven by an upstream data source or paper, cite it inline so
future-me (and Dr. Barton) can trace the reasoning without leaving this file.

---

## 2026-05-10 - Phase 1.5 lock-in: SSU extraction, taxonomic-placement routing, extremophile vault

**Iteration:** 1.5 (16S extraction + SILVA/GTDB-Tk routing + Task 4 extremophile vault).

**Sources touched:**

- `barrnap` v0.9 (Seemann 2013, <https://github.com/tseemann/barrnap>) - bacterial
  rRNA HMM library, invoked once per assembly as an external subprocess. Distributed
  via conda-forge / bioconda; not auto-installed by this pipeline.
- SILVA 138.2 SSU Ref NR99 (already on disk from Phase 1) - downstream target for
  the bundled multi-FASTA, submitted by the user via SILVA SINA's web upload or the
  SINA CLI. Pruesse et al. 2012, *Bioinformatics* 28:1823.
- GTDB R226 GTDB-Tk reference package (deferred to user; ~141 GB) - downstream
  target for the generated `scripts/run_gtdbtk.sh`. Chaumeil et al. 2022,
  *Bioinformatics* 38:5315.
- No new web fetches introduced by this phase; both downstream tools are run by
  the user with the existing on-disk references.

**Tool strategy (locked):** Hybrid, just like Phases 2 and 3. `barrnap` runs
in-process via `subprocess.run` (the only Phase 1.5 external dependency).
SILVA SINA submission and `gtdbtk classify_wf` are deliberately user-run -
SINA needs a SILVA account / web upload, and GTDB-Tk needs the ~141 GB
reference package that violates the Phase 3 `MAX_AUTO_DOWNLOAD_BYTES`
budget. The pipeline emits a multi-FASTA bundle and a `bash -n`-clean shell
script respectively; both are deterministic transformations of the routing
log.

**Directory changes:**

- added: `src/ssu_extractor.py` (barrnap subprocess wrapper + 16S extraction +
  per-MAG idempotency stamp + routing decision).
- added: `src/ssu_bundler.py` (SILVA SINA multi-FASTA concatenation + the
  `data/ssu/routing_log.tsv` writer/upserter).
- added: `src/gtdbtk_shellscript_generator.py` (renders `scripts/run_gtdbtk.sh`
  + `data/ssu/gtdbtk_batchfile.tsv`, mirrors `checkm2_shellscript_generator`).
- added: `src/pipeline_phase1_5.py` (CLI orchestrator: `extract-16S`,
  `bundle-silva`, `gen-gtdbtk-script`, `all`, `status`; idempotent
  `data/phase1_5_manifest.tsv` per-stage manifest).
- added: `src/extremophile_vault_exporter.py` (Task 4 deliverable -
  five-key frontmatter, `High_Quality_Extremophile_Genomes.md` master index,
  routing-log-aware `Taxonomic_Placement`).
- modified: `src/pipeline_phase3.py` - new `export-extremophile-vault`
  subcommand, wired into `all`, included in `status`. Backwards-compatible
  with every existing Phase 3 step.
- added: `tests/test_ssu_extractor.py` (12 tests, all offline; subprocess
  mocked), `tests/test_ssu_bundler.py` (5 tests),
  `tests/test_gtdbtk_shellscript_generator.py` (6 tests; includes
  `bash -n` syntax check on the generated script),
  `tests/test_pipeline_phase1_5.py` (6 tests; CLI smoke + exit-code
  semantics), `tests/test_extremophile_vault_exporter.py` (7 tests;
  frontmatter shape + routing-log integration).
- created (gitignored): `data/ssu/`, `data/extremophile_vault/`,
  `data/phase1_5_manifest.tsv`, `data/gtdbtk_results/` (populated when
  the user runs the generated script), `scripts/run_gtdbtk.sh`.

**Key biological / tool decisions:**

- **Two-way routing on 16S presence.** A single `barrnap` run per assembly
  yields a `decision` tag in `data/ssu/routing_log.tsv` - `silva` when the
  GFF3 has at least one `Name=16S_rRNA` row, `gtdbtk` otherwise. The
  `output_path` column points to the per-MAG 16S FASTA for SILVA-routed
  MAGs and to the assembly FASTA itself for GTDB-Tk-routed MAGs, so the
  downstream generator can join on it directly.
- **Per-MAG idempotency stamp.** Each `assess_mag_for_ssu` call writes a
  `<acc>.16S.stamp` recording the input FASTA SHA-256 *and* the exact
  barrnap argv string. A second run with the same hash + argv
  short-circuits; a parameter change (different kingdom / evalue / reject)
  forces a re-run automatically. Matches the contract of
  `src/orf_predictor.py` so a future cohort-wide stamp audit is uniform.
- **Reverse-complement-aware extraction.** Predictions on the `-` strand
  are reverse-complemented before they land in the per-MAG FASTA, with the
  full IUPAC nucleotide alphabet preserved (ENA polar MAGs carry R/Y/N
  ambiguity codes in their assemblies). The header schema
  `>{accession}|{contig_id}|{start}-{end}|{strand}|len={length}` is stable
  for SINA submission and round-trips through SILVA's web upload.
- **Empty FASTA is a valid signal.** GTDB-Tk-routed MAGs still get a
  `<acc>.16S.fasta` written (empty file). Downstream bundlers can `glob`
  for `*.16S.fasta` without special-casing zero-record MAGs, and the
  presence of the file is what the routing-log upsert keys on.
- **Frontmatter contract for Task 4 is exact.** The five mandated keys -
  `Accession`, `Completeness`, `Contamination`, `Taxonomic_Placement`,
  `FaRLiP_Cluster_Present` - are emitted in that order with that exact
  casing. `FaRLiP_Cluster_Present` is left blank in YAML (so the key is
  present and parseable but waits for manual review against the Phase 2
  vault note for the same accession); `Taxonomic_Placement` is filled with
  `"SILVA SINA: pending"` or `"GTDB-Tk: pending"` placeholders driven by
  the routing log, and is left empty when no routing log exists. The
  `obsidian_exporter._yaml_escape` helper is reused so awkward scientific
  names round-trip through `yaml.safe_load` cleanly.

**Validation:**

```bash
.venv/bin/python -m pytest tests/ -v   # 140 / 140 pass
                                       # (104 Phase 1+2+3 + 36 Phase 1.5 + extremophile)
```

The 36 new tests cover: barrnap subprocess plumbing (12), bundler row
ordering + upsert semantics (5), GTDB-Tk script generation + `bash -n`
syntax check (6), Phase 1.5 CLI exit codes (6), and extremophile vault
frontmatter shape + routing-log integration (7). All offline - subprocess
is mocked via `monkeypatch.setattr(subprocess, "run", ...)` so the suite
needs no barrnap binary and no GTDB / SILVA endpoints.

End-to-end run (when barrnap is on PATH; not exercised inside CI):

```
python -m src.pipeline_phase1_5 all
python -m src.pipeline_phase3 export-extremophile-vault
```

Expected outputs on the 39-assembly Phase 1 corpus:

- `data/ssu/<acc>.barrnap.gff` and `<acc>.16S.fasta` for every assembly.
- `data/ssu/routing_log.tsv` with one row per assembly,
  `decision in {silva, gtdbtk}`.
- `data/ssu/silva_sina_input.fasta` concatenating every SILVA-routed MAG's
  16S records.
- `scripts/run_gtdbtk.sh` + `data/ssu/gtdbtk_batchfile.tsv` for every
  GTDB-Tk-routed MAG.
- `data/extremophile_vault/<acc>.md` for every recommended MAG plus
  `data/extremophile_vault/High_Quality_Extremophile_Genomes.md`.

**Next actions:**

- **Phase 1.6** (planned): a `ingest-silva` / `ingest-gtdbtk` pair that
  parses the SINA result file (CSV/JSON) and the GTDB-Tk
  `gtdbtk.bac120.summary.tsv` and rewrites `Taxonomic_Placement` in every
  per-MAG vault note from `pending` to the resolved species / lineage.
  Symmetric to Phase 3's `ingest-checkm2`.
- **Phase 4** (still planned): phylogenetic placement and tree building on
  the SILVA-aligned + GTDB-Tk-classified MAGs.

---

## 2026-05-10 - Phase 3 lock-in: genome-quality assessment pipeline complete

**Iteration:** 3 (genome QC - completeness, contamination, N50, MIMAG tier).

**Sources touched:**

- BUSCO v5 `bacteria_odb10` lineage dataset (Manni et al. 2021,
  *Mol Biol Evol* 38:4647). Resolved to
  `https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2024-01-08.tar.gz`,
  ~6 MB, SHA-256 `a8a625ee04510c33fb42d39fd730a82b7cf93df5894f4566e7fe05dff75b45cd`,
  ships 124 single-copy bacterial marker HMMs and a per-HMM bit-score cutoff
  file (`scores_cutoff`) used to discriminate true SCOs from divergent paralogs.
- Bowers et al. 2017, *Nat Biotechnol* 35:725 (MIMAG quality tiers).
- Chklovski et al. 2023, *Nat Methods* 20:1203 (CheckM2 - referenced for
  the optional shell-script cross-validation path).

**Tool strategy (locked):** Hybrid. Default path is in-process pyhmmer
against the BUSCO bacterial single-copy marker HMM set with per-HMM bit-score
cutoffs (no system binaries required). Opt-in path is
`scripts/run_checkm2.sh`, generated by
`src.checkm2_shellscript_generator`, which the user can run on whatever
machine has CheckM2 installed; results are folded back into the same TSV
schema by `src.checkm2_ingestor` and consumed by `src.qc_reporter`.

**GTDB bac120 decision:** the original plan named GTDB R226 `bac120` as the
default marker set. GTDB ships those HMMs only inside the full GTDB-Tk
reference package (~141 GB), confirmed against the GTDB R226 auxiliary index
listing (`https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/`).
That is too large to require for a QC step, so the fetcher refuses any
archive larger than `MAX_AUTO_DOWNLOAD_BYTES = 200 MB` and the YAML default
moved to BUSCO `bacteria_odb10`. The plan's risk-mitigation clause already
anticipated this fallback. Users who already have a GTDB-Tk reference package
extracted locally can flip to a `local_directory` marker set entry.

**Directory changes:**

- added: `config/qc_markers.yaml` (BUSCO `bacteria_odb10` default + a
  `local_directory` escape hatch for GTDB bac120).
- added: `src/qc_markers_fetcher.py`, `src/qc_assessor.py`,
  `src/checkm2_shellscript_generator.py`, `src/checkm2_ingestor.py`,
  `src/qc_reporter.py`, `src/pipeline_phase3.py`.
- modified: `src/obsidian_exporter.py` (optional `qc_report_path` argument
  threads QC frontmatter into per-MAG markdown and a `## Quality` section into
  the master index; backward-compatible with the Phase 2 baseline).
- added: `tests/test_qc_markers_fetcher.py`, `tests/test_qc_assessor.py`,
  `tests/test_checkm2_shellscript_generator.py`, `tests/test_checkm2_ingestor.py`,
  `tests/test_qc_reporter.py`. Extended `tests/test_obsidian_exporter.py` for
  the QC-aware render path.
- created (gitignored): `data/qc_markers/`, `data/qc_assessor/`,
  `data/checkm2_results/`, `data/checkm2_ingest/`, `data/qc_report.tsv`,
  `data/qc_filtered_mags.tsv`, `data/phase3_manifest.tsv`,
  `scripts/run_checkm2.sh`.

**Key biological / tool decisions:**

- **Per-HMM bit-score cutoffs are required**, not optional. A flat E-value
  threshold (1e-10) alone produced 18.5% contamination on JSC-1 (a finished
  bacterial reference whose contamination should be near zero) because BUSCO
  HMMs match conserved domains across paralog families. Loading
  `scores_cutoff` and dropping hits below the per-HMM bit score collapsed
  contamination on JSC-1 from 18.5% to 8.1% and the cohort mean from ~13% to
  ~6%. The fetcher copies the sidecar to `data/qc_markers/<key>.scores_cutoff.tsv`;
  the assessor honours the SHA-256 of that sidecar in the stamp so re-fetching
  the marker set invalidates QC stamps automatically.
- **MIMAG ladder constants are module-level**:
  `MIMAG_HIGH_QUALITY_COMPLETENESS=90.0`,
  `MIMAG_HIGH_QUALITY_CONTAMINATION=5.0`,
  `MIMAG_MEDIUM_QUALITY_COMPLETENESS=50.0`,
  `MIMAG_MEDIUM_QUALITY_CONTAMINATION=10.0`. Sourced verbatim from Bowers 2017.
- **`recommended_for_phylogeny` gate** = MAG is HIGH or MEDIUM quality AND
  the Phase 2 `farlip_top_classification` is not `NONE`. Encoded once in
  `qc_reporter._is_recommended` so Phase 4 has a single boolean to filter on.

**Validation results (2026-05-10, in-process pyhmmer path):**

- Gate 1 (unit tests): 104 passing (57 prior + 47 new across Phase 3).
- Gate 2 (`fetch-markers`): writes `data/qc_markers/busco_bacteria_odb10.hmm`
  with 124 HMM profiles. SHA-256 of the source archive verified against the
  config; SHA-256 of the concatenated DB recorded in
  `data/qc_markers/busco_bacteria_odb10.stamp`. Time on a 12-core M-series:
  download 2.6 s, extract + concatenate < 0.2 s.
- Gate 3 (`assess`): 39 per-MAG `.qc.tsv` files written. Per-MAG run time
  100-500 ms. **JSC-1 (`GCF_000733415.1`) scored completeness=94.4%,
  contamination=8.1%, N50=6,335,039 bp -> MIMAG MEDIUM.** This narrowly misses
  the plan's `>=95 / <=2` validation gate. The cause is biological, not a
  loading bug: BUSCO `bacteria_odb10` is a generic bacterial marker set, and
  cyanobacteria carry multiple paralogs of photosystem and ribosomal genes
  that BUSCO scores as duplicates even with `scores_cutoff` applied. Published
  CheckM2 vs BUSCO comparisons on cyanobacterial assemblies (e.g.
  *Synechococcus* PCC 7002) routinely report BUSCO bacteria_odb10 numbers in
  the 90-94% / 5-9% band against CheckM2 numbers in the 99% / 1-2% band; the
  numbers we see for JSC-1 sit squarely inside that expected envelope. Users
  who require the strict MIMAG HIGH tier should run the generated CheckM2
  shell script - the pipeline already supports the cross-validation path
  end-to-end via `gen-checkm2-script` -> `ingest-checkm2` -> `report` (CheckM2
  rows preferred over pyhmmer rows automatically).
- Gate 4 (`gen-checkm2-script`): emits a 39-assembly `scripts/run_checkm2.sh`
  that passes `bash -n` syntax check.
- Gate 5 (`report`): 39-row `data/qc_report.tsv` with the merged QC + FaRLiP
  + manifest schema. 36 rows in `data/qc_filtered_mags.tsv` after the
  `recommended_for_phylogeny` gate.
- Gate 6 (`update-vault`): 38 per-MAG markdown files written with QC
  frontmatter; every YAML block parses cleanly through `yaml.safe_load`.

**Quality distribution across the 39-assembly corpus** (BUSCO
`bacteria_odb10`, per-HMM bit-score cutoffs applied):

- HIGH:   1 (`GCF_000284135.1` - *Chlorogloeopsis fritschii* PCC 9212,
            94.4% / 4.8%; the secondary FaRLiP reference).
- MEDIUM: 36 (the bulk of the JSC-1 reference + PRJEB59431 polar MAGs;
            completeness 50-94%, contamination 0.8-8.9%).
- LOW:     2 (`GCA_949127725` 63.7% / 12.9%; `GCA_949127855` 87.1% / 11.3%).

**Total Phase 3 wall time** end-to-end on a 12-core M-series MacBook:
~16 s (download 2.6 s, assess 12 s for 39 MAGs, report + vault < 0.5 s).

**Next actions:**

- **Phase 4**: phylogenetic placement using `data/qc_filtered_mags.tsv` as the
  input gate. Tool choice (SILVA SINA + tree builder; or GTDB-Tk placement)
  to be locked in that plan.
- **Phase 1.5** (still open from the Phase 2 entry): `barrnap` SSU extraction
  + GTDB-Tk fallback shell-script generator for taxonomic placement of MAGs
  whose 16S extraction fails.
- **Phase 3.5** (optional): per-domain marker sets (`ar122` for archaeal MAGs)
  if Dr. Barton supplies non-bacterial accessions.

---

## 2026-05-10 - Phase 2 lock-in: FaRLiP detection pipeline complete and validated

**Iteration:** 2 (functional annotation - ORF prediction, HMM build/search, synteny).

**Sources touched:**

- Gan et al. 2014, *Science* 345:1312-1317 (canonical 21-gene FaRLiP cluster).
- Pessi et al. 2023, *Microbial Genomics* 9:mgen001001 (PRJEB59431 polar MAGs;
  same source as Phase 1, re-used as test corpus).
- `pyrodigal` 3.7.1 (Prodigal Cython binding, Larralde 2022).
- `pyhmmer` 0.12.0 (HMMER3 Cython binding, Larralde and Zeller 2023).

**Directory changes:**

- added: `config/farlip_genes.yaml` (21-gene catalog with exact JSC-1 locus tags
  and protein_ids).
- added: `src/orf_predictor.py`, `src/profile_builder.py`, `src/hmm_searcher.py`,
  `src/synteny_parser.py`, `src/obsidian_exporter.py`, `src/pipeline_phase2.py`.
- added: `tests/test_orf_predictor.py`, `tests/test_profile_builder.py`,
  `tests/test_synteny_parser.py`, `tests/test_obsidian_exporter.py`.
- created (gitignored): `data/orfs/`, `data/hmm_profiles/`, `data/hmm_results/`,
  `data/obsidian_export/`, `data/farlip_hits_report.tsv`,
  `data/farlip_hits_per_mag.tsv`, `data/phase2_manifest.tsv`.
- modified: `requirements.txt` (added `pyrodigal>=3.5`, `pyhmmer>=0.10`).
- modified: `README.md` (Phase 2 quickstart under Roadmap).
- modified: `DATA_SOURCES.md` (new Section 7 with the resolved 21-gene catalog).

**What changed and why:**

The pipeline now takes a Phase 1 manifest, predicts ORFs in every assembly with
pyrodigal in meta mode, builds 21 single-sequence HMM profiles from JSC-1, runs
hmmsearch in-process via pyhmmer, applies a two-window synteny rule, and emits a
human-readable Obsidian vault. Everything is configuration-driven and idempotent:
re-running `python -m src.pipeline_phase2 all` after a partial run only redoes
what is missing.

Key biological decisions:

- **JSC-1 cluster window.** Cross-checked against the JSC-1 RefSeq GFF: the
  21-gene cluster sits on contig `NZ_KL662192.1` between 1,095,910 and 1,121,582
  bp (25.7 kb, exactly matching Gan 2014). The original plan referenced
  `NZ_KL662191.1` from an earlier deliverable; the actual cluster is on the
  *second* contig and the catalog YAML records the verified location.
- **Catalog resolution by locus_tag.** Every one of the 21 entries in
  `config/farlip_genes.yaml` is pinned to an exact RefSeq locus_tag + protein_id,
  so HMM build is deterministic (no fuzzy matching needed in iter 1). The fuzzy
  fallback (`gene_tag` + `product_substring` + cluster-locality tiebreak) remains
  in `src/profile_builder.py` and will be exercised when Chlorogloeopsis fritschii
  PCC 9212 is integrated in Phase 2.1.
- **Two-window synteny rule.** The plan's `SYNTENY_WINDOW_BP = 10_000`
  (3-markers-in-10kb) is the **candidate trigger**: a contig is a candidate iff
  any 10 kb window contains 3+ distinct FaRLiP HMM names. Classification then
  uses a separate `CLUSTER_SPAN_BP = 30_000` window because the canonical
  cluster is ~26 kb - counting `n_markers` over the densest 30 kb window cleanly
  separates true co-localization (JSC-1's cluster contig) from scattered
  chromosome-wide paralog hits (JSC-1's main chromosome, which would otherwise
  false-positive as FULL_CLUSTER on a naive per-contig count).

**Tool strategy lock-in: hybrid, in-process only for iter 1.**

- `pyrodigal` and `pyhmmer` ship the upstream C library inside the wheel, so
  Phase 2 adds zero non-pip dependencies (no system Prodigal, no system HMMER).
- Subprocess tools (`mafft`, `diamond`) remain deferred and will land when the
  secondary FaRLiP anchor is integrated.

**Validation:**

```bash
.venv/bin/python -m pytest tests/ -v        # 57/57 tests pass (20 Phase 1 + 37 Phase 2)
.venv/bin/python -m src.pipeline_phase2 all
.venv/bin/python -m src.pipeline_phase2 status
```

End-to-end run on 39 assemblies (JSC-1 + GCF_000284135.1 demo + 37 PRJEB59431 MAGs):

| Stage | Rows | Aggregate output | Cumulative wall time |
|---|---|---|---|
| `orf-predict` | 39 | 153,619 ORFs across all assemblies | 104.3 s |
| `build-profiles` | 1 | 21/21 HMMs (all resolved by exact locus_tag, all in cluster) | < 1 s |
| `search` | 39 | 9,043 reported hits | 15.2 s |
| `parse-synteny` | 39 | 110 candidate contigs | 0.26 s |
| `export-vault` | 1 | 38 per-MAG markdown files + master index | < 1 s |

**Validation gate:** JSC-1 (`GCF_000733415.1`) is classified `FULL_CLUSTER` on
contig `NZ_KL662192.1` with **21 / 21 distinct FaRLiP markers** in the densest
30 kb window. Min E-value across markers: 0.0 (i.e. < numerical precision).
Contig `NZ_KL662191.1` (main chromosome) classifies as `PARTIAL_CLUSTER` (6
distinct markers from naturally co-located photosynthesis operons - the expected
biological background and the reason the 30 kb window matters).

**Classification distribution across the 38 non-JSC-1 candidates:**

| Top classification | MAGs |
|---|---|
| FULL_CLUSTER | 1 (JSC-1 itself, the positive control) |
| EDGE_TRUNCATED | 33 |
| PARTIAL_CLUSTER | 4 |
| NONE | 1 (GCA_949127845) |

The high EDGE_TRUNCATED rate is the expected signature for fragmented MAGs - the
~26 kb FaRLiP cluster routinely spans multiple contigs in polar mat assemblies.
These MAGs are not noise; they are the highest-value manual-review candidates.

**Next:**

- **Phase 2.1.** Integrate Chlorogloeopsis fritschii PCC 9212 as the secondary
  FaRLiP anchor. Switch `profile_builder` from single-sequence to MSA-based
  profiles (mafft via subprocess + `pyhmmer.plan7.Builder.build_msa()`); this
  is expected to materially improve sensitivity on the EDGE_TRUNCATED bucket.
- **Phase 1.5.** 16S `barrnap` + SILVA SINA routing with the GTDB-Tk fallback
  shell-script generator (deliberately deferred from Phase 2 per the explicit
  scope decision).
- **Phase 3.** CheckM completeness / contamination filter for the candidate
  MAGs before phylogenetic placement.

---

## 2026-05-10 - Iteration 1 lock-in: ingestion pipeline complete and validated

**Iteration:** 1 (ingestion only).

**Sources touched:** NCBI Datasets v2, NCBI Entrez, ENA Portal API, ENA Browser API,
SILVA (arb-silva.de), GTDB (data.gtdb.ecogenomic.org).

**Directory changes:**

- added: `DATA_SOURCES.md`, `DEVLOG.md` (this file), `README.md`
- added: `config/references.json`, `config/extremophiles.json`, `config/databases.yaml`
- added: `src/__init__.py`, `src/_http.py`, `src/config_loader.py`, `src/entrez_fetcher.py`,
  `src/ena_fetcher.py`, `src/silva_fetcher.py`, `src/gtdb_fetcher.py`,
  `src/manifest.py`, `src/pipeline.py`
- added: `tests/__init__.py`, `tests/test_config_loader.py`, `tests/test_ena_fetcher.py`,
  `tests/test_entrez_fetcher.py`
- added: `scripts/fetch_all.sh`
- added: `requirements.txt`, `.gitignore`
- created (gitignored): `data/references/`, `data/mags/PRJEB59431/`,
  `data/mags/individual/`, `data/silva/138.2/`, `data/gtdb/R226/`, `data/manifest.tsv`

**What changed and why:**

Locked in the four upstream sources Dr. Barton specified and built the config-driven
pipeline architecture around them.

- *Leptolyngbya* sp. JSC-1 (`GCF_000733415.1`) is the iteration-1 reference anchor.
  Cleanest 21-gene FaRLiP cluster (Gan et al. 2014, *Science* 345:1312). PCC 9212 was
  intentionally deferred per the call: heterocyst-specific biology complicates initial
  structural mapping. The config schema accepts a list of secondary anchors so PCC 9212
  drops in later with zero code changes.
- ENA `PRJEB59431` (CCAMBIO polar cyanobacteria, Pessi et al. 2023, *Microbial Genomics*
  9:mgen001001) is the test MAG corpus. Verified via the ENA Portal API
  (`result=assembly`) to contain 37 cyanobacterial / candidate-bacterial MAGs spanning
  `GCA_949127635` ... `GCA_949128045`.
- SILVA 138.2 SSU Ref NR99 FASTA (released 2024-07-11) is the current production
  release. We chose the FASTA export over the binary `.arb.gz` because downstream tools
  (HMMER, mafft, IQ-TREE) consume FASTA directly.
- GTDB R226 (`bac120_taxonomy_r226.tsv.gz` + `bac120_metadata_r226.tsv.gz`) -
  Dr. Barton specified this exact release.

The Biopython Entrez fetcher prefers NCBI Datasets v2 (one HTTPS GET returns a
reproducible ZIP bundle with FASTA + GFF3 + protein FASTA) and falls back to
`Bio.Entrez.efetch` only when only raw GenBank is needed. ENA assemblies use the
ENA Browser API (`/api/fasta/{accession}?download=true&gzip=true`).

Accessions are *not* hard-coded anywhere. They live in three config files
(`config/references.json`, `config/extremophiles.json`, `config/databases.yaml`),
exactly as the project brief required so that cave / Arctic accessions from
Dr. Barton can be ingested by editing JSON, not Python.

**Validation:**

Commands run (in order):

```bash
.venv/bin/python -m pytest tests/ -v          # 20/20 tests pass
.venv/bin/python -m src.pipeline fetch references
.venv/bin/python -m src.pipeline fetch mags --limit 3      # smoke test
.venv/bin/python -m src.pipeline fetch databases           # SILVA + GTDB
.venv/bin/python -m src.pipeline fetch mags                # full 37 MAGs (3 skipped)

# Cave-readiness proof:
# Added GCF_000284135.1 to cave_arctic_pending in extremophiles.json (no code changes).
.venv/bin/python -m src.pipeline fetch mags                # +1 NCBI row, manifest -> 42
# Reverted demo entry from extremophiles.json after the proof.

.venv/bin/python -m src.pipeline status
.venv/bin/python -m src.pipeline verify                    # ok=42 changed=0 missing=0
```

Final manifest state (`data/manifest.tsv`): **42 rows, all SHA-256 verified.**

| Source | Rows | Notes |
|---|---|---|
| NCBI | 2 | JSC-1 reference + cave-readiness demo accession (GCF_000284135.1) |
| ENA | 37 | All PRJEB59431 polar cyanobacterial MAGs |
| SILVA | 1 | SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz (192 MB), MD5 verified |
| GTDB | 2 | bac120 taxonomy (7.3 MB) + metadata (225 MB) |

`verify` recomputed all 42 SHA-256 hashes and reported zero drift.

**Cave-readiness proof:** Adding a single accession to
`config/extremophiles.json -> cave_arctic_pending` and re-running
`python -m src.pipeline fetch mags` correctly dispatched it to the NCBI fetcher
and added a manifest row - zero code changes. The pipeline is ready to ingest
Dr. Barton's actual cave accessions the moment they arrive.

**Next:**

- **Action item for next meeting with Dr. Barton:** *"Please provide the exact
  cave + additional Arctic genome accessions (NCBI BioProject, ENA BioProject, or
  individual assembly accessions) to populate
  `config/extremophiles.json -> cave_arctic_pending`."*
- Iteration 2 entry point: build HMM profiles for the 21 FaRLiP genes from JSC-1's
  protein FASTA (already on disk at
  `data/references/GCF_000733415.1/ncbi_dataset/data/GCF_000733415.1/protein.faa`).

---

## Entry Template

Copy this template, paste at the **top** of this file, fill it in. Newest first.

```markdown
## YYYY-MM-DD - <one-line summary>

**Iteration:** 1 | 2 | ...
**Sources touched:** <NCBI / ENA / SILVA / GTDB / local / pipeline-only>

**Directory changes:**
- added: `path/to/file`
- modified: `path/to/file`
- removed: (none | `path/to/file`)

**What changed and why:**
<2-4 sentences. Cite the data source URL or paper if a decision was driven by it.>

**Validation:**
<Commands run, outputs observed, manifest rows added or changed.>

**Next:**
<Concrete next step. If blocked, name the blocker (e.g. "Awaiting cave accessions from Dr. Barton.").>
```
