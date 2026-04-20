# ECHO Test Suite

Tests live in `echo-nextflow/test/` and are split into two levels:

- **Unit tests** — run instantly with no pipeline, no MMseqs2, no NCBI DB required.
- **Pipeline tests** — validate real Nextflow output; require a pipeline run first.

All tests use pytest. Logs are written automatically to `echo-nextflow/test/logs/` at the
end of every session (`info_<timestamp>.log` and `debug_<timestamp>.log`).

---

## Prerequisites

```bash
pip install pytest pytest-repeat pandas duckdb
```

### NCBI Taxonomy database (`taxa.sqlite`)

The pipeline uses [ete3](https://etetoolkit.org/) to resolve taxonomic distances.
ete3 requires a local copy of the NCBI Taxonomy database as a SQLite file.

**Format:** SQLite file (`taxa.sqlite`), created and managed by ete3.  
**Version:** ete3 fetches the current NCBI taxonomy dump at the time of first use — there is no explicit version number; the file reflects the NCBI taxonomy tree as it was when you downloaded it.  
**Default location:** `~/.etetoolkit/taxa.sqlite`

To pre-download the database before your first pipeline run:

```bash
python -c "from ete3 import NCBITaxa; NCBITaxa()"
```

This downloads ~50 MB from the NCBI FTP and writes `taxa.sqlite` to `~/.etetoolkit/`.

**Using the database inside the container**

- **Singularity** (recommended on the cluster): `singularity.autoMounts = true` is set in `nextflow.config`, so Singularity automatically bind-mounts `$HOME` from the host into the container. The default path `~/.etetoolkit/taxa.sqlite` is therefore accessible inside the container without any extra configuration.
- **Docker** (local use): Docker does not auto-mount the host home directory. Either bind-mount the directory explicitly (`-v $HOME/.etetoolkit:/root/.etetoolkit`) or set `ncbi_taxa_db` in your params YAML to an absolute path and mount it as a Docker volume. The simplest approach is to set `ncbi_taxa_db` to an absolute host path and add `--volume` to `docker.runOptions` in `nextflow.config`.

---

## 1. Selection logic — unit tests (`test_selection.py`)

No pipeline required. Runs in a few seconds.

```bash
pytest echo-nextflow/test/test_selection.py -v --log-cli-level=INFO
```

**What is tested:**

| Test | Description |
|---|---|
| `test_more_taxa_than_n_selects_exactly_n` (9 cases) | Cluster with more taxa than N → exactly N selected by lowest distance |
| `test_fewer_taxa_than_n_retains_all_parametrized` (9 cases) | Cluster with fewer taxa than N → all retained, no cap applied |
| `test_equal_taxa_as_n` (5 cases) | Cluster with exactly N taxa → all retained |
| `test_no_matching_taxa_produces_empty_outputs` | No ranked taxa match → empty FASTA, zero-row manifest with correct column schema |
| `test_best_protein_per_taxon_selected` (5 cases) | Three-tier tiebreak: (1) lowest distance, (2) longer sequence, (3) first row |

Optional flag: `--cluster-id-mode random` assigns each protein a distinct cluster ID
(default is `same`, all proteins share one cluster).

---

## 2. Deduplication — unit tests (`test_dedup.py`)

No pipeline required. Runs in a few seconds.

```bash
pytest echo-nextflow/test/test_dedup.py -v --log-cli-level=INFO -k "not pipeline_dedup"
```

**What is tested:**

- Duplicate header removed; first occurrence kept; unique headers all retained.
- Sequence dedup disabled by default; exact duplicates removed when `--dedup_sequences` is passed.
- `dedup_report.tsv` written only when sequence duplicates exist; correct schema and content.
- No report written when there are no sequence duplicates.
- Cross-file header deduplication and cross-file sequence deduplication.
- Real pipeline sequences: header dup, sequence dup, combined header+sequence dups, report content.

---

## 3. End-to-end pipeline test (`test_e2e.py`)

Requires Nextflow, the NCBI taxonomy DB (see [Prerequisites](#prerequisites) above), and either the Singularity SIF or Docker image for the container.

**Step 1** — Run the standard pipeline (once, or whenever test data changes):

```bash
cd echo-nextflow

# On the cluster (Singularity + SLURM):
./nextflow run main.nf -params-file test/params.test.yaml -profile singularity,slurm
./nextflow run main.nf -params-file test/params.test.no_singletons.yaml -profile singularity,slurm

# Locally with Docker:
./nextflow run main.nf -params-file test/params.test.yaml -profile docker,local
./nextflow run main.nf -params-file test/params.test.no_singletons.yaml -profile docker,local
```

These produce `nextflow_test/` and `nextflow_test_no_singletons/` respectively.

**Step 2** — Run the e2e tests:

```bash
pytest echo-nextflow/test/test_e2e.py -v --log-cli-level=INFO
```

**What is tested:**

| Class | Tests |
|---|---|
| `TestClustering` | Total cluster count (7), total protein count (36), remaining cluster count (6), singleton count (1), singleton absent from remaining_clusters |
| `TestPerQueryOutputs` | FASTA exists and has 30 sequences, manifest exists with all 10 required columns and 30 rows, query species not self-selected, `selection_rank` ≤ `num_of_rel`, `query_tax_id` consistent, `selection_source` values correct |
| `TestSummaryReports` | `cluster_summary.txt` and `echo_pipeline_summary.txt` content checks |
| `TestDiagnostics` | `diagnostics_out/` directory and all 7 files (PDF, summary, 5 plots) exist and are non-empty |
| `TestNoSingletons` | Same checks against the `with_singletons=false` run: `singleton_cluster_summary.tsv` written, `singleton_manifest.tsv` absent, `discarded_singletons.fa` present, per-query counts reduced by 1 |

**Synthetic dataset properties (fixed-seed, deterministic):**

- 6 species, 36 proteins, 6 protein families
- 7 clusters after MMseqs2: 6 multi-member + 1 singleton
- 2 query species: `saccharomyces_cerevisiae` (4932), `candida_albicans` (5476)
- `num_of_rel = 5`; expected 30 relatives per query (29 ranked-cluster + 1 singleton)

---

## 4. Deduplication pipeline test (`test_dedup.py` — pipeline section)

Uses a variant of the test dataset where two deliberate sequence duplicate pairs have been
introduced into the input FASTAs:

- `neu03` and `sch03` share the same sequence (length 51)
- `yar02` and `asp02` share the same sequence (length 61)

**Step 1** — Run the dedup pipeline:

```bash
cd echo-nextflow

# On the cluster (Singularity + SLURM):
./nextflow run main.nf -params-file test/params.test.dedup.yaml -profile singularity,slurm

# Locally with Docker:
./nextflow run main.nf -params-file test/params.test.dedup.yaml -profile docker,local
```

This produces `nextflow_test_dedup/`.

**Step 2** — Run the full dedup test file (unit + pipeline):

```bash
pytest echo-nextflow/test/test_dedup.py -v --log-cli-level=INFO
```

**What is tested (pipeline section):**

| Test | Description |
|---|---|
| `test_pipeline_dedup_fasta_exists` | Output FASTA present for each query |
| `test_pipeline_dedup_sequence_count` | 28 sequences per query (30 base − 2 duplicate pairs removed) |
| `test_pipeline_dedup_output_is_clean` | Re-running dedup on pipeline output finds zero residual duplicates |
| `test_pipeline_dedup_report_exists` | `dedup_report.tsv` published to output directory |
| `test_pipeline_dedup_report_content` | 2 rows: pairs (neu03, sch03, len=51) and (yar02, asp02, len=61) |

---

## 5. Restart mode test (`test_restart.py`)

Tests the `existing_clusters_dir` parameter — running a new query species against
pre-computed clusters from a previous full run, without repeating MMseqs2 clustering.

**Step 1** — Ensure `nextflow_test/` exists (from step 3 above).

**Step 2** — Run the restart pipeline with a new query species (`schizosaccharomyces_pombe`):

```bash
cd echo-nextflow

# On the cluster (Singularity + SLURM):
./nextflow run main.nf -params-file test/params.test.restart.yaml -profile singularity,slurm

# Locally with Docker:
./nextflow run main.nf -params-file test/params.test.restart.yaml -profile docker,local
```

This produces `nextflow_test_restart/`.

**Step 3** — Run the restart tests:

```bash
pytest echo-nextflow/test/test_restart.py -v --log-cli-level=INFO
```

**What is tested:**

| Test | Description |
|---|---|
| `test_new_query_all_relatives_fasta_exists / _nonempty` | FASTA for `schizosaccharomyces_pombe` produced and non-empty |
| `test_new_query_manifest_exists / _columns / _row_count_matches_fasta` | Manifest produced with all 10 required columns; rows match FASTA sequence count |
| `test_clustering_artifacts_not_republished` | Heavy clustering outputs absent from restart outdir (they remain in the base dir) |
| `test_original_query_outputs_absent` | Original query FASTAs not re-created in restart outdir |
| `test_base_dir_clustering_artifacts_intact` | Base run outputs untouched by the restart |
| `test_diagnostics_files_present` | Full diagnostics bundle (7 files) produced for the new query |

---

## 6. Full suite in one command

Run all tests (unit tests pass immediately; pipeline tests skip gracefully if the relevant
output directory does not exist):

```bash
pytest echo-nextflow/test/ -v --log-cli-level=INFO
```

---

## Logs

Every pytest session writes two log files to `echo-nextflow/test/logs/`:

- `info_<timestamp>.log` — INFO-level summary (mirrors console output)
- `debug_<timestamp>.log` — full DEBUG trace across all tests

These are useful for diagnosing failures without re-running the pipeline.
