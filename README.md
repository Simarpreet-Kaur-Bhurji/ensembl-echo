# ECHO Nextflow Pipeline

ECHO is a Nextflow (DSL2) pipeline that:
1. Combines per-species protein FASTA files into a single FASTA.
2. Builds metadata (`processed_input.tsv` + `processed_input.parquet`).
3. Runs **MMseqs2** clustering.
4. Parses cluster output into `clusters.parquet` and separates singleton clusters.
5. Filters clusters and computes ranked taxonomic distances (`ranked_taxa.tsv`).
6. Finds the closest relatives per query species (chunked + parallel).
7. Produces per-query FASTA outputs (`*_all_relatives.fa`) and provenance manifests (`*_manifest.tsv`).
8. Optionally deduplicates identical sequences across the merged relatives FASTA.
9. Generates diagnostics plots and summary reports.

Steps 1–4 can be skipped for new query species by providing `existing_clusters_dir` (see [Restart mode](#restart-mode-existing_clusters_dir)).

---

## Repository layout

- `main.nf` – pipeline entrypoint
- `workflows/echo.nf` – main ECHO workflow DAG
- `modules/` – individual DSL2 process modules
- `bin/` – Python + shell wrappers executed by processes
- `test/` – test dataset, example params, and test scripts
- `nextflow.config` – profiles + resource configuration
- `nextflow_schema.json` – parameter schema (see [Parameter schema](#parameter-schema) below)
- `params.yaml` – typical run parameters
- `Dockerfile` – container definition for the pipeline runtime environment
- `requirements.txt` – Python package dependencies

---

## Requirements

### 1) Nextflow
You need Nextflow installed and available on `PATH`.

Check:
```bash
nextflow -version
```

### 2) Container runtime

All pipeline tools (Python packages, MMseqs2, DuckDB) are packaged in a single Docker image published to the GitHub Container Registry:

```
ghcr.io/simarpreet-kaur-bhurji/echo-container:latest
```

The image is defined by the `Dockerfile` at the repository root and built from `requirements.txt`.

**On the cluster (Singularity):**
Pull the image once as a `.sif` file, then use `-profile singularity,slurm`:

```bash
singularity pull echo-container_latest.sif docker://ghcr.io/simarpreet-kaur-bhurji/echo-container:latest
```

**Locally (Docker):**
Pull the image and use `-profile docker,local`:

```bash
docker pull ghcr.io/simarpreet-kaur-bhurji/echo-container:latest
```

> **Note (macOS Apple Silicon / ARM64):** Local Docker testing is not fully supported on
> macOS ARM64. MMseqs2's `linux/arm64` binary does not include all process modes available
> in the AVX2 x86-64 build, causing certain pipeline steps to fail. Use the Singularity
> image on an x86-64 cluster for production runs.

### 3) Running without a container (pyenv)

If you cannot or do not want to use Docker/Singularity, install dependencies into a
pyenv virtual environment and run the pipeline with the `slurm` profile (which has no
container directive and uses whatever Python and `mmseqs` are on `PATH`).

```bash
# Install pyenv if not already available
curl https://pyenv.run | bash

# Add to your shell profile (~/.bashrc or ~/.zshrc), then restart your shell:
export PYENV_ROOT="$HOME/.pyenv"
export PATH="$PYENV_ROOT/bin:$PATH"
eval "$(pyenv init -)"

# Install Python 3.10 (matches the container)
pyenv install 3.10.14
pyenv local 3.10.14

# Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate

# Install Python dependencies
pip install -r requirements.txt

# Download the NCBI Taxonomy database (one-time)
python -c "from ete3 import NCBITaxa; NCBITaxa()"
```

Install MMseqs2 natively — download the pre-built binary for your platform from
the [MMseqs2 releases page](https://github.com/soedinglab/MMseqs2#installation) and
ensure `mmseqs` is on your `PATH`.

Run on the cluster without a container:

```bash
time nextflow run main.nf \
  -profile slurm \
  -params-file params.yaml \
  -with-report \
  -resume
```

### 4) NCBI Taxonomy database

The taxonomy-ranking step requires a local copy of the NCBI Taxonomy database managed by [ete3](https://etetoolkit.org/).

**Default location:** `~/.etetoolkit/taxa.sqlite`

Download it once with:

```bash
python -c "from ete3 import NCBITaxa; NCBITaxa()"
```

On the cluster, Singularity automatically bind-mounts `$HOME` into the container (`singularity.autoMounts = true` in `nextflow.config`), so the default path works without extra configuration. For Docker, pass the path explicitly via `ncbi_taxa_db` in your params YAML and bind-mount the directory (see [Parameters](#parameters-paramsyaml) below).

---

## Inputs

### FASTA directory (optional: `params.input_fasta_dir`)
An optional base directory used to resolve relative `file_path` entries in the metadata TSV.
If all `file_path` values in the metadata TSV are absolute, this parameter can be omitted.

### Metadata TSV (`params.metadata_tsv`)
4-column TSV with a header row:

| Column | Description |
|--------|-------------|
| `sps_name` | Normalised species name (e.g. `homo_sapiens`) |
| `taxon_id` | NCBI taxon ID |
| `gca` | Genome assembly accession, or `NA` |
| `file_path` | Path to the per-species protein FASTA — absolute, or a filename/relative path resolved against `input_fasta_dir` |

Example:
```
sps_name	taxon_id	gca	file_path
homo_sapiens	9606	GCA_000001405.29	homo_sapiens.fa
mus_musculus	10090	NA	/absolute/path/mus_musculus.fa
```

### Query species TSV (`params.query_species`)
2-column TSV with a header row listing the species for which relatives will be found.

| Column | Description |
|--------|-------------|
| `sps_name` | Normalised species name matching entries in the metadata TSV |
| `tax_id` | NCBI taxon ID |

Example:
```
sps_name	tax_id
homo_sapiens	9606
mus_musculus	10090
```

Each query species produces its own `*_all_relatives.fa` and `*_manifest.tsv` in `params.outdir`.

---

## Configuration

### Parameters (`params.yaml`)

**Required:**
- `metadata_tsv` – path to the metadata TSV
- `query_species` – path to the query species TSV
- `outdir` – output directory

**Clustering (MMseqs2):**
- `min_seq_id` – minimum sequence identity (default: `0.75`)
- `coverage` – minimum alignment coverage (default: `0.8`)
- `cov_mode` – coverage mode: `0`=bidirectional, `1`=target, `2`=query, `3`=target+query (default: `0`)

**Relatives:**
- `num_of_rel` – maximum number of closest relatives per query species per cluster (default: `5`). Selection uses a three-tier tiebreak: (1) lowest taxonomic distance, (2) longer sequence if distance is equal, (3) first row if both are equal. Clusters with fewer distinct taxa than `num_of_rel` return all available taxa.
- `clusters_per_chunk` – clusters per parallel chunk (default: `10000`)

**Singleton policy (`with_singletons`):**
- `with_singletons: true` *(default)* — singletons are included as low-support but retained coverage evidence. They are appended to each query's `_all_relatives.fa` and appear in `*_manifest.tsv` with `selection_source=singleton`.
- `with_singletons: false` — singletons are excluded from all per-query outputs. They are written to `discarded_singletons.fa` and summarised in `singleton_cluster_summary.tsv` for reference.

**Sequence deduplication (`dedup_sequences`):**
- `dedup_sequences: false` *(default)* — no sequence-level deduplication. Duplicate headers are always removed (first occurrence kept).
- `dedup_sequences: true` — after merging relatives, sequences that are identical across different headers are deduplicated. The first occurrence is kept. A `dedup_report.tsv` is published to `outdir` listing every dropped/kept pair with its sequence length. Enable this when input FASTAs may contain identical proteins from different species or assemblies.

**FASTA output:**
- `fasta_line_width` – characters per sequence line in `combined_input_fasta.fa` (default: `60`). Standard FASTA is 60 or 80; adjust if downstream tools require a specific wrap length.

**Taxonomy:**
- `ncbi_taxa_db` – path to a pre-built NCBITaxa sqlite database. If omitted, ete3 uses its default (`~/.etetoolkit/taxa.sqlite`), downloading it on first use if absent. On a shared cluster, point this to a central copy to avoid repeated downloads and `~` filesystem pressure.

**Optional / advanced:**
- `input_fasta_dir` – base directory for resolving relative `file_path` entries in the metadata TSV
- `existing_clusters_dir` – path to a previous ECHO output directory to enable restart mode (see [Restart mode](#restart-mode-existing_clusters_dir))
- `workdir` – custom Nextflow work directory (default: `work` in the run location)

### Example `params.yaml`

```yaml
outdir: "results"
metadata_tsv: "/path/to/metadata.tsv"
query_species: "/path/to/query_species.tsv"
input_fasta_dir: "/path/to/input_fastas"   # optional if file_path in TSV is absolute

num_of_rel: 5
clusters_per_chunk: 10000

min_seq_id: 0.75
coverage: 0.8
cov_mode: 1

with_singletons: true
dedup_sequences: false
fasta_line_width: 60

# Point to a shared NCBITaxa DB to avoid ete3 re-downloading it
ncbi_taxa_db: "/shared/path/taxa.sqlite"

# existing_clusters_dir: "/path/to/previous/outdir"
# workdir: "/path/to/work/dir"
```

---

## Running

### Standard run

**On the cluster (Singularity + SLURM):**
```bash
time nextflow run main.nf \
  -profile singularity,slurm \
  -params-file params.yaml \
  -with-report \
  -resume
```

**Locally (Docker):**
```bash
time nextflow run main.nf \
  -profile docker,local \
  -params-file params.yaml \
  -with-report \
  -resume
```

The `time` prefix prints total wall-clock time when the run completes.
The `-with-report` flag produces a Nextflow HTML report (`report.html`) in the run directory — open it in a browser to see per-process CPU, memory, and wall-clock time for each task.
The `-resume` flag reuses cached task outputs from previous runs so only changed or failed steps are re-executed.

### Slurm (recommended for cluster)

Submit using a batch script like:

`run_nf.sh`
```bash
#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/echo_nf.out
#SBATCH --error=logs/echo_nf.err
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --job-name=echo_nf

set -euo pipefail

time nextflow run main.nf \
  -profile singularity,slurm \
  -params-file params.yaml \
  -with-report \
  -resume
```

Submit:

```bash
mkdir -p logs
sbatch run_nf.sh
```

### Restart mode (`existing_clusters_dir`)

If you need to add a new query species to a dataset that has already been clustered, use
restart mode to skip the expensive MMseqs2 steps. Set `existing_clusters_dir` to the
`outdir` of a previous full run, update `query_species` to list only the new query, and
point `outdir` to a new directory for the restart outputs.

The pipeline will read `processed_input.parquet`, `clusters.parquet`,
`remaining_clusters.parquet`, and `combined_input_fasta.fa` from the existing directory
and run only the ranked-taxa, closest-relatives, manifest, and diagnostics steps.

Example `params.restart.yaml`:
```yaml
existing_clusters_dir: "/path/to/previous/outdir"
outdir: "/path/to/restart/outdir"
query_species: "/path/to/new_query.tsv"

num_of_rel: 5
clusters_per_chunk: 10000
with_singletons: true
dedup_sequences: false
ncbi_taxa_db: "/shared/path/taxa.sqlite"
fasta_line_width: 60
```

Run:
```bash
# On the cluster:
nextflow run main.nf -profile singularity,slurm -params-file params.restart.yaml -with-report

# Locally:
nextflow run main.nf -profile docker,local -params-file params.restart.yaml -with-report
```

The clustering artifacts (`combined_input_fasta.fa`, `clusters.parquet`, etc.) remain
in the original output directory and are not re-published to the restart outdir.

---

## Parameter schema

`nextflow_schema.json` is a JSON Schema file that documents every pipeline parameter with its type, default, and description.
It can be used in several ways:

- **Validation** – Nextflow natively validates your `params.yaml` against the schema at startup and reports missing required parameters or type mismatches before any jobs run.
- **IDE auto-complete** – editors such as VS Code (with the nf-core schema extension) use it to provide parameter hints when editing `params.yaml`.
- **Documentation** – tools like `nf-core schema docs` can render the schema as a human-readable parameter table.

The schema is kept in sync with `nextflow.config`; if you add a new parameter, update both files.

---

## Outputs

All outputs are written to `params.outdir`.

### Clustering outputs
| File | Description |
|------|-------------|
| `combined_input_fasta.fa` | All input proteins concatenated into one FASTA |
| `processed_input.tsv` / `processed_input.parquet` | Per-protein metadata table used throughout the pipeline |
| `mmseqs_results_cluster.tsv` | Raw MMseqs2 cluster assignment output |
| `clusters.parquet` | Parsed cluster table with protein sequences and metadata |
| `remaining_clusters.parquet` | Multi-member clusters eligible for relative selection |

### Singleton outputs
Depends on `with_singletons`:

| File | `with_singletons: true` | `with_singletons: false` |
|------|------------------------|--------------------------|
| `singleton_manifest.tsv` | Written — one row per singleton protein | Not written |
| `discarded_singletons.fa` | Not written | Written — singleton sequences excluded from query outputs |
| `singleton_cluster_summary.tsv` | Not written | Written — summary table of discarded singleton clusters |

### Per-query outputs
For each query species, two files are produced:

**`<query_name>_all_relatives.fa`**
FASTA of all selected relative proteins. Contains up to `num_of_rel` relatives per cluster from `remaining_clusters`, plus singleton proteins if `with_singletons=true`. The query species' own proteins are never included.

**`<query_name>_manifest.tsv`**
Provenance handoff file for Genebuild. One row per selected relative with the following columns:

| Column | Description |
|--------|-------------|
| `query_tax_id` | NCBI taxon ID of the query species |
| `query_name` | Normalised name of the query species |
| `cluster_id` | Cluster the relative was drawn from |
| `protein_header` | FASTA header of the selected protein |
| `source_tax_id` | NCBI taxon ID of the species the protein belongs to |
| `distance` | Taxonomic distance between query and source species |
| `selection_rank` | Rank of this relative within its cluster (1 = closest) |
| `cluster_size` | Total number of proteins in the cluster |
| `unique_tax_ids` | Number of distinct taxa in the cluster |
| `selection_source` | `ranked_cluster` or `singleton` |

Row count in the manifest always matches the sequence count in the corresponding FASTA.

### Deduplication output
| File | When present |
|------|-------------|
| `dedup_report.tsv` | Only when `dedup_sequences=true` and at least one duplicate sequence is found. Columns: `dropped_header`, `kept_header`, `sequence_length`. |

If `dedup_sequences=true` but no duplicate sequences exist, no report is written.

### Summary reports
| File | Description |
|------|-------------|
| `cluster_summary.txt` | Total cluster count, singleton count, remaining cluster count |
| `echo_pipeline_summary.txt` | Total input sequences, query species processed, relative counts |

### Diagnostics
| File | Description |
|------|-------------|
| `diagnostics_out/diagnostics.pdf` | All plots compiled into a single PDF |
| `diagnostics_out/diagnostics_summary.txt` | Text summary of diagnostics values |
| `diagnostics_out/plot1_cluster_size_distribution.png` | Distribution of cluster sizes |
| `diagnostics_out/plot2_unique_taxids_distribution.png` | Distribution of unique taxa per cluster |
| `diagnostics_out/plot3_cluster_retention_summary.png` | Clusters retained vs filtered |
| `diagnostics_out/plot4_taxonomic_distance_distribution.png` | Distribution of selected relative distances |
| `diagnostics_out/plot5_size_vs_taxid_scatter.png` | Cluster size vs unique taxon count |

### Nextflow runtime report
- `report.html` — written to the run directory when `-with-report` is used. Shows per-process CPU, memory, and wall-clock time.

### Notes on reproducibility
The main reproducibility dependencies are the container image and the NCBI taxonomy database. For consistent runs use:

- a pinned image digest instead of `:latest` (e.g. `ghcr.io/simarpreet-kaur-bhurji/echo-container@sha256:<digest>`)
- the same `taxa.sqlite` snapshot across runs, or a shared central copy
- fixed Nextflow version (optional but recommended)

---

## Testing

Unit tests run without Nextflow or MMseqs2. From the repository root:

```bash
# Install test dependencies (one-time) if required
pip install pytest pandas duckdb

# Selection logic unit tests
pytest test/test_selection.py -v --log-cli-level=INFO

# Deduplication unit tests
pytest test/test_dedup.py -v --log-cli-level=INFO -k "not pipeline_dedup"
```

Pipeline-level tests (end-to-end, dedup pipeline, restart mode) require a completed Nextflow run first. See [test/README.md](test/README.md) for step-by-step instructions.

---

## Troubleshooting

If a process fails, inspect:
- `.nextflow.log`
- the process work directory printed by Nextflow

Run the failed command manually:
```bash
cd <work/xx/xxxx>
bash .command.run
```

If outputs look stale, rerun with `-resume` to pick up from where the pipeline left off:

```bash
nextflow run main.nf -profile singularity,slurm -params-file params.yaml -resume
```
