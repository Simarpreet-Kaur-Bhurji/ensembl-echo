# ECHO Nextflow Pipeline

ECHO is a Nextflow (DSL2) pipeline that:
1. Combines per-species protein FASTA files into a single FASTA.
2. Builds metadata (`processed_input.tsv` + `processed_input.parquet`).
3. Runs **MMseqs2** clustering.
4. Parses cluster output into `clusters.parquet`.
5. Filters clusters (singletons, low taxon diversity, remaining clusters).
6. Computes ranked taxonomic distances (`ranked_taxa.tsv`).
7. Finds closest relatives per query species (chunked + parallel).
8. Produces merged logs, per-query FASTA outputs, diagnostics and summary reports.

---

## Repository layout

- `main.nf` – pipeline entrypoint
- `workflows/echo.nf` – main ECHO workflow DAG
- `modules/` – individual DSL2 process modules
- `bin/` – Python + shell wrappers executed by processes
- `test/` – test dataset + example params
- `nextflow.config` – profiles + resource configuration
- `nextflow_schema.json` – parameter schema (see [Parameter schema](#parameter-schema) below)
- `params.yaml` – typical run parameters

---

## Requirements

### 1) Nextflow
You need Nextflow installed and available on `PATH`.

Check:
```bash
nextflow -version
```

### 2) Python environment (venv + requirements.txt)
Create and activate a virtual environment, then install dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -r requirements.txt
```

The pipeline runs the scripts in bin/ using your active Python environment on the compute nodes.
Make sure your Slurm job environment loads the same Python/venv (or module) consistently.

### 3) MMseqs2 container (Singularity)
The pipeline runs MMseqs2 via Singularity. You must have:

- Singularity available on compute nodes
- A valid MMseqs2 image path set via `params.mmseqs_singularity_image`

### Inputs

#### FASTA directory (optional: `params.input_fasta_dir`)
An optional base directory used to resolve relative `file_path` entries in the metadata TSV.
If all `file_path` values in the metadata TSV are absolute, this parameter can be omitted.

#### Metadata TSV (`params.metadata_tsv`)
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

#### Query species TSV (`params.query_species`)
TSV with columns: `tax_id`, `sps_name`

Example:
```
tax_id  sps_name
78070   Platismatia glauca
...
```

### Configuration

#### Parameters (`params.yaml`)

**Required:**
- `metadata_tsv` – path to the metadata TSV
- `query_species` – path to the query species TSV
- `outdir` – output directory
- `mmseqs_singularity_image` – path to the MMseqs2 `.sif` image

**Clustering (MMseqs2):**
- `min_seq_id` – minimum sequence identity (default: `0.75`)
- `coverage` – minimum alignment coverage (default: `0.8`)
- `cov_mode` – coverage mode: `0`=bidirectional, `1`=target, `2`=query, `3`=target+query (default: `0`)
- `mmseqs_threads` – CPU threads for MMseqs2 (default: `32`)

**Relatives:**
- `num_of_rel` – number of closest relatives per query species (default: `5`)
- `clusters_per_chunk` – clusters per parallel chunk (default: `10000`)
- `with_singletons` – include singleton clusters in `_all_relatives.fa` output (default: `false`)

**FASTA output:**
- `fasta_line_width` – characters per sequence line in `combined_input_fasta.fa` (default: `60`). Standard FASTA is 60 or 80; adjust if downstream tools require a specific wrap length.

**Taxonomy:**
- `ncbi_taxa_db` – path to a pre-built NCBITaxa sqlite database. If omitted, ete3 uses its default (`~/.etetoolkit/taxa.sqlite`), downloading it on first use if absent. On a shared cluster, point this to a central copy to avoid repeated downloads and `~` filesystem pressure.

**Optional / advanced:**
- `input_fasta_dir` – base directory for resolving relative `file_path` entries in the metadata TSV
- `existing_clusters_dir` – path to a previous ECHO output directory; when set, steps 1–4 (FASTA parsing, MMseqs2, cluster filtering) are skipped and their outputs are read from this directory
- `workdir` – custom Nextflow work directory (default: `work` in the run location)

Example `params.yaml`:

```yaml
outdir: "results"
metadata_tsv: "/path/to/metadata.tsv"
query_species: "/path/to/query_species.tsv"
input_fasta_dir: "/path/to/input_fastas"   # optional if file_path in TSV is absolute

num_of_rel: 5
clusters_per_chunk: 10000

mmseqs_singularity_image: "/path/to/mmseqs2_latest.sif"
min_seq_id: 0.75
coverage: 0.8
cov_mode: 1
mmseqs_threads: 16

with_singletons: false
fasta_line_width: 60

# Point to a shared NCBITaxa DB to avoid ete3 re-downloading it
ncbi_taxa_db: "/shared/path/taxa.sqlite"

# existing_clusters_dir: "/path/to/previous/outdir"
# workdir: "/path/to/work/dir"
```

---

### Running

```bash
nextflow run main.nf -profile slurm -params-file params.yaml -with-report
```

#### Note on run-time reporting
The `-with-report` flag produces a Nextflow HTML report (e.g. `report.html`) in the run directory.
Open it in a browser to see per-process CPU, memory, and **wall-clock time** for each task.
This is the easiest way to check total run time and identify bottlenecks without parsing logs manually.

### Slurm (recommended)
Submit using a batch script like:

`run_nf.sh`
```bash
#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/vgp_set_nf.out
#SBATCH --error=logs/vgp_set_nf.err
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --job-name=echo_nf

set -euo pipefail

# Activate your python environment (edit as needed)
source /path/to/echo-nextflow/.venv/bin/activate

time nextflow run main.nf \
  -profile slurm \
  -params-file params.yaml \
  -with-report \
  -resume
```

Submit:

```bash
mkdir -p logs
sbatch run_nf.sh
```

---

### Parameter schema

`nextflow_schema.json` is a JSON Schema file that documents every pipeline parameter with its type, default, and description.
It can be used in several ways:

- **Validation** – Nextflow natively validates your `params.yaml` against the schema at startup and reports missing required parameters or type mismatches before any jobs run.
- **IDE auto-complete** – editors such as VS Code (with the nf-core schema extension) use it to provide parameter hints when editing `params.yaml`.
- **Documentation** – tools like `nf-core schema docs` can render the schema as a human-readable parameter table.

The schema is kept in sync with `nextflow.config`; if you add a new parameter, update both files.

---

### Outputs (in `params.outdir`)
Common key outputs:

- `combined_input_fasta.fa`
- `processed_input.tsv`, `processed_input.parquet`
- `mmseqs_results_cluster.tsv`
- `clusters.parquet`
- `remaining_clusters.parquet`
- `discarded_singletons.fa`
- `clusters_with_fewer_tax_ids.fa`
- `ranked_taxa.tsv`
- `closest_relatives_log.tsv`
- Per-query FASTA:
  - `<query_name>_all_relatives.fa` (adds common clusters-with-fewer-taxids and optionally singletons)
- Diagnostics:
  - `diagnostics_out/diagnostics.pdf`
  - `diagnostics_out/*.png`
- Nextflow report:
  - `report.html` (written where you run Nextflow unless configured otherwise)

#### Notes on reproducibility
Nextflow DAG + scripts are fully version-controlled.

The main reproducibility dependency is the Python environment and the MMseqs container.

For consistent runs, use:

- pinned `requirements.txt`
- fixed MMseqs `.sif` path/version
- fixed Nextflow version (optional but recommended)

#### Test run
A small test dataset is included under `test/`.
Example:

```bash
nextflow run main.nf -profile slurm -params-file test/params.test.yaml
```

### Troubleshooting
If a process fails, inspect:
- `.nextflow.log`
- the process work directory printed by Nextflow

Run the command manually:
```bash
cd <work/xx/xxxx>
bash .command.run
```

If outputs look stale, rerun without cache:

```bash
nextflow run main.nf -profile slurm -params-file params.yaml -resume
```
