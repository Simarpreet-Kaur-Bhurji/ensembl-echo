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

singularity available on compute nodes

a valid MMseqs image path in params (e.g. params.singularity_image.sif)

### Inputs
FASTA directory (params.input_fasta_dir)
A directory containing per-species protein FASTA files (*.fa), e.g.

```
input_fastas/
  species1.prots.fa
  species2.prots.fa
  ...
```

#### Metadata TSV (params.metadata_tsv)
TSV header: Scientific name  Species taxon_id

#### Query species TSV (params.query_species)
TSV header: tax_id  sps_name

Example:
```
tax_id  sps_name
78070   Platismatia glauca
...
```
### Configuration

Parameters (params.yaml)
You typically set:
- input_fasta_dir
- metadata_tsv
- query_species
- outdir
- num_of_rel
- chunk_size
- singularity_image
- MMseqs settings: min_seq_id, coverage, cov_mode, threads
- optional: with_singletons

Example snippet:

```
outdir: "results"
input_fasta_dir: "/path/to/input_fastas"
metadata_tsv: "/path/to/metadata.tsv"
query_species: "/path/to/query_species.tsv"

num_of_rel: 5
chunk_size: 2500

singularity_image: "/path/to/mmseqs2_latest.sif"
min_seq_id: 0.75
coverage: 0.8
cov_mode: 1
threads: 16

with_singletons: false
```

### Running

Local (small tests)
```
nextflow run main.nf -profile local -params-file params.yaml -with-report
```

### Slurm (recommended)
Submit using a batch script like:

run_nf.sh
```
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

time
```

Submit:

```
mkdir -p logs
sbatch run_nf.sh
```

### Outputs (in params.outdir)
Common key outputs:

- combined_input_fasta.fa
- processed_input.tsv, processed_input.parquet
- mmseqs_results_cluster.tsv
- clusters.parquet
- remaining_clusters.parquet
- discarded_singletons.fa
- clusters_with_fewer_tax_ids.fa
- ranked_taxa.tsv
- closest_relatives_log.tsv
- Per-query FASTA:
  <query_name>_relatives.fa
  <query_name>_all_relatives.fa (adds common clusters-with-fewer-taxids and optionally singletons)
- Diagnostics:
  diagnostics_out/diagnostics.html
  diagnostics_out/*.png
- Nextflow report:
  report.html (written where you run Nextflow unless configured otherwise)

#### Notes on reproducibility
Nextflow DAG + scripts are fully version-controlled.

The main reproducibility dependency is the Python environment and the MMseqs container.

For consistent runs, use:

- pinned requirements.txt
- fixed MMseqs .sif path/version
- fixed Nextflow version (optional but recommended)

#### Test run
A small test dataset is included under test/.
Example:

```
nextflow run main.nf -profile local -params-file test/params.test.yaml
```

Troubleshooting
If a process fails, inspect:
.nextflow.log

the process work directory printed by Nextflow

run the command manually:
```
cd <work/xx/xxxx>
bash .command.run
```
If outputs look stale, rerun without cache:

```
nextflow run main.nf -profile slurm -params-file params.yaml -resume
```

