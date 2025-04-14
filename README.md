# Ensembl Conserved set of HOmologues (ECHO)

🧬 ECHO – Evolutionary Clustering for High-confidence Orthologues
ECHO is a pipeline designed to identify the closest taxonomic relatives of target species using high-confidence protein clusters. It leverages taxonomy-aware ranking and efficient database querying to extract the most relevant sequences for comparative genomics.

## 🚀 Example Command

time isrun -m 24gb -t 4d python main.py \
    --num_of_rel 5 \
    --species_path path/to/input.csv \
    --cluster_path path/to/clusters.txt \
    --input_hcp_fasta path/to/hcp.fa \
    --output_dir output
    
## Required Parameters
--num_of_rel: Number of closest relatives to retrieve per target species.

--species_path: CSV file with target species information (tax ID, species name, production name).

--cluster_path: File containing cluster assignments (protein groupings).

--input_hcp_fasta: FASTA file of high-confidence protein sequences.

--output_dir: Directory where output FASTA files and results will be saved.

## 📦 Module Overview
The pipeline is organized into modular components, each handling a specific task:

1. process_clusters_module.py
Parses the cluster file (clusters.txt) and returns a dictionary:

{ cluster_number: [protein_header1, protein_header2, ...], ... }

3. process_input_species_module.py
Parses the input species metadata CSV and returns a dictionary:

{ tax_id: (species_name, production_name), ... }

3. unique_taxon_ids_from_hcp_module.py
Parses the input HCP FASTA file to extract unique taxonomic IDs.
Also creates a DuckDB table from the FASTA for fast and efficient protein sequence lookups.

4. rank_by_taxon_module.py
Generates all possible combinations of target species and HCP taxon IDs.
Calculates taxonomic distances between each species pairs. 

6. closest_taxonomic_relatives_module.py
Identifies the closest taxonomic relatives within clusters for input species.
Queries the DuckDB FASTA database to extract matching protein sequences.
Outputs results as individual FASTA files for each target species for downstream analysis.

## 📁 Output
The pipeline writes:

Ranked closest relative FASTA files per target species to --output_dir.

Intermediate data (e.g., taxonomic distances) as serialized files or logs.

## 🧪 Example Input Files

input.csv
tax_id	species_name	production_name
9606	Homo sapiens	homo_sapiens
10090	Mus musculus	mus_musculus

clusters.txt

#Cluster1
protein1
protein2
#Cluster2
protein3
protein4
...

## 📌 Notes

The taxonomy resolution step requires access to a local or cached NCBI taxonomy database.

The pipeline supports parallel and resource-managed execution via isrun, but it can be adapted to other batch environments.
