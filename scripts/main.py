"""
main.py

This script processes the input target species, clusters and high confidence protein fasta file
to calculate taxonomic distances and identify closest taxonomic relatives within clusters. It uses several 
modules to handle specific tasks such as parsing input files, processing clusters, parsing hcp fasta file, 
ranking the taxons and obtaining the closest taxonomic relatives from every cluster.

Modules:
- rank_by_taxon_module: Takes as input a list of all pairs of taxonomy ids and calculates taxonomic distance.
- process_input_species_module: Processes input species data, these are the target species for which we want annotations.
- process_clusters_module: Parses cluster file.
- unique_taxon_ids_from_hcp_module: Parses HCP FASTA file and stores as a DuckDB table for easy access.
- closest_taxonomic_relatives_module: Identifies closest taxonomic relatives within clusters for every target species.

Example Command:
time isrun -m 24gb -t 4d python main.py \
    --num_of_rel 5 \
    --species_path /hps/software/users/ensembl/compara/sbhurji/modenv/ensembl/main/ensembl-echo/scripts/birds_clustering/birds_query.csv \
    --cluster_path /hps/software/users/ensembl/compara/sbhurji/modenv/ensembl/main/ensembl-echo/scripts/birds_clustering/alfatclust_birds_output/birds_cluster.txt \
    --input_hcp_fasta /hps/software/users/ensembl/compara/sbhurji/modenv/ensembl/main/ensembl-echo/scripts/birds_clustering/birds_hcp_replaced_headers.fa \
    --output_dir ranking_all
"""


import argparse
import duckdb
import pickle

from rank_by_taxon_module import get_all_input_species_combinations, calculate_taxonomic_distance
from process_input_species_module import get_input_sps
from process_clusters_module import parse_cluster_file
from unique_taxon_ids_from_hcp_module import parse_hcp_fasta, get_hcp_tax_ids, create_hcp_table
from closest_taxonomic_relatives_module import get_closest_rel_within_cluster


def main():
    """
    Main function to process input data and calculate taxonomic relationships.

    This function:
    - Parses command-line arguments.
    - Reads and processes cluster and species input files.
    - Parses HCP FASTA files.
    - Creates a DuckDB table for HCP taxonomic data for easy access.
    - Calculates the taxonomic distances and ranks the taxons.
    - Identifies the closest taxonomic relatives within clusters for each target species.

    Command-line Arguments:
    --num_of_rel: (int) Number of closest relatives to identify.
    --species_path: (str) Path to the input species file.
    --cluster_path: (str) Path to the cluster file.
    --input_hcp_fasta: (str) Path to the HCP FASTA file.
    --output_dir: (str) Directory to store output files (default: "output").
    """
    parser = argparse.ArgumentParser(
        description='Processes input species, clusters, and HCP data to calculate taxonomic relationships.'
        )
    parser.add_argument(
        '--num_of_rel', 
        type=int, 
        required=True,
        help='Number of closest taxonomic relatives to identify for each target species.'
        )
    parser.add_argument(
        '--species_path', 
        type=str, 
        required=True,
        help='Path to the input file containing target species taxonomy id, species name and  production name.'
        )
    parser.add_argument(
        '--cluster_path', 
        type=str, 
        required=True,
        help='Path to the cluster file obtained from alfatclust.'
        )
    parser.add_argument(
        '--input_hcp_fasta', 
        type=str, 
        required=True,
        help='Path to the High Confidence Protein FASTA file. ' \
        'This script assumes that the taxonomy id is added to the header of every sequence at the end.'
        )
    parser.add_argument(
        '--output_dir', 
        type=str, 
        default="output",
        help='Directory to store the output files. Defaults to "output".'
        )


    args = parser.parse_args()
    #print(args)
    if(args.num_of_rel):
        clusters = parse_cluster_file(args.cluster_path)
        target_sps = get_input_sps(args.species_path)
        con = duckdb.connect()
        hcp_tax_ids = parse_hcp_fasta(args.input_hcp_fasta)
        hcp_table = create_hcp_table(hcp_tax_ids, con)
        taxon_ids = get_hcp_tax_ids(con)
        #print(len(taxon_ids))
        combinations = get_all_input_species_combinations(target_sps, taxon_ids)
        #print(combinations)
        # Uncomment the following line to calculate taxonomic distance.
        ranked_taxa = calculate_taxonomic_distance(combinations, args.output_dir)

        # Once calculated, you can load the ranked taxa from a file for any following runs.
        #with open("ranked_taxa.json", "rb") as f:
        #    ranked_taxa = pickle.load(f)

        closest_neighbours = get_closest_rel_within_cluster(args.num_of_rel, target_sps, clusters, ranked_taxa, con, args.output_dir)

main()