"""
main.py

This script processes input protein sequences from either a High Confidence Protein (HCP) FASTA file 
or a directory of FASTA files, generates clusters using MMseqs2, and identifies the closest 
taxonomic relatives for a set of query species within each cluster.

Functionality includes:
- Combining multiple FASTA files into a single file with standardized headers.
- Generating a metadata TSV and Parquet file with sequence information, scientific names, and tax IDs.
- Running MMseqs2 to cluster proteins and parsing cluster outputs into dictionaries.
- Calculating ranked taxa based on taxonomic distances.
- Mapping query species to the closest relatives within protein clusters.
- Appending clusters with fewer tax IDs to all relative output files.

Command-line arguments:
- --num_of_rel: Number of closest taxonomic relatives to identify per species.
- --query_species: Path to query species file (tax_id, species name, production name).
- --output_dir: Directory to store outputs (default: "output").
- --input_hcp_fasta / --input_fasta_dir: Input HCP FASTA file or path to directory of FASTA files.
- --metadata_tsv: Required if using a FASTA directory; provides species metadata (Scientific name, Taxon id).
- MMseqs2 options: --min_seq_id, --coverage, --cov_mode, --threads, --singularity_image

This script relies on several helper modules for processing metadata, parsing FASTA files, 
generating clusters, and calculating taxonomic distances.
"""

import argparse
import duckdb
import pickle
import os
import time

from rank_by_taxon_module import (
    get_all_input_species_combinations,
    calculate_taxonomic_distance,
)
from process_input_species_module import get_input_sps
from generate_input_tsv import *
from parse_hcp_fasta import parse_hcp_fasta
from parse_metadata_parquet import create_hcp_table, get_hcp_tax_ids
from generate_clusters import run_mmseqs, parse_cluster_file
from closest_taxonomic_relatives_module import get_closest_rel_within_cluster
from diagnostics import get_diagnostic_stats_and_plots


def validate_args(args):
    """
    Validate command-line arguments and ensure required dependencies exist.

    Raises:
        ValueError: If --input_fasta_dir is provided without --metadata_tsv.
    """
    if args.input_fasta_dir and not args.metadata_tsv:
        raise ValueError("--metadata_tsv is required if --input_fasta_dir is provided")
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)


def process_fasta_dir(args):
    """
    Process a directory of FASTA files:
    - Combine FASTA files into one.
    - Generate metadata TSV with protein info.

    Returns:
        tuple: (combined_fasta_path, metadata_tsv_path)
    """
    combine_fasta_name = os.path.join(args.output_dir, "combined_input_fasta.fa")
    combined = combine_fastas(args.input_fasta_dir, combine_fasta_name)
    species_map = load_species_info(args.metadata_tsv)
    processed_input_tsv = os.path.join(args.output_dir, "processed_input.tsv")
    write_protein_metadata(combined, species_map, processed_input_tsv)
    return combine_fasta_name, processed_input_tsv


def process_hcp_fasta(args):
    """
    Process a High Confidence Protein (HCP) FASTA file:
    - Generate metadata TSV and Parquet from the HCP FASTA.

    Returns:
        str: Path to the processed metadata TSV.
    """
    output_tsv_name = os.path.join(args.output_dir, "processed_input.tsv")
    parse_hcp = parse_hcp_fasta(args.input_hcp_fasta, output_tsv_name)
    return output_tsv_name


def write_elapsed_time(start_time, end_time, output_dir):
    elapsed_time = end_time - start_time
    minutes, seconds = divmod(int(elapsed_time), 60)

    runtime_message = (
        f"Closest relatives for all queries were computed in {minutes} min {seconds} sec.\n"
        f"Results, diagnostic stats and plots are stored in: {output_dir}\n"
    )

    echo_file = os.path.join(output_dir, "echo_pipeline_summary.txt")

    # Append to cluster_summary.txt
    mode = "a" if os.path.exists(echo_file) else "w"
    with open(echo_file, mode) as f:
        if mode == "w":
            f.write("ECHO Pipeline Summary\n")
            f.write("======================\n\n")
        f.write(runtime_message)

    return runtime_message.strip()


def load_or_calculate_ranked_taxa(combinations, output_dir):
    """
    Load previously calculated ranked taxa from JSON if it exists,
    otherwise calculate and return taxonomic distances.

    Returns:
        dict: Ranked taxa mapping for input species combinations.
    """
    rank_taxon_file = os.path.join(output_dir, "ranked_taxa.pkl")
    if os.path.exists(rank_taxon_file):
        with open(rank_taxon_file, "r") as f:
            return pickle.load(f)
    return calculate_taxonomic_distance(combinations, output_dir)


def append_to_relatives(output_dir):
    """
    1. For each '*_relatives.fa', create a new '*_all_relatives.fa' file
       that contains the relatives plus clusters_with_fewer_tax_ids.fa.
    2. Create one combined_all_relatives.fa file containing all *_relatives.fa
       and clusters_with_fewer_tax_ids.fa (only once).
    3. Make hidden all original *_relatives.fa files and clusters_with_fewer_tax_ids.fa
       by prefixing them with a dot, but leave newly created files visible.
    """
    clusters_file = os.path.join(output_dir, "clusters_with_fewer_tax_ids.fa")
    if not os.path.exists(clusters_file):
        raise FileNotFoundError(f"{clusters_file} not found!")

    # Read clusters_with_fewer_tax_ids content once
    with open(clusters_file, "r") as cf:
        cluster_content = cf.read()

    # Prepare all_queries_combined_relatives.fa
    combined_path = os.path.join(output_dir, "all_queries_relatives_combined.fa")
    with open(combined_path, "w") as combined_out:
        for fname in os.listdir(output_dir):
            if fname.endswith("_relatives.fa") and not fname.endswith(
                "_all_relatives.fa"
            ):
                fpath = os.path.join(output_dir, fname)

                # Read original relatives content
                with open(fpath, "r") as rf:
                    relatives_content = rf.read()

                # Write *_all_relatives.fa
                all_relatives_path = os.path.join(
                    output_dir, fname.replace("_relatives.fa", "_all_relatives.fa")
                )
                with open(all_relatives_path, "w") as arf:
                    arf.write(
                        relatives_content.strip()
                        + "\n"
                        + cluster_content.strip()
                        + "\n"
                    )
                print(f"Created {all_relatives_path}")

                # Append to combined_all_relatives.fa
                combined_out.write(relatives_content.strip() + "\n")

        # Finally, append clusters_with_fewer_tax_ids content once
        combined_out.write(cluster_content.strip() + "\n")

    print(f"Created combined_all_relatives.fa: {combined_path}")

    # Hide only the original *_relatives.fa files and clusters_with_fewer_tax_ids.fa
    for fname in os.listdir(output_dir):
        fpath = os.path.join(output_dir, fname)
        if (
            (
                (
                    fname.endswith("_relatives.fa")
                    and not fname.endswith("_all_relatives.fa")
                )
                or fname == "clusters_with_fewer_tax_ids.fa"
            )
            and not fname.startswith(".")
            and fname != "combined_input_fasta.fa"
        ):
            hidden_path = os.path.join(output_dir, "." + fname)
            os.rename(fpath, hidden_path)
            print(f"Hidden {fname} -> {os.path.basename(hidden_path)}")


def main():
    """
    Main pipeline for processing input species data:
    1. Parse command-line arguments.
    2. Process either HCP FASTA or a directory of FASTA files.
    3. Run MMseqs2 clustering if needed.
    4. Create DuckDB table from metadata.
    5. Calculate ranked taxa for closest relative search.
    6. Identify closest taxonomic relatives within clusters.
    7. Append clusters with fewer tax IDs to relative files.
    """
    parser = argparse.ArgumentParser(
        description="Processes input species, clusters, and HCP data to calculate taxonomic relationships."
    )
    parser.add_argument(
        "--num_of_rel",
        type=int,
        required=True,
        help="Number of closest taxonomic relatives to identify for each target species.",
    )
    parser.add_argument(
        "--query_species",
        type=str,
        required=True,
        help="Path to the file containing query species taxonomy id, species name and production name.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="output",
        help='Directory to store the output files. Defaults to "output".',
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--input_hcp_fasta",
        type=str,
        help="Path to High Confidence Protein FASTA file (taxonomy id in headers).",
    )
    group.add_argument(
        "--input_fasta_dir",
        type=str,
        help="Directory containing FASTA files for processing. Provide metadata_tsv if using this option.",
    )
    parser.add_argument(
        "--metadata_tsv",
        type=str,
        help="Species metadata TSV (required if --input_fasta_dir is used).",
    )
    parser.add_argument(
        "--min_seq_id", type=float, default=0.75, help="Minimum sequence identity"
    )
    parser.add_argument(
        "--coverage", type=float, default=0.8, help="Coverage threshold"
    )
    parser.add_argument("--cov_mode", type=int, default=1, help="Coverage mode")
    parser.add_argument("--threads", type=int, default=16, help="Number of threads")
    parser.add_argument(
        "--singularity_image",
        default="/hps/nobackup/flicek/ensembl/compara/jitender/containers/mmseqs2_latest.sif",
        help="Path to MMseqs Singularity image",
    )

    args = parser.parse_args()

    # Validate arguments
    validate_args(args)

    # Process input data
    if args.input_fasta_dir:
        combine_fasta_name, processed_input_tsv = process_fasta_dir(args)
        cluster_file = os.path.join(args.output_dir, "mmseqs_results_cluster.tsv")
        if os.path.exists(cluster_file):
            clusters = parse_cluster_file(cluster_file)
        else:
            cluster_file = run_mmseqs(
                combine_fasta_name,
                args.output_dir,
                min_seq_id=args.min_seq_id,
                coverage=args.coverage,
                cov_mode=args.cov_mode,
                threads=args.threads,
                singularity_image=args.singularity_image,
            )
            clusters = parse_cluster_file(cluster_file)
    elif args.input_hcp_fasta:
        processed_input_tsv = process_hcp_fasta(args)
        cluster_file = os.path.join(args.output_dir, "mmseqs_results_cluster.tsv")
        if os.path.exists(cluster_file):
            clusters = parse_cluster_file(cluster_file)
        else:
            cluster_file = run_mmseqs(
                args.input_hcp_fasta,
                args.output_dir,
                min_seq_id=args.min_seq_id,
                coverage=args.coverage,
                cov_mode=args.cov_mode,
                threads=args.threads,
                singularity_image=args.singularity_image,
            )
            clusters = parse_cluster_file(cluster_file)

    # Connect to DuckDB and process metadata
    start_time = time.time()
    con = duckdb.connect()
    metadata_pq_file = next(
        (
            os.path.join(args.output_dir, f)
            for f in os.listdir(args.output_dir)
            if f.endswith(".parquet")
        ),
        None,
    )
    create_hcp_table(metadata_pq_file, con)
    query_sps = get_input_sps(args.query_species)
    taxon_ids = get_hcp_tax_ids(con)

    # Generate combinations and calculate ranked taxa
    combinations = get_all_input_species_combinations(query_sps, taxon_ids)
    ranked_taxa = load_or_calculate_ranked_taxa(combinations, args.output_dir)

    # Get closest relatives
    get_closest_rel_within_cluster(
        args.num_of_rel, query_sps, clusters, ranked_taxa, con, args.output_dir
    )
    append_to_relatives(args.output_dir)

    if args.input_hcp_fasta:
        get_diagnostic_stats_and_plots(args.input_hcp_fasta, args.output_dir)
    elif args.input_fasta_dir:
        files = os.listdir(args.output_dir)
        created_fasta_file = os.path.join(args.output_dir, "combined_input_fasta.fa")
        get_diagnostic_stats_and_plots(created_fasta_file, args.output_dir)

    end_time = time.time()
    write_elapsed_time(start_time, end_time, args.output_dir)


if __name__ == "__main__":
    main()
