"""
closest_taxonomic_relatives_module.py

This module provides functionality to identify the closest taxonomic relatives within clusters for target species 
and retrieve their protein sequences. It interacts with a DuckDB table to query protein data and writes the results 
to FASTA files for downstream analysis.

Functions:
- get_sequences_from_pep: Retrieves the sequence for a given protein ID from the DuckDB table and writes it to an output file.
- get_sequences_for_fasta: Retrieves sequences for a list of protein IDs and writes them to a FASTA file.
- get_distance_for_all_proteins: Calculates the taxonomic distances for all proteins in a cluster relative to a target species.
- get_closest_rel_within_cluster: Identifies the closest taxonomic relatives within clusters for each target species 
  and retrieves their protein sequences.

Dependencies:
- pandas: Used for handling data retrieved from the DuckDB table.
- duckdb: Used for querying the DuckDB table containing protein data.
- os: Used for file and directory management.

Output:
- FASTA files containing the protein sequences of the closest taxonomic relatives for each target species.
- Another FASTA file containing the sequences of all proteins in clusters with fewer than the specified number of relatives.
"""


import os
import pandas as pd


def get_sequences_from_pep(prot_id, con, outfile):
    """
    Get the sequence for a given protein ID from the DuckDB table and write it to an output file.

    Args:
        prot_id (str): The protein ID.
        con (duckdb.DuckDBPyConnection): The DuckDB connection object.
        outfile (file object): The output file object.

    Returns:
        tuple: The header and sequence that match the input protein ID.
    """
    # Protein from the input csv file
    row = con.execute(
        "SELECT * FROM hcp_table WHERE protein_id = ?", (prot_id,)
    ).fetchdf()
    if not row.empty:
        # Access values in the row
        header = row["header"].values[0]
        sequence = row["sequence"].values[0]
        print(header, sequence)
        outfile.write(f">{header}\n{sequence}\n")

        return header, sequence
    else:
        print(f"No match found for protein ID: {prot_id}")
        return None


def get_sequences_for_fasta(proteins, con, output_dir, file_name):
    """
    Retrieve sequences for the given proteins and write them to the output file.

    Args:
        proteins (list): List of proteins.
        con (duckdb.DuckDBPyConnection): The DuckDB connection object.
        output_dir (str): Path to the output directory.

    Returns:
        None
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_file = os.path.join(output_dir, file_name)

    with open(output_file, "a") as outfile:
        for i in proteins:
            split_header = i.split("|")
            prot_id = split_header[0]
            # print(prot_id)
            seq = get_sequences_from_pep(prot_id, con, outfile)


def get_distance_for_all_proteins(
    unique_cluster_protein_tax_ids, target_tax_id, ranked_taxa
):
    """
    Retrieve the taxonomy IDs of the ranked closest relatives for a target taxon
    within a cluster based on precomputed taxonomic distances.

    Args:
        unique_cluster_protein_tax_ids (list): A list of unique taxon IDs for proteins in the cluster.
        target_tax_id (int): The taxonomy ID of the target species.
        ranked_taxa (dict): A dictionary where keys are tuples of taxon ID pairs (target_tax_id, other_tax_id),
                            and values are tuples containing:
                            - Distance (int): The taxonomic distance between the two taxa.
                            - LCA (int): The last common ancestor taxon ID.

    Returns:
        list: A list of taxonomy IDs representing the ranked closest relatives for the target taxon,
              sorted by taxonomic distance in ascending order.

    Example:
        Input:
            unique_cluster_protein_tax_ids = [9606, 10090, 10116]
            target_tax_id = 9606
            ranked_taxa = {
                (9606, 10090): (3, 9605),
                (9606, 10116): (5, 9604)
            }

    """
    distances = []
    for taxon in unique_cluster_protein_tax_ids:
        if int(taxon) != int(target_tax_id):
            for k, v in ranked_taxa.items():
                if taxon in k and target_tax_id in k:
                    distances.append((int(v[0]), v[1], taxon))

    sorted_distances = sorted(distances, key=lambda x: (x[0], int(x[2])))
    ranked_closest_relatives = [int(i[2]) for i in sorted_distances]
    print("Ranked closest relatives for target taxon", ranked_closest_relatives)
    return sorted_distances


def get_proteins_from_taxid(tax_id, proteins, query_name, output_path, con):
    """
    Get proteins from a given taxonomy ID and write their sequences to an output file.

    Args:
        tax_id (str): The taxonomy ID.
        proteins (list): List of proteins.
        target (str): The target identifier.
        output_path (str): Path to the output directory.
        con (duckdb.DuckDBPyConnection): The DuckDB connection object.

    Returns:
        None
    """
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    query = query_name.lower().replace(" ", "_")
    file_name = f"{query}_relatives.fa"
    file_path = os.path.join(output_path, file_name)
    mode = "a" if os.path.isfile(file_path) else "w"

    with open(file_path, mode) as file:
        for j in proteins:
            # Query hcp_table for this protein and tax_id
            split_header = j.split("|")
            prot_id = split_header[0]
            query = """
                SELECT header, sequence 
                FROM hcp_table 
                WHERE protein_id = ? AND tax_id = ?
            """
            result = con.execute(query, (prot_id, tax_id)).fetchdf()

            if not result.empty:
                header = result["header"].values[0]
                sequence = result["sequence"].values[0]
                file.write(f">{header}\n")
                # Wrap sequence at 80 characters per line
                for i in range(0, len(sequence), 80):
                    file.write(sequence[i : i + 80] + "\n")
            # seq = get_sequences_from_pep(prot_id, con, file)


def break_taxonomic_ties(sorted_distances, proteins, con):
    """
    Break taxonomic ties by sorting based on distance, sequence length, and optionally confidence level.

    Args:
        sorted_distances (list): A list of tuples containing distance, LCA, and tax_id.
        proteins (list): A list of protein headers.
        con (duckdb.DuckDBPyConnection): The DuckDB connection object.

    Returns:
        pd.DataFrame: A DataFrame containing the resolved ties.
    """
    distances_df = pd.DataFrame(sorted_distances, columns=["distance", "lca", "tax_id"])
    distances_df["tax_id"] = distances_df["tax_id"].astype(str)

    prot_ids = [p.split("|")[0].strip() for p in proteins]
    placeholders = ",".join("?" for _ in prot_ids)
    query = f"""
        SELECT protein_id, tax_id, confidence_level, seq_len
        FROM hcp_table
        WHERE protein_id IN ({placeholders})
    """

    protein_df = con.execute(query, prot_ids).fetchdf()
    protein_df["tax_id"] = protein_df["tax_id"].astype(str)

    # Merge the protein data with distances
    merged_df = pd.merge(protein_df, distances_df, on="tax_id", how="inner")

    # Handle confidence_level as optional
    if "confidence_level" in merged_df.columns:
        merged_df["confidence_level"] = merged_df["confidence_level"].fillna("NA")
        merged_df["confidence_level"] = pd.Categorical(
            merged_df["confidence_level"],
            categories=["high1", "high2", "high3", "NA"],
            ordered=True,
        )
        # Check if all confidence_level values are "NA"
        if all(merged_df["confidence_level"] == "NA"):
            # Sort by distance and sequence length only
            merged_df = merged_df.sort_values(
                by=["distance", "seq_len"], ascending=[True, False]
            )
        else:
            # Sort by distance, sequence length, and confidence level
            merged_df = merged_df.sort_values(
                by=["distance", "seq_len", "confidence_level"],
                ascending=[True, False, True],
            )
    else:
        # Sort by distance and sequence length only
        merged_df = merged_df.sort_values(
            by=["distance", "seq_len"], ascending=[True, False]
        )

    # Drop duplicates based on tax_id, keeping the first occurrence
    final_df = merged_df.drop_duplicates(subset="tax_id", keep="first")

    return final_df[
        ["protein_id", "tax_id", "distance", "seq_len", "confidence_level"]
    ]


def get_taxids_from_proteins(con, proteins):
    tax_ids = []
    for p in proteins:
        split_header = p.split("|")
        prot_id = split_header[0]
        tax_id = con.execute(
            "SELECT tax_id FROM hcp_table WHERE protein_id = ?", (prot_id,)
        ).fetchdf()
        tax_ids.extend(tax_id["tax_id"].tolist())

    return list(set(map(int, tax_ids)))


def get_closest_rel_within_cluster(
    number_of_relatives, target_species, clusters_dict, ranked_taxa, con, output_dir
):
    """
    Identify the closest taxonomic relatives within clusters for each target species and retrieve their protein sequences.

    Args:
        number_of_relatives (int): The number of closest relatives to identify for each target species.
        target_species (dict): A dictionary where keys are target species taxon IDs, and values are tuples containing:
                               - Species name (str)
                               - Production name (str)
        clusters_dict (dict): A dictionary where keys are cluster IDs, and values are lists of protein identifiers
                              (e.g., "protein_id_name_tax_id").
        ranked_taxa (dict): A dictionary where keys are target taxon IDs, and values are dictionaries containing
                            precomputed taxonomic distances for each target taxon and other taxa.
        con (duckdb.DuckDBPyConnection): The DuckDB connection object for querying protein sequences.
        output_dir (str): The directory where the output FASTA files will be saved.

    Returns:
        None: The function writes the closest relatives' protein sequences to FASTA files in the specified output directory.

    Workflow:
        1. For each cluster, extract the unique taxonomic IDs of the proteins in the cluster.
        2. If the number of unique taxonomic IDs is less than or equal to `number_of_relatives`, retrieve all protein sequences.
        3. Otherwise, for each target species:
            - Calculate the distances to all proteins in the cluster using `get_distance_for_all_proteins`.
            - Identify the closest relatives based on the specified `number_of_relatives`.
            - Retrieve the protein sequences for the closest relatives using `get_proteins_from_taxid`.
    """
    log_data = []

    for cluster_id, proteins in clusters_dict.items():
        # tax_id_of_protein = [int(p.split("|")[4]) for p in proteins]
        # unique_cluster_protein_tax_ids = list(dict.fromkeys(tax_id_of_protein))
        unique_cluster_protein_tax_ids = get_taxids_from_proteins(con, proteins)
        print("Cluster No.", cluster_id.split(" ")[1])
        if len(proteins) == 1:
            # print(proteins)
            file_name_singleton = "discarded_singletons.fa"
            write_singletons = get_sequences_for_fasta(
                proteins, con, output_dir, file_name_singleton
            )

        elif len(unique_cluster_protein_tax_ids) < number_of_relatives:
            file_name_fewer = "clusters_with_fewer_tax_ids.fa"
            to_write = get_sequences_for_fasta(
                proteins, con, output_dir, file_name_fewer
            )

        else:
            for target_tax_id, target_name in target_species.items():
                sorted_distances = get_distance_for_all_proteins(
                    unique_cluster_protein_tax_ids,
                    target_tax_id,
                    ranked_taxa[target_tax_id],
                )
                resolved_ties = break_taxonomic_ties(sorted_distances, proteins, con)
                print(resolved_ties)
                ranked_closest_relatives = resolved_ties["tax_id"].tolist()
                print(
                    f"{number_of_relatives} closest relatives of target {target_tax_id}, {target_name}: {ranked_closest_relatives[:number_of_relatives]}"
                )
                closest_relatives = ranked_closest_relatives[:number_of_relatives]
                for i in closest_relatives:
                    get_proteins_from_taxid(
                        i, proteins, target_name, output_dir, con
                    )

                log_data.append(
                    {
                        "cluster_id": cluster_id,
                        "target_species_name": target_name,
                        "target_taxonomy_id": target_tax_id,
                        "closest_relatives": closest_relatives,
                    }
                )

    # Write log data to a TSV file
    tsv_log_path = os.path.join(output_dir, "closest_relatives_log.tsv")
    with open(tsv_log_path, "w") as tsv_file:
        tsv_file.write(
            "Cluster_ID\tTarget_Species_Name\tTarget_Taxonomy_ID\tClosest_Relatives\n"
        )
        for entry in log_data:
            tsv_file.write(
                f"{entry['cluster_id']}\t{entry['target_species_name']}\t{entry['target_taxonomy_id']}\t{','.join(map(str, entry['closest_relatives']))}\n"
            )
