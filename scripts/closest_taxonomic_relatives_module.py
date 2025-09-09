"""
closest_taxonomic_relatives_module.py

This module identifies the closest taxonomic relatives within clusters for target species 
and retrieves their protein sequences. It interacts with a DuckDB table created from a 
Parquet file and writes results to FASTA files for downstream analysis.

Improvements:
- Uses *batched queries per cluster* instead of one query per protein → much faster.
- Keeps Pandas for convenient deduplication and sorting.
"""

import os
import pandas as pd


def get_sequences_from_pep(header, con, outfile):
    """
    Get the sequence for a given protein header from the DuckDB table and write it to an output file.

    Args:
        header (str): The full protein header from the cluster.
        con (duckdb.DuckDBPyConnection): The DuckDB connection object.
        outfile (file object): The output file object.

    Returns:
        tuple: The header and sequence that match the input header.
    """
    row = con.execute(
        "SELECT header, sequence FROM hcp_table WHERE header = ? LIMIT 1", (header,)
    ).fetchdf()

    if not row.empty:
        sequence = row["sequence"].values[0]
        outfile.write(f">{header}\n{sequence}\n")
        return header, sequence
    else:
        print(f"No match found for header: {header}")
        return None


def get_sequences_for_fasta(protein_headers, con, output_dir, file_name):
    """
    Retrieve sequences for the given protein headers and write them to a FASTA file.

    Args:
        proteins (list): List of full protein headers from the cluster.
        con (duckdb.DuckDBPyConnection): The DuckDB connection object.
        output_dir (str): Path to the output directory.
        file_name (str): Name of the FASTA file to write.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_file = os.path.join(output_dir, file_name)

    with open(output_file, "a") as outfile:
        for header in protein_headers:
            get_sequences_from_pep(header.strip(), con, outfile)


def get_distance_for_all_proteins(unique_cluster_tax_ids, query_tax_id, ranked_taxa):
    """
    Calculate distances for all proteins in a cluster relative to a target taxon.
    """
    distances = []
    for taxon in unique_cluster_tax_ids:
        if int(taxon) != int(query_tax_id):
            for k, v in ranked_taxa.items():
                if taxon in k and query_tax_id in k:
                    distances.append((int(v[0]), v[1], int(taxon)))

    # Do I need to sort by tax id here? Yes to make the code reproducible
    sorted_distances = sorted(distances, key=lambda x: (x[0], int(x[2])))
    # sorted_distances = sorted(distances, key=lambda x: x[0])
    return sorted_distances


def write_closest_relatives_fasta(closest_df, query_name, output_dir):
    print(f"Writing {len(closest_df)} sequences for query '{query_name}'")
    query_file_name = query_name.lower().replace(" ", "_")
    file_name = f"{query_file_name}_relatives.fa"
    file_path = os.path.join(output_dir, file_name)
    os.makedirs(output_dir, exist_ok=True)
    mode = "a" if os.path.isfile(file_path) else "w"

    with open(file_path, mode) as fasta_file:
        for _, row in closest_df.iterrows():
            header = row["header"]
            sequence = row["sequence"]
            fasta_file.write(f">{header}\n")
            # Wrap sequence at 80 characters per line
            for i in range(0, len(sequence), 80):
                fasta_file.write(sequence[i : i + 80] + "\n")


def break_taxonomic_ties(sorted_distances, protein_headers, con):
    """
    Break taxonomic ties using DuckDB with a temporary table for distances.
    Efficient for large clusters and avoids Pandas merge.

    Args:
        sorted_distances (list): List of tuples (distance, lca, tax_id)
        proteins (list): List of full protein headers in the cluster
        con (duckdb.DuckDBPyConnection): DuckDB connection object

    Returns:
        pd.DataFrame: DataFrame with resolved ties containing columns
                      ["header", "tax_id", "distance", "seq_len", "confidence_level"]
    """
    if not sorted_distances or not protein_headers:
        return pd.DataFrame(
            columns=["header", "tax_id", "distance", "seq_len", "confidence_level"]
        )

    # Create a temporary table for distances
    con.execute(
        "CREATE TEMP TABLE distance_table(distance INTEGER, lca INTEGER, tax_id INTEGER)"
    )
    con.executemany("INSERT INTO distance_table VALUES (?, ?, ?)", sorted_distances)

    # Prepare placeholders for headers
    headers_placeholders = ",".join("?" for _ in protein_headers)

    # SQL query: join hcp_table with temp distance_table
    sql = f"""
    SELECT h.header, h.tax_id, h.seq_len,
           COALESCE(h.confidence_level, 'NA') AS confidence_level,
           d.distance, h.sequence
    FROM hcp_table h
    JOIN distance_table d USING(tax_id)
    WHERE h.header IN ({headers_placeholders})
    ORDER BY d.distance ASC, h.seq_len DESC,
             CASE COALESCE(h.confidence_level, 'NA')
                 WHEN 'high1' THEN 1
                 WHEN 'high2' THEN 2
                 WHEN 'high3' THEN 3
                 ELSE 99
             END ASC
    """

    # Execute query
    merged_df = con.execute(sql, protein_headers).fetchdf()
    # print("Merged and sorted DataFrame:\n", merged_df)

    # print("Duplicates", merged_df[merged_df.duplicated(subset=["header"])])
    # Drop duplicates per tax_id (keep first)
    final_df = merged_df.drop_duplicates(subset="tax_id", keep="first")

    # Drop temporary table
    con.execute("DROP TABLE distance_table")

    return final_df[
        ["header", "tax_id", "seq_len", "confidence_level", "distance", "sequence"]
    ]


def get_taxids_from_proteins(protein_headers, con):
    """
    Get unique tax_ids for a list of proteins using full headers (not just protein IDs).
    """
    headers = [p.strip() for p in protein_headers]
    placeholders = ",".join("?" for _ in headers)
    query = f"SELECT DISTINCT tax_id FROM hcp_table WHERE header IN ({placeholders})"
    df = con.execute(query, headers).fetchdf()
    # print("Distinct tax ids in cluster:", df)
    return df["tax_id"].astype(int).tolist()


def get_closest_rel_within_cluster(
    number_of_relatives, query_species, clusters_dict, ranked_taxa, con, output_dir
):
    """
    Identify the closest taxonomic relatives within clusters and write results to FASTA + log.
    """
    log_data = []
    singleton_cluster_summary = []
    less_than_req_num_rel_cluster_summary = []

    # Counters
    singletons_count = 0
    small_clusters_count = 0
    valid_clusters_count = 0

    for cluster_id, protein_headers in clusters_dict.items():
        unique_cluster_tax_ids = get_taxids_from_proteins(protein_headers, con)
        print(
            f"Processing {cluster_id} with {len(protein_headers)} proteins and {len(unique_cluster_tax_ids)} tax ids."
        )

        if len(protein_headers) == 1:
            singletons_count += 1
            print(f"Skipping {cluster_id}: only one protein in cluster.")
            get_sequences_for_fasta(
                protein_headers, con, output_dir, "discarded_singletons.fa"
            )

            singleton_cluster_summary.append(
                {
                    "cluster_id": cluster_id,
                    "prot_id": protein_headers[0],
                    "tax_id": list(unique_cluster_tax_ids)[0]
                    if unique_cluster_tax_ids
                    else None,
                }
            )

        elif (
            len(unique_cluster_tax_ids) < number_of_relatives
            and len(protein_headers) > 1
        ):
            small_clusters_count += 1
            print(
                f" {cluster_id} has fewer than {number_of_relatives} taxon ids; saving all proteins."
            )
            get_sequences_for_fasta(
                protein_headers, con, output_dir, "clusters_with_fewer_tax_ids.fa"
            )

            less_than_req_num_rel_cluster_summary.append(
                {
                    "cluster_id": cluster_id,
                    "proteins": ",".join(protein_headers),
                    "tax_ids": ",".join(map(str, unique_cluster_tax_ids)),
                    "total_prot_ids": len(protein_headers),
                    "total_tax_ids": len(unique_cluster_tax_ids),
                }
            )

        elif len(unique_cluster_tax_ids) >= number_of_relatives:
            valid_clusters_count += 1
            print(
                f" {cluster_id} has more than {number_of_relatives} taxon ids; selecting {number_of_relatives} closest from them.."
            )
            for query_tax_id, query_name in query_species.items():
                print("Query tax id:", query_tax_id, query_name)
                sorted_distances = get_distance_for_all_proteins(
                    unique_cluster_tax_ids, query_tax_id, ranked_taxa[query_tax_id]
                )

                resolved_ties = break_taxonomic_ties(
                    sorted_distances, protein_headers, con
                )
                # print(" Resolved ties:\n", resolved_ties)

                # Take top N closest proteins
                closest_df = resolved_ties.head(number_of_relatives)
                print(f"Closest relatives DataFrame: closest_df of length {closest_df}")

                # Extract closest headers, sequences, tax_ids, and distances
                closest_proteins = closest_df["header"].tolist()
                closest_seqs = closest_df["sequence"].tolist()
                closest_relatives = closest_df["tax_id"].tolist()
                closest_distances = closest_df["distance"].tolist()

                print(
                    f"{number_of_relatives} closest relatives for query {query_tax_id} ({query_name}) are {closest_relatives} and their distances to the query are {closest_distances}"
                )
                write_closest_relatives_fasta(closest_df, query_name, output_dir)

                log_data.append(
                    {
                        "cluster_id": cluster_id,
                        "query_species_name": query_name,
                        "query_taxonomy_id": query_tax_id,
                        "closest_relatives": closest_relatives,
                        "closest_proteins": closest_proteins,
                        "closest_distances": closest_distances,
                        "num_proteins": len(protein_headers),
                        "num_unique_tax_ids": len(unique_cluster_tax_ids),
                    }
                )

    # Write log file

    tsv_log_path = os.path.join(output_dir, "closest_relatives_log.tsv")
    with open(tsv_log_path, "w") as tsv_file:
        # Updated header
        tsv_file.write(
            "Cluster_ID\tQuery_Species_Name\tQuery_Taxonomy_ID\t"
            "Closest_Relatives\tClosest_Proteins\tClosest_Distances\t"
            "Num_Proteins\tNum_Unique_TaxIDs\n"
        )
        for entry in log_data:
            tsv_file.write(
                f"{entry['cluster_id']}\t"
                f"{entry['query_species_name']}\t"
                f"{entry['query_taxonomy_id']}\t"
                f"{','.join(map(str, entry['closest_relatives']))}\t"
                f"{','.join(entry['closest_proteins'])}\t"
                f"{','.join(map(str, entry['closest_distances']))}\t"
                f"{entry['num_proteins']}\t"
                f"{entry['num_unique_tax_ids']}\n"
            )

    # --- Write singleton cluster summary ---
    singleton_summary_path = os.path.join(output_dir, "singleton_cluster_summary.tsv")
    with open(singleton_summary_path, "w") as f:
        f.write("Cluster_ID\tProtein_ID\tTax_ID\n")
        for entry in singleton_cluster_summary:
            f.write(f"{entry['cluster_id']}\t{entry['prot_id']}\t{entry['tax_id']}\n")

    print(f"Singleton cluster summary written to {singleton_summary_path}")

    # --- Write small clusters (less than required relatives) summary ---
    small_clusters_summary_path = os.path.join(
        output_dir, "clusters_with_fewer_taxids_summary.tsv"
    )
    with open(small_clusters_summary_path, "w") as f:
        f.write("Cluster_ID\tProteins\tTax_ids\tTotal_Proteins\tTotal_TaxIDs\n")
        for entry in less_than_req_num_rel_cluster_summary:
            f.write(
                f"{entry['cluster_id']}\t{entry['proteins']}\t{entry['tax_ids']}\t{entry['total_prot_ids']}\t{entry['total_tax_ids']}\n"
            )

    print(
        f"Clusters with fewer than required relatives summary written to {small_clusters_summary_path}"
    )

    print("\n=== Cluster Summary ===")
    print(f"Singleton clusters: {singletons_count}")
    print(f"Clusters with < {number_of_relatives} tax IDs : {small_clusters_count}")
    print(f"Clusters with >= {number_of_relatives} tax IDs: {valid_clusters_count}")
