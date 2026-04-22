# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
closest_taxonomic_relatives_module.py

This module identifies the closest taxonomic relatives within clusters for target species 
and retrieves their protein sequences. It interacts with a DuckDB table created from a 
Parquet file and writes results to FASTA files for downstream analysis.

Improvements:
- Uses *batched queries per cluster* instead of one query per protein → much faster.
- Keeps Pandas for convenient deduplication and sorting.
"""

import contextlib
import os
import pandas as pd


# -----------------------------
# SQL-based tie breaking
# -----------------------------
def break_taxonomic_ties(con, cluster_id, query_tax_id, number_of_relatives):
    """
    For a given cluster + query_tax_id, return closest relatives using SQL.
    Tie-breaking happens with ROW_NUMBER().
    """
    sql = f"""
    SELECT header, tax_id, seq_len,
           COALESCE(confidence_level, 'NA') AS confidence_level,
           sequence, distance
    FROM (
        SELECT c.Cluster_ID, c.header, c.tax_id, c.seq_len,
               c.confidence_level, c.sequence,
               r.distance,
               ROW_NUMBER() OVER (
                   PARTITION BY c.tax_id 
                   ORDER BY r.distance ASC, c.seq_len DESC,
                            CASE COALESCE(c.confidence_level, 'NA')
                                WHEN 'high1' THEN 1
                                WHEN 'high2' THEN 2
                                WHEN 'high3' THEN 3
                                ELSE 99
                            END ASC
               ) AS rn
        FROM remaining_clusters c
        JOIN ranked_taxa r
          ON (c.tax_id = r.input_tid AND r.query_tax_id = {query_tax_id})
        WHERE c.Cluster_ID = '{cluster_id}'
    )
    WHERE rn = 1
    ORDER BY distance ASC, seq_len DESC
    LIMIT {number_of_relatives}
    """
    return con.execute(sql).fetchdf()


def write_closest_sequences(closest_df, fasta_file):
    """
    Batch write all sequences from a DataFrame into FASTA.
    """
    lines = []
    for _, row in closest_df.iterrows():
        header = row["header"]
        seq = row["sequence"]
        lines.append(f">{header}\n")
        lines.append("\n".join(seq[i : i + 80] for i in range(0, len(seq), 80)))
        lines.append("\n")
    fasta_file.write("".join(lines))


def get_closest_rel_within_cluster(
    remaining_clusters_df, number_of_relatives, query_species, con, output_dir
):
    """
    Identify closest relatives cluster by cluster.
    Cluster-first loop (your design), but SQL handles distances & ties.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Register input
    con.register("remaining_clusters", remaining_clusters_df)

    # Group clusters once
    grouped = con.execute(
        """
        SELECT 
            Cluster_ID, 
            array_agg(header) AS proteins,
            array_agg(DISTINCT CAST(tax_id AS INT)) AS unique_tax_ids
        FROM remaining_clusters
        GROUP BY Cluster_ID
    """
    ).fetchdf()

    log_data = []

    # Open one fasta file per query using ExitStack so all are properly closed
    with contextlib.ExitStack() as stack:
        fasta_files = {}
        for query_name in query_species.values():
            query_file_name = query_name.lower().replace(" ", "_")
            file_name = f"{query_file_name}_relatives.fa"
            file_path = os.path.join(output_dir, file_name)
            fasta_files[query_name] = stack.enter_context(
                open(file_path, "w", encoding="utf-8")
            )

        # Cluster-first loop
        for _, row in grouped.iterrows():
            cluster_id = row["Cluster_ID"]
            protein_headers = row["proteins"].tolist()
            unique_cluster_tax_ids = row["unique_tax_ids"].tolist()

            for query_tax_id, query_name in query_species.items():
                closest_df = break_taxonomic_ties(
                    con, cluster_id, query_tax_id, number_of_relatives
                )

                # Write sequences to correct FASTA
                if not closest_df.empty:
                    write_closest_sequences(closest_df, fasta_files[query_name])

                    # Append log info
                    log_data.append(
                        {
                            "cluster_id": cluster_id,
                            "query_species_name": query_name,
                            "query_taxonomy_id": query_tax_id,
                            "closest_relatives": closest_df["tax_id"].tolist(),
                            "closest_proteins": closest_df["header"].tolist(),
                            "closest_distances": closest_df["distance"].tolist(),
                            "num_proteins": len(protein_headers),
                            "num_unique_tax_ids": len(unique_cluster_tax_ids),
                        }
                    )

    # Write log file
    pd.DataFrame(log_data).to_csv(
        os.path.join(output_dir, "closest_relatives_log.tsv"), sep="\t", index=False
    )
    print(f"Sequences and log written to {output_dir}")
