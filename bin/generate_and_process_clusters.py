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
generate_and_process_clusters.py
---------------------------------
Library functions for ECHO pipeline steps 3–4 (clustering and cluster processing).

Functions used by the Nextflow pipeline:
  - parse_cluster_file()          – called by echo_parse_clusters.py (P3)
  - process_clusters()            – called by echo_process_clusters.py (P4)
  - annotate_clusters()           – internal helper for process_clusters()
  - write_all_sequences_to_fasta()– internal helper
  - get_singleton_sequences()     – internal helper
  - get_clusters_with_fewer_taxids() – internal helper

Note: The Nextflow pipeline drives MMseqs2 via echo_run_mmseqs.sh.
"""

import glob
import os
from collections import defaultdict

import pandas as pd
import duckdb



# ---------------------------------------------------------------------------
# Cluster parsing (P3)
# ---------------------------------------------------------------------------


def parse_cluster_file(raw_tsv, sequences_parquet, con, output_dir):
    """
    Parse the MMseqs2 cluster TSV and join with the processed-input parquet.

    MMseqs2 easy-cluster produces a two-column TSV: (seed, member).
    The seed (LHS) is included in the member list on the RHS, so we collect
    only RHS values to avoid double-counting.

    Steps:
      1. Build a seed → members mapping from the TSV.
      2. Flatten to a (Cluster_ID, header) dataframe.
      3. LEFT JOIN with processed_input.parquet to enrich with tax_id, sequence, etc.
      4. Write clusters.parquet to output_dir.

    Returns:
        str: path to the written clusters.parquet
    """
    output_parquet = os.path.join(output_dir, "clusters.parquet")
    cluster_pool = defaultdict(set)

    # Step 1: read TSV, group members under each seed
    with open(raw_tsv, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            lhs, rhs = parts[0], parts[1]
            cluster_pool[lhs].add(rhs)

    # Step 2: flatten to (Cluster_ID, header) records
    records = []
    for i, (_seed, members) in enumerate(cluster_pool.items(), start=1):
        cluster_id = f"Cluster_{i}"
        for protein in sorted(members):
            records.append({"Cluster_ID": cluster_id, "header": protein})

    df_flat = pd.DataFrame.from_records(records)

    # Step 3: enrich with metadata from processed_input.parquet (joined on header)
    con.register("clusters_flat", df_flat)
    df_final = con.execute(
        f"""
        SELECT c.Cluster_ID, c.header,
               s.tax_id, s.sequence, s.name, s.seq_len
        FROM clusters_flat AS c
        LEFT JOIN read_parquet('{sequences_parquet}') AS s
          ON c.header = s.header
        """
    ).fetchdf()

    # Step 4: write output — use duckdb to avoid pyarrow dependency
    print(f"[parse_cluster_file] clusters parquet head:\n{df_final.head(10)}")
    con.register("clusters_final", df_final)
    con.execute(f"COPY clusters_final TO '{output_parquet}' (FORMAT PARQUET)")
    print(f"[parse_cluster_file] saved {output_parquet}")
    return output_parquet


# ---------------------------------------------------------------------------
# FASTA output helpers
# ---------------------------------------------------------------------------


def write_all_sequences_to_fasta(df, fasta_file):
    """
    Write all sequences from df to a FASTA file (80-char line wrap).

    Args:
        df         : DataFrame with 'header' and 'sequence' columns
        fasta_file : output file path (str)
    """
    with open(fasta_file, "w", encoding="utf-8") as out:
        for _, row in df.iterrows():
            header = row["header"]
            sequence = row["sequence"]
            if sequence is not None:
                out.write(f">{header}\n")
                for i in range(0, len(sequence), 80):
                    out.write(sequence[i : i + 80] + "\n")
            else:
                print(
                    f"[write_all_sequences_to_fasta] WARNING: no sequence for {header}"
                )


def get_singleton_sequences(df, output_dir, with_singletons):
    """
    Write clusters of size 1 (singletons) to FASTA and either a manifest TSV
    (when with_singletons=True) or a summary TSV (when with_singletons=False).

    with_singletons=True  — singletons are included in the final per-query FASTA
      and manifest. singleton_manifest.tsv carries the same column schema as the
      closest-relatives manifest; query_tax_id and query_name are written as 'NA'
      placeholders that MERGE_MANIFESTS_PER_QUERY fills in per query.
      distance and selection_rank are 'NA' — singletons are not ranked by distance.

    with_singletons=False — singletons are discarded; singleton_cluster_summary.tsv
      records what was removed for diagnostics.
    """
    singletons = df[df["cluster_size"] == 1]
    print(f"[get_singleton_sequences] {len(singletons)} singleton proteins")

    fasta_file = os.path.join(output_dir, "discarded_singletons.fa")
    write_all_sequences_to_fasta(singletons, fasta_file)

    if with_singletons:
        manifest_rows = [
            {
                "query_tax_id":   "NA",
                "query_name":     "NA",
                "cluster_id":     row["Cluster_ID"],
                "protein_header": row["header"],
                "source_tax_id":  str(row["tax_id"]),
                "distance":       "NA",
                "selection_rank": "NA",
                "cluster_size":   1,
                "unique_tax_ids": 1,
                "selection_source": "singleton",
            }
            for _, row in singletons.iterrows()
        ]
        out_file = os.path.join(output_dir, "singleton_manifest.tsv")
        pd.DataFrame(manifest_rows).to_csv(out_file, sep="\t", index=False)
    else:
        out_file = os.path.join(output_dir, "singleton_cluster_summary.tsv")
        with open(out_file, "w", encoding="utf-8") as log:
            log.write("Cluster_ID\tProtein_id\tTax_ID\n")
            for _, row in singletons.iterrows():
                log.write(f"{row['Cluster_ID']}\t{row['header']}\t{row['tax_id']}\n")

    print(f"[get_singleton_sequences] wrote {fasta_file} and {out_file}")



# ---------------------------------------------------------------------------
# Cluster annotation and filtering (P4)
# ---------------------------------------------------------------------------


def annotate_clusters(con, clusters_table):
    """
    Add cluster_size and unique_tax_ids window columns to the clusters table.

    Args:
        con            : open DuckDB connection
        clusters_table : name of the registered table to annotate

    Returns:
        DataFrame with extra columns: cluster_size, unique_tax_ids
    """
    query = f"""
        WITH annotated AS (
            SELECT
                *,
                COUNT(*)           OVER (PARTITION BY Cluster_ID) AS cluster_size,
                COUNT(DISTINCT tax_id) OVER (PARTITION BY Cluster_ID) AS unique_tax_ids
            FROM {clusters_table}
        )
        SELECT * FROM annotated
    """
    return con.execute(query).fetchdf()


def process_clusters(output_dir, num_relatives, with_singletons):
    """
    Annotate clusters, write singleton side outputs, and return all multi-member
    clusters eligible for the relatives search.

    Expects clusters*.parquet to exist in output_dir (written by parse_cluster_file).
    Writes to output_dir:
      - discarded_singletons.fa / singleton_cluster_summary.tsv

    Returns:
        DataFrame: all clusters with cluster_size > 1 (SQL caps selection at
        LEAST(unique_tax_ids, num_relatives) per cluster)
    """
    parquet_files = glob.glob(os.path.join(output_dir, "clusters*.parquet"))
    if not parquet_files:
        raise FileNotFoundError(f"No cluster parquet file found in {output_dir}")

    output_cluster_name = parquet_files[0]
    print(f"[process_clusters] using cluster file: {output_cluster_name}")

    con = duckdb.connect()
    con.execute(
        f"CREATE OR REPLACE TABLE clusters AS SELECT * FROM read_parquet('{output_cluster_name}')"
    )

    df = annotate_clusters(con, "clusters")
    print(f"[process_clusters] annotated clusters head:\n{df.head()}")

    get_singleton_sequences(df, output_dir, with_singletons=with_singletons)

    # all multi-member clusters enter the SQL ranking path;
    # LEAST(unique_tax_ids, num_relatives) in the query caps selection per cluster
    remaining_clusters = df[df["cluster_size"] > 1]
    print(
    f" [process_clusters] remaining clusters: "
    f" {remaining_clusters['Cluster_ID'].nunique()} clusters, {len(remaining_clusters)} proteins"
    )
    return remaining_clusters
