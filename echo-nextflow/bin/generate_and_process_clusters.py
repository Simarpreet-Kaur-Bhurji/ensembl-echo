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

Note: run_mmseqs() is a standalone utility kept for reference; the Nextflow
pipeline drives MMseqs2 via echo_run_mmseqs.sh instead.
"""

import glob
import os
import subprocess
import time
from collections import defaultdict

import pandas as pd
import duckdb



# ---------------------------------------------------------------------------
# MMseqs2 wrapper (standalone utility; not called by NF pipeline)
# ---------------------------------------------------------------------------


def run_mmseqs(
    input_fasta,
    output_dir,
    min_seq_id=0.75,
    coverage=0.8,
    cov_mode=1,
    threads=16,
    singularity_image="mmseqs2_latest.sif",
):
    """
    Run MMseqs2 easy-cluster via Singularity and return the cluster TSV path.

    Not called by the Nextflow pipeline (which uses echo_run_mmseqs.sh).
    Kept as a standalone utility for ad-hoc use.
    """
    os.makedirs(output_dir, exist_ok=True)

    # MMseqs2 easy-cluster requires: input, output prefix, tmp dir
    output_prefix = os.path.join(output_dir, "mmseqs_results")
    tmp_dir = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    cmd = [
        "singularity",
        "exec",
        singularity_image,
        "mmseqs",
        "easy-cluster",
        input_fasta,
        output_prefix,
        tmp_dir,
        "--min-seq-id",
        str(min_seq_id),
        "-c",
        str(coverage),
        "--cov-mode",
        str(cov_mode),
        "--threads",
        str(threads),
    ]

    print("Running command:")
    print(" ".join(cmd))

    start_time = time.time()
    subprocess.run(cmd, check=True)
    elapsed = time.time() - start_time
    minutes, seconds = divmod(int(elapsed), 60)

    runtime_message = (
        f"MMseqs2 clustering completed in {minutes} min {seconds} sec.\n"
        f"Results stored in: {output_dir}\n"
    )
    summary_file = os.path.join(output_dir, "cluster_summary.txt")
    mode = "a" if os.path.exists(summary_file) else "w"
    with open(summary_file, mode, encoding="utf-8") as f:
        f.write("\nCluster Summary Report\n======================\n\n")
        f.write(runtime_message)

    return f"{output_prefix}_cluster.tsv"


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
               s.tax_id, s.sequence, s.name,
               s.confidence_score, s.confidence_level, s.seq_len
        FROM clusters_flat AS c
        LEFT JOIN read_parquet('{sequences_parquet}') AS s
          ON c.header = s.header
        """
    ).fetchdf()

    # Step 4: write output
    print(f"[parse_cluster_file] clusters parquet head:\n{df_final.head(10)}")
    df_final.to_parquet(output_parquet, index=False)
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


def get_singleton_sequences(df, output_dir):
    """
    Write clusters of size 1 (singletons) to FASTA and a summary TSV.

    Singletons are discarded from the main relatives search but optionally
    included in the final output via params.with_singletons.
    """
    singletons = df[df["cluster_size"] == 1]
    print(f"[get_singleton_sequences] {len(singletons)} singleton proteins")

    fasta_file = os.path.join(output_dir, "discarded_singletons.fa")
    write_all_sequences_to_fasta(singletons, fasta_file)

    log_file = os.path.join(output_dir, "singleton_cluster_summary.tsv")
    with open(log_file, "w", encoding="utf-8") as log:
        log.write("Cluster_ID\tProtein_id\tTax_ID\n")
        for _, row in singletons.iterrows():
            log.write(f"{row['Cluster_ID']}\t{row['header']}\t{row['tax_id']}\n")

    print(f"[get_singleton_sequences] wrote {fasta_file} and {log_file}")


def get_clusters_with_fewer_taxids(df, num_relatives, output_dir):
    """
    Write multi-protein clusters that have fewer unique tax_ids than num_relatives.

    These clusters can't yield the requested number of relatives, so they are
    pulled out separately and optionally added back to final output.
    """
    few_taxid = df[(df["cluster_size"] > 1) & (df["unique_tax_ids"] < num_relatives)]
    print(
        f"[get_clusters_with_fewer_taxids] {len(few_taxid)} proteins in few-taxid clusters"
    )

    fasta_file = os.path.join(output_dir, "clusters_with_fewer_tax_ids.fa")
    write_all_sequences_to_fasta(few_taxid, fasta_file)

    # one row per cluster summarising its members
    few_taxid_grouped = (
        few_taxid.groupby("Cluster_ID")
        .agg(
            {
                "header": ",".join,
                "tax_id": lambda x: ",".join(map(str, x)),
                "cluster_size": "first",
                "unique_tax_ids": "first",
            }
        )
        .reset_index()
    )

    log_file = os.path.join(output_dir, "clusters_with_fewer_taxids_summary.tsv")
    with open(log_file, "w", encoding="utf-8") as log:
        log.write("Cluster_ID\tProteins\tTax_ids\tTotal_Proteins\tTotal_TaxIDs\n")
        for _, row in few_taxid_grouped.iterrows():
            log.write(
                f"{row['Cluster_ID']}\t{row['header']}\t{row['tax_id']}\t"
                f"{row['cluster_size']}\t{row['unique_tax_ids']}\n"
            )

    print(f"[get_clusters_with_fewer_taxids] wrote {fasta_file} and {log_file}")


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


def process_clusters(output_dir, num_relatives):
    """
    Annotate clusters, write singleton / few-taxid side outputs, and return
    the remaining clusters that are eligible for the relatives search.

    Expects clusters*.parquet to exist in output_dir (written by parse_cluster_file).
    Writes to output_dir: 
      - discarded_singletons.fa / singleton_cluster_summary.tsv
      - clusters_with_fewer_tax_ids.fa / clusters_with_fewer_taxids_summary.tsv

    Returns:
        DataFrame: clusters with cluster_size > 1 AND unique_tax_ids >= num_relatives
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

    get_singleton_sequences(df, output_dir)
    get_clusters_with_fewer_taxids(df, num_relatives, output_dir)

    # clusters eligible for relatives search: multi-member with enough distinct species
    remaining_clusters = df[
        (df["cluster_size"] > 1) & (df["unique_tax_ids"] >= num_relatives)
    ]
    print(
    f" [process_clusters] remaining clusters: "
    f" {remaining_clusters['Cluster_ID'].nunique()} clusters, {len(remaining_clusters)} proteins"
    )
    return remaining_clusters
