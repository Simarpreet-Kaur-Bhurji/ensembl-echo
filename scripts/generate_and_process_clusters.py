import pandas as pd
import os
import time
import subprocess
import duckdb
from collections import defaultdict


def run_mmseqs(
    input_fasta,
    output_dir,
    min_seq_id=0.75,
    coverage=0.8,
    cov_mode=1,
    threads=16,
    singularity_image="/hps/nobackup/flicek/ensembl/compara/jitender/containers/mmseqs2_latest.sif",
):
    """
    Run MMseqs2 easy-cluster using Singularity and return the cluster file path.
    """
    os.makedirs(output_dir, exist_ok=True)

    # MMseqs2 easy-cluster requires three arguments: input, output prefix, tmp
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
    end_time = time.time()

    elapsed_time = end_time - start_time
    minutes, seconds = divmod(int(elapsed_time), 60)

    runtime_message = (
        f"MMseqs2 clustering completed in {minutes} min {seconds} sec.\n"
        f"Results stored in: {output_dir}\n"
    )
    summary_file = os.path.join(output_dir, "cluster_summary.txt")

    # Append to cluster_summary.txt
    mode = "a" if os.path.exists(summary_file) else "w"
    with open(summary_file, mode) as f:
        f.write("\n")
        f.write("Cluster Summary Report\n")
        f.write("======================\n\n")
        f.write(runtime_message)
    # Return the full path to the cluster file
    cluster_file = f"{output_prefix}_cluster.tsv"
    return cluster_file


def parse_cluster_file(raw_tsv, sequences_parquet, con, output_dir):
    """
    Parse raw TSV using LHS as cluster seed.
    Each cluster includes all unique RHS proteins (LHS already included in RHS).
    Join with sequences_parquet to get tax_id.
    """
    output_parquet = os.path.join(output_dir, "clusters.parquet")
    cluster_pool = defaultdict(set)

    # Step 1: collect clusters
    with open(raw_tsv, "r") as f:
        for i, line in enumerate(f, start=1):
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            lhs, rhs = parts[0], parts[1]
            cluster_pool[lhs].add(rhs)  # only add RHS

    # Step 2: build flat dataframe
    records = []
    for i, (seed, members) in enumerate(cluster_pool.items(), start=1):
        cluster_id = f"#Cluster_{i}"
        for protein in sorted(members):
            records.append({"Cluster_ID": cluster_id, "header": protein})

    df_flat = pd.DataFrame.from_records(records)

    # Step 3: join with sequences_parquet to get tax_id
    con.register("clusters_flat", df_flat)
    df_final = con.execute(
        f"""
        SELECT c.Cluster_ID, c.header, s.tax_id, s.sequence, s.name, s.confidence_score, s.confidence_level, seq_len
        FROM clusters_flat AS c
        LEFT JOIN read_parquet('{sequences_parquet}') AS s
        ON c.header = s.header
    """
    ).fetchdf()

    # Step 4: save parquet
    print("Clusters parquet head", df_final.head(10))
    df_final.to_parquet(output_parquet, index=False)
    # con.close()
    print(f"Saved {output_parquet}")

    return output_parquet


def write_all_sequences_to_fasta(df, fasta_file):
    """
    Write sequences from df to FASTA file.
    """

    with open(fasta_file, "w") as fasta_file:
        for _, row in df.iterrows():
            header = row["header"]
            sequence = row["sequence"]
            fasta_file.write(f">{header}\n")
            for i in range(0, len(sequence), 80):
                fasta_file.write(sequence[i : i + 80] + "\n")


def get_singleton_sequences(df, output_dir):
    singletons = df[df["cluster_size"] == 1]
    print("Singletons", singletons.head())
    fasta_file = os.path.join(output_dir, "discarded_singletons.fa")
    write_all_sequences_to_fasta(singletons, fasta_file)

    log_file = os.path.join(output_dir, "singleton_cluster_summary.tsv")
    with open(log_file, "w") as log:
        log.write("Cluster_ID\tProtein_id\tTax_ID\n")
        for _, row in singletons.iterrows():
            log.write(f"{row['Cluster_ID']}\t{row['header']}\t{row['tax_id']}\n")

    print(f"Wrote {len(singletons)} singletons to {fasta_file} and {log_file}")


def get_clusters_with_fewer_taxids(df, num_relatives, output_dir):
    few_taxid = df[(df["cluster_size"] > 1) & (df["unique_tax_ids"] < num_relatives)]
    print("Few taxid clusters", few_taxid.head())
    fasta_file = os.path.join(output_dir, "clusters_with_fewer_tax_ids.fa")
    write_all_sequences_to_fasta(few_taxid, fasta_file)

    log_file = os.path.join(output_dir, "clusters_with_fewer_taxids_summary.tsv")

    few_taxid_grouped = (
        few_taxid.groupby("Cluster_ID")
        .agg(
            {
                "header": lambda x: ",".join(x),
                "tax_id": lambda x: ",".join(map(str, x)),
                "cluster_size": "first",
                "unique_tax_ids": "first",
            }
        )
        .reset_index()
    )

    with open(log_file, "w") as log:
        log.write("Cluster_ID\tProteins\tTax_ids\tTotal_Proteins\tTotal_TaxIDs\n")
        for _, row in few_taxid_grouped.iterrows():
            log.write(
                f"{row['Cluster_ID']}\t{row['header']}\t{row['tax_id']}\t{row['cluster_size']}\t{row['unique_tax_ids']}\n"
            )

    print(
        f"Wrote {len(few_taxid)} proteins from few-taxid clusters to {fasta_file} and {log_file}"
    )


def annotate_clusters(con, clusters_table):
    """Add cluster size and unique tax_id count."""
    query = f"""
        WITH annotated AS (
            SELECT 
                *,
                COUNT(*) OVER (PARTITION BY Cluster_ID) AS cluster_size,
                COUNT(DISTINCT tax_id) OVER (PARTITION BY Cluster_ID) AS unique_tax_ids
            FROM {clusters_table}
        )
        SELECT * FROM annotated
    """

    return con.execute(query).fetchdf()


def process_clusters(output_dir, num_relatives):
    output_cluster_name = os.path.join(output_dir, "clusters.parquet")
    con = duckdb.connect()
    con.execute(
        f"CREATE OR REPLACE TABLE clusters AS SELECT * FROM read_parquet('{output_cluster_name}')"
    )

    df = annotate_clusters(con, "clusters")
    print("Annotated clusters", df.head())

    get_singleton_sequences(df, output_dir)
    get_clusters_with_fewer_taxids(df, num_relatives, output_dir)

    remaining_clusters = df[
        (df["cluster_size"] > 1) & (df["unique_tax_ids"] >= num_relatives)
    ]

    return remaining_clusters
