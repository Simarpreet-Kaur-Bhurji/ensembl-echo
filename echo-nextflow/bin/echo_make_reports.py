#!/usr/bin/env python3
import argparse
import os
import pandas as pd


def count_fasta_headers(path: str) -> int:
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def count_tsv_rows(tsv_path: str) -> int:
    """Count data rows (excluding header) in a TSV."""
    with open(tsv_path) as fh:
        # skip header
        next(fh, None)
        return sum(1 for _ in fh)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters_parquet", required=True)
    ap.add_argument("--remaining_clusters_parquet", required=True)
    ap.add_argument("--clusters_with_fewer_taxids_summary_tsv", required=True)
    ap.add_argument("--discarded_singletons_fa", required=True)
    ap.add_argument("--clusters_with_fewer_tax_ids_fa", required=True)

    ap.add_argument("--all_relatives_fastas", nargs="+", required=True)
    ap.add_argument("--input_fasta", required=True)  # combined_input_fasta.fa

    ap.add_argument("--with_singletons", action="store_true")

    ap.add_argument("--out_cluster_summary", required=True)
    ap.add_argument("--out_pipeline_summary", required=True)
    args = ap.parse_args()

    # -------------------------
    # Cluster summary
    # -------------------------
    df = pd.read_parquet(args.clusters_parquet)
    total_clusters = int(df["Cluster_ID"].nunique())

    cluster_sizes = df.groupby("Cluster_ID").size()
    singleton_clusters = int((cluster_sizes == 1).sum())

    fewer_taxid_clusters = int(count_tsv_rows(args.clusters_with_fewer_taxids_summary_tsv))

    # Load remaining clusters to get actual closest relatives clusters
    rdf = pd.read_parquet(args.remaining_clusters_parquet)
    closest_rel_clusters = int(rdf["Cluster_ID"].nunique())

    with open(args.out_cluster_summary, "w") as f:
        f.write("Cluster Summary Report\n")
        f.write("======================\n\n")
        f.write(f"Total number of clusters obtained from input fasta: {total_clusters}\n")
        f.write(f"Singleton clusters (clusters with only one protein): {singleton_clusters}\n")
        f.write(f"Number of clusters with fewer than the required tax IDs: {fewer_taxid_clusters}\n")
        f.write(f"Closest relatives clusters (clusters with required tax IDs or more): {closest_rel_clusters}\n")

    # -------------------------
    # ECHO pipeline summary
    # -------------------------
    total_input_seqs = count_fasta_headers(args.input_fasta)

    # Common seq counts
    singletons_seq_count = count_fasta_headers(args.discarded_singletons_fa)
    fewer_taxids_seq_count = count_fasta_headers(args.clusters_with_fewer_tax_ids_fa)

    # Per-query retained counts
    retained_counts = {}
    for fa in args.all_relatives_fastas:
        base = os.path.basename(fa)
        # expects: <query>_all_relatives.fa
        q = base.replace("_all_relatives.fa", "")
        retained_counts[q] = count_fasta_headers(fa)

    with open(args.out_pipeline_summary, "w") as f:
        f.write("ECHO Pipeline Summary\n")
        f.write("====================\n\n")
        f.write(f"Total sequences in input FASTA: {total_input_seqs}\n\n")
        f.write("Common files for all queries:\n")
        f.write(f"  Discarded singletons: {singletons_seq_count} sequences\n")
        f.write(f"  Clusters with fewer tax IDs: {fewer_taxids_seq_count} sequences\n\n")
        f.write("Sequences retained per query:\n")
        for q in sorted(retained_counts.keys()):
            count = retained_counts[q]
            pct_ret = (count / total_input_seqs) * 100 if total_input_seqs else 0.0
            pct_dis = 100.0 - pct_ret
            f.write(f"  {q} : {count} sequences ({pct_ret:.2f}% retained, {pct_dis:.2f}% discarded)\n")
        f.write("\nNote: Each query_all_relatives.fa file includes sequences from clusters_with_fewer_tax_ids.fa.\n")

if __name__ == "__main__":
    main()

