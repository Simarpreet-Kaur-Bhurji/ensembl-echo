#!/usr/bin/env python3

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

import argparse
import os
import pandas as pd


def count_fasta_headers(path: str) -> int:
    n = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters_parquet", required=True)
    ap.add_argument("--remaining_clusters_parquet", required=True)
    ap.add_argument("--discarded_singletons_fa", required=True)

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

    # Load remaining clusters (all multi-member clusters eligible for relatives search)
    rdf = pd.read_parquet(args.remaining_clusters_parquet)
    multi_member_clusters = int(rdf["Cluster_ID"].nunique())

    with open(args.out_cluster_summary, "w", encoding="utf-8") as f:
        f.write("Cluster Summary Report\n")
        f.write("======================\n\n")
        f.write(
            f"Total number of clusters obtained from input fasta: {total_clusters}\n"
        )
        f.write(
            f"Singleton clusters (clusters with only one protein): {singleton_clusters}\n"
        )
        f.write(
            f"Multi-member clusters eligible for relatives search: {multi_member_clusters}\n"
        )

    # -------------------------
    # ECHO pipeline summary
    # -------------------------
    total_input_seqs = count_fasta_headers(args.input_fasta)

    # Singleton seq count
    singletons_seq_count = count_fasta_headers(args.discarded_singletons_fa)

    # Per-query retained counts
    retained_counts = {}
    for fa in args.all_relatives_fastas:
        base = os.path.basename(fa)
        # expects: <query>_all_relatives.fa
        q = base.replace("_all_relatives.fa", "")
        retained_counts[q] = count_fasta_headers(fa)

    with open(args.out_pipeline_summary, "w", encoding="utf-8") as f:
        f.write("ECHO Pipeline Summary\n")
        f.write("====================\n\n")
        f.write(f"Total sequences in input FASTA: {total_input_seqs}\n")
        f.write(f"Discarded singletons: {singletons_seq_count} sequences\n\n")
        f.write("Sequences retained per query:\n")
        for q in sorted(retained_counts.keys()):
            count = retained_counts[q]
            pct_ret = (count / total_input_seqs) * 100 if total_input_seqs else 0.0
            pct_dis = 100.0 - pct_ret
            f.write(
                f"  {q} : {count} sequences ({pct_ret:.2f}% retained, {pct_dis:.2f}% discarded)\n"
            )


if __name__ == "__main__":
    main()
