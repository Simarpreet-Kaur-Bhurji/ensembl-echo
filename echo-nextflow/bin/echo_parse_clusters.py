#!/usr/bin/env python3
import argparse
import duckdb

from generate_and_process_clusters import parse_cluster_file


def main():
    p = argparse.ArgumentParser(description="ECHO P3: parse mmseqs cluster TSV -> clusters.parquet")
    p.add_argument("--cluster_tsv", required=True)
    p.add_argument("--processed_input_parquet", required=True)
    p.add_argument("--out_parquet", default="clusters.parquet")
    args = p.parse_args()

    con = duckdb.connect()

    # parse_cluster_file writes to output_dir/clusters.parquet
    # so we set output_dir='.' and then rename if needed
    parse_cluster_file(
        raw_tsv=args.cluster_tsv,
        sequences_parquet=args.processed_input_parquet,
        con=con,
        output_dir="."
    )

    # parse_cluster_file always writes ./clusters.parquet
    # If user asked for a different name, rename it
    if args.out_parquet != "clusters.parquet":
        import os
        os.rename("clusters.parquet", args.out_parquet)

    print(f"[OK] wrote {args.out_parquet}")


if __name__ == "__main__":
    main()

