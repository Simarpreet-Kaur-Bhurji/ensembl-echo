#!/usr/bin/env python3
import argparse
import pandas as pd

from generate_and_process_clusters import process_clusters


def main():
    p = argparse.ArgumentParser(description="ECHO P4: annotate/filter clusters and write remaining_clusters.parquet")
    p.add_argument("--num_rel", type=int, required=True)
    p.add_argument("--out_remaining", default="remaining_clusters.parquet")
    args = p.parse_args()

    # process_clusters expects output_dir containing clusters*.parquet
    # We run in the task work dir where clusters.parquet is staged, so output_dir='.'
    remaining_df = process_clusters(output_dir=".", num_relatives=args.num_rel)

    # Make remaining clusters explicit for Nextflow
    remaining_df.to_parquet(args.out_remaining, index=False)

    print(f"[OK] wrote {args.out_remaining}")


if __name__ == "__main__":
    main()

