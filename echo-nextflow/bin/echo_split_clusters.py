#!/usr/bin/env python3
import argparse
import os
import duckdb

def main():
    ap = argparse.ArgumentParser(
        description="Split remaining_clusters.parquet into chunk parquet files by Cluster_ID"
    )
    ap.add_argument("--in_parquet", required=True, help="remaining_clusters.parquet")
    ap.add_argument("--out_dir", default="cluster_chunks", help="Output directory")
    ap.add_argument("--clusters_per_chunk", type=int, default=5000, help="How many Cluster_IDs per chunk")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    con = duckdb.connect()

    # Collect unique Cluster_IDs
    cluster_ids = con.execute(
        f"SELECT DISTINCT Cluster_ID FROM read_parquet('{args.in_parquet}') ORDER BY Cluster_ID"
    ).fetchall()
    cluster_ids = [x[0] for x in cluster_ids]

    chunk_size = args.clusters_per_chunk
    if chunk_size <= 0:
        raise ValueError("--clusters_per_chunk must be > 0")

    for idx, start in enumerate(range(0, len(cluster_ids), chunk_size)):
        ids = cluster_ids[start : start + chunk_size]
        out_file = os.path.join(args.out_dir, f"chunk_{idx:03d}.parquet")

        placeholders = ",".join(["?"] * len(ids))
        con.execute(
            f"""
            COPY (
              SELECT *
              FROM read_parquet('{args.in_parquet}')
              WHERE Cluster_ID IN ({placeholders})
            )
            TO '{out_file}' (FORMAT PARQUET)
            """,
            ids,
        )
        print(f"[OK] wrote {out_file} (clusters={len(ids)})")

if __name__ == "__main__":
    main()

