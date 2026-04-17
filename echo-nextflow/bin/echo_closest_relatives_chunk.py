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
import duckdb
import pandas as pd


def write_fasta(df: pd.DataFrame, fasta_path: str) -> None:
    """
    df: rows to write (already selected). Writes in current order.
    """
    with open(fasta_path, "w", encoding="utf-8") as fh:
        for _, row in df.iterrows():
            header = row["header"]
            seq = row["sequence"]
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Compute closest relatives for one (chunk_parquet, query_tax_id) using one SQL query."
    )
    ap.add_argument("--chunk_parquet", required=True)
    ap.add_argument("--ranked_taxa_tsv", required=True)
    ap.add_argument("--query_tax_id", required=True)
    ap.add_argument("--query_name", required=True)
    ap.add_argument("--num_of_rel", type=int, required=True)
    ap.add_argument("--out_fasta", required=True)
    ap.add_argument("--out_manifest", required=True)  # replaces --out_log; flat per-protein manifest instead of per-cluster list rows
    args = ap.parse_args()

    qtid = int(args.query_tax_id)
    nrel = int(args.num_of_rel)

    con = duckdb.connect()

    # Load inputs
    con.execute(
    f"CREATE OR REPLACE TABLE remaining_clusters AS SELECT * FROM read_parquet('{args.chunk_parquet}')"
    )
    con.execute(
    f"CREATE OR REPLACE TABLE ranked_taxa AS SELECT * FROM read_csv_auto('{args.ranked_taxa_tsv}', sep='\\t')"
    )

    # Efficient: single query returns ALL selected proteins for ALL clusters in this chunk for this query
    selected = con.execute(
        f"""
      WITH joined AS (
        SELECT
          c.Cluster_ID,
          c.header,
          c.tax_id,
          c.seq_len,
          c.sequence,
          r.distance,
          COUNT(DISTINCT c.tax_id) OVER (PARTITION BY c.Cluster_ID) AS cluster_unique_tax_ids,

          -- pick best protein per (Cluster_ID, tax_id)
          -- tiebreak: distance → seq_len → header (alphabetical) for reproducibility
          ROW_NUMBER() OVER (
            PARTITION BY c.Cluster_ID, c.tax_id
            ORDER BY r.distance ASC, c.seq_len DESC, c.header ASC
          ) AS rn_taxon
        FROM remaining_clusters c
        JOIN ranked_taxa r
          ON CAST(c.tax_id AS INT) = CAST(r.input_tid AS INT)
         AND CAST(r.query_tax_id AS INT) = {qtid}
      ),
      best_per_taxon AS (
        SELECT *
        FROM joined
        WHERE rn_taxon = 1
      ),
      topn AS (
        SELECT *,
          -- pick up to N taxa per cluster; sparse clusters get all their taxa
          -- tiebreak: distance → seq_len → tax_id (ascending) for reproducibility
          ROW_NUMBER() OVER (
            PARTITION BY Cluster_ID
            ORDER BY distance ASC, seq_len DESC, tax_id ASC
          ) AS rn_cluster
        FROM best_per_taxon
      )
      SELECT
        Cluster_ID, header, tax_id, distance, seq_len, sequence,
        rn_cluster, cluster_unique_tax_ids  -- exposed so manifest rows carry selection rank and diversity count
      FROM topn
      WHERE rn_cluster <= LEAST(cluster_unique_tax_ids, {nrel})
      ORDER BY Cluster_ID, distance ASC, seq_len DESC, tax_id ASC
    """
    ).fetchdf()

    # always create outputs (even if empty)
    if selected.empty:
        open(args.out_fasta, "w", encoding="utf-8").close()
        # write empty manifest with correct columns so downstream concat always has a header
        pd.DataFrame(columns=[
            "query_tax_id", "query_name", "cluster_id", "protein_header",
            "source_tax_id", "distance", "selection_rank", "cluster_size",
            "unique_tax_ids", "selection_source",
        ]).to_csv(args.out_manifest, sep="\t", index=False)
        return

    # Write FASTA for this query+chunk
    write_fasta(selected, args.out_fasta)

    # Total proteins per cluster (all proteins in the chunk, before selection);
    # used in the manifest so Genebuild can see how large the cluster was
    stats = con.execute(
        """
      SELECT Cluster_ID, COUNT(*) AS cluster_size
      FROM remaining_clusters
      GROUP BY Cluster_ID
    """
    ).fetchdf()
    cluster_size_map = stats.set_index("Cluster_ID")["cluster_size"].to_dict()

    # Flat manifest: one row per selected protein rather than one row per cluster with embedded lists;
    # this makes the output directly queryable by protein_header, source_tax_id, or distance
    manifest_rows = []
    for _, row in selected.iterrows():
        manifest_rows.append(
            {
                "query_tax_id": str(qtid),
                "query_name": args.query_name,
                "cluster_id": row["Cluster_ID"],
                "protein_header": row["header"],
                "source_tax_id": str(row["tax_id"]),
                "distance": row["distance"],
                "selection_rank": int(row["rn_cluster"]),
                "cluster_size": int(cluster_size_map.get(row["Cluster_ID"], 0)),
                "unique_tax_ids": int(row["cluster_unique_tax_ids"]),
                "selection_source": "ranked_cluster",
            }
        )

    pd.DataFrame(manifest_rows).to_csv(args.out_manifest, sep="\t", index=False)


if __name__ == "__main__":
    main()
