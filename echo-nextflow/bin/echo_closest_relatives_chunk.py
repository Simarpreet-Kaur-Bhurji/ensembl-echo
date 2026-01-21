#!/usr/bin/env python3
import argparse
import duckdb
import pandas as pd


def write_fasta(df: pd.DataFrame, fasta_path: str) -> None:
    """
    df: rows to write (already selected). Writes in current order.
    """
    with open(fasta_path, "w") as fh:
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
    ap.add_argument("--out_log", required=True)
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
    selected = con.execute(f"""
      WITH joined AS (
        SELECT
          c.Cluster_ID,
          c.header,
          c.tax_id,
          c.seq_len,
          COALESCE(c.confidence_level, 'NA') AS confidence_level,
          c.sequence,
          r.distance,

          -- pick best protein per (Cluster_ID, tax_id)
          ROW_NUMBER() OVER (
            PARTITION BY c.Cluster_ID, c.tax_id
            ORDER BY
              r.distance ASC,
              c.seq_len DESC,
              CASE COALESCE(c.confidence_level, 'NA')
                WHEN 'high1' THEN 1
                WHEN 'high2' THEN 2
                WHEN 'high3' THEN 3
                ELSE 99
              END ASC
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
          -- pick top N taxa per cluster; ties on distance broken by seq_len DESC
          ROW_NUMBER() OVER (
            PARTITION BY Cluster_ID
            ORDER BY distance ASC, seq_len DESC
          ) AS rn_cluster
        FROM best_per_taxon
      )
      SELECT
        Cluster_ID, header, tax_id, distance, seq_len, confidence_level, sequence
      FROM topn
      WHERE rn_cluster <= {nrel}
      ORDER BY Cluster_ID, distance ASC, seq_len DESC
    """).fetchdf()

    # always create outputs (even if empty)
    if selected.empty:
        open(args.out_fasta, "w").close()
        pd.DataFrame([]).to_csv(args.out_log, sep="\t", index=False)
        return

    # Write FASTA for this query+chunk
    write_fasta(selected, args.out_fasta)

    # Stats per cluster (same chunk)
    stats = con.execute("""
      SELECT Cluster_ID,
             COUNT(*) AS num_proteins,
             COUNT(DISTINCT CAST(tax_id AS INT)) AS num_unique_tax_ids
      FROM remaining_clusters
      GROUP BY Cluster_ID
    """).fetchdf()

    # Build one log row per Cluster_ID like your original code
    log_rows = []
    for cid, sub in selected.groupby("Cluster_ID", sort=False):
        st = stats[stats["Cluster_ID"] == cid].iloc[0]
        log_rows.append(
            {
                "cluster_id": cid,
                "query_species_name": args.query_name,
                "query_taxonomy_id": str(qtid),
                "closest_relatives": sub["tax_id"].tolist(),
                "closest_proteins": sub["header"].tolist(),
                "closest_distances": sub["distance"].tolist(),
                "num_proteins": int(st["num_proteins"]),
                "num_unique_tax_ids": int(st["num_unique_tax_ids"]),
            }
        )

    pd.DataFrame(log_rows).to_csv(args.out_log, sep="\t", index=False)


if __name__ == "__main__":
    main()

