#!/usr/bin/env python3
import argparse
import duckdb
import os

from process_input_species_module import get_input_sps
from rank_by_taxon_module import (
    get_all_input_species_combinations,
    calculate_taxonomic_distance,
)

def main():
    ap = argparse.ArgumentParser(
        description="Compute ranked_taxa.tsv from query_species and processed_input.parquet"
    )
    ap.add_argument("--query_species", required=True, help="TSV with tax_id and sps_name columns")
    ap.add_argument("--processed_input_parquet", required=True, help="processed_input.parquet")
    ap.add_argument("--out_tsv", default="ranked_taxa.tsv", help="Output ranked taxa TSV")
    args = ap.parse_args()

    # query species dict: {tax_id: species_name}
    query_sps = get_input_sps(args.query_species)

    # get all distinct tax_ids present in processed_input
    con = duckdb.connect()
    taxon_ids = con.execute(
        f"SELECT DISTINCT tax_id FROM read_parquet('{args.processed_input_parquet}')"
    ).fetchdf()["tax_id"].tolist()

    # compute combinations + distances (your function writes ranked_taxa.tsv to output_dir)
    combos = get_all_input_species_combinations(query_sps, taxon_ids)
    calculate_taxonomic_distance(combos, output_dir=".")

    # standardize output filename if needed
    if args.out_tsv != "ranked_taxa.tsv":
        if os.path.exists(args.out_tsv):
            os.remove(args.out_tsv)
        os.rename("ranked_taxa.tsv", args.out_tsv)

    print(f"[OK] wrote {args.out_tsv}")

if __name__ == "__main__":
    main()

