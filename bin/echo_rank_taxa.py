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

import duckdb

from process_input_species_module import get_input_sps
from rank_by_taxon_module import (
    init_ncbi,
    get_all_input_species_combinations,
    calculate_taxonomic_distance,
)


def main():
    ap = argparse.ArgumentParser(
        description="Compute ranked_taxa.tsv from query_species and processed_input.parquet"
    )
    ap.add_argument(
        "--query_species", required=True, help="TSV with tax_id and sps_name columns"
    )
    ap.add_argument(
        "--processed_input_parquet", required=True, help="processed_input.parquet"
    )
    ap.add_argument(
        "--out_tsv", default="ranked_taxa.tsv", help="Output ranked taxa TSV"
    )
    ap.add_argument(
        "--ncbi_taxa_db",
        default=None,
        help=(
            "Path to NCBITaxa sqlite database. "
            "Defaults to ete3's own default (~/.etetoolkit/taxa.sqlite). "
            "Set params.ncbi_taxa_db in your Nextflow params YAML to use a shared copy."
        ),
    )
    args = ap.parse_args()

    # Override NCBITaxa db path before any taxonomy lookups
    if args.ncbi_taxa_db:
        print(f"[echo_rank_taxa] using custom NCBITaxa db: {args.ncbi_taxa_db}")
        init_ncbi(dbfile=args.ncbi_taxa_db)

    # query species dict: {tax_id: species_name}
    query_sps = get_input_sps(args.query_species)

    # get all distinct tax_ids present in processed_input
    con = duckdb.connect()
    taxon_ids = (
        con.execute(
            f"SELECT DISTINCT tax_id FROM read_parquet('{args.processed_input_parquet}')"
        )
        .fetchdf()["tax_id"]
        .tolist()
    )

    # compute combinations + distances (your function writes ranked_taxa.tsv to output_dir)
    combos = get_all_input_species_combinations(query_sps, taxon_ids)
    calculate_taxonomic_distance(combos, output_dir=".", query_sps=query_sps)

    # standardize output filename if needed
    if args.out_tsv != "ranked_taxa.tsv":
        if os.path.exists(args.out_tsv):
            os.remove(args.out_tsv)
        os.rename("ranked_taxa.tsv", args.out_tsv)

    print(f"[OK] wrote {args.out_tsv}")


if __name__ == "__main__":
    main()
