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

"""
echo_parse_clusters.py  –  ECHO pipeline step 3 (P3)

Parses the MMseqs2 cluster TSV and joins it with processed_input.parquet to
produce clusters.parquet: one row per (Cluster_ID, protein) enriched with
tax_id, sequence, name, and confidence fields.

Called by: modules/parse_clusters.nf
"""
import argparse
import os

import duckdb

from generate_and_process_clusters import parse_cluster_file


def main():
    p = argparse.ArgumentParser(
        description="ECHO P3: parse mmseqs cluster TSV -> clusters.parquet"
    )
    p.add_argument(
        "--cluster_tsv", required=True, help="MMseqs2 two-column cluster TSV"
    )
    p.add_argument(
        "--processed_input_parquet",
        required=True,
        help="processed_input.parquet from P1",
    )
    p.add_argument(
        "--out_parquet", default="clusters.parquet", help="Output parquet path"
    )
    args = p.parse_args()

    print(f"[echo_parse_clusters] cluster_tsv:             {args.cluster_tsv}")
    print(
        f"[echo_parse_clusters] processed_input_parquet: {args.processed_input_parquet}"
    )
    print(f"[echo_parse_clusters] out_parquet:             {args.out_parquet}")

    con = duckdb.connect()

    # parse_cluster_file always writes ./clusters.parquet; rename afterwards if needed
    parse_cluster_file(
        raw_tsv=args.cluster_tsv,
        sequences_parquet=args.processed_input_parquet,
        con=con,
        output_dir=".",
    )

    # rename if the caller requested a non-default output path
    if args.out_parquet != "clusters.parquet":
        os.rename("clusters.parquet", args.out_parquet)

    print(f"[OK] wrote {args.out_parquet}")


if __name__ == "__main__":
    main()
