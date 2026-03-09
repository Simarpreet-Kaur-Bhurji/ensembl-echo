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
echo_process_clusters.py  –  ECHO pipeline step 4 (P4)

Annotates clusters with cluster_size and unique_tax_ids, writes singleton and
few-taxid side outputs, and emits remaining_clusters.parquet containing only
the clusters eligible for the relatives search (multi-member, enough species).

Called by: modules/process_clusters.nf
Expects:   clusters.parquet staged in the Nextflow work directory (output_dir='.')
Writes to work dir:
  - remaining_clusters.parquet
  - discarded_singletons.fa / singleton_cluster_summary.tsv
  - clusters_with_fewer_tax_ids.fa / clusters_with_fewer_taxids_summary.tsv
"""
import argparse

from generate_and_process_clusters import process_clusters


def main():
    p = argparse.ArgumentParser(
        description="ECHO P4: annotate/filter clusters -> remaining_clusters.parquet"
    )
    p.add_argument(
        "--num_rel",
        type=int,
        required=True,
        help="Minimum unique tax_ids for a cluster to be eligible",
    )
    p.add_argument(
        "--out_remaining",
        default="remaining_clusters.parquet",
        help="Output path for remaining clusters",
    )
    args = p.parse_args()

    print(f"[echo_process_clusters] num_rel:       {args.num_rel}")
    print(f"[echo_process_clusters] out_remaining: {args.out_remaining}")

    # process_clusters reads clusters*.parquet from CWD (Nextflow work dir)
    remaining_df = process_clusters(output_dir=".", num_relatives=args.num_rel)

    # write remaining clusters as an explicit Nextflow output file
    remaining_df.to_parquet(args.out_remaining, index=False)
    print(
    f" [OK] wrote {args.out_remaining}  ({remaining_df['Cluster_ID'].nunique()} clusters,"
    f" {len(remaining_df)} proteins)"
    )


if __name__ == "__main__":
    main()
