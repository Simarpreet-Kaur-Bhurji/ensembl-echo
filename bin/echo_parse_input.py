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
from generate_input_tsv import (
    combine_fastas_from_map,
    load_species_info,
    write_protein_metadata,
)


def main():
    p = argparse.ArgumentParser(
        description="ECHO P1: Combine FASTAs and generate processed_input TSV + Parquet"
    )
    p.add_argument(
        "--input_fasta_dir",
        required=False,
        default=None,
        help=(
            "Optional base directory for resolving relative file_path entries in the "
            "metadata TSV. If file_path in the TSV is absolute this argument is ignored. "
            "If omitted, relative file_path values are resolved against the current "
            "working directory."
        ),
    )
    p.add_argument(
        "--metadata_tsv",
        required=True,
        help="Metadata TSV with columns: sps_name, taxon_id, gca, file_path",
    )
    p.add_argument(
        "--out_fasta",
        required=True,
        help="Output combined FASTA (e.g. combined_input_fasta.fa)",
    )
    p.add_argument(
        "--out_tsv", required=True, help="Output TSV (e.g. processed_input.tsv)"
    )
    p.add_argument(
        "--out_parquet",
        required=True,
        help="Output Parquet (e.g. processed_input.parquet)",
    )
    p.add_argument(
        "--fasta_line_width",
        type=int,
        default=60,
        help="Characters per sequence line in the output combined FASTA (default: 60)",
    )
    args = p.parse_args()

    print(f"[echo_parse_input] metadata_tsv:      {args.metadata_tsv}")
    print(f"[echo_parse_input] input_fasta_dir:   {args.input_fasta_dir!r}")
    print(f"[echo_parse_input] out_fasta:         {args.out_fasta}")
    print(f"[echo_parse_input] out_tsv:           {args.out_tsv}")
    print(f"[echo_parse_input] out_parquet:       {args.out_parquet}")
    print(f"[echo_parse_input] fasta_line_width:  {args.fasta_line_width}")

    # 1) Load species metadata (includes file_path per species)
    species_map = load_species_info(args.metadata_tsv)

    # 2) Combine FASTAs using file_path from the TSV; resolve relative paths against input_fasta_dir
    combined = combine_fastas_from_map(
        species_map, args.out_fasta, args.input_fasta_dir, args.fasta_line_width
    )

    # 3) Write processed TSV + Parquet
    write_protein_metadata(combined, species_map, args.out_tsv, args.out_parquet)

    print(f"[OK] wrote: {args.out_fasta}")
    print(f"[OK] wrote: {args.out_tsv}")
    print(f"[OK] wrote: {args.out_parquet}")


if __name__ == "__main__":
    main()
