#!/usr/bin/env python3

import argparse
from generate_input_tsv import combine_fastas, load_species_info, write_protein_metadata


def main():
    p = argparse.ArgumentParser(
        description="ECHO P1: Combine FASTAs and generate processed_input TSV + Parquet"
    )
    p.add_argument("--input_fasta_dir", required=True, help="Directory containing FASTA files")
    p.add_argument("--metadata_tsv", required=True, help="Metadata TSV (Scientific name, Species taxon_id, etc.)")
    p.add_argument("--out_fasta", required=True, help="Output combined FASTA (e.g. combined_input_fasta.fa)")
    p.add_argument("--out_tsv", required=True, help="Output TSV (e.g. processed_input.tsv)")
    p.add_argument("--out_parquet", required=True, help="Output Parquet (e.g. processed_input.parquet)")
    args = p.parse_args()

    # 1) Combine all FASTAs into one FASTA + return list of (header, seq)
    combined = combine_fastas(args.input_fasta_dir, args.out_fasta)

    # 2) Load species metadata mapping
    species_map = load_species_info(args.metadata_tsv)

    # 3) Write TSV + Parquet with explicit output paths
    write_protein_metadata(combined, species_map, args.out_tsv, args.out_parquet)

    print(f"[OK] wrote: {args.out_fasta}")
    print(f"[OK] wrote: {args.out_tsv}")
    print(f"[OK] wrote: {args.out_parquet}")


if __name__ == "__main__":
    main()

