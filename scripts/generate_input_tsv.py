import os
import argparse
import csv
import pandas as pd


def write_sequence(out_f, header, seq, combined):
    """
    Writes a sequence to the combined FASTA file with header length.
    Appends the sequence to the combined list.
    """
    seq_len = len(seq)
    header_with_len = f"{header}|{seq_len}"
    combined.append((header_with_len, seq))
    out_f.write(f">{header_with_len}\n")
    for i in range(0, seq_len, 80):
        out_f.write(seq[i : i + 80] + "\n")


def process_fasta_file(filepath, prefix, out_f, combined):
    """
    Reads a single FASTA file, prepends the prefix to headers,
    and writes sequences to the combined file.
    """
    header = None
    seq_lines = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and seq_lines:
                    seq = "".join(seq_lines)
                    write_sequence(out_f, header, seq, combined)
                # Normalize header: take first token before any space
                raw_id = line[1:].strip().split()[0]
                header = f"{raw_id}|{prefix}"
                seq_lines = []
            else:
                seq_lines.append(line)
        # Write last sequence in the file
        if header and seq_lines:
            seq = "".join(seq_lines)
            write_sequence(out_f, header, seq, combined)


def combine_fastas(input_fasta_dir, combined_fasta_file):
    """
    Combine all FASTA files in input_fasta_dir into a single FASTA file.
    Returns a list of tuples (header, sequence).
    """
    combined = []
    with open(combined_fasta_file, "w") as out_f:
        for fname in os.listdir(input_fasta_dir):
            if fname.endswith(".fa"):
                prefix = fname.replace("_canonical_proteins.fa", "")
                filepath = os.path.join(input_fasta_dir, fname)
                process_fasta_file(filepath, prefix, out_f, combined)
    return combined


def load_species_info(tsv_file):
    mapping = {}
    with open(tsv_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames
        for row in reader:
            sci_name = row["Scientific name"].strip()
            norm_name = sci_name.lower().replace(
                " ", "_"
            )  # e.g. "Octopus rubescens" -> "octopus_rubescens"
            mapping[norm_name] = {
                "scientific_name": sci_name,
                "tax_id": row.get("Species taxon_id", "NA"),
                "confidence_score": row.get("confidence score", "NA")
                if "confidence score" in headers
                else "NA",
                "confidence_level": row.get("confidence level", "NA")
                if "confidence level" in headers
                else "NA",
            }
    return mapping


def header_to_species_key(header):
    """
    Convert FASTA header species part to TSV scientific name key for lookup.
    E.g., 'archivesica_marissinica_gca014843695v1' -> 'archivesica_marissinica'
    """
    try:
        species_part = header.split("|")[1]  # second field is species prefix
    except IndexError:
        species_part = "unknown"

    # Split on '_gca' if present, otherwise take the full string
    species_base = (
        species_part.split("_gca")[0] if "_gca" in species_part else species_part
    )

    parts = species_base.split("_")
    return "_".join(parts[:2]).lower()


def write_protein_metadata(combined_fasta, species_map, output_tsv):
    """
    Write protein metadata TSV from combined FASTA sequences.
    - Protein ID from header
    - Scientific name and taxon ID from species_map
    - Sequence length from header
    - Confidence score and level if present
    """
    rows = []

    with open(output_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(
            [
                "protein_id",
                "scientific_name",
                "sequence_length",
                "tax_id",
                "confidence_score",
                "confidence_level",
            ]
        )

        for header, seq in combined_fasta:
            parts = header.split("|")
            if len(parts) < 3:
                raise ValueError(f"Header not in expected format: {header}")

            protein_id = parts[0]
            seq_len_str = parts[-1]
            if not seq_len_str.isdigit():
                seq_len = len(seq)
            else:
                seq_len = int(seq_len_str)

            # Normalize species for lookup
            species_key = header_to_species_key(header)
            info = species_map.get(
                species_key,
                {
                    "scientific_name": "NA",
                    "tax_id": "NA",
                    "confidence_score": "NA",
                    "confidence_level": "NA",
                },
            )

            # Write TSV row
            writer.writerow(
                [
                    protein_id,
                    info["scientific_name"],
                    seq_len,
                    info["tax_id"],
                    info["confidence_score"],
                    info["confidence_level"],
                ]
            )

            # Collect row for Parquet
            rows.append(
                [
                    header,
                    seq,
                    protein_id,
                    info["scientific_name"],
                    info["tax_id"],
                    info["confidence_score"],
                    info["confidence_level"],
                    seq_len,
                ]
            )

    # Write Parquet with sequence included
    parquet_file = os.path.splitext(output_tsv)[0] + ".parquet"
    df = pd.DataFrame(
        rows,
        columns=[
            "header",
            "sequence",
            "protein_id",
            "scientific_name",
            "tax_id",
            "confidence_score",
            "confidence_level",
            "seq_len",
        ],
    )
    df.to_parquet(parquet_file, index=False)
    print(f"Parquet file (with sequences) written to: {parquet_file}")
