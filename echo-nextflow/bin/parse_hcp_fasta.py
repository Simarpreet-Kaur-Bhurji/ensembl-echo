import pysam
import csv
import os
import pandas as pd


def parse_hcp_fasta(input_pep_files, outdir):
    """
    Parse a FASTA file with pysam and write sequence metadata.

    Args:
        input_pep_files (str): Path to input FASTA file.
        outdir (str): Path to output TSV file (Parquet is also created).

    Returns:
        list: List of tuples with header, sequence, protein_id, name,
              tax_id, confidence_score, confidence_level, and seq_len.
    """
    sequences = []
    fasta = pysam.FastaFile(input_pep_files)

    for header in fasta.references:
        sequence = fasta.fetch(header)
        seq_len = len(sequence)
        split_header = header.split("|")
        protein_id = split_header[0]
        name = split_header[2].replace("_", "")
        tax_id = split_header[4]
        # gene_quality = split_header[-1]
        confidence_level = split_header[-1]
        confidence_score = "NA"
        # sequences.append((header, sequence, protein_id, name, tax_id, gene_quality, seq_len))
        sequences.append(
            (
                header,
                sequence,
                protein_id,
                name,
                tax_id,
                confidence_score,
                confidence_level,
                seq_len,
            )
        )

    fasta.close()

    with open(outdir, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(
            [
                "header",
                "sequence",
                "protein_id",
                "name",
                "tax_id",
                "confidence_score",
                "confidence_level",
                "seq_len",
            ]
        )
        for row in sequences:
            writer.writerow(row)

    print(f"Metadata TSV written to: {outdir}")

    # Parquet
    parquet_file = os.path.splitext(outdir)[0] + ".parquet"
    df = pd.DataFrame(
        sequences,
        columns=[
            "header",
            "sequence",
            "protein_id",
            "name",
            "tax_id",
            "confidence_score",
            "confidence_level",
            "seq_len",
        ],
    )
    df.to_parquet(parquet_file, index=False)
    print(f"Metadata Parquet written to: {parquet_file}")

    return sequences
