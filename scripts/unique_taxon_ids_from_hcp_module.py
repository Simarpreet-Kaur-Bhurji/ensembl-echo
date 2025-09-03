"""
unique_taxon_ids_from_hcp_module.py

This module provides functionality to parse a FASTA file containing high-confidence protein (HCP) data 
and extract unique taxonomic IDs. It also allows creating a DuckDB table for efficient querying of the data.

Functions:
- parse_hcp_fasta: Reads a FASTA file and extracts relevant information such as header, sequence, 
  protein ID, name, and taxonomic ID.
- create_hcp_table: Creates a DuckDB table from the extracted sequences for efficient data storage and querying.

Dependencies:
- pysam: Used for reading and parsing the FASTA file.

Input File Format:
The headers of the input Fasta file should have taxon id appended to their name like the following example:
><transcript_stable_id>|<gene_stable_id>|<species_name>|<GCA>|<taxon_id>|<gene_quality>
"""

import pandas as pd

def create_hcp_table(metadata_pq_file, con):
    # Create a table in DuckDB

    # con.execute("""
    # CREATE TABLE hcp_table (
    #     header TEXT,
    #     sequence TEXT,
    #     protein_id TEXT,
    #     name TEXT,
    #     tax_id TEXT,
    #     confidence_score TEXT,
    #     confidence_level TEXT,
    #     seq_len INTEGER
    # )
    # """)
    con.execute(f"""
    CREATE TABLE hcp_table AS SELECT * FROM '{metadata_pq_file}'
    """)
    
    #print(f"DuckDB table 'hcp_table' created from Parquet file: {metadata_pq_file}")
    
    
    # Fetch and print the data to verify
    #df = con.execute("SELECT * FROM fasta_table").fetchdf()
    #print(df)
    return

def get_hcp_tax_ids(con):
    """
    Get the unique tax IDs from the DuckDB table.

    Args:
        con (duckdb.DuckDBPyConnection): The DuckDB connection object.

    Returns:
        list: A list of unique tax IDs.
    """
    # Fetch unique tax IDs from the DuckDB table
    tax_ids = con.execute("SELECT DISTINCT tax_id FROM hcp_table").fetchdf()
    return tax_ids['tax_id'].tolist()




# def parse_hcp_fasta(input_pep_files, outdir):
#     """
#     Read a FASTA file using pysam and extract sequences.

#     Args:
#         input_pep_files (str): Path to the FASTA file.

#     Returns:
#         list: A list of tuples, where each tuple contains:
#               - header (str): The full header of the FASTA entry.
#               - sequence (str): The protein sequence.
#               - protein_id (str): The unique protein identifier extracted from the header.
#               - name (str): The name extracted from the header.
#               - tax_id (str): The taxonomic ID extracted from the header.
#               - gene_quality (str): The gene quality extracted from the header.

#     """
#     sequences = []
#     fasta = pysam.FastaFile(input_pep_files)
    
#     for header in fasta.references:
#         sequence = fasta.fetch(header)
#         seq_len = len(sequence)
#         split_header = header.split("|")
#         protein_id = split_header[0]
#         name = split_header[2]
#         tax_id = split_header[4]
#         #gene_quality = split_header[-1]
#         confidence_level = split_header[-1]
#         confidence_score = "NA"
#         #sequences.append((header, sequence, protein_id, name, tax_id, gene_quality, seq_len))
#         sequences.append((header, sequence, protein_id, name, tax_id, confidence_score, confidence_level, seq_len))
    
#     fasta.close()

#     with open(outdir, "w", newline="") as out:
#             writer = csv.writer(out, delimiter="\t")
#             writer.writerow([
#                 "header", "sequence", "protein_id",
#                 "name", "tax_id", "confidence_score",
#                 "confidence_level", "seq_len"
#             ])
#             for row in sequences:
#                 writer.writerow(row)
        
#     print(f"Metadata TSV written to: {outdir}")

#         # Parquet
#     parquet_file = os.path.splitext(outdir)[0] + ".parquet"
#     df = pd.DataFrame(sequences, columns=[
#             "header", "sequence", "protein_id",
#             "name", "tax_id", "confidence_score",
#             "confidence_level", "seq_len"
#         ])
#     df.to_parquet(parquet_file, index=False)
#     print(f"Metadata Parquet written to: {parquet_file}")

#     return sequences

