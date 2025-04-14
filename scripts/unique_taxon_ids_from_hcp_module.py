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
>ENSGALG00000008145_gallus_gallus_gca000002315v5_9031, where 9031 is the taxonomy id.
"""
import pysam

def parse_hcp_fasta(input_pep_files):
    """
    Read a FASTA file using pysam and extract sequences.

    Args:
        input_pep_files (str): Path to the FASTA file.

    Returns:
        list: A list of tuples, where each tuple contains:
              - header (str): The full header of the FASTA entry.
              - sequence (str): The protein sequence.
              - protein_id (str): The unique protein identifier extracted from the header.
              - name (str): The name extracted from the header.
              - tax_id (str): The taxonomic ID extracted from the header.

    """
    sequences = []
    fasta = pysam.FastaFile(input_pep_files)
    
    for header in fasta.references:
        sequence = fasta.fetch(header)
        split_header = header.split("_")
        protein_id = split_header[0]
        name = '_'.join(split_header[1:-1])
        tax_id = split_header[-1]
        sequences.append((header, sequence, protein_id, name, tax_id))
    
    fasta.close()
    return sequences

def create_hcp_table(sequences, con):
    """
    Create a DuckDB table from the given sequences and insert the data.

    Args:
        sequences (list): A list of tuples containing sequence data. Each tuple should have:
                          - header (str): The full header of the FASTA entry.
                          - sequence (str): The protein sequence.
                          - protein_id (str): The unique protein identifier.
                          - name (str): The name associated with the protein.
                          - tax_id (str): The taxonomic ID associated with the protein.
        con (duckdb.DuckDBPyConnection): A DuckDB connection object.

    Returns:
        None: The function creates a table named "hcp_table" in the DuckDB connection and inserts the data.
    """
    # Create a table in DuckDB
    con.execute("""
    CREATE TABLE hcp_table (
        header TEXT,
        sequence TEXT,
        protein_id TEXT,
        name TEXT,
        tax_id TEXT
    )
    """)
    # Insert the sequences into the DuckDB table
    con.executemany("INSERT INTO hcp_table VALUES (?, ?, ?, ?, ?)", sequences)
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

