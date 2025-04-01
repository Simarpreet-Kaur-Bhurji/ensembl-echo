import pandas as pd
from collections import defaultdict
import re
import requests
import time
import pytaxonkit
import os
import pysam
import duckdb


def process_query_sps(file_path):
    """
    Process the query species file and create a dictionary.

    Args:
        file_path (str): Path to the query species file.

    Returns:
        dict: A dictionary with the taxonomy id as the key and the species name and production name
        as values in a tuple.
    """
    query_sps = {}
    df = pd.read_csv(file_path, delimiter=",")
    print(df)
    for _, row in df.iterrows():
        key = row[0]
        value = (row[1], row[2])
        query_sps[key] = value

    return query_sps
    

def parse_cluster_file(cluster_file):
    """
    Parse the cluster file and create a dictionary of clusters.

    Args:
        cluster_file (str): Path to the cluster file.

    Returns:
        dict: A dictionary where the keys are cluster numbers and the values are lists of cluster members.
    """
    clusters = {}
    current_cluster = None

    with open(cluster_file, "r") as f:
        for line in f:
            line = line.strip()  # Remove leading/trailing whitespaces
            if line.startswith("#Cluster"):
                # Extract cluster number
                cluster_number = line.split()[1]
                current_cluster = cluster_number
                clusters.setdefault(cluster_number, set())
            elif current_cluster and line:
                # Append line to current cluster
                clusters[current_cluster].add(line)

    # Convert sets to lists before returning
    for cluster_number, values_set in clusters.items():
        clusters[cluster_number] = list(values_set)

    return clusters

def get_parents_pytaxon(tax_id):
    """
    Retrieve the lineage of a given taxonomy ID using pytaxonkit.

    Args:
        tax_id (str): The taxonomy ID.

    Returns:
        list: A list of ancestor taxonomy IDs.
    """
    result = pytaxonkit.lineage([tax_id])
    lineage_string = result["FullLineageTaxIDs"].tolist()
    for i in lineage_string:
        val = i.split(";")
        return list(map(int,val))


def read_fasta(input_pep_files):
    """
    Read a FASTA file using pysam and extract sequences.

    Args:
        input_pep_files (str): Path to the FASTA file.

    Returns:
        list: A list of tuples containing header, sequence, protein_id, name, and tax_id.
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

def create_fasta_table(sequences, con):
    """
    Create a DuckDB table from the given sequences and insert the data.

    Args:
        sequences (list): A list of tuples containing header, sequence, protein_id, name, and tax_id.
        con (duckdb.DuckDBPyConnection): The DuckDB connection object.

    Returns:
        None
    """
    # Create a table in DuckDB
    con.execute("""
    CREATE TABLE fasta_table (
        header TEXT,
        sequence TEXT,
        protein_id TEXT,
        name TEXT,
        tax_id TEXT
    )
    """)
    # Insert the sequences into the DuckDB table
    con.executemany("INSERT INTO fasta_table VALUES (?, ?, ?, ?, ?)", sequences)
    # Fetch and print the data to verify
    #df = con.execute("SELECT * FROM fasta_table").fetchdf()
    #print(df)
    return

def get_sequences_from_pep(prot_id, con, outfile):
    """
    Get the sequence for a given protein ID from the DuckDB table and write it to an output file.

    Args:
        prot_id (str): The protein ID.
        con (duckdb.DuckDBPyConnection): The DuckDB connection object.
        outfile (file object): The output file object.

    Returns:
        tuple: The header and sequence that match the input protein ID.
    """
    # Protein from the input csv file
    row = con.execute("SELECT * FROM fasta_table WHERE protein_id = ?", (prot_id,)).fetchdf()
    if not row.empty:
        # Access values in the row
        header = row['header'].values[0]
        sequence = row['sequence'].values[0]
        outfile.write(f">{header}\n{sequence}\n")
        
        return header, sequence
    else:
        print(f"No match found for protein ID: {prot_id}")
        return None

def get_proteins_from_taxid(tax_id, proteins, target, output_path, con):
    """
    Get proteins from a given taxonomy ID and write their sequences to an output file.

    Args:
        tax_id (str): The taxonomy ID.
        proteins (list): List of proteins.
        target (str): The target identifier.
        output_path (str): Path to the output directory.
        con (duckdb.DuckDBPyConnection): The DuckDB connection object.

    Returns:
        None
    """
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    file_name = target[1] + '_relatives.fa'
    file_path = os.path.join(output_path, file_name)
    mode = "a" if os.path.isfile(file_path) else "w"
    
    with open(file_path, mode) as file:
        for j in proteins:
            split_header = j.split("_")
            if str(tax_id) == split_header[-1]:
                print(j)
                prot_id = split_header[0]
                seq = get_sequences_from_pep(prot_id, con, file)

def match_input_and_cluster_proteins(unique_tax_ids, desc, relative, number_of_relatives, proteins, target, tax_id, output, con):
    """
    Match input and cluster proteins based on unique tax IDs and write their sequences to the output file.

    Args:
        unique_tax_ids (list): List of unique taxonomy IDs.
        desc (list): List of descendant taxonomy IDs.
        relative (list): List of relative taxonomy IDs.
        number_of_relatives (int): Number of relatives to match.
        proteins (list): List of proteins.
        target (str): The target identifier.
        tax_id (str): The taxonomy ID.
        output (str): Path to the output directory.
        con (duckdb.DuckDBPyConnection): The DuckDB connection object.

    Returns:
        None
    """
    for i in desc:
        if i in unique_tax_ids and i not in relative and i != tax_id:
            relative.append(i)
            prot_ids = get_proteins_from_taxid(i, proteins, target, output, con)
            if len(relative) == number_of_relatives:
                break

def get_all_desc(tax_id, desc):
    """
    Retrieve all descendant taxa for a given tax ID.

    Args:
        tax_id (str): The taxonomy ID.
        desc (list): List to store descendant taxonomy IDs.

    Returns:
        list: List of descendant taxonomy IDs.
    """
    result = pytaxonkit.list([tax_id])
    for taxon, tree in result:
        subtaxa = [t for t in tree.traverse]
        for j in subtaxa:
            if j[1] == "species" or j[1] == "subspecies":
                desc.append(j[0])
    return desc

def get_lca(tax_ids):
    """
    Retrieve the least common ancestor for a list of tax IDs.

    Args:
        tax_ids (list): List of taxonomy IDs.

    Returns:
        str: The lowest common ancestor.
    """
    lca = pytaxonkit.lca(tax_ids)
    return lca

def get_rel_args(file_path_arg, clusters_arg, input_pep_files):
    """
    Process the query species file, parse the cluster file, retrieve the lowest common ancestor, and create the FASTA table.

    Args:
        file_path_arg (str): Path to the query species file.
        clusters_arg (str): Path to the cluster file.
        input_pep_files (str): Path to the input peptide files.

    Returns:
        tuple: The query species dictionary, clusters dictionary, lowest common ancestor, and DuckDB connection object.
    """
    file_path = file_path_arg
    query_sps = process_query_sps(file_path)
    print("Number of input species", len(query_sps))
    clusters = parse_cluster_file(clusters_arg)
    lca = get_lca(list(query_sps.keys()))
    sequences = read_fasta(input_pep_files)
    con = duckdb.connect()
    input_pep_df = create_fasta_table(sequences, con)
    return query_sps, clusters, lca, con

def get_sequences_for_fasta(proteins, con, output_dir):
    """
    Retrieve sequences for the given proteins and write them to the output file.

    Args:
        proteins (list): List of proteins.
        con (duckdb.DuckDBPyConnection): The DuckDB connection object.
        output_dir (str): Path to the output directory.

    Returns:
        None
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_file = os.path.join(output_dir, "clusters_with_less_proteins.fa")
    
    with open(output_file, "a") as outfile:
        for i in proteins:
            split_header = i.split("_")
            prot_id = split_header[0]
            seq = get_sequences_from_pep(prot_id, con, outfile)

def get_closest_rel_within_cluster(
    number_of_relatives, file_path_arg, clusters_arg, input_pep_files, output
):
    """
    Process clusters to find the closest relatives within each cluster and write their sequences to the output file.

    Args:
        number_of_relatives (int): Number of relatives to match.
        file_path_arg (str): Path to the query species file.
        clusters_arg (str): Path to the cluster file.
        input_pep_files (str): Path to the input peptide files.
        output (str): Path to the output directory.

    Returns:
        None
    """
    query_sps, clusters, lca, con = get_rel_args(
        file_path_arg, clusters_arg, input_pep_files
    )
    for cluster_id, proteins in clusters.items():
        tax_id_of_protein = [int(p.split("_")[-1]) for p in proteins]
        unique_tax_ids = list(dict.fromkeys(tax_id_of_protein))
        print("Cluster No.", cluster_id)
        if len(unique_tax_ids) <= number_of_relatives:
            to_write = get_sequences_for_fasta(proteins, con, output)
        else:
            for tax_id, target in query_sps.items():
                relative = []
                print("tax_id", tax_id, "target", target)
                parents = get_parents_pytaxon(tax_id)
                print(f"Parents of {tax_id}", parents)
                while len(relative) < number_of_relatives:
                    descendants = []
                    # Initialising the list in a loop because this will help reduce the search space.
                    # A great grand parent will have all the child nodes that a grand parent or a parent has.
                    parents.pop()
                    gp = parents[-1]
                    all_desc = get_all_desc(gp, descendants)
                    print(len(descendants))
                    get_gp_rel = match_input_and_cluster_proteins(
                        unique_tax_ids, descendants, relative, number_of_relatives, proteins, target, tax_id, output, con
                    )
                    print(
                        "The tax ids of the genomes that are present in the cluster and in ensembl (gp)",
                        relative,
                    )
                    print("Relatives:",target,":",cluster_id,relative)
                    if gp == lca:
                        print("LCA reached")
                        print("Relatives:",target,":",cluster_id,relative)
                        if len(relative) <  number_of_relatives:
                            print(f"LCA reached and less than {number_of_relatives} found in this cluster")
                            print( "Relatives:", target, ":", cluster_id, relative)
                            break
                        
    return
