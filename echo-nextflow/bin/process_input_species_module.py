"""
process_input_species_module.py

This script provides functionality to process the input species file and extract relevant data.

Functions:
- get_input_sps: Reads a CSV file containing species information and converts it into a dictionary.

Dependencies:
- pandas: Used for reading and processing the CSV file.

Input File Format:
The input species file should be a CSV file with the following columns:
1. Taxonomy ID (unique identifier for each species)
2. Species Name
3. Production Name
"""
import pandas as pd


def get_input_sps(file_path):
    """
    Process the input species file and create a dictionary.

    Args:
        file_path (str): Path to the query species file.

    Returns:
        dict: A dictionary with the taxonomy id as the key and the species name
        as values.
    """
    all_species = {}
    df = pd.read_csv(file_path, delimiter="\t", header = 0)
    print(df)
    for _, row in df.iterrows():
        key = row["tax_id"]
        value = row["sps_name"]
        all_species[key] = value

    return all_species
