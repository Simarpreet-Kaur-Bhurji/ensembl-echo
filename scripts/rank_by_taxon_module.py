"""
rank_by_taxon_module.py

This module provides functionality to calculate taxonomic distances between species 
and generate all possible combinations of target species and HCP taxon IDs.

Functions:
- get_all_input_species_combinations: Generates all combinations of target species and HCP taxon IDs.
- calculate_taxonomic_distance: Calculates taxonomic distances between species pairs using NCBI TaxonHub.

Dependencies:
- NCBITaxonHub: Used to access taxonomic data and calculate distances.
- collections.defaultdict: Used to store taxonomic distances in a nested dictionary.
- pickle, os: Used for file handling and local cache management.
"""

from prototype_material.hub import NCBITaxonHub
from collections import defaultdict
import pickle
import os

def get_all_input_species_combinations(target_sps, hcp_taxon_ids):
    """
    Generate all combinations of target species and HCP taxon IDs.

    Args:
        target_sps (dict): A dictionary where keys are target species taxon IDs and 
                    values are a tuple of speceis name and production name.
        hcp_taxon_ids (list): A list of unique taxon IDs present in high confidence protein file.

    Returns:
        dict: A dictionary where keys are target species taxon IDs, and values are 
              lists of tuples representing all pairs of combinations of target species and HCP taxon IDs.
    """
    target_hcp_combination = dict()
    target_taxon_ids = list(target_sps.keys())
    print("Target Taxon IDs:", target_taxon_ids)
    for i in target_taxon_ids:
        target_hcp_combination[int(i)] = [(int(i), int(qid)) for qid in hcp_taxon_ids]
        
    return target_hcp_combination

local_cache_dir = "ncbi_taxa/archives"
ncbi_taxa_hub = NCBITaxonHub(local_cache_dir)

def calculate_taxonomic_distance(target_hcp_combinations, output_dir):
    """
    Calculate taxonomic distances between species pairs using NCBI TaxonHub.

    Args:
        target_hcp_combinations (dict): A dictionary where keys are target species taxon IDs, 
                                        and values are lists of tuples representing species pairs.
        output_dir (str): Path to the directory where the output will be saved.

    Returns:
        defaultdict: A sorted nested dictionary where the outer keys are target species taxon IDs, 
                     inner keys are the species pair taxon IDs, and values are taxonomic distances.

        Example:
        {
            9606: { #This is the taxon id of target species
                (9606, 10090): (3, 9605), #First tuple is the pair of taxon ids for which distance is being calculated 
                (9606, 10116): (4, 9605), #Second tuple has the first element as taxonomic distance and second element 
                as the last common ancestor taxon id
                ...
            },
           } 
    """
    
    taxa_distance = defaultdict(dict)
    
    for k, v in target_hcp_combinations.items():
        for pair in v:
            species_1 = pair[0]
            species_2 = pair[1]
            if species_1 != species_2:
                with ncbi_taxa_hub.taxon_db('2025-02-01') as ncbi_taxon_db:
                    species_1_lineage = ncbi_taxon_db.lineage_ids(species_1)
                    species_2_lineage = ncbi_taxon_db.lineage_ids(species_2)
                    common_lineage = ncbi_taxon_db.common_lineage_ids([species_1, species_2])
                    lca = ncbi_taxon_db.last_common_taxon_id([species_1, species_2])
                    common_lineage_length = len(common_lineage)
                    distance_species_1 = len(species_1_lineage) - common_lineage_length
                    distance_species_2 = len(species_2_lineage) - common_lineage_length
                    total_distance = distance_species_1 + distance_species_2
                    taxa_distance[k][pair] = (total_distance, lca)

    #print(taxa_distance)           
    # Sort each inner dict by distance
    sorted_taxa = {
        k: dict(sorted(inner.items(), key=lambda item: item[1][0]))
        for k, inner in taxa_distance.items()
    }

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_file = os.path.join(output_dir, "ranked_taxa.json")
    
    with open(output_file, 'wb') as f:
        pickle.dump(sorted_taxa, f)

    print(f"Saved calculated taxonomic distances for every target species to {output_file}")
    return sorted_taxa
