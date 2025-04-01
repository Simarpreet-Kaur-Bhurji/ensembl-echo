from prototype_material.hub import NCBITaxonHub
import json

local_cache_dir = "/hps/software/users/ensembl/compara/sbhurji/modenv/ensembl/main/ensembl-echo/scripts/prototype-material/ncbi_taxa/archives"
ncbi_taxa_hub = NCBITaxonHub(local_cache_dir)


def calculate_taxonomic_distance(query_taxon_id, target_taxon_ids):
    """
    Calculate the taxonomic distance between a query taxon and a list of target taxa.
    
    Args:
        query_taxon_id (int): The NCBI taxonomy ID of the query taxon.
        target_taxon_ids (list of int): A list of NCBI taxonomy IDs for target taxa.
    
    Returns:
        dict: A dictionary where keys are target taxonomy IDs and values are tuples of 
              (total taxonomic distance, last common ancestor taxon ID), sorted by distance.
    """
    
    taxa_distance = dict()
    
    with ncbi_taxa_hub.taxon_db('2025-02-01') as ncbi_taxon_db:
        # Get the lineage for the query taxon
        query_lineage = ncbi_taxon_db.lineage_ids(query_taxon_id)
        #print("Query Lineage:", query_lineage)
        
        for target_taxon_id in target_taxon_ids:
            # Get the lineage for the target taxon
            target_lineage = ncbi_taxon_db.lineage_ids(target_taxon_id)
            #print("Target Lineage:", target_lineage)
            
            # Find common lineage IDs and last common ancestor (LCA)
            common_lineage = ncbi_taxon_db.common_lineage_ids([target_taxon_id, query_taxon_id])
            lca = ncbi_taxon_db.last_common_taxon_id([target_taxon_id, query_taxon_id])
            
            # Calculate taxonomic distances
            common_lineage_length = len(common_lineage)
            distance_target = len(target_lineage) - common_lineage_length
            distance_query = len(query_lineage) - common_lineage_length
            total_distance = distance_target + distance_query
            
            #print(f"Total Distance to {target_taxon_id}: {total_distance}, LCA: {lca}")
            
            # Store the result
            taxa_distance[target_taxon_id] = (total_distance, lca)
        
    # Sort taxa by total distance (closer taxa first)
    sorted_taxa = dict(sorted(taxa_distance.items(), key=lambda item: item[1][0]))
    print("Sorted Taxa by Distance:", sorted_taxa)
    
    return sorted_taxa


#print(json.dumps(results1, indent=4))
# Example usage:
if __name__ == "__main__":
    # Example taxon IDs: mouse (10090), chicken (9031), fly (7227) compared to target human (9606)
    target_taxon = [48155, 30419, 57412, 8969]
    query_taxon = [9031]
    ranked_taxa = calculate_taxonomic_distance(query_taxon[0], target_taxon)
    for taxon, (distance, lca) in ranked_taxa.items():
        print(f"Taxon {taxon} is at a distance {distance} from query {query_taxon[0]} (LCA: {lca})")

