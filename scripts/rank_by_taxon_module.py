from ete3 import NCBITaxa
from collections import defaultdict
import pickle
import os

ncbi = NCBITaxa(dbfile="~/.etetoolkit/taxa.sqlite")
# ncbi.update_taxonomy_database()


def get_all_input_species_combinations(query_sps, hcp_taxon_ids):
    query_hcp_combination = {}
    query_taxon_ids = list(query_sps.keys())
    print("Query Taxon IDs:", query_taxon_ids)

    for i in query_taxon_ids:
        query_hcp_combination[int(i)] = [(int(i), int(qid)) for qid in hcp_taxon_ids]

    return query_hcp_combination


def is_valid_taxid(taxid):
    taxid = int(taxid)
    try:
        valid = taxid in ncbi.get_taxid_translator([taxid])
    except Exception:
        valid = False
    return valid


def get_lca_from_lineages(lineage_1, lineage_2):
    common = set(lineage_1) & set(lineage_2)
    if not common:
        return None
    # the deepest (closest to leaves) is the one with the largest index in lineage
    return max(
        common, key=lambda taxid: max(lineage_1.index(taxid), lineage_2.index(taxid))
    )


def calculate_taxonomic_distance(query_hcp_combinations, output_dir):
    """
    Calculate taxonomic distances between species pairs using ete3.NCBITaxa.

    Args:
        query_hcp_combinations (dict): Keys are query species taxon IDs,
                                        values are lists of tuples (query_taxid, hcp_taxid)
        output_dir (str): Directory where output pickle will be saved.

    Returns:
        defaultdict: Sorted nested dictionary:
        {
            query_taxid: {
                (query_taxid, hcp_taxid): (total_distance, lca_taxid),
                ...
            },
            ...
        }
    """
    taxa_distance = defaultdict(dict)

    for query_taxid, pairs in query_hcp_combinations.items():
        if not is_valid_taxid(query_taxid):
            print(f"Skipping taxid {query_taxid} because it is invalid")
            continue

        for species_1, species_2 in pairs:
            if not is_valid_taxid(species_1) or not is_valid_taxid(species_2):
                print(
                    f"Skipping taxids {species_1}, {species_2} because one is invalid"
                )
                continue

            if species_1 == species_2:
                continue

            try:
                lineage_1 = ncbi.get_lineage(species_1)
                lineage_2 = ncbi.get_lineage(species_2)
            except Exception:
                print(
                    f"Skipping taxids {species_1}, {species_2} due to error fetching lineage"
                )
                continue

            if not lineage_1 or not lineage_2:
                print(
                    f"Skipping taxids {species_1}, {species_2} because lineage is empty"
                )
                continue

            lca = get_lca_from_lineages(lineage_1, lineage_2)
            common_lineage_length = len(set(lineage_1) & set(lineage_2))
            distance_species_1 = len(lineage_1) - common_lineage_length
            distance_species_2 = len(lineage_2) - common_lineage_length
            total_distance = distance_species_1 + distance_species_2

            if total_distance > 0:
                taxa_distance[query_taxid][(species_1, species_2)] = (
                    total_distance,
                    lca,
                )

        # Sort each inner dict by distance
        taxa_distance[query_taxid] = dict(
            sorted(taxa_distance[query_taxid].items(), key=lambda item: item[1][0])
        )

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_file = os.path.join(output_dir, "ranked_taxa.pkl")
    with open(output_file, "wb") as f:
        pickle.dump(taxa_distance, f)

    print(
        f"Saved calculated taxonomic distances for every query species to {output_file}"
    )
    return taxa_distance
