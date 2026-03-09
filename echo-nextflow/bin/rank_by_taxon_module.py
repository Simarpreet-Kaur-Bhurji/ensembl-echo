# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
rank_by_taxon_module.py
-----------------------
Computes pairwise taxonomic distances between query species and all species
present in the input protein set, using the NCBI taxonomy tree (ete3.NCBITaxa).

Key functions:
  - get_all_input_species_combinations() : build (query, input) taxid pairs
  - calculate_taxonomic_distance()       : compute distances and write ranked_taxa.tsv

Note: NCBITaxa is initialised at module import time using ete3's default path.
      Call init_ncbi(dbfile=...) from echo_rank_taxa.py to override the path.
"""

import os

from ete3 import NCBITaxa
import pandas as pd

# NCBITaxa database — initialised with ete3's default path (~/.etetoolkit/taxa.sqlite).
# To use a custom path pass --ncbi_taxa_db to echo_rank_taxa.py, which calls init_ncbi().
# Previous hardcoded path (kept for reference): /homes/sbhurji/.etetoolkit/taxa.sqlite
ncbi = NCBITaxa()


def init_ncbi(dbfile=None):
    """
    (Re-)initialise the module-level NCBITaxa instance.

    Call this from echo_rank_taxa.py before invoking any other function if you
    need a non-default sqlite path (e.g. a shared cluster copy of the DB).

    Args:
        dbfile : absolute path to taxa.sqlite, or None to use ete3's default.
    """
    global ncbi
    ncbi = NCBITaxa(dbfile=dbfile) if dbfile else NCBITaxa()


def get_all_input_species_combinations(query_sps, hcp_taxon_ids):
    """
    Build all (query_taxid, input_taxid) pairs for taxonomic distance calculation.

    Args:
        query_sps     : dict {tax_id: species_name} for the query species
        hcp_taxon_ids : list of tax_ids present in the input protein set

    Returns:
        dict {query_taxid (int): [(query_taxid, input_taxid), ...]}
    """
    query_hcp_combination = {}
    query_taxon_ids = list(query_sps.keys())
    print(f"[get_all_input_species_combinations] query taxon IDs: {query_taxon_ids}")
    print(
        f"[get_all_input_species_combinations] input species count: {len(hcp_taxon_ids)}"
    )

    for i in query_taxon_ids:
        query_hcp_combination[int(i)] = [(int(i), int(qid)) for qid in hcp_taxon_ids]

    return query_hcp_combination


def is_valid_taxid(taxid):
    """Return True if taxid is recognised in the local NCBI taxonomy database."""
    taxid = int(taxid)
    try:
        return taxid in ncbi.get_taxid_translator([taxid])
    except Exception:
        return False


def get_lca_from_lineages(lineage_1, lineage_2):
    """
    Return the deepest (lowest common ancestor) taxid shared by two lineages.

    Lineages are ordered root→leaf; the deepest common node has the largest
    index in either lineage list.
    """
    common = set(lineage_1) & set(lineage_2)
    if not common:
        return None
    return max(
        common, key=lambda taxid: max(lineage_1.index(taxid), lineage_2.index(taxid))
    )


def _ncbi_db_path():
    """Return the path to the NCBITaxa sqlite DB currently in use."""
    return getattr(ncbi, "dbfile", "unknown path")


def _taxid_label(taxid, name_lookup=None):
    """
    Return a human-readable label for a taxid: 'name (taxid)' if a name is
    available, otherwise just the taxid.

    Args:
        taxid       : integer taxid
        name_lookup : optional dict {tax_id (str or int): species_name}
                      (e.g. query_sps from echo_rank_taxa.py)
    """
    if name_lookup:
        name = name_lookup.get(str(taxid)) or name_lookup.get(int(taxid))
        if name:
            return f"{name} (taxid {taxid})"
    # fall back to NCBI DB lookup
    try:
        names = ncbi.get_taxid_translator([taxid])
        if taxid in names:
            return f"{names[taxid]} (taxid {taxid})"
    except Exception:
        pass
    return f"taxid {taxid}"


def calculate_taxonomic_distance(query_hcp_combinations, output_dir, query_sps=None):
    """
    Calculate taxonomic distances between all query–input species pairs.

    Distance = (lineage_depth_1 - common_depth) + (lineage_depth_2 - common_depth),
    i.e. the total number of nodes traversed from each species to their LCA.
    Pairs where distance == 0 (same species) are skipped.

    Args:
        query_hcp_combinations : dict from get_all_input_species_combinations()
        output_dir             : directory to write ranked_taxa.tsv
        query_sps              : optional dict {tax_id: species_name} used to
                                 enrich warning messages with species names

    Returns:
        DataFrame with columns: query_tax_id, input_tid, distance, lca
    """
    rows = []
    skipped_queries = []
    skipped_input_taxids = set()
    db_path = _ncbi_db_path()

    for query_taxid, pairs in query_hcp_combinations.items():
        if not is_valid_taxid(query_taxid):
            label = _taxid_label(query_taxid, query_sps)
            print(
                f"[WARNING] Query species {label} was not found in the NCBI taxonomy database.\n"
                f"  No relatives will be computed for this species.\n"
                f"  DB in use : {db_path}\n"
                f"  Action    : Verify the tax ID at https://www.ncbi.nlm.nih.gov/taxonomy, "
                f" or update the DB with: "
                f' python -c "from ete3 import NCBITaxa; NCBITaxa().update_taxonomy_database()"'
            )
            skipped_queries.append(label)
            continue

        for species_1, species_2 in pairs:
            if not is_valid_taxid(species_2) and species_2 not in skipped_input_taxids:
                label = _taxid_label(species_2, query_sps)
                print(
                    f"[WARNING] Input species {label} was not found in the NCBI taxonomy database.\n"
                    f"  All pairs involving this tax ID will be skipped.\n"
                    f"  DB in use : {db_path}\n"
                    f"  Action    : Verify the tax ID at https://www.ncbi.nlm.nih.gov/taxonomy, "
                    f" or update the DB with: "
                    f' python -c "from ete3 import NCBITaxa; NCBITaxa().update_taxonomy_database()"'
                )
                skipped_input_taxids.add(species_2)

            if not is_valid_taxid(species_2):
                continue

            if species_1 == species_2:
                continue  # distance 0 — query is in the input set

            try:
                lineage_1 = ncbi.get_lineage(species_1)
                lineage_2 = ncbi.get_lineage(species_2)
            except Exception:
                print(
                f" [calculate_taxonomic_distance] error fetching lineage for "
                f" ({species_1}, {species_2}), skipping "
                )
                continue

            if not lineage_1 or not lineage_2:
                print(
                f"[calculate_taxonomic_distance] empty lineage for ({species_1}, {species_2}), skipping"
                )
                continue

            lca = get_lca_from_lineages(lineage_1, lineage_2)
            common_len = len(set(lineage_1) & set(lineage_2))
            total_distance = (len(lineage_1) - common_len) + (
                len(lineage_2) - common_len
            )

            if total_distance > 0:
                rows.append(
                    {
                        "query_tax_id": query_taxid,
                        "input_tid": species_2,
                        "distance": total_distance,
                        "lca": lca,
                    }
                )

    if skipped_queries:
        print(
        f"\n[SUMMARY] {len(skipped_queries)} query species were skipped due to unrecognised tax IDs "
        f"(DB: {db_path}):"
        )
        for label in skipped_queries:
            print(f"  - {label}")

    if skipped_input_taxids:
        print(
        f"\n[SUMMARY] {len(skipped_input_taxids)} input tax IDs were skipped due to unrecognised tax IDs "
        f"(DB: {db_path}):"
        )
        for tid in sorted(skipped_input_taxids):
            print(f"  - {_taxid_label(tid, query_sps)}")

    df = pd.DataFrame(rows)
    os.makedirs(output_dir, exist_ok=True)

    tsv_file = os.path.join(output_dir, "ranked_taxa.tsv")
    df.to_csv(tsv_file, sep="\t", index=False)
    print(f"[calculate_taxonomic_distance] wrote {len(df)} rows to {tsv_file}")
    return df
