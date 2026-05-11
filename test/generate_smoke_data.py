#!/usr/bin/env python3

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
generate_smoke_data.py
Generate synthetic input data for the ECHO end-to-end smoke test.

Outputs (relative to this script's location, i.e. test/):
  data/input_fastas/
    saccharomyces_cerevisiae.fa
    schizosaccharomyces_pombe.fa
    aspergillus_niger.fa
    neurospora_crassa.fa
    candida_albicans.fa
    yarrowia_lipolytica.fa
  data/test_metadata.tsv   — 4-column pipeline metadata (sps_name, taxon_id, gca, file_path)
  data/test_query.tsv      — 2-row query TSV (sps_name, tax_id)

Sequence design
---------------
6 families × 6 species = 36 proteins total (6 per species FASTA).

Family bases: each is a random sequence of a different length (40–62 AA), generated with
a fixed random seed, so sequences across families are unrelated (~5% pairwise identity).

Per-species variants: 5 internal positions mutated (90% identity to base) + 0–8 random
AA appended at the C-terminus.  Pairwise identity between any two variants of the same
family: worst case = base_len - 10 / base_len ≥ 75% → safely above MMseqs2 min_seq_id.
Coverage (cov_mode=1, shorter seq): base is always fully covered by longer variants.

MMseqs2 parameters used (from params.test.yaml):
  min_seq_id=0.75, coverage=0.8, cov_mode=1

Expected clustering result: 7 clusters — 6 multi-member (one per family, 6 proteins each)
+ 1 singleton (aspergillus_niger family 5, asp05) whose sequence is fully random and
shares <10% identity with the family-5 base, so MMseqs2 never groups it.

Deliberate sequence duplicates (for dedup smoke test):
  sch03 ≡ neu03  (schizosaccharomyces_pombe family 3 == neurospora_crassa family 3)
  asp02 ≡ yar02  (aspergillus_niger family 2 == yarrowia_lipolytica family 2)
These appear as relative pairs within separate query outputs so that
params.test.dedup.yaml (dedup_sequences=true) removes exactly 2 sequences.

Tax IDs are real NCBI taxonomy IDs (fungi) so ete3/NCBITaxa can compute distances.

Run:
  python test/generate_smoke_data.py
"""

import csv
import os
import random

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, "data")
FASTA_DIR  = os.path.join(DATA_DIR, "input_fastas")

AA = "ACDEFGHIKLMNPQRSTVWY"  # 20 standard amino acids
N_MUTATIONS = 5   # internal mutations per variant; 90% identity to base, ≥80% pairwise

# Different base length per family; per-species variants extend by 0–8 AA at the C-terminus.
# Coverage is computed on the shorter sequence (cov_mode=1), so even the longest variant
# (base + 8 AA) covers the base at 100% → always above the 0.8 threshold.
FAMILY_LENGTHS = [40, 55, 48, 62, 45, 58]   # base lengths for families 0–5
MAX_EXTENSION  = 8                            # max extra AA appended per species variant

# 6 input species: (sps_name, taxon_id, gca)
# Real NCBI fungal taxonomy IDs so ete3 tree lookups succeed.
SPECIES = [
    ("saccharomyces_cerevisiae",  4932, "GCA_000146045.2"),
    ("schizosaccharomyces_pombe", 4896, "GCA_000002945.2"),
    ("aspergillus_niger",         5061, "GCA_000002655.1"),
    ("neurospora_crassa",         5141, "GCA_000182925.2"),
    ("candida_albicans",          5476, "GCA_000182965.3"),
    ("yarrowia_lipolytica",       4952, "GCA_000002525.1"),
]

# 2 query species (subset of SPECIES above so ete3 distances are meaningful)
QUERIES = [
    (4932, "saccharomyces_cerevisiae"),
    (5476, "candida_albicans"),
]

N_FAMILIES = 6

# Deliberate sequence duplicates: (dst_sps_idx, dst_fam_idx) -> (src_sps_idx, src_fam_idx)
# dst protein gets the identical sequence as src protein, enabling the dedup smoke test.
DUPLICATES = {
    (1, 2): (3, 2),  # sch03 ≡ neu03
    (2, 1): (5, 1),  # asp02 ≡ yar02
}


# ---------------------------------------------------------------------------
# Sequence generation (deterministic via fixed seeds)
# ---------------------------------------------------------------------------

def make_base(family_idx: int) -> str:
    """Random base sequence for one family; length set by FAMILY_LENGTHS[family_idx]."""
    rng = random.Random(42 + family_idx * 997)
    return "".join(rng.choices(AA, k=FAMILY_LENGTHS[family_idx]))


def make_singleton(family_idx: int) -> str:
    """Fully random sequence used for the asp05 singleton — shares ~5% identity with
    the family-5 base by chance, well below min_seq_id=0.75."""
    rng = random.Random(99999 + family_idx * 7)
    return "".join(rng.choices(AA, k=FAMILY_LENGTHS[family_idx]))


def make_variant(base: str, species_idx: int, family_idx: int) -> str:
    """
    Mutate N_MUTATIONS internal positions and append 0–MAX_EXTENSION random AA at the
    C-terminus.  Seed combines species_idx and family_idx so every (species, family)
    pair is unique.  Species 0 (saccharomyces_cerevisiae) receives the unmodified base.
    """
    if species_idx == 0:
        return base
    rng = random.Random(42 + species_idx * 113 + family_idx * 997)
    seq = list(base)
    positions = rng.sample(range(len(base)), N_MUTATIONS)
    for pos in positions:
        alts = [aa for aa in AA if aa != seq[pos]]
        seq[pos] = rng.choice(alts)
    # C-terminal extension: 0 to MAX_EXTENSION extra residues
    extension_len = rng.randint(0, MAX_EXTENSION)
    seq += rng.choices(AA, k=extension_len)
    return "".join(seq)


# ---------------------------------------------------------------------------
# File writers
# ---------------------------------------------------------------------------

def write_fasta(path: str, records: list) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


def write_tsv(path: str, header: list, rows: list) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    os.makedirs(FASTA_DIR, exist_ok=True)

    bases = [make_base(i) for i in range(N_FAMILIES)]

    print(f"Writing to: {DATA_DIR}\n")

    # --- per-species FASTA files ---
    for sps_idx, (sps_name, _tax_id, _gca) in enumerate(SPECIES):
        prefix = sps_name[:3]  # e.g. "sac", "sch", "asp", "neu", "can", "yar"
        records = []
        for fam_idx in range(N_FAMILIES):
            # asp05: aspergillus_niger (sps_idx=2), family 5 (fam_idx=4) is the singleton
            if sps_idx == 2 and fam_idx == 4:
                seq = make_singleton(fam_idx)
            elif (sps_idx, fam_idx) in DUPLICATES:
                src_sps, src_fam = DUPLICATES[(sps_idx, fam_idx)]
                seq = make_variant(bases[src_fam], src_sps, src_fam)
            else:
                seq = make_variant(bases[fam_idx], sps_idx, fam_idx)
            records.append((f"{prefix}{fam_idx + 1:02d}", seq))
        path = os.path.join(FASTA_DIR, f"{sps_name}.fa")
        write_fasta(path, records)
        print(f"  {path}  ({len(records)} proteins)")

    # --- metadata TSV ---
    metadata_path = os.path.join(DATA_DIR, "test_metadata.tsv")
    write_tsv(
        metadata_path,
        ["sps_name", "taxon_id", "gca", "file_path"],
        [(sps, tid, gca, f"{sps}.fa") for sps, tid, gca in SPECIES],
    )
    print(f"\n  {metadata_path}  ({len(SPECIES)} species)")

    # --- query TSV ---
    query_path = os.path.join(DATA_DIR, "test_query.tsv")
    write_tsv(
        query_path,
        ["sps_name", "tax_id"],
        [(sps, tid) for tid, sps in QUERIES],
    )
    print(f"  {query_path}  ({len(QUERIES)} queries)")

    print("\nDone.")
    print("\nExpected clustering outcome: 7 clusters.")
    print("  6 multi-member clusters (one per family, 6 proteins each).")
    print("  1 singleton: asp05 (aspergillus_niger, family 5).")
    print("\nDeliberate sequence duplicates introduced:")
    print("  sch03 ≡ neu03  (schizosaccharomyces_pombe fam3 == neurospora_crassa fam3)")
    print("  asp02 ≡ yar02  (aspergillus_niger fam2 == yarrowia_lipolytica fam2)")
    print("  With dedup_sequences=true, 2 sequences are removed from the final FASTAs.")
    print("\nRun the pipeline with:")
    print("  nextflow run workflows/echo.nf -params-file test/params.test.yaml")


if __name__ == "__main__":
    main()
