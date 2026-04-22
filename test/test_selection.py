"""
test_selection.py
-----------------
Unit tests for echo_closest_relatives_chunk.py selection logic.

Tests the three cases the reviewer specified:
  a) cluster with more than N taxa → exactly N selected
  b) cluster with fewer than N taxa (but > 1) → all taxa retained, ranked by distance
  c) empty chunk (no ranked taxa match) → empty FASTA and valid manifest written

Run with:
  pytest test/test_selection.py -v
"""

import logging
import os
import subprocess
import sys
import random

import duckdb
import pandas as pd
import pytest

log = logging.getLogger("echo.test.selection")


def _test_dir(log_dir: str, test_name: str) -> str:
    """Create and return a named subdirectory inside log_dir for this test's artifacts."""
    path = os.path.join(log_dir, test_name)
    os.makedirs(path, exist_ok=True)
    return path

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "bin", "echo_closest_relatives_chunk.py")
QUERY_TAX_ID = "9999"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_chunk_parquet(path, rows):
    """
    rows: list of dicts with keys Cluster_ID, header, tax_id, sequence, seq_len.
    Writes a parquet file using duckdb.
    """
    df = pd.DataFrame(rows)
    con = duckdb.connect()
    con.register("df", df)
    con.execute(f"COPY df TO '{path}' (FORMAT PARQUET)")


def make_ranked_taxa(path, taxa_distances):
    """
    taxa_distances: list of (input_tid, distance).
    Writes a ranked_taxa.tsv for QUERY_TAX_ID. Always writes the header row
    so duckdb read_csv_auto can infer column names even when there are no data rows.
    """
    rows = [
        {"query_tax_id": int(QUERY_TAX_ID), "input_tid": tid, "distance": dist, "lca": 10}
        for tid, dist in taxa_distances
    ]
    # Always write header — matches real pipeline output; avoids duckdb column inference failure
    pd.DataFrame(rows, columns=["query_tax_id", "input_tid", "distance", "lca"]).to_csv(
        str(path), sep="\t", index=False
    )


def run_chunk_script(work_dir, chunk_parquet, ranked_taxa, num_of_rel):
    out_fasta    = os.path.join(str(work_dir), "out.fa")
    out_manifest = os.path.join(str(work_dir), "out_manifest.tsv")
    result = subprocess.run(
        [
            sys.executable, SCRIPT,
            "--chunk_parquet",   str(chunk_parquet),
            "--ranked_taxa_tsv", str(ranked_taxa),
            "--query_tax_id",    QUERY_TAX_ID,
            "--query_name",      "test_species",
            "--num_of_rel",      str(num_of_rel),
            "--out_fasta",       str(out_fasta),
            "--out_manifest",    str(out_manifest),
        ],
        capture_output=True, text=True
    )
    assert result.returncode == 0, f"Script failed:\num_relatives{result.stderr}"
    return out_fasta, out_manifest


def count_fasta_headers(fasta_path):
    with open(fasta_path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


# ---------------------------------------------------------------------------
# Test a: cluster with more than N taxa → exactly N selected
# ---------------------------------------------------------------------------
def _make_selection_data(work_dir, cluster_id_mode, num_tax_id_per_cluster):
    """
    Create synthetic chunk.parquet and ranked.tsv for selection tests.

    Returns (chunk_parquet_path, ranked_taxa_path, ranked_pairs) where
    ranked_pairs is [(tax_id, distance), ...] sorted by distance ascending.
    """
    assert cluster_id_mode in {"same", "random"}

    tax_ids      = random.sample(range(1, 9001), num_tax_id_per_cluster)
    distances    = random.sample(range(1, 36),   num_tax_id_per_cluster)
    ranked_pairs = sorted(zip(tax_ids, distances), key=lambda x: x[1])
    log.debug("ranked pairs (tax_id, distance): %s", ranked_pairs)

    shared_cluster_id = f"Cluster_{random.randint(1, 10000)}"
    def make_cluster_id():
        return f"Cluster_{random.randint(1, 10000)}" if cluster_id_mode == "random" else shared_cluster_id

    chunk_rows = [
        {"Cluster_ID": make_cluster_id(), "header": f"prot_{tid}", "tax_id": tid,
         "sequence": "ACDEFGHIKL", "seq_len": 10}
        for tid, _ in ranked_pairs
    ]
    log.debug("cluster IDs assigned: %s", [r["Cluster_ID"] for r in chunk_rows])

    chunk_parquet = os.path.join(work_dir, "chunk.parquet")
    ranked_taxa   = os.path.join(work_dir, "ranked.tsv")

    log.debug("writing chunk parquet: %d proteins -> %s", len(chunk_rows), chunk_parquet)
    make_chunk_parquet(chunk_parquet, chunk_rows)
    log.debug("writing ranked taxa tsv (sorted by distance asc) -> %s", ranked_taxa)
    make_ranked_taxa(ranked_taxa, ranked_pairs)

    return chunk_parquet, ranked_taxa, ranked_pairs


def _run_taxa_selection(log_dir, cluster_id_mode, num_relatives, num_tax_id_per_cluster,
                        expected_count, test_label):
    """
    Shared core for taxa selection tests.

    Delegates data creation to _make_selection_data, runs the selection script,
    and asserts that exactly expected_count tax_ids are selected and ascending selection ranks.
    """
    work_dir = _test_dir(log_dir, f"{test_label}_n{num_relatives}_tax_ids{num_tax_id_per_cluster}_clusterid{cluster_id_mode}")

    chunk_parquet, ranked_taxa, ranked_pairs = _make_selection_data(
        work_dir, cluster_id_mode, num_tax_id_per_cluster
    )

    log.debug("running selection script (num_of_rel=%d)", num_relatives)
    out_fasta, out_manifest = run_chunk_script(work_dir, chunk_parquet, ranked_taxa, num_relatives)

    manifest = pd.read_csv(out_manifest, sep="\t")
    log.debug("manifest returned %d rows", len(manifest))
    log.debug("selected tax_ids: %s", sorted(manifest["source_tax_id"].astype(int).tolist()))
    log.debug("selection ranks:  %s", sorted(manifest["selection_rank"].tolist()))
    log.debug("artifacts saved to: %s", work_dir)

    log.info("num_tax_ids=%d  num_relatives=%d  expected_selected=%d",
             num_tax_id_per_cluster, num_relatives, expected_count)

    assert len(manifest) == expected_count, \
        f"Expected {expected_count} rows, got {len(manifest)}"
    assert count_fasta_headers(out_fasta) == expected_count

    selected_tax_ids = set(manifest["source_tax_id"].astype(int))
    expected_tax_ids = {tid for tid, _ in ranked_pairs[:expected_count]}
    assert selected_tax_ids == expected_tax_ids, \
        f"Wrong taxa selected: got {selected_tax_ids}, expected {expected_tax_ids}"

    assert sorted(manifest["selection_rank"].tolist()) == list(range(1, expected_count + 1))

    log.info("PASS: %d proteins selected, tax_ids=%s, ranks=1..%d",
             expected_count, sorted(selected_tax_ids), expected_count)


# cluster with more taxa than N → exactly N selected
@pytest.mark.parametrize("num_relatives,num_tax_id_per_cluster",
    [(n, tids) for n in range(2, 5) for tids in range(n + 1, 7)]
)
def test_more_taxa_than_n_selects_exactly_n(log_dir, cluster_id_mode,
                                            num_relatives, num_tax_id_per_cluster):
    """
    Cluster has more taxa than N — expects exactly N proteins selected (the N closest).
    num_tax_id_per_cluster > num_relatives in all parametrize cases.
    """
    _run_taxa_selection(log_dir, cluster_id_mode,
                        num_relatives, num_tax_id_per_cluster,
                        expected_count=num_relatives,
                        test_label="test_more_taxa")


# cluster with fewer taxa than N → all taxa retained
@pytest.mark.parametrize("num_relatives,num_tax_id_per_cluster",
    [(n, tids) for n in range(3, 7) for tids in range(2, n)]
)
def test_fewer_taxa_than_n_retains_all_parametrized(log_dir, cluster_id_mode,
                                                    num_relatives, num_tax_id_per_cluster):
    """
    Cluster has fewer taxa than N — expects all proteins retained (no artificial cap).
    num_tax_id_per_cluster < num_relatives in all parametrize cases.
    """
    _run_taxa_selection(log_dir, cluster_id_mode,
                        num_relatives, num_tax_id_per_cluster,
                        expected_count=num_tax_id_per_cluster,
                        test_label="test_fewer_taxa")
    
# cluster with equal taxa as N → all taxa retained
@pytest.mark.parametrize("num_relatives,num_tax_id_per_cluster",
    [(n,n) for n in range(2, 7)]
)
def test_equal_taxa_as_n(log_dir, cluster_id_mode,
                         num_relatives, num_tax_id_per_cluster):
    """
    Cluster has equal taxa as N — expects all proteins retained (no artificial cap).
    num_tax_id_per_cluster == num_relatives in all parametrize cases.
    """
    _run_taxa_selection(log_dir, cluster_id_mode,
                        num_relatives, num_tax_id_per_cluster,
                        expected_count=num_tax_id_per_cluster,
                        test_label="test_equal_taxa")


# ---------------------------------------------------------------------------
# Test c: no taxa in ranked_taxa match → empty outputs with valid schema
# ---------------------------------------------------------------------------

def test_no_matching_taxa_produces_empty_outputs(tmp_path):
    """
    Chunk has proteins from taxa not present in ranked_taxa (e.g. query has no
    relatives in this chunk). Expects an empty FASTA and a manifest with correct
    column headers but zero data rows.
    """
    chunk_rows = [
        {"Cluster_ID": "Cluster_1", "header": "prot_999", "tax_id": 999,
         "sequence": "ACDEFGHIKL", "seq_len": 10}
    ]
    chunk_parquet = tmp_path / "chunk.parquet"
    ranked_taxa   = tmp_path / "ranked.tsv"
    make_chunk_parquet(chunk_parquet, chunk_rows)
    # ranked_taxa is empty — no distances for any taxon
    make_ranked_taxa(ranked_taxa, [])
    log.debug("chunk: 1 protein (tax_id=999); ranked_taxa: empty — no matches expected")

    out_fasta, out_manifest = run_chunk_script(tmp_path, chunk_parquet, ranked_taxa, 5)

    n_seqs = count_fasta_headers(out_fasta)
    log.debug("output FASTA sequences: %d (expected 0)", n_seqs)
    assert n_seqs == 0

    manifest = pd.read_csv(out_manifest, sep="\t")
    log.debug("manifest rows: %d (expected 0) | columns: %s", len(manifest), list(manifest.columns))
    assert len(manifest) == 0
    expected_cols = {
        "query_tax_id", "query_name", "cluster_id", "protein_header",
        "source_tax_id", "distance", "selection_rank", "cluster_size",
        "unique_tax_ids", "selection_source",
    }
    missing = expected_cols - set(manifest.columns)
    assert not missing, f"Missing columns: {missing}"
    log.info("PASS: empty chunk produces 0 FASTA seqs and 0 manifest rows with correct schema")


# ---------------------------------------------------------------------------
# Test: best protein per taxon is selected when a taxon has multiple proteins
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("chunk_rows,ranked_pairs,num_relatives,expected_selections", [
    pytest.param(
        # distance decides — two taxa, same seq_len, different distances
        # N=1; tax 301 (dist 5) beats tax 302 (dist 15)
        [
            {"Cluster_ID": "Cluster_1", "header": "prot_a", "tax_id": 301, "sequence": "A" * 15, "seq_len": 15},
            {"Cluster_ID": "Cluster_1", "header": "prot_b", "tax_id": 302, "sequence": "B" * 15, "seq_len": 15},
        ],
        [(301, 5), (302, 15)],
        1,
        {301: "prot_a"},
        id="distance_decides",
    ),
    pytest.param(
        # seq_len decides between taxa — same distance, different seq_len
        # N=1; tax 302 (len 25) beats tax 301 (len 15) since both at dist 5
        [
            {"Cluster_ID": "Cluster_1", "header": "prot_a", "tax_id": 301, "sequence": "A" * 15, "seq_len": 15},
            {"Cluster_ID": "Cluster_1", "header": "prot_b", "tax_id": 302, "sequence": "B" * 25, "seq_len": 25},
        ],
        [(301, 5), (302, 5)],
        1,
        {302: "prot_b"},
        id="seqlen_decides_between_taxa",
    ),
    pytest.param(
        # tax_id decides — same distance, same seq_len across taxa
        # N=1; tax 301 beats tax 302 because 301 < 302 (tax_id ASC tiebreaker)
        [
            {"Cluster_ID": "Cluster_1", "header": "prot_a", "tax_id": 301, "sequence": "A" * 15, "seq_len": 15},
            {"Cluster_ID": "Cluster_1", "header": "prot_b", "tax_id": 302, "sequence": "B" * 15, "seq_len": 15},
        ],
        [(301, 5), (302, 5)],
        2,
        {301: "prot_a", 302: "prot_b"},
        id="taxid_decides_when_distance_and_seqlen_tied",
    ),
    pytest.param(
        # seq_len decides within a taxon — tax 301 has 3 proteins, best (len 15) is prot_a
        # tax 302 also present; both at same distance; N=2 so both taxa are selected
        [
            {"Cluster_ID": "Cluster_1", "header": "prot_a",  "tax_id": 301, "sequence": "A" * 15, "seq_len": 15},
            {"Cluster_ID": "Cluster_1", "header": "prot_a1", "tax_id": 301, "sequence": "A" * 12, "seq_len": 12},
            {"Cluster_ID": "Cluster_1", "header": "prot_b1", "tax_id": 301, "sequence": "B" * 12, "seq_len": 12},
            {"Cluster_ID": "Cluster_1", "header": "prot_b",  "tax_id": 302, "sequence": "B" * 25, "seq_len": 25},
        ],
        [(301, 5), (302, 5)],
        2,
        {301: "prot_a", 302: "prot_b"},
        id="seqlen_decides_within_taxon",
    ),
    pytest.param(
        # first row decides — one taxon, two proteins, distance and seq_len both tied
        # prot_a (inserted first) beats prot_b
        [
            {"Cluster_ID": "Cluster_1", "header": "prot_a", "tax_id": 301, "sequence": "A" * 12, "seq_len": 12},
            {"Cluster_ID": "Cluster_1", "header": "prot_b", "tax_id": 301, "sequence": "B" * 12, "seq_len": 12},
        ],
        [(301, 5)],
        5,
        {301: "prot_a"},
        id="first_row_decides",
    ),
])
def test_best_protein_per_taxon_selected(tmp_path, chunk_rows, ranked_pairs,
                                         num_relatives, expected_selections):
    """
    Best-protein selection follows a three-tier tiebreak:
      1. lowest distance wins
      2. if distance equal → longer sequence wins
      3. if distance and seq_len equal → first row (insertion order) wins
    Each parametrize case isolates exactly one tier.
    """
    chunk_parquet = tmp_path / "chunk.parquet"
    ranked_taxa   = tmp_path / "ranked.tsv"
    make_chunk_parquet(chunk_parquet, chunk_rows)
    make_ranked_taxa(ranked_taxa, ranked_pairs)

    log.debug("input chunk (%d proteins):", len(chunk_rows))
    for r in chunk_rows:
        log.debug("  tax_id=%-5s  seq_len=%-4s  header=%s", r["tax_id"], r["seq_len"], r["header"])
    log.debug("ranked_taxa: %s", ranked_pairs)
    log.debug("num_relatives=%d  expected=%s", num_relatives, expected_selections)

    _, out_manifest = run_chunk_script(tmp_path, chunk_parquet, ranked_taxa, num_relatives)

    manifest = pd.read_csv(out_manifest, sep="\t")
    log.debug("manifest (%d rows):", len(manifest))
    for _, row in manifest.iterrows():
        log.debug("  tax_id=%-5s  header=%-12s  distance=%-4s  rank=%s",
                  row["source_tax_id"], row["protein_header"], row["distance"], row["selection_rank"])

    assert len(manifest) == len(expected_selections), \
        f"Expected {len(expected_selections)} row(s), got {len(manifest)}"

    for tax_id, expected_header in expected_selections.items():
        row = manifest[manifest["source_tax_id"].astype(int) == tax_id]
        actual = row["protein_header"].values[0]
        log.info("tax_id=%s  expected=%s  got=%s  %s",
                 tax_id, expected_header, actual, "PASS" if actual == expected_header else "FAIL")
        assert actual == expected_header, \
            f"tax_id {tax_id}: expected {expected_header!r}, got {actual!r}"
