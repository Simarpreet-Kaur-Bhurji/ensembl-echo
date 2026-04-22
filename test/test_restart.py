"""
test_restart.py
---------------
Tests for existing_clusters_dir (restart) mode.

Verifies that a new query species can be run against pre-computed clusters
from a previous full run, skipping the expensive clustering steps.

Prerequisites:
  1. Run the standard pipeline first (nextflow_test/ must exist):
       cd repo root
       ./nextflow run main.nf -params-file test/params.test.yaml

  2. Run the restart pipeline:
       cd repo root
       ./nextflow run main.nf -params-file test/params.test.restart.yaml

  3. Run this test:
       pytest test/test_restart.py -v
"""

import logging
import os

import pandas as pd
import pytest

log = logging.getLogger("echo.test.restart")

REPO_ROOT       = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
NF_DIR          = REPO_ROOT
OUTDIR_BASE     = os.path.join(NF_DIR, "nextflow_test")
OUTDIR_RESTART  = os.path.join(NF_DIR, "nextflow_test_restart")

NEW_QUERY       = "schizosaccharomyces_pombe"

# Files that belong to clustering (steps 1-4) and must NOT be re-published
# in the restart outdir (they stay in the base dir).
CLUSTERING_ARTIFACTS = [
    "combined_input_fasta.fa",
    "clusters.parquet",
    "mmseqs_results_cluster.tsv",
    "processed_input.parquet",
    "processed_input.tsv",
]

# Queries from the original run — must NOT appear in the restart outdir.
ORIGINAL_QUERIES = [
    "saccharomyces_cerevisiae",
    "candida_albicans",
]


@pytest.fixture(scope="module")
def base_outdir():
    if not os.path.isdir(OUTDIR_BASE):
        pytest.skip(
            "nextflow_test/ not found — run the standard pipeline first:\n"
            "  ./nextflow run main.nf -params-file test/params.test.yaml"
        )
    log.info("base outdir: %s", OUTDIR_BASE)
    return OUTDIR_BASE


@pytest.fixture(scope="module")
def restart_outdir(base_outdir):
    if not os.path.isdir(OUTDIR_RESTART):
        pytest.skip(
            "nextflow_test_restart/ not found — run the restart pipeline first:\n"
            "  ./nextflow run main.nf -params-file test/params.test.restart.yaml"
        )
    log.info("restart outdir: %s", OUTDIR_RESTART)
    return OUTDIR_RESTART


# ---------------------------------------------------------------------------
# New query outputs exist
# ---------------------------------------------------------------------------

def test_new_query_all_relatives_fasta_exists(restart_outdir):
    path = os.path.join(restart_outdir, f"{NEW_QUERY}_all_relatives.fa")
    log.debug("checking FASTA: %s", path)
    assert os.path.isfile(path), f"Missing: {path}"
    log.info("PASS: %s_all_relatives.fa found", NEW_QUERY)


def test_new_query_all_relatives_fasta_nonempty(restart_outdir):
    path = os.path.join(restart_outdir, f"{NEW_QUERY}_all_relatives.fa")
    with open(path, encoding="utf-8") as fh:
        seqs = [l for l in fh if l.startswith(">")]
    log.debug("%s_all_relatives.fa: %d sequences", NEW_QUERY, len(seqs))
    assert len(seqs) > 0, f"{NEW_QUERY}_all_relatives.fa is empty"
    log.info("PASS: %s_all_relatives.fa has %d sequences", NEW_QUERY, len(seqs))


def test_new_query_manifest_exists(restart_outdir):
    path = os.path.join(restart_outdir, f"{NEW_QUERY}_manifest.tsv")
    log.debug("checking manifest: %s", path)
    assert os.path.isfile(path), f"Missing: {path}"
    log.info("PASS: %s_manifest.tsv found", NEW_QUERY)


def test_new_query_manifest_columns(restart_outdir):
    path = os.path.join(restart_outdir, f"{NEW_QUERY}_manifest.tsv")
    df = pd.read_csv(path, sep="\t")
    expected = {
        "query_tax_id", "query_name", "cluster_id", "protein_header",
        "source_tax_id", "distance", "selection_rank", "cluster_size",
        "unique_tax_ids", "selection_source",
    }
    missing = expected - set(df.columns)
    log.debug("manifest columns: %s", sorted(df.columns.tolist()))
    log.debug("expected columns: %s", sorted(expected))
    assert not missing, f"Missing columns: {missing}"
    log.info("PASS: all %d expected manifest columns present", len(expected))


def test_new_query_manifest_row_count_matches_fasta(restart_outdir):
    """Every sequence in the FASTA has a row in the manifest."""
    fa_path  = os.path.join(restart_outdir, f"{NEW_QUERY}_all_relatives.fa")
    tsv_path = os.path.join(restart_outdir, f"{NEW_QUERY}_manifest.tsv")
    with open(fa_path, encoding="utf-8") as fh:
        n_seqs = sum(1 for l in fh if l.startswith(">"))
    n_rows   = len(pd.read_csv(tsv_path, sep="\t"))
    log.debug("FASTA sequences: %d  manifest rows: %d", n_seqs, n_rows)
    assert n_rows == n_seqs, f"Manifest rows ({n_rows}) != FASTA sequences ({n_seqs})"
    log.info("PASS: manifest row count (%d) matches FASTA sequence count", n_rows)


# ---------------------------------------------------------------------------
# Clustering was NOT repeated — heavy artifacts absent from restart outdir
# ---------------------------------------------------------------------------

def test_clustering_artifacts_not_republished(restart_outdir):
    """Clustering outputs must not be re-published; they stay in existing_clusters_dir."""
    for artifact in CLUSTERING_ARTIFACTS:
        path = os.path.join(restart_outdir, artifact)
        exists = os.path.isfile(path)
        log.debug("clustering artifact '%s' in restart dir: %s (expected: absent)", artifact, exists)
        assert not exists, \
            f"{artifact} should not be in restart outdir (it lives in the base dir)"
    log.info("PASS: none of %d clustering artifacts re-published in restart dir",
             len(CLUSTERING_ARTIFACTS))


# ---------------------------------------------------------------------------
# Original queries are not processed again
# ---------------------------------------------------------------------------

def test_original_query_outputs_absent(restart_outdir):
    """The restart run only processes the new query; original query FASTAs are not re-created."""
    for query in ORIGINAL_QUERIES:
        path = os.path.join(restart_outdir, f"{query}_all_relatives.fa")
        exists = os.path.isfile(path)
        log.debug("original query '%s' FASTA in restart dir: %s (expected: absent)", query, exists)
        assert not exists, \
            f"{query}_all_relatives.fa should not appear in restart outdir"
    log.info("PASS: original query outputs absent from restart dir (%s)",
             ", ".join(ORIGINAL_QUERIES))


# ---------------------------------------------------------------------------
# Clustering artifacts are still present in the base dir (unchanged)
# ---------------------------------------------------------------------------

def test_base_dir_clustering_artifacts_intact(base_outdir):
    """The base run's clustering outputs are untouched by the restart."""
    for artifact in CLUSTERING_ARTIFACTS:
        path = os.path.join(base_outdir, artifact)
        exists = os.path.isfile(path)
        log.debug("base dir artifact '%s': %s (expected: present)", artifact, exists)
        assert exists, \
            f"{artifact} missing from base outdir — base run may be incomplete"
    log.info("PASS: all %d clustering artifacts intact in base dir", len(CLUSTERING_ARTIFACTS))


# ---------------------------------------------------------------------------
# Diagnostics produced for the new query
# ---------------------------------------------------------------------------

DIAGNOSTICS_FILES = {
    "diagnostics.pdf",
    "diagnostics_summary.txt",
    "plot1_cluster_size_distribution.png",
    "plot2_unique_taxids_distribution.png",
    "plot3_cluster_retention_summary.png",
    "plot4_taxonomic_distance_distribution.png",
    "plot5_size_vs_taxid_scatter.png",
}


def test_diagnostics_files_present(restart_outdir):
    diag_dir = os.path.join(restart_outdir, "diagnostics_out")
    log.debug("checking diagnostics dir: %s", diag_dir)
    assert os.path.isdir(diag_dir), "diagnostics_out/ directory missing"
    for fname in DIAGNOSTICS_FILES:
        path = os.path.join(diag_dir, fname)
        exists = os.path.isfile(path)
        log.debug("diagnostics file '%s': %s", fname, "present" if exists else "MISSING")
        assert exists, f"Missing diagnostics file: {fname}"
    log.info("PASS: all %d diagnostics files present in restart outdir", len(DIAGNOSTICS_FILES))
