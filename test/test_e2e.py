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
test_e2e.py — End-to-end smoke test for the ECHO Nextflow pipeline.

Validates all four required output types against known properties of the
synthetic test dataset (6 species, 36 proteins, 6 families):

  1. Per-query FASTA       (*_all_relatives.fa)
  2. Per-query manifest    (*_manifest.tsv)
  3. Summary reports       (cluster_summary.txt, echo_pipeline_summary.txt)
  4. Diagnostics bundle    (diagnostics_out/)

Synthetic dataset properties (deterministic — fixed-seed sequences):
  - 36 input proteins across 6 species
  - 7 clusters after MMseqs2: 6 multi-member + 1 singleton
  - 1 singleton protein (aspergillus_niger, Cluster_5)
  - 2 query species: saccharomyces_cerevisiae (4932), candida_albicans (5476)
  - num_of_rel = 5 (from params.test.yaml)
  - Per query: 29 relatives selected (query's own protein excluded per cluster),
               + 1 singleton appended → 30 proteins in *_all_relatives.fa

Generate test data:
  python test/generate_smoke_data.py

Run pipeline (once, or whenever test data changes):
  cd repo root
  nextflow run workflows/echo.nf -params-file test/params.test.yaml

Run this test against existing output (default nextflow_test/):
  pytest test/test_e2e.py -v

Run this test against a specific output directory:
  pytest test/test_e2e.py -v --outdir /path/to/your/run

Run this test AND execute the pipeline:
  pytest test/test_e2e.py -v --run-pipeline
"""

import logging
import os
import subprocess
import duckdb
import pandas as pd
import pytest
import yaml

log = logging.getLogger("echo.test.e2e")


def _outdir_from_params(params_file: str) -> str:
    """Read the outdir value from a params YAML file, resolved relative to NF_DIR."""
    with open(params_file, encoding="utf-8") as fh:
        params = yaml.safe_load(fh)
    outdir = params.get("outdir")
    if not outdir:
        raise ValueError(f"'outdir' not set in {params_file}")
    return outdir if os.path.isabs(outdir) else os.path.join(NF_DIR, outdir)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO_ROOT              = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
NF_DIR                 = REPO_ROOT
PARAMS_FILE            = os.path.join(NF_DIR, "test", "params.test.yaml")
PARAMS_FILE_NO_SINGLETONS = os.path.join(NF_DIR, "test", "params.test.no_singletons.yaml")
OUTDIR                 = os.path.join(NF_DIR, "nextflow_test")
OUTDIR_NO_SINGLETONS   = os.path.join(NF_DIR, "nextflow_test_no_singletons")

MANIFEST_COLS = {
    "query_tax_id", "query_name", "cluster_id", "protein_header",
    "source_tax_id", "distance", "selection_rank", "cluster_size",
    "unique_tax_ids", "selection_source",
}

DIAGNOSTICS_FILES = {
    "diagnostics.pdf",
    "diagnostics_summary.txt",
    "plot1_cluster_size_distribution.png",
    "plot2_unique_taxids_distribution.png",
    "plot3_cluster_retention_summary.png",
    "plot4_taxonomic_distance_distribution.png",
    "plot5_size_vs_taxid_scatter.png",
}

# Properties of the synthetic dataset (from generate_smoke_data.py)
TOTAL_INPUT_PROTEINS  = 36
TOTAL_CLUSTERS        = 7
SINGLETON_COUNT       = 1
REMAINING_CLUSTERS    = 6         # multi-member clusters eligible for selection
QUERY_SPECIES = {
    "saccharomyces_cerevisiae": 4932,
    "candida_albicans":          5476,
}
NUM_OF_REL            = 5         # matches params.test.yaml
EXPECTED_RELATIVES    = 30        # per query: 29 ranked_cluster rows + 1 singleton row
EXPECTED_ALL_RELATIVES = 30       # relatives + 1 singleton


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def outdir(request):
    """
    Returns the pipeline output directory.
    If --outdir is passed, uses that path directly (no pipeline execution).
    If --run-pipeline is passed, executes the pipeline first and fails fast on error.
    Falls back to the default nextflow_test/ location.
    """
    custom_outdir = request.config.getoption("--outdir")
    if custom_outdir:
        resolved = os.path.abspath(custom_outdir)
        if not os.path.isdir(resolved):
            pytest.skip(f"--outdir path not found: {resolved}")
        log.info("Using custom output directory: %s", resolved)
        return resolved

    if request.config.getoption("--run-pipeline"):
        log.info("Running Nextflow pipeline (this may take several minutes)...")
        result = subprocess.run(
            ["nextflow", "run", "workflows/echo.nf", "-params-file", PARAMS_FILE],
            cwd=NF_DIR,
            capture_output=False,   # stream output so progress is visible
            text=True,
            check=False,
        )
        assert result.returncode == 0, "Nextflow pipeline failed — check output above"
        log.info("Pipeline complete.")
        resolved = _outdir_from_params(PARAMS_FILE)
        log.info("Output directory (from params): %s", resolved)
        return resolved
    log.info("Skipping pipeline execution (pass --run-pipeline to run it).")

    if not os.path.isdir(OUTDIR):
        pytest.skip(
            f"Output directory not found: {OUTDIR}\n"
            "Run the pipeline first:  cd repo root && "
            "nextflow run workflows/echo.nf -params-file test/params.test.yaml"
        )
    return OUTDIR


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _count_fasta_headers(path: str) -> int:
    with open(path, encoding="utf-8") as fh:
        return sum(1 for line in fh if line.startswith(">"))


# ---------------------------------------------------------------------------
# 1. Clustering outputs
# ---------------------------------------------------------------------------

class TestClustering:
    """Validates clusters.parquet and remaining_clusters.parquet."""

    def test_total_clusters(self, outdir):
        con = duckdb.connect()
        n = con.execute(
            f"SELECT COUNT(DISTINCT Cluster_ID) FROM read_parquet('{outdir}/clusters.parquet')"
        ).fetchone()[0]
        log.info("total clusters: %d", n)
        assert n == TOTAL_CLUSTERS, f"Expected {TOTAL_CLUSTERS} clusters, got {n}"

    def test_total_proteins(self, outdir):
        con = duckdb.connect()
        n = con.execute(
            f"SELECT COUNT(*) FROM read_parquet('{outdir}/clusters.parquet')"
        ).fetchone()[0]
        log.info("total proteins in clusters: %d", n)
        assert n == TOTAL_INPUT_PROTEINS

    def test_remaining_clusters_count(self, outdir):
        con = duckdb.connect()
        n = con.execute(
            f"SELECT COUNT(DISTINCT Cluster_ID) FROM read_parquet('{outdir}/remaining_clusters.parquet')"
        ).fetchone()[0]
        log.info("remaining (multi-member) clusters: %d", n)
        assert n == REMAINING_CLUSTERS

    def test_singleton_count(self, outdir):
        """Singleton count derived from clusters.parquet minus remaining_clusters.parquet.
        singleton_manifest.tsv is an internal channel only and is never published to outdir."""
        con = duckdb.connect()
        n = con.execute(f"""
            SELECT COUNT(DISTINCT Cluster_ID)
            FROM read_parquet('{outdir}/clusters.parquet')
            WHERE Cluster_ID NOT IN (
                SELECT DISTINCT Cluster_ID
                FROM read_parquet('{outdir}/remaining_clusters.parquet')
            )
        """).fetchone()[0]
        log.info("singleton clusters (derived from parquet): %d", n)
        assert n == SINGLETON_COUNT, f"Expected {SINGLETON_COUNT} singleton cluster(s), got {n}"

    def test_singleton_not_in_remaining(self, outdir):
        """Clusters with exactly 1 protein must not appear in remaining_clusters."""
        con = duckdb.connect()
        singleton_ids = set(
            con.execute(f"""
                SELECT Cluster_ID FROM read_parquet('{outdir}/clusters.parquet')
                GROUP BY Cluster_ID HAVING COUNT(*) = 1
            """).fetchdf()["Cluster_ID"].tolist()
        )
        remaining_ids = set(
            con.execute(
                f"SELECT DISTINCT Cluster_ID FROM read_parquet('{outdir}/remaining_clusters.parquet')"
            ).fetchdf()["Cluster_ID"].tolist()
        )
        log.debug("singleton cluster IDs: %s", singleton_ids)
        log.debug("remaining cluster IDs (sample): %s", list(remaining_ids)[:5])
        overlap = singleton_ids & remaining_ids
        assert not overlap, f"Singleton cluster(s) found in remaining_clusters: {overlap}"
        log.info("PASS: %d singleton cluster(s) correctly absent from remaining_clusters",
                 len(singleton_ids))


# ---------------------------------------------------------------------------
# 2. Per-query FASTA and manifest
# ---------------------------------------------------------------------------

class TestPerQueryOutputs:
    """Validates *_all_relatives.fa and *_manifest.tsv for each query."""

    @pytest.mark.parametrize("query_name,query_tax_id", QUERY_SPECIES.items())
    def test_all_relatives_fasta_exists(self, outdir, query_name, query_tax_id):
        path = os.path.join(outdir, f"{query_name}_all_relatives.fa")
        assert os.path.isfile(path), f"Missing: {path}"

    @pytest.mark.parametrize("query_name,query_tax_id", QUERY_SPECIES.items())
    def test_all_relatives_fasta_sequence_count(self, outdir, query_name, query_tax_id):
        path = os.path.join(outdir, f"{query_name}_all_relatives.fa")
        n = _count_fasta_headers(path)
        log.info("%s_all_relatives.fa: %d sequences", query_name, n)
        assert n == EXPECTED_ALL_RELATIVES, \
            f"{query_name}: expected {EXPECTED_ALL_RELATIVES} sequences, got {n}"

    @pytest.mark.parametrize("query_name,query_tax_id", QUERY_SPECIES.items())
    def test_manifest_exists(self, outdir, query_name, query_tax_id):
        path = os.path.join(outdir, f"{query_name}_manifest.tsv")
        assert os.path.isfile(path), f"Missing: {path}"

    @pytest.mark.parametrize("query_name,query_tax_id", QUERY_SPECIES.items())
    def test_manifest_columns(self, outdir, query_name, query_tax_id):
        path = os.path.join(outdir, f"{query_name}_manifest.tsv")
        df = pd.read_csv(path, sep="\t")
        missing = MANIFEST_COLS - set(df.columns)
        assert not missing, f"{query_name} manifest missing columns: {missing}"

    @pytest.mark.parametrize("query_name,query_tax_id", QUERY_SPECIES.items())
    def test_manifest_row_count(self, outdir, query_name, query_tax_id):
        path = os.path.join(outdir, f"{query_name}_manifest.tsv")
        df = pd.read_csv(path, sep="\t")
        log.info("%s manifest: %d rows", query_name, len(df))
        assert len(df) == EXPECTED_RELATIVES, \
            f"{query_name}: expected {EXPECTED_RELATIVES} manifest rows, got {len(df)}"

    @pytest.mark.parametrize("query_name,query_tax_id", QUERY_SPECIES.items())
    def test_query_not_selected_for_itself(self, outdir, query_name, query_tax_id):
        """The query species' own proteins must never appear as selected relatives."""
        path = os.path.join(outdir, f"{query_name}_manifest.tsv")
        df = pd.read_csv(path, sep="\t")
        self_selected = df[df["source_tax_id"].astype(int) == query_tax_id]
        assert len(self_selected) == 0, \
            f"{query_name}: query's own proteins appear in manifest:\n{self_selected}"

    @pytest.mark.parametrize("query_name,query_tax_id", QUERY_SPECIES.items())
    def test_selection_rank_bounded_by_num_of_rel(self, outdir, query_name, query_tax_id):
        """selection_rank must never exceed num_of_rel (= 5)."""
        path = os.path.join(outdir, f"{query_name}_manifest.tsv")
        df = pd.read_csv(path, sep="\t")
        over_limit = df[df["selection_rank"] > NUM_OF_REL]
        assert len(over_limit) == 0, \
            f"{query_name}: {len(over_limit)} rows have selection_rank > {NUM_OF_REL}"

    @pytest.mark.parametrize("query_name,query_tax_id", QUERY_SPECIES.items())
    def test_manifestquery_tax_id_consistent(self, outdir, query_name, query_tax_id):
        """All rows in the manifest must carry the correct query_tax_id."""
        path = os.path.join(outdir, f"{query_name}_manifest.tsv")
        df = pd.read_csv(path, sep="\t")
        wrong = df[df["query_tax_id"].astype(int) != query_tax_id]
        assert len(wrong) == 0, \
            f"{query_name}: rows with wrong query_tax_id:\n{wrong}"

    @pytest.mark.parametrize("query_name,query_tax_id", QUERY_SPECIES.items())
    def test_selection_source_values(self, outdir, query_name, query_tax_id):
        """
        Relatives must have selection_source='ranked_cluster'.
        Singletons (with_singletons=true) must have selection_source='singleton'.
        No other values are permitted.
        """
        path = os.path.join(outdir, f"{query_name}_manifest.tsv")
        df = pd.read_csv(path, sep="\t")
        allowed = {"ranked_cluster", "singleton"}
        unexpected = set(df["selection_source"].unique()) - allowed
        assert not unexpected, \
            f"{query_name}: unexpected selection_source values: {unexpected}"
        ranked = df[df["selection_source"] == "ranked_cluster"]
        singletons = df[df["selection_source"] == "singleton"]
        log.info("%s: ranked_cluster=%d  singleton=%d", query_name, len(ranked), len(singletons))
        assert len(ranked) == EXPECTED_RELATIVES - SINGLETON_COUNT
        assert len(singletons) == SINGLETON_COUNT


# ---------------------------------------------------------------------------
# 3. Summary reports
# ---------------------------------------------------------------------------

class TestSummaryReports:

    def test_cluster_summary_exists(self, outdir):
        path = os.path.join(outdir, "cluster_summary.txt")
        log.debug("checking cluster_summary.txt: %s", path)
        assert os.path.isfile(path)
        log.info("PASS: cluster_summary.txt found")

    def test_cluster_summary_total_clusters(self, outdir):
        with open(os.path.join(outdir, "cluster_summary.txt"), encoding="utf-8") as fh:
            text = fh.read()
        expected = f"Total number of clusters obtained from input fasta: {TOTAL_CLUSTERS}"
        log.debug("looking for: %r", expected)
        assert expected in text, f"Expected line not found in cluster_summary.txt:\n  {expected}"
        log.info("PASS: cluster_summary total_clusters=%d confirmed", TOTAL_CLUSTERS)

    def test_cluster_summary_singleton_count(self, outdir):
        with open(os.path.join(outdir, "cluster_summary.txt"), encoding="utf-8") as fh:
            text = fh.read()
        expected = f"Singleton clusters (clusters with only one protein): {SINGLETON_COUNT}"
        log.debug("looking for: %r", expected)
        assert expected in text
        log.info("PASS: cluster_summary singleton_count=%d confirmed", SINGLETON_COUNT)

    def test_pipeline_summary_exists(self, outdir):
        path = os.path.join(outdir, "echo_pipeline_summary.txt")
        log.debug("checking echo_pipeline_summary.txt: %s", path)
        assert os.path.isfile(path)
        log.info("PASS: echo_pipeline_summary.txt found")

    def test_pipeline_summary_total_sequences(self, outdir):
        with open(os.path.join(outdir, "echo_pipeline_summary.txt"), encoding="utf-8") as fh:
            text = fh.read()
        expected = f"Total sequences in input FASTA: {TOTAL_INPUT_PROTEINS}"
        log.debug("looking for: %r", expected)
        assert expected in text, \
            f"Expected line not found in echo_pipeline_summary.txt:\n  {expected}"
        log.info("PASS: pipeline_summary total_sequences=%d confirmed", TOTAL_INPUT_PROTEINS)

    def test_pipeline_summary_both_queries_present(self, outdir):
        with open(os.path.join(outdir, "echo_pipeline_summary.txt"), encoding="utf-8") as fh:
            text = fh.read()
        for query_name in QUERY_SPECIES:
            log.debug("checking query '%s' mentioned in pipeline summary", query_name)
            assert query_name in text, \
                f"Query '{query_name}' not mentioned in echo_pipeline_summary.txt"
        log.info("PASS: all %d queries mentioned in pipeline summary", len(QUERY_SPECIES))


# ---------------------------------------------------------------------------
# 4. Diagnostics bundle
# ---------------------------------------------------------------------------

class TestDiagnostics:

    def test_diagnostics_dir_exists(self, outdir):
        path = os.path.join(outdir, "diagnostics_out")
        log.debug("checking diagnostics_out/: %s", path)
        assert os.path.isdir(path), "diagnostics_out/ directory not found"
        log.info("PASS: diagnostics_out/ directory found")

    @pytest.mark.parametrize("fname", sorted(DIAGNOSTICS_FILES))
    def test_diagnostics_file_exists(self, outdir, fname):
        path = os.path.join(outdir, "diagnostics_out", fname)
        log.debug("checking diagnostics file: %s", path)
        assert os.path.isfile(path), f"Missing diagnostics file: {fname}"
        log.info("PASS: diagnostics file present — %s", fname)

    @pytest.mark.parametrize("fname", sorted(DIAGNOSTICS_FILES))
    def test_diagnostics_file_nonempty(self, outdir, fname):
        path = os.path.join(outdir, "diagnostics_out", fname)
        size = os.path.getsize(path)
        log.debug("diagnostics file size: %s = %d bytes", fname, size)
        assert size > 0, f"Diagnostics file is empty: {fname}"
        log.info("PASS: diagnostics file non-empty — %s (%d bytes)", fname, size)


# ---------------------------------------------------------------------------
# 5. with_singletons=false run (nextflow_test_no_singletons)
# ---------------------------------------------------------------------------

# Expected counts for the no-singletons run (same synthetic dataset)
NS_EXPECTED_RELATIVES    = 29  # only ranked_cluster rows; singletons not added to manifest
NS_EXPECTED_ALL_RELATIVES = 29  # FASTA also excludes singletons


@pytest.fixture(scope="module")
def outdir_no_singletons(request):
    """
    Returns the pipeline output directory for the with_singletons=false run.
    If --run-pipeline is passed, executes the pipeline first.
    """
    if request.config.getoption("--run-pipeline"):
        log.info("Running Nextflow pipeline (no-singletons)...")
        result = subprocess.run(
            ["nextflow", "run", "workflows/echo.nf", "-params-file", PARAMS_FILE_NO_SINGLETONS],
            cwd=NF_DIR,
            capture_output=False,
            text=True,
            check=False,
        )
        assert result.returncode == 0, "Nextflow no-singletons pipeline failed"
        log.info("Pipeline (no-singletons) complete.")
    else:
        log.info("Skipping no-singletons pipeline execution (pass --run-pipeline to run it).")

    if not os.path.isdir(OUTDIR_NO_SINGLETONS):
        pytest.skip(
            f"Output directory not found: {OUTDIR_NO_SINGLETONS}\n"
            "Run the pipeline first:  cd repo root && "
            "nextflow run workflows/echo.nf -params-file test/params.test.no_singletons.yaml"
        )
    return OUTDIR_NO_SINGLETONS


class TestNoSingletons:
    """Validates with_singletons=false run output (nextflow_test_no_singletons)."""

    # --- clustering ---

    def test_total_clusters(self, outdir_no_singletons):
        con = duckdb.connect()
        n = con.execute(
            f"SELECT COUNT(DISTINCT Cluster_ID) FROM read_parquet('{outdir_no_singletons}/clusters.parquet')"
        ).fetchone()[0]
        log.info("no-singletons total clusters: %d", n)
        assert n == TOTAL_CLUSTERS

    def test_singleton_cluster_summary_exists(self, outdir_no_singletons):
        path = os.path.join(outdir_no_singletons, "singleton_cluster_summary.tsv")
        log.debug("checking singleton_cluster_summary.tsv: %s", path)
        assert os.path.isfile(path), \
            "singleton_cluster_summary.tsv not found (expected when with_singletons=false)"
        log.info("PASS: singleton_cluster_summary.tsv found")

    def test_singleton_cluster_summary_count(self, outdir_no_singletons):
        df = pd.read_csv(
            os.path.join(outdir_no_singletons, "singleton_cluster_summary.tsv"), sep="\t"
        )
        log.debug("singleton_cluster_summary rows: %d (expected %d)", len(df), SINGLETON_COUNT)
        assert len(df) == SINGLETON_COUNT
        log.info("PASS: singleton_cluster_summary count=%d", len(df))

    def test_singleton_manifest_not_published(self, outdir_no_singletons):
        """singleton_manifest.tsv must never appear in the output dir."""
        path = os.path.join(outdir_no_singletons, "singleton_manifest.tsv")
        exists = os.path.isfile(path)
        log.debug("singleton_manifest.tsv present: %s (expected: absent)", exists)
        assert not exists, "singleton_manifest.tsv should not be published"
        log.info("PASS: singleton_manifest.tsv correctly absent when with_singletons=false")

    def test_discarded_singletons_fa_exists(self, outdir_no_singletons):
        path = os.path.join(outdir_no_singletons, "discarded_singletons.fa")
        log.debug("checking discarded_singletons.fa: %s", path)
        assert os.path.isfile(path)
        log.info("PASS: discarded_singletons.fa found")

    def test_discarded_singletons_fa_count(self, outdir_no_singletons):
        n = _count_fasta_headers(os.path.join(outdir_no_singletons, "discarded_singletons.fa"))
        log.debug("discarded_singletons.fa sequences: %d (expected %d)", n, SINGLETON_COUNT)
        assert n == SINGLETON_COUNT
        log.info("PASS: discarded_singletons.fa has %d sequences", n)

    # --- per-query outputs ---

    @pytest.mark.parametrize("query_name,query_tax_id", QUERY_SPECIES.items())
    def test_all_relatives_fasta_sequence_count(self, outdir_no_singletons, query_name, query_tax_id):
        path = os.path.join(outdir_no_singletons, f"{query_name}_all_relatives.fa")
        assert os.path.isfile(path), f"Missing: {path}"
        n = _count_fasta_headers(path)
        assert n == NS_EXPECTED_ALL_RELATIVES, \
            f"{query_name}: expected {NS_EXPECTED_ALL_RELATIVES} sequences, got {n}"

    @pytest.mark.parametrize("query_name,query_tax_id", QUERY_SPECIES.items())
    def test_manifest_row_count(self, outdir_no_singletons, query_name, query_tax_id):
        path = os.path.join(outdir_no_singletons, f"{query_name}_manifest.tsv")
        df = pd.read_csv(path, sep="\t")
        assert len(df) == NS_EXPECTED_RELATIVES, \
            f"{query_name}: expected {NS_EXPECTED_RELATIVES} rows, got {len(df)}"

    @pytest.mark.parametrize("query_name,query_tax_id", QUERY_SPECIES.items())
    def test_manifest_only_ranked_cluster(self, outdir_no_singletons, query_name, query_tax_id):
        """With singletons=false, manifest must contain only ranked_cluster rows."""
        path = os.path.join(outdir_no_singletons, f"{query_name}_manifest.tsv")
        df = pd.read_csv(path, sep="\t")
        unexpected = set(df["selection_source"].unique()) - {"ranked_cluster"}
        assert not unexpected, \
            f"{query_name}: unexpected selection_source values: {unexpected}"

    @pytest.mark.parametrize("query_name,query_tax_id", QUERY_SPECIES.items())
    def test_query_not_selected_for_itself(self, outdir_no_singletons, query_name, query_tax_id):
        path = os.path.join(outdir_no_singletons, f"{query_name}_manifest.tsv")
        df = pd.read_csv(path, sep="\t")
        self_selected = df[df["source_tax_id"].astype(int) == query_tax_id]
        assert len(self_selected) == 0

    # --- diagnostics ---

    def test_diagnostics_dir_exists(self, outdir_no_singletons):
        path = os.path.join(outdir_no_singletons, "diagnostics_out")
        log.debug("checking diagnostics_out/: %s", path)
        assert os.path.isdir(path)
        log.info("PASS: diagnostics_out/ found in no-singletons outdir")

    @pytest.mark.parametrize("fname", sorted(DIAGNOSTICS_FILES))
    def test_diagnostics_file_exists(self, outdir_no_singletons, fname):
        path = os.path.join(outdir_no_singletons, "diagnostics_out", fname)
        log.debug("checking diagnostics file (no-singletons): %s", path)
        assert os.path.isfile(path), f"Missing diagnostics file: {fname}"
        log.info("PASS: diagnostics file present (no-singletons) — %s", fname)
