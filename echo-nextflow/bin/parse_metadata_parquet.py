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
parse_metadata_parquet.py
--------------------------
DuckDB helpers for loading processed_input.parquet into an in-memory table
and querying unique tax_ids from it.

Used by the rank_taxa step to determine which species are present in the
input protein set before computing taxonomic distances.
"""


def create_hcp_table(metadata_pq_file, con):
    """
    Load processed_input.parquet into a DuckDB table named 'hcp_table'.

    Explicit column selection guards against schema changes in the parquet
    unexpectedly breaking downstream queries.

    Args:
        metadata_pq_file : path to processed_input.parquet
        con              : open DuckDB connection
    """
    con.execute(
        f"""
        CREATE OR REPLACE TABLE hcp_table AS
        SELECT
            header,
            sequence,
            protein_id,
            name,
            tax_id,
            confidence_score,
            confidence_level,
            seq_len
        FROM read_parquet('{metadata_pq_file}')
        """
    )


def get_hcp_tax_ids(con):
    """
    Return the list of unique tax_ids present in 'hcp_table'.

    Args:
        con : open DuckDB connection with hcp_table loaded

    Returns:
        list of tax_id values (type matches what is stored in the parquet)
    """
    tax_ids = con.execute("SELECT DISTINCT tax_id FROM hcp_table").fetchdf()
    return tax_ids["tax_id"].tolist()
