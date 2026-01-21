import pandas as pd
import duckdb


def create_hcp_table(metadata_pq_file, con):
    """
    Create a DuckDB table 'hcp_table' from a Parquet file.
    """
    # con.execute(
    #    f"""
    # CREATE TABLE hcp_table AS SELECT * FROM '{metadata_pq_file}'
    # """
    # )

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
    # print(f"DuckDB table 'hcp_table' created from Parquet file: {metadata_pq_file}")

    # Fetch and print the data to verify
    #df = con.execute("SELECT * FROM hcp_table").fetchdf()
    #print(df['tax_id'].head())
    return


def get_hcp_tax_ids(con):
    """
    Get unique tax IDs from the 'hcp_table'.
    """

    # Fetch unique tax IDs from the DuckDB table
    tax_ids = con.execute("SELECT DISTINCT tax_id FROM hcp_table").fetchdf()

    return tax_ids["tax_id"].tolist()

