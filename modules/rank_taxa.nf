/* Step 5: Compute pairwise taxonomic distances between query species and all input species.
 * Input:  processed_input.parquet, query_species TSV
 * Output: ranked_taxa.tsv — sorted by ascending distance for each query taxon
 */
process RANK_TAXA {
  tag "rank_taxa"
  publishDir params.outdir, mode: 'copy'

  input:
    path processed_input_parquet
    path query_species_tsv

  output:
    path "ranked_taxa.tsv", emit: ranked_taxa_tsv

  script:
  // ncbi_taxa_db is optional; omit the flag to use ete3's default ~/.etetoolkit/taxa.sqlite
  def ncbi_db_arg = params.ncbi_taxa_db
    ? "--ncbi_taxa_db ${params.ncbi_taxa_db}"
    : ""
  """
  echo_rank_taxa.py \
    --query_species ${query_species_tsv} \
    --processed_input_parquet ${processed_input_parquet} \
    --out_tsv ranked_taxa.tsv \
    ${ncbi_db_arg}
  """
}

