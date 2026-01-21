process RANK_TAXA {
  tag "rank_taxa"
  publishDir params.outdir, mode: 'copy'

  input:
    path processed_input_parquet
    path query_species_tsv

  output:
    path "ranked_taxa.tsv", emit: ranked_taxa_tsv

  script:
  """
  echo_rank_taxa.py \
    --query_species ${query_species_tsv} \
    --processed_input_parquet ${processed_input_parquet} \
    --out_tsv ranked_taxa.tsv
  """
}

