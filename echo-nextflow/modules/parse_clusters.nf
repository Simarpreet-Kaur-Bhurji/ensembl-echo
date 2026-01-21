process PARSE_CLUSTERS {

  tag "parse_clusters"
  publishDir params.outdir, mode: 'copy'

  input:
    path cluster_tsv
    path processed_input_parquet

  output:
    path "clusters.parquet", emit: clusters_parquet

  script:
  """
  echo_parse_clusters.py \
    --cluster_tsv ${cluster_tsv} \
    --processed_input_parquet ${processed_input_parquet} \
    --out_parquet clusters.parquet
  """
}

