process SPLIT_REMAINING_CLUSTERS {
  tag "split_clusters"

  input:
    path remaining_clusters_parquet

  output:
    path "cluster_chunks/chunk_*.parquet", emit: chunks

  script:
  """
  echo_split_clusters.py \
    --in_parquet ${remaining_clusters_parquet} \
    --out_dir cluster_chunks \
    --clusters_per_chunk ${params.clusters_per_chunk}
  """
}

