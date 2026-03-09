/* Step 6: Partition remaining_clusters.parquet into fixed-size chunk files.
 * Input:  remaining_clusters.parquet
 * Output: cluster_chunks/chunk_*.parquet (params.clusters_per_chunk clusters each)
 * Chunking enables embarrassingly parallel closest-relative searches.
 */
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

