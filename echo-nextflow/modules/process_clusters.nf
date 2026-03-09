/* Step 4: Filter clusters by taxon diversity and singleton status.
 * Input:  clusters.parquet
 * Output: remaining_clusters.parquet (pass-through clusters)
 *         discarded_singletons.fa + summary (optional, if any singletons)
 *         clusters_with_fewer_tax_ids.fa + summary (optional, if any under-represented)
 */
process PROCESS_CLUSTERS {

  tag "process_clusters"
  publishDir params.outdir, mode: 'copy'

  input:
    path clusters_parquet

  output:
    path "remaining_clusters.parquet",              emit: remaining_clusters
    path "discarded_singletons.fa",                 optional: true, emit: singletons_fa
    path "singleton_cluster_summary.tsv",           optional: true, emit: singletons_summary
    path "clusters_with_fewer_tax_ids.fa",          optional: true, emit: fewer_tax_fa
    path "clusters_with_fewer_taxids_summary.tsv",  optional: true, emit: few_taxids_summary

  script:
  """
  echo_process_clusters.py \
    --num_rel ${params.num_of_rel} \
    --out_remaining remaining_clusters.parquet
  """
}

