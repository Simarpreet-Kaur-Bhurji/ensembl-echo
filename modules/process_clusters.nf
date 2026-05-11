/* Step 4: Filter clusters by singleton status and write remaining.
 * Input:  clusters.parquet
 * Output: remaining_clusters.parquet (all multi-member clusters)
 *         discarded_singletons.fa + summary (optional, if any singletons)
 * Publishing rules:
 *   - discarded_singletons.fa is only published when with_singletons=false
 *     (when included in the final FASTA, publishing separately is misleading)
 *   - remaining_clusters.parquet is suppressed in restart mode (already exists in base dir)
 */
process PROCESS_CLUSTERS {

  tag "process_clusters"
  publishDir params.outdir, mode: 'copy', saveAs: { filename ->
    if (filename == "remaining_clusters.parquet" && params.existing_clusters_dir) return null
    if (filename == "discarded_singletons.fa" && params.with_singletons) return null
    if (filename == "singleton_manifest.tsv") return null
    return filename
  }

  input:
    path clusters_parquet

  output:
    path "remaining_clusters.parquet",              emit: remaining_clusters
    path "discarded_singletons.fa",                 optional: true, emit: singletons_fa
    path "singleton_manifest.tsv",                  optional: true, emit: singletons_manifest
    path "singleton_cluster_summary.tsv",           optional: true, emit: singletons_summary

  script:
  def with_singletons_flag = params.with_singletons ? "--with_singletons" : ""
  """
  echo_process_clusters.py \
    --out_remaining remaining_clusters.parquet \
    ${with_singletons_flag}
  """
}

