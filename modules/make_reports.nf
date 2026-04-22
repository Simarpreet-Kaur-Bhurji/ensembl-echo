/* Step 10a: Generate human-readable summary reports.
 * Input:  clusters_parquet, remaining_clusters, singletons_fa,
 *         all_relatives_fastas, input_fasta
 * Output: cluster_summary.txt, echo_pipeline_summary.txt
 */
process MAKE_REPORTS {
  tag "make_reports"
  publishDir params.outdir, mode: 'copy'

  input:
    path clusters_parquet
    path remaining_clusters
    path singletons_fa
    path all_relatives_fastas
    path input_fasta

  output:
    path "cluster_summary.txt",        emit: cluster_summary
    path "echo_pipeline_summary.txt",  emit: pipeline_summary

  script:
  """
  echo_make_reports.py \
    --clusters_parquet ${clusters_parquet} \
    --remaining_clusters_parquet ${remaining_clusters} \
    --discarded_singletons_fa ${singletons_fa} \
    --all_relatives_fastas ${all_relatives_fastas} \
    --input_fasta ${input_fasta} \
    ${params.with_singletons ? "--with_singletons" : ""} \
    --out_cluster_summary cluster_summary.txt \
    --out_pipeline_summary echo_pipeline_summary.txt
  """
}

