/* Step 10b: Generate QC diagnostics PDF and summary.
 * Input:  clusters_parquet, remaining_clusters, manifest_tsvs,
 *         input_fasta, cluster_summary_txt, pipeline_summary_txt
 * Output: diagnostics_out/ directory, diagnostics.pdf, diagnostics_summary.txt
 * Note:   singleton count derived from clusters_parquet (cluster_size==1) —
 *         no separate singleton file needed.
 */

process MAKE_DIAGNOSTICS {
  tag "diagnostics"
  publishDir params.outdir, mode: 'copy'
  cache false

  input:
    path clusters_parquet
    path remaining_clusters
    path manifest_tsvs  // all per-query manifests collected; replaces the single closest_relatives_log.tsv
    path input_fasta
    path cluster_summary_txt
    path pipeline_summary_txt

  output:
    path "diagnostics_out", emit: diagnostics_dir
    path "diagnostics_out/diagnostics.pdf", emit: diagnostics_pdf
    path "diagnostics_out/diagnostics_summary.txt", emit: diagnostics_summary

  script:
  """
  mkdir -p diagnostics_out

  echo_make_diagnostics.py \
    --clusters_parquet ${clusters_parquet} \
    --remaining_clusters_parquet ${remaining_clusters} \
    --manifest_tsvs ${manifest_tsvs} \
    --input_fasta ${input_fasta} \
    --cluster_summary_txt ${cluster_summary_txt} \
    --pipeline_summary_txt ${pipeline_summary_txt} \
    ${params.with_singletons ? "--with_singletons" : ""} \
    --outdir diagnostics_out
  """
}
