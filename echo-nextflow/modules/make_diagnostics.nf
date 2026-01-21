process MAKE_DIAGNOSTICS {
  tag "diagnostics"
  publishDir params.outdir, mode: 'copy'
  cache false

  input:
    path clusters_parquet
    path singleton_summary_tsv
    path remaining_clusters
    path few_taxids_summary_tsv
    path closest_log_tsv
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
    --singleton_summary_tsv ${singleton_summary_tsv} \
    --remaining_clusters_parquet ${remaining_clusters} \
    --few_taxids_summary_tsv ${few_taxids_summary_tsv} \
    --closest_log_tsv ${closest_log_tsv} \
    --input_fasta ${input_fasta} \
    --cluster_summary_txt ${cluster_summary_txt} \
    --pipeline_summary_txt ${pipeline_summary_txt} \
    ${params.with_singletons ? "--with_singletons" : ""} \
    --outdir diagnostics_out
  """
}

