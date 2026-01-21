process MAKE_REPORTS {
  tag "make_reports"
  publishDir params.outdir, mode: 'copy'

  input:
    path clusters_parquet
    path remaining_clusters
    path few_taxids_summary_tsv
    path singletons_fa
    path fewer_tax_fa
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
    --clusters_with_fewer_taxids_summary_tsv ${few_taxids_summary_tsv} \
    --discarded_singletons_fa ${singletons_fa} \
    --clusters_with_fewer_tax_ids_fa ${fewer_tax_fa} \
    --all_relatives_fastas ${all_relatives_fastas} \
    --input_fasta ${input_fasta} \
    ${params.with_singletons ? "--with_singletons" : ""} \
    --out_cluster_summary cluster_summary.txt \
    --out_pipeline_summary echo_pipeline_summary.txt
  """
}

