process RUN_MMSEQS {

  tag "mmseqs"
  publishDir params.outdir, mode: 'copy'

  cpus params.mmseqs_threads

  input:
    path combined_fasta

  output:
    path "mmseqs_results_cluster.tsv", emit: cluster_tsv

  script:
  """
  echo_run_mmseqs.sh \
    --input_fasta ${combined_fasta} \
    --out_prefix mmseqs_results \
    --tmp_dir tmp \
    --singularity_image ${params.mmseqs_singularity_image} \
    --min_seq_id ${params.min_seq_id} \
    --coverage ${params.coverage} \
    --cov_mode ${params.cov_mode} \
    --threads ${task.cpus}
  """
}

