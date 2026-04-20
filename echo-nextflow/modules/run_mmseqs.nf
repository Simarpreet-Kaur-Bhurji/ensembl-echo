/* Step 2: cluster all input proteins with MMseqs2 easy-cluster via container.
 * Input:  combined_input_fasta.fa  (from PARSE_INPUT_FASTA)
 * Output: mmseqs_results_cluster.tsv — two-column (seed, member) TSV
 * Threads are set per-profile in nextflow.config (e.g., 16 for slurm, 1 for local).
 */
process RUN_MMSEQS {

  tag "mmseqs"
  publishDir params.outdir, mode: 'copy'

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
    --min_seq_id ${params.min_seq_id} \
    --coverage ${params.coverage} \
    --cov_mode ${params.cov_mode} \
    --threads ${task.cpus}
  """
}

