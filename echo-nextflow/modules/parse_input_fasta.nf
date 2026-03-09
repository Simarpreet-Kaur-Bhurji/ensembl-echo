/* Step 1: Parse the metadata TSV and combine per-species FASTAs into one file.
 * Input:  metadata_tsv (4-col: sps_name, taxon_id, gca, file_path)
 * Output: combined_input_fasta.fa, processed_input.tsv, processed_input.parquet
 * Reads file_path entries; resolves relative paths against params.input_fasta_dir.
 */
process PARSE_INPUT_FASTA {

  tag "parse_input_fasta"
  publishDir params.outdir, mode: 'copy'

  input:
    path metadata_tsv

  output:
    path "combined_input_fasta.fa", emit: combined_fasta
    path "processed_input.tsv",     emit: processed_tsv
    path "processed_input.parquet", emit: processed_parquet

  script:
  // input_fasta_dir is optional: used as a base directory for relative file_path
  // entries in the metadata TSV.  Absolute file_path values in the TSV are used as-is.
  def base_dir_arg = params.input_fasta_dir
    ? "--input_fasta_dir ${file(params.input_fasta_dir).toAbsolutePath()}"
    : ""
  """
  echo_parse_input.py \
    ${base_dir_arg} \
    --metadata_tsv ${metadata_tsv} \
    --out_fasta combined_input_fasta.fa \
    --out_tsv processed_input.tsv \
    --out_parquet processed_input.parquet \
    --fasta_line_width ${params.fasta_line_width}
  """
}

