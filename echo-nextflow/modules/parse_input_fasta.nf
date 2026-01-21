process PARSE_INPUT_FASTA {

  tag "parse_input_fasta"
  publishDir params.outdir, mode: 'copy'

  input:
    path input_fasta_dir
    path metadata_tsv

  output:
    path "combined_input_fasta.fa", emit: combined_fasta
    path "processed_input.tsv",     emit: processed_tsv
    path "processed_input.parquet", emit: processed_parquet

  script:
  """
  echo_parse_input.py \
    --input_fasta_dir ${input_fasta_dir} \
    --metadata_tsv ${metadata_tsv} \
    --out_fasta combined_input_fasta.fa \
    --out_tsv processed_input.tsv \
    --out_parquet processed_input.parquet
  """
}

