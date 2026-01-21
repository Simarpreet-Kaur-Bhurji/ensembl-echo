process MERGE_FASTAS_PER_QUERY {
  tag { query_name }
  publishDir params.outdir, mode: 'copy'

  input:
    tuple val(query_tax_id), val(query_name), path(fasta_parts)

  output:
    path "${query_name}_relatives.fa", emit: relatives_fa

  script:
  """
  # fasta_parts is a list of files; cat works for 1 or many
  cat ${fasta_parts.join(' ')} > ${query_name}_relatives.fa
  """
}

