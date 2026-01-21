process CLOSEST_RELATIVES_CHUNK {
  tag { "${query_tax_id}_${chunk_file.simpleName}" }

  input:
    tuple path(chunk_file), path(ranked_taxa_tsv), val(query_tax_id), val(query_name)

  output:
    path "${query_tax_id}_${chunk_file.simpleName}.fa",  emit: partial_fastas
    path "${query_tax_id}_${chunk_file.simpleName}.tsv", emit: partial_logs

  script:
  """
  echo_closest_relatives_chunk.py \
    --chunk_parquet ${chunk_file} \
    --ranked_taxa_tsv ${ranked_taxa_tsv} \
    --query_tax_id ${query_tax_id} \
    --query_name '${query_name}' \
    --num_of_rel ${params.num_of_rel} \
    --out_fasta ${query_tax_id}_${chunk_file.simpleName}.fa \
    --out_log ${query_tax_id}_${chunk_file.simpleName}.tsv
  """
}

