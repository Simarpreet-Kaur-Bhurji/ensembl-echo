/* Step 9: Assemble the final all_relatives FASTA for each query.
 * Input:  relatives_fa (from MERGE_FASTAS_PER_QUERY), fewer_tax_fa, singletons_fa
 * Output: {query_name}_all_relatives.fa
 * Singletons are included only when params.with_singletons is true.
 */
process MAKE_ALL_RELATIVES {
  tag { query_name }
  publishDir params.outdir, mode: 'copy'

  input:
    tuple val(query_name), path(relatives_fa), path(fewer_tax_fa), path(singletons_fa)

  output:
    path "${query_name}_all_relatives.fa", emit: all_relatives_fa

  script:
  """
  if ${params.with_singletons} ; then
    cat ${relatives_fa} ${fewer_tax_fa} ${singletons_fa} > ${query_name}_all_relatives.fa
  else
    cat ${relatives_fa} ${fewer_tax_fa} > ${query_name}_all_relatives.fa
  fi
  """
}

