/* Step 9: Assemble the final all_relatives FASTA for each query.
 * Input:  relatives_fa (from MERGE_FASTAS_PER_QUERY), singletons_fa
 * Output: {query_name}_all_relatives.fa, optional dedup_report.tsv
 * Singletons are included by default (params.with_singletons = true).
 * Header dedup always runs. Exact-sequence dedup runs when params.dedup_sequences = true,
 * producing a dedup_report.tsv per query listing dropped headers and their duplicates.
 */
process MAKE_ALL_RELATIVES {
  tag { query_name }
  publishDir params.outdir, mode: 'copy'

  input:
    tuple val(query_name), path(relatives_fa), path(singletons_fa)

  output:
    path "${query_name}_all_relatives.fa", emit: all_relatives_fa
    path "dedup_report.tsv",               optional: true, emit: dedup_report

  script:
  """
  if ${params.with_singletons} ; then
    echo_dedup_fasta.py \
      --inputs ${relatives_fa} ${singletons_fa} \
      --output ${query_name}_all_relatives.fa \
      ${params.dedup_sequences ? '--dedup_sequences' : ''}
  else
    echo_dedup_fasta.py \
      --inputs ${relatives_fa} \
      --output ${query_name}_all_relatives.fa \
      ${params.dedup_sequences ? '--dedup_sequences' : ''}
  fi
  """
}

