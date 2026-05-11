/* Step 9: Assemble the final all_relatives FASTA for each query.
 * Input:  relatives_fa (from MERGE_FASTAS_PER_QUERY), singletons_fa, manifest_tsv
 * Output: {query_name}_all_relatives.fa, {query_name}_manifest.tsv, optional dedup_report.tsv
 * Singletons are included by default (params.with_singletons = true).
 * Header dedup always runs. Exact-sequence dedup runs when params.dedup_sequences = true,
 * producing a dedup_report.tsv per query listing dropped headers and their duplicates.
 * When dedup_report.tsv is produced, manifest rows for dropped headers are marked
 * selection_source=dropped_duplicate to keep the manifest in sync with the final FASTA.
 */
process MAKE_ALL_RELATIVES {
  tag { query_name }
  publishDir params.outdir, mode: 'copy'

  input:
    tuple val(query_name), path(relatives_fa), path(singletons_fa), path(manifest_tsv)

  output:
    path "${query_name}_all_relatives.fa", emit: all_relatives_fa
    path "${query_name}_manifest.tsv",     emit: query_manifest
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

  python3 -c "
import pandas as pd, os
dropped = set()
if os.path.exists('dedup_report.tsv'):
    report = pd.read_csv('dedup_report.tsv', sep='\t', dtype=str)
    dropped = set(report['dropped_header'])
df = pd.read_csv('${manifest_tsv}', sep='\t', dtype=str)
if dropped:
    df.loc[df['protein_header'].isin(dropped), 'selection_source'] = 'dropped_duplicate'
out = '${query_name}_manifest.tsv'
if os.path.islink(out):
    os.unlink(out)
df.to_csv(out, sep='\t', index=False)
"
  """
}

