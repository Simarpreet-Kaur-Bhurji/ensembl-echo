/* Step 8b (manifest): Merge all per-chunk manifest TSVs for a single query into one file.
 * Mirrors MERGE_FASTAS_PER_QUERY: groups by query, produces one manifest per query species.
 * Input:  list of partial manifest TSVs for one query (all chunks)
 * Output: {query_name}_manifest.tsv — flat per-protein provenance file, published to outdir
 *         Genebuild uses this to trace every protein in the FASTA back to its cluster and source taxon.
 */
process MERGE_MANIFESTS_PER_QUERY {
  tag { query_name }

  input:
    tuple val(query_tax_id), val(query_name), path(manifest_parts), path(singleton_manifest)

  output:
    path "${query_name}_manifest.tsv", emit: query_manifest

  script:
  """
  set -euo pipefail

  out="${query_name}_manifest.tsv"
  files=( ${manifest_parts.join(' ')} )

  # write header from first file, then data rows from all files (skip header on each)
  head -n 1 "\${files[0]}" > "\$out"
  for f in "\${files[@]}"; do
    tail -n +2 "\$f" >> "\$out"
  done

  # append singleton rows if the singleton manifest has data (with_singletons=true)
  if [ -s "${singleton_manifest}" ] && [ \$(wc -l < "${singleton_manifest}") -gt 1 ]; then
    python3 -c "
import pandas as pd, sys
df = pd.read_csv('${singleton_manifest}', sep='\t', dtype=str)
df['query_tax_id'] = '${query_tax_id}'
df['query_name']   = '${query_name}'
df.to_csv(sys.stdout, sep='\t', index=False, header=False)
" >> "\$out"
  fi
  """
}
