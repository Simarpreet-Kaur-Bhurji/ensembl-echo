/* Step 8a: Concatenate all per-chunk closest-relative log TSVs into one file.
 * Input:  all partial log TSVs from CLOSEST_RELATIVES_CHUNK (collected)
 * Output: closest_relatives_log.tsv — single log covering all queries and chunks
 */
process MERGE_ALL_LOGS {
  tag "merge_all_logs"
  publishDir params.outdir, mode: 'copy'

  input:
    path logs

  output:
    path "closest_relatives_log.tsv", emit: merged_log

  script:
  """
  set -euo pipefail

  out="closest_relatives_log.tsv"
  : > "\$out"

  files=( ${logs} )

  if [[ \${#files[@]} -eq 0 ]]; then
    exit 0
  fi

  head -n 1 "\${files[0]}" > "\$out" || true

  for f in "\${files[@]}"; do
    tail -n +2 "\$f" >> "\$out" || true
  done
  """
}

