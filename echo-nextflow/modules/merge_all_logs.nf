process MERGE_ALL_LOGS {
  tag "merge_all_logs"
  publishDir params.outdir, mode: 'copy'

  input:
    path logs

  output:
    path "closest_relatives_log.tsv", emit: merged_log

  script:
  """
  # keep only first header
  first=\$(ls -1 ${logs} | head -n 1)
  head -n 1 "\$first" > closest_relatives_log.tsv || true
  for f in ${logs}; do
    tail -n +2 "\$f" >> closest_relatives_log.tsv || true
  done
  """
}

