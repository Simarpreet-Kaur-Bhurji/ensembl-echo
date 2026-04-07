nextflow.enable.dsl = 2

include { PARSE_INPUT_FASTA }        from '../modules/parse_input_fasta.nf'
include { RUN_MMSEQS }               from '../modules/run_mmseqs.nf'
include { PARSE_CLUSTERS }           from '../modules/parse_clusters.nf'
include { PROCESS_CLUSTERS }         from '../modules/process_clusters.nf'

include { RANK_TAXA }                from '../modules/rank_taxa.nf'
include { SPLIT_REMAINING_CLUSTERS } from '../modules/split_remaining_clusters.nf'
include { CLOSEST_RELATIVES_CHUNK }    from '../modules/closest_relatives_chunk.nf'
include { MERGE_FASTAS_PER_QUERY }     from '../modules/merge_fastas_per_query.nf'
include { MERGE_MANIFESTS_PER_QUERY }  from '../modules/merge_manifests_per_query.nf'  // new: per-query manifest merge; replaces MERGE_ALL_LOGS
include { MAKE_ALL_RELATIVES }       from '../modules/make_all_relatives.nf'
include { MAKE_REPORTS }             from '../modules/make_reports.nf'
include { MAKE_DIAGNOSTICS }         from '../modules/make_diagnostics.nf'

workflow ECHO {

  main:
  /*
   * Decide whether to reuse existing clustering
   */
  def use_existing = params.existing_clusters_dir != null
  
  /*
   * Define inputs depending on mode
   */
  if (use_existing) {

    def base = params.existing_clusters_dir

    // Hard prerequisites — cannot be derived from anything else; fail fast if missing
    assert file("${base}/processed_input.parquet").exists()
    assert file("${base}/clusters.parquet").exists()
    assert file("${base}/remaining_clusters.parquet").exists()
    assert file("${base}/combined_input_fasta.fa").exists()

    processed_parquet  = channel.value(file("${base}/processed_input.parquet"))
    clusters_parquet   = channel.value(file("${base}/clusters.parquet"))
    remaining_clusters = channel.value(file("${base}/remaining_clusters.parquet"))
    combined_fasta     = channel.value(file("${base}/combined_input_fasta.fa"))

    // Singleton files are derived outputs of clusters.parquet — fast to regenerate
    // (DuckDB + pandas, no MMseqs2). May be absent if the previous run used
    // with_singletons=false or output was partially cleaned up.
    def singletons_fa_file  = file("${base}/discarded_singletons.fa")
    def singletons_tsv_file = file("${base}/singleton_cluster_summary.tsv")

    if (singletons_fa_file.exists() && singletons_tsv_file.exists()) {
      singletons_fa       = channel.value(singletons_fa_file)
      singletons_summary  = channel.value(singletons_tsv_file)
      singletons_manifest = channel.empty()
    } else {
      // Regenerate from clusters.parquet — PROCESS_CLUSTERS also emits
      // remaining_clusters but we use the one from base dir above
      log.info "[ECHO] Singleton files missing from ${base} — regenerating from clusters.parquet"
      def regen           = PROCESS_CLUSTERS(clusters_parquet)
      singletons_fa       = regen.singletons_fa
      singletons_summary  = regen.singletons_summary
      singletons_manifest = regen.singletons_manifest
    }

  } else {

    /*
     * Full pipeline from scratch
     */
    parsed = PARSE_INPUT_FASTA(
      file(params.metadata_tsv)
    )

    mm = RUN_MMSEQS(parsed.combined_fasta)

    cl = PARSE_CLUSTERS(
      mm.cluster_tsv,
      parsed.processed_parquet
    )

    filtered = PROCESS_CLUSTERS(cl.clusters_parquet)

    processed_parquet     = parsed.processed_parquet
    clusters_parquet      = cl.clusters_parquet
    remaining_clusters    = filtered.remaining_clusters
    singletons_fa         = filtered.singletons_fa
    singletons_summary    = filtered.singletons_summary
    singletons_manifest   = filtered.singletons_manifest
    combined_fasta        = parsed.combined_fasta
  }

  /*
   * Rank taxa (always recomputed)
   */
  ranked = RANK_TAXA(
    processed_parquet,
    file(params.query_species)
  )

  /*
   * Split clusters into chunks
   */
  chunks = SPLIT_REMAINING_CLUSTERS(remaining_clusters)
  chunk_files = chunks.chunks.flatten()

  /*
   * Load queries
   */
  queries = channel
    .fromPath(params.query_species)
    .splitCsv(sep: '\t', header: true)
    .map { row -> tuple(row.tax_id.toString(), row.sps_name.toString()) }
    .distinct { item -> item[0] }

  /*
   * Cartesian product: (chunk × query × ranked_taxa)
   */
  closest_inputs = chunk_files
    .combine(queries)
    .combine(ranked.ranked_taxa_tsv)
    .map { it ->
      tuple(it[0], it[3], it[1], it[2])
    }

  /*
   * Closest relatives per chunk
   */
  partial = CLOSEST_RELATIVES_CHUNK(closest_inputs)

  /*
   * Merge FASTAs per query
   */
  fastas_grouped = partial.partial_fastas
    .map { fa ->
      def tid = fa.baseName.tokenize('_')[0]
      tuple(tid, fa)
    }
    .groupTuple()

  merged_fastas = fastas_grouped
    .combine(queries)
    .filter { tid1, _fas, tid2, _name -> tid1 == tid2 }
    .map { tid, fas, _unused, name ->
      tuple(tid, name.toLowerCase().replaceAll(' ', '_'), fas)
    }
    | MERGE_FASTAS_PER_QUERY

  /*
   * Merge manifest TSVs per query; mirrors the FASTA merge above.
   * Each {query_name}_manifest.tsv is the provenance handoff file for Genebuild,
   * published alongside the FASTA.
   */
  manifests_grouped = partial.partial_manifests
    .map { tsv ->
      def tid = tsv.baseName.tokenize('_')[0]
      tuple(tid, tsv)
    }
    .groupTuple()

  merged_manifests = manifests_grouped
    .combine(queries)
    .filter { tid1, _tsvs, tid2, _name -> tid1 == tid2 }
    .map { tid, tsvs, _unused, name ->
      tuple(tid, name.toLowerCase().replaceAll(' ', '_'), tsvs)
    }
    .combine(singletons_manifest.ifEmpty([file("NO_SINGLETON_MANIFEST")]))
    | MERGE_MANIFESTS_PER_QUERY

  /*
   * Build *_all_relatives.fa
   */
  all_inputs = merged_fastas.relatives_fa
    .map { fa ->
      def qname = fa.baseName.replaceFirst(/_relatives$/, '')
      tuple(qname, fa)
    }
    .combine(singletons_fa)
    .map { it ->
      tuple(it[0], it[1], it[2])
    }

  all_rel = MAKE_ALL_RELATIVES(all_inputs)

  /*
   * Diagnostics & reports
   */

  reports = MAKE_REPORTS(
    clusters_parquet,
    remaining_clusters,
    singletons_fa,
    all_rel.all_relatives_fa.collect(),
    combined_fasta
  )

  diagnostics = MAKE_DIAGNOSTICS(
    clusters_parquet,
    remaining_clusters,
    merged_manifests.query_manifest.collect(),  // all per-query manifests collected; diagnostics reads distances from them directly instead of a separate combined log
    combined_fasta,
    reports.cluster_summary,
    reports.pipeline_summary
  )

  emit:
    clusters_parquet          = clusters_parquet
    remaining_clusters        = remaining_clusters
    ranked_taxa_tsv           = ranked.ranked_taxa_tsv
    query_manifests           = merged_manifests.query_manifest  // per-query provenance TSVs; handoff contract for Genebuild
    relatives_fastas          = merged_fastas.relatives_fa
    all_relatives_fa          = all_rel.all_relatives_fa
    cluster_summary_txt       = reports.cluster_summary
    echo_pipeline_summary_txt = reports.pipeline_summary
    diagnostics_pdf           = diagnostics.diagnostics_pdf
    diagnostics_summary       = diagnostics.diagnostics_summary
    diagnostics_dir           = diagnostics.diagnostics_dir

}
