nextflow.enable.dsl = 2

include { PARSE_INPUT_FASTA }        from '../modules/parse_input_fasta.nf'
include { RUN_MMSEQS }               from '../modules/run_mmseqs.nf'
include { PARSE_CLUSTERS }           from '../modules/parse_clusters.nf'
include { PROCESS_CLUSTERS }         from '../modules/process_clusters.nf'

include { RANK_TAXA }                from '../modules/rank_taxa.nf'
include { SPLIT_REMAINING_CLUSTERS } from '../modules/split_remaining_clusters.nf'
include { CLOSEST_RELATIVES_CHUNK }  from '../modules/closest_relatives_chunk.nf'
include { MERGE_ALL_LOGS }           from '../modules/merge_all_logs.nf'
include { MERGE_FASTAS_PER_QUERY }   from '../modules/merge_fastas_per_query.nf'
include { MAKE_ALL_RELATIVES }       from '../modules/make_all_relatives.nf'
include { MAKE_REPORTS }             from '../modules/make_reports.nf'
include { MAKE_DIAGNOSTICS }         from '../modules/make_diagnostics.nf'

workflow ECHO {

  /*
   * Decide whether to reuse existing clustering
   */
  boolean use_existing = params.existing_clusters_dir != null
  
  /*
   * Define inputs depending on mode
   */
  if (use_existing) {

    def base = params.existing_clusters_dir
    // Fail fast if required files are missing
    assert file("${base}/processed_input.parquet").exists()
    assert file("${base}/clusters.parquet").exists()
    assert file("${base}/remaining_clusters.parquet").exists()
    assert file("${base}/clusters_with_fewer_tax_ids.fa").exists()
    assert file("${base}/clusters_with_fewer_taxids_summary.tsv").exists()
    assert file("${base}/discarded_singletons.fa").exists()
    assert file("${base}/singleton_cluster_summary.tsv").exists()
    assert file("${base}/combined_input_fasta.fa").exists()

    processed_parquet     = Channel.value(file("${base}/processed_input.parquet"))
    clusters_parquet      = Channel.value(file("${base}/clusters.parquet"))
    remaining_clusters    = Channel.value(file("${base}/remaining_clusters.parquet"))
    fewer_tax_fa          = Channel.value(file("${base}/clusters_with_fewer_tax_ids.fa"))
    few_taxids_summary    = Channel.value(file("${base}/clusters_with_fewer_taxids_summary.tsv"))
    singletons_fa         = Channel.value(file("${base}/discarded_singletons.fa"))
    singletons_summary    = Channel.value(file("${base}/singleton_cluster_summary.tsv"))
    combined_fasta        = Channel.value(file("${base}/combined_input_fasta.fa"))

  } else {

    /*
     * Full pipeline from scratch
     */
    parsed = PARSE_INPUT_FASTA(
      file(params.input_fasta_dir),
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
    fewer_tax_fa          = filtered.fewer_tax_fa
    few_taxids_summary    = filtered.few_taxids_summary
    singletons_fa         = filtered.singletons_fa
    singletons_summary    = filtered.singletons_summary
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
  queries = Channel
    .fromPath(params.query_species)
    .splitCsv(sep: '\t', header: true)
    .map { row -> tuple(row.tax_id.toString(), row.sps_name.toString()) }
    .distinct { it[0] }

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
   * Merge logs
   */
  merged_logs = MERGE_ALL_LOGS(partial.partial_logs.collect())

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
    .filter { tid1, fas, tid2, name -> tid1 == tid2 }
    .map { tid, fas, _, name ->
      tuple(tid, name.toLowerCase().replaceAll(' ', '_'), fas)
    }
    | MERGE_FASTAS_PER_QUERY

  /*
   * Build *_all_relatives.fa
   */
  all_inputs = merged_fastas.relatives_fa
    .map { fa ->
      def qname = fa.baseName.replaceFirst(/_relatives$/, '')
      tuple(qname, fa)
    }
    .combine(fewer_tax_fa)
    .combine(singletons_fa)
    .map { it ->
      tuple(it[0], it[1], it[2], it[3])
    }

  all_rel = MAKE_ALL_RELATIVES(all_inputs)

  /*
   * Diagnostics & reports
   */

  reports = MAKE_REPORTS(
    clusters_parquet,
    remaining_clusters,
    few_taxids_summary,
    singletons_fa,
    fewer_tax_fa,
    all_rel.all_relatives_fa.collect(),
    combined_fasta
  )

  diagnostics = MAKE_DIAGNOSTICS(
    clusters_parquet,
    singletons_summary,
    remaining_clusters,
    few_taxids_summary,
    merged_logs.merged_log,
    combined_fasta,
    reports.cluster_summary,
    reports.pipeline_summary
  )

  emit:
    clusters_parquet          = clusters_parquet
    remaining_clusters        = remaining_clusters
    ranked_taxa_tsv           = ranked.ranked_taxa_tsv
    closest_relatives_log     = merged_logs.merged_log
    relatives_fastas          = merged_fastas.relatives_fa
    all_relatives_fa          = all_rel.all_relatives_fa
    cluster_summary_txt       = reports.cluster_summary
    echo_pipeline_summary_txt = reports.pipeline_summary
    diagnostics_pdf           = diagnostics.diagnostics_pdf
    diagnostics_summary       = diagnostics.diagnostics_summary
    diagnostics_dir           = diagnostics.diagnostics_dir

}
