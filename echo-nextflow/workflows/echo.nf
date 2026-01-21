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
   * P1 – parse input FASTAs + metadata
   */
  parsed = PARSE_INPUT_FASTA(
    file(params.input_fasta_dir),
    file(params.metadata_tsv)
  )

  /*
   * P2 – MMseqs clustering
   */
  mm = RUN_MMSEQS(parsed.combined_fasta)

  /*
   * P3 – parse MMseqs clusters
   */
  cl = PARSE_CLUSTERS(
    mm.cluster_tsv,
    parsed.processed_parquet
  )

  /*
   * P4 – filter clusters (singletons / fewer tax IDs)
   */
  filtered = PROCESS_CLUSTERS(cl.clusters_parquet)

  /*
   * P5 – rank taxa (query × input species distances)
   */
  ranked = RANK_TAXA(
    parsed.processed_parquet,
    file(params.query_species)
  )

  /*
   * P6 – split remaining clusters into chunk parquet files
   */
  chunks = SPLIT_REMAINING_CLUSTERS(filtered.remaining_clusters)
  chunk_files = chunks.chunks.flatten()

  /*
   * Load queries: (tax_id, species_name)
   */
  queries = Channel
    .fromPath(params.query_species)
    .splitCsv(sep: '\t', header: true)
    .map { row -> tuple(row.tax_id.toString(), row.sps_name.toString()) }

  /*
   * Build (chunk × query) inputs
   */

  closest_inputs = chunk_files
    .combine(queries)
    .combine(ranked.ranked_taxa_tsv)     // (chunk_file, tax_id, query_name, ranked_taxa_tsv)
    .map { it ->
        tuple(it[0], it[3], it[1], it[2])
    }

  // chunk_files.combine(queries).view { "PAIR=" + it }

  /*
   * P7 – closest relatives per (chunk × query)
   */
  partial = CLOSEST_RELATIVES_CHUNK(closest_inputs)

  /*
   * Merge all chunk logs → single TSV
   */
  merged_logs = MERGE_ALL_LOGS(partial.partial_logs.collect())

  /*
   * P8 – merge chunk FASTAs per query
   */
  fastas_grouped = partial.partial_fastas
    .map { fa ->
      def tax_id = fa.baseName.tokenize('_')[0]   // 78070_chunk_003.fa → 78070
      tuple(tax_id, fa)
    }
    .groupTuple()

  merge_inputs = fastas_grouped
    .combine(queries)
    .filter { tid1, fas, tid2, name -> tid1 == tid2 }
    .map { tid, fas, _, name ->
        tuple(tid, name.toLowerCase().replaceAll(' ', '_'), fas)
    }

  merged_fastas = MERGE_FASTAS_PER_QUERY(merge_inputs)

  /*
   * P9 – build *_all_relatives.fa
   * always includes clusters_with_fewer_tax_ids
   * optionally includes discarded_singletons
   */
  
  all_inputs = merged_fastas.relatives_fa
    .map { fa ->
      def qname = fa.baseName.replaceFirst(/_relatives$/, '')
      tuple(qname, fa)
    }
    .combine(filtered.fewer_tax_fa)
    .combine(filtered.singletons_fa)
    .map { it ->
          // it = [ qname, relatives_fa, fewer_tax_fa, singletons_fa ]
         tuple(it[0], it[1], it[2], it[3])
    }
  
  all_rel = MAKE_ALL_RELATIVES(all_inputs)

  reports = MAKE_REPORTS(
  cl.clusters_parquet,
  filtered.remaining_clusters,
  filtered.few_taxids_summary,
  filtered.singletons_fa,
  filtered.fewer_tax_fa,
  all_rel.all_relatives_fa.collect(),
  parsed.combined_fasta
  )
  
  diagn = MAKE_DIAGNOSTICS(
  cl.clusters_parquet,
  filtered.singletons_summary,
  filtered.remaining_clusters,
  filtered.few_taxids_summary,
  merged_logs.merged_log,
  parsed.combined_fasta,
  reports.cluster_summary,
  reports.pipeline_summary
  )

  /*
   * Final outputs
   */
  emit:
    clusters_parquet         = cl.clusters_parquet
    remaining_clusters       = filtered.remaining_clusters
    ranked_taxa_tsv          = ranked.ranked_taxa_tsv
    closest_relatives_log    = merged_logs.merged_log
    relatives_fastas         = merged_fastas.relatives_fa
    all_relatives_fastas     = all_rel.all_relatives_fa
    cluster_summary_txt       = reports.cluster_summary
    echo_pipeline_summary_txt = reports.pipeline_summary
    diagnostics_pdf     = diagn.diagnostics_pdf
    diagnostics_summary  = diagn.diagnostics_summary
    diagnostics_dir      = diagn.diagnostics_dir

}

