//
// NOVELTY subworkflow  (contig-first path, pluggable backend)
//
// Reference-free / open-set novelty branch. Classifies the de novo assembly built from ALL
// post-QC reads (aligned or not) -- the same contigs handed to the annotation branch -- rather
// than only the unmapped residual. This lets the search report organisms that DID align to a
// pulled reference (e.g. targeted Ebola) alongside anything novel.
//
// The backend is chosen by --novelty (params.novelty). By default ALL backends consume the de novo
// contigs DIRECTLY. Optionally (--novelty_gene) a Pyrodigal ORF-prediction step is inserted first
// and the predicted nucleotide CDS (.fna) become the query instead of the raw contigs. This trades
// a little runtime for finer granularity: one whole-genome contig becomes one query per gene, so a
// single near-complete genome contributes several LCA rows (counts) instead of one -- closer to a
// per-gene BLAST view and friendlier to the contig-count-based novelty score.
//
//   --novelty mmseqs2 : contigs (or genes) -> MMSEQS_TAXONOMY (translated LCA; mmseqs extracts ORFs
//                       internally, so nucleotide input is fine -- predicted genes are just cleaner
//                       and faster than whole contigs)
//   --novelty kaiju   : contigs (or genes) -> KAIJU (greedy translated search; self-translates)
//   --novelty bracken : contigs (or genes) -> KRAKEN2_NOVELTY -> BRACKEN (count-weighted abundance)
//
// --novelty_gene OFF (default) : query = de novo contigs (current behaviour, unchanged)
// --novelty_gene ON            : query = PYRODIGAL predicted genes (.fna) over those same contigs
//
// Every backend emits the SAME contract for the score: a per-query (or count-weighted) LCA tsv,
// an optional best-hit table with %identity, and a kraken-style report.
//
// Shared-cluster path (--microbert): the MicrobeRT MMSEQS_EASYCLUSTER reps are used as the query
// instead of the raw contigs. They are still nucleotide fasta, so the same backends apply.
//
// EXTRACT_UNMAPPED produces only the read-accounting counts that NOVELTY_SCORE needs as
// denominators -- it does not extract a read stream for taxonomy.
//

include { EXTRACT_UNMAPPED } from '../../modules/local/extract_unmapped'
include { PYRODIGAL        } from '../../modules/local/pyrodigal'
include { MMSEQS_TAXONOMY  } from '../../modules/local/mmseqs_taxonomy'
include { KAIJU            } from '../../modules/local/kaiju'
include { KRAKEN2_KRAKEN2 as KRAKEN2_NOVELTY } from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN          } from '../../modules/local/bracken'
include { NOVELTY_SCORE    } from '../../modules/local/novelty_score'

workflow NOVELTY {
    take:
        ch_bams            // [meta, bam, csi]   ALIGNMENT.out.bams (reference-merged per sample)
        ch_kraken2_report  // [meta, kreport]    CLASSIFIER.out.ch_kraken2_report (or NO_FILE)
        ch_denovo          // [meta, contigs]    ASSEMBLY.out.ch_denovo_assembly (may be empty)
        ch_novelty_db      // staged novelty db DIRECTORY (resolved local-or-downloaded by caller)
        ch_precluster_reps // [meta, rep_seq.fasta] shared MicrobeRT cluster reps (empty unless --microbert)

    main:
        ch_versions = Channel.empty()
        ch_empty    = file("$projectDir/assets/NO_FILE")
        def method  = params.novelty   // 'mmseqs2' | 'kaiju' | 'bracken'

        // Skip controls -- a flagged negative control is a process signal, not a sample call.
        ch_bams_nc = ch_bams.filter { meta, bam, csi -> !(meta.control == true) }

        // Attach the kreport. INNER join: keeps only samples with both a BAM and a kreport.
        ch_extract_in = ch_bams_nc
            .join(ch_kraken2_report)
            .map { meta, bam, csi, kr -> [meta, bam, csi, kr] }

        // Read accounting only (no read extraction): supplies NOVELTY_SCORE's denominators.
        EXTRACT_UNMAPPED(ch_extract_in)
        ch_versions = ch_versions.mix(EXTRACT_UNMAPPED.out.versions)

        // Fold the counts onto meta so NOVELTY_SCORE has its denominators.
        ch_meta_counts = EXTRACT_UNMAPPED.out.counts.map { meta, total, mapped, unmapped, k2 ->
            def m = meta + [ total_reads   : (total    as Integer),
                             ref_aligned   : (mapped   as Integer),
                             ref_unaligned : (unmapped as Integer),
                             k2_classified : (k2       as Integer) ]
            [m.id, m]
        }

        // ---------------------------------------------------------------------------------
        // Query source. Two mutually-exclusive paths, each with a SINGLE map so the DSL2
        // graph stays valid:
        //   --microbert ON  : reuse the shared MicrobeRT cluster reps.
        //   --microbert OFF : the full de novo contigs (all post-QC reads, aligned or not).
        // The folded-counts meta is attached so it rides through to NOVELTY_SCORE.
        // ---------------------------------------------------------------------------------
        if (params.microbert) {
            ch_query = ch_precluster_reps
                .map { meta, reps -> [meta.id, reps] }
                .join(ch_meta_counts)                                      // [id, reps, meta]
                .map { id, reps, meta -> [meta, reps] }
        } else {
            ch_query = ch_denovo
                .map { meta, c -> [meta.id, c] }
                .join(ch_meta_counts)                                      // [id, contigs, meta]
                .filter { id, c, meta -> c != null && c.name != 'NO_FILE' } // require real contigs
                .map { id, c, meta -> [meta, c] }
        }

        // ---------------------------------------------------------------------------------
        // Optional gene-prediction step (--novelty_gene). When on, Pyrodigal predicts ORFs on
        // the selected query (contigs or MicrobeRT reps) and the predicted nucleotide CDS (.fna)
        // replace the contigs as the per-query unit for every backend. meta is preserved so the
        // folded read-accounting counts still ride through to NOVELTY_SCORE.
        // ---------------------------------------------------------------------------------
        if (params.novelty_gene) {
            PYRODIGAL(ch_query)
            ch_versions = ch_versions.mix(PYRODIGAL.out.versions)
            ch_query = PYRODIGAL.out.genes            // [meta, *.fna] predicted nucleotide CDS
        }

        // ---------------------------------------------------------------------------------
        // Backend dispatch. Each branch ends by setting ch_lca / ch_tophit / ch_report so the
        // scoring + emit below has a single call site per process.
        // ---------------------------------------------------------------------------------
        ch_lca    = Channel.empty()
        ch_tophit = Channel.empty()
        ch_report = Channel.empty()

        if (method == 'kaiju') {
            KAIJU(ch_query, ch_novelty_db)
            ch_versions = ch_versions.mix(KAIJU.out.versions)
            ch_lca    = KAIJU.out.lca
            ch_tophit = KAIJU.out.tophit
            ch_report = KAIJU.out.report
        }
        else if (method == 'bracken') {
            // Kraken2 over the contigs first, then bracken re-estimates abundance.
            // Contigs are always a single FASTA -- force single_end:true so the module
            // does not add --paired (which requires two files and would cause Kraken2 to fail).
            ch_query_se = ch_query.map { meta, contigs -> [ meta + [single_end: true], contigs ] }
            KRAKEN2_NOVELTY(ch_query_se, ch_novelty_db, false, false)
            ch_versions = ch_versions.mix(KRAKEN2_NOVELTY.out.versions)
            BRACKEN(KRAKEN2_NOVELTY.out.report, ch_novelty_db)
            ch_versions = ch_versions.mix(BRACKEN.out.versions)
            ch_lca    = BRACKEN.out.lca
            ch_tophit = BRACKEN.out.tophit
            ch_report = BRACKEN.out.report
        }
        else {
            // default / 'mmseqs2': translated LCA straight on the contigs.
            MMSEQS_TAXONOMY(ch_query, ch_novelty_db)
            ch_versions = ch_versions.mix(MMSEQS_TAXONOMY.out.versions)
            ch_lca    = MMSEQS_TAXONOMY.out.lca
            ch_tophit = MMSEQS_TAXONOMY.out.tophit
            ch_report = MMSEQS_TAXONOMY.out.report
        }

        // Single-pass scoring. The run_summaries input is already plumbed for a future
        // two-pass across-sample z-score mode. NOVELTY_SCORE reads params.novelty to decide
        // whether the LCA carries a count column (bracken) or one row per query (mmseqs/kaiju).
        ch_score_in = ch_lca.join(ch_tophit)
        NOVELTY_SCORE(ch_score_in, ch_empty)
        ch_versions = ch_versions.mix(NOVELTY_SCORE.out.versions)

    emit:
        summary    = NOVELTY_SCORE.out.summary      // [meta, *.novelty.summary.tsv]
        candidates = NOVELTY_SCORE.out.candidates   // [meta, *.novelty.candidates.tsv]
        report     = ch_report                      // kraken-style, reusable by krona/kreport
        versions   = ch_versions
}
