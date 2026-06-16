//
// NOVELTY subworkflow  (fully wired)
//
// Reference-free / open-set novelty branch. Consumes the closed-set residual: reads that
// aligned to NONE of the pulled references, optionally pooled with de novo contigs. Runs a
// translated-search LCA (protein space survives genus/family divergence) and turns it into
// (a) a per-sample novelty signal and (b) genus-or-higher candidate calls.
//
// Self-contained: it extracts the unmapped reads and the read accounting itself, so the only
// thing the caller has to provide is the standard ALIGNMENT/CLASSIFIER outputs + the DB.
//
include { EXTRACT_UNMAPPED } from '../../modules/local/extract_unmapped'
include { MMSEQS_TAXONOMY  } from '../../modules/local/mmseqs_taxonomy'
include { NOVELTY_SCORE    } from '../../modules/local/novelty_score'

workflow NOVELTY {
    take:
        ch_bams            // [meta, bam, csi]   ALIGNMENT.out.bams (reference-merged per sample)
        ch_kraken2_report  // [meta, kreport]    CLASSIFIER.out.ch_kraken2_report (or NO_FILE)
        ch_denovo          // [meta, contigs]    ASSEMBLY.out.ch_denovo_assembly (may be empty)
        ch_seqtaxdb        // path to prebuilt mmseqs seqTaxDB (params.novelty_db)

    main:
        ch_versions = Channel.empty()
        ch_empty    = file("$projectDir/assets/NO_FILE")

        // Skip controls -- a flagged negative control is a process signal, not a sample call.
        ch_bams_nc = ch_bams.filter { meta, bam, csi -> !(meta.control == true) }

        // Attach the kreport (left join so samples without one still proceed with NO_FILE).
        ch_extract_in = ch_bams_nc
            .join(ch_kraken2_report, remainder: true)
            .map { meta, bam, csi, kr -> [meta, bam, csi, kr ?: ch_empty] }
            .filter { meta, bam, csi, kr -> bam }   // drop the kreport-only remainder rows

        EXTRACT_UNMAPPED(ch_extract_in)
        ch_versions = ch_versions.mix(EXTRACT_UNMAPPED.out.versions)

        // Fold the counts (env strings) onto meta so NOVELTY_SCORE has its denominators.
        ch_meta_counts = EXTRACT_UNMAPPED.out.counts.map { meta, total, mapped, k2 ->
            def m = meta + [ total_reads   : (total as Integer),
                             ref_aligned   : (mapped as Integer),
                             k2_classified : (k2 as Integer) ]
            [m.id, m]
        }

        // Query = unmapped reads (+ contigs when assembly ran). Everything is keyed by id
        // (string) for joins, because the enriched meta map is not a stable join key.
        ch_query = EXTRACT_UNMAPPED.out.reads.map { meta, reads -> [meta.id, reads] }
            .join(ch_meta_counts)                                   // [id, reads, meta+counts]
            .join(ch_denovo.map { meta, c -> [meta.id, c] }, remainder: true)  // [id, reads, m, contigs]
            .map { id, reads, m, contigs ->
                def parts = [reads, contigs].flatten().findAll { it && it.name != 'NO_FILE' }
                [m, parts]
            }
            .filter { m, parts -> m && parts }

        MMSEQS_TAXONOMY(ch_query, ch_seqtaxdb)
        ch_versions = ch_versions.mix(MMSEQS_TAXONOMY.out.versions)

        // Single-pass scoring (in-sample fallback baseline). For a true across-sample z-score,
        // see the two-pass note in README -- the run_summaries input is already plumbed.
        ch_score_in = MMSEQS_TAXONOMY.out.lca.join(MMSEQS_TAXONOMY.out.tophit)

        NOVELTY_SCORE(ch_score_in, ch_empty)
        ch_versions = ch_versions.mix(NOVELTY_SCORE.out.versions)

    emit:
        summary    = NOVELTY_SCORE.out.summary      // [meta, *.novelty.summary.tsv]
        candidates = NOVELTY_SCORE.out.candidates   // [meta, *.novelty.candidates.tsv]
        report     = MMSEQS_TAXONOMY.out.report     // kraken-style, reusable by krona/kreport tooling
        versions   = ch_versions
}
