//
// NOVELTY subworkflow  (ORF-first path)
//
// Reference-free / open-set novelty branch. Consumes the closed-set residual: reads that
// aligned to NONE of the pulled references, optionally pooled with de novo contigs.
//
// Shared-cluster path (when --microbert is set):
//   MicrobeRT MMSEQS_EASYCLUSTER reps (passed in) → PYRODIGAL → MMSEQS_TAXONOMY
//   The clustering is computed once upstream and reused by both the MicrobeRT classifier
//   and this branch, so de novo contigs / MMSEQS_EASYCLUSTER_NOVELTY are bypassed here.
//
// Preferred path (no --microbert, contigs available via --use_denovo):
//   contigs → PYRODIGAL (metagenome ORF prediction)
//           → MMSEQS_EASYCLUSTER_NOVELTY (protein deduplication, ~90 % id)
//           → MMSEQS_TAXONOMY (protein-protein LCA against seqTaxDB)
//
// This is orders-of-magnitude faster than feeding raw reads to mmseqs taxonomy because:
//   1. Contigs are longer → Pyrodigal predicts actual ORFs rather than 6-frame guessing
//   2. Clustering removes near-identical proteins before the expensive DB search
//   3. The search is protein-protein rather than nucleotide → internal 6-frame per read
//
// Fallback path (no contigs — --use_denovo not set):
//   unmapped reads → MMSEQS_TAXONOMY (nucleotide translated-search, existing behaviour)
//
// Self-contained: it extracts the unmapped reads and the read accounting itself, so the only
// thing the caller has to provide is the standard ALIGNMENT/CLASSIFIER outputs + the DB.
//

include { EXTRACT_UNMAPPED                          } from '../../modules/local/extract_unmapped'
include { PYRODIGAL                                 } from '../../modules/local/pyrodigal'
include { MMSEQS_EASYCLUSTER as MMSEQS_EASYCLUSTER_NOVELTY } from '../../modules/local/mmseqs2_easycluster'
include { MMSEQS_TAXONOMY                           } from '../../modules/local/mmseqs_taxonomy'
include { NOVELTY_SCORE                             } from '../../modules/local/novelty_score'

workflow NOVELTY {
    take:
        ch_bams            // [meta, bam, csi]   ALIGNMENT.out.bams (reference-merged per sample)
        ch_kraken2_report  // [meta, kreport]    CLASSIFIER.out.ch_kraken2_report (or NO_FILE)
        ch_denovo          // [meta, contigs]    ASSEMBLY.out.ch_denovo_assembly (may be empty)
        ch_seqtaxdb        // staged seqTaxDB directory (resolved local-or-downloaded by caller)
        ch_precluster_reps // [meta, rep_seq.fasta] shared MicrobeRT cluster reps (empty unless --microbert)

    main:
        ch_versions = Channel.empty()
        ch_empty    = file("$projectDir/assets/NO_FILE")

        // Skip controls -- a flagged negative control is a process signal, not a sample call.
        ch_bams_nc = ch_bams.filter { meta, bam, csi -> !(meta.control == true) }

        // Attach the kreport. INNER join: keeps only samples with both a BAM and a kreport.
        ch_extract_in = ch_bams_nc
            .join(ch_kraken2_report)
            .map { meta, bam, csi, kr -> [meta, bam, csi, kr] }

        EXTRACT_UNMAPPED(ch_extract_in)
        ch_versions = ch_versions.mix(EXTRACT_UNMAPPED.out.versions)

        // Fold the counts onto meta so NOVELTY_SCORE has its denominators.
        ch_meta_counts = EXTRACT_UNMAPPED.out.counts.map { meta, total, mapped, k2 ->
            def m = meta + [ total_reads   : (total as Integer),
                             ref_aligned   : (mapped as Integer),
                             k2_classified : (k2 as Integer) ]
            [m.id, m]
        }

        // ---------------------------------------------------------------------------------
        // Choose the Pyrodigal input depending on whether the shared MicrobeRT cluster is
        // available (--microbert). Two mutually-exclusive paths, each with a SINGLE call site
        // per process so the DSL2 graph stays valid:
        //
        //   --microbert ON  (shared-cluster path):
        //       MMSEQS_EASYCLUSTER reps (nucleotide) -> PYRODIGAL -> MMSEQS_TAXONOMY
        //       The clustering already happened upstream and is reused here, so the de novo
        //       contigs and MMSEQS_EASYCLUSTER_NOVELTY are bypassed entirely.
        //
        //   --microbert OFF (fallback path, unchanged):
        //       de novo contigs -> PYRODIGAL -> MMSEQS_EASYCLUSTER_NOVELTY -> MMSEQS_TAXONOMY,
        //       with a reads-only fallback for samples that produced no contigs.
        // ---------------------------------------------------------------------------------
        ch_reads_query = Channel.empty()
        if (params.microbert) {
            // Attach the count-enriched meta (controls drop out via the inner join with
            // ch_meta_counts, which is derived from the non-control BAMs).
            ch_pyrodigal_in = ch_precluster_reps
                .map { meta, reps -> [meta.id, reps] }
                .join(ch_meta_counts)                                      // [id, reps, meta]
                .map { id, reps, meta -> [meta, reps] }
        } else {
            ch_reads_keyed  = EXTRACT_UNMAPPED.out.reads.map { meta, r -> [meta.id, r] }
            ch_denovo_keyed = ch_denovo.map                  { meta, c -> [meta.id, c] }

            ch_joined = ch_reads_keyed
                .join(ch_meta_counts)                                      // [id, reads, meta]
                .join(ch_denovo_keyed, remainder: true)                   // [id, reads, meta, contigs?]
                .filter { row -> row[1] != null && row[2] != null }       // require reads + meta

            ch_joined.branch { row ->
                with_contigs : row.size() > 3 && row[3] != null && row[3].name != 'NO_FILE'
                reads_only   : true
            }.set { ch_branched }

            ch_pyrodigal_in = ch_branched.with_contigs.map { row -> [ row[2], row[3] ] }
            // Reads-only fallback (no contigs): translated-search the raw reads directly.
            ch_reads_query  = ch_branched.reads_only.map  { row -> [ row[2], [row[1]] ] }
        }

        // Predict ORFs in metagenome mode; outputs .faa protein sequences. Single call site.
        PYRODIGAL(ch_pyrodigal_in)
        ch_versions = ch_versions.mix(PYRODIGAL.out.versions)

        if (params.microbert) {
            // Reps were already clustered upstream, so go straight to taxonomy.
            ch_taxonomy_input = PYRODIGAL.out.proteins
        } else {
            // Cluster predicted proteins to remove redundancy before the DB search.
            // Aliased MMSEQS_EASYCLUSTER so modules.config can give it protein identity
            // thresholds without touching the annotation/MicrobeRT clustering settings.
            MMSEQS_EASYCLUSTER_NOVELTY(PYRODIGAL.out.proteins)
            ch_versions = ch_versions.mix(MMSEQS_EASYCLUSTER_NOVELTY.out.versions)

            // Merge ORF reps with the reads-only fallback before taxonomy.
            ch_taxonomy_input = MMSEQS_EASYCLUSTER_NOVELTY.out.representatives.mix(ch_reads_query)
        }

        MMSEQS_TAXONOMY(ch_taxonomy_input, ch_seqtaxdb)
        ch_versions = ch_versions.mix(MMSEQS_TAXONOMY.out.versions)

        // Single-pass scoring. The run_summaries input is already plumbed for a future
        // two-pass across-sample z-score mode.
        ch_score_in = MMSEQS_TAXONOMY.out.lca.join(MMSEQS_TAXONOMY.out.tophit)

        NOVELTY_SCORE(ch_score_in, ch_empty)
        ch_versions = ch_versions.mix(NOVELTY_SCORE.out.versions)

    emit:
        summary    = NOVELTY_SCORE.out.summary      // [meta, *.novelty.summary.tsv]
        candidates = NOVELTY_SCORE.out.candidates   // [meta, *.novelty.candidates.tsv]
        report     = MMSEQS_TAXONOMY.out.report     // kraken-style, reusable by krona/kreport
        versions   = ch_versions
}
