//
// Annotate de novo assemblies with protein-level DIAMOND BLASTx
//
include { DIAMOND_MAKEDB as ANNOTATE_DIAMOND_MAKEDB } from '../../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTX as ANNOTATE_DIAMOND_BLASTX } from '../../modules/nf-core/diamond/blastx/main'
include { ANNOTATE_REPORT                           } from '../../modules/local/annotate_report'

workflow PROTEINS {
    take:
        ch_denovo  // [meta, fasta]  -- de novo assembly per sample

    main:
        ch_versions        = Channel.empty()
        ch_diamond_output  = Channel.empty()
        ch_annotate_report = Channel.empty()
        ch_proteins_fasta = params.annotate_proteins ? Channel.fromPath(params.annotate_proteins, checkIfExists: true) : Channel.fromPath("$projectDir/assets/bvbrc_specialty_genes_with_sequences_taxids_and_sites.faa", checkIfExists: true)
            .map { fasta -> [ [id: 'annotate_proteins_db'], fasta ] }
        ch_annotate_meta = params.annotate_meta ? Channel.fromPath(params.annotate_meta, checkIfExists: true) : Channel.empty()


        if (params.annotate_proteins) {
            // Build a single DIAMOND protein DB from the reference fasta
            ANNOTATE_DIAMOND_MAKEDB(
                ch_proteins_fasta,
                [],
                [],
                [],
            )
            ch_versions = ch_versions.mix(ANNOTATE_DIAMOND_MAKEDB.out.versions)

            // Pair every assembly with the single DB file (drop the DB meta)
            ch_prep = ch_denovo.combine(
                ANNOTATE_DIAMOND_MAKEDB.out.db.map { _meta, db -> db }
            )

            // Branch: only run DIAMOND on non-control samples
            ch_prep_branch = ch_prep.branch {
                control:    it[0].control == true
                noncontrol: true
            }

            ANNOTATE_DIAMOND_BLASTX(
                ch_prep_branch.noncontrol.map { meta, fasta, db -> [meta, fasta] },
                ch_prep_branch.noncontrol.map { meta, fasta, db -> [[id: 'annotate_proteins_db'], db] },
                'txt',
                false,
            )
            ch_versions = ch_versions.mix(ANNOTATE_DIAMOND_BLASTX.out.versions)

            // For controls, emit an empty placeholder file
            ch_control_empty = ch_prep_branch.control.map { meta, fasta, db ->
                [ meta, file("$projectDir/assets/NO_FILE") ]
            }
            ch_diamond_output = ANNOTATE_DIAMOND_BLASTX.out.txt
                .mix(ch_control_empty)

            // ── Annotate matched hits against metadata and produce xlsx report ──
            if (params.annotate_meta) {

                ANNOTATE_REPORT(
                    ANNOTATE_DIAMOND_BLASTX.out.txt,
                    ch_annotate_meta.collect(),
                )
                ch_versions        = ch_versions.mix(ANNOTATE_REPORT.out.versions)
                ch_annotate_report = ANNOTATE_REPORT.out.xlsx
            }
        }

    emit:
        diamond_output  = ch_diamond_output    // [meta, txt]
        annotate_report = ch_annotate_report   // [meta, xlsx]
        versions        = ch_versions
}
