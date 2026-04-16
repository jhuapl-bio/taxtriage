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

        // Resolve paths: use user-supplied value or fall back to bundled project assets
        def proteins_fasta_path = params.annotate_proteins
            ?: "$projectDir/assets/bvbrc_specialty_genes_with_sequences_taxids_and_sites.faa"
        def annotate_meta_path = params.annotate_meta
            ?: "$projectDir/assets/bvbrc_specialty_genes_with_sequences_taxids_and_sites.tsv"

        ch_proteins_fasta = Channel.fromPath(proteins_fasta_path, checkIfExists: true)
            .map { fasta -> [ [id: 'annotate_proteins_db'], fasta ] }

        ch_annotate_meta = Channel.fromPath(annotate_meta_path, checkIfExists: true)

        ANNOTATE_DIAMOND_MAKEDB(
            ch_proteins_fasta,
            [],
            [],
            [],
        )
        ch_versions = ch_versions.mix(ANNOTATE_DIAMOND_MAKEDB.out.versions)

        ch_prep = ch_denovo.combine(
            ANNOTATE_DIAMOND_MAKEDB.out.db.map { _meta, db -> db }
        )

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

        ch_control_empty = ch_prep_branch.control.map { meta, fasta, db ->
            [ meta, file("$projectDir/assets/NO_FILE") ]
        }

        ch_diamond_output = ANNOTATE_DIAMOND_BLASTX.out.txt
            .mix(ch_control_empty)

        ANNOTATE_REPORT(
            ANNOTATE_DIAMOND_BLASTX.out.txt,
            ch_annotate_meta.collect(),
        )
        ch_versions        = ch_versions.mix(ANNOTATE_REPORT.out.versions)
        ch_annotate_report = ANNOTATE_REPORT.out.xlsx

    emit:
        diamond_output  = ch_diamond_output    // [meta, txt]
        annotate_report = ch_annotate_report   // [meta, xlsx]
        versions        = ch_versions
}
