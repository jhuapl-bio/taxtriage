//
// Annotate de novo assemblies with protein-level DIAMOND BLASTx
//
include { DIAMOND_MAKEDB } from '../../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTX } from '../../modules/nf-core/diamond/blastx/main'

workflow PROTEINS {
    take:
        ch_denovo  // [meta, fasta]  — de novo assembly per sample

    main:
        ch_versions      = Channel.empty()
        ch_diamond_output = Channel.empty()

        if (params.annotate_proteins) {
            // Build a single DIAMOND protein DB from the reference fasta
            ch_proteins_fasta = Channel.fromPath(params.annotate_proteins, checkIfExists: true)
                .map { fasta -> [ [id: 'annotate_proteins_db'], fasta ] }

            DIAMOND_MAKEDB(
                ch_proteins_fasta,
                [],
                [],
                [],
            )
            ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions)

            // Pair every assembly with the single DB file (drop the DB meta)
            ch_prep = ch_denovo.combine(
                DIAMOND_MAKEDB.out.db.map { _meta, db -> db }
            )
            // ch_prep: [meta, fasta, db]

            DIAMOND_BLASTX(
                ch_prep.map { meta, fasta, db -> [meta, fasta] },
                ch_prep.map { meta, fasta, db -> [[id: 'annotate_proteins_db'], db] },
                'txt',
                false,
            )
            ch_versions      = ch_versions.mix(DIAMOND_BLASTX.out.versions)
            ch_diamond_output = DIAMOND_BLASTX.out.txt
        }

    emit:
        diamond_output = ch_diamond_output   // [meta, txt]
        versions       = ch_versions
}
