//
// Check input samplesheet and get read channels
//

include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/modules/kraken2/kraken2/main'
include { MOVE_FILES } from '../../modules/local/moveFiles.nf'

workflow FILTER_READS {
    take:
        fastq_reads
        db

    main:
        KRAKEN2_KRAKEN2(
            fastq_reads,
            db,
            true,
            true
        )
        MOVE_FILES(
            KRAKEN2_KRAKEN2.out.unclassified_reads_fastq,
            "filtered_",
            false,
            []
        )
        ch_unclassified_reads = MOVE_FILES.out.files


    emit:
        reads = ch_unclassified_reads
        versions = KRAKEN2_KRAKEN2.out.versions
}