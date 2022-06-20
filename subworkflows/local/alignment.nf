//
// Check input samplesheet and get read channels
//

include { BWAMEM2_INDEX } from '../../modules/nf-core/modules/bwamem2/index/main'
include { BWAMEM2_MEM } from '../../modules/nf-core/modules/bwamem2/mem/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/modules/minimap2/align/main'
include { MINIMAP2_INDEX } from '../../modules/nf-core/modules/minimap2/index/main'
include { BAM_TO_PAF } from "../../modules/local/bam_to_paf"
workflow ALIGNMENT {
    take:
    fastq_reads
    fasta

    main:
    fastq_reads.map{ meta, record -> meta.platform }.set{ ch_platform }
    fastq_reads.view()
    ch_platform.view()
    println"____"
    if (ch_platform =~ 'ILLUMINA'){
        println "running minimap2"
        MINIMAP2_ALIGN (
            fastq_reads,
            fasta,
            true,
            true,
            true
        )
    } else {
        println "running bwamem2"
        BWAMEM2_INDEX (
            fasta
        )
        BWAMEM2_MEM(
            fastq_reads,
            BWAMEM2_INDEX.out.index, 
            true
        )
    }



    emit:
    bam  = ch_platform =~ 'ILLUMINA' ? MINIMAP2_ALIGN.out.bam : BWAMEM2_MEM.out.bam // channel: [ val(meta), [ bamfile ] ] ]
    paf  = ch_platform =~ 'ILLUMINA'  ? MINIMAP2_ALIGN.out.paf : '' // channel: [ val(meta), [ bamfile ] ] ]
    versions = ch_platform =~ 'ILLUMINA' ? MINIMAP2_ALIGN.out.versions.first() : BWAMEM2_MEM.out.versions.first()
}