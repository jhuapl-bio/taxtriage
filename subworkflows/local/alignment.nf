//
// Check input samplesheet and get read channels
//

include { BWA_INDEX } from '../../modules/nf-core/modules/bwa/index/main'
include { BWA_MEM } from '../../modules/nf-core/modules/bwa/mem/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/modules/minimap2/align/main'
include { MINIMAP2_INDEX } from '../../modules/nf-core/modules/minimap2/index/main'
include { BAM_TO_SAM } from "../../modules/local/bam_to_sam"
include { RSEQC_BAMSTAT } from '../../modules/nf-core/modules/rseqc/bamstat/main'

workflow ALIGNMENT {
    take:
    fastq_reads
    

    main:
    ch_bams = Channel.empty()
    ch_pafs = Channel.empty()

    ch_aligners = fastq_reads
    .branch{
        minimap2: it[0].platform =~ 'OXFORD'
            
        bwamem2: it[0].platform =~ 'ILLUMINA'
    }
    ch_versions = 1
    ch_aligners.bwamem2.view()
    println ch_aligners.bwamem2
    BWA_INDEX (
        ch_aligners.bwamem2.map{ meta, fastq, fasta -> fasta }
    )
    fastq_reads.view()
    BWA_INDEX.out.index.view()
    BWA_MEM(
        ch_aligners.bwamem2.map{ meta, fastq, fasta -> [ meta, fastq ] },
        BWA_INDEX.out.index,
        true
    )
    ch_bams = ch_bams.mix(BWA_MEM.out.bam)
    
    MINIMAP2_ALIGN (
        ch_aligners.minimap2.map{ meta, fastq, fasta -> [ meta, fastq ] },
        ch_aligners.minimap2.map{ meta, fastq, fasta -> fasta },
        true,
        true,
        true
    )
    ch_bams = ch_bams.mix(MINIMAP2_ALIGN.out.bam)
    BAM_TO_SAM(
        ch_bams
    )
    ch_sams=BAM_TO_SAM.out.sam




    // RSEQC_BAMSTAT(
    //     ch_bams
    // )    
    // RSEQC_BAMSTAT.out.txt.set{ ch_stats }

    emit:
    sam  = ch_sams // channel: [ val(meta), [ paffile ] ] ]
    bams  = ch_bams // channel: bamfile
    // stats = ch_stats
    versions = ch_versions
}