//
// Check input samplesheet and get read channels
//

include { BWA_INDEX } from '../../modules/nf-core/modules/bwa/index/main'
include { BWA_MEM } from '../../modules/nf-core/modules/bwa/mem/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/modules/minimap2/align/main'
include { MINIMAP2_INDEX } from '../../modules/nf-core/modules/minimap2/index/main'
include { BAM_TO_SAM } from "../../modules/local/bam_to_sam"
include { SPLIT_READS } from '../../modules/local/split_reads'
include { SAMTOOLS_MERGE } from '../../modules/nf-core/modules/samtools/merge/main'

workflow ALIGNMENT {
    take:
    fastq_reads
    

    main:
    ch_bams = Channel.empty()
    ch_pafs = Channel.empty()
    
    ch_meta = fastq_reads.first().map{
        meta, fastq, fasta ->
        return  meta
    }.first()


    ch_aligners = fastq_reads
    .branch{
        minimap2: it[0].platform =~ 'OXFORD'
        bwamem2: it[0].platform =~ 'ILLUMINA'
    }
    ch_versions = 1
    ch_aligners.bwamem2
        .transpose(by:[2])
        .map{
            m, fastq, fasta ->
                def basename = fasta.baseName
                def id = basename
                return [ [id:id, platform:m.platform, base: m.id] , fastq, fasta ]
        }
        .set{ transposed_fastas_illumina }
    ch_aligners.minimap2
        .transpose(by:[2])
        .map{
            m, fastq, fasta ->
                def basename = fasta.baseName
                def id = basename
                return [ [id:id, platform:m.platform, base: m.id], fastq, fasta ]
        }
        .set{ transposed_fastas_oxford }

    BWA_INDEX (
        transposed_fastas_illumina.map{ m, fastq, fasta -> 
            return fasta 
        }
    )
    // println "____"
    BWA_MEM(
        transposed_fastas_illumina.map{ m, fastq, fasta -> 
            return [ m, fastq ] 
        },
        BWA_INDEX.out.index,
        true
    )
    ch_bams = ch_bams.mix(BWA_MEM.out.bam)
    MINIMAP2_ALIGN (
        transposed_fastas_oxford.map{ m, fastq, fasta -> [ m, fastq ] },
        transposed_fastas_oxford.map{ meta, fastq, fasta -> fasta },
        true,
        true,
        true
    )
    ch_bams = ch_bams.mix(MINIMAP2_ALIGN.out.bam)

    ch_bams.map{
        m, bams ->
        return bams

    }.collect().set{ collected_bams }
    collected_bams = ch_meta.combine(collected_bams.map{
        m -> 
        return [m]
    })
    
    

    SAMTOOLS_MERGE(
        collected_bams
    )

    BAM_TO_SAM(
        SAMTOOLS_MERGE.out.bam
    )
    ch_sams=BAM_TO_SAM.out.sam
    ch_pileups=BAM_TO_SAM.out.mpileup

    

    



    emit:
    sam  = ch_sams // channel: [ val(meta), [ paffile ] ] ]
    bams  = SAMTOOLS_MERGE.out.bam // channel: bamfile
    mpileup = ch_pileups
    // // stats = ch_stats
    versions = ch_versions
}