//
// Check input samplesheet and get read channels
//

include { BWA_INDEX } from '../../modules/nf-core/modules/bwa/index/main'
include { BWA_MEM } from '../../modules/nf-core/modules/bwa/mem/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/modules/minimap2/align/main'
include { MINIMAP2_INDEX } from '../../modules/nf-core/modules/minimap2/index/main'
include { BAM_TO_SAM } from "../../modules/local/bam_to_sam"
include { SAMTOOLS_MERGE } from '../../modules/nf-core/modules/samtools/merge/main'
include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/modules/bcftools/consensus/main'
include { BCFTOOLS_MPILEUP as BCFTOOLS_MPILEUP_ILLUMINA } from '../../modules/nf-core/modules/bcftools/mpileup/main'
include { BCFTOOLS_MPILEUP as BCFTOOLS_MPILEUP_OXFORD } from '../../modules/nf-core/modules/bcftools/mpileup/main'
include { BCFTOOLS_STATS } from '../../modules/nf-core/modules/bcftools/stats/main'
include { MERGE_FASTA } from '../../modules//local/merge_fasta'
include { SPLIT_READS } from '../../modules/local/split_reads'


workflow ALIGNMENT {
    take:
    fastq_reads
    

    main:
    ch_bams = Channel.empty()
    ch_fasta = Channel.empty()
    ch_pileups = Channel.empty()
    ch_stats = Channel.empty()
    ch_pafs = Channel.empty()
    ch_sams = Channel.empty()
    ch_aligners = Channel.empty()
    collected_bams  = Channel.empty()
    
    ch_merged_fasta = Channel.empty()

    
    fastq_reads
    .branch{
        minimap2: it[0].platform =~ 'OXFORD'
        bwamem2: it[0].platform =~ 'ILLUMINA'
    }.set { ch_aligners }

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
                return [ [id:id, single_end: true, platform:m.platform, base: m.id], fastq, fasta ]
        }
        .set{ transposed_fastas_oxford }

    BWA_INDEX (
        transposed_fastas_illumina.map{ m, fastq, fasta -> 
            return fasta 
        }
    )
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
    ch_bams = ch_bams.mix(BWA_MEM.out.bam)
    ch_bams.map{
        m, bams ->
        return [m.base,m, bams]

    }.groupTuple(by:0).map{
        base,meta,fasta -> [meta.first(), fasta]
    }.set{ collected_bams }
    fastq_reads.map{
        m, reads, fasta ->
            [ m.id, m, reads ]
        
    }
    .join(
    collected_bams
    .map{
        m,bams -> 
        return [ m.base, bams ]
    }
    , by:0)
    .map{
        id, m, fastqs, bams -> [ m, bams ]
    }.set { collected_bams }

    if (!params.skip_variants){
        BCFTOOLS_MPILEUP_OXFORD(
            MINIMAP2_ALIGN.out.bam, 
            transposed_fastas_oxford.map{m,fastq,fasta -> fasta },
            false
        )
        ch_merged_mpileup_oxford = BCFTOOLS_MPILEUP_OXFORD.out.vcf.join(BCFTOOLS_MPILEUP_OXFORD.out.tbi)
        ch_merged_mpileup_oxford = ch_merged_mpileup_oxford.join(transposed_fastas_oxford.map{m,fastq,fasta -> [m,fasta] })
        
        BCFTOOLS_MPILEUP_ILLUMINA(
            BWA_MEM.out.bam, 
            transposed_fastas_illumina.map{m,fastq,fasta -> fasta },
            false
        )    
        ch_merged_mpileup_illumina = BCFTOOLS_MPILEUP_ILLUMINA.out.vcf.join(BCFTOOLS_MPILEUP_ILLUMINA.out.tbi)
        ch_merged_mpileup_illumina= ch_merged_mpileup_illumina.join(transposed_fastas_illumina.map{m,fastq,fasta -> [m,fasta] })

        ch_merged_mpileup = ch_merged_mpileup_illumina.mix(ch_merged_mpileup_oxford)
        BCFTOOLS_STATS(
            ch_merged_mpileup.map{
                m, vcf, tbi, fasta -> [m, vcf, tbi]
                
            },
            [],
            [],
            []
        )
        if (!params.skip_consensus){
            BCFTOOLS_CONSENSUS (
                ch_merged_mpileup
            )
            ch_fasta = BCFTOOLS_CONSENSUS.out.fasta
            BCFTOOLS_CONSENSUS.out.fasta.map{
                m,fasta-> [m.base,m,fasta]
            }.groupTuple(by:0).map{
                base,meta,fasta -> [meta.first(), fasta]
            }.set{ ch_fasta }
            MERGE_FASTA(
                ch_fasta
            )
            ch_fasta = MERGE_FASTA.out.fasta
            collected_bams.map{
                m, bams->
                    [ m.id, m, bams ]
                
            }.join(
            ch_fasta
            .map{
                m,fasta -> 
                return [ m.base, fasta ]
            }
            , by:0)
            .map{
                id, m, bams, fasta -> [ m, fasta ]
            }.set { ch_merged_fasta }
        }
        ch_stats = BCFTOOLS_STATS.out.stats
    }


    SAMTOOLS_MERGE(
        collected_bams
    )
    BAM_TO_SAM(
        SAMTOOLS_MERGE.out.bam
    )
    ch_sams=BAM_TO_SAM.out.sam
    ch_pileups=BAM_TO_SAM.out.mpileup
    ch_bams = SAMTOOLS_MERGE.out.bam
    
    
    

    



    emit:
    sam  = ch_sams // channel: [ val(meta), [ paffile ] ] ]
    mpileup = ch_pileups
    bams = ch_bams
    fasta  = ch_merged_fasta
    stats = ch_stats
    versions = ch_versions
}