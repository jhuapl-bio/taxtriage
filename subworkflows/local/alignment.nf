//
// Check input samplesheet and get read channels
//

include { BWA_INDEX } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM } from '../../modules/nf-core/bwa/mem/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
// include { BAM_TO_SAM } from "../../modules/local/bam_to_sam"
include { SAMTOOLS_DEPTH } from '../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_COVERAGE } from '../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_HIST_COVERAGE  }  from '../../modules/local/samtools_hist_coverage'
include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/bcftools/consensus/main'
include { BCFTOOLS_MPILEUP as BCFTOOLS_MPILEUP_ILLUMINA } from '../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_INDEX  } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_MPILEUP as BCFTOOLS_MPILEUP_OXFORD } from '../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_STATS } from '../../modules/nf-core/bcftools/stats/main'
include { MERGE_FASTA } from '../../modules//local/merge_fasta'
include { SPLIT_VCF } from '../../modules/local/split_vcf'
include { BOWTIE2_ALIGN  } from '../../modules/nf-core/bowtie2/align/main'
include { RSEQC_BAMSTAT } from '../../modules/nf-core/rseqc/bamstat/main'
include { BEDTOOLS_DEPTHCOVERAGE } from '../../modules/local/bedtools_coverage'

workflow ALIGNMENT {
    take:
    fastq_reads

    main:
    ch_bams = Channel.empty()
    ch_fasta = Channel.empty()
    ch_pileups = Channel.empty()
    ch_stats = Channel.empty()
    ch_depths = Channel.empty()
    ch_pafs = Channel.empty()
    ch_sams = Channel.empty()
    ch_aligners = Channel.empty()
    collected_bams  = Channel.empty()
    ch_bamstats = Channel.empty()
    ch_merged_fasta = Channel.empty()
    ch_merged_mpileup = Channel.empty()




    fastq_reads
        .branch{
            longreads: it[0].platform =~ 'OXFORD'
            shortreads: it[0].platform =~ 'ILLUMINA'
    }.set { ch_aligners }

    ch_versions = 1

    def idx = 0

    ch_aligners.shortreads
        .flatMap { meta, fastq, fastas, _ ->
            def size = fastas.size()
            fastas.collect{ fasta ->

                def id = "${meta.id}"
                if (size > 1){
                    def basen = fasta[0].getBaseName()
                    id = "${id}.${basen}"
                }
                def mm = [id: id, oid: meta.id, single_end: meta.single_end  ]
                idx++

                return [ mm, fastq, fasta]
            }
        }
        .set { ch_fasta_shortreads_files_for_alignment }
    ch_aligners.longreads
        .flatMap { meta, fastq, fastas, _ ->
            def size = fastas.size()
            fastas.collect{ fasta ->
                def id = "${meta.id}"
                if (size > 1){
                    def basen = fasta[0].getBaseName()
                    id = "${id}.${basen}"
                }
                def mm = [id: id,  oid: meta.id, single_end: meta.single_end  ]
                idx++
                return [ mm, fastq, fasta[0]]
            }
        }
        .set { ch_fasta_longreads_files_for_alignment }

    MINIMAP2_ALIGN(
        ch_fasta_longreads_files_for_alignment,
        true,
        true,
        true,
        params.minmapq
    )

    BOWTIE2_ALIGN(
        ch_fasta_shortreads_files_for_alignment,
        true,
        true,
        params.minmapq
    )

    collected_bams = BOWTIE2_ALIGN.out.aligned
        .mix(MINIMAP2_ALIGN.out.bam)




    sorted_bams = collected_bams


    sorted_bams
        .map { meta, file -> [meta.oid, file] } // Extract oid and file
        .groupTuple() // Group by oid
        .map { oid, files -> [[id:oid], files] } // Replace oid with id in the metadata
        .set{ merged_bams_channel }


    if (params.get_variants || params.reference_assembly){
    //     // // // branch out the samtools_sort output to nanopore and illumina


        BCFTOOLS_MPILEUP_OXFORD(
            ch_fasta_longreads_files_for_alignment.map{ m, bam, fasta -> [m, bam] },
            ch_fasta_longreads_files_for_alignment.map{ m, bam, fasta -> fasta },
            false
        )
        ch_merged_mpileup_oxford = BCFTOOLS_MPILEUP_OXFORD.out.vcf.join(BCFTOOLS_MPILEUP_OXFORD.out.tbi)
        ch_merged_mpileup_oxford = ch_merged_mpileup_oxford.join(ch_fasta_longreads_files_for_alignment.map{ m, bam, fasta -> [m, fasta] })
        BCFTOOLS_MPILEUP_ILLUMINA(
            ch_fasta_shortreads_files_for_alignment.map{ m, bam, fasta -> [m, bam] },
            ch_fasta_shortreads_files_for_alignment.map{ m, bam, fasta -> fasta[0] },
            false
        )
        .set{ transposed_fastas_oxford }

        ch_merged_mpileup_illumina = BCFTOOLS_MPILEUP_ILLUMINA.out.vcf.join(BCFTOOLS_MPILEUP_ILLUMINA.out.tbi)
        ch_merged_mpileup_illumina = ch_merged_mpileup_illumina.join(ch_fasta_shortreads_files_for_alignment.map{ m, bam, fasta -> [m, fasta[0]] })

        ch_merged_mpileup = ch_merged_mpileup_illumina.mix(ch_merged_mpileup_oxford)

        if (params.reference_assembly){
            BCFTOOLS_CONSENSUS(
                ch_merged_mpileup
            )
            ch_fasta = BCFTOOLS_CONSENSUS.out.fasta
        }
    }
    // Split the channel based on the condition
    merged_bams_channel
        .branch {
            mergeNeeded: it[1].size() > 1
            noMergeNeeded: it[1].size() == 1
        }
        .set { branchedChannels }


    SAMTOOLS_MERGE(
        branchedChannels.mergeNeeded
    )
    // sort the multi bams from the merge
    SAMTOOLS_SORT(
        SAMTOOLS_MERGE.out.bam
    )

     // // Unified channel from both merged and non-merged BAMs
    SAMTOOLS_SORT.out.bam.mix( branchedChannels.noMergeNeeded ).set { collected_bams }



    // // Join the channels on 'id' and append the BAM files to the fastq_reads entries
    fastq_reads.map { item -> [item[0].id, item] } // Map to [id, originalItem]
        .join(collected_bams.map { item -> [item[0].id, item] }, by: 0)
        .map { joinedItems ->
            def item1 = joinedItems[1] // Original item from Channel1
            def item2 = joinedItems[2] // Original item from Channel2
            // Now you can merge item1 and item2 as needed
            return [item1[0], item2[1]] // Adjust based on how you need the merged items
        }.set{ collected_bams }

    // // Example to view the output
    SAMTOOLS_DEPTH(
        collected_bams
    )
    SAMTOOLS_INDEX(
        collected_bams
    )
    collected_bams.join(SAMTOOLS_INDEX.out.bai).set{ sorted_bams_with_index }


    SAMTOOLS_COVERAGE(
        sorted_bams_with_index
    )
    ch_stats = SAMTOOLS_COVERAGE.out.coverage

    gcf_with_bam = collected_bams.join(fastq_reads.map{ m, fastq, fasta, map -> return [m, map] })

    SAMTOOLS_HIST_COVERAGE (
        gcf_with_bam
    )

    sorted_bams_with_index.join(fastq_reads.map{ m, fastq, fasta, map -> return [m, fasta] }).set{ sorted_bams_with_index_fasta }

    // RSEQC_BAMSTAT(
    //     collected_bams
    // )
    // ch_bamstats = RSEQC_BAMSTAT.out.txt
    ch_bams =  sorted_bams_with_index
    ch_depths = SAMTOOLS_DEPTH.out.tsv
    emit:
        depth = ch_depths
        mpileup = ch_merged_mpileup
        bams = ch_bams
        fasta  = ch_merged_fasta
        stats = ch_stats
        // bamstats = ch_bamstats
        versions = ch_versions
    }
