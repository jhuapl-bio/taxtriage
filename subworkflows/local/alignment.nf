//
// Check input samplesheet and get read channels
//

include { BWA_INDEX } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM } from '../../modules/nf-core/bwa/mem/main'
include { HISAT2_ALIGN } from '../../modules/nf-core/hisat2/align/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
// include { BAM_TO_SAM } from "../../modules/local/bam_to_sam"
include { SAMTOOLS_DEPTH } from '../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_COVERAGE } from '../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_HIST_COVERAGE  }  from '../../modules/local/samtools_hist_coverage'
include { MERGE_FASTA } from '../../modules//local/merge_fasta'
include { SPLIT_VCF } from '../../modules/local/split_vcf'
include { BOWTIE2_ALIGN  } from '../../modules/nf-core/bowtie2/align/main'
include { RSEQC_BAMSTAT } from '../../modules/nf-core/rseqc/bamstat/main'
include { BEDTOOLS_DEPTHCOVERAGE } from '../../modules/local/bedtools_coverage'
include { REFERENCE_ASSEMBLY } from './reference_assembly'

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


    ch_default_aligner = params.default_aligner
    fastq_reads
        .branch{
            bowtie2: params.default_aligner =~ /(?i)bowtie2/ && it[0].platform =~ /(?i)illumina/ || it[0].aligner =~ /(?i)bowtie2/
            hisat2: params.default_aligner =~ /(?i)hisat2/ || it[0].aligner =~ /(?i)hisat2/
            minimap2: true
    }.set { ch_aligners }



    ch_versions = 1

    def idx = 0

    ch_aligners.bowtie2
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
        .set { ch_fasta_bowtie2_files_for_alignment }
    ch_aligners.minimap2
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
                return [ mm, fastq, fasta]
            }
        }
        .set { ch_fasta_minimap2_files_for_alignment }
    ch_aligners.hisat2
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
                return [ mm, fastq, fasta]
            }
        }
        .set { ch_fasta_hisat2_files_for_alignment }

    // make empty channel
    ch_aligned_output = Channel.empty()

    MINIMAP2_ALIGN(
        ch_fasta_minimap2_files_for_alignment,
        true,
        true,
        true,
        params.minmapq
    )
    ch_fasta_minimap2_files_for_alignment = ch_fasta_minimap2_files_for_alignment
        .map{ m, fastq, fasta -> { return [ m, fastq, fasta[0] ] }  }
        .join(MINIMAP2_ALIGN.out.bam)



    BOWTIE2_ALIGN(
        ch_fasta_bowtie2_files_for_alignment,
        true,
        true,
        params.minmapq
    )
    ch_fasta_bowtie2_files_for_alignment = ch_fasta_bowtie2_files_for_alignment
        .map{ m, fastq, fasta -> { return [ m, fastq, fasta[0] ] }  }
        .join(BOWTIE2_ALIGN.out.aligned)

    HISAT2_ALIGN(
        ch_fasta_hisat2_files_for_alignment.map{ m, fastq, fasta -> [m, fastq] },
        ch_fasta_hisat2_files_for_alignment.map{ m, fastq, fasta -> [m, fasta] },
        ch_fasta_hisat2_files_for_alignment.map{ m, fastq, fasta -> [m, null] },
    )
    ch_fasta_hisat2_files_for_alignment = ch_fasta_hisat2_files_for_alignment
        .map{ m, fastq, fasta -> [m, fastq, fasta[0]] }
        .join(HISAT2_ALIGN.out.bam)


    // Merging all outputs here to a single channel for further processing
    ch_aligned_output = ch_fasta_minimap2_files_for_alignment
        .mix(ch_fasta_bowtie2_files_for_alignment)
        .mix(ch_fasta_hisat2_files_for_alignment)



    if (params.get_variants || params.reference_assembly){
        REFERENCE_ASSEMBLY(
            ch_aligned_output
        )
    }
    collected_bams = ch_aligned_output.map { meta, fastqs, fasta, bam -> [meta.oid, bam] } // Extract oid and file
        .groupTuple() // Group by oid
        .map { oid, files -> [[id:oid], files ] } // Replace oid with id in the metadata
        .set{ merged_bams_channel }


    // // Split the channel based on the condition
    merged_bams_channel
        .branch {
            mergeNeeded: it[1].size() > 1
            noMergeNeeded: it[1].size() == 1
        }
        .set { branchedChannels }

    SAMTOOLS_MERGE(
        branchedChannels.mergeNeeded
    )

    // // sort the multi bams from the merge
    SAMTOOLS_SORT(
        SAMTOOLS_MERGE.out.bam
    )
    //  // // Unified channel from both merged and non-merged BAMs
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
    // // // Example to view the output
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

    collected_bams.view()

    gcf_with_bam = collected_bams.join(fastq_reads.map{ m, fastq, fasta, map -> return [m, map] })

    SAMTOOLS_HIST_COVERAGE(
        gcf_with_bam
    )

    sorted_bams_with_index.join(fastq_reads.map{ m, fastq, fasta, map -> return [m, fasta] }).set{ sorted_bams_with_index_fasta }

    // // RSEQC_BAMSTAT(
    // //     collected_bams
    // // )
    // // ch_bamstats = RSEQC_BAMSTAT.out.txt
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
