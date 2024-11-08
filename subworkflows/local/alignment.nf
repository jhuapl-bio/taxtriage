//
// Check input samplesheet and get read channels
//

include { BWA_INDEX } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM } from '../../modules/nf-core/bwa/mem/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
include { HISAT2_ALIGN } from '../../modules/nf-core/hisat2/align/main'
include { SAMTOOLS_DEPTH } from '../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_COVERAGE } from '../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_HIST_COVERAGE  }  from '../../modules/local/samtools_hist_coverage'
include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/bcftools/consensus/main'
include { BCFTOOLS_MPILEUP } from '../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_INDEX  } from '../../modules/nf-core/bcftools/index/main'
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
    ch_versions = 1

    fastq_reads.set { ch_aligners }
    def idx = 0

    ch_aligners
        .flatMap { meta, fastq, fastas, _ ->
            // Print 'fastas' for debugging purposes

            def outputs = []

            // Flatten 'fastas' by one level if it's nested
            def flattenedFastas = fastas.collectMany { it }

            flattenedFastas.each { fastaItem -> {
                    // if fastaItem is a list
                    if (fastaItem instanceof List) {
                        // If the item is a list of files, return each file separately
                        def fasta = fastaItem[0]
                        def id = "${meta.id}.${fasta.getBaseName()}"
                        def mm = [id: id,  oid: meta.id, single_end: meta.single_end, platform: meta.platform   ]
                        outputs << [mm, fastq, fastaItem]
                    } else {
                        // If the item is a single file, return it as is
                        def id = "${meta.id}.${fastaItem.getBaseName()}"
                        def mm = [id: id,  oid: meta.id, single_end: meta.single_end, platform: meta.platform   ]
                        outputs.add([mm, fastq, fastaItem])
                    }
                }
            }
            // Return the collected outputs
            return outputs
        }.set { ch_fasta_files_for_alignment }



    if (params.use_bt2){
        BOWTIE2_ALIGN(
            ch_fasta_files_for_alignment,
            true,
            true
        )

        collected_bams = BOWTIE2_ALIGN.out.bam

    } else if (params.use_hisat2) {
        // Set null for hisat2 splicesites
        HISAT2_ALIGN(
            ch_fasta_files_for_alignment.map{ m, fastq, fasta -> [m, fastq] },
            ch_fasta_files_for_alignment.map{ m, fastq, fasta -> [m, fasta] },
            ch_fasta_files_for_alignment.map{ m, fastq, fasta -> [m, null ] },
        )

        collected_bams = HISAT2_ALIGN.out.bam
    } else {
        MINIMAP2_ALIGN(
            ch_fasta_files_for_alignment,
            true,
            true,
            true
        )

        collected_bams = MINIMAP2_ALIGN.out.bam
    }

    sorted_bams = collected_bams

    sorted_bams
        .map { meta, file -> [meta.oid, file] } // Extract oid and file
        .groupTuple() // Group by oid
        .map { oid, files -> [[id:oid], files] } // Replace oid with id in the metadata
        .set{ merged_bams_channel }

    if (params.get_variants || params.reference_assembly){
        //     // // // branch out the samtools_sort output to nanopore and illumina
        BCFTOOLS_MPILEUP(
            ch_fasta_files_for_alignment.map{ m, bam, fasta -> [m, bam] },
            ch_fasta_files_for_alignment.map{ m, bam, fasta -> fasta },
            false
        )
        ch_merged_mpileup = BCFTOOLS_MPILEUP.out.vcf.join(BCFTOOLS_MPILEUP.out.tbi)
        ch_merged_mpileup = ch_merged_mpileup.join(ch_fasta_files_for_alignment.map{ m, bam, fasta -> [m, fasta] })
        BCFTOOLS_MPILEUP(
            ch_fasta_files_for_alignment.map{ m, bam, fasta -> [m, bam] },
            ch_fasta_files_for_alignment.map{ m, bam, fasta -> fasta[0] },
            false
        )
        .set{ transposed_fastas }


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
    //  // // Unified channel from both merged and non-merged BAMs
    SAMTOOLS_MERGE.out.bam.mix ( branchedChannels.noMergeNeeded ).set { collected_bams }

    // sort the multi bams from the merge

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
    collected_bams.join(SAMTOOLS_INDEX.out.csi).set{ sorted_bams_with_index }

    SAMTOOLS_COVERAGE(
        sorted_bams_with_index
    )
    ch_stats = SAMTOOLS_COVERAGE.out.coverage

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
