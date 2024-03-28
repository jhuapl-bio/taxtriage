//
// Check input samplesheet and get read channels
//
// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # All rights reserved.
// # Permission is hereby granted, free of charge, to any person obtaining a copy of this
// # software and associated documentation files (the "Software"), to deal in the Software
// # without restriction, including without limitation the rights to use, copy, modify,
// # merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// # permit persons to whom the Software is furnished to do so.
// #
// # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// # INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// # PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// # LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// # TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
// # OR OTHER DEALINGS IN THE SOFTWARE.
// #

include { MINIMAP2_ALIGN as FILTER_MINIMAP2 } from '../../modules/nf-core/minimap2/align/main'
include { BOWTIE2_ALIGN as FILTER_BOWTIE2 } from '../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD  as FILTER_BOWTIE2_IDX } from '../../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FASTQ } from '../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_INDEX as FILTERED_SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS as FILTERED_STATS } from '../../modules/nf-core/samtools/stats/main'

workflow HOST_REMOVAL {
    take:
        ch_reads
        genome

    main:
        // // Remove the human reads first
        ch_bt2_index = Channel.empty()
        ch_filt_illumina = Channel.empty()
        ch_filt_oxfo = Channel.empty()
        ch_filtered_stats = Channel.empty()
        ch_reference_fasta = Channel.empty()
        ch_reads.branch{
            longreads: it[0].platform =~ /(?i)OXFORD/
            shortreads: it[0].platform =~ /(?i)ILLUMINA/
        }.set { ch_aligned_for_filter }


        if (params.remove_reference_file){
            ch_reference_fasta_removal =  file(params.remove_reference_file, checkIfExists: true)
        } else if (params.genome) {
            ch_reference_fasta_removal = params.genomes[params.genome]['fasta']
        }

        if (params.remove_reference_file || params.genome){
            // Run minimap2 module on all OXFORD platform reads and Bowtie2 on ILLUMINA  reads
            // if ch_aligned_for_filter.shorteads is not empty
            // then run bowtie2 on it
            // else run minimap2 on ch_aligned_for_filter.longreads
            if (ch_aligned_for_filter.shortreads){
                ch_meta_reference_fasta = [ [id: 'filterreadsbt2'] , ch_reference_fasta_removal]
                FILTER_BOWTIE2_IDX(
                    ch_meta_reference_fasta
                )
                ch_bt2_index = FILTER_BOWTIE2_IDX.out.index

            }
            FILTER_BOWTIE2(
                ch_aligned_for_filter.shortreads.map{ m, fastq -> return [m, fastq] }.combine(ch_bt2_index.map{ m, idx -> return idx }),
                true,
                true
            )
            FILTER_MINIMAP2(
                ch_aligned_for_filter.longreads.map{ m, fastq -> return [m, fastq, ch_reference_fasta_removal] },
                true,
                true,
                true
            )


            SAMTOOLS_VIEW (
                FILTER_MINIMAP2.out.bam.map{ m, bam ->
                    return [m, bam, []
                ]},
                [ [],[] ],
                []
            )
            SAMTOOLS_FASTQ ( SAMTOOLS_VIEW.out.bam, false )
            ch_reads = SAMTOOLS_FASTQ.out.other.mix(
                FILTER_BOWTIE2.out.fastq
            )
            ch_all_Bams = FILTER_MINIMAP2.out.bam.mix(FILTER_BOWTIE2.out.aligned)

            FILTERED_SAMTOOLS_INDEX ( ch_all_Bams )
            ch_bai_files = ch_all_Bams.join(FILTERED_SAMTOOLS_INDEX.out.bai)
            FILTERED_STATS (
                ch_bai_files,
                [ [], file(params.remove_reference_file, checkIfExists: true) ]
            )
            ch_filtered_stats = FILTERED_STATS.out.stats.collect{it[1]}.ifEmpty([])

        }
        // if (params.pre_remove_taxids){
        //     ch_pre_remove_taxids = Channel.of(params.pre_remove_taxids)
        // }


    emit:
        unclassified_reads = ch_reads
        stats_filtered = ch_filtered_stats
}
