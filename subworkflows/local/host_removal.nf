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
include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main'
include { REMOVE_HOSTREADS } from '../../modules/local/remove_unaligned'
include { SAMTOOLS_INDEX as FILTERED_SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS as FILTERED_STATS } from '../../modules/nf-core/samtools/stats/main'
include { CHECK_GZIPPED_READS } from '../../modules/local/check_reads_exist'

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

        if (params.remove_reference_file){
            ch_reference_fasta_removal =  file(params.remove_reference_file, checkIfExists: true)
        } else if (params.genome) {
            ch_reference_fasta_removal = params.genomes[params.genome]['fasta']
        }

        if (params.remove_reference_file || params.genome){
            // Run minimap2 module on all OXFORD platform reads and Bowtie2 on ILLUMINA  reads
            // if ch_aligned_for_filter.shorteads is not empty
            // Run minimap2 on all for host removal - as host removal outperforms bowtie2 for host false negative rate https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9040843/

            FILTER_MINIMAP2(
                ch_reads.map{ m, fastq -> return [m, fastq, ch_reference_fasta_removal] },
                true,
                true,
                true,
                0
            )

            ch_bam_hosts = FILTER_MINIMAP2.out.bam
            REMOVE_HOSTREADS(ch_bam_hosts)
            ch_filtered_reads = REMOVE_HOSTREADS.out.reads

            // Check the filtered output and fallback to original reads if filtered reads are empty
            CHECK_GZIPPED_READS(ch_filtered_reads, 4)
            ch_valid_reads = CHECK_GZIPPED_READS.out.check_result
            ch_orig_reads = ch_valid_reads.filter({
                it[1].name == 'emptyfile.txt'
            }).join(ch_reads).map({
                meta, result, reads -> return [meta, reads]
            })
            ch_filtered_reads = ch_valid_reads.filter({
                it[1].name == 'minimum_reads_check.txt'
            }).join(ch_filtered_reads).map({
                meta, result, reads -> return [meta, reads]
            })
            ch_reads = ch_orig_reads.mix(ch_filtered_reads)

            // Continue processing the final reads
            FILTERED_SAMTOOLS_INDEX(ch_bam_hosts)

            ch_bai_files = ch_bam_hosts.join(FILTERED_SAMTOOLS_INDEX.out.bai)
            FILTERED_STATS (
                ch_bai_files,
                [ [], file(params.remove_reference_file) ]
            )
            ch_filtered_stats = FILTERED_STATS.out.stats.collect{it[1]}.ifEmpty([])

        }
        // filter out all ch_reads fastq files that are empty

    emit:
        unclassified_reads = ch_reads
        stats_filtered = ch_filtered_stats
}
