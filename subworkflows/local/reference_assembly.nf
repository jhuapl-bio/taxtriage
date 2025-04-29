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

include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/bcftools/consensus/main'
include { BCFTOOLS_MPILEUP } from '../../modules/nf-core/bcftools/mpileup/main'

workflow REFERENCE_ASSEMBLY {
    take:
        ch_aligned_output

    main:
        ch_assemblies = Channel.empty()
        ch_versions = Channel.empty()
        //// // // branch out the samtools_sort output to nanopore and illumina

        BCFTOOLS_MPILEUP(
            ch_aligned_output.map{ m, fastq, fasta, bam -> [m, bam] },
            ch_aligned_output.map{ m, fastq, fasta, bam -> fasta },
            false
        )
        ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)
        ch_merged_mpileup = BCFTOOLS_MPILEUP.out.vcf.join(BCFTOOLS_MPILEUP.out.tbi)
        ch_merged_mpileup = ch_merged_mpileup.join(ch_aligned_output.map{ m, fastq, fasta, bam -> [m, fasta] })


        if (params.reference_assembly){
            BCFTOOLS_CONSENSUS(
                ch_merged_mpileup
            )
            ch_versions = BCFTOOLS_CONSENSUS.out.versions
            ch_fasta = BCFTOOLS_CONSENSUS.out.fasta
        }

    emit:
        ch_assembly = ch_assemblies
        versions = ch_versions
}
