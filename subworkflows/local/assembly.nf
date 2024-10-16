//
// Cross-check all upstream steps for final confidence report
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
//
include { MEGAHIT } from '../../modules/nf-core/megahit/main'
include { MEGAHIT as MEGAHIT_LONG } from '../../modules/nf-core/megahit/main'
include { DIAMOND_MAKEDB } from '../../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTX } from '../../modules/nf-core/diamond/blastx/main'
include { BEDTOOLS_COVERAGE } from '../../modules/nf-core/bedtools/coverage/main'
include { FEATURES_MAP } from '../../modules/local/features_map'
include { MAP_PROT_ASSEMBLY } from '../../modules/local/map_prot_assembly'
include { FLYE } from '../../modules/nf-core/flye/main'

workflow ASSEMBLY {
    take:
        postalignmentfiles
        assemblyfile
    main:
        ch_empty_file = file("$projectDir/assets/NO_FILE2")
        ch_assembly_report = Channel.empty()
        ch_diamond_output = Channel.empty()
        ch_diamond_analysis = Channel.empty()
        ch_assembly_alignment = Channel.empty()
        ch_diamond_analysis = postalignmentfiles.map{ [it[0], ch_empty_file] }
        if (!params.skip_denovo){
            // branch long and short reads
            postalignmentfiles.branch {
                shortreads: it[0].platform =~ 'ILLUMINA'
                longreads: it[0].platform != 'ILLUMINA'
            }.set { branchedChannels }

            ch_longreads_assembled = Channel.empty()
            if (params.use_megahit_longreads){
                MEGAHIT_LONG(
                    branchedChannels.longreads.map{ meta, bam, bai, mapping, bed, cds, mapcd,  reads -> [meta, reads] }
                )
                ch_longreads_assembled = MEGAHIT_LONG.out.contigs
            } else {
                FLYE(
                    branchedChannels.longreads.map{ meta, bam, bai, mapping, bed, cds, mapcd,  reads -> [meta, reads] },
                    '--nano-raw'
                )
                ch_longreads_assembled = FLYE.out.fasta
            }
            MEGAHIT(
                branchedChannels.shortreads.map{ meta, bam, bai, mapping, bed, cds, mapcd,  reads -> [meta, reads] }
            )
            ch_assembled_files = MEGAHIT.out.contigs.mix(ch_longreads_assembled)

            BEDTOOLS_COVERAGE(
                postalignmentfiles.map{ meta, bam, bai, mapping, bed, cds, mapcd, reads -> [meta, bed, bam] }
            )
            ch_bedout = BEDTOOLS_COVERAGE.out.bed.join(
                postalignmentfiles.map{ meta, bam, bai, mapping, bed, cds, mapcd, reads -> [meta, mapping] }
            )
            FEATURES_MAP(
                ch_bedout
            )

            try {

                valid_aligners  = postalignmentfiles.filter{
                    return it[5] != []
                }
                DIAMOND_MAKEDB(
                    valid_aligners.map{ meta, bam, bai, mapping, bed, cds, mapcd, reads -> [meta, cds] }
                )
                DIAMOND_BLASTX(
                    ch_assembled_files.join(DIAMOND_MAKEDB.out.db),
                    'txt',
                    false
                )
                MAP_PROT_ASSEMBLY(
                    DIAMOND_BLASTX.out.txt.join(valid_aligners.map{meta, bam, bai, mapping, bed, cds, mapcd, reads -> [meta, mapcd]}),
                    assemblyfile
                )
                ch_diamond_analysis = MAP_PROT_ASSEMBLY.out.promap
                // Run minimap2 on the contigs against reference fasta files
                ch_diamond_output = DIAMOND_BLASTX.out.txt
            } catch (Exception e){
                ch_diamond_analysis = postalignmentfiles.map{ [it[0], ch_empty_file] }
            }
        } else {
            postalignmentfiles.map{ meta, bam, bai, mapping, bed, cds, mapcd, reads  -> {
                    return [ meta,  ch_empty_file]
            }
            }.set{ ch_diamond_output }
        }
    emit:
        ch_diamond_output
        ch_diamond_analysis
}
