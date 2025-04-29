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
        ch_versions = Channel.empty()

        ch_diamond_analysis = postalignmentfiles.map{ [it[0], ch_empty_file] }
        if (params.use_denovo || params.use_diamond){
            // branch long and short reads
            postalignmentfiles.branch {
                shortreads: it[0].platform =~ 'ILLUMINA'
                longreads: it[0].platform != 'ILLUMINA'
            }.set { branchedChannels }

            ch_longreads_assembled = Channel.empty()
            if (params.use_megahit_longreads){
                MEGAHIT_LONG(
                    branchedChannels.longreads.map{ meta, bam, bai, mapping, bed, cds, features, mapcd,  reads -> [meta, reads, []] }
                )
                ch_versions = ch_versions.mix(MEGAHIT_LONG.out.versions)
                ch_longreads_assembled = MEGAHIT_LONG.out.contigs
            } else {
                FLYE(
                    branchedChannels.longreads.map{ meta, bam, bai, mapping, bed, cds, features, mapcd,  reads -> [meta, reads] },
                    '--nano-raw'
                )
                ch_versions = ch_versions.mix(FLYE.out.versions)
                ch_longreads_assembled = FLYE.out.fasta
            }
            // if paired end then make a reads1 and reads2, else do single reads input
            MEGAHIT(
                branchedChannels.shortreads.map{ meta, bam, bai, mapping, bed, cds, features, mapcd,  reads -> {
                        if (meta.single_end) {
                            return [meta, reads, []]
                        } else {
                            return [meta, reads[0], reads[1]]
                        }
                    }
                }
            )
            ch_versions = ch_versions.mix(MEGAHIT.out.versions)
            ch_assembled_files = MEGAHIT.out.contigs.mix(ch_longreads_assembled)
            if ((params.use_diamond)){
                BEDTOOLS_COVERAGE(
                    postalignmentfiles.map{ meta, bam, bai, mapping, bed, cds, features, mapcd, reads -> [meta, bed, bam] },
                    postalignmentfiles.map{ meta, bam, bai, mapping, bed, cds, features, mapcd, reads -> [] },
                )
                ch_bedout = BEDTOOLS_COVERAGE.out.bed.join(
                    postalignmentfiles.map{ meta, bam, bai, mapping, bed, cds, features, mapcd, reads -> [meta, mapping] }
                )
                FEATURES_MAP(
                    ch_bedout
                )
            }

            try {
                valid_aligners  = postalignmentfiles.filter{
                    return it[5] != []
                }
                DIAMOND_MAKEDB(
                    valid_aligners.map{ meta, bam, bai, mapping, bed, cds, features, mapcd, reads -> [meta, cds] },
                    [],
                    [],
                    [],
                )
                ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions)

                ch_prep_dmndout = ch_assembled_files.join(DIAMOND_MAKEDB.out.db)

                DIAMOND_BLASTX(
                    ch_prep_dmndout.map{  m , fasta, dmnd -> [m, fasta] },
                    ch_prep_dmndout.map{  m , fasta, dmnd -> [m, dmnd] },
                    'txt',
                    false
                )
                ch_versions = ch_versions.mix(DIAMOND_BLASTX.out.versions)
                MAP_PROT_ASSEMBLY(
                    DIAMOND_BLASTX.out.txt.join(valid_aligners.map{meta, bam, bai, mapping, bed, cds, features, mapcd, reads -> [meta, mapcd, features, mapping] }),
                    assemblyfile
                )
                ch_versions = ch_versions.mix(MAP_PROT_ASSEMBLY.out.versions)
                ch_diamond_analysis = ch_diamond_analysis.join(MAP_PROT_ASSEMBLY.out.promap, remainder: true)
                ch_diamond_analysis = ch_diamond_analysis.map{
                    meta, nullfile, promap-> {
                        if (promap == null){
                            return [meta, ch_empty_file]
                        } else {
                            return [meta, promap]
                        }
                    }
                }
                // for any meta.id missing from diamond analysis, add empty file

                // Run minimap2 on the contigs against reference fasta files
                ch_diamond_output = DIAMOND_BLASTX.out.txt
            } catch (Exception e){
                ch_diamond_analysis = postalignmentfiles.map{ [it[0], ch_empty_file] }
            }
        } else {
            postalignmentfiles.map{ meta, bam, bai, mapping, bed, cds, features, mapcd, reads  -> {
                    return [ meta,  ch_empty_file]
            }
            }.set{ ch_diamond_output }

        }
    emit:
        ch_diamond_output
        ch_diamond_analysis
        versions = ch_versions
}
