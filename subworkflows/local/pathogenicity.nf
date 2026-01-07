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

// include minimap2
include { DIAMOND_MAKEDB as PATHOGENICITY_DMND_MAKEDB } from '../../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTX as PATHOGENICITY_DMND_BLASTX } from '../../modules/nf-core/diamond/blastx/main'


workflow PATHOGENICITY {
    take:
        assemblyfiles
    main:
        // ch_empty_file = file("$projectDir/assets/NO_FILE2")
        ch_versions = Channel.empty()
        ch_dmnd_pathogenicity = Channel.empty()
        // assets/vfs_metadata.faa.gz
        ch_amr_fasta = Channel.fromPath("$projectDir/assets/vfs_metadata.faa.gz")
        println "PATHOGENICITY WORKFLOW: Starting pathogenicity subworkflow..."
        // add ch_amr_fasta to all assemblyfiles
        assemblyfiles.combine(ch_amr_fasta).set { assembly_with_amr }
        PATHOGENICITY_DMND_MAKEDB(
            assembly_with_amr.map{ meta, fasta, cds -> [meta, cds] },
            [],
            [],
            [],
        )
        ch_versions = ch_versions.mix(PATHOGENICITY_DMND_MAKEDB.out.versions)

        ch_prep_dmndout = assembly_with_amr.join(PATHOGENICITY_DMND_MAKEDB.out.db)

        PATHOGENICITY_DMND_BLASTX(
            ch_prep_dmndout.map{  m , fasta, cds, dmnd -> [m, fasta] },
            ch_prep_dmndout.map{  m , fasta, cds, dmnd -> [m, dmnd] },
            'txt',
            false
        )
        ch_dmnd_pathogenicity = PATHOGENICITY_DMND_BLASTX.out.txt
    emit:
        versions = ch_versions
        dmnd = ch_dmnd_pathogenicity
}
