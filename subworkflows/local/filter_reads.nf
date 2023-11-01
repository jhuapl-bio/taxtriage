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

include { KRAKEN2_KRAKEN2 as FILTERREADSKRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { MOVE_FILES } from '../../modules/local/moveFiles.nf'

workflow READSFILTER {
    take:
        fastq_reads
        db

    main:
        FILTERREADSKRAKEN2(
            fastq_reads,
            db,
            true,
            true
        )
        MOVE_FILES(
            FILTERREADSKRAKEN2.out.unclassified_reads_fastq,
            "filtered_",
            false,
            []
        )
        ch_unclassified_reads = MOVE_FILES.out.files


    emit:
        reads = ch_unclassified_reads
        versions = FILTERREADSKRAKEN2.out.versions
}
