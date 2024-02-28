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

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
        reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.platform = row.platform ? row.platform : 'ILLUMINA'
    meta.fastq_1 = row.fastq_1
    meta.fastq_2 = row.fastq_2
    // if meta.fastq_2 it is not single end, set meta.single_end as true else meta.single_end is false
    meta.single_end = row.fastq_2  ? false : true
    if (row.trim && row.trim.toLowerCase() == "true"){
        meta.trim = true
    } else if (!row.trim  || (row.trim && row.trim.toLowerCase() == "false")){
        meta.trim = false
    }
    meta.type = row.type
    meta.directory = row.directory ?  row.directory.toBoolean() : null
    meta.sequencing_summary = row.sequencing_summary ? file(row.sequencing_summary) : null
    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (meta.directory ){
        if (meta.platform == 'OXFORD'){
            fastq_meta = [ meta, [ file(meta.fastq_1) ]  ]
        } else {
            exit 1, "ERROR: Please check input samplesheet -> the platform is not specified as OXFORD \n${meta.sample}"
        }
    } else {
        if (!file(meta.fastq_1).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${meta.fastq_1}"
        }
        if (meta.single_end) {
            fastq_meta = [ meta, [ file(meta.fastq_1) ] ]
        } else {
            if (meta.fastq_2 && !file(meta.fastq_2).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${meta.fastq_2}"
            }
            fastq_meta = [ meta, [ file(meta.fastq_1), file(meta.fastq_2) ] ]
        }
    }

    return fastq_meta
}
