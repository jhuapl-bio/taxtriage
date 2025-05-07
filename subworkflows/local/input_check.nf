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
include { GENERATE_SAMPLESHEET } from '../../modules/local/samplesheet_generate'

workflow INPUT_CHECK {


    main:


    if (params.fastq_1) {
        // Create a synthetic samplesheet from fastq params
        GENERATE_SAMPLESHEET(
            [
                sampleName: params.sample ?: 'sample',
                platform  : params.platform ?: 'ILLUMINA',
                fastq_1   : params.fastq_1,
                fastq_2   : params.fastq_2,
                seq_summary: params.seq_summary,
                trim      : params.trim,
                type      : params.type
            ]
        )

        GENERATE_SAMPLESHEET.out.csv
            .splitCsv(header: true, sep: ',')
            .map { create_fastq_channel(it) }
            .set { reads }
        versions = GENERATE_SAMPLESHEET.out.versions
    } else if (params.input) {
        // Use the provided samplesheet
        SAMPLESHEET_CHECK(file(params.input))
            .csv
            .splitCsv(header: true, sep: ',')
            .map { create_fastq_channel(it) }
            .set { reads }
            versions = SAMPLESHEET_CHECK.out.versions
        }
    else {
        error "ERROR: Must specify either --input or --fastq_1"
    }

    emit:
        reads                                     // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    // if fastq_2 is not a column then set it as null for all rows


    meta.id         = row.sample
    meta.platform = row.platform ? row.platform : 'ILLUMINA'
    meta.fastq_1 = row.fastq_1
    // Check if 'fastq_2' exists in 'row'
    if (row.containsKey('fastq_2')) {
        meta.fastq_2 = row.fastq_2
    } else {
        meta.fastq_2 = null
    }
    meta.needscompressing = row.needscompressing ? row.needscompressing : null

    // if meta.needscompressing is null or false AND the filename ends with .fastq or .fq then set to true
    // if (!meta.needscompressing && (meta.fastq_1.endsWith('.fastq') || meta.fastq_1.endsWith('.fq'))) {
    //     meta.needscompressing = true
    // }
    // if meta.fastq_2 it is not single end, set meta.single_end as true else meta.single_end is false
    meta.single_end = row.fastq_2  ? false : true
    meta.aligner  = row.aligner ? row.aligner : 'minimap2'
    // if meta.aligner is not minimap2, hisat2, or bowtie2 then exit and send error
    if (meta.aligner != 'minimap2' && meta.aligner != 'hisat2' && meta.aligner != 'bowtie2') {
        exit 1, "ERROR: Please check input samplesheet -> aligner is not specified as minimap2, hisat2, or bowtie2 \n${meta.sample}"
    }
    // Check if fastq_1 exists if not then error out and print error
    if (!file(meta.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${meta.fastq_1}"
    }
    if (!meta.single_end && !file(meta.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${meta.fastq_2}"
    }

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
        if (meta.platform == 'OXFORD' || meta.platform == "PACBIO"){
            fastq_meta = [ meta, [ file(meta.fastq_1) ]  ]
        } else {
            exit 1, "ERROR: Please check input samplesheet -> the platform is not specified as OXFORD or PACBIO \n${meta.sample}"
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
