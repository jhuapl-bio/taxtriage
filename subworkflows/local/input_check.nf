//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    println samplesheet
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
    meta.single_end = row.single_end.toBoolean()
    meta.platform = row.platform
    meta.barcode = row.barcode.toBoolean()
    meta.trim = row.trim.toBoolean()
    meta.sequencing_summary = row.sequencing_summary ? file(row.sequencing_summary) : null
    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    
    if (meta.barcode){
        if (row.platform == 'OXFORD'){
            meta.bc = row.bc
            meta.from = row.from
            fastq_meta = [ meta, [ file(row.from) ]  ]
        } else {
            exit 1, "ERROR: Please check input samplesheet -> the platform is not specified as OXFORD \n${row.sample}"
        }
    } else {
        if (!file(row.fastq_1).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
        }
        if (meta.single_end) {
            fastq_meta = [ meta, [ file(row.fastq_1) ] ]
        } else {
            if (!file(row.fastq_2).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
            }
            fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
        }
    }
    
    return fastq_meta
}
