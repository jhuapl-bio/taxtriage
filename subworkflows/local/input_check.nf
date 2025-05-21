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
    ch_infile = Channel.empty()

    if (params.fastq_1) {
        GENERATE_SAMPLESHEET(
            [
                sampleName: params.sample ?: 'sample',
                platform  : params.platform ?: 'ILLUMINA',
                fastq_1   : params.fastq_1,
                fastq_2   : params.fastq_2,
                seq_summary: params.seq_summary,
                trim      : params.trim,
                type      : params.type,
                batch     : params.batch,
            ]
        )

        ch_infile = GENERATE_SAMPLESHEET.out.csv
        versions = GENERATE_SAMPLESHEET.out.versions
    } else if (params.input) {
        println "INFO: Using input file: ${params.input}"
        ch_infile = file(params.input)
    } else {
        error "ERROR: Must specify either --input or --fastq_1 or --batch"
    }
    // Use the provided samplesheet
    SAMPLESHEET_CHECK(ch_infile)
        .csv
        .splitCsv(header: true, sep: ',')
        .flatMap { row -> create_fastq_channel(row) }
        .set { reads }

        reads.view()

        versions = SAMPLESHEET_CHECK.out.versions


    emit:
        reads                                     // channel: [ val(meta), [ reads ] ]
}


/**
 * @param row  A LinkedHashMap containing:
 *             sample, platform, fastq_1, fastq_2, sequencing_summary,
 *             trim, type, batch
 * @return     A List of [ metaMap, List<File> ] tuples
 */
 def create_fastq_channel( LinkedHashMap row ) {
    def results = []

    // ─── declare all flags / inputs up front ────────────────────────────────────
    boolean isBatch   = row.batch?.toString()?.toLowerCase() == 'true'
    boolean doTrim    = row.trim?.toString()?.toLowerCase()   == 'true'
    String  samp      = row.sample
    String  plat      = row.platform ?: 'ILLUMINA'
    String  seqSum    = row.sequencing_summary ?: ''
    String  type      = row.type ?: 'UNKNOWN'

    // ─── build our "inPath" and verify it exists ───────────────────────────────
    def inPath = new File(row.fastq_1)
    if( ! inPath.exists() ) {
        error "ERROR: Path does not exist: ${inPath}"
    }
    println "INFO: Input path: ${row}"
    def dir = inPath.isDirectory() ? inPath : inPath.parentFile
    if( isBatch && dir ) {
        // ─── batch mode: scan only the top level for .fastq/.fq (with or without .gz)
        def fastqs = dir.listFiles().findAll {
            it.name ==~ /(?i).+\.(fastq|fq)(\.gz)?$/
        }

        // if Oxford (or PacBio), skip R1/R2 pairing entirely:
        if( plat.equalsIgnoreCase('OXFORD') ) {
            fastqs.each { f ->
                // derive a sample ID from the filename (strip ext)
                def sampleId = f.name.replaceAll(/(?i)\.(fastq|fq)(?:\.gz)?$/, '')
                def needsCompress = f.name ==~ /(?i).+\.(fastq|fq)$/

                def meta = [
                    id                 : sampleId,
                    platform           : plat,
                    type               : row.type   ?: 'UNKNOWN',
                    trim               : row.trim.toBoolean() ?: false,
                    sequencing_summary : row.sequencing_summary ?: null,
                    single_end         : true,
                    directory          : false,
                    batch              : true,
                    fastq_1            : f1,
                    fastq_2            : null,
                    needscompressing   : needsCompress
                ]
                results << [ meta, [ file(f1) ] ]
            }
        }
        else {
            // non-Oxford: do your normal R1/R2 grouping
            def pattern = ~/(?i)(.+?)(?:[_.]R?([12]))?\.(?:fastq|fq)(?:\.gz)?$/
            def groups  = fastqs.groupBy { f ->
                def m = (f.name =~ pattern)
                m ? m[0][1] : f.name
            }
            groups.each { base, files ->
                def f1 = files.find{ it.name =~ /(?i)[_.]R?1\./ } ?: files[0]
                def f2 = files.find{ it.name =~ /(?i)[_.]R?2\./ }
                def needsCompress = f1.name ==~ /(?i).+\.(fastq|fq)$/

                // strip off the R1 suffix for your sample ID
                def sampleId = f1.name.replaceAll(/(?i)[_.]R?1\.(?:fastq|fq)(?:\.gz)?$/, '')

                def meta = [
                    id                 : sampleId,
                    sample             : sampleId,
                    platform           : plat,
                    type               : row.type   ?: 'UNKNOWN',
                    trim               : row.trim?.toString().toLowerCase()=='true',
                    sequencing_summary : row.sequencing_summary ?: null,
                    single_end         : (f2==null),
                    directory          : false,
                    batch              : true,
                    needscompressing   : needsCompress
                ]
                def fileList = [ file(f1) ] + (f2 ? [ file(f2) ] : [])
                results << [ meta, fileList ]
            }
        }
    } else {
        // ─── single‐sample mode: emit exactly one tuple ──────────────────────────
        def f1 = file(row.fastq_1)
        def f2 = row.fastq_2 ? file(row.fastq_2) : null
        if( f2 && ! f2.exists() ) {
            error "ERROR: Read2 does not exist: ${f2}"
        }
        def meta = [
            id                 : row.sample ?: 'Unknown_Sample',
            platform           : row.platform ?: 'ILLUMINA',
            type               : row.type ?: 'UNKNOWN',
            trim               : row.trim.toBoolean() ?: false,
            sequencing_summary : row.sequencing_summary ?: null,
            single_end         : row.single_end.toBoolean() ?: false,
            directory          : row.directory.toBoolean() ?: false,
            batch              : row.batch.toBoolean() ?: false,
            needscompressing   : row.needscompressing.toBoolean() ?: false

        ]


        def fileList = [ f1 ] + (f2 ? [ f2 ] : [])
        results << [ meta, fileList ]
    }

    return results
}
