/* groovylint-disable DuplicateMapLiteral */
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowTaxtriage.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

if (params.minq) {
    ch_minq_illumina = params.minq
    ch_minq_oxford = params.minq
} else {
    ch_minq_illumina = 20
    ch_minq_oxford = 7
    println 'Min Quality set to default'
}

ch_assembly_txt = null
ch_kraken_reference = false

if (!params.assembly) {
    println 'No assembly file given, downloading the standard ncbi one'
    ch_assembly_txt = null
} else {
    println "Assembly file present, using it to pull genomes from... ${params.assembly}"
    ch_assembly_txt = file(params.assembly, checkIfExists: true)
}

if (!params.assembly_file_type) {
    ch_assembly_file_type = 'ncbi'
} else {
    ch_assembly_file_type = params.assembly_file_type
}

if (params.assembly && ch_assembly_txt.isEmpty()) {
    exit 1, "File provided with --assembly is empty: ${ch_assembly_txt.getName()}!"
}  else if (params.assembly) {
    println "Assembly file present, using it to pull genomes from ncbi... ${params.assembly}"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_css       = file("$projectDir/assets/mqc.css", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
ch_merged_table_config        = Channel.fromPath("$projectDir/assets/table_explanation_mqc.yml", checkIfExists: true)
ch_alignment_stats = Channel.empty()
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { ALIGNMENT } from '../subworkflows/local/alignment'
include { READSFILTER } from '../subworkflows/local/filter_reads'
include { KRONA_KTUPDATETAXONOMY  } from '../modules/nf-core/krona/ktupdatetaxonomy/main'
include { KRONA_KTIMPORTTEXT  } from '../modules/nf-core/krona/ktimporttext/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { DOWNLOAD_DB } from '../modules/local/download_db'
include { DOWNLOAD_TAXTAB } from '../modules/local/download_taxtab'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { PYCOQC                      } from '../modules/nf-core/pycoqc/main'
include { FASTP } from '../modules/nf-core/fastp/main'
include { KRAKEN2_KRAKEN2                      } from '../modules/nf-core/kraken2/kraken2/main'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'
include { ARTIC_GUPPYPLEX } from '../modules/nf-core/artic/guppyplex/main'
include { MOVE_FILES } from '../modules/local/moveFiles.nf'
include { MOVE_NANOPLOT } from '../modules/local/move_nanoplot.nf'
include { PORECHOP } from '../modules/nf-core/porechop/main'
include { SEQTK_SAMPLE } from '../modules/nf-core/seqtk/sample/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { FLYE                     } from '../modules/nf-core/flye/main'
include { KRAKENTOOLS_COMBINEKREPORTS   } from '../modules/nf-core/krakentools/combinekreports/main'
include { KRONA   } from '../modules/local/krona.nf'
include { SPADES as SPADES_ILLUMINA } from '../modules/nf-core/spades/main'
include { SPADES as SPADES_OXFORD } from '../modules/nf-core/spades/main'
include { NANOPLOT                     } from '../modules/nf-core/nanoplot/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CONFIDENCE_METRIC } from '../modules/local/confidence'
include { CONVERT_CONFIDENCE } from '../modules/local/convert_confidence'
include { PULL_TAXID } from '../modules/local/pull_taxid'
include { REFERENCE } from '../modules/local/download_reference'
include { PULL_FASTA } from '../modules/local/pullFASTA'
include { TOP_HITS } from '../modules/local/top_hits'
include { GET_ASSEMBLIES } from '../modules/local/get_assembly_refs'
include { REMOVETAXIDSCLASSIFICATION } from '../modules/local/remove_taxids.nf'
include { KRAKENREPORT } from '../modules/local/krakenreport'
include { MERGEDKRAKENREPORT } from '../modules/local/merged_krakenreport'
include { FILTERKRAKEN } from '../modules/local/filter_krakenreport'
include { MERGE_CONFIDENCE } from '../modules/local/merge_confidence'
include { VISUALIZE_REPORTS } from '../subworkflows/local/visualize_reports'
include { HOST_REMOVAL } from '../subworkflows/local/host_removal'
include { KREPORT_TO_KRONATXT } from '../modules/local/generate_krona_txtfile'
include { NCBIGENOMEDOWNLOAD }  from '../modules/nf-core/ncbigenomedownload/main'
include { DOWNLOAD_ASSEMBLY } from '../modules/local/download_assembly'
include { NCBIGENOMEDOWNLOAD_FEATURES } from '../modules/local/get_feature_tables'
include { FEATURES_TO_BED } from '../modules/local/convert_features_to_bed'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow TAXTRIAGE {
    supported_dbs = [
        'flukraken2': [
            'url': 'https://media.githubusercontent.com/media/jhuapl-bio/mytax/master/databases/flukraken2.tar.gz',
            'checksum': '9d388703b1fa7c2e269bb63acf1043dbec7bb62da0a57c4fb1c41d8ab7f9c953',
            'size': '180M'
        ],
        'minikraken2': [
            'url': 'ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz',
            'checksum': 'a184ae5c1e382abfff34574e135ceaaace4ac27605b205f4fb83dca11cfa42ac',
            'size': '7.5G'
        ],
        'standard8': [
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz',
            'checksum': 'a184ae5c1e382abfff34574e135ceaaace4ac27605b205f4fb83dca11cfa42ac',
            'size': '7.5G'
        ],
        'viral': [
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20230605.tar.gz',
            'checksum': 'adf5deba8a62f995609592aa86e2f7aac7e49162e995e132a765b96edb456f99',
            'size': '553M'
        ],
        'test': [
            'url': 'https://github.com/jhuapl-bio/datasets/raw/main/databases/kraken2/test_metagenome.tar.gz',
            'checksum': 'c7d50ca4f46885ce7d342a06e748f9390cf3f4157a54c995d34ecdabbc83e1b8',
            'size': '112M'
        ],

    ]
    // if the download_db params is called AND the --db is not existient as a path
    // then download the db
    if (params.download_db) {
        if (supported_dbs.containsKey(params.db)) {
            println "Kraken db ${params.db} will be downloaded if it cannot be found. This requires ${supported_dbs[params.db]['size']} of space."
            DOWNLOAD_DB(
                params.db,
                supported_dbs[params.db]['url'],
                supported_dbs[params.db]['checksum']
            )
            /* groovylint-disable-next-line UnnecessaryGetter */
            ch_db = DOWNLOAD_DB.out.k2d.map { file -> file.getParent() }

            println('_____________')
        } else if (params.db)  {
            ch_db = file(params.db, checkIfExists: true)
        } else {
            println "Database ${params.db} not found in download list. Currently supported databases are ${supported_dbs.keySet()}. If this database has already been downloaded, indicate it with --db <exact path>"
        }
    } else {
        if (params.db) {
            file(params.db, checkIfExists: true)
            ch_db = params.db
        }
    }

    if (!ch_assembly_txt) {
        println 'empty'
        GET_ASSEMBLIES()
        GET_ASSEMBLIES.out.assembly.map {  record -> record }.set { ch_assembly_txt }
    }

    ch_versions = Channel.empty()
    ch_mergedtsv = Channel.empty()
    // make an empty path channel
    ch_accession_mapping  = Channel.empty()
    ch_empty_file = file("$projectDir/assets/NO_FILE")

    ch_hit_to_kraken_report = Channel.empty()
    // // // //
    // // // // MODULE: MultiQC
    // // // //
    workflow_summary    = WorkflowTaxtriage.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_css))
    // // // //
    // // // //
    // // // //

    // //
    // // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    // //
    INPUT_CHECK(
        ch_input
    )
    ch_reads = INPUT_CHECK.out.reads
    ARTIC_GUPPYPLEX(
        ch_reads.filter { it[0].directory   }
    )
    ch_reads = ARTIC_GUPPYPLEX.out.fastq
    ch_reads = ch_reads.mix(INPUT_CHECK.out.reads.filter { !it[0].directory   })

    if (params.subsample && params.subsample > 0) {
        ch_subsample  = params.subsample
        SEQTK_SAMPLE(
            ch_reads,
            ch_subsample
        )
        ch_reads = SEQTK_SAMPLE.out.reads
    }

    PYCOQC(
        ch_reads.filter { it[0].platform == 'OXFORD' && it[0].sequencing_summary != null }.map {
            meta, reads -> meta.sequencing_summary
}
    )

    // // // //

    // // // // MODULE: Run FastQC or Porechop, Trimgalore
    // // //
    ch_porechop_out = Channel.empty()
    if (params.trim) {
        nontrimmed_reads = ch_reads.filter { !it[0].trim }
        TRIMGALORE(
            ch_reads.filter { it[0].platform == 'ILLUMINA' && it[0].trim }
        )
        PORECHOP(
            ch_reads.filter { it[0].platform == 'OXFORD' && it[0].trim  }
        )
        ch_porechop_out  = PORECHOP.out.reads

        trimmed_reads = TRIMGALORE.out.reads.mix(PORECHOP.out.reads)
        ch_reads = nontrimmed_reads.mix(trimmed_reads)
    }
    ch_fastp_reads = Channel.empty()
    ch_fastp_html = Channel.empty()
    if (!params.skip_fastp) {
        FASTP(
            ch_reads,
            [],
            false,
            false
        )
        ch_reads = FASTP.out.reads
        ch_fastp_reads = FASTP.out.json
        ch_fastp_html = FASTP.out.html
    }

    HOST_REMOVAL(
        ch_reads,
        params.genome
    )
    ch_reads = HOST_REMOVAL.out.unclassified_reads
    ch_multiqc_files = ch_multiqc_files.mix(HOST_REMOVAL.out.stats_filtered)

    if (!params.skip_plots) {
        FASTQC(
            ch_reads.filter { it[0].platform =~ /(?i)ILLUMINA/ }
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        NANOPLOT(
            ch_reads.filter { it[0].platform  =~ /(?i)OXFORD/ }
        )

    }

    // // // // // //
    // // // // // // MODULE: Run Kraken2
    // // // // // //
    KRAKEN2_KRAKEN2(
        ch_reads,
        ch_db,
        true,
        true
    )

    ch_kraken2_report = KRAKEN2_KRAKEN2.out.report

    KREPORT_TO_KRONATXT(
        ch_kraken2_report
    )
    ch_krona_txt = KREPORT_TO_KRONATXT.out.txt
    ch_combined = ch_krona_txt
                .map{ it[1] }        // Get the file path
                .collect()            // Collect all file parts into a list
                .map { files ->
                    // Join the files with single quotes and space
                    String joinedFiles = files.collect { "'$it'" }.join(' ')
                    [[id:'combined_krona_kreports'], joinedFiles]  // Combine with new ID
                }

    KRONA_KTIMPORTTEXT(
        ch_combined
    )

    if (params.remove_taxids) {
        remove_input = ch_kraken2_report.map {
            meta, report -> [
                meta, report, params.remove_taxids
            ]
        }
        REMOVETAXIDSCLASSIFICATION(
            remove_input
        )
        ch_kraken2_report = REMOVETAXIDSCLASSIFICATION.out.report
    }

    TOP_HITS(
        ch_kraken2_report
    )
    MERGEDKRAKENREPORT(
        TOP_HITS.out.krakenreport.map { meta, file ->  file }.collect()
    )

    FILTERKRAKEN(
        MERGEDKRAKENREPORT.out.krakenreport
    )
    ch_filtered_reads = KRAKEN2_KRAKEN2.out.classified_reads_fastq.map { m, r-> [m, r.findAll { it =~ /.*\.classified.*(fq|fastq)(\.gz)?/  }] }
    ch_hit_to_kraken_report = TOP_HITS.out.tops.join(ch_filtered_reads)
    ch_accessions = Channel.empty()
    ch_bedfiles = Channel.empty()
    ch_hit_to_kraken_report_merged = Channel.empty()
    ch_bedfiles_or_default = Channel.empty()
    if (params.reference_fasta) { //
        // format of the FASTA file MUST be "kraken:taxid|<taxidnumber>" in each reference accession
        ch_reference_fasta = params.reference_fasta ? Channel.fromPath(params.reference_fasta, checkIfExists: true) : Channel.empty()
        ch_hit_to_kraken_report = ch_hit_to_kraken_report.combine(
            ch_reference_fasta
        )
    } else {
        DOWNLOAD_ASSEMBLY(
            ch_hit_to_kraken_report.map {
                meta, report, readsclass ->  [ meta, report ]
            },
            ch_assembly_txt
        )
        ch_hit_to_kraken_report = ch_hit_to_kraken_report.join(
            DOWNLOAD_ASSEMBLY.out.fasta
        )
        ch_accessions = DOWNLOAD_ASSEMBLY.out.accessions
        ch_accession_mapping = DOWNLOAD_ASSEMBLY.out.mappings

        if (!params.skip_features){
            NCBIGENOMEDOWNLOAD_FEATURES(
                ch_accessions
            )

            FEATURES_TO_BED(
                NCBIGENOMEDOWNLOAD_FEATURES.out.features
            )

            ch_bedfiles = FEATURES_TO_BED.out.bed
        }
    }


    if (!params.skip_realignment) {
        // Combine ch_hit_to_kraken_report with ch_bedfiles_or_default
        if (params.skip_features){
            ch_hit_to_kraken_report = ch_hit_to_kraken_report.map {
                meta, report, readsclass, fasta -> [ meta, report, readsclass, fasta, null ]
            }
        } else {
            ch_hit_to_kraken_report = ch_hit_to_kraken_report.join(ch_bedfiles)
        }

        ALIGNMENT(
            ch_hit_to_kraken_report
        )

        ch_alignment_stats = ALIGNMENT.out.stats
        ch_bamstats = ALIGNMENT.out.bamstats
        ch_depth = ALIGNMENT.out.depth

        ch_alignment_outmerg = ALIGNMENT.out.bams.join(ALIGNMENT.out.depth)
        // ch_alignment_outmerg = ch_alignment_outmerg.join(Channel.empty(), remainder: true)

        ch_combined = ch_alignment_outmerg
            .join(ch_accession_mapping, by: 0, remainder: true)
            .map { meta, bam, depth, mapping ->
                // If mapping is not present, replace it with null or an empty placeholder
                return [meta, bam, depth, mapping ?: ch_empty_file]
            }
        CONFIDENCE_METRIC(
           ch_combined
        )

        ch_joined_confidence_report = KRAKEN2_KRAKEN2.out.report.join(
            CONFIDENCE_METRIC.out.tsv
        )


        if (!params.skip_confidence) {
            CONVERT_CONFIDENCE(
                ch_joined_confidence_report
            )

            MERGE_CONFIDENCE(
                CONVERT_CONFIDENCE.out.tsv.map {  file ->  file }.collect()
            )
            ch_mergedtsv = MERGE_CONFIDENCE.out.confidence_report
        }
    }

    if (!params.skip_assembly) {
        illumina_reads = ch_hit_to_kraken_report.filter {
            it[0].platform == 'ILLUMINA'
        }.map {
            meta, reads -> [meta, reads, [], []]
        }
        SPADES_ILLUMINA(
            illumina_reads
        )

        println('____________________________')

        nanopore_reads = ch_hit_to_kraken_report.filter {
            it[0].platform == 'OXFORD'
        }.map {
            meta, reads -> [meta, [], [], [], reads]
        }

        FLYE(
            nanopore_reads,
            '--nano-raw'
        )
    }

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    // // //
    // // // MODULE: MultiQC Pt 2
    // // //

    ch_multiqc_files = ch_multiqc_files.mix(MERGEDKRAKENREPORT.out.krakenreport.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FILTERKRAKEN.out.reports.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TOP_HITS.out.krakenreport.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_alignment_stats.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_porechop_out.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_merged_table_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_fastp_html.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_kraken2_report.collect { it[1] }.ifEmpty([]))
    if (params.trim) {
        ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.reads.collect { it[1] }.ifEmpty([]))
    }
    if (!params.skip_plots) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect { it[1] }.ifEmpty([]))
    }
    ch_multiqc_files = ch_multiqc_files.mix(ch_mergedtsv.collect().ifEmpty([]))

    // Unused or Incomplete
    // if (params.blastdb && !params.remoteblast){
    //     ch_multiqc_files = ch_multiqc_files.mix(BLAST_BLASTN.out.txt.collect{it[1]}.ifEmpty([]))
    // } else if (params.blastdb && params.remoteblast){
    //     ch_multiqc_files = ch_multiqc_files.mix(REMOTE_BLASTN.out.txt.collect{it[1]}.ifEmpty([]))
    // }
    if (!params.skip_multiqc){
        MULTIQC(
            ch_multiqc_files.collect()
        )
        multiqc_report = MULTIQC.out.report.toList()
        ch_versions    = ch_versions.mix(MULTIQC.out.versions)
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
