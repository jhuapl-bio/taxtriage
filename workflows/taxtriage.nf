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
// println "Initialising taxtriage workflow with parameters: ${workflow}"
// def wfsummary = NfcoreSchema.generateJSONSchema(workflow)

// print wfsummary to json file to working directory
// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
]

if (workflow.containerEngine !== 'singularity' && workflow.containerEngine !== 'docker'){
    exit 1 , "Neither Docker or Singularity was selected as the container engine. Please specify with `-profile docker` or `-profile singularity`. Exiting..."
}

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// check that the params.classifiers is either kraken2 or centrifuge or metaphlan4
if (params.classifier != 'kraken2' && params.classifier != 'centrifuge' && params.classifier != 'metaphlan') {
    exit 1, "Classifier must be either kraken2, centrifuge or metaphlan"
}



// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

if (params.minq) {
    ch_minq_shortreads = params.minq
    ch_minq_longreads = params.minq
} else {
    ch_minq_shortreads = 20
    ch_minq_longreads = 7
    println 'Min Quality set to default'
}

// if params.save_fastq_classified then set ch_save_fastq_classified to true
// else set ch_save_fastq_classified to false
ch_save_fastq_classified = params.skip_classified_fastq ? false : true
ch_assembly_txt = null
ch_kraken_reference = false


def validateBt2Scoremin(String scoremin) {
    def pattern = /^(G|L),-?\d+(\.\d+)?,-?\d+(\.\d+)?$/
    if (!scoremin || !pattern.matcher(scoremin).matches()) {
        error "ERROR: The parameter 'bt2_scoremin' is in an incorrect format or not provided. It should be 'G' or 'L' followed by two comma-separated values (e.g., 'G,-10,-2')."
    }
}
String value = "G,-10,-2"
boolean matches = value.matches('^(G|L),-?\\d+(\\.\\d+)?,-?\\d+(\\.\\d+)?$')
ch_empty_file = file("$projectDir/assets/NO_FILE")

if (matches) {
    println("The value matches the pattern.")
} else {
    println("The value does not match the pattern.")
}
// if (params.bt2_scoremin) {
//     // Call the validation function early in the script
//     validateBt2Scoremin(params.bt2_scoremin)
// }

// if skip_kraken2 and reference_fasta is empty AND organisms is empty and organisms_file is empty print and exit that organisms is required
if (params.skip_kraken2 && !params.reference_fasta && !params.get_pathogens && !params.organisms && !params.organisms_file) {
    exit 1, "If you are skipping kraken2, you must provide a reference fasta, --get_pathogens to pull the pathogens file, organisms, or organisms_file"
}

// if params.pathogens, check if file ends with .tsv or .txt
if (params.pathogens) {
    if (params.pathogens.endsWith('.csv') || params.pathogens.endsWith('.txt')) {
        ch_pathogens = Channel.fromPath(params.pathogens, checkIfExists: true)
    } else {
        exit 1, "Pathogens file must end with .csv or .txt i.e. it is a .csv (comma-delimited) file!"
    }
} else {
    ch_pathogens = ch_empty_file
}
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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// // // //
// // // // MODULE: MultiQC
// // // //
workflow_summary    = WorkflowTaxtriage.paramsSummaryMultiqc(workflow, summary_params)
ch_workflow_summary = Channel.value(workflow_summary)
ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_css       = file("$projectDir/assets/mqc.css", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
ch_multiqc_files = Channel.empty()
ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_css))
ch_merged_table_config        = Channel.fromPath("$projectDir/assets/table_explanation_mqc.yml", checkIfExists: true)
ch_multiqc_files = ch_multiqc_files.mix(ch_merged_table_config.collect().ifEmpty([]))

// // // //
// // // // MODULE:  Pathogens
// // // //
//// Get Pathogen sheet by default
ch_pathogens = Channel.fromPath("$projectDir/assets/pathogen_sheet.csv", checkIfExists: true)
// // // //
// // // //

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
include { REPORT } from '../subworkflows/local/report'
include { READSFILTER } from '../subworkflows/local/filter_reads'
include { HOST_REMOVAL } from '../subworkflows/local/host_removal'
include { REFERENCE_PREP } from '../subworkflows/local/reference_prep'
include { ASSEMBLY } from '../subworkflows/local/assembly'
include { CLASSIFIER } from '../subworkflows/local/classifier'

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
include { DOWNLOAD_PATHOGENS } from '../modules/local/download_pathogens'
include { DOWNLOAD_TAXDUMP } from '../modules/local/download_taxdump'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { PYCOQC                      } from '../modules/nf-core/pycoqc/main'
include { PIGZ_COMPRESS } from '../modules/nf-core/pigz/compress/main'
include { FASTP } from '../modules/nf-core/fastp/main'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'
include { ARTIC_GUPPYPLEX } from '../modules/nf-core/artic/guppyplex/main'
include { MOVE_FILES } from '../modules/local/moveFiles.nf'
include { MOVE_NANOPLOT } from '../modules/local/move_nanoplot.nf'
include { PORECHOP } from '../modules/nf-core/porechop/main'
include { SEQTK_SAMPLE } from '../modules/nf-core/seqtk/sample/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { NANOPLOT                     } from '../modules/nf-core/nanoplot/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CONFIDENCE_METRIC } from '../modules/local/confidence'
include { CONVERT_CONFIDENCE } from '../modules/local/convert_confidence'
include { PULL_TAXID } from '../modules/local/pull_taxid'
include { REFERENCE } from '../modules/local/download_reference'
include { GET_ASSEMBLIES } from '../modules/local/get_assembly_refs'
include { PULL_FASTA } from '../modules/local/pullFASTA'
include { MERGE_CONFIDENCE } from '../modules/local/merge_confidence'
include { NCBIGENOMEDOWNLOAD }  from '../modules/nf-core/ncbigenomedownload/main'
include { NCBIGENOMEDOWNLOAD_FEATURES } from '../modules/local/get_feature_tables'
include { CONFIDENCE_MERGE } from '../modules/local/merge_confidence_contigs'
include { MAP_GCF } from '../modules/local/map_gcfs'
include {  REFERENCE_REHEADER } from '../modules/local/reheader'

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
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240605.tar.gz',
            'checksum': 'a184ae5c1e382abfff34574e135ceaaace4ac27605b205f4fb83dca11cfa42ac',
            'size': '7.5G'
        ],
        'standard': [
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240605.tar.gz',
            'checksum': 'a184ae5c1e382abfff34574e135ceaaace4ac27605b205f4fb83dca11cfa42ac',
            'size': '78G'
        ],
        'viral': [
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240605.tar.gz',
            'checksum': 'adf5deba8a62f995609592aa86e2f7aac7e49162e995e132a765b96edb456f99',
            'size': '553M'
        ],
        'pluspf': [
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240605.tar.gz',
            'checksum': 'adf5deba8a62f995609592aa86e2f7aac7e49162e995e132a765b96edb456f99',
            'size': '77G'
        ],
        'pluspfp16': [
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20240605.tar.gz',
            'checksum': 'adf5deba8a62f995609592aa86e2f7aac7e49162e995e132a765b96edb456f99',
            'size': '16G'
        ],
        'pluspf8': [
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08gb_20240605.tar.gz',
            'checksum': 'adf5deba8a62f995609592aa86e2f7aac7e49162e995e132a765b96edb456f99',
            'size': '7.5G'
        ],
        'eupath': [
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_eupathdb48_20230407.tar.gz',
            'checksum': 'adf5deba8a62f995609592aa86e2f7aac7e49162e995e132a765b96edb456f99',
            'size': '11G'
        ],
        'test': [
            'url': 'https://github.com/jhuapl-bio/datasets/raw/main/databases/kraken2/test_metagenome.tar.gz',
            'checksum': 'c7d50ca4f46885ce7d342a06e748f9390cf3f4157a54c995d34ecdabbc83e1b8',
            'size': '112M'
        ],

    ]
    ch_reference_fasta = params.reference_fasta ? Channel.fromPath(params.reference_fasta, checkIfExists: true) : Channel.empty()

    if (params.get_pathogens){
        DOWNLOAD_PATHOGENS()
        ch_reference_fasta = DOWNLOAD_PATHOGENS.out.fasta
    }
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
        } else if (params.db)  {
            ch_db = file(params.db, checkIfExists: true)
        } else {
            println "Database ${params.db} not found in download list. Currently supported databases are ${supported_dbs.keySet()}. If this database has already been downloaded, indicate it with --db <exact path>. You may also retrieve them, locally, from https://benlangmead.github.io/aws-indexes/k2. Make sure to download and decompress/untar the .tar.gz file which contains a folder you can specify with --db <localpath>. "
        }
    } else {
        if (params.db) {
            file(params.db, checkIfExists: true)
            ch_db = params.db
        }
    }

    ch_taxdump_dir = Channel.empty()
    if (params.classifier){
        // split params.classifier on command and optional space assign to list channel
        ch_classifier = params.classifiers.split(",\\s*")
    } else {
        ch_classifier = ['kraken2']
    }
    if (!ch_assembly_txt) {
        GET_ASSEMBLIES()
        GET_ASSEMBLIES.out.assembly.map {  record -> record }.set { ch_assembly_txt }
    }

    ch_versions = Channel.empty()
    ch_mergedtsv = Channel.empty()
    // make an empty path channel
    ch_accession_mapping  = Channel.empty()
    ch_organisms = Channel.empty()
    ch_pass_files = Channel.empty()
    // // // //
    // // // //
    // // // //

    // //
    // // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    // //
    INPUT_CHECK(
        ch_input
    )

    // Example nextflow.config file or within the script
    println "Nextflow version: ${workflow.nextflow.version}"
    // if commitId is null then set it to "local"
    if (!workflow.commitId) {
        workflow.commitId = 'local'
    }
    println "Commit ID: ${workflow.commitId}"
    if (!workflow.repository) {
        workflow.repository = 'local'
    }
    println "Repository URL: ${workflow.repository}"
    // get the date and time of the run
    def date = new Date()
    println "Date: ${date}"

    ch_reads = INPUT_CHECK.out.reads
    ch_pass_files = ch_reads.map{ meta, reads -> {
            return [ meta ]
        }
    }

    ARTIC_GUPPYPLEX(
        ch_reads.filter { it[0].directory   }
    )
    ch_reads = ch_reads.filter({ !it[0].directory   }).mix(ARTIC_GUPPYPLEX.out.fastq)

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
    ch_fastp_reads = Channel.empty()
    ch_fastp_html = Channel.empty()
    ch_mapaa_taxid = Channel.empty()
    ch_diamond_output = Channel.empty()
    params.pathogens ? ch_pathogens = Channel.fromPath(params.pathogens, checkIfExists: true) : ''


    // if (params.trim) {
    nontrimmed_reads = ch_reads.filter { !it[0].trim }
    TRIMGALORE(
        ch_reads.filter { it[0].platform == 'ILLUMINA' && it[0].trim }
    )

    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.reads.collect { it[1] }.ifEmpty([]))

    PORECHOP(
        ch_reads.filter { (it[0].platform == 'OXFORD' || it[0].platform == "PACBIO") && it[0].trim  }
    )
    ch_porechop_out  = PORECHOP.out.reads
    trimmed_reads = TRIMGALORE.out.reads.mix(PORECHOP.out.reads)
    ch_reads = nontrimmed_reads.mix(trimmed_reads)
    ch_multiqc_files = ch_multiqc_files.mix(ch_porechop_out.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_fastp_html.collect { it[1] }.ifEmpty([]))
//
    // }
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

    // test to make sure that fastq files are not empty files
    ch_multiqc_files = ch_multiqc_files.mix(HOST_REMOVAL.out.stats_filtered)

    if (!params.skip_plots) {
        FASTQC(
            ch_reads.filter { it[0].platform =~ /(?i)ILLUMINA/ }
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        NANOPLOT(
            ch_reads.filter { it[0].platform =~ /(?i)OXFORD/ || it[0].platform =~ /(?i)PACBIO/ }
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect { it[1] }.ifEmpty([]))
    }
    ch_filtered_reads = ch_reads
    ch_profile = Channel.empty()
    ch_preppedfiles = Channel.empty()
    ch_organisms_to_download = ch_filtered_reads.map { meta, reads -> return [meta, []] }

    def empty_organism_file = false
    if (params.unknown_sample){
        distributions = Channel.fromPath(ch_empty_file)
    } else if (!params.distributions){
        distributions = Channel.fromPath("$projectDir/assets/taxid_abundance_stats.hmp.tsv.gz", checkIfExists: true)
    } else{
        distributions = Channel.fromPath(params.distributions)
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    CLASSIFIER(
        ch_filtered_reads,
        ch_db,
        ch_save_fastq_classified,
        distributions,
        ch_pathogens,
        ch_organisms_to_download
    )
    ch_kraken2_report = CLASSIFIER.out.ch_kraken2_report
    ch_reads = CLASSIFIER.out.ch_reads
    ch_pass_files = ch_pass_files.join(ch_kraken2_report)
    ch_organisms_to_download = CLASSIFIER.out.ch_organisms_to_download
    ////////////////////////////////////////////////////////////////////////////////////////////////
    REFERENCE_PREP(
        ch_organisms_to_download,
        ch_reference_fasta,
        ch_assembly_txt
    )
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // todo - add alignment for contigs
    ch_mapped_assemblies = Channel.empty()
    ch_preppedfiles = REFERENCE_PREP.out.ch_preppedfiles
    ch_mapped_assemblies = ch_preppedfiles.map{
        meta, fastas, map, gcfids -> {
            return [meta, map]
        }
    }
    ch_accessions = Channel.empty()
    ch_bedfiles = Channel.empty()
    ch_bedfiles_or_default = Channel.empty()
    ch_alignment_stats = Channel.empty()
    ch_assembly_analysis = Channel.empty()

    // If you use a local genome Refseq FASTA file
    // if ch_refernece_fasta is empty

    if (!params.skip_realignment) {
        ch_prepfiles = ch_reads.join(ch_preppedfiles.map{ meta, fastas, map, gcfids -> {
                return [meta, fastas, map]
            }
        })
        ////////////////////////////////////////////////////////////////////////////////////////////////
        ALIGNMENT(
            ch_prepfiles
        )
        ch_postalignmentfiles = ch_reads.map{
            meta, reads -> {
                return [meta, null, null, null, null, null, null,  []]
            }
        }
        ch_depthfiles = ALIGNMENT.out.depth
        ch_covfiles = ALIGNMENT.out.stats
        ch_alignment_stats = ALIGNMENT.out.stats
        ch_multiqc_files = ch_multiqc_files.mix(ch_alignment_stats.collect { it[1] }.ifEmpty([]))

        ch_depth = ALIGNMENT.out.depth

        ch_alignment_outmerg = ALIGNMENT.out.bams.join(ALIGNMENT.out.depth)

        ch_combined = ch_alignment_outmerg
            .join(ch_mapped_assemblies, by: 0, remainder: true)
            .map { meta, bam, bai, depth, mapping ->
                // If mapping is not present, replace it with null or an empty placeholder
                return [meta, bam, bai, depth, mapping ?: ch_empty_file]
            }
        ch_bedfiles = REFERENCE_PREP.out.ch_bedfiles


        ch_postalignmentfiles = ch_combined.map {
            meta, bam, bai,  depth, mapping ->  return [ meta, bam, bai, mapping ]
        }.join(ch_bedfiles)
        .join(REFERENCE_PREP.out.ch_reference_cds)
        .join(REFERENCE_PREP.out.ch_cds_to_taxids)
        .join(
            ch_filtered_reads.map{
                meta, reads -> {
                    return [meta, reads]
                }
            }
        )

        ////////////////////////////////////////////////////////////////////////////////////////////////
        ASSEMBLY(
            ch_postalignmentfiles,
            ch_assembly_txt
        )
        ch_diamond_output = ASSEMBLY.out.ch_diamond_output
        ch_assembly_analysis = ASSEMBLY.out.ch_diamond_analysis
        // ch_assembly_alignment = ASSEMBLY.out.ch_assembly_alignment
        ////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if (!params.skip_report){
            // if ch_kraken2_report is empty join on empty
            // Define a channel that emits a placeholder value if ch_kraken2_report is empty

            input_alignment_files = ALIGNMENT.out.bams
                .join(ch_mapped_assemblies)
                .join(ch_depthfiles)
                .join(ch_covfiles)
                .join(ch_kraken2_report)
                .join(ch_assembly_analysis)
            REPORT(
                input_alignment_files,
                ch_pathogens,
                distributions,
                ch_assembly_txt
            )
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////

        if (!params.skip_confidence) {
            CONFIDENCE_METRIC(
                ch_combined
            )

            CONFIDENCE_MERGE(
                CONFIDENCE_METRIC.out.tsv
            )
            CONVERT_CONFIDENCE(
                CONFIDENCE_MERGE.out.confidence
            )

            MERGE_CONFIDENCE(
                CONVERT_CONFIDENCE.out.tsv.map {  file ->  file }.collect()
            )

            ch_mergedtsv = MERGE_CONFIDENCE.out.confidence_report
            ch_multiqc_files = ch_multiqc_files.mix(ch_mergedtsv.collect().ifEmpty([]))
        }
    }

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    // // //
    // // // MODULE: MultiQC Pt 2
    // // //
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
