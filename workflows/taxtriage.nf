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

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.FileVisitOption

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


// check that the params.classifiers is either kraken2 or centrifuge or metaphlan4
if (params.classifier != 'kraken2' && params.classifier != 'centrifuge' && params.classifier != 'metaphlan') {
    exit 1, "Classifier must be either kraken2, centrifuge or metaphlan"
}

println "Working Directory: ${workflow.workDir}"

// Either use an existing samplesheet file or build one from fastq parameters
if (params.fastq_1) {
    // check that fastq_1 exists
    if (!file(params.fastq_1).exists()) {
        exit 1, "ERROR: fastq_1 file does not exist: ${params.fastq_1}"
    }
    if (params.fastq_2){
        // check that fastq_2 exists
        if (!file(params.fastq_2).exists()) {
            exit 1, "ERROR: fastq_2 file does not exist: ${params.fastq_2}"
        }
    }
} else if (params.input) {
    if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not available or non-existent!' }
} else {
    ex
    it 1, 'ERROR: Please specify either an input samplesheet (--input) or at least a fastq_1 file (--fastq_1)!'
}

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
ch_save_fastq_classified = params.save_classified_fastq ? true : false
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
ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_css       = Channel.fromPath("$projectDir/assets/mqc.css", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
ch_multiqc_files = Channel.empty()
ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

ch_merged_table_config        = Channel.fromPath("$projectDir/assets/table_explanation_mqc.yml", checkIfExists: true)

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
include { COUNT_READS  } from '../modules/local/count_reads'
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
include { METRIC_ALIGNMENT } from '../modules/local/confidence'
include { CONVERT_METRICS } from '../modules/local/convert_confidence'
include { PULL_TAXID } from '../modules/local/pull_taxid'
include { REFERENCE } from '../modules/local/download_reference'
include { GET_ASSEMBLIES } from '../modules/local/get_assembly_refs'
include { PULL_FASTA } from '../modules/local/pullFASTA'
include { MERGE_ALIGNMENT_MERGES } from '../modules/local/merge_confidence'
include { NCBIGENOMEDOWNLOAD }  from '../modules/nf-core/ncbigenomedownload/main'
include { NCBIGENOMEDOWNLOAD_FEATURES } from '../modules/local/get_feature_tables'
include { METRIC_MERGE } from '../modules/local/merge_confidence_contigs'
include { MAP_GCF } from '../modules/local/map_gcfs'
include { REFERENCE_REHEADER } from '../modules/local/reheader'
include { BBMAP_BBNORM } from '../modules/nf-core/bbmap/bbnorm/main'

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
    // ch_reference_fasta = params.reference_fasta ? Channel.fromPath(params.reference_fasta, checkIfExists: true) : Channel.empty()

    ch_reference_fasta = params.reference_fasta ? Channel.from(params.reference_fasta.split(" ").collect { it  }) : Channel.empty()
    ch_reference_fasta
    .flatMap { fasta ->
        // Replace tilde with home directory
        def normalizedPath = fasta.replaceFirst('^~', System.getProperty('user.home'))
        // Convert to Path object
        def path = file(normalizedPath)

        if (path.isDirectory()) {
            // Use Files.walk to traverse the directory recursively
            def fastaFiles = []
            if (params.recursive_reference){
                Files.walk(path)
                    .filter { p ->
                        Files.isRegularFile(p) &&
                        (p.fileName.toString().toLowerCase().endsWith('.fa') || p.fileName.toString().toLowerCase().endsWith('.fasta'))
                    }
                    .forEach { p -> fastaFiles << file(p.toString()) } // Use Nextflow's 'file' for consistency
            } else {
                def filesArray = path.listFiles()
                if (!filesArray) {
                    println "Warning: The directory '${normalizedPath}' is empty or inaccessible."
                    return []
                }
                fastaFiles = filesArray.toList().findAll { file ->
                    file.name.toLowerCase().endsWith('.fa') || file.name.toLowerCase().endsWith('.fasta')
                }
            }
            if (fastaFiles.isEmpty()) {
                println "Warning: No .fa or .fasta files found in directory '${normalizedPath}'."
                return []
            }
            // Return the list of files
            return fastaFiles
        } else if (path.isFile()) {
            // Return the file itself as a single-item list
            return [path]
        } else {
            // Warn about invalid paths
            println "Warning: The path '${normalizedPath}' is not a valid file or directory and will be skipped."
            return []
        }
    }
    .set { ch_reference_fasta }

    if (params.get_pathogens){
        DOWNLOAD_PATHOGENS()
        ch_reference_fasta = DOWNLOAD_PATHOGENS.out.fasta
    }
    ch_taxdump_dir = Channel.empty()
    ch_taxdump_nodes = Channel.empty()
    ch_taxdump_names = Channel.empty()

    if (params.download_taxdump || (!params.taxdump && params.metaphlan)){
        DOWNLOAD_TAXDUMP()
        ch_taxdump_nodes = DOWNLOAD_TAXDUMP.out.nodes
        ch_taxdump_dir = DOWNLOAD_TAXDUMP.out.nodes.parent
    }
    else if (params.taxdump){
        // set ch_taxdump_nodes BUT add nodes.dmp to end as a file
        ch_taxdump_nodes = file("$params.taxdump/nodes.dmp", checkIfExists: true)
        ch_taxdump_dir = file(params.taxdump, checkIfExists: true)
    } else {
        // set to empty_file
        ch_taxdump_nodes = ch_empty_file
        ch_taxdump_dir = ch_empty_file.parent
    }

    // if the download_db params is called AND the --db is not existient as a path
    // then download the db
    ch_db = Channel.empty()
    if (params.download_db && !params.skip_kraken2) {
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
        if (params.db  && !params.skip_kraken2) {
            file(params.db, checkIfExists: true)
            ch_db = params.db
        }
    }

    if (params.classifiers){
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
    INPUT_CHECK()

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


    // compress reads if needed

    ch_reads.branch {
        needsCompress: it[0].needscompressing
        noCompress: !it[0].needscompressing
    }.set{ split_compressing }

    PIGZ_COMPRESS(
        split_compressing.needsCompress
    )
    ch_reads = split_compressing.noCompress.mix(PIGZ_COMPRESS.out.archive)

    // // // //
    // // // // MODULE: Run FastQC or Porechop, Trimgalore
    // // //
    ch_porechop_out = Channel.empty()
    ch_fastp_reads = Channel.empty()
    ch_fastp_html = Channel.empty()
    ch_mapaa_taxid = Channel.empty()
    ch_diamond_output = Channel.empty()
    params.pathogens ? ch_pathogens = Channel.fromPath(params.pathogens, checkIfExists: true) : ''

    nontrimmed_reads = ch_reads.filter { !it[0].trim }
    TRIMGALORE(
        ch_reads.filter { it[0].platform == 'ILLUMINA' && it[0].trim }
    )

    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.reads.collect { it[1] }.ifEmpty([]) )

    PORECHOP(
        ch_reads.filter { (it[0].platform == 'OXFORD' || it[0].platform == "PACBIO") && it[0].trim  }
    )
    ch_porechop_out  = PORECHOP.out.reads
    trimmed_reads = TRIMGALORE.out.reads.mix(PORECHOP.out.reads)
    ch_reads = nontrimmed_reads.mix(trimmed_reads)
    ch_multiqc_files = ch_multiqc_files.mix(ch_porechop_out.collect { it[1] }.ifEmpty([]))
    // Create an empty file if se_reads is null
    // When calling the module, pass the empty file instead of null:
    if (params.downsample) {
        BBMAP_BBNORM(
            ch_reads
        )
        ch_reads = BBMAP_BBNORM.out.fastq
        ch_versions = ch_versions.mix(BBMAP_BBNORM.out.versions)
    }
    COUNT_READS(ch_reads)
    readCountChannel = COUNT_READS.out.count
    // Update the meta with the read count by reading the file content
    readCountChannel.map { meta, countFile, reads ->
        def count = countFile.text.trim().toInteger()
        // Update meta map by adding a new key 'read_count'
        meta.read_count = count
        return [meta, reads]
    }.set{ ch_reads }

    //////////////////// RUN PYCOQC on any seq summary file ////////////////////

    PYCOQC(
        ch_reads.filter { it[0].platform == 'OXFORD' && it[0].sequencing_summary != null }.map {
            meta, reads -> meta.sequencing_summary
        }
    )
    //////////////////// RUN FASTP to get qc plots and output reads ////////////////////
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
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { it[1] }.ifEmpty([]))
    }

    //////////////////// RUN ALIGNEMNT to filter out host reads ////////////////////
    HOST_REMOVAL(
        ch_reads,
        params.genome
    )
    //////////////////// RUN OPTIONAL SEQTK to subsample arbitrarily ////////////////////

    ch_reads = HOST_REMOVAL.out.unclassified_reads
    if (params.subsample && params.subsample > 0) {
        ch_subsample  = params.subsample
        SEQTK_SAMPLE(
            ch_reads,
            ch_subsample
        )
        ch_reads = SEQTK_SAMPLE.out.reads
    }
    // test to make sure that fastq files are not empty files
    ch_multiqc_files = ch_multiqc_files.mix(HOST_REMOVAL.out.stats_filtered)

    //////////////////// RUN OPTIONAL FASTQC to get qc plots for multiqc  ////////////////////
    if (!params.skip_plots) {
        FASTQC(
            ch_reads.filter { it[0].platform =~ /(?i)ILLUMINA/ }
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        NANOPLOT(
            ch_reads.filter { it[0].platform =~ /(?i)OXFORD/ || it[0].platform =~ /(?i)PACBIO/ }
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] }.ifEmpty([]) )
        ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect { it[1] }.ifEmpty([]) )
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
    //////////////////////////////RUN CLASSIFIER(S) for top hits calculations//////////////////////////////////////////////////////////////////
    CLASSIFIER(
        ch_filtered_reads,
        ch_db,
        ch_save_fastq_classified,
        distributions,
        ch_pathogens,
        ch_organisms_to_download,
        ch_taxdump_dir
    )
    ch_kraken2_report = CLASSIFIER.out.ch_kraken2_report
    ch_reads = CLASSIFIER.out.ch_reads
    ch_krona = CLASSIFIER.out.ch_krona_plot
    // ch_multiqc_files = ch_multiqc_files.mix(ch_krona.collect { it[1] }.ifEmpty([]))
    ch_krakenreport = CLASSIFIER.out.ch_tops
    ch_pass_files = ch_pass_files.join(ch_kraken2_report)
    // add ch_kraken2_report to ch_multiqc, only unique names
    ch_multiqc_files = ch_multiqc_files.mix(
    ch_kraken2_report
            .map { it[1] } // Correctly map to the second element of each tuple
            .ifEmpty(Channel.empty()) // Handle empty channels appropriately
            .distinct() // Remove duplicates if necessary
    )
    ch_organisms_to_download = CLASSIFIER.out.ch_organisms_to_download
    ch_multiqc_files = ch_multiqc_files.mix(ch_krakenreport.collect { it[1] }.ifEmpty([]))

    ///////////////////////////RUN REFERENCE pull or processing/////////////////////////////////////////////////////////////////////
    REFERENCE_PREP(
        ch_organisms_to_download,
        ch_reference_fasta,
        ch_assembly_txt,
        ch_pathogens

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
    ch_fastas = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////
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
        ch_covfiles = ALIGNMENT.out.stats
        ch_alignment_stats = ALIGNMENT.out.stats
        ch_bedgraphs = ALIGNMENT.out.bedgraphs
        ch_versions = ch_versions.mix(ALIGNMENT.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(ch_alignment_stats.collect { it[1] }.ifEmpty([]))

        ch_alignment_outmerg = ALIGNMENT.out.bams

        ch_alignment_outmerg
            .join(ch_mapped_assemblies, by: 0, remainder: true)
            .filter{
                it[1]
            }
            .map { meta, bam, bai, mapping ->
                // If mapping is not present, replace it with null or an empty placeholder
                return [meta, bam, bai, mapping ?: ch_empty_file]
            }.set{ ch_combined }

        ch_bedfiles = REFERENCE_PREP.out.ch_bedfiles
        ch_fastas = REFERENCE_PREP.out.fastas

        ch_postalignmentfiles = ch_combined.map {
            meta, bam, bai, mapping ->  return [ meta, bam, bai, mapping ]
        }.filter{
            it[1]
        }
        ch_postalignmentfiles = ch_combined.map {
            meta, bam, bai, mapping ->  return [ meta, bam, bai, mapping ]
        }.filter{
            it[1]
        }
        .join(ch_bedfiles)
        .join(REFERENCE_PREP.out.ch_reference_cds)
        .join(REFERENCE_PREP.out.features)
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
        ch_versions = ch_versions.mix(ASSEMBLY.out.versions)

        ch_assembly_analysis = ASSEMBLY.out.ch_diamond_analysis
        ch_assembly_analysis_opt = ch_assembly_analysis.ifEmpty {
            Channel.value(null)
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////


        if (!params.skip_report){
            // if ch_kraken2_report is empty join on empty
            // Define a channel that emits a placeholder value if ch_kraken2_report is empty
            input_alignment_files = ALIGNMENT.out.bams
                .join(ch_mapped_assemblies)
                .join(ch_bedgraphs)
                .join(ch_covfiles)
                .join(ch_kraken2_report)
                .join(ch_assembly_analysis)
                .join(ch_fastas)

            ////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////
            all_samples = ch_pass_files.map{ it[0].id }.collect().flatten().toSortedList()

            REPORT(
                input_alignment_files,
                ch_pathogens,
                distributions,
                ch_assembly_txt,
                ch_taxdump_nodes,
                all_samples
            )
            ch_multiqc_files = ch_multiqc_files.mix(REPORT.out.merged_report_txt.collect { it }.ifEmpty([]))
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////

        // if (!params.skip_confidence) {
        //     METRIC_ALIGNMENT(
        //         ch_combined
        //     )

        //     METRIC_MERGE(
        //         METRIC_ALIGNMENT.out.tsv
        //     )
        //     CONVERT_METRICS(
        //         METRIC_MERGE.out.metrics
        //     )

        //     MERGE_ALIGNMENT_MERGES(
        //         CONVERT_METRICS.out.tsv.map {  file ->  file }.collect()
        //     )

        //     ch_mergedtsv = MERGE_ALIGNMENT_MERGES.out.metrics_report
        //     ch_multiqc_files = ch_multiqc_files.mix(ch_mergedtsv.collect().ifEmpty([]))
        // }
    }
    ch_collated_versions = ch_versions.unique().collectFile(name: 'all_mqc_versions.yml')
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_collated_versions
    )
    // // //
    // // // MODULE: MultiQC Pt 2
    // // //
    // Unused or Incomplete

    if (!params.skip_multiqc){
        MULTIQC(
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_multiqc_logo.toList(),
            [],
            []
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
