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

// NF v26 does not support Groovy `import` declarations; fully-qualified names are used inline below.

// (All former top-level statements moved inside workflow TAXTRIAGE below — NF v25+ requirement)

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
include { INSILICO } from '../subworkflows/local/insilico'
include { PROTEINS } from '../subworkflows/local/proteins'
include { NOVELTY } from '../subworkflows/local/novelty'
// Shared MicrobeRT clustering, lifted up to the workflow level so its output can feed BOTH
// the MicrobeRT classifier (in REPORT) and the Pyrodigal->mmseqs taxonomy novelty branch.
include { CLUSTER_ALIGNMENT } from '../modules/local/cluster_alignment'
include { MMSEQS_EASYCLUSTER } from '../modules/local/mmseqs2_easycluster'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { DOWNLOAD_DB } from '../modules/local/download_db'
include { MMSEQS_DOWNLOADDB } from '../modules/local/mmseqs_downloaddb'
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

workflow TAXTRIAGE {
    // ── Initialisation (moved from top level for NF v25+ compatibility) ──────────
    def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
    WorkflowTaxtriage.initialise(params, log)

    def checkPathParamList = [ params.input ]

    if (workflow.containerEngine != 'singularity' && workflow.containerEngine != 'docker') {
        exit 1, "Neither Docker or Singularity was selected as the container engine. Please specify with `-profile docker` or `-profile singularity`. Exiting..."
    }

    if (params.classifier != 'kraken2' && params.classifier != 'centrifuge' && params.classifier != 'metaphlan') {
        exit 1, "Classifier must be either kraken2, centrifuge or metaphlan"
    }

    println "Working Directory: ${workflow.workDir}"

    if (params.fastq_1) {
        if (!file(params.fastq_1).exists()) {
            exit 1, "ERROR: fastq_1 file does not exist: ${params.fastq_1}"
        }
        if (params.fastq_2) {
            if (!file(params.fastq_2).exists()) {
                exit 1, "ERROR: fastq_2 file does not exist: ${params.fastq_2}"
            }
        }
    } else if (params.input) {
        if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not available or non-existent!' }
    } else {
        exit 1, 'ERROR: Please specify either an input samplesheet (--input) or at least a fastq_1 file (--fastq_1)!'
    }

    if (params.minq) {
        ch_minq_shortreads = params.minq
        ch_minq_longreads  = params.minq
    } else {
        ch_minq_shortreads = 20
        ch_minq_longreads  = 7
        println 'Min Quality set to default'
    }

    ch_save_fastq_classified = params.save_classified_fastq ? true : false
    ch_assembly_txt          = null
    ch_kraken_reference      = false

    String  value   = 'G,-10,-2'
    boolean matches = value.matches('^(G|L),-?\\d+(\\.\\d+)?,-?\\d+(\\.\\d+)?$')
    ch_empty_file   = file("$projectDir/assets/NO_FILE")

    if (matches) { println('The value matches the pattern.') }
    else          { println('The value does not match the pattern.') }

    // Require Kraken2 DB unless Kraken2 is skipped
    if (!params.skip_kraken2 && !params.db && !params.download_db) {
        exit 1, "If --skip_kraken2 is false, you must provide --db or --download_db"
    }

    if (params.skip_kraken2 && !params.reference_fasta && !params.get_pathogens && !params.organisms && !params.organisms_file) {
        exit 1, "If you are skipping kraken2, you must provide a reference fasta, --get_pathogens to pull the pathogens file, organisms, or organisms_file"
    }

    ch_pathogens = Channel.fromPath("$projectDir/assets/pathogen_sheet.csv", checkIfExists: true)
    if (params.pathogens) {
        if (params.pathogens.endsWith('.csv') || params.pathogens.endsWith('.txt')) {
            ch_pathogens = Channel.fromPath(params.pathogens, checkIfExists: true)
        } else {
            exit 1, "Pathogens file must end with .csv or .txt i.e. it is a .csv (comma-delimited) file!"
        }
    }

    if (!params.assembly) {
        println 'No assembly file given, downloading the standard NCBI RefSeq summary' + (params.enable_genbank ? ' and GenBank summary (--enable_genbank)' : ' (GenBank pulling disabled; enable with --enable_genbank)')
        ch_assembly_txt = null
    } else {
        println "Assembly file present, using it to pull genomes from... ${params.assembly}"
        def _assembly_files = [file(params.assembly, checkIfExists: true)]
        if (params.assembly_summary_genbank) {
            println "GenBank assembly file also provided: ${params.assembly_summary_genbank}"
            _assembly_files << file(params.assembly_summary_genbank, checkIfExists: true)
        }
        ch_assembly_txt = _assembly_files.size() == 1 ? _assembly_files[0] : _assembly_files
    }

    if (!params.assembly_file_type) {
        ch_assembly_file_type = 'ncbi'
    } else {
        ch_assembly_file_type = params.assembly_file_type
    }

    workflow_summary          = WorkflowTaxtriage.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary       = Channel.value(workflow_summary)
    ch_multiqc_config         = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_css            = Channel.fromPath("$projectDir/assets/mqc.css", checkIfExists: true)
    ch_multiqc_custom_config  = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo           = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo,   checkIfExists: true) : Channel.empty()
    ch_multiqc_files          = Channel.empty()
    ch_multiqc_files          = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_merged_table_config    = Channel.fromPath("$projectDir/assets/table_explanation_mqc.yml", checkIfExists: true)
    // ── End initialisation ───────────────────────────────────────────────────────

    // Info required for completion email and summary
    def multiqc_report = []

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
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260226.tar.gz',
            'checksum': '34cff72c6bd67a0892709c7b472c8351',
            'size': '7.5G'
        ],
        'standard': [
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20260226.tar.gz',
            'checksum': 'bbff202fcbc7280e8ea22077e57e3b25',
            'size': '78G'
        ],
        'viral': [
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20260226.tar.gz',
            'checksum': 'cd3c624c6fd774ea71dc164350b83783',
            'size': '553M'
        ],
        'pluspf': [
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20260226.tar.gz',
            'checksum': '37294169dc75fd21da324bb99e1e0d85',
            'size': '77G'
        ],
        'pluspfp16': [
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16_GB_20260226.tar.gz',
            'checksum': 'cf80ea5ad50b3b6276132fda43bfe714',
            'size': '16G'
        ],
        'pluspf8': [
            'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08_GB_20260226.tar.gz',
            'checksum': '4a0e671232c516b4d017a5165288ede0',
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
                java.nio.file.Files.walk(path)
                    .filter { p ->
                        java.nio.file.Files.isRegularFile(p) &&
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
    ch_taxdump_nodes = ch_empty_file
    ch_taxdump_dir = ch_empty_file.parent

    if (!params.taxdump && (params.download_taxdump || (!params.taxdump && params.metaphlan) ) ){
        println "No taxdump provided, downloading the latest taxdump from NCBI"
        DOWNLOAD_TAXDUMP()
        ch_taxdump_nodes = DOWNLOAD_TAXDUMP.out.nodes
        ch_taxdump_dir = DOWNLOAD_TAXDUMP.out.nodes.parent
    }
    else if (params.taxdump){
        // set ch_taxdump_nodes BUT add nodes.dmp to end as a file
        ch_taxdump_nodes = file("$params.taxdump/nodes.dmp", checkIfExists: true)
        ch_taxdump_dir = file(params.taxdump, checkIfExists: true)
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
    ch_meta = Channel.empty()
    ch_reads = INPUT_CHECK.out.reads
    ch_meta = INPUT_CHECK.out.ch_meta
    // ── Run-level metadata: --meta CSV ───────────────────────────────────────
    // Supports two formats:
    //   a) CSV with 'sample' column → each row maps to one sample by name
    //   b) CSV with only 'run_id' column → applied to all samples whose
    //      meta.run_id (from the samplesheet) matches the run_id value
    // Fields merged into meta: run_id, latitude, longitude, depth, salinity,
    //   collection_time, location.  Samplesheet columns (if present) are
    //   used as defaults; --meta CSV values override them.
    // ── Meta CSV channel (passed directly to ALIGNMENT_PER_SAMPLE via REPORT) ─
    // The CSV file is passed as a Nextflow file input to the process, where
    // match_paths.py reads it directly — no Groovy channel-combination needed.
    if (params.meta) {
        // Use .first() to convert the single-file queue channel into a value channel
        // so it is broadcast to ALL ALIGNMENT_PER_SAMPLE invocations, not just one.
        ch_meta_csv = Channel.fromPath(params.meta, checkIfExists: true).first()
    } else {
        // Use the meta channel emitted by INPUT_CHECK (from samplesheet extra columns)
        // if it produced a file; otherwise fall back to the empty sentinel.
        ch_meta_csv = ch_meta
            .ifEmpty(file("$projectDir/assets/NO_FILE_meta_csv"))
    }

    // Apply --positive / --negative CLI param overrides to force a sample as a control
    if (params.negative || params.positive) {
        ch_reads = ch_reads.map { meta, reads ->
            if (params.negative && meta.id == params.negative) {
                meta.control = true
                meta.control_type = 'negative'
            }
            if (params.positive && meta.id == params.positive) {
                meta.control = true
                meta.control_type = 'positive'
            }
            [meta, reads]
        }
    }

    // ── FASTA inputs: branch off before any QC / trimming / host-removal ────────
    // Samples whose fastq_1 column contains a FASTA file (.fa/.fasta/.fna or .gz
    // variants) set meta.is_fasta = true in input_check.nf.  They skip every
    // read-quality step and are rejoined into the pipeline just before the
    // classifier so that Kraken2 and alignment still run on them normally.
    ch_reads.branch {
        fasta: it[0].is_fasta == true
        fastq: !(it[0].is_fasta == true)
    }.set { reads_by_type }
    ch_fasta_reads = reads_by_type.fasta
    ch_reads       = reads_by_type.fastq

    // ch_pass_files tracks every sample (FASTQ + FASTA) through to REPORT
    ch_pass_files = ch_reads.mix(ch_fasta_reads).map{ meta, reads -> [ meta ] }

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
    ch_reads = ch_reads
        .map { meta, reads -> [meta.id, meta, reads] }
        .join(readCountChannel.map { meta, countFile -> [meta.id, countFile] }, by: 0)
        .map { id, meta, reads, countFile ->
            def count = countFile.text.trim().toInteger()
            meta.read_count = count
            return [meta, reads]
        }


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
    // Force singleton removal when de novo assembly or diamond is enabled,
    // as singletons can cause issues with assemblers
    if ((params.use_denovo || params.use_diamond) && params.include_singletons_removal) {
        println "WARNING: --include_singletons_removal has been overwritten to false because --use_denovo or --use_diamond was specified. Singletons will be removed from paired-end host-removed reads."
        params.include_singletons_removal = false
    }

    // Re-join FASTA inputs into the main read channel now so they benefit from
    // host removal.  minimap2 handles FASTA queries natively; if host removal
    // runs, REMOVE_HOSTREADS converts them to FASTQ via samtools-fastq (dummy
    // quality scores).  If host removal is not configured, FASTA files pass
    // through unchanged.  Either way, meta.is_fasta stays true and gates the
    // QC visualisation steps below.
    ch_fasta_reads = ch_fasta_reads.map { meta, reads ->
        meta.read_count = 0   // COUNT_READS was skipped for these samples
        [meta, reads]
    }
    ch_reads = ch_reads.mix(ch_fasta_reads)

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
    ch_multiqc_files = ch_multiqc_files.mix(HOST_REMOVAL.out.host_removal_stats)

    //////////////////// RUN OPTIONAL FASTQC to get qc plots for multiqc  ////////////////////
    // Skip QC plots for FASTA inputs: if host removal ran, those reads were
    // converted to FASTQ with dummy quality scores by samtools-fastq, so the
    // plots would be meaningless.  If host removal did not run, the files are
    // still FASTA and FastQC / NanoPlot would error.
    if (!params.skip_plots) {
        FASTQC(
            ch_reads.filter { it[0].platform =~ /(?i)ILLUMINA/ && !it[0].is_fasta }
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        NANOPLOT(
            ch_reads.filter { (it[0].platform =~ /(?i)OXFORD/ || it[0].platform =~ /(?i)PACBIO/) && !it[0].is_fasta }
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] }.ifEmpty([]) )
        ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect { it[1] }.ifEmpty([]) )
    }
    // ch_reads now contains both processed FASTQ samples and FASTA samples
    // (the latter having passed through host-removal but skipped all QC plots)
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
    ch_mapped_assemblies = ch_preppedfiles.map{ meta, fastas, map, gcfids -> [meta, map] }
    ch_accessions = Channel.empty()
    ch_bedfiles = Channel.empty()
    ch_bedfiles_or_default = Channel.empty()
    ch_alignment_stats = Channel.empty()
    ch_assembly_analysis = Channel.empty()
    ch_fastas = Channel.empty()

    ////////////////////////////// OPTIONAL: IN-SILICO READ SIMULATION //////////////////////////////
    // Simulate reads (ISS for Illumina, NanoSim for ONT) from Kraken2 abundance profiles.
    // Simulated reads are emitted as new samples tagged with insilico meta.  They are then
    // injected into the normal ALIGNMENT → REPORT pipeline.  In REPORT, insilico samples are
    // split out, processed through ALIGNMENT_PER_SAMPLE, and their JSONs are used as insilico
    // controls for the real (non-control) samples.
    ch_insilico_reads = Channel.empty()
    if (params.generate_iss || params.generate_nanosim) {
        // Extract merged_taxid map from REFERENCE_PREP prepped files (non-controls only)
        ch_sim_merged_taxid = ch_preppedfiles
            .filter { !it[0].control }
            .map { meta, fastas, mergedmap, mergedids ->
                [meta, mergedmap]
            }

        // Extract reference FASTAs from REFERENCE_PREP (non-controls only)
        ch_sim_fastas = REFERENCE_PREP.out.fastas
            .filter { !it[0].control }

        // Top hits from CLASSIFIER (non-controls only)
        ch_sim_tops = CLASSIFIER.out.ch_tops
            .filter { !it[0].control }

        INSILICO(
            ch_sim_tops,
            ch_sim_merged_taxid,
            ch_sim_fastas
        )
        ch_versions = ch_versions.mix(INSILICO.out.versions)
        ch_insilico_reads = INSILICO.out.insilico_reads

        // ── Inject insilico reads as new samples into the pipeline channels ──
        // Mix insilico reads into ch_reads so they flow through ALIGNMENT
        ch_reads = ch_reads.mix(ch_insilico_reads)

        // Create ch_preppedfiles entries for insilico samples by cloning the
        // parent sample's reference prep data with the insilico meta.
        // INSILICO.out.insilico_reads has tuple(insilico_meta, reads) where
        // insilico_meta.parent_id == original sample id.
        ch_insilico_prepfiles = ch_insilico_reads
            .map { meta, reads -> [meta.parent_id, meta] }
            .combine(
                ch_preppedfiles.map { meta, fastas, map, gcfids -> [meta.id, fastas, map, gcfids] },
                by: 0
            )
            .map { parent_id, insilico_meta, fastas, map, gcfids ->
                [insilico_meta, fastas, map, gcfids]
            }
        ch_preppedfiles = ch_preppedfiles.mix(ch_insilico_prepfiles)

        // Update ch_mapped_assemblies to include insilico entries
        ch_mapped_assemblies = ch_mapped_assemblies.mix(
            ch_insilico_prepfiles.map { meta, fastas, map, gcfids -> [meta, map] }
        )

        // Create placeholder ch_kraken2_report entries for insilico samples
        ch_kraken2_report = ch_kraken2_report.mix(
            ch_insilico_reads.map { meta, reads ->
                [meta, file("$projectDir/assets/NO_FILE")]
            }
        )

        // Add insilico entries to ch_fastas (REFERENCE_PREP.out.fastas)
        ch_insilico_fastas = ch_insilico_reads
            .map { meta, reads -> [meta.parent_id, meta] }
            .combine(
                REFERENCE_PREP.out.fastas.map { meta, fastas -> [meta.id, fastas] },
                by: 0
            )
            .map { parent_id, insilico_meta, fastas ->
                [insilico_meta, fastas]
            }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////
    if (!params.skip_realignment) {
        ch_prepfiles = ch_reads.join(ch_preppedfiles.map{ meta, fastas, map, gcfids -> [meta, fastas, map] })
        ////////////////////////////////////////////////////////////////////////////////////////////////
        ALIGNMENT(
            ch_prepfiles
        )
        ch_postalignmentfiles = ch_reads.map{ meta, reads -> [meta, null, null, null, null, null, null, []] }
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

        // Add insilico fastas if simulation was run
        if (params.generate_iss || params.generate_nanosim) {
            ch_fastas = ch_fastas.mix(ch_insilico_fastas)
        }

        ch_postalignmentfiles = ch_combined.map {
            meta, bam, bai, mapping ->  return [ meta, bam, bai, mapping ]
        }.filter{
            it[1]
        }
        // Insilico samples will be naturally filtered out here because they
        // don't have entries in ch_bedfiles, ch_reference_cds, etc.
        // This is correct — insilico samples skip ASSEMBLY.
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
            ch_filtered_reads.map{ meta, reads -> [meta, reads] }
        )
        ////////////////////////////////////////////////////////////////////////////////////////////////
        ASSEMBLY(
            ch_postalignmentfiles,
            ch_assembly_txt
        )
        ch_diamond_output = ASSEMBLY.out.ch_diamond_output
        ch_versions = ch_versions.mix(ASSEMBLY.out.versions)

        ch_assembly_analysis = ASSEMBLY.out.ch_diamond_analysis
        ch_denovo = ASSEMBLY.out.ch_denovo_assembly

        // Seed a per-sample placeholder for EVERY sample coming out of ALIGNMENT
        // (covers insilico and control samples that skip ASSEMBLY/PROTEINS)
        ch_annotate_report_tsv = ALIGNMENT.out.bams.map { meta, bam, bai ->
            [meta, file("$projectDir/assets/NO_FILE_annotate_report")]
        }

        if (params.annotate) {
            // if params.annotate_proteins is null then set to 'assets/bvbrc_specialty_genes_with_sequences_taxids_and_sites.faa'
            PROTEINS(
                ch_denovo,
            )


            // Overlay the real per-sample xlsx where available; keep placeholder elsewhere
            ch_annotate_report_tsv = ch_annotate_report_tsv
                .join(PROTEINS.out.annotate_report, remainder: true)
                .map { meta, placeholder, real_file ->
                    [meta, real_file ?: placeholder]
            }

        }

        ////////////////////////////////////////////////////////////////////////////////////////////////
        // SHARED MicrobeRT CLUSTERING
        // One CLUSTER_ALIGNMENT -> MMSEQS_EASYCLUSTER pass, computed here (before both NOVELTY and
        // REPORT) so its representatives feed BOTH consumers:
        //   * REPORT  -> MICROBERT_PREDICT / MICROBERT_PARSE   (the MicrobeRT classifier)
        //   * NOVELTY -> PYRODIGAL -> MMSEQS_TAXONOMY          (the novelty taxonomy branch)
        // Only runs under --microbert; otherwise the channels stay empty and NOVELTY falls back
        // to its de-novo-contig path while REPORT emits the MicrobeRT placeholder.
        ////////////////////////////////////////////////////////////////////////////////////////////////
        ch_microbert_reps     = Channel.empty()
        ch_microbert_clusters = Channel.empty()
        if (params.microbert) {
            CLUSTER_ALIGNMENT(
                ALIGNMENT.out.bams.map { meta, bam, csi -> [meta, bam] }
            )
            MMSEQS_EASYCLUSTER(
                CLUSTER_ALIGNMENT.out.fasta
            )
            ch_microbert_reps     = MMSEQS_EASYCLUSTER.out.representatives
            ch_microbert_clusters = MMSEQS_EASYCLUSTER.out.tsv
            ch_versions = ch_versions.mix(MMSEQS_EASYCLUSTER.out.versions)
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////
        // NOVELTY: reference-free / open-set detection on the closed-set residual
        // (reads that aligned to NO reference + de novo contigs) via translated-search LCA.
        ////////////////////////////////////////////////////////////////////////////////////////////////
        // Per-sample novelty outputs, surfaced to REPORT so they can be collected into the
        // interactive HTML report (panel + downloadable JSON/XLSX). Empty unless --detect_novelty.
        ch_novelty_summary    = Channel.empty()
        ch_novelty_candidates = Channel.empty()
        // Embedding files (*.umap.tsv, *.cluster_summary.tsv, *.clusters.tsv) from NOVEL_HOMOLOGS.
        // Collected into a flat list so CREATE_COMPARISON_REPORT can stage them all at once.
        // Populated inside the detect_novelty block below when a NOVEL_HOMOLOGS module is wired up.
        ch_embedding_files = Channel.value(file("$projectDir/assets/NO_FILE_embedding"))
        if (params.detect_novelty) {
            // Resolve the seqTaxDB local-first (like --db); otherwise download + cache it
            // once via `mmseqs databases` (storeDir-cached, reused across runs and on -resume).
            def ch_novelty_db
            if (params.novelty_db && file(params.novelty_db).exists()) {
                println "Novelty: using local mmseqs seqTaxDB at ${params.novelty_db}"
                ch_novelty_db = Channel.value(file(params.novelty_db, checkIfExists: true))
            } else {
                def novelty_dbname = params.novelty_db ?: 'UniProtKB'
                println "Novelty: seqTaxDB '${novelty_dbname}' not found locally; will download via " +
                        "'mmseqs databases' (cached at ${params.novelty_db_cache})."
                MMSEQS_DOWNLOADDB(novelty_dbname)
                ch_versions = ch_versions.mix(MMSEQS_DOWNLOADDB.out.versions)
                ch_novelty_db = MMSEQS_DOWNLOADDB.out.db.first()
            }

            NOVELTY(
                ALIGNMENT.out.bams,     // [meta, bam, csi] reference-merged per sample
                ch_kraken2_report,      // [meta, kreport]  from CLASSIFIER
                ch_denovo,              // [meta, contigs]  ASSEMBLY.out.ch_denovo_assembly
                ch_novelty_db,
                ch_microbert_reps       // [meta, rep_seq.fasta] shared MicrobeRT cluster reps (empty unless --microbert)
            )
            ch_versions = ch_versions.mix(NOVELTY.out.versions)
            ch_novelty_summary    = NOVELTY.out.summary       // [meta, *.novelty.summary.tsv]
            ch_novelty_candidates = NOVELTY.out.candidates    // [meta, *.novelty.candidates.tsv]
            // TODO: when NOVEL_HOMOLOGS gets a proper NF module, replace the placeholder below with:
            //   ch_embedding_files = NOVELTY.out.embedding_files.flatten().collect()
        }

        // Add placeholder assembly analysis entries for insilico samples so
        // they are not filtered out by the inner join in input_alignment_files
        if (params.generate_iss || params.generate_nanosim) {
            ch_assembly_analysis = ch_assembly_analysis.mix(
                ch_insilico_reads.map { meta, reads ->
                    [meta, file("$projectDir/assets/NO_FILE2")]
                }
            )
        }

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
                .join(ch_annotate_report_tsv)

            ////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////
            all_samples = ch_pass_files.map{ it[0].id }.collect().flatten().toSortedList()


            REPORT(
                input_alignment_files,
                ch_pathogens,
                distributions,
                ch_assembly_txt,
                ch_taxdump_dir,
                all_samples,
                ch_meta_csv,
                ch_microbert_reps,      // [meta, rep_seq.fasta] shared MicrobeRT cluster reps
                ch_microbert_clusters,  // [meta, *.tsv]          shared MicrobeRT cluster membership
                ch_novelty_summary,     // [meta, *.novelty.summary.tsv]
                ch_novelty_candidates,  // [meta, *.novelty.candidates.tsv]
                ch_embedding_files      // flat: *.umap.tsv, *.cluster_summary.tsv, *.clusters.tsv
            )
            ch_multiqc_files = ch_multiqc_files.mix(REPORT.out.merged_report_txt.collect { it }.ifEmpty([]))
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////
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

    // ── Completion handler (moved inside workflow for NF v25+ compatibility) ──
    // Capture references explicitly to avoid delegate-resolution surprises in onComplete.
    def wf_meta             = workflow
    def run_params          = params
    def run_summary_params  = summary_params
    def run_multiqc_report  = multiqc_report

    workflow.onComplete {
        if (run_params?.email || run_params?.email_on_fail) {
            NfcoreTemplate.email(wf_meta, run_params, run_summary_params, projectDir, log, run_multiqc_report)
        }
        NfcoreTemplate.summary(wf_meta, run_params, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
