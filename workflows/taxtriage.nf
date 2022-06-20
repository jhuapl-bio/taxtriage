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
// def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]

def checkPathParamList = [  params.db, params.reference ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }



if (params.assembly) {
    ch_assembly_txt = file(params.assembly, checkIfExists: true)
    ch_kraken_reference=false
    if (ch_assembly_txt.isEmpty()) {exit 1, "File provided with --assembly is empty: ${ch_assembly_txt.getName()}!"}
} else if (params.genomes) {
    ch_assembly_genomes = file(params.genomes, checkIfExists: true)
    ch_kraken_reference=true
    if (ch_assembly_genomes.isEmpty()) {exit 1, "File provided with --genomes is empty: ${ch_assembly_genomes.getName()}!"}
} else {
    ch_assembly_txt = false
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()



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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { KRAKEN2_KRAKEN2                      } from '../modules/nf-core/modules/kraken2/kraken2/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { CONFIDENCE_METRIC } from '../modules/local/confidence'
include { PULL_TAXID } from '../modules/local/pull_taxid'
include { REFERENCE } from '../modules/local/download_reference'
include { DOWNLOAD_ASSEMBLY } from '../modules/local/download_assembly'
include { PULL_FASTA } from '../modules/local/pullFASTA'
include { TOP_HITS } from '../modules/local/top_hits'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []



workflow TAXTRIAGE {
    ch_db = params.db
    
    ch_versions = Channel.empty()
    
    // //
    // // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    // //
    INPUT_CHECK (
        ch_input
    )
    INPUT_CHECK.out.reads.view()
    
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )


    // //
    // // MODULE: Run Kraken2
    // //
    KRAKEN2_KRAKEN2(
        INPUT_CHECK.out.reads,
        ch_db,
        params.save_output_fastqs,
        params.save_reads_assignment
    )
    // // ch_reference = Channel.fromPath(params.reference, type: 'file')

    // ch_bam_format = false
    // ch_cigar_paf_format =  true
    // ch_cigar_bam = true

    TOP_HITS (
        KRAKEN2_KRAKEN2.out.report
    )
    TOP_HITS.out.tops
        .splitCsv(header: true, sep: '\t')
        .map {
            meta, record -> [ meta,  record.taxid ]
        }.groupTuple( by: [0,0] ).first().set{ ch_hits }

    if (params.assembly){
        DOWNLOAD_ASSEMBLY (
            ch_assembly_txt
        )
        DOWNLOAD_ASSEMBLY.out.assembly.view()
        ch_assembly_txt_2 = Channel.fromPath(ch_assembly_txt)

        ch_assembly_txt_2
        .splitCsv( skip:1, sep: '\t', header: true  ).first()
        .view { row -> "${row.wgs_master} - ${row.taxid} - ${row.bioproject}" }
    } else {

        KRAKEN2_KRAKEN2.out.classified_reads_fastq
        .map { meta, record -> 
            record.findAll{ it =~ /.*\.classified.*(fq|fastq)(\.gz)?/  }
        }
        .set { classified_reads_fastq }

        KRAKEN2_KRAKEN2.out.classified_reads_assignment.map {
            meta, file  -> file
        }
        .set { ch_kraken2_assignments }
        

        PULL_FASTA (
            ch_hits,
            Channel.fromPath(ch_assembly_genomes),
            classified_reads_fastq,
            ch_kraken2_assignments
        )
    }
    ALIGNMENT(
        PULL_FASTA.out.fastq, 
        PULL_FASTA.out.fasta
    )
    


    // CONFIDENCE_METRIC (
    //     MINIMAP2_ALIGN.out.paf,
    //     PULL_TAXID.out.genome.first()
    // )




    
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )




    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowTaxtriage.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )

    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {
//     if (params.email || params.email_on_fail) {
//         NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//     }
//     NfcoreTemplate.summary(workflow, params, log)
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
