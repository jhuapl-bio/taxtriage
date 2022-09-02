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
// def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]

def checkPathParamList = [  params.db, params.reference ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.top_hits_count) { 
    ch_top_hits_count = params.top_hits_count 
} else { 
    ch_top_hits_count=2
    println 'Top hits not specified, defaulting to 10 per rank level in taxonomy tree for database for kraken2' 
}




ch_assembly_txt=null
ch_kraken_reference=false
if (!params.assembly){
    println "No assembly file given, downloading the standard ncbi one"
    ch_assembly_txt=null
} else {
    println "Assembly file present, using it to pull genomes from... ${params.assembly}"
    ch_assembly_txt=file(params.assembly, checkIfExists: true)
}
if (!params.assembly_file_type){
    ch_assembly_file_type = 'ncbi'
} else {
    ch_assembly_file_type = params.assembly_file_type
}
if (params.assembly && ch_assembly_txt.isEmpty() ) {
    exit 1, "File provided with --assembly is empty: ${ch_assembly_txt.getName()}!"
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
include { READSFILTER } from '../subworkflows/local/filter_reads'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { PYCOQC                      } from '../modules/nf-core/modules/pycoqc/main'
include { KRAKEN2_KRAKEN2                      } from '../modules/nf-core/modules/kraken2/kraken2/main'
include { TRIMGALORE } from '../modules/nf-core/modules/trimgalore/main'
include { ARTIC_GUPPYPLEX } from '../modules/nf-core/modules/artic/guppyplex/main'
include { MOVE_FILES } from '../modules/local/moveFiles.nf'
include { MOVE_NANOPLOT } from '../modules/local/move_nanoplot.nf'
include { PORECHOP } from '../modules/nf-core/modules/porechop/main'
include { SEQTK_SAMPLE } from '../modules/nf-core/modules/seqtk/sample/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { FLYE                     } from '../modules/nf-core/modules/flye/main'
include { SPADES as SPADES_ILLUMINA } from '../modules/nf-core/modules/spades/main'
include { SPADES as SPADES_OXFORD } from '../modules/nf-core/modules/spades/main'
include { NANOPLOT                     } from '../modules/nf-core/modules/nanoplot/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { CONFIDENCE_METRIC } from '../modules/local/confidence'
include { CONVERT_CONFIDENCE } from '../modules/local/convert_confidence'
include { PULL_TAXID } from '../modules/local/pull_taxid'
include { REFERENCE } from '../modules/local/download_reference'
include { DOWNLOAD_ASSEMBLY } from '../modules/local/download_assembly'
include { PULL_FASTA } from '../modules/local/pullFASTA'
include { TOP_HITS } from '../modules/local/top_hits'
include { GET_ASSEMBLIES } from '../modules/local/get_assembly_refs'
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
    ch_reads = INPUT_CHECK.out.reads
    if (params.subsample && params.subsample > 0){
        ch_subsample  = params.subsample
        SEQTK_SAMPLE (
            ch_reads,
            ch_subsample
        )
        ch_reads = SEQTK_SAMPLE.out.reads
        ch_reads.view()
    }
    if (params.demux){
        ARTIC_GUPPYPLEX(
            ch_reads.filter{ it[0].barcode }
        )
        ch_reads = ARTIC_GUPPYPLEX.out.fastq
        ch_reads = ch_reads.mix(INPUT_CHECK.out.reads.filter{ !it[0].barcode })
    } else {
        ch_reads = INPUT_CHECK.out.reads
    }

    PYCOQC(
        ch_reads.filter { it[0].platform == 'OXFORD' && it[0].sequencing_summary != null }.map{
            meta, reads -> meta.sequencing_summary
        }
    )
    if (!ch_assembly_txt){
        println "empty"
        GET_ASSEMBLIES(
            ch_reads
        )
        GET_ASSEMBLIES.out.assembly.map{ meta, record -> record }.set{ ch_assembly_txt }
        
    }
    // // //
    
    // // // MODULE: Run FastQC
    // // //
    if (params.trim){
        nontrimmed_reads = ch_reads.filter { !it[0].trim }
        TRIMGALORE(
            ch_reads.filter { it[0].platform == 'ILLUMINA' && it[0].trim }
        )
        PORECHOP(
            ch_reads.filter { it[0].platform == 'OXFORD' && it[0].trim  }
        )

        trimmed_reads = TRIMGALORE.out.reads.mix(PORECHOP.out.reads)
        ch_reads=nontrimmed_reads.mix(trimmed_reads)
    } 
    
    if (params.filter){
        ch_filter_db = file(params.filter)
        println "${ch_filter_db} <-- filtering reads on this db"
        READSFILTER(
            ch_reads,
            ch_filter_db
        )
        ch_reads = READSFILTER.out.reads
    }
    if (!params.skip_plots){
        FASTQC (
            ch_reads.filter { it[0].platform == 'ILLUMINA'}
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        NANOPLOT (
            ch_reads.filter { it[0].platform == 'OXFORD'}
        )
        ch_nanoplot_files_reformatted = NANOPLOT.out.html.map{
            meta, record -> [ meta, record.findAll{ !( it =~ /.*NanoPlot-report.html/) }  ]
        }
        MOVE_NANOPLOT(
            ch_nanoplot_files_reformatted
        )
    }
    
    // // // // //
    // // // // // MODULE: Run Kraken2
    // // // // //
    KRAKEN2_KRAKEN2(
        ch_reads,
        ch_db,
        true,
        true
    )
    TOP_HITS (
        KRAKEN2_KRAKEN2.out.report,
        ch_top_hits_count
    )
    if (!params.skip_realignment){
        ch_hit_to_kraken_report = TOP_HITS.out.tops.join(
            KRAKEN2_KRAKEN2.out.classified_reads_fastq
        )
        ch_hit_to_kraken_report = ch_hit_to_kraken_report.join(
            KRAKEN2_KRAKEN2.out.classified_reads_assignment
        )
    
        if (ch_assembly_file_type == 'ncbi' ){
            
            DOWNLOAD_ASSEMBLY (
                ch_hit_to_kraken_report,
                ch_assembly_txt
            )
            PULL_FASTA (
                DOWNLOAD_ASSEMBLY.out.fasta
            )
        } else {
            ch_hit_to_kraken_report = ch_hit_to_kraken_report.map{
                meta, report, classified_fastqs, reads_class -> [ meta, report, classified_fastqs, reads_class, ch_assembly_txt]
            }
            PULL_FASTA (
                ch_hit_to_kraken_report
            )
        }
        ALIGNMENT(
            PULL_FASTA.out.fastq
        )

        CONFIDENCE_METRIC (
            ALIGNMENT.out.sam,
        )
        ch_joined_confidence_report = KRAKEN2_KRAKEN2.out.report.join(
            CONFIDENCE_METRIC.out.tsv
        )
        CONVERT_CONFIDENCE (
            ch_joined_confidence_report
        )

        CONVERT_CONFIDENCE.out.tsv.collectFile(name: 'merged_mqc.tsv', keepHeader: true, storeDir: 'merged_mqc',  newLine: true)
        .set{ mergedtsv }
    }
    if (!params.spades_hmm ){
        ch_spades_hmm = []
    } else {
        ch_spades_hmm = params.spades_hmm
    }
     
    if (!params.skip_assembly){
        println "spades"
        illumina_reads =  PULL_FASTA.out.fastq.filter { it[0].platform == 'ILLUMINA'  }.map{
            meta, reads, ref -> 
            [meta, reads, [], []]
        }
        SPADES_ILLUMINA(
           illumina_reads,
           ch_spades_hmm
        )
        println "nanopore"
        nanopore_reads = PULL_FASTA.out.fastq.filter { it[0].platform == 'OXFORD'  }.map{
            meta, reads, ref -> 
            [meta, [], [], reads]
        }
	    // SPADES_OXFORD(
        //    nanopore_reads,
        //    ch_spades_hmm
        // )
        
        FLYE(
           nanopore_reads,
           "--nano-raw"
        )
    
    }


    


    
    

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
    ch_multiqc_files = ch_multiqc_files.mix(mergedtsv.collect().ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(CONVERT_CONFIDENCE.out.tsv.collect())
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))
    // // ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.stats.collect{it[1]}.ifEmpty([]))
    if (params.trim){
        ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.reads.collect{it[1]}.ifEmpty([]))
    }
    if (!params.skip_plots){
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MOVE_NANOPLOT.out.html.collect{it[1]}.ifEmpty([]))
    }
    
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
