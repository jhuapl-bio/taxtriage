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
process CREATE_COMPARISON_REPORT {
    tag "comparison_full_report"
    label 'process_medium'
    publishDir "${params.outdir}/report", mode: 'copy'

    conda (params.enable_conda ? "bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/reportlab-pdf:4.0.9' :
        'jhuaplbio/reportlab-pdf:4.0.9' }"

    input:
    // One or more .paths.json files from ALIGNMENT_PER_SAMPLE (collected across all samples)
    path(json_files)
    // HTML template (heatmap.html)
    path(template)
    // Optional protein-annotation XLSX files from ORGANISM_MERGE_REPORT --output_annot_xlsx
    // Pass a NO_FILE placeholder when protein annotations are not available
    path(protein_annotations)

    output:
        path "versions.yml"           , emit: versions
        path("*html")                 , optional: true, emit: html

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    def output_html = "all.comparison.report.html"

    // Build the list of JSON input files (filter out any NO_FILE placeholders)
    def json_inputs = json_files instanceof List
        ? json_files.findAll { it.name != 'NO_FILE' && it.name.endsWith('.json') }.join(' ')
        : (json_files.name != 'NO_FILE' && json_files.name.endsWith('.json') ? json_files.toString() : '')

    // Build optional protein annotations argument
    def prot_arg = ''
    if (protein_annotations) {
        def prot_files = protein_annotations instanceof List
            ? protein_annotations.findAll { it.name != 'NO_FILE' && !it.name.startsWith('NO_FILE') }.join(' ')
            : (protein_annotations.name != 'NO_FILE' && !protein_annotations.name.startsWith('NO_FILE') ? protein_annotations.toString() : '')
        if (prot_files) {
            prot_arg = "-p ${prot_files}"
        }
    }

    """
    make_report.py -i ${json_inputs} \\
        -t ${template} \\
        -o ${output_html} \\
        ${prot_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS

    """
}


