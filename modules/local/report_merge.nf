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
process ORGANISM_MERGE_REPORT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/report", mode: 'copy'

    conda (params.enable_conda ? "bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/reportlab-pdf:4.0.8' :
        'jhuaplbio/reportlab-pdf:4.0.8' }"

    input:
    tuple val(meta), file(files_of_pathogens), file(distributions)
    val(missing_samples)
    path(nodes)

    output:
        path "versions.yml"           , emit: versions
        path("*organisms.report.txt")    , optional: true, emit: report
        path("*organisms.report.pdf")    , optional: false, emit: pdf

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/

    def output_txt = "${meta.id}.organisms.report.txt"
    def output_pdf = "${meta.id}.organisms.report.pdf"
    def distribution_arg = distributions.name != "NO_FILE" ? " -d $distributions " : ""
    distribution_arg = ""
    def min_conf = params.min_conf || params.min_conf == 0 ? " -c $params.min_conf " : ""

    def missing_arg = ''
    if (missing_samples) {
        missing_arg = "-m \"${missing_samples.join(' ')}\""
    }
    def show_potentials = params.show_potentials ? " --show_potentials " : ""
    def show_commensals = params.show_commensals ? " --show_commensals " : ""
    def show_unidentified = params.show_unidentified ? " --show_unidentified " : ""
    def taxdump = nodes.name != "NO_FILE" ? " --taxdump $nodes " : ""
    """

    create_report.py -i $files_of_pathogens -u $output_txt  \\
        -o $output_pdf \\
        $show_potentials \\
        $show_commensals \\
        $show_unidentified \\
        $distribution_arg \\
        $min_conf $taxdump \\
        $missing_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1)
    END_VERSIONS

    """
}


