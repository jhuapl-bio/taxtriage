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
process PATHOGENS_MERGE_REPORT {
    tag "Pathogens_Report"
    label 'process_medium'
    publishDir "${params.outdir}/report", mode: 'copy'

    conda (params.enable_conda ? "bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://bmerritt1762/jhuaplbio/reportlab-pdf:4.0.7' :
        'jhuaplbio/reportlab-pdf:4.0.7' }"

    input:
    tuple file(files_of_pathogens)

    output:
        path "versions.yml"           , emit: versions
        path("pathogens.report.txt")    , optional: false, emit: report
        path("pathogens.report.pdf")    , optional: false, emit: pdf


    when:
    task.ext.when == null || task.ext.when




    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/

    def output_txt = "pathogens.report.txt"
    def output_pdf = "pathogens.report.pdf"
    def distributions = params.distributions ? " -d ${file(params.distributions)} " : " "

    """

    awk 'NR==1{print; next} FNR>1' $files_of_pathogens > $output_txt

    create_report.py -i $output_txt  -o $output_pdf $distributions

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1)
    END_VERSIONS

    """
}


