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
process TOP_HITS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(report) 

    output:
    path "versions.yml"           , emit: versions
    tuple val(meta), path("*top_report.tsv"), optional:false, emit: tops
    tuple val(meta), path("*.krakenreport_mqc.tsv"), optional:false, emit: krakenreport



    when:
    task.ext.when == null || task.ext.when


    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    def id = "${meta.id}"
    ch_top_per_taxa = ""
    def top_per_taxa  = params.top_per_taxa ? " -s ${params.top_per_taxa} " : ''
    def top_hits_count = params.top_hits_count ? " -t ${params.top_hits_count}" : ' -t 10 '
    """
    echo ${meta.id} "-----------------META variable------------------"
    get_top_hits.py \\
        -i \"$report\" \\
        -o ${id}.top_report.tsv   \\
        $top_hits_count  $top_per_taxa
    
    awk -F '\\t' -v id=${id} \\
        'BEGIN{OFS=\"\\t\"} { if (NR==1){ print \"Sample_Taxid\",  \$1, \$4, \$6} else { \$5 = id\"_\"\$5;  print \$5, \$1, \$4, \$6  }}'  ${id}.top_report.tsv > ${id}.krakenreport_mqc.tsv 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS

    """
}

