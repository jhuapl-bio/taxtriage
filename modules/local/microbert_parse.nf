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
process MICROBERT_PARSE {
    label 'process_low'
    tag "${meta.id}"

    conda (params.enable_conda ? "bioconda::pandas" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(meta), path(reps), path(predictions), path(clusters), val(modelname)


    output:
    tuple val(meta), path("*microbert.annotations.tsv"), optional: false, emit: annotations
    tuple val(meta), path("*microbert.report.tsv"), optional: false, emit: report

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline for maine coon

    outfile = "${meta.id}.microbert.annotations.tsv"
    outreport = "${meta.id}.microbert.report.tsv"

    """
        map_clusters_to_taxa.py \\
            --fasta  ${reps} \\
            --json ${predictions} \\
            --taxa-report ${outreport} \\
            --clusters ${clusters} \\
            --out ${outfile} \\
            --modelname ${modelname} \\
            --include-prob
    """
}
