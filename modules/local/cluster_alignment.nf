// ##############################################################################################
// # Copyright 2022-25 The Johns Hopkins University Applied Physics Laboratory LLC
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
process CLUSTER_ALIGNMENT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pandas" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(meta), path(bamfile)


    output:
    tuple val(meta), path("*fasta"), optional: false, emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline for maine coon

    outfile = "${meta.id}.cluster.fasta"

    """

        cluster_alignments.py \\
            --bam $bamfile \\
            --out $outfile \\
            --window 25 \\
            --top-unique-pct 5 \\
            --min-unique-per-ref 3 \\
            --min-mapq 10 \\
    """
}
