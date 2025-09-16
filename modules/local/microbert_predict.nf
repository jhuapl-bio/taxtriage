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
process MICROBERT_PREDICT {
    label 'process_high'
    tag "${meta.id}"

    conda (params.enable_conda ? "R" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jhuaplbio/microbert-classify:1.0.0' :
        'microbert-classify:1.0.0' }"          // Fallback Docker image

    input:
    tuple val(meta), path(input), path(model), val(modelname)


    output:
    tuple val(meta), path("*json"), optional: false, emit: predictions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline for maine coon

    outfile = "${meta.id}_microbert_predictions.json"

    """
        ln -s $model /analysis/data

        python3 /analysis/analysis/experiment/test_sequences.py -i ${input} \\
            -o ${outfile} \\
            -b 50 \\
            -m ${modelname}
    """
}
