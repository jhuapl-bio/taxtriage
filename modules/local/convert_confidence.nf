process CONVERT_CONFIDENCE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'pegi3s/biopython:latest' }"

    input:
    tuple val(meta),  path(kraken_report), path(tsv)

    output:
    path("*confidences.merged.tsv"), optional: false, emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when


    

    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    
    def output_parsed = "${meta.id}.confidences.merged.tsv"
    // println "${meta.id} ${kraken_report} ${tsv}"

    """

    mergeConfidence.py \\
        -i $tsv \\
        -o $output_parsed \\
        -s ${meta.id} \\
        -k $kraken_report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version )
    END_VERSIONS

    """
}

