process FEATURES_DOWNLOAD {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://pegi3s/biopython:latest' :
        'biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(gcfs)
    path(assembly)

    output:
    tuple val(meta), path("**/*feature_table.txt"), emit: features, optional: false
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outputdir = "${meta.id}_features"
    def output = "${meta.id}_features.txt"

    """

    download_features.py \\
        -i ${gcfs} \\
        -a  ${assembly} \\
        -o  ${outputdir} \\
        -d

    for f in \$(find ${outputdir} -name "*feature_table.txt"); do
        cat \$f
    done > ${output}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version )
    END_VERSIONS
    """
}
