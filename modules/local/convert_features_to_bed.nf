process FEATURES_TO_BED {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:4.2.0' :
        'biocontainers/gawk:4.2.0' }"

    input:
    tuple val(meta), path(features)

    output:
    tuple val(meta),  path("*.features.bed"), optional: false, emit: bed
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    """



    awk -F '\t' 'NR > 1 { if (\$1 == "CDS" && \$8 < \$9 ){print \$7\"\t\"\$8\"\t\"\$9\"\t\"\$14 }}' \\
        ${features} > ${meta.id}.features.bed


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version )
    END_VERSIONS
    """
}
