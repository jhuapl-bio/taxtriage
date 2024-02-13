process REFERENCE_REHEADER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:4.2.0' :
        'biocontainers/gawk:4.2.0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta),  path("*.reformat.fasta"), optional: false, emit: fasta
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}.reformat.fasta"

    """

    reheader_fasta.sh -i ${fasta} -o ${output} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$( awk --version )
    END_VERSIONS
    """
}
