process NANOSIM_SIMULATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanosim:3.2.3--hdfd78af_0' :
        'biocontainers/nanosim:3.2.3--hdfd78af_2' }"

    input:
    tuple val(meta), path(genome_list), path(abundance_file), path(nanosim_training)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: simulated
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    """

    simulator.py metagenome \
          -gl "$genome_list" \
          -a "$$abundance_file" \
          -c "$nanosim_training" \
          -o "ont_simreads" \
          --fastq \


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanosim: \$(echo \$(nanosim --version 2>&1) | sed 's/^.*nanosim v//' ))
    END_VERSIONS
    """
}
