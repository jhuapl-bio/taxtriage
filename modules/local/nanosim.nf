process NANOSIM_SIMULATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanosim:3.2.3--hdfd78af_0' :
        'biocontainers/nanosim:3.2.3--hdfd78af_2' }"

    input:
    tuple val(meta), path(genome_list), path(size_file), path(nanosim_training)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.nanosim"

    """
    simulator.py metagenome \
        -gl ${genome_list} \
        -a ${size_file} \
        -c ${nanosim_training} \
        -o ${prefix} \
        --fastq \
        -t ${task.cpus} \
        ${args}

    # Compress output FASTQ files
    for f in ${prefix}*.fastq ${prefix}*.fq; do
        [ -f "\$f" ] && gzip "\$f"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanosim: \$(simulator.py --version 2>&1 | head -1 || echo "3.2.3")
    END_VERSIONS
    """
}
