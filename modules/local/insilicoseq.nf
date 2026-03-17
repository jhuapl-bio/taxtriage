process INSILICOSEQ_SIMULATE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/insilicoseq:2.0.1--pyh7cba7a3_0' :
        'biocontainers/insilicoseq:2.0.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(reference_fasta), path(abundance_file)
    val(num_reads)
    val(model)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.iss"
    def iss_model = model ?: 'miseq'
    def mode = task.ext.mode ?: 'kde'

    """
    iss generate \
        --genomes ${reference_fasta} \
        --model ${iss_model} \
        --output ${prefix} \
        --mode ${mode} \
        --abundance_file ${abundance_file} \
        -n ${num_reads} \
        --cpus ${task.cpus} \
        ${args}

    # Compress FASTQ files if not already gzipped
    for f in ${prefix}_R*.fastq; do
        [ -f "\$f" ] && gzip "\$f"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        insilicoseq: \$(iss --version 2>&1 | sed 's/^.*iss //')
    END_VERSIONS
    """
}
