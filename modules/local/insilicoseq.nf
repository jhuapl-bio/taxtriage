process INSILICOSEQ_SIMULATE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/insilicoseq:1.3.5--py_0' :
        'biocontainers/insilicoseq:1.3.5--py_0' }"

    input:
    tuple val(meta), path(genome_list), path(abundance_file)
    val(reference_fasta)
    val(num_reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: simulated
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.simulated"
    def mode = "kde" // default to KDE mode, but allow override with task.ext.mode
    if (task.ext.mode) {
        if (["kde", "empirical"].contains(task.ext.mode)) {
            mode = task.ext.mode
        } else {
            log.warn "Invalid mode specified in task.ext.mode: ${task.ext.mode}. Defaulting to 'kde'."
        }
    }



    """

    iss generate --genomes $reference_fasta \
        --model $iss_model \
        --output $prefix \
        --mode $mode \
        --abundance_file $abundance_file \
        -n $num_reads \
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        insilicoseq: \$(echo \$(iss --version 2>&1) | sed 's/^.*iss v//' ))
    END_VERSIONS
    """
}
