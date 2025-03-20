process KHMER_NORMALIZEBYMEDIAN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/khmer:3.0.0a3--py37haa7609a_2' :
        'biocontainers/khmer:3.0.0a3--py37haa7609a_2' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: reads
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def m = task.memory * 0.8
    def mode = meta.single_end ? ' ' : ' --force_single '


    """
    normalize-by-median.py \\
        -M ${m.toMega()}e6 \\
        --gzip $args \\
        -o ${meta.id}.fastq.gz \\
        $mode \\
        $files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer: \$( normalize-by-median.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """
}
