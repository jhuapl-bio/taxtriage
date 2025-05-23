process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(bam)
    val(minmapq)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    tuple val(meta), path("*.csi"), emit: csi, optional: true
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def minmapq = minmapq ? " -q ${minmapq} " :  ""

    def I_value = "${(task.memory.toMega() * Math.min(0.85 / task.cpus, 0.85)).longValue()}M"



    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    samtools sort \\
        $args \\
        -@ $task.cpus \\
        -m $I_value \\
        $bam | \\
        samtools view $minmapq -@ ${task.cpus} -b -h -o ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
