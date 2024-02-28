process BEDTOOLS_DEPTHCOVERAGE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bed), path(bam)

    output:
    tuple val(meta), path("*.features.coverage"), emit: coverage
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def reference = genome_file ? "-g ${genome_file} -sorted" : ""


    """

    mkdir -p split_bed

    awk '{print > "split_bed/"\$1".bed"}'  $bed


    for bed_file in split_bed/*.bed; do

        output_file="\${bed_file%.bed}_coverage.txt"

        bedtools coverage $args  \\
            -a "\$bed_file" \\
            -b $bam -hist;

    done >  ${meta.id}.features.coverage

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//' ))
    END_VERSIONS
    """
}
