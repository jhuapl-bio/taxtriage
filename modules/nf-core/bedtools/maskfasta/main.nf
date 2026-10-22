process BEDTOOLS_MASKFASTA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bed)
    path  fasta

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // bedtools maskfasta cannot read gzip/bgzip FASTA (arq5x/bedtools2#985): it aborts with
    // std::out_of_range. Stream it through `zcat -f`, which also passes plain FASTA through unchanged.
    """
    bedtools \\
        maskfasta \\
        $args \\
        -fi <(zcat -f $fasta) \\
        -bed $bed \\
        -fo ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
