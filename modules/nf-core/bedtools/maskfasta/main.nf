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
    // std::out_of_range. Stream it through a pipe, decompressing only when the file is gzip
    // (magic bytes) -- portable, unlike relying on `zcat -f` to pass plain FASTA through.
    """
    bedtools \\
        maskfasta \\
        $args \\
        -fi <(if [ "\$(head -c 2 $fasta | od -An -tx1 | tr -d ' \\n')" = "1f8b" ]; then gzip -dc $fasta; else cat $fasta; fi) \\
        -bed $bed \\
        -fo ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
