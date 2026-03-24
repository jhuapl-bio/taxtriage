process BEDTOOLS_GENOMECOVERAGE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bedtools=2.31.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/6397750e9730a3fbcc5b4c43f14bd141c64c723fd7dad80e47921a68a7c3cd21/data'
        : 'community.wave.seqera.io/library/bedtools_coreutils:a623c13f66d5262b'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bedgraph"), emit: bedgraph
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def reference = genome_file ? "-g ${genome_file} -sorted" : ""
    def bedgraph = "${meta.id}.bedgraph"
    def minmapq = params.minmapq ? "-q ${params.minmapq}" : ""

    """

    bedtools genomecov -ibam $bam -bg > $bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//' ))
    END_VERSIONS
    """
}
