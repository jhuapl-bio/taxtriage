//
// Propagate per-representative taxonomy (LCA + best-hit) back to every linclust cluster member,
// so NOVELTY_SCORE sees member-level abundance even though the expensive search ran on reps only.
// Bundled script: bin/expand_clusters.py.
//
process EXPAND_CLUSTERS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(lca), path(tophit), path(clusters)

    output:
    tuple val(meta), path("*.lca.expanded.tsv")   , emit: lca
    tuple val(meta), path("*.tophit.expanded.tsv"), emit: tophit
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.novelty"
    """
    expand_clusters.py \\
        --clusters ${clusters} \\
        --lca ${lca} \\
        --tophit ${tophit} \\
        --out-lca ${prefix}.lca.expanded.tsv \\
        --out-tophit ${prefix}.tophit.expanded.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.novelty"
    """
    touch ${prefix}.lca.expanded.tsv ${prefix}.tophit.expanded.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
