//
// Linear-time protein clustering (mmseqs easy-linclust) used as a SPEED dedup in front of the
// novelty translated-search. linclust is near-linear (vs the cascaded easy-cluster), so it
// collapses near-identical ORFs cheaply; only the cluster representatives go to MMSEQS_TAXONOMY.
// The cluster TSV (representative<TAB>member) is emitted so EXPAND_CLUSTERS can propagate each
// representative's LCA call back to every member -- nothing is dropped from the downstream score.
//
process MMSEQS_LINCLUST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:18.8cc5c--hd6d6fdc_0':
        'staphb/mmseqs2' }"

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("*_rep_seq.fasta"), emit: representatives
    tuple val(meta), path("*_cluster.tsv")  , emit: clusters
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '--min-seq-id 0.9 -c 0.8 --cov-mode 1'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mmseqs \\
        easy-linclust \\
        ${sequence} \\
        ${prefix} \\
        tmp_lc \\
        $args \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs version 2>/dev/null || mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rep_seq.fasta ${prefix}_cluster.tsv
    echo '"${task.process}": {mmseqs: stub}' > versions.yml
    """
}
