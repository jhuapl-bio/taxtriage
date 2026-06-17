// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # (header trimmed for the sketch -- copy the full APL header from modules/local/confidence.nf)
// ##############################################################################################
//
// Compute the per-sample novelty score from the translated-search LCA output, the best-hit
// identity table, and the read accounting already produced upstream (K2 + alignment).
// Bundled script lives in bin/novelty_score.py (same pattern as sam_to_confidence.py).
//
process NOVELTY_SCORE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.75' }"

    input:
    // meta carries the read-accounting ints folded on by EXTRACT_UNMAPPED (see subworkflow)
    tuple val(meta), path(lca), path(tophit)
    path  run_summaries        // optional sibling-sample baseline; NO_FILE when first pass

    output:
    tuple val(meta), path("*.novelty.summary.tsv")   , emit: summary
    tuple val(meta), path("*.novelty.candidates.tsv"), emit: candidates
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def baseline  = run_summaries.name != 'NO_FILE' ? "--run-summaries ${run_summaries}" : ''
    def flag_z    = params.novelty_flag_z   ?: 2.0
    def weights   = params.novelty_weights  ?: '0.5,0.3,0.2'
    def idnt_cut  = params.novelty_idnt_cut ?: 50.0
    """
    novelty_score.py \\
        -s ${meta.id} \\
        --lca ${lca} \\
        --tophit ${tophit} \\
        --total-reads ${meta.total_reads ?: 0} \\
        --k2-classified ${meta.k2_classified ?: 0} \\
        --ref-aligned ${meta.ref_aligned ?: 0} \\
        --flag-threshold ${flag_z} \\
        --weights ${weights} \\
        --idnt-cut ${idnt_cut} \\
        ${baseline} ${args} \\
        -o ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.novelty.summary.tsv ${meta.id}.novelty.candidates.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
