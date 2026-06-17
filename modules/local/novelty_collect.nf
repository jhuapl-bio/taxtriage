// ##############################################################################################
// # Copyright 2022-25 The Johns Hopkins University Applied Physics Laboratory LLC
// # (header trimmed for the sketch -- copy the full APL header from the other modules)
// ##############################################################################################
//
// Collect the per-sample NOVELTY_SCORE outputs into downloadable JSON/XLSX artifacts and the
// combined all.novelty.json feed that make_report.py bakes into the interactive HTML report.
//
process NOVELTY_COLLECT {
    tag "novelty_collect"
    label 'process_low'
    publishDir "${params.outdir}/report", mode: 'copy'

    conda (params.enable_conda ? "conda-forge::pandas conda-forge::openpyxl" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/reportlab-pdf:4.0.9' :
        'jhuaplbio/reportlab-pdf:4.0.9' }"

    input:
    // Collected across all samples (each may be a NO_FILE placeholder list when novelty is off)
    path(summaries)
    path(candidates)

    output:
    path "all.novelty.json"          , emit: combined_json
    path "*.novelty.json"            , emit: json_files     // per-sample + combined
    path "*.novelty.xlsx"            , optional: true, emit: xlsx_files
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sum_args  = (summaries instanceof List ? summaries : [summaries])
        .findAll { it.name != 'NO_FILE' && it.name.endsWith('.tsv') }.join(' ')
    def cand_args = (candidates instanceof List ? candidates : [candidates])
        .findAll { it.name != 'NO_FILE' && it.name.endsWith('.tsv') }.join(' ')
    """
    novelty_collect.py \\
        --summaries ${sum_args} \\
        --candidates ${cand_args} \\
        --outdir . \\
        --combined-prefix all.novelty

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"taxtriage_novelty": true, "version": "1.0", "samples": {}}' > all.novelty.json
    touch all.novelty.xlsx
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
