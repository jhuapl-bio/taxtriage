// Normalise the user's spike-in sheet (csv/tsv/xlsx) into the TSV the rest of the
// spike-in path consumes, plus the distinct accession list to fetch references for.
process PARSE_SPIKEIN {
    tag "spikein_sheet"
    label 'process_single'
    publishDir "${params.outdir}/spikein", mode: 'copy'

    conda (params.enable_conda ? "conda-forge::openpyxl" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    path(sheet)
    val(default_replicates)

    output:
    path("spikein.tsv")        , emit: tsv
    path("spikein_accessions.txt"), emit: accessions
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reps = default_replicates ?: 1
    """
    parse_spikein_sheet.py \\
        --input ${sheet} \\
        --output spikein.tsv \\
        --accessions spikein_accessions.txt \\
        --default-replicates ${reps}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
