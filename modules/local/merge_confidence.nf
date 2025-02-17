process MERGE_ALIGNMENT_MERGES {
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.75' }"


    input:
    file metrics

    output:
    path("*metrics*.merged_mqc.tsv"), optional:false, emit: metrics_report

    script:
    """
    merge_tsvs.py \\
        -o ./metrics.merged_mqc.tsv \\
        -i $metrics



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
