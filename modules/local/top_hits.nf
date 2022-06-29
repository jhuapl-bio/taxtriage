process TOP_HITS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(report)
    val(top_hits_count)

    output:
    path "versions.yml"           , emit: versions
    tuple val(meta), path("*_top_report.tsv"), optional:false, emit: tops



    when:
    task.ext.when == null || task.ext.when


    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    def id = "${meta.id}"
    """
    echo ${meta.id} "-----------------META variable------------------"
    get_top_hits.py \\
        -i $report \\
        -o $id \\
        -t $top_hits_count
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS

    """
}


