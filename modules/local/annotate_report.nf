process ANNOTATE_REPORT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3 conda-forge::pandas conda-forge::openpyxl" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(meta), path(diamond_txt)
    path(annotate_meta)

    output:
    tuple val(meta), path("*.xlsx"), emit: xlsx
    tuple val(meta), path("*.tsv"),  emit: tsv
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    annotate_report.py \\
        --diamond ${diamond_txt} \\
        --meta ${annotate_meta} \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        openpyxl: \$(python3 -c "import openpyxl; print(openpyxl.__version__)")
    END_VERSIONS
    """
}
