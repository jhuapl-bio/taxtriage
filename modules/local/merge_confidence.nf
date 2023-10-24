process MERGE_CONFIDENCE {
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.75' }"

        
    input:
    file confidences

    output:
    path("*confidences*.merged_mqc.tsv"), optional:false, emit: confidence_report

    script:
    """
    merge_tsvs.py \\
        -o ./confidences.merged_mqc.tsv \\
        -i $confidences 

    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}