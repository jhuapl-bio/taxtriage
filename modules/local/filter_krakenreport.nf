process FILTERKRAKEN {
    label 'process_medium'

    /* conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'pegi3s/biopython:latest' }" */

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.8.3' }"


    input:
    path krakenreport

    output:
    path("*_filtered_mqc.tsv"), optional:true, emit: reports
    /* path("species_filtered.tsv"), optional:true, emit: speciesreport
    path("genus_filtered.tsv"), optional:true, emit: genusreport
    path("family_filtered.tsv"), optional:true, emit: familyreport
    path("order_filtered.tsv"), optional:true, emit: orderreport
    path("class_filtered.tsv"), optional:true, emit: classreport
    path("phylum_filtered.tsv"), optional:true, emit: phylumreport */

    script:

    """
    filter_tsv.py \\
        -o ./ \\
        -i $krakenreport

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
