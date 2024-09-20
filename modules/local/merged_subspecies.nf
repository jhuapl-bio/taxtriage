process MERGEDSUBSPECIES {
    label 'process_medium'
    tag "MergeTopSpeciesAllSamples"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.75' }"


    input:
    file(reports)
    file(pathogens)

    output:
    path("*complete.krakenreport.merged.hierarchy.csv"), optional:false, emit: mergedhierarchy
    path("*.single.csv"), optional:true, emit: individualhierarchy


    script:

    pathogens_Var = " -p ${pathogens} "

    """

    merge_subspecies.py \\
        -o complete.krakenreport.merged.hierarchy.csv \\
        -i $reports \\
        -e S $pathogens_Var



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
