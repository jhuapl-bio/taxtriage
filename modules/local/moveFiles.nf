process MOVE_FILES {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:4.2.0' :
        'cicirello/gnu-on-alpine' }"

    input:
    tuple val(meta), path(files)
    val(pattern)
    val(append_id)
    val(replaces)

    output:
    tuple val(meta), path("*")     , optional:true, emit: files



    when:
    task.ext.when == null || task.ext.when


    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    

    """
    mv ${files} "${pattern}${files}"
    

    """
}


