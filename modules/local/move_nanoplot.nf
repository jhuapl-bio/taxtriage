process MOVE_NANOPLOT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:4.2.0' :
        'cicirello/gnu-on-alpine' }"

    input:
    tuple val(meta), path(html)

    output:
    tuple val(meta), path("*_mqc.html"), optional: false, emit: html
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when


    

    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    
    def output = "${meta.id}_mqc.html"

    """

    array=( $html )
    for i in "\${array[@]}"
    do
        mv \$i ${meta.id}_\${i}_mqc.html
    done


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version )
    END_VERSIONS

    """
}
