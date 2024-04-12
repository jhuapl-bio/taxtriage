process COMBINE_MAPFILES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:4.2.0' :
        'biocontainers/gawk:4.2.0' }"

    input:
    tuple val(meta), path(gcfmaps), path(gcfids)

    output:
    tuple val(meta), path("*.combined.gcfmap.tsv"), path("*.combined.gcfids.txt"), optional: false, emit: mergefiles
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """


    cat ${gcfids} | sort | uniq > ${meta.id}.combined.gcfids.txt
    cat ${gcfmaps} | sort | uniq  > ${meta.id}.combined.gcfmap.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version 2>&1)
    END_VERSIONS
    """
}
