process NCBIGENOMEDOWNLOAD_FEATURES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-genome-download:0.3.3--pyh7cba7a3_0' :
        'biocontainers/ncbi-genome-download:0.3.3--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(gcf_file)

    output:
    tuple val(meta),  path("*.features.txt"), optional: false, emit: features
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    ncbi-genome-download \\
        --assembly-accessions ${gcf_file} \\
        --output-folder ./ \\
        --flat-output -F features \\
        $args all


    find . -type f -name "*_feature_table.txt.gz" -exec gzip -d {}    \\;


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbigenomedownload: \$( ncbi-genome-download --version )
    END_VERSIONS
    """
}
