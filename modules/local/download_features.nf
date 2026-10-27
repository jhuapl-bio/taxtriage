process FEATURES_DOWNLOAD {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython%3A1.75' :
        'biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(gcfs)
    path(assembly)
    val(get_aas)
    val(product_types)   // --mask_products: feature types for <id>.products.bed ('' = no BED)

    output:
    tuple val(meta), path("${meta.id}_features.txt"), emit: features, optional: false
    tuple val(meta), path("${meta.id}_proteins.faa"), emit: proteins, optional: true
    tuple val(meta), path("${meta.id}.map.txt"), emit: mapfile, optional: true
    tuple val(meta), path("${meta.id}.products.bed"), emit: products, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outputdir = "${meta.id}_features"
    def output = "${meta.id}_features.txt"
    def aa_output = "${meta.id}_proteins.faa"
    def mapfile = "${meta.id}.map.txt"
    def geta = get_aas ? "-p" : ""
    def bed  = product_types ? "-b ${meta.id}.products.bed -t '${product_types}'" : ""

    """
    download_features.py \\
        -i ${gcfs} \\
        -a ${assembly} \\
        -o ${outputdir} \\
        -m ${mapfile} \\
        -d ${geta} ${bed}

    for f in \$(find ${outputdir} -name "*feature_table.txt"); do
        cat \$f
    done > ${output}

    for f in \$(find ${outputdir} -name "*protein.faa"); do
        cat \$f
    done > ${aa_output}

    # remove all faa files that is not full.faa
    find ${outputdir} -name "*.faa" -not -name "${aa_output}" -exec rm {} \\;
    find ${outputdir} -name "*_feature_table.txt" -not -name "${output}" -exec rm {} \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version )
    END_VERSIONS
    """
}
