process KRONA { //https://github.com/nf-core/mag/blob/66cf53aff834d2a254b78b94fc54cd656b8b7b57/modules/local/krona.nf
    tag "${meta.classifier}-${meta.id}"

    conda "bioconda::krona=2.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krona:2.7.1--pl526_5' :
        'biocontainers/krona:2.7.1--pl526_5' }"

    input:
    tuple val(meta), path(report)

    output:
    path "*.html"       , emit: html
    path "versions.yml" , emit: versions

    script:
    """
    ktImportTaxonomy "$report" -tax taxonomy

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ktImportTaxonomy: \$(ktImportTaxonomy 2>&1 | sed -n '/KronaTools /p' | sed 's/^.*KronaTools //; s/ - ktImportTaxonomy.*//')
    END_VERSIONS
    """
}
