process PULL_TAXID {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:4.2.0' :
        'cicirello/gnu-on-alpine' }"

    input:
    tuple val(meta), val(assembly_hits)
    val row


    output:
    tuple val(meta), path("*.fasta"), optional: false, emit: genome
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when


    

    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    def output = "."
    println assembly_hits
    """
    mkdir -p $output
    bash refseq_download_single.sh -i $assembly_hits -o $output
    gzip -d $output/*.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(gawk --version 2>&1)
    END_VERSIONS

    """
}
