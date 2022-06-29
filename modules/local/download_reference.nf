process REFERENCE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:4.2.0' :
        'cicirello/gnu-on-alpine' }"

    


    input:
    tuple val(meta), val(taxid)
    path(assembly)
    


    output:
    path("*.ftmp"), optional: false, emit: assembly_hits
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when


     

    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    
    // def files = tops.splitCsv(header: true, sep="\t")
    def id = "${meta.id}"

    """

    bash taxid_to_reflist.sh -i $taxid -o $id".ftmp" -a $assembly
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version 2>&1)
    END_VERSIONS

    """
}
