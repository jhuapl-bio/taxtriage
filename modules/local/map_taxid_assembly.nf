process MAP_TAXID_ASSEMBLY {
    label 'process_medium'
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://bmerritt1762/jhuaplbio/reportlab-pdf:4.0.7' :
        'jhuaplbio/reportlab-pdf:4.0.7' }"


    input:
    tuple val(meta), file(gcfmapping)
    file(assembly)

    output:
    tuple val(meta), path("*merged.taxid.tsv"), optional:false, emit: taxidmerged



    script:


    """

   append_taxid.py \\
        -i $gcfmapping \\
        -r $assembly \\
        -o ${meta.id}.merged.taxid.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
