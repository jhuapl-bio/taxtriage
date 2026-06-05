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
    path(custom_map, stageAs: "custom_accession_map.tsv")

    output:
    tuple val(meta), path("*merged.taxid.tsv"), optional:false, emit: taxidmerged



    script:
    def email          = params.email ? " -e ${params.email}" : ""
    def custom_map_arg = custom_map.name != "NO_FILE" ? " --custom-map ${custom_map}" : ""

    """

   append_taxid.py \\
        -i $gcfmapping \\
        -r $assembly \\
        --ncbi-backup${email}${custom_map_arg} \\
        -o ${meta.id}.merged.taxid.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
