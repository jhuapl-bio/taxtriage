process MAP_TAXID_ASSEMBLY {
    label 'process_medium'
    tag "$meta.id"

    // The NCBI backup lookup hits E-utilities, which rate limits per IP
    // (3 req/s without an API key, 10 with). Serialise the tasks so parallel
    // samples don't collectively blow past the limit and get HTTP 429s.
    maxForks (params.ncbi_api_key ? 3 : 1)

    // Rate limiting is transient - give the task a couple of chances.
    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2

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
    def api_key        = params.ncbi_api_key ? " --api-key ${params.ncbi_api_key}" : ""
    def custom_map_arg = custom_map.name != "NO_FILE" ? " --custom-map ${custom_map}" : ""

    """

   append_taxid.py \\
        --ncbi-batch-size ${params.ncbi_batch_size} \\
        -i $gcfmapping \\
        -r $assembly \\
        --ncbi-backup${email}${api_key}${custom_map_arg} \\
        -o ${meta.id}.merged.taxid.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
