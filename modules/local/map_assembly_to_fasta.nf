process MAP_LOCAL_ASSEMBLY_TO_FASTA {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://pegi3s/biopython:latest' :
        'biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(fasta)
    path(assembly)
    path(pathogens_file)

    output:
    tuple val(meta), path("*localmap.tsv"), emit: map, optional: false
    tuple val(meta), path("*output.gcfids.txt"), emit: accessions, optional: false
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${meta.id}.localmap.tsv"
    def pathogens  = pathogens_file  ? "-p ${pathogens_file}" : ""

    """

    fuzzy_match_assembly.py  \\
        -i ${fasta} \\
        -a ${assembly} \\
        -o ${output} $pathogens

    cut -f 2 ${output} > ${meta.id}.output.gcfids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version )
    END_VERSIONS
    """
}
