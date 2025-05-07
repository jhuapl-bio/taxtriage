process GENERATE_SAMPLESHEET {
    tag "generate_temp_samplesheet"
    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.8.3' }"

    input:
    val(meta) // a map with all fields

    output:
    path 'temp_samplesheet.csv', emit: csv
    path "versions.yml", emit: versions

    script:
    def sampleName = meta.sampleName ?: 'sample'
    def platform   = meta.platform   ?: 'ILLUMINA'
    def fastq_1    = meta.fastq_1
    def fastq_2    = meta.fastq_2    ?: ''
    def seq_sum    = meta.seq_summary ?: ''
    def trim       = meta.trim       ?: 'false'
    def type       = meta.type       ?: 'UNKNOWN'

    """
    echo "sample,platform,fastq_1,fastq_2,sequencing_summary,trim,type" > temp_samplesheet.csv
    echo "${sampleName},${platform},${fastq_1},${fastq_2},${seq_sum},${trim},${type}" >> temp_samplesheet.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS

    """
}
