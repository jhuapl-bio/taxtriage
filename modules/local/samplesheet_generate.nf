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
    // Deliberately NOT defaulted to ILLUMINA here.  When fastq_1 is an SRA/ENA
    // accession, input_check.nf's merge_sra_row() treats a non-empty platform
    // column as "the user declared this" and skips the instrument_platform that
    // ENA/NCBI reported -- so defaulting here would silently process an ONT or
    // PacBio accession as Illumina.  The ILLUMINA fallback lives in
    // create_fastq_channel(), which still applies to local-file rows.
    def platform   = meta.platform ?: ''
    def fastq_1    = meta.fastq_1    ?: ''
    def fastq_2    = meta.fastq_2    ?: ''
    def bam        = meta.bam        ?: ''
    def seq_sum    = meta.seq_summary ?: ''
    def trim       = meta.trim       ?: 'false'
    def type       = meta.type       ?: 'UNKNOWN'

    """
    echo "sample,platform,fastq_1,fastq_2,bam,sequencing_summary,trim,type" > temp_samplesheet.csv
    echo "${sampleName},${platform},${fastq_1},${fastq_2},${bam},${seq_sum},${trim},${type}" >> temp_samplesheet.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS

    """
}
