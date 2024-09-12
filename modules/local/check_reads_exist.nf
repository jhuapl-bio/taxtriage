process CHECK_GZIPPED_READS {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(fastq_files)
    val(minimum_reads_check)

    output:
    tuple val(meta), path("{emptyfile,minimum_reads_check}.txt"), optional: true, emit: check_result


    script:
    // Create a Bash command string that checks each file
    def checks = fastq_files.collect { fq ->
        // Use gzip -dc to decompress and wc -l to count lines
        return "if [ \$(gzip -dc ${fq} | wc -l) -ge ${minimum_reads_check} ]; then echo '${fq} contains fewer than ${minimum_reads_check} reads'; fi"
    }.join('; ')

    """
    echo 'Checking if FASTQ.gz files contain more than ${ minimum_reads_check } reads...' > minimum_reads_check.tmp
    $checks > minimum_reads_check.tmp
    if [ -s minimum_reads_check.tmp ]; then
        mv minimum_reads_check.tmp minimum_reads_check.txt
    else
        rm minimum_reads_check.tmp
        touch emptyfile.txt
    fi
    """
}
