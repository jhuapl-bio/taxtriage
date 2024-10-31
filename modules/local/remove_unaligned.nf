process REMOVE_HOSTREADS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.hostremoved.{fastq,fq}.gz"), optional: true, emit: reads
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def paired = meta.paired ? "-1 ${input[0]} -2 ${input[1]}" : "-U ${input}"
    def flag = !meta.single_end ? "-f 12" : "-f 4"
    def cmd = !meta.single_end ? \
        " -1 ${meta.id}_1.hostremoved.fastq.gz -2 ${meta.id}_2.hostremoved.fastq.gz" : \
        meta.platform == "OXFORD" ? "-0 ${meta.id}.hostremoved.fastq.gz -s /dev/null" : \
        " -o ${meta.id}.hostremoved.fastq.gz -s /dev/null -0 /dev/null"
    """
    samtools view -b ${flag} ${input} | \\
        samtools fastq \\
        ${cmd} - \\
        ${args} -n

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
