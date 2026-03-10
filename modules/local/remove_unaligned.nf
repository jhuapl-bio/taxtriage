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
    // if params.use_denovo or params.us_diamond is true and there are paired-end reads, then force ignore include_singletons_hostremoval
    def flag = !params.include_singletons_hostremoval ? "-f 12 -F 0x900" : "-f 4 -F 0x900"
    if ((params.use_denovo || params.use_diamond) && paired && params.include_singletons_hostremoval) {
        println "ALERT: paired unmapped only (include_singletons_hostremoval parameter) has been ignored because --use_denovo or --use_diamond was specified with paired-end reads. Both reads in a pair must be unmapped to be retained as Megahit fails without equal read counts."
        flag = "-f 12 -F 0x900"
    }

    def singleton_out = params.include_singletons_removal ? "${meta.id}.hostremoved.singletons.fastq.gz" : "/dev/null"
    def cmd = !meta.single_end ? \
        " -1 ${meta.id}_1.hostremoved.fastq.gz -2 ${meta.id}_2.hostremoved.fastq.gz -0 ${singleton_out}" : \
        meta.platform == "OXFORD" ? "-0 ${meta.id}.hostremoved.fastq.gz -s /dev/null" : \
        " -o ${meta.id}.hostremoved.fastq.gz -s /dev/null -0 /dev/null"
    """
    samtools view -b ${flag} ${input} | \\
        samtools fastq -n ${args} \\
        ${cmd} -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
