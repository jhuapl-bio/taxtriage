process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(input_files)
    // path fasta
    // path fai

    output:
    tuple val(meta), path("${prefix}.bam") , optional:true, emit: bam
    tuple val(meta), path("${prefix}.cram"), optional:true, emit: cram
    path  "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def file_type = input_files[0].getExtension()
    // def reference = fasta ? "--reference ${fasta}" : ""
    def cpus = task.cpus > 1 ? task.cpus - 1 : 1
    def S_value = "${(task.memory.toMega() * Math.min(0.15 / task.cpus, 0.15)).longValue()}M"

    """
    samtools \\
        merge -f \\
        $args \\
        --threads ${cpus} \\
        -u - \\
        $input_files | samtools sort -@ ${cpus} -m $S_value -o ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
