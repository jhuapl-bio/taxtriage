process UNMAPPED_MINIMAP2 {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::minimap2=2.21 bioconda::samtools=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' }"

    input:
    tuple val(meta), path(reads), path(reference)

    output:
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // define mapx if paired end reads or single end illumina or single end oxford, match with regex ignore case
    def mapx = ''
    if (meta.platform =~ /(?i)illumina/) {
        mapx = '-ax sr'
    } else if (meta.platform =~ /(?i)pacbio/) {
        mapx = '-ax map-hifi'
    } else {
        mapx = '-ax map-ont'
    }
    def input_reads = reads.findAll { it != null }.join(' ')
    // def minmapq = minmapq ? " -q ${minmapq} " :  ""
    def I_value = "${(task.memory.toMega() * Math.min(0.8 / task.cpus, 0.8)).longValue()}M"
    def S_value = "${(task.memory.toMega() * Math.min(0.15 / task.cpus, 0.15)).longValue()}M"
    // if input is illumina then use -ax sr else use -ax map-ont

    """

    minimap2 \\
        $args $mapx \\
        -a -t $task.cpus -I $I_value \\
        $reference \\
        $input_reads | \\
        samtools sort -@ ${task.cpus} -m $S_value | \\
        samtools view -@ ${task.cpus} -b -h -o ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
