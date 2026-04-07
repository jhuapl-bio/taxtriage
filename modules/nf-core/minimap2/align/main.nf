process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_standard'

    conda (params.enable_conda ? 'bioconda::minimap2=2.21 bioconda::samtools=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' }"

    input:
    tuple val(meta), path(reads), path(reference)
    val bam_format
    val cigar_paf_format
    val cigar_bam

    output:
    tuple val(meta), path("*.paf"), optional: true, emit: paf
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def mapx = ''
    if (meta.platform =~ /(?i)illumina/) {
        mapx = '-ax sr'
    } else if (meta.platform =~ /(?i)pacbio/) {
        mapx = '-ax map-hifi'
    } else {
        mapx = '-ax map-ont'
    }

    def input_reads = reads.findAll { it != null }.join(' ')

    // Reserve most CPUs for minimap2; keep sort small
    def sort_threads = task.cpus > 4 ? 2 : 1
    def mm2_threads  = Math.max(task.cpus - sort_threads, 1)

    // Sort memory is PER THREAD, so set it explicitly and conservatively
    def sort_mem_total_mb      = (task.memory.toMega() * 0.20).longValue()
    def sort_mem_per_thread_mb = Math.max((sort_mem_total_mb / sort_threads) as long, 768L)
    sort_mem_per_thread_mb     = Math.min(sort_mem_per_thread_mb, 4096L)
    def S_value = "${sort_mem_per_thread_mb}M"

    // Treat minimap2 memory knobs as tuning params, not hard caps
    def I_value = params.mmap2_I ?: '8G'
    def K_value = params.mmap2_K ?: '100M'

    def cigar_paf     = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    def mmap2_window  = params.mmap2_window ? "-w ${params.mmap2_window}" : ''
    def mmap2_fraction_filter = params.mmap2_fraction_filter ? "-f ${params.mmap2_fraction_filter}" : ''
    def split_prefix  = params.no_split_prefix ? "" : "--split-prefix ${meta.id}.prefix"

    def bam_output = bam_format
        ? "-a | samtools sort -@ ${sort_threads} -m ${S_value} -T ${prefix}.tmp -O BAM -o ${prefix}.bam -"
        : "-o ${prefix}.paf"

    """
    minimap2 \\
        $args $mapx \\
        -t ${mm2_threads} -I ${I_value} -K ${K_value} ${split_prefix} \\
        ${reference} \\
        ${input_reads} \\
        ${cigar_paf} ${mmap2_window} ${mmap2_fraction_filter} \\
        ${set_cigar_bam} \\
        ${bam_output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}