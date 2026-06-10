process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_standard'

    conda (params.enable_conda ? 'bioconda::minimap2=2.28 bioconda::samtools=1.20' : null)
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37671219cfd244eb9b33db9345d3543ffd83037419a1c57f4648aace493ec2c2/data' :
        'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce' }"

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

    // Resolve the minimap2 preset.
    // Priority: meta.minimap2_preset (user-supplied per-sample) >
    //           meta.platform-derived default.
    // Valid preset names (passed without the leading -ax):
    //   map-ont, map-pb, map-pb (CLR), map-hifi (HiFi/CCS), lr:hq (ONT Q20),
    //   sr (short reads), splice (spliced long reads), splice:hq (Iso-seq),
    //   splice:sr (short-read RNA-seq), asm5/asm10/asm20 (assembly alignment),
    //   ava-pb, ava-ont (read overlap).
    def mapx = ''
    if (meta.minimap2_preset) {
        mapx = "-ax ${meta.minimap2_preset}"
    } else if (meta.platform =~ /(?i)illumina/) {
        mapx = '-ax sr'
    } else if (meta.platform =~ /(?i)pacbio/) {
        mapx = '-ax map-hifi'
    } else {
        // Default: Oxford Nanopore (also used for unspecified / FASTA inputs)
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

    // Derive -I and -K from 80% of the container's allocated RAM
    def ram_80pct_mb = (task.memory.toMega() * 0.80).longValue()
    def I_gb         = Math.max((ram_80pct_mb / 1024) as long, 1L)
    // Cap the derived query minibatch at 2G: a large -K is the main OOM driver
    // (it loads that many query bases into RAM on top of the index).
    def K_mb         = Math.min(Math.max((ram_80pct_mb / 8) as long, 64L), 2048L)
    def I_value      = params.mmap2_I ?: "${I_gb}G"
    def K_value      = params.mmap2_K ?: "${K_mb}M"

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