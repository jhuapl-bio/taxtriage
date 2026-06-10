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

    // minimap2 and samtools now run SEQUENTIALLY (minimap2 writes an intermediate
    // file, then samtools sorts it), so each stage gets the FULL cpu count and up to
    // 70% of container RAM -- they no longer run concurrently and compete for memory.
    def threads = task.cpus

    // -I / -K are NOT derived automatically anymore: when unset, minimap2 picks its own
    // defaults and manages memory itself. Only pass them when explicitly overridden via
    // params.mmap2_I / params.mmap2_K.
    def I_value = params.mmap2_I ? "-I ${params.mmap2_I}" : ""
    def K_value = params.mmap2_K ? "-K ${params.mmap2_K}" : ""

    // samtools sort -m is PER THREAD; total (-m * threads) stays within 70% of task RAM.
    def ram_70pct_mb = (task.memory.toMega() * 0.70).longValue()
    def sort_mem_per_thread_mb = Math.max((ram_70pct_mb / threads) as long, 768L)
    def S_value = "${sort_mem_per_thread_mb}M"

    def cigar_paf     = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    def mmap2_window  = params.mmap2_window ? "-w ${params.mmap2_window}" : ''
    def mmap2_fraction_filter = params.mmap2_fraction_filter ? "-f ${params.mmap2_fraction_filter}" : ''
    // Off by default; enable --split-prefix only for references too large to index in RAM as one block
    def split_prefix  = params.split_prefix ? "--split-prefix ${meta.id}.prefix" : ""

    // Write minimap2 output to an intermediate file instead of piping into samtools.
    // Decoupling the stages avoids broken-pipe / "truncated file" failures and lets
    // sort use the full CPU/RAM budget. The intermediate SAM is removed afterwards.
    def mm2_out   = bam_format ? "-a -o ${prefix}.unsorted.sam" : "-o ${prefix}.paf"
    def sort_step = bam_format ? """
    samtools sort -@ ${threads} -m ${S_value} -T ${prefix}.tmp -O BAM -o ${prefix}.bam ${prefix}.unsorted.sam""" : ""

    """
    # Abort on any command failure (and on a failure anywhere in a pipe) so a failed
    # minimap2 can't silently fall through to samtools sort.
    set -e -o pipefail

    # cleanup() ALWAYS runs on exit (success, error, or signal) via the EXIT trap,
    # so intermediate/temp files are removed no matter how the task ends. It preserves
    # the original exit code so Nextflow's errorStrategy still sees the real status.
    cleanup() {
        rc=\$?
        rm -f ${prefix}.unsorted.sam ${prefix}.tmp*.bam ${meta.id}.prefix.* 2>/dev/null || true
        exit \$rc
    }
    trap cleanup EXIT

    minimap2 \\
        $args $mapx \\
        -t ${threads} ${I_value} ${K_value} ${split_prefix} \\
        ${reference} \\
        ${input_reads} \\
        ${cigar_paf} ${mmap2_window} ${mmap2_fraction_filter} \\
        ${set_cigar_bam} \\
        ${mm2_out}
    ${sort_step}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}