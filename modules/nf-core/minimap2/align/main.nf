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

    // Sort memory is PER THREAD, so set it explicitly and conservatively.
    // Keep the sort buffer small: on Fusion/S3, the host reference is cached
    // by the Fusion agent (~reference size) AND loaded by minimap2 as an index.
    // A generous sort buffer can push the process over its cgroup limit, causing
    // samtools sort to fork() a merge helper and fail with "Cannot allocate memory".
    // 8% of task memory, capped at 1 GB per thread, is safe across all machine sizes.
    def sort_mem_total_mb      = (task.memory.toMega() * 0.08).longValue()
    def sort_mem_per_thread_mb = Math.max((sort_mem_total_mb / sort_threads) as long, 256L)
    sort_mem_per_thread_mb     = Math.min(sort_mem_per_thread_mb, 1024L)
    def S_value = "${sort_mem_per_thread_mb}M"

    // Treat minimap2 memory knobs as tuning params, not hard caps.
    // On retries, shrink chunk sizes so minimap2 uses less RAM (more passes, but survives OOM).
    // Start at 4G (-I) so a 7 GB reference is processed in 2 passes rather than 1,
    // halving minimap2's peak index footprint.
    // User-supplied params.mmap2_I / mmap2_K always take priority.
    def I_defaults = ['4G', '2G', '1G', '512M']
    def K_defaults = ['500M', '200M', '100M', '50M']
    def attempt_idx = Math.min(task.attempt - 1, I_defaults.size() - 1)
    def I_value = params.mmap2_I ?: I_defaults[attempt_idx]
    def K_value = params.mmap2_K ?: K_defaults[attempt_idx]

    def cigar_paf     = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    def mmap2_window  = params.mmap2_window ? "-w ${params.mmap2_window}" : ''
    def mmap2_fraction_filter = params.mmap2_fraction_filter ? "-f ${params.mmap2_fraction_filter}" : ''
    def split_prefix  = params.no_split_prefix ? "" : "--split-prefix ${meta.id}.prefix"

    // BAM output pipeline.
    //
    // Why name-sort → fixmate → coordinate-sort instead of a single sort?
    //
    // --split-prefix causes minimap2 to process the reference index in
    // multiple chunks and emit reads in multiple passes.  Each pass ends with
    // unmapped reads for that chunk, so the SAM stream delivered to the pipe
    // is not globally ordered: RNAME=* (unmapped) records appear between
    // mapped records from different index splits.  A single coordinate sort
    // with limited per-thread memory (-m) creates temp files and merges them;
    // when the merge is also memory-constrained the final BAM can contain
    // RNAME=* records before some mapped records, which samtools index rejects
    // with "NO_COOR reads not in a single block at the end".
    //
    // The name-sort → fixmate pipeline fixes this:
    //   1. Name sort groups mates together so fixmate can run.
    //   2. fixmate normalises FLAG bits and RNEXT/PNEXT/TLEN for every read,
    //      including setting RNAME=* / POS=0 for reads whose mate is also
    //      unmapped.  This ensures truly unmapped reads carry no residual
    //      reference coordinate that would misplace them during coordinate sort.
    //   3. The final coordinate sort produces a BAM where all RNAME=* records
    //      are in a single block at the end, which samtools index requires.
    //
    // For single-end and long-read samples fixmate is a no-op (passes reads
    // through unchanged) so there is no correctness cost, only a small
    // additional sort pass.
    def bam_output = bam_format
        ? """-a \\
        | samtools sort -n -@ ${sort_threads} -m ${S_value} -T ${prefix}.nsort.tmp \\
        | samtools fixmate -m -@ ${sort_threads} - - \\
        | samtools sort    -@ ${sort_threads} -m ${S_value} -T ${prefix}.csort.tmp -O BAM -o ${prefix}.bam"""
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