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
    // Priority: meta.minimap2_preset user-supplied per-sample >
    //           meta.platform-derived default.
    // Valid preset names passed without the leading -ax:
    //   map-ont, map-pb, map-hifi, lr:hq,
    //   sr, splice, splice:hq, splice:sr,
    //   asm5, asm10, asm20,
    //   ava-pb, ava-ont.
    def mapx = ''
    if (meta.minimap2_preset) {
        mapx = "-ax ${meta.minimap2_preset}"
    } else if (meta.platform =~ /(?i)illumina/) {
        mapx = '-ax sr'
    } else if (meta.platform =~ /(?i)pacbio/) {
        mapx = '-ax map-hifi'
    } else {
        // Default: Oxford Nanopore, also used for unspecified / FASTA inputs.
        mapx = '-ax map-ont'
    }

    def input_reads = reads.findAll { it != null }.join(' ')

    /*
     * Dynamic memory / thread model
     *
     * Keep 20% of task.memory unused.
     *
     * Of the remaining 80%:
     *   - samtools sort gets 10%
     *   - minimap2 tuning budget gets 70%
     *   - the remaining 20% is extra headroom for pipes, compression,
     *     process overhead, and memory not directly controlled by -I/-K/-m.
     *
     * Example with task.memory = 20 GB and task.cpus = 6:
     *
     *   usable memory         = 20 GB * 0.80 = 16 GB
     *   samtools sort budget  = 16 GB * 0.10 = 1.6 GB total
     *   minimap2 budget       = 16 GB * 0.70 = 11.2 GB total
     *   samtools threads      = 1
     *   minimap2 threads      = 5
     *
     * Because the BAM path has two samtools sort stages in the same pipe,
     * the samtools sort budget is divided across both sort commands.
     *
     * So in the 20 GB / 6 CPU example:
     *
     *   each samtools sort -m  = 1.6 GB / 2 = 800M
     *   minimap2 -I           = 11200M
     *   minimap2 -K           = 11200M / 5 = 2240M
     *
     * Note:
     *   samtools sort -m is per thread.
     *   minimap2 -I and -K are tuning parameters, not strict RAM caps.
     */
    def toDoubleOrDefault = { value, defaultValue ->
        value != null ? value.toString().toDouble() : defaultValue
    }

    def mem_reserve_fraction = toDoubleOrDefault(params.mmap2_mem_reserve_fraction, 0.20d)
    def sort_budget_fraction = toDoubleOrDefault(params.mmap2_sort_budget_fraction, 0.10d)
    def mm2_budget_fraction  = toDoubleOrDefault(params.mmap2_mm2_budget_fraction, 0.70d)

    def total_mem_mb  = task.memory.toMega().longValue()
    def usable_mem_mb = Math.max((total_mem_mb * (1.0d - mem_reserve_fraction)).longValue(), 1L)

    // minimap2 and samtools now run sequentially (file-based, not piped),
    // so each stage can use all available CPUs in turn.
    def sort_threads = bam_format ? Math.max(task.cpus, 1) : 0
    def mm2_threads  = Math.max(task.cpus, 1)

    // samtools sort -m is per thread.
    // The BAM pipeline has two sort stages:
    //   1. name sort
    //   2. coordinate sort
    //
    // Divide the total samtools sort budget across both stages.
    // Paired-end runs two coordinate-affecting sort stages (name + coord);
    // single-end runs only the single coordinate sort.
    def sort_stages = bam_format ? (meta.single_end ? 1 : 2) : 0

    def sort_budget_total_mb = bam_format
        ? Math.max((usable_mem_mb * sort_budget_fraction).longValue(), 1L)
        : 0L

    def sort_mem_per_thread_mb = bam_format
        ? Math.max((sort_budget_total_mb / (sort_stages * sort_threads)) as long, 1L)
        : 0L

    def S_value = "${sort_mem_per_thread_mb}M"

    // minimap2 memory tuning budget.
    // -I uses the total minimap2 budget.
    // -K is scaled across minimap2 threads.
    //
    // User-supplied params.mmap2_I or params.mmap2_K still override these
    // dynamic defaults.
    def mm2_budget_total_mb = Math.max((usable_mem_mb * mm2_budget_fraction).longValue(), 1L)

    def mm2_budget_per_thread_mb = Math.max(
        (mm2_budget_total_mb / mm2_threads) as long,
        1L
    )

    def I_value = params.mmap2_I ?: "${mm2_budget_total_mb}M"
    def K_value = params.mmap2_K ?: "${mm2_budget_per_thread_mb}M"

    def cigar_paf     = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    def mmap2_window  = params.mmap2_window ? "-w ${params.mmap2_window}" : ''
    def mmap2_fraction_filter = params.mmap2_fraction_filter ? "-f ${params.mmap2_fraction_filter}" : ''
    def split_prefix  = params.no_split_prefix ? "" : "--split-prefix ${meta.id}.prefix"

    // File-based BAM pipeline (no shell pipe).
    //
    // minimap2 writes its SAM to disk, then each samtools stage reads from a
    // file and writes the next file. Running the stages sequentially instead
    // of through a pipe keeps peak memory low: minimap2 releases its index
    // memory before samtools sort allocates its own, which avoids the OOM
    // kills that previously surfaced as "samtools sort: failed to read header
    // from -" (samtools receiving an empty stream from a dead minimap2).
    //
    // Why name-sort -> fixmate -> coordinate-sort instead of a single sort?
    //
    // --split-prefix causes minimap2 to process the reference index in
    // multiple chunks and emit reads in multiple passes. Each pass ends with
    // unmapped reads for that chunk, so the SAM is not globally ordered:
    // RNAME=* unmapped records appear between mapped records from different
    // index splits. A single coordinate sort with limited per-thread memory
    // -m creates temp files and merges them; when the merge is also
    // memory-constrained the final BAM can contain RNAME=* records before some
    // mapped records, which samtools index rejects with "NO_COOR reads not in
    // a single block at the end".
    //
    // The name-sort -> fixmate pipeline fixes this:
    //   1. Name sort groups mates together so fixmate can run.
    //   2. fixmate normalises FLAG bits and RNEXT/PNEXT/TLEN for every read,
    //      including setting RNAME=* / POS=0 for reads whose mate is also
    //      unmapped. This ensures truly unmapped reads carry no residual
    //      reference coordinate that would misplace them during coordinate sort.
    //   3. The final coordinate sort produces a BAM where all RNAME=* records
    //      are in a single block at the end, which samtools index requires.
    //
    // For single-end (and long-read) samples there are no mate pairs, so the
    // name-sort -> fixmate stages are pure no-ops. We skip them entirely and
    // run a single coordinate sort straight from the minimap2 SAM, which is
    // both faster and uses less disk. The name-sort -> fixmate pipeline only
    // runs for paired-end data, where it is required.
    def minimap2_out = bam_format ? "-a -o ${prefix}.aln.sam" : "-o ${prefix}.paf"

    def samtools_block
    if (!bam_format) {
        samtools_block = ""
    } else if (meta.single_end) {
        samtools_block = """
    samtools sort      -@ ${sort_threads} -m ${S_value} -T ${prefix}.csort.tmp -O BAM -o ${prefix}.bam ${prefix}.aln.sam"""
    } else {
        samtools_block = """
    samtools sort   -n -@ ${sort_threads} -m ${S_value} -T ${prefix}.nsort.tmp -O BAM -o ${prefix}.nsort.bam ${prefix}.aln.sam
    samtools fixmate -m -@ ${sort_threads} ${prefix}.nsort.bam ${prefix}.fixmate.bam
    samtools sort      -@ ${sort_threads} -m ${S_value} -T ${prefix}.csort.tmp -O BAM -o ${prefix}.bam ${prefix}.fixmate.bam"""
    }

    """
    set -euo pipefail

    # Always remove intermediate and temp files, whether the task succeeds or
    # fails. The EXIT trap preserves the original exit code, so a minimap2 or
    # samtools failure still propagates to Nextflow.
    cleanup() {
        rm -f ${prefix}.aln.sam ${prefix}.nsort.bam ${prefix}.fixmate.bam
        rm -f ${prefix}.nsort.tmp*.bam ${prefix}.csort.tmp*.bam
        rm -f ${meta.id}.prefix.* 2>/dev/null || true
    }
    trap cleanup EXIT

    minimap2 \\
        $args $mapx \\
        -t ${mm2_threads} -I ${I_value} -K ${K_value} ${split_prefix} \\
        ${reference} \\
        ${input_reads} \\
        ${cigar_paf} ${mmap2_window} ${mmap2_fraction_filter} \\
        ${set_cigar_bam} \\
        ${minimap2_out}
${samtools_block}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}