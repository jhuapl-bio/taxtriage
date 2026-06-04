process MINIMAP2_ALIGN {

    tag "$meta.id"
    label 'process_standard'

    /*
     * Retry likely OOM / signal failures with the larger memory
     * you define in config. Keep this here or move it to nextflow.config.
     */
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

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

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    /*
     * Shell quoting for staged Nextflow paths.
     */
    def shellQuote = { obj ->
        def s = obj.toString()
        return "'" + s.replace("'", "'\"'\"'") + "'"
    }

    /*
     * Reads may arrive as a single Path, a List<Path>, or a list with nulls.
     */
    def read_list = reads instanceof List ? reads.findAll { it != null } : [reads].findAll { it != null }
    def input_reads = read_list.collect { shellQuote(it) }.join(' ')

    /*
     * Reference may be a single FASTA/MMI or a staged list with sidecars.
     * Prefer .mmi if present, otherwise FASTA-like files, otherwise first path.
     */
    def ref_list = reference instanceof List ? reference.findAll { it != null } : [reference].findAll { it != null }
    def reference_file = ref_list.find { it.getName() ==~ /.*\.mmi$/ } ?:
                         ref_list.find { it.getName() ==~ /.*\.(fa|fasta|fna)(\.gz)?$/ } ?:
                         ref_list[0]
    def reference_arg = shellQuote(reference_file)
    def reference_is_mmi = reference_file.getName().endsWith('.mmi')

    /*
     * Preset selection.
     *
     * Important: use -x only here.
     * Add -a only in BAM/SAM mode.
     */
    def platform = (meta.platform ?: '').toString()
    def preset = meta.minimap2_preset ?: (
        (platform =~ /(?i)illumina/) ? 'sr' :
        (platform =~ /(?i)pacbio/)   ? 'map-hifi' :
                                       'map-ont'
    )
    def preset_arg = "-x ${preset}"

    /*
     * Paired-end Illumina gets mate repair.
     * ONT/PacBio/single-end skip fixmate to save memory and I/O.
     */
    def is_paired   = read_list.size() == 2
    def is_illumina = platform ? (platform =~ /(?i)illumina/) : false
    def needs_fixmate = bam_format && is_paired && is_illumina

    /*
     * CPU budgeting.
     *
     * minimap2 can use an extra I/O thread during mapping, and samtools -@
     * creates additional worker/compression threads. Keep overlapped streaming
     * stages conservative; let coordinate sort use a small number of threads
     * after minimap2 has exited.
     */
    def cpus = Math.max((task.cpus ?: 1) as int, 1)

    def stream_threads_req = params.samtools_stream_threads != null ? (params.samtools_stream_threads as int) : 0
    def stream_threads = Math.max(0, Math.min(stream_threads_req, Math.max(cpus - 2, 0)))

    def sort_threads_req = params.samtools_sort_threads != null ? (params.samtools_sort_threads as int) : Math.min(Math.max(cpus - 1, 0), 2)
    def sort_threads = Math.max(0, Math.min(sort_threads_req, Math.max(cpus - 1, 0)))

    /*
     * Reserve one CPU for the minimap2 I/O thread and one for the low-memory
     * samtools view process that runs in the pipe.
     */
    def mm2_threads = Math.max(1, cpus - stream_threads - 2)

    /*
     * task.memory must be set in config for this to be meaningful.
     * If not set, this fallback only protects the calculation; it does not
     * reserve memory from the scheduler.
     */
    def task_mem_mb = task.memory ? (task.memory.toMega() as long) : ((params.minimap2_fallback_memory_mb ?: 8192) as long)

    /*
     * samtools sort -m is approximately memory PER sort thread.
     * Keep total sort memory small because minimap2/index memory is the main risk.
     *
     * Expert override:
     *   --samtools_sort_mem 512M
     *
     * Safer numeric cap:
     *   --samtools_sort_mem_cap_mb 768
     */
    def sort_threads_for_mem = Math.max(sort_threads, 1)
    def sort_budget_mb = Math.max((task_mem_mb * 0.20) as long, 64L)
    def auto_sort_mem_mb = Math.floor((sort_budget_mb / sort_threads_for_mem) as double) as long
    auto_sort_mem_mb = Math.max(auto_sort_mem_mb, 64L)
    auto_sort_mem_mb = Math.min(auto_sort_mem_mb, ((params.samtools_sort_mem_cap_mb ?: 768) as long))

    def S_value = params.samtools_sort_mem ?: "${auto_sort_mem_mb}M"

    /*
     * minimap2 memory knobs.
     *
     * -I is target bases per index batch, not a hard memory cap.
     * Very small -I can create a multi-part index. That reduces peak memory but
     * can affect MAPQ. If the reference is an existing .mmi, minimap2 uses the
     * index settings stored in the .mmi, so -I/-k/-w cannot shrink it here.
     */
    def auto_I = task_mem_mb >= 24576 ? '8G' :
                 task_mem_mb >= 12288 ? '4G' :
                 task_mem_mb >=  6144 ? '2G' :
                                         '1G'

    def I_value = params.minimap2_I ?: params.mmap2_I ?: auto_I
    def K_value = params.minimap2_K ?: params.mmap2_K ?: '25M'

    def mm2_index_opts = reference_is_mmi ? "" : "-I ${I_value}"

    def window_value = params.minimap2_window ?: params.mmap2_window
    def mm2_window = window_value ? "-w ${window_value}" : ''

    def fraction_value = params.minimap2_fraction_filter ?: params.mmap2_fraction_filter
    def mm2_fraction_filter = fraction_value ? "-f ${fraction_value}" : ''

    def split_prefix = params.no_split_prefix ? "" : "--split-prefix \${TMPDIR}/${prefix}.mm2split"

    def paf_cigar = cigar_paf_format && !bam_format ? "-c" : ''
    def bam_long_cigar = cigar_bam && bam_format ? "-L" : ''

    def sort_level = params.samtools_sort_level
    def sort_level_opt = sort_level != null ? "-l ${sort_level}" : ""

    def collate_temp_files = params.samtools_collate_temp_files ?: 64

    /*
     * Keep the minimap2 command options before target/query.
     * Put safety-controlled options after task.ext.args so accidental -t/-K
     * in args is less likely to override the safe values.
     */
    def mm2_common = """
        minimap2 \\
            ${args} \\
            ${preset_arg} \\
            -t ${mm2_threads} \\
            ${mm2_index_opts} \\
            -K ${K_value} \\
            ${split_prefix} \\
            ${mm2_window} \\
            ${mm2_fraction_filter}
    """

    def run_alignment

    if (bam_format && needs_fixmate) {

        /*
         * OOM-safe paired-end path:
         *
         *   minimap2 | low-memory samtools view -> unsorted BAM
         *   samtools collate -> name-grouped BAM
         *   samtools fixmate -> mate-fixed BAM
         *   samtools sort -> final coordinate BAM
         *
         * This avoids the high-risk resident chain:
         *   minimap2 | sort -n | fixmate | sort
         */
        run_alignment = """
        ${mm2_common} \\
            -a \\
            ${bam_long_cigar} \\
            ${reference_arg} \\
            ${input_reads} \\
        | samtools view \\
            -@ ${stream_threads} \\
            -u \\
            -o "${prefix}.unsorted.bam" \\
            -

        samtools collate \\
            -@ ${stream_threads} \\
            -n ${collate_temp_files} \\
            -u \\
            -O \\
            "${prefix}.unsorted.bam" \\
            "\${TMPDIR}/${prefix}.collate.tmp" \\
        > "${prefix}.collate.bam"

        samtools fixmate \\
            -m \\
            -@ ${stream_threads} \\
            -u \\
            "${prefix}.collate.bam" \\
            "${prefix}.fixmate.bam"

        samtools sort \\
            ${sort_level_opt} \\
            -@ ${sort_threads} \\
            -m ${S_value} \\
            -T "\${TMPDIR}/${prefix}.csort.tmp" \\
            -O BAM \\
            -o "${prefix}.bam" \\
            "${prefix}.fixmate.bam"

        samtools quickcheck -v "${prefix}.bam"
        """

    } else if (bam_format) {

        /*
         * OOM-safe non-fixmate BAM path:
         *
         * Do not overlap minimap2 and coordinate sort.
         * That costs disk but prevents minimap2 index/mapping memory and
         * samtools sort memory from peaking at the same time.
         */
        run_alignment = """
        ${mm2_common} \\
            -a \\
            ${bam_long_cigar} \\
            ${reference_arg} \\
            ${input_reads} \\
        | samtools view \\
            -@ ${stream_threads} \\
            -u \\
            -o "${prefix}.unsorted.bam" \\
            -

        samtools sort \\
            ${sort_level_opt} \\
            -@ ${sort_threads} \\
            -m ${S_value} \\
            -T "\${TMPDIR}/${prefix}.csort.tmp" \\
            -O BAM \\
            -o "${prefix}.bam" \\
            "${prefix}.unsorted.bam"

        samtools quickcheck -v "${prefix}.bam"
        """

    } else {

        /*
         * PAF path:
         *
         * No -a here. The preset is -x only, so this really writes PAF.
         */
        run_alignment = """
        ${mm2_common} \\
            ${paf_cigar} \\
            -o "${prefix}.paf" \\
            ${reference_arg} \\
            ${input_reads}
        """
    }

    """
    set -euo pipefail

    export TMPDIR="\${TMPDIR:-\$PWD/tmp}"
    mkdir -p "\${TMPDIR}"

    # Reduces glibc arena growth in multi-threaded tools.
    export MALLOC_ARENA_MAX="${params.malloc_arena_max ?: 2}"

    cleanup() {
        rm -f "${prefix}.unsorted.bam" "${prefix}.collate.bam" "${prefix}.fixmate.bam" "${prefix}.name.bam" 2>/dev/null || true
        rm -f "\${TMPDIR}/${prefix}".*.tmp* "\${TMPDIR}/${prefix}.mm2split"* 2>/dev/null || true
    }
    trap cleanup EXIT

    ${run_alignment}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/^samtools //')
    END_VERSIONS
    """
}