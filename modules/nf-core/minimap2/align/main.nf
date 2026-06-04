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

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    /*
     * Resolve minimap2 preset.
     *
     * Priority:
     *   1. meta.minimap2_preset
     *   2. platform-derived default
     */
    def mapx = ''
    if (meta.minimap2_preset) {
        mapx = "-ax ${meta.minimap2_preset}"
    } else if (meta.platform =~ /(?i)illumina/) {
        mapx = '-ax sr'
    } else if (meta.platform =~ /(?i)pacbio/) {
        mapx = '-ax map-hifi'
    } else {
        mapx = '-ax map-ont'
    }

    /*
     * Make reads robust whether `reads` arrives as:
     *   - a List<Path>
     *   - a single Path
     *   - a list containing nulls
     */
    def read_list = reads instanceof List ? reads.findAll { it != null } : [reads].findAll { it != null }
    def input_reads = read_list.join(' ')

    /*
     * Paired-end detection.
     *
     * Only use name-sort -> fixmate -> coordinate-sort for paired-end Illumina.
     * For ONT, PacBio, FASTA, or single-end reads, skip fixmate to save memory.
     */
    def is_paired = read_list.size() == 2
    def is_illumina = meta.platform ? (meta.platform =~ /(?i)illumina/) : false
    def needs_fixmate = bam_format && is_paired && is_illumina

    /*
     * Threading strategy.
     *
     * Keep samtools sort small. minimap2 gets the remaining CPUs.
     * This avoids memory spikes from multiple high-threaded tools running together.
     */
    def sort_threads = 1
    def mm2_threads  = Math.max(task.cpus - sort_threads, 1)

    /*
     * samtools sort memory.
     *
     * `samtools sort -m` is memory PER THREAD.
     * Keep it conservative because minimap2 and samtools sort overlap during piping.
     *
     * Default behavior:
     *   - use about 10% of task memory for the active samtools sort
     *   - minimum 256M
     *   - maximum 1024M per sort thread
     *
     * Override with:
     *   --samtools_sort_mem 768M
     */
    def sort_mem_per_thread_mb = Math.max(((task.memory.toMega() * 0.10) / sort_threads) as long, 256L)
    sort_mem_per_thread_mb = Math.min(sort_mem_per_thread_mb, 1024L)
    def S_value = params.samtools_sort_mem ?: "${sort_mem_per_thread_mb}M"

    /*
     * minimap2 memory-related knobs.
     *
     * These are not hard caps, but lowering them reduces peak memory.
     *
     * Override with:
     *   --mmap2_I 4G
     *   --mmap2_K 50M
     */
    def I_value = params.mmap2_I ?: '2G'
    def K_value = params.mmap2_K ?: '25M'

    def cigar_paf = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    def mmap2_window = params.mmap2_window ? "-w ${params.mmap2_window}" : ''
    def mmap2_fraction_filter = params.mmap2_fraction_filter ? "-f ${params.mmap2_fraction_filter}" : ''
    def split_prefix = params.no_split_prefix ? "" : "--split-prefix ${prefix}.prefix"

    /*
     * Build output command.
     *
     * Important:
     *   The previous implementation streamed:
     *
     *     minimap2 | sort -n | fixmate | sort
     *
     *   That can have minimap2 plus two samtools sort processes resident at once.
     *
     *   This version writes intermediate BAMs for paired-end Illumina:
     *
     *     minimap2 | name-sort -> name.bam
     *     fixmate name.bam -> fixmate.bam
     *     coordinate-sort fixmate.bam -> final.bam
     *
     *   This lowers peak RAM substantially.
     */
    def run_alignment

    if (bam_format && needs_fixmate) {
        run_alignment = """
        minimap2 \\
            ${args} ${mapx} \\
            -t ${mm2_threads} -I ${I_value} -K ${K_value} ${split_prefix} \\
            ${reference} \\
            ${input_reads} \\
            ${mmap2_window} ${mmap2_fraction_filter} \\
            ${set_cigar_bam} \\
            -a \\
        | samtools sort \\
            -n \\
            -@ ${sort_threads} \\
            -m ${S_value} \\
            -T \${TMPDIR:-.}/${prefix}.nsort.tmp \\
            -O BAM \\
            -o ${prefix}.name.bam \\
            -

        samtools fixmate \\
            -m \\
            -@ 1 \\
            ${prefix}.name.bam \\
            ${prefix}.fixmate.bam

        samtools sort \\
            -@ ${sort_threads} \\
            -m ${S_value} \\
            -T \${TMPDIR:-.}/${prefix}.csort.tmp \\
            -O BAM \\
            -o ${prefix}.bam \\
            ${prefix}.fixmate.bam

        rm -f ${prefix}.name.bam ${prefix}.fixmate.bam
        """
    } else if (bam_format) {
        run_alignment = """
        minimap2 \\
            ${args} ${mapx} \\
            -t ${mm2_threads} -I ${I_value} -K ${K_value} ${split_prefix} \\
            ${reference} \\
            ${input_reads} \\
            ${mmap2_window} ${mmap2_fraction_filter} \\
            ${set_cigar_bam} \\
            -a \\
        | samtools sort \\
            -@ ${sort_threads} \\
            -m ${S_value} \\
            -T \${TMPDIR:-.}/${prefix}.csort.tmp \\
            -O BAM \\
            -o ${prefix}.bam \\
            -
        """
    } else {
        run_alignment = """
        minimap2 \\
            ${args} ${mapx} \\
            -t ${mm2_threads} -I ${I_value} -K ${K_value} ${split_prefix} \\
            ${reference} \\
            ${input_reads} \\
            ${cigar_paf} ${mmap2_window} ${mmap2_fraction_filter} \\
            -o ${prefix}.paf
        """
    }

    """
    set -euo pipefail

    export TMPDIR="\${TMPDIR:-.}"

    cleanup() {
        rm -f ${prefix}.name.bam ${prefix}.fixmate.bam || true
        rm -f \${TMPDIR}/${prefix}.nsort.tmp* || true
        rm -f \${TMPDIR}/${prefix}.csort.tmp* || true
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