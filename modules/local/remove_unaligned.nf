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
    tuple val(meta), path("*.host_removal_stats_mqc.tsv"), emit: stats
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def paired = meta.paired ? "-1 ${input[0]} -2 ${input[1]}" : "-U ${input}"
    // For single-end reads (including ONT), only check bit 4 (read unmapped).
    // Bit 8 (mate unmapped) is never set for unpaired alignments, so -f 12
    // would discard ALL reads and produce an empty fastq.
    // For paired-end, -f 12 requires both mates unmapped; -f 4 keeps singletons.
    def flag = ""
    if (meta.single_end) {
        flag = "-f 4 -F 0x900"
    } else {
        flag = !params.include_singletons_hostremoval ? "-f 12 -F 0x900" : "-f 4 -F 0x900"
        if ((params.use_denovo || params.use_diamond) && params.include_singletons_hostremoval) {
            println "ALERT: paired unmapped only (include_singletons_hostremoval parameter) has been ignored because --use_denovo or --use_diamond was specified with paired-end reads. Both reads in a pair must be unmapped to be retained as Megahit fails without equal read counts."
            flag = "-f 12 -F 0x900"
        }
    }

    def min_mapq = params.min_mapq_host ? params.min_mapq_host as int : 0
    def singleton_out = params.include_singletons_removal ? "${meta.id}.hostremoved.singletons.fastq.gz" : "/dev/null"
    def cmd = !meta.single_end ? \
        " -1 ${meta.id}_1.hostremoved.fastq.gz -2 ${meta.id}_2.hostremoved.fastq.gz -0 ${singleton_out}" \
        : "-0 ${meta.id}.hostremoved.fastq.gz -s /dev/null"
    // When min_mapq_host > 0, use expression-based filtering to also keep
    // reads with MAPQ below the threshold (not confidently aligned to host).
    // NOTE: use bare `mapq` (the SAM column), NOT `[mapq]` which would look
    // for an auxiliary tag and cause a samtools parse error.
    def view_filter = ""
    if (min_mapq > 0) {
        if (meta.single_end) {
            // Single-end mode: keep unmapped reads OR reads with MAPQ below threshold
            view_filter = "-F 0x900 -e 'flag.unmap || mapq < ${min_mapq}'"
        } else if (flag.contains("-f 12")) {
            // Paired mode: keep reads where both mates unmapped OR MAPQ below threshold
            view_filter = "-F 0x900 -e '(flag & 12) == 12 || mapq < ${min_mapq}'"
        } else {
            // Paired singleton mode: keep unmapped reads OR reads with MAPQ below threshold
            view_filter = "-F 0x900 -e 'flag.unmap || mapq < ${min_mapq}'"
        }
    } else {
        view_filter = flag
    }
    def prefix = meta.id
    """
    # Count total reads in the input BAM (exclude supplementary/secondary via -F 0x900)
    total_reads=\$(samtools view -c -F 0x900 ${input})

    # Extract non-host reads from BAM to fastq
    samtools view -b ${view_filter} ${input} | \\
        samtools fastq -n ${args} \\
        ${cmd} -

    # Count retained reads from the output fastq(s)
    retained_reads=0
    for fq in ${prefix}*.hostremoved.fastq.gz ${prefix}*.hostremoved.fq.gz; do
        if [ -f "\$fq" ]; then
            fq_lines=\$(gzip -dc "\$fq" | wc -l)
            retained_reads=\$(( retained_reads + fq_lines / 4 ))
        fi
    done

    removed_reads=\$(( total_reads - retained_reads ))
    if [ "\$total_reads" -gt 0 ]; then
        pct_removed=\$(awk "BEGIN {printf \\"%.2f\\", (\$removed_reads / \$total_reads) * 100}")
    else
        pct_removed="0.00"
    fi

    # Log to stdout
    echo "==========================================="
    echo "  HOST REMOVAL STATS: ${prefix}"
    echo "==========================================="
    echo "  Total input reads:    \$total_reads"
    echo "  Retained (non-host):  \$retained_reads"
    echo "  Removed (host):       \$removed_reads (\${pct_removed}%)"
    if [ "\$retained_reads" -eq 0 ] && [ "\$total_reads" -gt 0 ]; then
        echo ""
        echo "  WARNING: ALL reads were classified as host for sample '${prefix}'."
        echo "  The sample will fall back to its original (unfiltered) reads."
        echo "  Consider checking host reference or adjusting --min_mapq_host."
        echo ""
    fi
    echo "==========================================="

    # Write MultiQC-compatible TSV stats file
    printf "Sample\\tTotal Reads\\tRetained Reads\\tRemoved (Host) Reads\\tPercent Host\\tAll Removed\\n" > ${prefix}.host_removal_stats_mqc.tsv
    if [ "\$retained_reads" -eq 0 ] && [ "\$total_reads" -gt 0 ]; then
        all_removed="YES"
    else
        all_removed="NO"
    fi
    printf "${prefix}\\t\$total_reads\\t\$retained_reads\\t\$removed_reads\\t\${pct_removed}%%\\t\$all_removed\\n" >> ${prefix}.host_removal_stats_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
