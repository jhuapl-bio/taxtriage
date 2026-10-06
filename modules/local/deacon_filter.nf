// Host depletion with Deacon (https://github.com/bede/deacon).
//
// Drop-in alternative to the minimap2 + REMOVE_HOSTREADS route: reads sharing
// enough minimizers with the host index are discarded (`deacon filter -d`).
// Outputs use the same names (*.hostremoved.fastq.gz) and the same MultiQC
// stats table as REMOVE_HOSTREADS so everything downstream is unchanged.
process DEACON_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::deacon=0.18.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deacon:0.18.0--h3edb6b3_0' :
        'biocontainers/deacon:0.18.0--h3edb6b3_0' }"

    input:
    tuple val(meta), path(reads)
    path(index)

    output:
    tuple val(meta), path("*.hostremoved.fastq.gz")                , emit: reads
    tuple val(meta), path("*.host_removal_stats_mqc.tsv")         , emit: stats
    tuple val(meta), path("*.deacon.json")                         , emit: summary
    path  "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def abs_t     = params.deacon_abs_threshold != null ? "-a ${params.deacon_abs_threshold}" : ''
    def rel_t     = params.deacon_rel_threshold != null ? "-r ${params.deacon_rel_threshold}" : ''
    def read_list = reads instanceof List ? reads : [reads]
    def paired    = !meta.single_end && read_list.size() == 2
    def outs      = paired ? ["${prefix}_1.hostremoved.fastq.gz", "${prefix}_2.hostremoved.fastq.gz"]
                           : ["${prefix}.hostremoved.fastq.gz"]
    def in_args   = paired ? "${read_list[0]} ${read_list[1]}" : "${read_list[0]}"
    def out_args  = paired ? "-o ${outs[0]}.tmp -O ${outs[1]}.tmp" : "-o ${outs[0]}.tmp"
    """
    set -o pipefail

    deacon filter -d \\
        -t ${task.cpus} \\
        ${abs_t} ${rel_t} \\
        ${args} \\
        -s ${prefix}.deacon.json \\
        ${index} \\
        ${in_args} \\
        ${out_args}

    # Normalise every output to gzipped FASTQ.  Deacon keeps the input format, so
    # FASTA samples get placeholder qualities (the minimap2 route does the same).
    for OUT in ${outs.join(' ')}; do
        TMP=\${OUT}.tmp
        if gzip -t "\$TMP" 2>/dev/null; then READ="gzip -dc"; else READ="cat"; fi
        FIRST=\$(\$READ "\$TMP" | head -c 1 || true)
        if [ "\$FIRST" = ">" ]; then
            \$READ "\$TMP" | awk 'BEGIN{seq=""} /^>/{if(seq!=""){q=seq; gsub(/./,"I",q); print "@"name"\\n"seq"\\n+\\n"q}; name=substr(\$0,2); seq=""; next} {seq=seq\$0} END{if(seq!=""){q=seq; gsub(/./,"I",q); print "@"name"\\n"seq"\\n+\\n"q}}' | gzip -c > "\$OUT"
        elif [ "\$READ" = "cat" ]; then
            gzip -c "\$TMP" > "\$OUT"
        else
            mv "\$TMP" "\$OUT"
        fi
        rm -f "\$TMP"
    done

    # Read counts for the MultiQC table: prefer Deacon's own summary, fall back
    # to counting.  Paired-end counts include both mates, as in REMOVE_HOSTREADS.
    TOTAL=\$(grep -o '"seqs_in"[^0-9]*[0-9]*' ${prefix}.deacon.json | grep -o '[0-9]*\$' || true)
    RETAINED=\$(grep -o '"seqs_out"[^0-9]*[0-9]*' ${prefix}.deacon.json | grep -o '[0-9]*\$' || true)
    if [ -z "\$TOTAL" ]; then
        TOTAL=0
        for f in ${read_list.join(' ')}; do
            if gzip -t "\$f" 2>/dev/null; then R="gzip -dc"; else R="cat"; fi
            if [ "\$(\$R "\$f" | head -c 1)" = ">" ]; then
                n=\$(\$R "\$f" | grep -c '^>' || true)
            else
                n=\$(( \$(\$R "\$f" | wc -l) / 4 ))
            fi
            TOTAL=\$(( TOTAL + n ))
        done
    fi
    if [ -z "\$RETAINED" ]; then
        RETAINED=0
        for f in ${outs.join(' ')}; do
            RETAINED=\$(( RETAINED + \$(gzip -dc "\$f" | wc -l) / 4 ))
        done
    fi

    REMOVED=\$(( TOTAL - RETAINED ))
    PCT=\$(awk -v t="\$TOTAL" -v r="\$REMOVED" 'BEGIN{ if (t>0) printf "%.2f", (r/t)*100; else printf "0.00" }')
    ALL="NO"
    if [ "\$RETAINED" -eq 0 ] && [ "\$TOTAL" -gt 0 ]; then
        ALL="YES"
        echo "WARNING: Deacon classified ALL reads in '${prefix}' as host; the sample falls back to its original reads." >&2
    fi
    printf "Sample\\tTotal Reads\\tRetained Reads\\tRemoved (Host) Reads\\tPercent Host\\tAll Removed\\n" > ${prefix}.host_removal_stats_mqc.tsv
    printf "%s\\t%s\\t%s\\t%s\\t%s%%\\t%s\\n" "${prefix}" "\$TOTAL" "\$RETAINED" "\$REMOVED" "\$PCT" "\$ALL" >> ${prefix}.host_removal_stats_mqc.tsv

    # Empty outputs are kept on purpose: HOST_REMOVAL's CHECK_GZIPPED_READS sees
    # them and falls back to the original reads, as on the minimap2 route.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: \$(deacon --version 2>&1 | sed 's/^deacon //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.hostremoved.fastq.gz
    echo '{}' > ${prefix}.deacon.json
    printf "Sample\\tTotal Reads\\tRetained Reads\\tRemoved (Host) Reads\\tPercent Host\\tAll Removed\\n${prefix}\\t0\\t0\\t0\\t0.00%%\\tNO\\n" > ${prefix}.host_removal_stats_mqc.tsv
    touch versions.yml
    """
}
