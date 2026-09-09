// Build every dataset of a spike-in series: exact read counts drawn from each
// organism's simulated pool, concatenated onto the FIXED background.
//
// The background is byte-identical in every dataset — that is what makes this a
// spike-in series rather than a dilution — so it is concatenated in the shell
// (gzip members concatenate cleanly) instead of being read through Python. Memory
// stays flat no matter how deep the background is.
process SPIKE_INTO_BACKGROUND {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(meta), path(reads)          // the background sample's cleaned reads
    path(spikein_tsv)                     // normalised sheet from PARSE_SPIKEIN
    path(pools, stageAs: "pools/*")       // every accession's simulated pool
    val(mode)
    val(seed)

    output:
    tuple val(meta), path("datasets/*.fastq.gz"), emit: reads
    path("*_spikein_manifest.tsv")               , emit: manifest
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = !meta.single_end
    def paired_arg = paired ? "--paired" : ""
    def bg_depth = params.spikein_background_depth ?: 0
    // Identify the background mates the same way SUBSAMPLE_INSILICO does, so a
    // sample that reaches here through either path behaves identically.
    // Use the STAGED PATHS Nextflow gives us, not a glob. Cleaned reads reach this
    // process named by whatever produced them (fastp emits background_1.fastp.fastq.gz,
    // which matches neither *_R1*.fastq.gz nor *_1.fastq.gz), so globbing for mates is
    // guesswork that fails on exactly the inputs this process is meant to take.
    def _rl = (reads instanceof List) ? reads : [reads]
    def mate_prep = paired ? """\
    BG_R1=${_rl[0]}
    BG_R2=${_rl.size() > 1 ? _rl[1] : ''}
    if [ ! -s "\$BG_R1" ] || [ ! -s "\$BG_R2" ]; then
        echo "ERROR: paired background FASTQs missing or empty (\$BG_R1 / \$BG_R2)" >&2
        exit 1
    fi""" : """\
    cat ${reads} > background.se.fastq.gz
    BG_R1=background.se.fastq.gz
    BG_R2=''"""

    """
    mkdir -p datasets spike

    ${mate_prep}

    # Optionally trim the background to a fixed depth first, so every dataset in
    # the series sits on a background of known, comparable size.
    if [ "${bg_depth}" -gt 0 ] 2>/dev/null; then
        head -n \$(( ${bg_depth} * 4 )) <(zcat "\$BG_R1") | gzip > bg.trim_R1.fastq.gz
        BG_R1=bg.trim_R1.fastq.gz
        if [ -n "\$BG_R2" ] && [ "${paired}" = "true" ]; then
            head -n \$(( ${bg_depth} * 4 )) <(zcat "\$BG_R2") | gzip > bg.trim_R2.fastq.gz
            BG_R2=bg.trim_R2.fastq.gz
        fi
    fi

    BG_READS=\$(( \$(zcat "\$BG_R1" | wc -l) / 4 ))
    echo "[spike] background \$BG_R1 holds \$BG_READS records" >&2

    # Map each staged pool file back to its accession. Pool files are named
    # <accession>[_R1|_R2].fastq.gz by SPIKEIN_POOL.
    POOL_ARGS=""
    for acc in \$(cut -f3 ${spikein_tsv} | tail -n +2 | sort -u); do
        # printf, NOT echo: echo appends a newline, tr -c maps it (it is outside the
        # allowed set) to '_', and the name gains a trailing underscore — so
        # GCF_000859985.2 became GCF_000859985.2_ and never matched the pool file
        # SPIKEIN_POOL wrote from the Groovy-side sanitisation.
        safe=\$(printf '%s' "\$acc" | tr -c 'A-Za-z0-9._-' '_')
        if [ "${paired}" = "true" ]; then
            p1=pools/\${safe}_R1.fastq.gz
            p2=pools/\${safe}_R2.fastq.gz
            [ -f "\$p1" ] && [ -f "\$p2" ] || { echo "ERROR: missing pool for \$acc (\$p1 / \$p2)" >&2; exit 1; }
            POOL_ARGS="\$POOL_ARGS --pool \$acc=\$p1,\$p2"
        else
            p1=pools/\${safe}.fastq.gz
            [ -f "\$p1" ] || { echo "ERROR: missing pool for \$acc (\$p1)" >&2; exit 1; }
            POOL_ARGS="\$POOL_ARGS --pool \$acc=\$p1"
        fi
    done

    spike_into_background.py \\
        --spikein ${spikein_tsv} \\
        --parent ${prefix} \\
        --mode ${mode} \\
        --seed ${seed} \\
        ${paired_arg} \\
        \$POOL_ARGS \\
        --outdir spike \\
        --manifest ${prefix}_spikein_manifest.tsv \\
        --background-reads \$BG_READS \\
        --background-name ${meta.parent_id ?: meta.id}

    # Background + spike -> the dataset FASTQs the pipeline will analyse.
    for f in spike/*.spike_R1.fastq.gz; do
        [ -e "\$f" ] || continue
        ds=\$(basename "\$f" .spike_R1.fastq.gz)
        cat "\$BG_R1" "\$f" > datasets/\${ds}.spikein_R1.fastq.gz
        cat "\$BG_R2" "spike/\${ds}.spike_R2.fastq.gz" > datasets/\${ds}.spikein_R2.fastq.gz
    done
    for f in spike/*.spike.fastq.gz; do
        [ -e "\$f" ] || continue
        ds=\$(basename "\$f" .spike.fastq.gz)
        cat "\$BG_R1" "\$f" > datasets/\${ds}.spikein.fastq.gz
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
