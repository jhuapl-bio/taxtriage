process SUBSAMPLE_INSILICO {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(meta), path(reads)
    val(mode)
    val(counts)        // comma/space separated list of read counts
    val(replicates)
    val(seed)

    output:
    tuple val(meta), path("datasets/*.fastq.gz"), emit: reads
    tuple val(meta), path("indices/*.idx.gz")   , emit: indices
    path("*_subsample_manifest.tsv")            , emit: manifest
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // NOTE: FASTQs are always materialised here because downstream ALIGNMENT
    // consumes them. Whether they are *published* (kept on disk) is controlled
    // by sim_keep_subsampled_fastq in conf/modules.config. Index files always
    // allow reconstruction via bin/reconstruct_insilico_reads.py.
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = !meta.single_end

    // Build the master-mate preparation snippet in Groovy so we avoid nested
    // string-escaping pitfalls in the shell heredoc below.
    def mate_prep
    if (paired) {
        mate_prep = """\
        R1=\$( (ls *_R1*.fastq.gz *_1.fastq.gz 2>/dev/null || true) | head -n1 )
        R2=\$( (ls *_R2*.fastq.gz *_2.fastq.gz 2>/dev/null || true) | head -n1 )
        if [ -z "\$R1" ] || [ -z "\$R2" ]; then
            echo "ERROR: could not identify paired R1/R2 master FASTQs" >&2
            exit 1
        fi
        MATE_ARGS="--r1 \$R1 --r2 \$R2\""""
    } else {
        mate_prep = """\
        cat ${reads} > master.se.fastq.gz
        MATE_ARGS="--r1 master.se.fastq.gz\""""
    }

    """
    mkdir -p datasets indices

    ${mate_prep}

    subsample_insilico_reads.py \\
        \$MATE_ARGS \\
        --mode ${mode} \\
        --counts "${counts}" \\
        --replicates ${replicates} \\
        --seed ${seed} \\
        --parent ${prefix} \\
        --outdir datasets \\
        --index-dir indices \\
        --manifest ${prefix}_subsample_manifest.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
