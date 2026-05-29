process REMOVE_HOSTREADS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    // reads is the original fastp-trimmed FASTQ(s) — needed so the QNAME-based
    // fallback can filter the original files when the BAM lacks paired-end flags.
    // Single-end: a single path.  Paired-end: [R1, R2].
    tuple val(meta), path(bam), path(reads)

    output:
    tuple val(meta), path("*.hostremoved.{fastq,fq}.gz"), optional: true, emit: reads
    tuple val(meta), path("*.host_removal_stats_mqc.tsv"), emit: stats
    path  "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = meta.id
    def min_mapq     = params.min_mapq_host     ? params.min_mapq_host as int : 0
    def singleton_out = params.include_singletons_removal \
                            ? "${prefix}.hostremoved.singletons.fastq.gz" \
                            : "/dev/null"

    def single_end_flag        = meta.single_end                                ? "-e" : ""
    def include_singletons_flag = params.include_singletons_hostremoval         ? "-i" : ""
    def denovo_flag            = (params.use_denovo || params.use_diamond)      ? "-d" : ""
    def extra_args_flag        = args ? "-x '${args}'" : ""

    def reads_flags = meta.single_end ? "" : "-1 ${reads[0]} -2 ${reads[1]}"

    if ((params.use_denovo || params.use_diamond) && params.include_singletons_hostremoval) {
        log.warn "ALERT: --include_singletons_hostremoval ignored because --use_denovo or " +
                 "--use_diamond was set; both reads in a pair must be unmapped for Megahit."
    }

    """
    remove_hostreads.sh \\
        -b ${bam} \\
        ${reads_flags} \\
        -p ${prefix} \\
        -q ${min_mapq} \\
        -s ${singleton_out} \\
        ${single_end_flag} \\
        ${include_singletons_flag} \\
        ${denovo_flag} \\
        ${extra_args_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
