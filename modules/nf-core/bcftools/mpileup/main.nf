process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5acacb55c52bec97c61fd34ffa8721fce82ce823005793592e2a80bf71632cd0/data':
        'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11' }"

    input:
    tuple val(meta), path(bam), path(intervals)
    tuple val(meta2), path(fasta)
    val save_mpileup

    output:
    tuple val(meta), path("*vcf.gz")     , emit: vcf
    tuple val(meta), path("*vcf.gz.tbi") , emit: tbi
    tuple val(meta), path("*stats.txt")  , emit: stats
    tuple val(meta), path("*.mpileup.gz"), emit: mpileup, optional: true
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mpileup = save_mpileup ? "| tee ${prefix}.mpileup" : ""
    def bgzip_mpileup = save_mpileup ? "bgzip ${prefix}.mpileup" : ""
    def intervals_cmd = intervals ? "-T ${intervals}" : ""
    // --fasta-ref must be faidx-able (plain FASTA). A gzip/bgzip reference is
    // unpacked inside this task only (detected by magic bytes, not extension).
    """
    echo "${meta.id}" > sample_name.list

    REF="$fasta"
    if [ "\$(head -c 2 "\$REF" | od -An -tx1 | tr -d ' \\n')" = "1f8b" ]; then
        gzip -dc "\$REF" > ${prefix}.mpileup_ref.fa
        REF=${prefix}.mpileup_ref.fa
    fi

    bcftools \\
        mpileup \\
        --fasta-ref \$REF \\
        $args \\
        $bam \\
        $intervals_cmd \\
        $mpileup \\
        | bcftools call --output-type v $args2 \\
        | bcftools reheader --samples sample_name.list \\
        | bcftools view --output-file ${prefix}.vcf.gz --output-type z $args3

    $bgzip_mpileup

    tabix -p vcf -f ${prefix}.vcf.gz
    rm -f ${prefix}.mpileup_ref.fa ${prefix}.mpileup_ref.fa.fai

    bcftools stats ${prefix}.vcf.gz > ${prefix}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcftools_stats.txt
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.mpileup.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
