process GENERATE_ABUNDANCE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(meta), path(report), path(merged_taxid)
    val(sampletype)
    val (nsamples)

    output:
    tuple val(meta), path("*abu.txt"), emit: abundance
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.abu"



    """
    make_sample_abundance.py \
        --top_report "$report" \
        --merged_taxid "$merged_taxid" \
        --fasta "$prefix".abu.fasta \
        --outdir . \
        --nsamples "$nsamples" \
        --ranks S S1 S2 S3 \


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        insilicoseq: \$(echo \$(iss --version 2>&1) | sed 's/^.*iss v//' ))
    END_VERSIONS
    """
}
