process PREPARE_NANOSIM_INPUTS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::seqkit" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(meta), path(abundance_tsv), path(reference_fasta), path(merged_taxid)
    val(num_reads)

    output:
    tuple val(meta), path("genome_list.tsv"), path("size_file.tsv"), emit: nanosim_inputs
    tuple val(meta), path("genomes")     , emit: genomes
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    prepare_nanosim_inputs.py \
        --abundance ${abundance_tsv} \
        --reference ${reference_fasta} \
        --merged_taxid ${merged_taxid} \
        --num_reads ${num_reads} \
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
