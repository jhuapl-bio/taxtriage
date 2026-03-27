process MAKE_SIMULATED_SAMPLES {
    tag "$meta.id"
    label 'process_medium'

    // conda (params.enable_conda ? "bioconda::" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(meta), path(top_report), path(merged_taxid), path(fastas)
    val(nsamples)
    val(ranks)
    val(minreads)
    val(exclude_taxids)
    val(include_taxids)
    path(abundance_input)

    output:
    tuple val(meta), path("sample_*/abundance.tsv"), path("sample_*/reference.fasta"), emit: samples
    tuple val(meta), path("reference_sequences.fasta"), emit: shared_reference
    path("manifest.tsv")       , emit: manifest
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def exclude_arg = exclude_taxids ? "--exclude_taxids ${exclude_taxids}" : ""
    def include_arg = include_taxids ? "--include_taxids ${include_taxids}" : ""
    def ranks_arg = ranks ? "--ranks ${ranks}" : "--ranks S S1 S2 S3"
    def random_abu_arg = params.sim_random_abundance ? "--random_abundance" : ""
    def nreads_arg = params.sim_nreads ? "--nreads ${params.sim_nreads}" : ""
    // If custom abundance_input is provided (not the NO_FILE sentinel), use it instead of top_report
    def has_custom_abundance = abundance_input.name != 'NO_FILE'
    def abundance_input_arg = has_custom_abundance ? "--abundance_input ${abundance_input}" : ""
    def top_report_arg = has_custom_abundance ? "" : "--top_report ${top_report}"

    """
    # Concatenate all reference FASTAs into a single file
    cat ${fastas} > merged_reference.fasta

    make_simulated_samples.py \
        ${top_report_arg} \
        --merged_taxid ${merged_taxid} \
        --fasta merged_reference.fasta \
        --outdir . \
        --nsamples ${nsamples} \
        ${ranks_arg} \
        --minreads ${minreads} \
        ${exclude_arg} \
        ${include_arg} \
        ${abundance_input_arg} \
        ${random_abu_arg} \
        ${nreads_arg} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
