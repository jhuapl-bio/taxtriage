// Generate an in-silico control JSON directly from the simulated sample
// reference sequences — without running a full read simulator or aligner.
//
// Strategy: self-align each per-sample reference FASTA to the shared
// reference using minimap2 (very fast, one alignment per sequence) to
// produce a valid BAM, then run match_paths.py on that BAM to create
// the insilico control JSON in the standard format consumed by
// ALIGNMENT_PER_SAMPLE.
process GENERATE_INSILICO_JSON {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::minimap2 bioconda::samtools bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(meta), path(sample_reference), path(shared_reference), path(mapping)
    path(pathogens_list)
    file assembly
    val minmapq
    path(taxdump)

    output:
    tuple val(meta), path("*.insilico.json"), emit: json
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}.insilico.json"
    def minmapq_arg = minmapq ? " --minmapq ${minmapq} " : ""
    def assemblyi = assembly ? " -j ${assembly} " : " "
    def mapping_arg = mapping.name != "NO_FILE" ? " -m ${mapping} " : " "
    def taxonomy = taxdump.name != "NO_FILE" ? " --taxdump ${taxdump} " : " "
    def pathogens_arg = pathogens_list ? " -p ${pathogens_list} " : " "
    def platform = meta.platform ?: "ILLUMINA"
    def mm2_preset = platform == "OXFORD" ? "map-ont" : "sr"

    """
    # Index shared reference
    minimap2 -d ref.mmi ${shared_reference}

    # Self-align per-sample reference sequences to the shared reference.
    # This produces one near-perfect alignment per reference contig — a
    # lightweight stand-in for simulated reads that still gives
    # match_paths.py a valid BAM to work with.
    minimap2 -ax ${mm2_preset} -t ${task.cpus} ref.mmi ${sample_reference} \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam

    samtools index ${prefix}.sorted.bam

    # Run match_paths.py to produce the insilico control JSON
    match_paths.py \\
        -i ${prefix}.sorted.bam \\
        -o ${output} \\
        -s ${prefix}_insilico \\
        ${mapping_arg} ${minmapq_arg} ${pathogens_arg} ${assemblyi} \\
        --output_dir search_results \\
        --scaled 8000 \\
        --min_threshold 0.002 \\
        --fast \\
        --platform ${platform} \\
        ${taxonomy} \\
        --control_type positive

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
