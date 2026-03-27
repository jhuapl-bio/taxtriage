// Align in-silico simulated reads against the reference and produce a
// match_paths.py JSON that can be used as an "insilico control" for
// non-control samples — analogous to lab negative/positive controls.
process ALIGN_INSILICO_READS {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::minimap2 bioconda::samtools bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(meta), path(reads), path(reference_fasta), path(mapping)
    path(pathogens_list)
    file assembly
    val minmapq
    path(taxdump)

    output:
    tuple val(meta), path("*.insilico.json"), emit: json
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}.insilico.json"
    def minmapq_arg = minmapq ? " --minmapq ${minmapq} " : ""
    def assemblyi = assembly ? " -j ${assembly} " : " "
    def mapping_arg = mapping.name != "NO_FILE" ? " -m ${mapping} " : " "
    def taxonomy = taxdump ? " --taxdump ${taxdump} " : " "
    def pathogens_arg = pathogens_list ? " -p ${pathogens_list} " : " "
    def platform = meta.platform ?: "ILLUMINA"
    // Determine minimap2 preset based on platform
    def mm2_preset = platform == "OXFORD" ? "map-ont" : "sr"

    """
    # Index reference
    minimap2 -d ref.mmi ${reference_fasta}

    # Align simulated reads
    minimap2 -ax ${mm2_preset} -t ${task.cpus} ref.mmi ${reads} \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam

    samtools index ${prefix}.sorted.bam

    # Run match_paths.py to produce the JSON (same format as control JSONs)
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
