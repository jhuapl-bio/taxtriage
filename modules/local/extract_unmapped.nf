// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # (header trimmed for the sketch -- copy the full APL header from the other modules)
// ##############################################################################################
//
// Read accounting for the novelty score, in one place. Counts only -- the taxonomy branch
// now classifies the FULL de novo contigs (built from all post-QC reads), so we no longer
// extract the unmapped read stream here.
//
// Input BAM is the per-sample reference-merged alignment (ALIGNMENT.out.bams). We count
// primary mapped/unmapped (unmapped = reads that hit NONE of the pulled references = the
// closed-set residual), and (optionally) the Kraken2 classified total from the kreport, and
// surface them as `env` outputs so the subworkflow can fold them onto `meta` without any
// read-the-file-back gymnastics.
//
process EXTRACT_UNMAPPED {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::samtools=1.17' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(bam), path(csi), path(kreport)   // kreport may be NO_FILE

    output:
    tuple val(meta), env('TOTAL'), env('MAPPED'), env('UNMAPPED'), env('K2CLASS') , emit: counts
    path "versions.yml"                                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # --- read accounting (primary alignments only, no double counting) ---
    MAPPED=\$(samtools view -c -F 0x904 ${bam})           # primary, mapped
    UNMAPPED=\$(samtools view -c -f 0x4 -F 0x900 ${bam})  # primary, unmapped (closed-set residual)
    TOTAL=\$(( MAPPED + UNMAPPED ))

    # --- Kraken2 classified total from the root ("R") line, if a report was provided ---
    if [ "${kreport}" != "NO_FILE" ] && [ -s "${kreport}" ]; then
        K2CLASS=\$(awk -F'\\t' '\$4=="R"{gsub(/ /,"",\$2); print \$2; exit}' ${kreport})
        [ -z "\$K2CLASS" ] && K2CLASS=0
    else
        K2CLASS=0
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    TOTAL=0; MAPPED=0; UNMAPPED=0; K2CLASS=0
    echo '"${task.process}": {samtools: stub}' > versions.yml
    """
}
