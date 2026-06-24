// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # (header trimmed for the sketch -- copy the full APL header from the other modules)
// ##############################################################################################
//
// Pull the residual ("dark") reads + the read accounting the novelty score needs, in one place.
//
// Input BAM is the per-sample reference-merged alignment (ALIGNMENT.out.bams). Reads with the
// unmapped flag = reads that hit NONE of the pulled references = exactly the closed-set residual.
// We also count primary mapped/unmapped here, and (optionally) the Kraken2 classified total from
// the kreport, and surface them as `env` outputs so the subworkflow can fold them onto `meta`
// without any read-the-file-back gymnastics.
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
    tuple val(meta), path("*.unmapped.fastq.gz")                          , emit: reads
    tuple val(meta), env('TOTAL'), env('MAPPED'), env('K2CLASS')          , emit: counts
    path "versions.yml"                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # --- residual reads: anything that mapped to no reference ---
    # (translated search scores reads independently, so a single stream is fine for PE too)
    samtools fastq -f 0x4 -F 0x900 ${bam} 2>/dev/null | gzip -c > ${prefix}.unmapped.fastq.gz

    # --- read accounting (primary alignments only, no double counting) ---
    MAPPED=\$(samtools view -c -F 0x904 ${bam})        # primary, mapped
    UNMAPPED=\$(samtools view -c -f 0x4 -F 0x900 ${bam})  # primary, unmapped
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
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip -c > ${prefix}.unmapped.fastq.gz
    TOTAL=0; MAPPED=0; K2CLASS=0
    echo '"${task.process}": {samtools: stub}' > versions.yml
    """
}
