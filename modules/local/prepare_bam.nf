// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # All rights reserved.
// # Permission is hereby granted, free of charge, to any person obtaining a copy of this
// # software and associated documentation files (the "Software"), to deal in the Software
// # without restriction, including without limitation the rights to use, copy, modify,
// # merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// # permit persons to whom the Software is furnished to do so.
// #
// # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// # INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// # PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// # LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// # TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
// # OR OTHER DEALINGS IN THE SOFTWARE.
// #

//
// Normalise a user-supplied alignment file (BAM / CRAM / SAM) so it matches what
// the rest of the pipeline expects from ALIGNMENT:
//   * converted to BAM (CRAM/SAM inputs are decoded; CRAM needs --reference)
//   * coordinate sorted (skipped when the header already says SO:coordinate)
//   * CSI indexed (same index flavour SAMTOOLS_INDEX produces with `-c`)
//   * optionally MAPQ filtered, matching the aligner-side --minmapq behaviour
// Also reports the number of primary reads, which stands in for COUNT_READS since
// no FASTQ ever exists for these samples.
//
process PREPARE_BAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(alignment)
    path(reference)

    output:
    tuple val(meta), path("*.prepared.bam"), path("*.prepared.bam.csi"), emit: bam
    tuple val(meta), path("*.readcount.txt")                           , emit: count
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    // Only filter when the user actually asked for it; params.minmapq defaults to 5
    // for the aligner but a pre-aligned BAM should be taken as given unless
    // --bam_minmapq is set explicitly.
    def mapq    = params.bam_minmapq ? " -q ${params.bam_minmapq} " : ""
    def ref_arg = (reference && reference.name != 'NO_FILE') ? " --reference ${reference} " : ""

    """
    set -euo pipefail

    IN="${alignment}"

    # 1. Normalise the container to BAM.  CRAM is decoded here (needs the
    #    reference); SAM is compressed.  A BAM input is left untouched.
    case "\${IN}" in
        *.cram|*.CRAM|*.sam|*.SAM)
            samtools view ${ref_arg} -@ ${task.cpus} -b -h -o converted.bam "\${IN}"
            IN=converted.bam
            ;;
    esac

    # 2. Coordinate sort unless the header already says so.  Indexing and
    #    samtools coverage both require coordinate order.
    sort_order=\$(samtools view -H "\${IN}" | grep -m1 '^@HD' | tr '\\t' '\\n' | grep '^SO:' | cut -d: -f2 || true)
    if [ "\${sort_order:-unknown}" != "coordinate" ]; then
        echo "Sort order of ${alignment} is '\${sort_order:-unknown}'; coordinate sorting."
        samtools sort -@ ${task.cpus} $args -o sorted.bam "\${IN}"
        IN=sorted.bam
    fi

    # 3. Optional MAPQ filter; otherwise link the normalised file through
    #    without copying it again.
    if [ -n "${mapq}" ]; then
        samtools view -@ ${task.cpus} ${mapq} -b -h -o ${prefix}.prepared.bam "\${IN}"
    elif [ "\${IN}" = "${alignment}" ]; then
        ln -s "\$(readlink -f "\${IN}")" ${prefix}.prepared.bam
    else
        mv "\${IN}" ${prefix}.prepared.bam
    fi

    samtools index -c -@ ${task.cpus} ${prefix}.prepared.bam

    # Primary, non-supplementary reads — the closest analogue to the FASTQ read
    # count that COUNT_READS provides for read-based samples.
    samtools view -c -F 0x900 ${prefix}.prepared.bam > ${prefix}.readcount.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.prepared.bam
    touch ${prefix}.prepared.bam.csi
    echo 0 > ${prefix}.readcount.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}
