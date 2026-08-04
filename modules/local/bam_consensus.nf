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
// Reconstruct reference sequence from a pre-aligned BAM when the user did not
// supply --reference_fasta.
//
// match_paths.py can already fall back to the BAM header for reference LENGTHS,
// but every sketching step (sourmash region signatures, the shared-window
// conflict report, and the taxid-level ANI matrix) needs actual bases.  This
// module calls a per-reference consensus straight from the alignment so those
// steps have sequence to work with.
//
// CAVEATS (documented in the README as well):
//   * only positions covered by reads are recovered; everything else is N, so
//     ANI/containment is computed over the covered fraction of each reference.
//   * multi-mapping reads contribute to the consensus of EVERY reference they
//     were placed on, which inflates the apparent similarity between related
//     organisms and therefore makes conflict-driven read removal more
//     aggressive than it would be against true reference genomes.
// Supplying --reference_fasta is always preferable when it is available.
//
process BAM_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::samtools=1.19.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(csi)

    output:
    tuple val(meta), path("*.consensus.fasta"), emit: fasta
    tuple val(meta), path("*.consensus_stats.tsv"), emit: stats
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def min_depth = params.consensus_min_depth ?: 1
    def min_bases = params.consensus_min_bases ?: 500
    def min_mq    = params.consensus_min_mapq != null ? " --min-MQ ${params.consensus_min_mapq} " : ""
    def mode      = params.consensus_mode ?: 'simple'

    """
    set -euo pipefail

    # `-a` keeps every reference at full length with uncovered positions as N so
    # consensus coordinates stay identical to the original reference — the
    # shared-window sketching in conflict_regions.py slides fixed windows along
    # these sequences and would otherwise compare mismatched coordinates.
    samtools consensus \\
        -@ ${task.cpus} \\
        -m ${mode} \\
        -a \\
        -d ${min_depth} \\
        ${min_mq} \\
        -f fasta \\
        --show-ins no \\
        --show-del yes \\
        $args \\
        ${bam} \\
        > raw_consensus.fasta

    # Drop references whose consensus is (almost) all N — i.e. references present
    # in the BAM header that no read actually supports.  Keeping them would add
    # empty sketches and drag every ANI comparison toward zero.
    awk -v minb=${min_bases} '
        function flush_rec() {
            if (hdr != "" && nonN >= minb) { print hdr; print seq }
            if (hdr != "") { printf "%s\\t%d\\t%d\\n", substr(hdr, 2), len, nonN >> statsf }
        }
        /^>/ {
            flush_rec()
            hdr = \$0; seq = ""; nonN = 0; len = 0; next
        }
        {
            seq = seq \$0
            len += length(\$0)
            t = \$0
            gsub(/[Nn]/, "", t)
            nonN += length(t)
        }
        END { flush_rec() }
    ' statsf="${prefix}.consensus_stats.tsv" raw_consensus.fasta > ${prefix}.consensus.fasta

    # Always emit the files even when nothing passed the filter, so downstream
    # channel joins do not silently drop the sample.
    touch ${prefix}.consensus.fasta ${prefix}.consensus_stats.tsv

    n_refs=\$(grep -c '^>' ${prefix}.consensus.fasta || true)
    echo "Recovered consensus for \${n_refs} reference(s) from ${bam} (min depth ${min_depth}, min ${min_bases} non-N bases)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.consensus.fasta
    touch ${prefix}.consensus_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}
