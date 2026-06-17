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
// ##############################################################################################
//
// Predict open reading frames from assembled contigs using Pyrodigal in metagenome mode.
// Protein sequences (.faa) feed directly into MMSEQS_LINCLUST → MMSEQS_TAXONOMY for the
// novelty-detection branch. Using predicted proteins instead of raw reads avoids the internal
// 6-frame translation that makes mmseqs taxonomy so slow on read-level input.
//
// Pyrodigal docs: https://pyrodigal.readthedocs.io
//

process PYRODIGAL {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pyrodigal=3.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyrodigal:3.7.1--py313hd72fa03_0' :
        'quay.io/biocontainers/pyrodigal:3.7.1--py313hd72fa03_0' }"

    input:
    tuple val(meta), path(fasta)   // assembled contigs (fasta / fasta.gz)

    output:
    tuple val(meta), path("*.faa"), emit: proteins   // predicted protein sequences
    tuple val(meta), path("*.fna"), emit: genes       // predicted nucleotide CDS sequences
    tuple val(meta), path("*.gff"), emit: gff         // GFF3 annotation
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '-p meta'    // -p meta = metagenome mode (no training)
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Decompress transparently if gzipped; pyrodigal reads from stdin via -i
    def input_cmd = fasta.name.endsWith('.gz') ? "zcat ${fasta}" : "cat ${fasta}"
    """
    ${input_cmd} | pyrodigal \\
        -i /dev/stdin \\
        -a ${prefix}.faa \\
        -d ${prefix}.fna \\
        -f gff \\
        -o ${prefix}.gff \\
        ${args} \\
        -j ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyrodigal: \$(pyrodigal --version 2>&1 | sed 's/pyrodigal //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.faa ${prefix}.fna ${prefix}.gff
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyrodigal: stub
    END_VERSIONS
    """
}
