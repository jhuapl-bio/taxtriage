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

/**
 * DECOMPRESS_FASTA
 *
 * gzip / bgzip FASTA -> plain FASTA, for consumers that cannot read compressed
 * references (bowtie2-build, hisat2-build, bcftools, bedtools, CRAM decoding).
 *
 * The output name is the input minus its compression suffix, matched
 * case-insensitively (x.fna.gz -> x.fna, x.FASTA.GZ -> x.FASTA, x.fa.bgz -> x.fa),
 * and is guaranteed to differ from the input name: the stock GUNZIP module
 * strips only a literal '.gz', so for '.bgz' / '.GZ' it would write its output
 * through the staged input symlink and overwrite the user's original file.
 * Not published: the plain copy only lives in the work directory.
 */
process DECOMPRESS_FASTA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(archive, stageAs: 'compressed/*')

    output:
    tuple val(meta), path("${outname}"), emit: fasta
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def base = archive.toString().tokenize('/')[-1]
    outname  = base.replaceFirst(/(?i)\.(gz|bgz)$/, '')
    if (outname == base) { outname = base + '.fasta' }
    """
    gzip -cd "${archive}" > "${outname}"
    if ! grep -q -m1 '^>' "${outname}"; then
        echo "ERROR: ${base} does not decompress to a FASTA (no '>' header line found)" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
    END_VERSIONS
    """
}
