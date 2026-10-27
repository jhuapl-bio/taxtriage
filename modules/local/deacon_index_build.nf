// Build a Deacon minimizer index from one or more host FASTA files.
//
// Inputs can be plain or gzipped FASTA (a local --remove_reference_file list,
// an iGenomes `fasta`, or the genomes FETCH_HOST_REFS pulled for a named host
// target).  Multiple files are concatenated into a single index.  The index is
// cached through `storeDir`; the file name carries k and w so changing either
// triggers a rebuild.
process DEACON_INDEX_BUILD {
    tag "$name"
    label 'process_high'
    storeDir params.deacon_index_dir ?: "${params.host_reference_dir ?: "${params.outdir}/host_references"}/deacon"

    conda "bioconda::deacon=0.18.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deacon:0.18.0--h3edb6b3_0' :
        'biocontainers/deacon:0.18.0--h3edb6b3_0' }"

    input:
    tuple val(name), path(fastas, stageAs: 'ref/*')

    output:
    path("${name}.k${params.deacon_kmer}w${params.deacon_window}.idx"), emit: index

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def out  = "${name}.k${params.deacon_kmer}w${params.deacon_window}.idx"
    // deacon reads gzip/zstd/xz FASTA itself, so nothing is decompressed to disk.
    // Several references are streamed (decompressed on the fly) into stdin ('-').
    """
    set -o pipefail
    N=\$(ls ref/ | wc -l)
    if [ "\$N" -eq 1 ]; then
        deacon index build \\
            -k ${params.deacon_kmer} \\
            -w ${params.deacon_window} \\
            -t ${task.cpus} \\
            ${args} \\
            "\$(ls ref/*)" > ${out}
    else
        for f in ref/*; do
            if [ "\$(head -c 2 "\$f" | od -An -tx1 | tr -d ' \\n')" = "1f8b" ]; then
                gzip -dc "\$f"
            else
                cat "\$f"
            fi
        done | deacon index build \\
            -k ${params.deacon_kmer} \\
            -w ${params.deacon_window} \\
            -t ${task.cpus} \\
            ${args} \\
            - > ${out}
    fi
    """

    stub:
    """
    touch ${name}.k${params.deacon_kmer}w${params.deacon_window}.idx
    """
}
