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
    """
    N=\$(ls ref/ | wc -l)
    if [ "\$N" -eq 1 ]; then
        INPUT=\$(ls ref/*)
    else
        # Several references: concatenate (decompressing as needed) into one.
        for f in ref/*; do
            case "\$f" in
                *.gz) gzip -dc "\$f" ;;
                *)    cat "\$f" ;;
            esac
        done > combined_host.fa
        INPUT=combined_host.fa
    fi

    deacon index build \\
        -k ${params.deacon_kmer} \\
        -w ${params.deacon_window} \\
        -t ${task.cpus} \\
        ${args} \\
        "\$INPUT" > ${out}

    rm -f combined_host.fa
    """

    stub:
    """
    touch ${name}.k${params.deacon_kmer}w${params.deacon_window}.idx
    """
}
