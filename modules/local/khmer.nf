process KHMER_NORMALIZEBYMEDIAN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/khmer:3.0.0a3--py37haa7609a_2' :
        'biocontainers/khmer:3.0.0a3--py37haa7609a_2' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("${meta.id}.normalized_*.fastq.gz"), emit: reads
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // calculate the amount of memory for khmer
    def m = task.memory * 0.8
    // build two alternate commands
    def seCmd = """\
        normalize-by-median.py \\
            -M ${m.toMega()}e6 \\
            --gzip \\
            --force-single \\
            ${files} \\
            -o ${meta.id}.normalized.fastq.gz
        """.stripIndent()

    def peCmd = """\
        normalize-by-median.py \\
            -M ${m.toMega()}e6 \\
            --gzip \\
            --paired \\
            ${files} \\
            -o - \\
        | split-paired-reads.py --gzip \\
            -1 ${meta.id}.normalized_R1.fastq.gz \\
            -2 ${meta.id}.normalized_R2.fastq.gz \\
            -
        """.stripIndent()

    // choose the right one
    if ( meta.single_end ) {
        seCmd
    } else {
        peCmd
    }
    
    // and then write out your versions.yml as before
    """
    ${ meta.single_end ? seCmd : peCmd }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer: \$( normalize-by-median.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """
}
