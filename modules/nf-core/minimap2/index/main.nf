process MINIMAP2_INDEX {
    label 'process_high'
    tag "$meta.id"

    conda (params.enable_conda ? 'bioconda::minimap2=2.31' : null)
    // --use_alphaevolve: build the index with the same minimap2 that will align against it.
    container "${ params.use_alphaevolve ? params.alphaevolve_minimap2_container :
        workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.31--h118bc1c_0' :
        'biocontainers/minimap2:2.31--h118bc1c_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mmi"), emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // minimap2 reads gzip FASTA directly: no decompression. x.fa.gz -> x.mmi
    def stem = WorkflowTaxtriage.fastaStem(fasta.name)
    """
    minimap2 \\
        -t $task.cpus \\
        -d ${stem}.mmi \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
