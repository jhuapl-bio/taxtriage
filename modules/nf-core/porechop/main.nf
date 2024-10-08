process PORECHOP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::porechop=0.2.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39h7cff6ad_2' :
        'biocontainers/porechop:0.2.4--py39h7cff6ad_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
   
    porechop \\
        -i $reads \\
        -t $task.cpus \\
        $args \\
        -o ${prefix}.tmp.fastq.gz

    if [ -s ${prefix}.tmp.fastq.gz ]
    then
        mv ${prefix}.tmp.fastq.gz ${prefix}.fastq.gz
    else
        rm ${prefix}.tmp.fastq.gz
        cp $reads ${prefix}.fastq.gz
    fi
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop: \$( porechop --version )
    END_VERSIONS
    """
}
