process REMOTE_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::blast=2.12.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0' :
        'biocontainers/blast:2.12.0--pl5262h3289130_0' }"

    input:
    tuple val(meta), path(fasta)
    path  db

    output:
    tuple val(meta), path('*.blastn.txt'), emit: txt
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > ${prefix}.blastn.txt
    blastn \\
        \\
        -db $db -remote \\
        -query $fasta \\
        $args \\
        >>  ${prefix}.blastn.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
