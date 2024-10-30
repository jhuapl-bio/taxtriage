process MEGAHIT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0' :
        'biocontainers/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("megahit_out/*.contigs.fa.gz")                            , emit: contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.contigs.fa.gz")      , emit: k_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.addi.fa.gz")         , emit: addi_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.local.fa.gz")        , emit: local_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.final.contigs.fa.gz"), emit: kfinal_contigs
    path "versions.yml"                                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def decompress_files = ""
    // remove .gz from end of reads
    def remove_files = ""
    if (params.decompress_pre_megahit){
        if (!meta.single_end) {
            decompress_files = "gzip -dfc -k ${reads[0]} > 1.fastq ; gzip -dfc -k ${reads[1]} > 2.fastq"
            reads[0] = "1.fastq"
            reads[1] = "2.fastq"
            reads[0] = reads[0].toString().replace(".gz", "")
            reads[1] = reads[1].toString().replace(".gz", "")
            remove_files = "rm 1.fastq 2.fastq"
        } else {
            decompress_files = "gzip -dfc -k ${reads} > 1.fastq"
            reads = reads.toString().replace(".gz", "")
            remove_files = "rm 1.fastq"
            reads = "1.fastq"
        }
    }
    
    if (meta.single_end) {
        """
        ${decompress_files}
        megahit \\
            -r $reads \\
            -t $task.cpus \\
            $args \\
            --out-prefix $prefix; ${remove_files}

        pigz \\
            --no-name \\
            -p $task.cpus \\
            $args2 \\
            megahit_out/*.fa \\
            megahit_out/intermediate_contigs/*.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    } else {
        """
        ${decompress_files}
        megahit \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -t $task.cpus \\
            $args \\
            --out-prefix $prefix; ${remove_files}

        pigz \\
            --no-name \\
            -p $task.cpus \\
            $args2 \\
            megahit_out/*.fa \\
            megahit_out/intermediate_contigs/*.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    }
    
}
