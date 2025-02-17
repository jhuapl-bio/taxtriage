process KRAKEN2_KRAKEN2 {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::kraken2=2.1.2 conda-forge::pigz=2.6' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' :
        'biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' }"

    input:
    tuple val(meta), path(reads)
    path  db
    val save_output_fastqs
    val save_reads_assignment

    output:
    tuple val(meta), path('*classified*')     , optional:true, emit: classified_reads_fastq
    tuple val(meta), path('*unclassified*')   , optional:true, emit: unclassified_reads_fastq
    tuple val(meta), path('*classifiedreads*'), optional:true, emit: classified_reads_assignment
    tuple val(meta), path('*report.txt')                     , emit: report
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired       = meta.single_end ? "" : "--paired"
    // get basename of the db 
    def db_basename = db.getName()

    def classified   = meta.single_end ? "${prefix}.${db_basename}.classified.fastq"   : "${prefix}.${db_basename}.classified#.fastq"
    def unclassified = meta.single_end ? "${prefix}.${db_basename}.unclassified.fastq" : "${prefix}.${db_basename}.unclassified#.fastq"
    def classified_command = save_output_fastqs ? "--classified-out ${classified}" : ""
    def confidence = params.k2_confidence ? "--confidence ${params.k2_confidence}" : ""
    def unclassified_command = save_output_fastqs ? "--unclassified-out ${unclassified}" : ""
    def readclassification_command = save_reads_assignment ? "--output ${prefix}.kraken2.classifiedreads.txt" : ""
    def compress_reads_command = save_output_fastqs ? "pigz -p $task.cpus *.fastq" : ""
    def minimum_hit_groups = params.k2_minimum_hit_groups ? "--minimum-hit-groups ${params.k2_minimum_hit_groups}" : ""
    
    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.kraken2.report.txt \\
        --gzip-compressed \\
        $unclassified_command \\
        $classified_command \\
        $minimum_hit_groups \\
        $readclassification_command \\
        $paired $confidence \\
        $args  \\
        $reads > /dev/null 
        
    $compress_reads_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}