process DOWNLOAD_ASSEMBLY {
    label 'process_medium'


    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'pegi3s/biopython:latest' }"

    


    input:
    tuple val(meta), val(taxid_containing_file), val(classified_reads_fastq), val(classified_reads)
    val(assembly)


    output:
    tuple val(meta), val(taxid_containing_file), val(classified_reads_fastq), val(classified_reads), path("*.output.references.fasta"), optional: false, emit: fasta
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when




     

    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    def email = params.email ? params.email : "brian.merritt@jhuapl.edu"
    // def taxid_join = taxids.join(" ")


    """
    
    

    download_fastas.py \\
            -i "${taxid_containing_file}" \\
            -o ${meta.id}.output.references.fasta  \\
            -e ${email} -f file \\
            -t ${assembly} -k 


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS

    """
}