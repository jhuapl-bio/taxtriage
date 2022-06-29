process PULL_FASTA {
    
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'pegi3s/biopython:latest' }"

    input:
        tuple val(meta), val(taxid_file), val(classified_reads), val(classified_reads_assignment), val(genomes)

    output:
        tuple val(meta), path("*filtered.fastq"), path("*refs.fasta"), optional: false, emit: fastq



    
    
    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    def outfile = "${meta.id}_refs.fasta"
    def filtered = classified_reads.findAll{ it =~ /.*\.classified.*(fq|fastq)(\.gz)?/  }
    def classifieds = filtered.join(" ")

    """
        
        getfasta_refs.py  \\
            -i "$taxid_file" \\
            -o $outfile \\
            -t file \\
            -r $genomes \\
            -d $classifieds \\
            -a $classified_reads_assignment \\
            -q 
        
    
    """
}

 