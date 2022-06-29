process GET_ASSEMBLIES {
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:4.2.0' :
        'cicirello/gnu-on-alpine' }"

    


    input:
    tuple val(meta), val(inputs)


    output:
    tuple val(meta), path("assembly_summary_refseq.txt"), optional: true, emit: assembly
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when




     

    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/


    """
    if [[ ! -s 'assembly_summary_refseq.txt' ]] ; then
        echo "Downloading the assembly summary file from ncbi...."
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -O assembly_summary_refseq.txt
    else 
        echo "Assembly file exists"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(wget --version 2>&1)
    END_VERSIONS

    """
}
