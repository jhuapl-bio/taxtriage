// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # All rights reserved.
// # Permission is hereby granted, free of charge, to any person obtaining a copy of this
// # software and associated documentation files (the "Software"), to deal in the Software
// # without restriction, including without limitation the rights to use, copy, modify,
// # merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// # permit persons to whom the Software is furnished to do so.
// #
// # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// # INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// # PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// # LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// # TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
// # OR OTHER DEALINGS IN THE SOFTWARE.
// #
process GET_ASSEMBLIES {
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget%3A1.18--h7132678_6' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"


    output:
        path("assembly_summary_refseq.txt"), optional: false, emit: assembly
        path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when






    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/


    """
    if [[ ! -s 'assembly_summary_refseq.txt' ]] ; then
        echo "Downloading the RefSeq assembly summary file from ncbi...."
        wget --no-check-certificate https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt  -O assembly_summary_refseq.txt
    else
        echo "RefSeq assembly summary file exists"
    fi

    if [[ "${params.enable_genbank}" == "true" ]] ; then
        if [[ ! -s 'assembly_summary_genbank.txt' ]] ; then
            echo "Downloading the GenBank assembly summary file from ncbi...."
            wget --no-check-certificate https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt -O assembly_summary_genbank.txt
        else
            echo "GenBank assembly summary file exists"
        fi
    else
        echo "GenBank pulling disabled (enable with --enable_genbank)"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version 2>&1)
    END_VERSIONS

    """
}
