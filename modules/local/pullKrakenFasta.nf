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
process PULL_FASTA {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.75' }"

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

