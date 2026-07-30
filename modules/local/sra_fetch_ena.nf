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
//
// Fast path: pull the ready-made gzipped FASTQs that ENA hosts for a run.
//
// This is preferred over sra-tools because ENA serves the files already split into
// _1/_2 and already gzipped -- there is no .sra download and no fasterq-dump conversion,
// which is where most of the wall-clock time (and disk) goes on the NCBI path.
//
// Caching: `storeDir` writes the finished FASTQs to a persistent per-run folder and
// SKIPS the task entirely on any later run where those files already exist. That covers
// `-resume` within a run and fresh runs across sessions, so a given accession is
// downloaded exactly once. Point --sra_cache_dir at shared storage to share downloads
// between users/projects.
//
process SRA_FETCH_ENA {
    tag "$run_accession"
    label 'process_low'

    storeDir { params.sra_cache_dir ? "${params.sra_cache_dir}/${run_accession}"
                                    : "${params.outdir}/sra_downloads/${run_accession}" }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.8.3' }"

    errorStrategy { task.attempt <= 2 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(run_accession), val(single_end), val(url1), val(url2), val(md5_1), val(md5_2)

    output:
    tuple val(run_accession), path("${run_accession}*.fastq.gz"), emit: reads
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // A run is paired only when the manifest says so AND ENA actually listed a
    // second file; otherwise the single file is written as <run>.fastq.gz.
    def is_paired = single_end.toString() != 'true' && url2
    def out1 = is_paired ? "${run_accession}_1.fastq.gz" : "${run_accession}.fastq.gz"
    def sum1 = md5_1 ? "--md5 ${md5_1}" : ''
    def read2 = ''
    if (is_paired) {
        read2 = "--url2 '${url2}' --out2 ${run_accession}_2.fastq.gz"
        if (md5_2) {
            read2 += " --md52 ${md5_2}"
        }
    }
    """
    fetch_ena_fastq.py \\
        --url '${url1}' \\
        --out ${out1} \\
        ${sum1} \\
        ${read2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
