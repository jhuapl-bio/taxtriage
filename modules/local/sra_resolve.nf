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
// Turn user-supplied SRA/ENA accessions into a concrete per-run download manifest.
//
// Runs ONCE for the whole pipeline (all accessions are batched into a single request
// file) so ENA/NCBI see one query per accession, not one per retry of every sample.
// The manifest it emits is what drives SRA_FETCH_ENA / SRA_FETCH_SRATOOLS: it carries
// the resolved run accessions, the FASTQ URLs, the md5s, the instrument platform and —
// critically — whether each run is paired-end, which is detected from ENA's file
// listing rather than guessed from the accession.
//
process SRA_RESOLVE {
    tag "resolve_sra_accessions"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.8.3' }"

    // Network access is required; retry a couple of times before failing the run.
    errorStrategy { task.attempt <= 2 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    path(requests)      // CSV: sample,accession

    output:
    path("sra_manifest.csv"), emit: manifest
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def force  = params.sra_force_sratools ? '--force-sratools' : ''
    def apikey = params.ncbi_api_key ? "--api-key ${params.ncbi_api_key}" : ''
    """
    resolve_sra.py \\
        -i ${requests} \\
        -o sra_manifest.csv \\
        ${force} \\
        ${apikey}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
