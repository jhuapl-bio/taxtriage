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
// Fallback path: fetch a run straight from NCBI with prefetch + fasterq-dump.
//
// Used only when ENA has no FASTQ files for the run (recently released runs, runs ENA
// has not mirrored) or when --sra_force_sratools is set. Slower than the ENA path
// because it downloads the .sra archive and then converts it, but it works for anything
// public in SRA.
//
// Paired detection: --split-3 writes _1/_2 for paired runs and a single file for
// single-end ones. The manifest's `single_end` flag comes from ENA/NCBI metadata, but
// the file layout produced here is what is actually emitted, so a run mislabelled
// upstream still yields the correct read files. --split-3 also drops orphan mates into
// a bare <run>.fastq; that file is discarded for paired runs so downstream tools never
// see a stray third FASTQ.
//
// Caching: `storeDir` -- same contract as SRA_FETCH_ENA. See --sra_cache_dir.
//
process SRA_FETCH_SRATOOLS {
    tag "$run_accession"
    label 'process_medium'

    storeDir { params.sra_cache_dir ? "${params.sra_cache_dir}/${run_accession}"
                                    : "${params.outdir}/sra_downloads/${run_accession}" }

    conda (params.enable_conda ? "bioconda::sra-tools=3.0.8 conda-forge::pigz=2.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:3.0.8--h9f5acd7_0' :
        'quay.io/biocontainers/sra-tools:3.0.8--h9f5acd7_0' }"

    errorStrategy { task.attempt <= 2 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(run_accession), val(single_end)

    output:
    tuple val(run_accession), path("${run_accession}*.fastq.gz"), emit: reads
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: '--split-3 --skip-technical'
    def max_size = params.sra_max_size ?: 'u'
    """
    # sra-tools refuses to run without a settings file, and \$HOME is not always
    # writable inside the container -- point NCBI_SETTINGS at the task work dir.
    # printf (not a heredoc) because the mkfg parser rejects leading whitespace.
    export NCBI_SETTINGS="\$PWD/user-settings.mkfg"
    printf '/LIBS/GUID = "%s"\\n/libs/cloud/report_instance_identity = "true"\\n' \\
        "\$(cat /proc/sys/kernel/random/uuid 2>/dev/null || echo taxtriage-sra)" \\
        > "\$NCBI_SETTINGS"

    prefetch \\
        --max-size ${max_size} \\
        --progress \\
        -O sra_tmp \\
        ${run_accession}

    fasterq-dump \\
        ${args} \\
        --threads ${task.cpus} \\
        --outdir . \\
        sra_tmp/${run_accession}

    # --split-3 emits <run>.fastq alongside _1/_2 when a paired run has orphan mates.
    # Keep only the proper pair in that case.
    if [ -f "${run_accession}_1.fastq" ] && [ -f "${run_accession}_2.fastq" ]; then
        rm -f "${run_accession}.fastq"
    fi

    if ! ls ${run_accession}*.fastq > /dev/null 2>&1; then
        echo "ERROR: fasterq-dump produced no FASTQ files for ${run_accession}" >&2
        exit 1
    fi

    if command -v pigz > /dev/null 2>&1; then
        pigz -p ${task.cpus} -f ${run_accession}*.fastq
    else
        gzip -f ${run_accession}*.fastq
    fi

    rm -rf sra_tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(fasterq-dump --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n1)
    END_VERSIONS
    """
}
