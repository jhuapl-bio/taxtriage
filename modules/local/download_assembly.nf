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
process DOWNLOAD_ASSEMBLY {
    label 'process_medium'
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://pegi3s/biopython:latest' :
        'biocontainers/biopython:1.75' }"




    input:
    tuple val(meta), path(hits_containing_file)
    path assembly


    output:
    tuple val(meta),  path("*.dwnld.references.fasta"), optional: true, emit: fasta
    tuple val(meta),  path("*.dwnld.gcfids.txt"), optional: false, emit: gcfids
    tuple val(meta),  path("*.dwnld.gcfmapping.tsv"), optional: false, emit: mapfile
    tuple val(meta),  path("missing.txt"), optional: true, emit: missings

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when






    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    def email = params.email ? " -e ${params.email}" : ""
    def column = " -c 1 "
    def columnAssembly = params.fuzzy ? " -a 7 "  : " -a 5,6 "
    def matchcol = params.fuzzy ? " -a 7 "  : " -a 5,6 "
    def refresh_download = params.refresh_download ? " -r " : ""
    def type = hits_containing_file ? " -f file " : " -f list  "


    """



    download_fastas.py \\
            -i  ${hits_containing_file} \\
            -o ${meta.id}.dwnld.references.fasta ${refresh_download} -m ${meta.id}.missing.txt \\
            ${email} $type -g ${meta.id}.dwnld.gcfmapping.tsv \\
            -t ${assembly} -k  $column $columnAssembly -y 7 -r

    cut -f 2 ${meta.id}.dwnld.gcfmapping.tsv  | sort | uniq  > ${meta.id}.dwnld.gcfids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version )
    END_VERSIONS

    """
}
