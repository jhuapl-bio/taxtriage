// ##############################################################################################
// # Copyright 2025 The Johns Hopkins University Applied Physics Laboratory LLC
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

/**
 * COMBINE_SAMPLES_JSON
 *
 * Merges all per-sample *.paths.json files produced by ALIGNMENT_PER_SAMPLE
 * into a single all.samples.json and publishes it to the report folder.
 *
 * The combined file uses the format:
 *   { "taxtriage_combined": true, "version": "1.0", "samples": [...] }
 *
 * Users can drag this single file onto any TaxTriage heatmap.html report
 * instead of importing individual per-sample JSONs one by one.
 */
process COMBINE_SAMPLES_JSON {
    tag "all_samples_combined"
    label 'process_low'
    publishDir "${params.outdir}/report", mode: 'copy'

    conda (params.enable_conda ? "conda-forge::python" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/reportlab-pdf:4.0.9' :
        'jhuaplbio/reportlab-pdf:4.0.9' }"

    input:
    // All per-sample .paths.json files (collected across all sample types)
    path(json_files)

    output:
    path "all.odr.json", emit: combined_json
    path "versions.yml",     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def json_inputs = json_files instanceof List
        ? json_files.findAll { it.name != 'NO_FILE' && it.name.endsWith('.json') }.join(' ')
        : (json_files.name != 'NO_FILE' && json_files.name.endsWith('.json') ? json_files.toString() : '')

    """
    combine_samples_json.py -i ${json_inputs} -o all.odr.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}
