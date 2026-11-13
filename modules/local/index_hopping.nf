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
 * INDEX_HOPPING
 *
 * Flags detections that Illumina index hopping from another sample of the same
 * lane could explain (bin/index_hopping.py). Runs once over every real and
 * control sample's *.paths.json, after scoring and before any report is built.
 * Flag only: each JSON comes back with an `index_hop` block on every organism
 * node and `metadata.index_hopping`; TASS and calls are untouched.
 *
 * Inputs
 *   json_files      every per-sample *.paths.json (staged under in/)
 *   hop_files       the distinct index-hopping tables, staged as hop01/, hop02/...
 *                   in list order, so sample_map_json can point at them by index
 *   sample_map_json {"samples": {id: {"files": [index], "index_sample_id": str|null}}}
 */
process INDEX_HOPPING {
    tag "index_hopping"
    label 'process_low'
    publishDir "${params.outdir}/index_hopping", mode: 'copy',
        saveAs: { fn -> fn == 'versions.yml' ? null : (fn.startsWith('annotated/') ? "json/${fn.tokenize('/')[-1]}" : fn) }

    conda (params.enable_conda ? "conda-forge::python" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/reportlab-pdf:4.0.9' :
        'jhuaplbio/reportlab-pdf:4.0.9' }"

    input:
    path(json_files, stageAs: 'in/*')
    path(hop_files, stageAs: 'hop??/*')
    val(sample_map_json)

    output:
    path "annotated/*.json",   emit: jsons
    path "index_hopping.tsv",  emit: tsv
    path "index_hopping.json", emit: json
    path "versions.yml",       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def jsons = (json_files instanceof List ? json_files : [json_files]).join(' ')
    def hops  = (hop_files instanceof List ? hop_files : [hop_files]).join(' ')
    """
    cat > sample_map.json <<'END_MAP'
    ${sample_map_json}
    END_MAP

    index_hopping.py \\
        -i ${jsons} \\
        --hopping ${hops} \\
        --sample-map sample_map.json \\
        --rate ${params.index_hopping_rate} \\
        --factor ${params.index_hopping_factor} \\
        --alpha ${params.index_hopping_alpha} \\
        -o annotated \\
        --summary index_hopping.tsv \\
        --report index_hopping.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}
