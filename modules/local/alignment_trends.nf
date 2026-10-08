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
 * ALIGNMENT_TRENDS
 *
 * Cross-sample analysis of alignment depth (on by default; --alignment_trends false
 * turns it off). Lines up every
 * sample's per-reference depth_profile (written by match_paths.py) on a common
 * window grid and reports how often each window is zero-depth, low-depth or
 * high-depth across the samples, plus the merged recurrent regions. Recurrent
 * zero / low regions mark stretches of the reference assemblies this sample set
 * does not carry (deletions, divergent or novel loci, assembly artefacts);
 * recurrent high regions mark repeats, rRNA operons, mobile elements and the like.
 *
 * Published to <outdir>/alignment_trends/, and its JSON is embedded in the HTML
 * report (Trends > Alignment Trends > Pipeline results) next to a live version.
 */
process ALIGNMENT_TRENDS {
    tag "all_samples"
    label 'process_low'
    publishDir "${params.outdir}/alignment_trends", mode: 'copy'

    conda (params.enable_conda ? "bioconda::pysam conda-forge::matplotlib conda-forge::openpyxl" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.2' :
        'jhuaplbio/taxtriage_confidence:2.2' }"

    input:
    path(json_files)

    output:
    path "all.alignment_trends.*.tsv" , emit: tsv
    path "all.alignment_trends.json"  , emit: json
    path "all.alignment_trends.xlsx"  , optional: true, emit: xlsx
    path "all.alignment_trends.plots" , optional: true, emit: plots
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def jsons = (json_files instanceof List ? json_files : [json_files])
        .findAll { !it.name.startsWith('NO_FILE') && it.name.endsWith('.json') }
        .join(' ')
    def low_abs  = params.trend_low_abs  != null ? " --low_abs ${params.trend_low_abs} "   : ' '
    def high_abs = params.trend_high_abs != null ? " --high_abs ${params.trend_high_abs} " : ' '
    def plots    = params.trend_plots           ? " --plots ${params.trend_plots} "        : ' ' 
    def matrix   = params.trend_matrix ? ' --matrix ' : ' '
    """
    alignment_trends.py \\
        -i ${jsons} \\
        -o all.alignment_trends \\
        --min_samples ${params.trend_min_samples} \\
        --min_reads ${params.trend_min_reads} \\
        --min_reads_per_window ${params.trend_min_reads_per_window} \\
        --low_frac ${params.trend_low_frac} \\
        --high_frac ${params.trend_high_frac} \\
        --min_freq ${params.trend_min_freq} \\
        --min_region_windows ${params.trend_min_region_windows} \\
        ${low_abs} ${high_abs} ${plots} ${matrix} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}
