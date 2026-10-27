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

/**
 * MASK_REPORT
 *
 * Visual summary of a masked run (bin/mask_report.py): per-sample read
 * accounting, an organism table (masked vs unmasked TASS when --mask_compare is
 * on), coverage tracks with the masked / targeted spans drawn in, and the NCBI
 * FCS assembly table. Published to <outdir>/masking/.
 */
process MASK_REPORT {
    tag "all_samples"
    label 'process_low'
    publishDir "${params.outdir}/masking", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.equals('versions.yml') ? null : fn }

    conda (params.enable_conda ? "conda-forge::python" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.2' :
        'jhuaplbio/taxtriage_confidence:2.2' }"

    input:
    path(masked_jsons,   stageAs: 'masked/*')
    path(baseline_jsons, stageAs: 'unmasked/*')
    path(mask_files,     stageAs: 'mask/*')

    output:
    path "all.masking.html" , emit: html
    path "all.masking.tsv"  , emit: tsv
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // stageAs into sub-directories makes .name the staged relative path
    // (e.g. 'beds/1/x.bed'), so compare on the basename
    def base = { f -> f.toString().tokenize('/')[-1] }
    def real = { files -> (files instanceof List ? files : [files]).findAll { it && !base.call(it).startsWith('NO_FILE') } }
    def mj   = real.call(masked_jsons).join(' ')
    def bj   = real.call(baseline_jsons)
    def mf   = real.call(mask_files)
    def asm  = mf.find { base.call(it) == 'fcs_assemblies.tsv' }
    def pref = ((params.mask_mode ?: 'exclude') == 'target' || (params.mask_products && params.mask_products.toString() != 'false')) ? 'target' : 'exclude'
    def summ = mf.find { base.call(it) == "mask.${pref}.summary.tsv" } ?: mf.find { base.call(it).endsWith('.summary.tsv') }
    def bjs  = mf.find { base.call(it) == "mask.${pref}.build.json" } ?: mf.find { base.call(it).endsWith('.build.json') }
    def top  = params.mask_report_top ?: 24
    """
    mask_report.py --jsons ${mj} \\
        ${ bj ? '--baseline ' + bj.join(' ') : '' } \\
        ${ asm ? '--assembly_table ' + asm : '' } \\
        ${ summ ? '--mask_summary ' + summ : '' } \\
        ${ bjs ? '--mask_json ' + bjs : '' } \\
        --top ${top} \\
        --prefix all.masking

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}
