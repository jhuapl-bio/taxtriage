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
 * BUILD_MASK
 *
 * One region mask for the whole run, restricted to the references the samples
 * actually aligned to (read from the BAM headers), so the multi-GB NCBI FCS
 * details files are streamed once per run rather than once per sample.
 *
 * Emits (all optional; ALIGNMENT_PER_SAMPLE picks them up by file name):
 *   mask.exclude.bed     never-scored regions: FCS spans (+ --mask_bed in exclude mode)
 *   mask.target.bed      --mask_mode target: --mask_bed + annotation-selected features
 *   fcs_assemblies.tsv   NCBI FCS per-assembly summary for the run's assemblies
 *   mask.*.summary.tsv   per-reference mask contents
 *   mask.*.build.json    build statistics
 *   fcs_subset.txt       the matching FCS rows; pass back as --fcs_details to skip the big scan
 */
process BUILD_MASK {
    tag "region_mask"
    label 'process_medium'
    publishDir "${params.outdir}/masking", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.equals('versions.yml') ? null : fn }

    conda (params.enable_conda ? "bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.2' :
        'jhuaplbio/taxtriage_confidence:2.2' }"

    input:
    path(bams,          stageAs: 'bams/?/*')
    path(mappings,      stageAs: 'maps/?/*')
    path(fcs_details,   stageAs: 'fcs_details/?/*')
    path(fcs_summary,   stageAs: 'fcs_summary/?/*')
    path(fcs_gx,        stageAs: 'fcs_gx/?/*')
    path(fcs_adaptor,   stageAs: 'fcs_adaptor/?/*')
    path(user_beds,     stageAs: 'beds/?/*')
    path(feature_files, stageAs: 'features/?/*')

    output:
    path "mask.*.bed"         , optional: true, emit: beds
    path "fcs_assemblies.tsv" , optional: true, emit: assemblies
    path "*.summary.tsv"      , optional: true, emit: summary
    path "*.build.json"       , optional: true, emit: json
    path "fcs_subset.txt"     , optional: true, emit: subset
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // stageAs into sub-directories makes .name the staged relative path
    // (e.g. 'beds/1/x.bed'), so compare on the basename
    def base = { f -> f.toString().tokenize('/')[-1] }
    def real = { files -> (files instanceof List ? files : [files]).findAll { it && !base.call(it).startsWith('NO_FILE') } }
    def bam_list  = real.call(bams).join(' ')
    def map_list  = real.call(mappings).join(' ')
    def det_list  = real.call(fcs_details).join(' ')
    def sum_list  = real.call(fcs_summary).join(' ')
    def gx_list   = real.call(fcs_gx).join(' ')
    def ad_list   = real.call(fcs_adaptor).join(' ')
    def bed_list  = real.call(user_beds).join(' ')
    def feats     = real.call(feature_files)
    def gffs      = feats.findAll { base.call(it) ==~ /.*\.gff3?(\.gz)?$/ }.join(' ')
    def ftabs     = feats.findAll { !(base.call(it) ==~ /.*\.gff3?(\.gz)?$/) }.join(' ')
    def refs_arg  = bam_list ? "--refs ${bam_list}" : ''
    def map_arg   = map_list ? "--mapping ${map_list}" : ''
    def actions   = params.fcs_actions ?: 'EXCLUDE,TRIM,FIX'
    def div_arg   = params.fcs_divisions ? "--divisions '${params.fcs_divisions}'" : ''
    def minlen    = params.fcs_min_len ? "--min_contam_len ${params.fcs_min_len}" : ''
    def target    = (params.mask_mode ?: 'exclude') == 'target'
    def has_fcs   = det_list || gx_list || ad_list
    def sel_arg   = params.mask_feature_select ? "--feature_select '${params.mask_feature_select}'" : ''
    def ftype_arg = params.mask_feature_types ? "--feature_types '${params.mask_feature_types}'" : ''
    def pad_arg   = params.mask_pad ? "--pad ${params.mask_pad}" : ''

    // exclude BED: FCS spans, plus the user BEDs when the mode is exclude
    def excl_srcs = [
        det_list ? "--fcs_details ${det_list}" : '',
        gx_list  ? "--fcs_gx_report ${gx_list}" : '',
        ad_list  ? "--fcs_adaptor_report ${ad_list}" : '',
        (!target && bed_list) ? "--bed ${bed_list}" : '',
    ].findAll { it }.join(' ')
    // target BED: user BEDs + annotation features selected by --mask_feature_select
    def tgt_srcs = [
        bed_list ? "--bed ${bed_list}" : '',
        ftabs    ? "--feature_table ${ftabs}" : '',
        gffs     ? "--gff ${gffs}" : '',
    ].findAll { it }.join(' ')
    def sum_arg = sum_list ? "--fcs_summary ${sum_list} --assembly_out fcs_assemblies.tsv" : ''

    def excl_cmd = (excl_srcs || sum_arg) ? """
    build_mask.py ${excl_srcs} ${sum_arg} ${refs_arg} ${map_arg} ${pad_arg} \\
        --actions '${actions}' ${div_arg} ${minlen} \\
        --out mask.exclude.bed \\
        --summary_out mask.exclude.summary.tsv \\
        --json_out mask.exclude.build.json \\
        ${ has_fcs ? '--subset_out fcs_subset.txt' : '' }
    # nothing to exclude -> drop the empty BED so scoring stays unmasked
    if ! grep -qv '^#' mask.exclude.bed; then rm -f mask.exclude.bed; fi
    """ : ''

    def tgt_cmd = (target && tgt_srcs) ? """
    build_mask.py ${tgt_srcs} ${refs_arg} ${map_arg} ${sel_arg} ${ftype_arg} ${pad_arg} \\
        --out mask.target.bed \\
        --summary_out mask.target.summary.tsv \\
        --json_out mask.target.build.json
    if ! grep -qv '^#' mask.target.bed; then
        echo "ERROR: --mask_mode target, but no target interval matched any aligned reference." >&2
        echo "       Check --mask_bed / --mask_feature_select against the reference accessions." >&2
        exit 1
    fi
    """ : ''
    """
    ${excl_cmd}
    ${tgt_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}
