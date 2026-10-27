/**
 * BUILD_SAMPLE_MASK
 *
 * Region mask for ONE sample, from its samplesheet `mask_fasta` / `mask_bed`
 * columns (docs/samplesheet.md "Per-sample references and masks"), merged with
 * the run-wide mask from BUILD_MASK when there is one (--fcs, --mask_bed, ...).
 *
 * Output files are named exactly like BUILD_MASK's, so ALIGNMENT_PER_SAMPLE picks
 * them up unchanged:
 *   mask.exclude.bed     never-scored regions: masked bases of mask_fasta, mask_bed
 *                        (exclude mode) and the run-wide exclusions
 *   mask.products.bed    --mask_products: this sample's product regions (default CDS) from
 *                        its FEATURES_DOWNLOAD BED or --mask_products_features
 *   mask.target.bed      mask_mode target: the regions in mask_bed (or the run-wide
 *                        targets when the row sets mask_mode target without a BED)
 *   fcs_assemblies.tsv   passed through from the run-wide mask
 * A NO_FILE_sample_mask placeholder is always present, so a row whose mask turns
 * out empty still yields a (no-op) mask.
 *
 * Per-row overrides: meta.mask_mode, meta.mask_fasta_case (else --mask_mode,
 * --mask_fasta_case). Published to <outdir>/masking/samples/<sample>/.
 */
process BUILD_SAMPLE_MASK {
    tag "$meta.id"
    label 'process_low'
    // closure path: evaluated per task, where meta is defined (Nextflow strict syntax)
    publishDir path: { "${params.outdir}/masking/samples/${meta.id}" }, mode: params.publish_dir_mode,
        saveAs: { fn ->
            def n = fn.tokenize('/')[-1]
            (n.equals('versions.yml') || n.startsWith('NO_FILE')) ? null : n
        }

    conda (params.enable_conda ? "bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.2' :
        'jhuaplbio/taxtriage_confidence:2.2' }"

    input:
    tuple val(meta), path(bams, stageAs: 'bams/?/*'), path(sample_masks, stageAs: 'sample_mask/?/*'), path(products, stageAs: 'products/?/*')
    path(run_mask, stageAs: 'run_mask/*')

    output:
    tuple val(meta), path("out/*"), emit: files
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def base     = { f -> f.toString().tokenize('/')[-1] }
    def real     = { files -> (files instanceof List ? files : [files]).findAll { it && !base.call(it).startsWith('NO_FILE') } }
    def _is_fa   = { f -> base.call(f) ==~ /.*\.(fa|fasta|fna|fas|fsa)(\.gz)?$/ }
    def bam_list = real.call(bams).findAll { base.call(it) ==~ /.*\.(bam|cram)$/ }.join(' ')
    def srcs     = real.call(sample_masks)
    def fa_list  = srcs.findAll { _is_fa.call(it) }.join(' ')
    def bed_list = srcs.findAll { !_is_fa.call(it) }.join(' ')
    def run      = real.call(run_mask)
    def run_excl = run.find { base.call(it) == 'mask.exclude.bed' }
    def run_tgt  = run.find { base.call(it) == 'mask.target.bed' }
    def run_asm  = run.find { base.call(it) == 'fcs_assemblies.tsv' }
    def target   = (meta.mask_mode ?: params.mask_mode ?: 'exclude') == 'target'
    // --mask_products: this sample's product regions (BED from FEATURES_DOWNLOAD, or
    // --mask_products_features: BED / NCBI feature table / GFF3) become its targets
    def prod_on  = (params.mask_products && params.mask_products.toString() != 'false') ? true : false
    def prods    = real.call(products)
    def _is_gff  = { f -> base.call(f) ==~ /.*\.gff3?(\.gz)?$/ }
    def _is_bed  = { f -> base.call(f) ==~ /.*\.bed(\.gz)?$/ }
    def p_bed    = prods.findAll { _is_bed.call(it) }.join(' ')
    def p_gff    = prods.findAll { _is_gff.call(it) }.join(' ')
    def p_ftab   = prods.findAll { !_is_bed.call(it) && !_is_gff.call(it) }.join(' ')
    def p_srcs   = [p_bed ? "--bed ${p_bed}" : '', p_gff ? "--gff ${p_gff}" : '', p_ftab ? "--feature_table ${p_ftab}" : ''].findAll { it }.join(' ')
    def sel_arg  = params.mask_feature_select ? "--feature_select '${params.mask_feature_select}' --select_bed_names" : ''
    def fa_case  = meta.mask_fasta_case ?: params.mask_fasta_case ?: 'hard'
    def refs_arg = bam_list ? "--refs ${bam_list}" : ''
    def pad_arg  = params.mask_pad ? "--pad ${params.mask_pad}" : ''

    // build_mask.py --bed is nargs='+' (a repeated --bed REPLACES the earlier one),
    // so every exclusion BED goes after a single --bed
    def excl_beds = ([(!target && bed_list) ? bed_list : '', run_excl ? run_excl.toString() : ''].findAll { it }).join(' ')
    def excl_srcs = [
        fa_list ? "--masked_fasta ${fa_list} --masked_fasta_case ${fa_case} --masked_fasta_min_run ${params.mask_fasta_min_run ?: 1}" : '',
        excl_beds ? "--bed ${excl_beds}" : '',
    ].findAll { it }.join(' ')
    def tgt_beds = target ? (bed_list ?: (run_tgt ? run_tgt.toString() : '')) : ''

    def excl_cmd = excl_srcs ? """
    build_mask.py ${excl_srcs} ${refs_arg} ${pad_arg} \\
        --out out/mask.exclude.bed \\
        --summary_out out/mask.exclude.summary.tsv \\
        --json_out out/mask.exclude.build.json
    # nothing to exclude -> drop the empty BED so this sample stays unmasked
    if ! grep -qv '^#' out/mask.exclude.bed; then rm -f out/mask.exclude.bed; fi
    """ : ''
    // product regions: exact sequence-name matches only (--feature_exact), so annotation
    // of an assembly a contig was merely matched to is never placed on that contig
    def prod_cmd = (prod_on && p_srcs) ? """
    build_mask.py ${p_srcs} --feature_exact ${sel_arg} --feature_types '${params.mask_products_types ?: 'CDS'}' ${refs_arg} \\
        --out out/mask.products.bed \\
        --summary_out out/mask.products.summary.tsv \\
        --json_out out/mask.products.build.json
    if ! grep -qv '^#' out/mask.products.bed; then
        echo "WARNING: sample ${meta.id}: --mask_products: no product region matched an aligned reference by exact name; this sample is not restricted to products." >&2
        rm -f out/mask.products.bed
    fi
    """ : (prod_on ? """
    echo "WARNING: sample ${meta.id}: --mask_products: no features for this sample (FEATURES_DOWNLOAD found none, or no --mask_products_features)." >&2
    """ : '')
    // target = the row's own target BED (or the run-wide one) + the product regions
    def tgt_cmd = (target || prod_on) ? """
    TGT="${tgt_beds}"
    if [ -s out/mask.products.bed ]; then TGT="\$TGT out/mask.products.bed"; fi
    if [ -n "\$(echo \$TGT)" ]; then
        build_mask.py --bed \$TGT ${refs_arg} ${pad_arg} \\
            --out out/mask.target.bed \\
            --summary_out out/mask.target.summary.tsv \\
            --json_out out/mask.target.build.json
        if ! grep -qv '^#' out/mask.target.bed; then
            if [ -n "${target ? 'yes' : ''}" ]; then
                echo "ERROR: sample ${meta.id}: mask_mode target, but no target interval in mask_bed matched an aligned reference." >&2
                exit 1
            fi
            rm -f out/mask.target.bed
        fi
    fi
    """ : ''
    """
    mkdir -p out
    touch out/NO_FILE_sample_mask
    ${excl_cmd}
    ${prod_cmd}
    ${tgt_cmd}
    ${ run_asm ? "cp ${run_asm} out/fcs_assemblies.tsv" : '' }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}
