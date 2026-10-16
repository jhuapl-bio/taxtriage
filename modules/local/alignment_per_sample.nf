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
process ALIGNMENT_PER_SAMPLE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.2' :
        'jhuaplbio/taxtriage_confidence:2.2' }"

    input:
    tuple val(meta),
        path(bamfiles),
        path(bai),
        path(mapping),
        path(bedgraph),
        path(covfile),
        path(k2_report),
        path(ch_diamond_analysis),
        path(fastas),
        path(ch_annotate_report),
        path(microbert_report),
        path(pathogens_list),
        path(sampletype_thresholds_file)
    file assembly
    val minmapq
    path(taxdump)
    path(negative_control_jsons)
    path(positive_control_jsons)
    path(insilico_control_jsons)
    path(meta_csv)
    // Region mask from BUILD_MASK (mask.exclude.bed / mask.target.bed /
    // fcs_assemblies.tsv) or a NO_FILE_mask placeholder. See docs/masking.md.
    path(mask_files, stageAs: 'mask/*')

    output:
        path "versions.yml"           , emit: versions
        tuple val(meta), path("*.json")    , optional:false, emit: txt
        tuple val(meta), path("*_removal_stats_by_taxid.xlsx")    , optional:true, emit: removal_by_taxid
        // --mask_write_bam: the input BAM with masked alignments removed (+ .bai)
        tuple val(meta), path("*.masked.bam*")    , optional:true, emit: masked_bam
    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    // if params.sensitive is true then set breadth_weight to 0.3, gini_weight to 0.26 and minhash_weight to 0.44
    def minhash_weight = params.minhash_weight ? " --minhash_weight ${params.minhash_weight} " : " "
    def gini_weight = params.gini_weight ? " --gini_weight ${params.gini_weight} " : " "
    def breadth_weight = params.breadth_weight ? " --breadth_weight ${params.breadth_weight} " : " "
    def mapq_weight = params.mapq_weight ? " --mapq_weight ${params.mapq_weight} " : " "
    def hmp_weight = params.hmp_weight ? " --hmp_weight ${params.hmp_weight} " : " "
    def disparity_score_weight = params.disparity_score_weight ? " --disparity_weight ${params.disparity_score_weight} " : " "
    def auto_score_power = params.auto_score_power ? " --auto_score_power " : " "
    def score_power = params.score_power != null ? " --score_power ${params.score_power} " : " "
    def depth_concentration_power = params.depth_concentration_power != null ? " --depth_concentration_power ${params.depth_concentration_power} " : " "
    def mapq_breadth_power = params.mapq_breadth_power != null ? " --mapq_breadth_power ${params.mapq_breadth_power} " : " "
    def mapq_gini_power = params.mapq_gini_power != null ? " --mapq_gini_power ${params.mapq_gini_power} " : " "
    // Conflict read-removal controls. taxid_removal_stats defaults to true, which
    // preserves the previous behaviour of always emitting the per-taxid report.
    def no_read_removal = params.no_read_removal ? " --no_read_removal " : " "
    def taxid_removal_stats = params.taxid_removal_stats ? " --taxid_removal_stats " : " "
    // ext.suffix lets the --mask_compare unmasked baseline publish alongside
    // the masked run without a file-name clash (see conf/modules.config).
    def output = "${meta.id}${task.ext.suffix ?: ''}.paths.json"
    def id = meta.id
    def minmapq_arg = minmapq ? " --minmapq ${minmapq} " :  ""
    // ── Ambiguous-read (multimapper) rescue ───────────────────────────────────
    // Re-scores MAPQ-0 ties on their own alignment phred instead of dropping
    // them wholesale; see docs/multimapper-rescue.md.
    def rescue_multimapped = params.rescue_multimapped == false ? " --no_rescue_multimapped " : " "
    def rescue_min_aln_phred = params.rescue_min_aln_phred != null ? " --rescue_min_aln_phred ${params.rescue_min_aln_phred} " : " "
    def rescue_max_mapq = params.rescue_max_mapq != null ? " --rescue_max_mapq ${params.rescue_max_mapq} " : " "
    def rescue_min_aln_frac = params.rescue_min_aln_frac != null ? " --rescue_min_aln_frac ${params.rescue_min_aln_frac} " : " "
    def rescue_max_nm_rate = params.rescue_max_nm_rate != null ? " --rescue_max_nm_rate ${params.rescue_max_nm_rate} " : " "
    def rescue_min_aln_len = params.rescue_min_aln_len != null ? " --rescue_min_aln_len ${params.rescue_min_aln_len} " : " "
    def rescue_proper_pair = params.rescue_require_proper_pair ? " --rescue_require_proper_pair " : " "
    def rescue_as_highmapq = params.rescue_counts_as_highmapq ? " --rescue_counts_as_highmapq " : " "
    def rescue_args = rescue_multimapped + rescue_min_aln_phred + rescue_max_mapq + rescue_min_aln_frac +
                      rescue_max_nm_rate + rescue_min_aln_len + rescue_proper_pair + rescue_as_highmapq
    def type = meta.type ? " -t ${meta.type} " : " -t Unknown "
    def min_reads_align = params.min_reads_align  ? " -r ${params.min_reads_align} " : " -r 3 "
    // def assemblyi = assembly ? " -j ${assembly} " : " "
    def read_count = meta.read_count && meta.read_count > 0 ? " --readcount ${meta.read_count} " : " "
    // Input read file size, as supplied (reporting only: metadata.input_data)
    def _isz = WorkflowTaxtriage.inputSize(meta)
    def _inames   = _isz[1].size() <= 50 ? _isz[1].join(',').replace("'", "") : "${_isz[1].size()} files"
    def input_size = _isz[0] ? " --input_bytes ${_isz[0]} --input_files '${_inames}' " : " "
    def cpu_count = task.cpus ? " -X ${task.cpus} "  : ""
    def k2 = k2_report.name == "NO_FILE" ? " " : " --k2 ${k2_report} "
    def mapping_arg  = mapping.name  != "NO_FILE"         ? "-m $mapping"  : " "
    def bedgraph_arg = bedgraph.name != "NO_FILE_bedgraph" ? "-b $bedgraph" : " "
    def diamond_output = ch_diamond_analysis.name != "NO_FILE2" ? " --diamond $ch_diamond_analysis" : " "
    def output_dir = "search_results"
    def compress_species = params.rank ? " --rank ${params.rank} " : " "
    def pident = params.pident ? " --pident ${params.pident} " : " "
    // Windowed per-reference depth profile for the Alignment Trends tab / --alignment_trends
    def depth_profile = params.depth_profile_windows != null ? " --depth_profile_windows ${params.depth_profile_windows} " : " "

    /* groovylint-disable-next-line UnnecessaryCollectCall */
    def fastas_arg = fastas && fastas.size() > 0 ? " -f ${fastas} " : " "
    def reward_factor = params.reward_factor ? " --reward_factor ${params.reward_factor} " : "  "
    def dispersion_factor = params.dispersion_factor ? " --dispersion_factor ${params.dispersion_factor} " : " "
    def fast = params.fast ? "  " : " --compare_references "
    def gap_allowance = params.gap_allowance ? " --gap_allowance ${params.gap_allowance} " : " "
    def jump_threshold = params.jump_threshold ? " --jump_threshold ${params.jump_threshold} " : " "
    def mbert_report = microbert_report.name != "NO_FILEmicrobert" ? " --microbert ${microbert_report} " : " "
    def taxonomy  = taxdump ? " --taxdump ${taxdump} " : " "
    def alpha = params.alpha ? " --alpha ${params.alpha} " : " --alpha 1.0 "
    def enable_matrix = params.enable_matrix ? " --enable_matrix " : " "
    def ani_threshold = params.ani_threshold ? " --ani_threshold $params.ani_threshold " : ""
    def workflow_revision = workflow.revision ? " --workflow_revision ${workflow.revision} " : " --workflow_revision NA "
    def commitID = workflow.commitId ? " --commit_id ${workflow.commitId} " : " --commit_id NA "
    def platform = meta.platform ? " --platform ${meta.platform} " : " "
    def sampletype_thresholds = sampletype_thresholds_file.name != "NO_FILE_thresholds" ? " --thresholds_json ${sampletype_thresholds_file} " : " "
    def annotate_report_arg = ch_annotate_report.name != "NO_FILE_annotate_report" ? " --annotate_report ${ch_annotate_report} " : " "

    // ── Run-level metadata: pass the CSV file directly to match_paths.py ────────
    // Python does the sample lookup, which is simpler and avoids Nextflow
    // channel-combination quoting issues.
    def meta_csv_arg = meta_csv.name != "NO_FILE_meta_csv" ? " --meta_csv ${meta_csv} " : " "

    // ── Region masking (--mask_bed / --fcs / --mask_mode) ─────────────────────
    // stageAs 'mask/*' makes .name the staged relative path -> compare basenames
    def _mbase      = { f -> f.toString().tokenize('/')[-1] }
    def _mfiles     = (mask_files instanceof List ? mask_files : [mask_files]).findAll { it && !_mbase.call(it).startsWith('NO_FILE') }
    def _excl_bed   = _mfiles.find { _mbase.call(it) == 'mask.exclude.bed' }
    def _tgt_bed    = _mfiles.find { _mbase.call(it) == 'mask.target.bed' }
    def _fcs_asm    = _mfiles.find { _mbase.call(it) == 'fcs_assemblies.tsv' }
    def _is_target  = (params.mask_mode ?: 'exclude') == 'target' && _tgt_bed
    def mask_args   = ''
    if (_is_target) {
        mask_args = " --mask_bed ${_tgt_bed} --mask_mode target --mask_untargeted ${params.mask_untargeted ?: 'drop'} " +
                    (_excl_bed ? " --mask_exclude_bed ${_excl_bed} " : '')
    } else if (_excl_bed) {
        mask_args = " --mask_bed ${_excl_bed} --mask_mode exclude "
    }
    if (mask_args) {
        mask_args += (params.mask_read_overlap != null ? " --mask_read_overlap ${params.mask_read_overlap} " : '') +
                     (params.mask_penalty_power ? " --mask_penalty_power ${params.mask_penalty_power} " : '') +
                     (params.mask_min_region ? " --mask_min_region ${params.mask_min_region} " : '') +
                     (params.mask_write_bam ? " --mask_write_bam ${meta.id}.masked.bam " : '')
    }
    // the assembly-level FCS annotation is useful even when nothing is masked
    mask_args += _fcs_asm ? " --fcs_assembly_table ${_fcs_asm} " : ''

    // Control sample arguments
    def ctrl_type = meta.control_type ? " --control_type ${meta.control_type} " : " "
    def neg_ctrls = negative_control_jsons.name != "NO_FILE_neg_ctrl" ? " --negative_controls ${negative_control_jsons} " : " "
    def pos_ctrls = positive_control_jsons.name != "NO_FILE_pos_ctrl" ? " --positive_controls ${positive_control_jsons} " : " "
    def insilico_ctrls = insilico_control_jsons.name != "NO_FILE_insilico_ctrl" ? " --insilico_controls ${insilico_control_jsons} " : " "

    // Simulated datasets (ISS / NanoSim / spike-into-background) stamp their own
    // provenance into the output JSON, so every downstream consumer can tell them
    // apart from real samples without pattern-matching the sample id.
    def insilico_flag = meta.insilico ? " --insilico " : " "
    def parent_arg    = (meta.insilico && meta.parent_id) ? " --parent_sample ${meta.parent_id} " : " "

    """


    match_paths.py \\
        -i $bamfiles \\
        -o $output $bedgraph_arg \\
        -s $id \\
        $type $read_count $input_size \\
        --output_dir $output_dir $fastas_arg $cpu_count \\
        --scaled 8000 \\
        $alpha \\
        --min_threshold 0.002 \\
        -p $pathogens_list  $mapping_arg $k2 $gap_allowance $jump_threshold \\
        --min_similarity_comparable 0.8 $taxid_removal_stats $no_read_removal \\
        $breadth_weight $disparity_score_weight $gini_weight $minhash_weight $mapq_weight $hmp_weight \\
        $auto_score_power $score_power $depth_concentration_power \\
        --fast \\
        $min_reads_align $compress_species $mbert_report $minmapq_arg $fast $taxonomy $enable_matrix $ani_threshold \\
        $workflow_revision $commitID $platform $sampletype_thresholds \\
        $ctrl_type $neg_ctrls $pos_ctrls $insilico_ctrls $insilico_flag $parent_arg $reward_factor $dispersion_factor \\
        $mapq_breadth_power $mapq_gini_power $rescue_args \\
        $annotate_report_arg $pident $depth_profile \\
        $meta_csv_arg $mask_args

    cp search_results/removal_stats.xlsx "${meta.id}_removal_stats.xlsx" || true
    cp search_results/removal_stats_by_taxid.xlsx "${meta.id}_removal_stats_by_taxid.xlsx" || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version)
    END_VERSIONS

    """
}


