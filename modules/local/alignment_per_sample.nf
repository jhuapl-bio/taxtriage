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
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(meta), path(bamfiles), path(bai), path(mapping), path(bedgraph), path(covfile), path(k2_report), path(ch_diamond_analysis), path(fastas), path(microbert_report), path(pathogens_list), path(sampletype_thresholds_file)
    file assembly
    path(annotate_report)
    val minmapq
    path(taxdump)
    path(negative_control_jsons)
    path(positive_control_jsons)
    path(insilico_control_jsons)

    output:
        path "versions.yml"           , emit: versions
        tuple val(meta), path("*.json")    , optional:false, emit: txt
        tuple val(meta), path("*_removal_stats_by_taxid.xlsx")    , optional:true, emit: removal_by_taxid
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
    def output = "${meta.id}.paths.json"
    def id = meta.id
    def minmapq = minmapq ? " --minmapq ${minmapq} " :  ""
    def type = meta.type ? " -t ${meta.type} " : " -t Unknown "
    def min_reads_align = params.min_reads_align  ? " -r ${params.min_reads_align} " : " -r 3 "
    // def assemblyi = assembly ? " -j ${assembly} " : " "
    def read_count = meta.read_count && meta.read_count > 0 ? " --readcount ${meta.read_count} " : " "
    def cpu_count = task.cpus ? " -X ${task.cpus} "  : ""
    def k2 = k2_report.name == "NO_FILE" ? " " : " --k2 ${k2_report} "
    def mapping = mapping.name != "NO_FILE" ? "-m $mapping " : " "
    def bedgraph = bedgraph.name != "NO_FILE_bedgraph" ? "-b $bedgraph" :  " "
    def diamond_output = ch_diamond_analysis.name != "NO_FILE2" ? " --diamond $ch_diamond_analysis" : " "
    def output_dir = "search_results"
    def compress_species = params.rank ? " --rank ${params.rank} " : " "
    def pident = params.pident ? " --pident ${params.pident} " : " "

    /* groovylint-disable-next-line UnnecessaryCollectCall */
    def fastas = fastas && fastas.size() > 0 ? " -f ${fastas} " : " "
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
    def annotate_report_arg = annotate_report.name != "NO_FILE_annotate_report" ? " --annotate_report ${annotate_report} " : " "

    // Control sample arguments
    def ctrl_type = meta.control_type ? " --control_type ${meta.control_type} " : " "
    def neg_ctrls = negative_control_jsons.name != "NO_FILE_neg_ctrl" ? " --negative_controls ${negative_control_jsons} " : " "
    def pos_ctrls = positive_control_jsons.name != "NO_FILE_pos_ctrl" ? " --positive_controls ${positive_control_jsons} " : " "
    def insilico_ctrls = insilico_control_jsons.name != "NO_FILE_insilico_ctrl" ? " --insilico_controls ${insilico_control_jsons} " : " "

    """


    match_paths.py \\
        -i $bamfiles \\
        -o $output $bedgraph \\
        -s $id \\
        $type $read_count \\
        --output_dir $output_dir $fastas $cpu_count \\
        --scaled 8000 \\
        $alpha \\
        --min_threshold 0.002 \\
        -p $pathogens_list  $mapping $k2 $gap_allowance $jump_threshold \\
        --min_similarity_comparable 0.8 --taxid_removal_stats \\
        $breadth_weight $disparity_score_weight $gini_weight $minhash_weight $mapq_weight $hmp_weight \\
        $auto_score_power $score_power $depth_concentration_power \\
        --fast \\
        $min_reads_align $compress_species $mbert_report $minmapq $fast $taxonomy $enable_matrix $ani_threshold \\
        $workflow_revision $commitID $platform $sampletype_thresholds \\
        $ctrl_type $neg_ctrls $pos_ctrls $insilico_ctrls $reward_factor $dispersion_factor \\
        $mapq_breadth_power $mapq_gini_power \\
        $annotate_report_arg $pident

    cp search_results/removal_stats.xlsx "${meta.id}_removal_stats.xlsx" || true
    cp search_results/removal_stats_by_taxid.xlsx "${meta.id}_removal_stats_by_taxid.xlsx" || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version)
    END_VERSIONS

    """
}


