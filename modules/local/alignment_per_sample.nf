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
    label 'process_standard'

    conda (params.enable_conda ? "bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/taxtriage_confidence:2.1' :
        'jhuaplbio/taxtriage_confidence:2.1' }"

    input:
    tuple val(meta), path(bamfiles), path(bai), path(mapping), path(bedgraph), path(covfile), path(k2_report), path(ch_diamond_analysis), path(fastas), path(pathogens_list)
    file assembly

    output:
        path "versions.yml"           , emit: versions
        tuple val(meta), path("*.txt")    , optional:false, emit: txt

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/

    def output = "${meta.id}.paths.txt"
    def id = meta.id
    def type = meta.type ? " -t ${meta.type} " : " -t Unknown "
    def min_reads_align = params.min_reads_align  ? " -r ${params.min_reads_align} " : " -r 3 "
    def assemblyi = assembly ? " -j ${assembly} " : " "
    def read_count = meta.read_count && meta.read_count > 0 ? " --readcount ${meta.read_count} " : " "
    def cpu_count = task.cpus ? " -X ${task.cpus} "  : ""
    // if k2_report exists and is not null add --k2 flag
    // if k2_report != NO_FILE, add --k2 flag
    def k2 = k2_report.name == "NO_FILE" ? " " : " --k2 ${k2_report} "
    def mapping = mapping.name != "NO_FILE" ? "-m $mapping " : " "
    def bedgraph = bedgraph.name != "NO_FILE" ? "-b $bedgraph" :  " "
    def diamond_output = ch_diamond_analysis.name != "NO_FILE2" ? " --diamond $ch_diamond_analysis" : " "
    def ignore_alignment = params.ignore_missing ? " --ignore_missing_inputs " : " "
    def output_dir = "search_results"
    def compress_species = params.compress_species ? " --compress_species " : " "

    /* groovylint-disable-next-line UnnecessaryCollectCall */
    def fastas = fastas && fastas.size() > 0 ? " -f ${fastas} " : " "
    def minhash_weight = params.minhash_weight ? " --minhash_weight ${params.minhash_weight} " : " --minhash_weight 0.05 "
    def mapq_weight = params.mapq_weight ? " --mapq_weight ${params.mapq_weight} " : " --mapq_weight 0.0 "
    def hmp_weight = params.hmp_weight ? " --hmp_weight ${params.hmp_weight} " : " --hmp_weight 0.0 "
    def gini_weight = params.gini_weight ? " --gini_weight ${params.gini_weight} " : " --gini_weight 0.70 "
    def disparity_score_weight = params.disparity_score_weight ? " --disparity_score_weight ${params.disparity_score_weight} " : " "
    def breadth_weight = params.breadth_weight ? " --breadth_weight ${params.breadth_weight} " : " --breadth_weight 0.25 "
    def reward_factor = params.reward_factor ? " --reward_factor ${params.reward_factor} " : "  "
    def dispersion_factor = params.dispersion_factor ? " --dispersion_factor ${params.dispersion_factor} " : " "
    def sensitive = params.sensitive ? " --sensitive " : " "
    def gap_allowance = params.gap_allowance ? " --gap_allowance ${params.gap_allowance} " : " "
    def jump_threshold = params.jump_threshold ? " --jump_threshold ${params.jump_threshold} " : " "

    """

    match_paths.py \\
        -i $bamfiles \\
        -o $output $bedgraph \\
        -s $id $assemblyi \\
        $type $read_count \\
        --output_dir $output_dir $fastas $cpu_count \\
        --scaled 8000 \\
        --alpha 1.5 \\
        --min_threshold 0.002 \\
        -p $pathogens_list  $mapping $k2 $sensitive $gap_allowance $jump_threshold \\
        --min_similarity_comparable 0.8 \\
        $breadth_weight $disparity_score_weight $gini_weight $minhash_weight $mapq_weight $hmp_weight \\
        --fast \\
        $min_reads_align $ignore_alignment $compress_species

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version)
    END_VERSIONS

    """
}


