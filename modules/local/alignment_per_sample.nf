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
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.21.0--py39hcada746_1' :
        'jhuaplbio/taxtriage_confidence:2.0' }"

    input:
    tuple val(meta), path(bamfiles), path(bai), path(mapping), path(bedgraph), path(covfile), path(k2_report), path(ch_diamond_analysis), path(pathogens_list)
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

    // if k2_report exists and is not null add --k2 flag
    // if k2_report != NO_FILE, add --k2 flag

    def k2 = k2_report.name == "NO_FILE" ? " " : " --k2 ${k2_report} "
    def mapping = mapping.name != "NO_FILE" ? "-m $mapping " : " "
    def bedgraph = bedgraph.name != "NO_FILE" ? "-b $bedgraph" :  " "
    def diamond_output = ch_diamond_analysis.name != "NO_FILE2" ? " --diamond $ch_diamond_analysis" : " "
    def ignore_alignment = params.ignore_missing ? " --ignore_missing_inputs " : " "
    def output_dir = "search_results"

    """

    match_paths.py \\
        -i $bamfiles \\
        -o $output $bedgraph \\
        -s $id $assemblyi \\
        $type \\
        --output_dir $output_dir \\
        --scaled 50 \\
        --alpha 2.5 \\
        --min_threshold 0.002 \\
        -p $pathogens_list  $mapping \\
        --min_similarity_comparable 0.8 \\
        --gini_weight 0.75  \\
        --disparity_score_weight 0.0  \\
        --breadth_weight 0.15 --minhash_weight 0.1 --mapq_weight 0.0 --hmp_weight 0.0 \\
        --fast \\
        $min_reads_align $ignore_alignment

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version)
    END_VERSIONS

    """
}


