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
process TOP_HITS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://bmerritt1762/jhuaplbio/reportlab-pdf:4.0.7' :
        'jhuaplbio/reportlab-pdf:4.0.7' }"

    input:
    tuple val(meta), path(report), path(distributions), path(pathogens)

    output:
    path "versions.yml"           , emit: versions
    tuple val(meta), path("*top_report.tsv"), optional:true, emit: tops
    tuple val(meta), path("*toptaxids.txt"), optional:true, emit: taxids
    tuple val(meta), path("*topnames.txt"), optional:true, emit: names
    tuple val(meta), path("*.krakenreport_mqc.tsv"), optional:true, emit: krakenreport



    when:
    task.ext.when == null || task.ext.when


    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    def id = "${meta.id}"
    def zscore_accepted = params.zscore_accepted && params.zscore_accepted != "" ? " -z ${params.zscore_accepted} " : ' -z 1.5 '
    ch_top_per_taxa = ""
    def top_per_taxa  = params.top_per_taxa && params.top_per_taxa != "" ? " -s ${params.top_per_taxa} " : ''
    def top_hits_count = params.top_hits_count ? " -t ${params.top_hits_count}" : ' -t 10 '
    def distribution_arg = params.add_irregular_top_hits && distributions.name != "NO_FILE" ? " -d $distributions  $zscore_accepted " : ""
    def site = meta.type ? " -b ${meta.type} " : ""
    def pathogen_sheet = pathogens.name != "NO_FILE" ? " -p $pathogens  " : ""
    def remove_commensal = params.remove_commensal ? " --remove_commensals " : ""

    """
    echo ${meta.id} "-----------------META variable------------------"
    get_top_hits.py \\
        -i \"$report\" \\
        -o ${id}.top_report.tsv  ${distribution_arg} ${pathogen_sheet} ${site} \\
        $top_hits_count  $top_per_taxa $remove_commensal



    awk -F '\\t' -v id=${id} \\
        'BEGIN{OFS=\"\\t\"} { if (NR==1){ print \"Sample_Taxid\", \$2, \$1, \$4, \$6} else { \$5 = id\"_\"\$5;  print \$5, \$2, \$1, \$4, \$6  }}'  ${id}.top_report.tsv > ${id}.krakenreport_mqc.tsv

    awk -F '\\t' 'NR>1 {if (\$4 ~ "^S"){print \$6}}' ${id}.top_report.tsv | sort | uniq > ${id}.topnames.txt
    awk -F '\\t' 'NR>1 {if (\$4 ~ "^S"){print \$5}}' ${id}.top_report.tsv | sort | uniq > ${id}.toptaxids.txt

    if [ ! -s ${id}.topnames.txt ]; then
        echo "No top hits found"
        rm ${id}.top_report.tsv ${id}.krakenreport_mqc.tsv ${id}.topnames.txt ${id}.toptaxids.txt
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS

    """
}

