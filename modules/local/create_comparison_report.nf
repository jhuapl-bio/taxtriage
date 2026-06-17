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
process CREATE_COMPARISON_REPORT {
    tag "comparison_full_report"
    label 'process_medium'
    publishDir "${params.outdir}/report", mode: 'copy'

    conda (params.enable_conda ? "bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/jhuaplbio/reportlab-pdf:4.0.9' :
        'jhuaplbio/reportlab-pdf:4.0.9' }"

    input:
    // One or more .paths.json files from ALIGNMENT_PER_SAMPLE (collected across all samples)
    path(json_files)
    // HTML template (heatmap.html)
    path(template)
    // Optional protein-annotation XLSX files from ORGANISM_MERGE_REPORT --output_annot_xlsx
    // Pass a NO_FILE placeholder when protein annotations are not available
    path(protein_annotations)
    // Optional novelty files from NOVELTY_COLLECT: per-sample + combined JSON/XLSX. The combined
    // all.novelty.json within this set drives the Novelty panel; the rest become download links.
    // Pass a NO_FILE placeholder when novelty detection did not run. Single input (not two) so the
    // combined json isn't staged twice -> avoids an "input file name collision" on all.novelty.json.
    path(novelty_files)
    // Optional embedding files from NOVEL_HOMOLOGS: per-sample *.umap.tsv, *.cluster_summary.tsv,
    // and *.clusters.tsv. Pass a NO_FILE placeholder when the step did not run.
    path(embedding_files)

    output:
        path "versions.yml"           , emit: versions
        path("*odr.html")                 , optional: true, emit: html

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/
    def output_html = "all.odr.html"

    // Build the list of JSON input files (filter out any NO_FILE placeholders)
    def json_inputs = json_files instanceof List
        ? json_files.findAll { it.name != 'NO_FILE' && it.name.endsWith('.json') }.join(' ')
        : (json_files.name != 'NO_FILE' && json_files.name.endsWith('.json') ? json_files.toString() : '')

    // Build optional protein annotations argument.
    // Filter out NO_FILE placeholders AND Nextflow conflict-rename artifacts
    // (files staged as "~original_name" when two inputs share the same filename).
    def prot_arg = ''
    if (protein_annotations) {
        def prot_list = protein_annotations instanceof List ? protein_annotations : [protein_annotations]
        def valid_prot = prot_list.findAll { f ->
            !f.name.startsWith('NO_FILE') &&
            f.name != 'NO_FILE' &&
            !f.name.startsWith('~')   // drop Nextflow same-name conflict copies
        }
        // Deduplicate by basename to avoid passing the same XLSX twice
        def seen_names = [] as Set
        def deduped_prot = valid_prot.findAll { f -> seen_names.add(f.name) }
        def prot_files = deduped_prot.join(' ')
        if (prot_files) {
            prot_arg = "-p ${prot_files}"
        }
    }
    def pident = params.pident ? " --pident ${params.pident} " : " "
    def mintass = params.mintass ? " --mintass ${params.mintass} " : " "

    // ── Novelty panel feed + download links ───────────────────────────────────
    // One staged set of files: pick the combined all.novelty.json as the -n feed and expose the
    // whole set as download links. Filter NO_FILE placeholders and Nextflow '~' rename artifacts.
    def nov_list = novelty_files instanceof List ? novelty_files : [novelty_files]
    def nov_valid = nov_list.findAll { f ->
        f && !f.name.startsWith('NO_FILE') && !f.name.startsWith('~') &&
        (f.name.endsWith('.json') || f.name.endsWith('.xlsx'))
    }
    def nov_combined = nov_valid.find { it.name == 'all.novelty.json' }
    def nov_arg = nov_combined ? "-n ${nov_combined}" : ''
    def nov_dl_files = nov_valid.join(' ')
    def nov_dl_arg = nov_dl_files ? "--novelty-downloads ${nov_dl_files}" : ''

    // ── Embedding / UMAP cluster files (NOVEL_HOMOLOGS) ──────────────────────
    def emb_list = embedding_files instanceof List ? embedding_files : [embedding_files]
    def emb_valid = emb_list.findAll { f ->
        f && !f.name.startsWith('NO_FILE') && !f.name.startsWith('~') &&
        (f.name.endsWith('.umap.tsv') || f.name.endsWith('.cluster_summary.tsv') || f.name.endsWith('.clusters.tsv'))
    }
    def emb_arg = emb_valid ? "--embedding ${emb_valid.join(' ')}" : ''

    """
    make_report.py -i ${json_inputs} \\
        -t ${template} \\
        -o ${output_html} \\
        ${prot_arg} ${pident} ${mintass} \\
        ${nov_arg} ${nov_dl_arg} \\
        ${emb_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS

    """
}


