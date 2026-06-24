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
    // Pathogen reference sheet (assets/pathogen_sheet.csv or --pathogens override). Used to flag
    // listed pathogens that have NO reference alignment but appear in novelty / VF-AMR results.
    // Pass a NO_FILE placeholder to disable the cross-reference.
    path(pathogens_sheet)
    // bvbrc specialty-gene reference TSV (source_id -> taxids). Lets VF/AMR pathogen matching
    // key on canonical taxid instead of the merged sheet's Genus/Species text. NO_FILE to skip.
    path(vfamr_taxid_tsv)
    // Standalone per-sample annotate_report.tsv files (annotate_report.py output). Carry de-novo /
    // unaligned VF-AMR hits for samples with NO reference alignment, whose annotation is otherwise
    // dropped from the merged XLSX. Pass a NO_FILE placeholder to disable. Single collected channel.
    path(annotate_reports)

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

    // Pathogen sheet cross-reference (skip on NO_FILE placeholder / '~' rename artifact).
    def path_arg = ''
    if (pathogens_sheet && !pathogens_sheet.name.startsWith('NO_FILE') && !pathogens_sheet.name.startsWith('~')) {
        path_arg = "--pathogens ${pathogens_sheet}"
    }

    // bvbrc source-id -> taxids reference for VF/AMR pathogen matching by taxid.
    def vfamr_tax_arg = ''
    if (vfamr_taxid_tsv && !vfamr_taxid_tsv.name.startsWith('NO_FILE') && !vfamr_taxid_tsv.name.startsWith('~')) {
        vfamr_tax_arg = "--vfamr-taxids ${vfamr_taxid_tsv}"
    }

    // ── Standalone annotate_report.tsv files (de-novo / unaligned VF-AMR) ──────
    // Supplement annotation for samples with no reference alignment. Filter NO_FILE
    // placeholders and Nextflow '~' same-name conflict copies; dedupe by basename.
    def annot_arg = ''
    if (annotate_reports) {
        def annot_list = annotate_reports instanceof List ? annotate_reports : [annotate_reports]
        def seen_annot = [] as Set
        def valid_annot = annot_list.findAll { f ->
            f && !f.name.startsWith('NO_FILE') && !f.name.startsWith('~') &&
            (f.name.endsWith('.tsv') || f.name.endsWith('.xlsx')) &&
            seen_annot.add(f.name)
        }
        def annot_files = valid_annot.join(' ')
        if (annot_files) {
            annot_arg = "--annotate_reports ${annot_files}"
        }
    }

    """
    make_report.py -i ${json_inputs} \\
        -t ${template} \\
        -o ${output_html} \\
        ${prot_arg} ${pident} ${mintass} \\
        ${nov_arg} ${nov_dl_arg} ${path_arg} ${vfamr_tax_arg} ${annot_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS

    """
}


