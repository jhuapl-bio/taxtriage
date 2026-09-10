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
    // HTML template (heatmap.html) — a thin shell that references its CSS/JS as
    // external parts under assets/src/. make_report.py inlines them into one
    // self-contained report.
    path(template)
    // The assets/src directory (css/ + js/) staged next to the template so the
    // template's external references resolve inside the task workdir.
    path(src_assets)
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
    // Optional directory of local CDN library copies (d3, xlsx, jspdf, Leaflet,
    // Font Awesome + their fonts/images) for a fully offline report build. When a
    // real directory is staged (params.offline_report_files set) make_report.py
    // embeds those files inline; otherwise a NO_FILE placeholder is passed and the
    // report either downloads the libs at build time (params.offline_report) or
    // leaves them as CDN links (default). Single input.
    path(offline_report_files)
    // Optional JSON file holding a full sample-QC rule list (params.report_flag_rules).
    // Replaces every report_flag_* param when supplied. NO_FILE placeholder to skip.
    path(flag_rules_json)
    // Optional per-dataset .paths.json files for the in-silico subsample datasets (from
    // ALIGNMENT_PER_SAMPLE_INSILICO). Used ONLY to build the In-Silico suite tab — NOT added to
    // the main multi-run heatmap/table. Pass a NO_FILE placeholder when --sim_subsample is off.
    path(insilico_json)
    // Optional *_subsample_manifest.tsv file(s) from SUBSAMPLE_INSILICO. Provide the In-Silico
    // suite tab with authoritative target/actual read counts + master totals + seeds. NO_FILE to skip.
    path(insilico_manifests)
    // Optional insilico_params.json describing the subsampling run parameters (mode, series,
    // replicates, seed, sim_nreads, iss_model, ...). Populates the suite tab's provenance panel. NO_FILE to skip.
    path(insilico_params)

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
    // Propagate the same --min_conf used upstream (match_paths.py/create_report.py)
    // so the HTML report's global AND per-sample-type sliders default to it too,
    // instead of the auto-computed best_cutoffs recommendation.
    def min_conf_arg = params.min_conf != null ? " --min_conf ${params.min_conf} " : " "

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

    // ── Offline report embedding ──────────────────────────────────────────────
    // A staged directory (not a NO_FILE placeholder) -> embed those local library
    // copies inline. Else if params.offline_report -> download + embed at build
    // time. Else (default) -> leave CDN links so the report fetches libs on load.
    def offline_arg = ''
    if (offline_report_files && !offline_report_files.name.startsWith('NO_FILE')) {
        offline_arg = "--offline_report_files ${offline_report_files}"
    } else if (params.offline_report) {
        offline_arg = "--offline_report"
    }

    // ── Sample-QC flag defaults ───────────────────────────────────────────────
    // params.report_flag_* seed the report's whole-sample rule set. Nothing here
    // filters the RUN — make_report.py just bakes the rules into the HTML so the
    // report opens with them applied (and the user can edit them live).
    def flag_bits = []
    def has_flag_rules_file = flag_rules_json && !flag_rules_json.name.startsWith('NO_FILE') && !flag_rules_json.name.startsWith('~')
    if (has_flag_rules_file) {
        flag_bits << "--flag-rules ${flag_rules_json}"
    } else {
        if (params.report_flag_min_reads != null)         flag_bits << "--flag-min-reads ${params.report_flag_min_reads}"
        if (params.report_flag_min_aligned_reads != null) flag_bits << "--flag-min-aligned-reads ${params.report_flag_min_aligned_reads}"
        if (params.report_flag_min_organisms != null)     flag_bits << "--flag-min-organisms ${params.report_flag_min_organisms}"
        if (params.report_flag_organism_tass != null)     flag_bits << "--flag-organism-tass ${params.report_flag_organism_tass}"
        if (params.report_flag_min_detections != null)    flag_bits << "--flag-min-detections ${params.report_flag_min_detections}"
        if (params.report_flag_metadata)                  flag_bits << "--flag-metadata '${params.report_flag_metadata}'"
    }
    // Only meaningful alongside a criterion; emitting them on their own made a
    // run with NO rules look configured in .command.sh.
    if (flag_bits) {
        if (params.report_flag_logic)  flag_bits << "--flag-logic ${params.report_flag_logic}"
        if (params.report_flag_action) flag_bits << "--flag-action ${params.report_flag_action}"
        if (params.report_flag_missing) flag_bits << "--flag-missing"
        if (params.report_flag_view && params.report_flag_view != 'all') flag_bits << "--flag-view ${params.report_flag_view}"
        // An empty string is a real setting ("count host too"), so test for null
        // rather than truthiness -- otherwise it would silently fall back to 9606.
        if (params.report_flag_exclude_taxids != null) flag_bits << "--flag-exclude-taxids '${params.report_flag_exclude_taxids}'"
    }
    def flag_arg = flag_bits.join(' ')
    // ── In-silico subsampling suite feed ──────────────────────────────────────
    // Subsample dataset JSON(s) — used only to build the suite tab (not the heatmap).
    def insil_json_arg = ''
    if (insilico_json) {
        def ij_list = insilico_json instanceof List ? insilico_json : [insilico_json]
        def seen_ij = [] as Set
        def valid_ij = ij_list.findAll { f ->
            f && !f.name.startsWith('NO_FILE') && !f.name.startsWith('~') &&
            f.name.endsWith('.json') && seen_ij.add(f.name)
        }
        def ij_files = valid_ij.join(' ')
        if (ij_files) {
            insil_json_arg = "--insilico_json ${ij_files}"
        }
    }
    // Subsample manifest(s): drop NO_FILE placeholders + Nextflow '~' rename copies,
    // dedupe by basename.
    def insil_manifest_arg = ''
    if (insilico_manifests) {
        def man_list = insilico_manifests instanceof List ? insilico_manifests : [insilico_manifests]
        def seen_man = [] as Set
        def valid_man = man_list.findAll { f ->
            f && !f.name.startsWith('NO_FILE') && !f.name.startsWith('~') &&
            f.name.endsWith('.tsv') && seen_man.add(f.name)
        }
        def man_files = valid_man.join(' ')
        if (man_files) {
            insil_manifest_arg = "--insilico_manifests ${man_files}"
        }
    }
    // Params JSON (single file).
    def insil_params_arg = ''
    if (insilico_params && !insilico_params.name.startsWith('NO_FILE') && !insilico_params.name.startsWith('~')) {
        insil_params_arg = "--insilico_params ${insilico_params}"
    }

    """
    make_report.py -i ${json_inputs} \\
        -t ${template} \\
        -o ${output_html} \\
        ${prot_arg} ${pident} ${mintass} ${min_conf_arg} \\
        ${nov_arg} ${nov_dl_arg} ${path_arg} ${vfamr_tax_arg} ${annot_arg} ${offline_arg} \\
        ${flag_arg} \\
        ${insil_json_arg} ${insil_manifest_arg} ${insil_params_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS

    """
}


