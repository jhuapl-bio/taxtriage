//
// Cross-check all upstream steps for final confidence report
//
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

include { ALIGNMENT_PER_SAMPLE as ALIGNMENT_PER_SAMPLE_CONTROLS } from '../../modules/local/alignment_per_sample'
include { ALIGNMENT_PER_SAMPLE } from '../../modules/local/alignment_per_sample'
include { ORGANISM_MERGE_REPORT } from '../../modules/local/report_merge'
include { ORGANISM_MERGE_REPORT as SINGLE_REPORT } from '../../modules/local/report_merge'
include { CREATE_COMPARISON_REPORT } from '../../modules/local/create_comparison_report'
include { MICROBERT_PREDICT } from '../../modules/local/microbert_predict'
include { CLUSTER_ALIGNMENT } from '../../modules/local/cluster_alignment'
include { MMSEQS_EASYCLUSTER } from '../../modules/local/mmseqs2_easycluster'
include { MICROBERT_PARSE } from '../../modules/local/microbert_parse'

workflow REPORT {
    take:
        alignments
        pathogens_list
        distributions
        assemblyfile
        ch_taxdump_dir
        all_samples
    main:
        ch_pathogens_report = Channel.empty()
        ch_pathognes_list = Channel.empty()
        // make the ch_report_microbert an empty path channel
        ch_report_microbert = Channel.empty()
        ch_empty_file3 = Channel.fromPath("$projectDir/assets/NO_FILE3", checkIfExists: true)
        ch_sampletype_thresholds = params.disable_auto_weights
            ? Channel.value(file("$projectDir/assets/NO_FILE_thresholds"))
            : Channel.fromPath("$projectDir/assets/sampletype_best_thresholds.json", checkIfExists: false)
                .ifEmpty { Channel.value(file("$projectDir/assets/NO_FILE_thresholds")) }

        // Placeholder files for when no control JSONs are available
        ch_no_neg_ctrl = Channel.value(file("$projectDir/assets/NO_FILE_neg_ctrl"))
        ch_no_pos_ctrl = Channel.value(file("$projectDir/assets/NO_FILE_pos_ctrl"))

        // get the list of meta.id from alignments
        // and assign it to the variable accepted_list
        accepted_list = alignments.map { it[0].id }.collect()
        accepted_list = accepted_list.flatten().toSortedList()
        // get just the BAM files from alignments
        ch_bams = alignments.map { [it[0], it[1]] }
        if (params.microbert){
            // ch_microbert_model = Channel.fromPath(params.microbert)
            // get the parent path of the model params.microbert
            ch_microbert_model = Channel.fromPath(params.microbert, checkIfExists: true)
            // get the parent path of ch_microbert_model
            // get the basename of the model
            ch_basename_microbert = ch_microbert_model.map { path -> path.getName() }
            // ch_microbert_model = ch_microbert_model.map { path -> path.parent }
            CLUSTER_ALIGNMENT(
                ch_bams
            )
            MMSEQS_EASYCLUSTER(
                CLUSTER_ALIGNMENT.out.fasta,
            )

            MICROBERT_PREDICT(
                MMSEQS_EASYCLUSTER.out.representatives.combine(ch_microbert_model)
            )
            ch_parse_files = MMSEQS_EASYCLUSTER.out.representatives.join(MICROBERT_PREDICT.out.predictions).join(MMSEQS_EASYCLUSTER.out.tsv)
            MICROBERT_PARSE(
                ch_parse_files.combine(ch_basename_microbert)
            )
            ch_report_microbert = MICROBERT_PARSE.out.report
        } else {
            ch_report_microbert = alignments.map { [ it[0], file("$projectDir/assets/NO_FILEmicrobert") ] }
        }
        alignments = alignments.join(ch_report_microbert)

        // Perform the difference operation
        missing_samples = all_samples - accepted_list
        missing_samples = missing_samples.flatten().toSortedList()
        if (!pathogens_list){
            println ("No pathogens list provided, skipping pathogen detection")
        } else{
            // ── Split alignments into control and non-control samples ──────────
            alignments.branch {
                control: it[0].control == true
                noncontrol: true
            }.set { split_alns }
            // alignments.view { println "All samples: ${it[0].id}, control: ${it[0].control}, negative: ${it[0].negative}, positive: ${it[0].positive}" }
            // split_alns.control.view { println "Control sample: ${it[0].id}" }
            // split_alns.noncontrol.view { println "Non-control sample: ${it[0].id}" }
            // ── Step 1: Run control samples FIRST ──────────────────────────────
            // Controls get NO_FILE placeholders for control JSON inputs and
            // their control_type is read from meta inside the process script.
            // errorStrategy 'ignore' is set so that control failures do not
            // block non-control samples from proceeding.
            ALIGNMENT_PER_SAMPLE_CONTROLS(
                split_alns.control
                    .combine(pathogens_list)
                    .combine(ch_sampletype_thresholds),
                assemblyfile,
                params.minmapq,
                ch_taxdump_dir,
                ch_no_neg_ctrl,
                ch_no_pos_ctrl,
            )

            // ── Step 2: Collect control JSON outputs into a value channel ──────
            // Build a single value [neg_jsons_list, pos_jsons_list] that the
            // non-control branch will .combine() with. Using toList() ensures
            // we wait for ALL controls to finish before non-controls start.
            control_json_map = ALIGNMENT_PER_SAMPLE_CONTROLS.out.txt
                .map { meta, json -> [meta.id, meta.control_type, json] }
                .toList()
                .map { entries ->
                    def neg_jsons = entries.findAll { it[1] == 'negative' }.collect { it[2] }
                    def pos_jsons = entries.findAll { it[1] == 'positive' }.collect { it[2] }
                    [neg_jsons, pos_jsons]
                }
                .ifEmpty { [[[], []]] }

            // ── Step 3: Run non-control samples with control JSONs ─────────────
            // For each non-control sample, resolve its specific negative/positive
            // control JSONs from the collected map based on meta.negative and
            // meta.positive sample name references.
            noncontrol_prepped = split_alns.noncontrol
                .combine(pathogens_list)
                .combine(ch_sampletype_thresholds)

            // Combine non-control alignment data with the collected control map
            // The control_json_map is a single-element value channel, so
            // .combine() will pair it with every non-control sample.
            noncontrol_with_ctrls = noncontrol_prepped
                .combine(control_json_map)

            // Now split into the tuple for ALIGNMENT_PER_SAMPLE and the
            // resolved neg/pos control JSON file channels.
            // The combined channel has shape:
            //   [meta, bam, bai, mapping, bedgraph, cov, k2, diamond, fastas,
            //    microbert, pathogens, thresholds, [neg_jsons], [pos_jsons]]
            // We need to extract the last two elements and resolve per-sample.

            noncontrol_tuple = noncontrol_with_ctrls.map { items ->
                // items[-1] = [pos_jsons], items[-2] = [neg_jsons]
                def all_neg = items[-2]
                def all_pos = items[-1]
                def meta = items[0]

                // Resolve this sample's specific negative control JSON
                def neg_json = null
                if (meta.negative && all_neg) {
                    neg_json = all_neg.find { it.name.startsWith(meta.negative.replaceAll(/\s+/, '_')) ||
                                               it.name.startsWith(meta.negative) ||
                                               it.name.contains(meta.negative.replaceAll(/\s+/, '_')) }
                }

                // Resolve this sample's specific positive control JSON
                def pos_json = null
                if (meta.positive && all_pos) {
                    pos_json = all_pos.find { it.name.startsWith(meta.positive.replaceAll(/\s+/, '_')) ||
                                               it.name.startsWith(meta.positive) ||
                                               it.name.contains(meta.positive.replaceAll(/\s+/, '_')) }
                }

                // Return: alignment tuple (first 12 items), neg_json, pos_json
                def aln_tuple = items[0..11]
                return aln_tuple + [neg_json, pos_json]
            }

            // Split the resolved channel into the process inputs
            noncontrol_aln_input = noncontrol_tuple.map { it[0..11] }
            noncontrol_neg_json = noncontrol_tuple.map { it[12] ?: file("$projectDir/assets/NO_FILE_neg_ctrl") }
            noncontrol_pos_json = noncontrol_tuple.map { it[13] ?: file("$projectDir/assets/NO_FILE_pos_ctrl") }
            ALIGNMENT_PER_SAMPLE(
                noncontrol_aln_input,
                assemblyfile,
                params.minmapq,
                ch_taxdump_dir,
                noncontrol_neg_json,
                noncontrol_pos_json,
            )

            // ── Step 4: Merge outputs from both control and non-control runs ───
            all_alignment_outputs = ALIGNMENT_PER_SAMPLE_CONTROLS.out.txt
                .mix(ALIGNMENT_PER_SAMPLE.out.txt)

            // collect all outputs into a single channel
            all_alignment_outputs.map{ m, txt -> txt }.collect().map{
                [[id: "all"], it]
            }.set{ full_list_pathogen_files }

            // merge
            SINGLE_REPORT(
                all_alignment_outputs.combine(distributions),
                false,
            )

            ORGANISM_MERGE_REPORT(
                full_list_pathogen_files.combine(distributions),
                missing_samples,
            )

            ch_template = Channel.fromPath("$projectDir/assets/heatmap.html", checkIfExists: true)

            CREATE_COMPARISON_REPORT(
                ORGANISM_MERGE_REPORT.out.report,
                ch_template
            )

            ch_pathogens_report = ORGANISM_MERGE_REPORT.out.report
        }
    emit:
        merged_report_txt = ch_pathogens_report
}
