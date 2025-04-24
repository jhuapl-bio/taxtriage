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

include { ALIGNMENT_PER_SAMPLE } from '../../modules/local/alignment_per_sample'
include { ORGANISM_MERGE_REPORT } from '../../modules/local/report_merge'
include { ORGANISM_MERGE_REPORT as SINGLE_REPORT } from '../../modules/local/report_merge'

workflow REPORT {
    take:
        alignments
        pathogens_list
        distributions
        assemblyfile
        ch_taxdump_nodes
        all_samples
    main:
        ch_pathogens_report = Channel.empty()
        ch_pathognes_list = Channel.empty()
        // get the list of meta.id from alignments
        // and assign it to the variable accepted_list
        accepted_list = alignments.map { it[0].id }.collect()
        accepted_list = accepted_list.flatten().toSortedList()

        // Perform the difference operation
        missing_samples = all_samples - accepted_list
        missing_samples = missing_samples.flatten().toSortedList()
        if (!pathogens_list){
            println ("No pathogens list provided, skipping pathogen detection")
        } else{
            ALIGNMENT_PER_SAMPLE(
                alignments.combine(pathogens_list),
                assemblyfile
            )

            // collect all outputs FIND_PATHOGENS.out.txt into a single channel
            // and assign it to the variable pathogens_list
            ALIGNMENT_PER_SAMPLE.out.txt.map{ m, txt ->txt }.collect().map{
                [[id: "all"], it]
            }.set{ full_list_pathogen_files }

            SINGLE_REPORT(
                ALIGNMENT_PER_SAMPLE.out.txt.combine(distributions),
                false,
                ch_taxdump_nodes
            )

            ORGANISM_MERGE_REPORT(
                full_list_pathogen_files.combine(distributions),
                missing_samples,
                ch_taxdump_nodes
            )
            ch_pathogens_report = ORGANISM_MERGE_REPORT.out.report
        }
    emit:
        merged_report_txt = ch_pathogens_report
}
