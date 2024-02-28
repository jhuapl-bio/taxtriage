//
// Cross-check alignment abundances with pathogen list
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


include { PATHOGENS_FIND_SAMPLE } from '../../modules/local/pathogens_find'
include { PATHOGENS_MERGE_REPORT } from '../../modules/local/pathogens_merge'

workflow PATHOGENS {
    take:
        alignments
        pathogens_list
    main:
        PATHOGENS_FIND_SAMPLE(
            alignments,
            pathogens_list
        )
        // collect all outputs FIND_PATHOGENS.out.txt into a single channel
        // and assign it to the variable pathogens_list
        full_list_pathogen_files = PATHOGENS_FIND_SAMPLE.out.txt.map{m, txt -> txt}.collect()
        full_list_pathogen_files.view()
        PATHOGENS_MERGE_REPORT(
            full_list_pathogen_files
        )
        ch_pathogens_report = PATHOGENS_MERGE_REPORT.out.report

        println("Checking for pathogens")
    emit:
        pathogens_list = ch_pathogens_report
}
