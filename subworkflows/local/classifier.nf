//
// Run upstream classifier using Kraken2 and, in future, metaphlan or centrifuge
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


include { METAPHLAN_METAPHLAN } from '../../modules/nf-core/metaphlan/metaphlan/main'
include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { TOP_HITS } from '../../modules/local/top_hits'
include { REMOVETAXIDSCLASSIFICATION } from '../../modules/local/remove_taxids.nf'
include { KRAKENREPORT } from '../../modules/local/krakenreport'
include { TAXPASTA_STANDARDISE } from '../../modules/nf-core/taxpasta/standardise/main'
include { TAXPASTA_MERGE } from '../../modules/nf-core/taxpasta/merge/main'
include { MERGEDKRAKENREPORT } from '../../modules/local/merged_krakenreport'
include { FILTERKRAKEN } from '../../modules/local/filter_krakenreport'
include { KREPORT_TO_KRONATXT } from '../../modules/local/generate_krona_txtfile'
include { KREPORT_TO_TAXONOMY } from '../../modules/local/kreport_to_taxonomy'
include { KRONA   } from '../../modules/local/krona.nf'
include { EXTRACT_TOP_SEQS   } from '../../modules/local/extract_top_seqs.nf'
include { KRONA_KTIMPORTTEXT  } from '../../modules/nf-core/krona/ktimporttext/main'
include { KRAKENTOOLS_COMBINEKREPORTS   } from '../../modules/nf-core/krakentools/combinekreports/main'
include { MERGEDSUBSPECIES } from '../../modules/local/merged_subspecies'

//include { CENTRIFUGE_CENTRIFUGE } from '../../modules/nf-core/centrifuge/centrifuge/main'
//include { CENTRIFUGE_KREPORT } from '../../modules/nf-core/centrifuge/kreport/main'

workflow CLASSIFIER {
    take:
        ch_reads
        ch_db
        ch_save_fastq_classified
        distributions
        ch_pathogens
        ch_organisms_to_download
        ch_taxdump_dir

    main:
        ch_kraken2_report = Channel.empty()
        ch_metaphlan_report = Channel.empty()
        ch_tops = Channel.empty()
        ch_krona_plot = Channel.empty()
        ch_empty_file = file("$projectDir/assets/NO_FILE")
        if (!params.skip_kraken2){
            // // // // // //
            // // // // // // MODULE: Run Kraken2
            // // // // // //

            // // // // // // //
            // // // // // // // MODULE: Run Kraken2
            // // // // // // //
            KRAKEN2_KRAKEN2(
                ch_reads,
                ch_db,
                ch_save_fastq_classified,
                false
            )

            ch_kraken2_report = KRAKEN2_KRAKEN2.out.report

            // KREPORT_TO_TAXONOMY(
            //     ch_kraken2_report
            // )
            // ch_taxnames = KREPORT_TO_TAXONOMY.out.names
            // ch_taxnodes = KREPORT_TO_TAXONOMY.out.nodes

            KREPORT_TO_KRONATXT(
                ch_kraken2_report
            )

            ch_krona_txt = KREPORT_TO_KRONATXT.out.txt

            ch_combined = ch_krona_txt
                        .map{ it[1] }        // Get the file path
                        .collect()            // Collect all file parts into a list
                        .map { files ->
                            // Join the files with single quotes and space
                            // String joinedFiles = files.collect { "'$it'" }.join(' ')
                            // if single file then make it [files] otherwise just files
                            [[id:'combined_krona_kreports'], files instanceof List ? files : [files]]  // Combine with new ID
                        }
            KRONA_KTIMPORTTEXT(
                ch_combined
            )
            ch_krona_plot = KRONA_KTIMPORTTEXT.out.html

            if (params.remove_taxids) {
                remove_input = ch_kraken2_report.map {
                    meta, report -> [
                        meta, report, params.remove_taxids
                    ]
                }
                REMOVETAXIDSCLASSIFICATION(
                    remove_input
                )
                ch_kraken2_report = REMOVETAXIDSCLASSIFICATION.out.report
            }

            TOP_HITS(
                ch_kraken2_report.combine(distributions).combine(ch_pathogens)
            )
            ch_tops = TOP_HITS.out.tops
            MERGEDSUBSPECIES(
                ch_kraken2_report.map{
                    meta, report -> report
                }.collect(),
                ch_pathogens
            )

            MERGEDKRAKENREPORT(
                TOP_HITS.out.krakenreport.map { meta, file ->  file }.collect()
            )
            FILTERKRAKEN(
                MERGEDKRAKENREPORT.out.krakenreport
            )

            if (ch_save_fastq_classified){
                ch_reads = KRAKEN2_KRAKEN2.out.classified_reads_fastq.map { m, r-> [m, r.findAll { it =~ /.*\.classified.*(fq|fastq)(\.gz)?/  }] }
             
            }

            if (params.fuzzy){
                ch_organisms = TOP_HITS.out.names
            } else {
                ch_organisms = TOP_HITS.out.taxids
            }
            // mix ch_organisms_to_download with ch_organisms 2nd index list
            ch_organisms_to_download = ch_organisms_to_download.join(
                ch_organisms
            ).map{
                meta, report, organisms -> {
                    report.add(organisms)
                    return [meta, report]
                }
            }
        } else {
            // set ch_kraken2_report to meta, null
            ch_kraken2_report = ch_reads.map{ meta, reads -> {
                    return [ meta,  ch_empty_file]
                }
            }
        }
        if (params.metaphlan) {
            METAPHLAN_METAPHLAN(
                ch_reads,
                params.metaphlan
            )
            ch_metaphlan_report = METAPHLAN_METAPHLAN.out.profile.map{ meta, file -> {
                    return [ meta, file, 'metaphlan' ]
                }
            }

            TAXPASTA_STANDARDISE(
                ch_metaphlan_report,
                ch_taxdump_dir
            )
            ch_standardized = TAXPASTA_STANDARDISE.out.standardised_profile
        }
    emit:
        ch_kraken2_report
        ch_metaphlan_report
        ch_reads
        ch_organisms_to_download
        ch_tops
        ch_krona_plot
}
