include { MAKE_FILE } from '../../modules/local/make_file'
include { DOWNLOAD_ASSEMBLY } from '../../modules/local/download_assembly'
include { MAP_LOCAL_ASSEMBLY_TO_FASTA } from '../../modules/local/map_assembly_to_fasta'
include { FEATURES_TO_BED } from '../../modules/local/convert_features_to_bed'


workflow  REFERENCE_PREP {
    take:
    ch_samples
    ch_reference_fasta
    ch_assembly_txt


    main:
    ch_versions = Channel.empty()
    ch_accessions = Channel.empty()
    ch_mapped_assemblies = Channel.empty()

    ch_reports_to_download = ch_samples.map{ meta, reads, report -> {
            if (report){
                return [meta,  report  ]
            } else {
                return [meta, []]
            }

        }
    }
    if (params.organisms_file){
        // check if params.organisms is a file or a string
        ch_organisms = Channel.fromPath(params.organisms_file, checkIfExists: true)
        ch_reports_to_download = ch_reports_to_download.combine(ch_organisms).map{
            meta, report, organisms -> {
                report.add(organisms)
                return [meta, report]
            }
        }
    }
    if (params.organisms) {
        ch_organisms_taxids = Channel.from(params.organisms)
        // print params.organisms as a tsv, separated by space per
        MAKE_FILE(
            ch_organisms_taxids
        )
        ch_organisms = MAKE_FILE.out.file

        ch_reports_to_download = ch_reports_to_download.combine(
            ch_organisms
        ).map{
            meta, report, organisms -> {
                report.add(organisms)
                return [meta, report]
            }
        }
    }
    ch_prepfiles = ch_reports_to_download.map{
        meta, reports -> {
            return [meta, reports, [], [], [] ]
        }
    }
    if (params.reference_fasta || params.get_pathogens) { //
        // format of the FASTA file MUST be "kraken:taxid|<taxidnumber>" in each reference accession
        // merge ch_reference_fasta on all of the krakenreports. single channel merged to multiple
        MAP_LOCAL_ASSEMBLY_TO_FASTA(
            ch_reference_fasta.map {  fasta ->  {
                    // get basename of fasta path
                    def basen = fasta.getName()
                    return [ [id: basen ], fasta ]
                }
            },
            ch_assembly_txt
        )



        // MAP_LOCAL_ASSEMBLY_TO_FASTA.out.map.map {  meta, mapfile ->  return mapfile  }.set{ch_mapped_assemblies }

        // MAP_LOCAL_ASSEMBLY_TO_FASTA.out.accessions.map {  meta, accessions ->  return accessions  }.set{ch_accessions }

        ch_accessions.view()
        ch_mapped_assemblies.view()

        ch_prepfiles = ch_prepfiles.combine(
            MAP_LOCAL_ASSEMBLY_TO_FASTA.out.map.map {  meta, mapfile ->  return mapfile  }
        ).combine(
            MAP_LOCAL_ASSEMBLY_TO_FASTA.out.accessions.map {  meta, accessions ->  return accessions  }
        ).combine(ch_reference_fasta)

        ch_prepfiles.map { meta, reports, listfasta, listmaps, listids, accessions, ids,  fasta -> {
                listfasta.add(fasta)
                listmaps.add(accessions)
                listids.add(ids)
                return [meta, reports, listfasta, listmaps, listids]
            }
        }.set{ch_prepfiles}

        ch_prepfiles.view()
    }
    DOWNLOAD_ASSEMBLY(
        ch_organisms_to_download.map {
            meta, report ->  return [ meta, report ]
        },
        ch_assembly_txt
    )
    // ch_filtered_reads = ch_filtered_reads.join(DOWNLOAD_ASSEMBLY.out.fasta)

    // ch_accessions = DOWNLOAD_ASSEMBLY.out.accessions
    // ch_mapped_assemblies = DOWNLOAD_ASSEMBLY.out.mappings


    // if (params.get_features){

    //     FEATURES_DOWNLOAD(
    //         ch_accessions,
    //         ch_assembly_txt
    //     )

    //     FEATURES_TO_BED(
    //         FEATURES_DOWNLOAD.out.features
    //     )

    //     ch_bedfiles = FEATURES_TO_BED.out.bed
    //     fastq_reads  = fastq_reads.join(ch_bedfiles, remainder: true)
    // }
    // else {
    //     fastq_reads = fastq_reads.map {
    //         meta, reads, fasta -> [ meta, reads, fasta, null ]
    //     }
    // }

    // fastq_reads = fastq_reads.join(ch_mapped_assemblies)

    emit:
        versions = ch_versions
}
