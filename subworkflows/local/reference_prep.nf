include { MAKE_FILE } from '../../modules/local/make_file'
include { DOWNLOAD_ASSEMBLY } from '../../modules/local/download_assembly'
include { MAP_LOCAL_ASSEMBLY_TO_FASTA } from '../../modules/local/map_assembly_to_fasta'
include { FEATURES_DOWNLOAD } from '../../modules/local/download_features'
include { FEATURES_TO_BED } from '../../modules/local/convert_features_to_bed'
include { COMBINE_MAPFILES } from '../../modules/local/combine_mapfiles'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_LOCAL } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_DWNLD } from '../../modules/nf-core/bowtie2/build/main'

workflow  REFERENCE_PREP {
    take:
    ch_samples
    ch_reference_fasta
    ch_assembly_txt


    main:
    ch_versions = Channel.empty()
    ch_accessions = Channel.empty()
    ch_prepfiles = Channel.empty()
    ch_empty_file = file("$projectDir/assets/NO_FILE")


    ch_mapped_assemblies = ch_samples.map{ meta, report -> {
            return [meta, [], [], []  ]
        }
    }


    ch_bedfiles = Channel.empty()

    ch_reports_to_download = ch_samples.map{ meta, report -> {
            if (report){
                return [meta,  report  ]
            } else {
                return [meta, []]
            }

        }
    }
//     if (params.organisms_file){
//         // check if params.organisms is a file or a string
//         ch_organisms = Channel.fromPath(params.organisms_file, checkIfExists: true)
//         ch_reports_to_download = ch_reports_to_download.combine(ch_organisms).map{
//             meta, report, organisms -> {
//                 report.add(organisms)
//                 return [meta, report]
//             }
//         }
//     }
//     if (params.organisms) {
//         ch_organisms_taxids = Channel.from(params.organisms)
//         // print params.organisms as a tsv, separated by space per
//         MAKE_FILE(
//             ch_organisms_taxids
//         )
//         ch_organisms = MAKE_FILE.out.file

//         ch_reports_to_download = ch_reports_to_download.combine(
//             ch_organisms
//         ).map{
//             meta, report, organisms -> {
//                 report.add(organisms)
//                 return [meta, report]
//             }
//         }
//     }

//     if (params.reference_fasta || params.get_pathogens) { //
//         // format of the FASTA file MUST be "kraken:taxid|<taxidnumber>" in each reference accession
//         // merge ch_reference_fasta on all of the krakenreports. single channel merged to multiple
//         MAP_LOCAL_ASSEMBLY_TO_FASTA(
//             ch_reference_fasta.map {  fasta ->  {
//                     // get basename of fasta path
//                     def basen = fasta.getName()
//                     return [ [id: basen ], fasta ]
//             }
//             },
//             ch_assembly_txt
//         )

//         if (!params.skip_realignment) {
//             if (params.bt2_indices) {
//                 // If bt2_indices parameter is provided, create a channel from the provided path
//                 ch_bt2indices = Channel.fromPath(params.bt2_indices)
//                 ch_mapped_assemblies.map { meta, fastas, listmaps, listids -> {
//                         fastas.add([file(ch_reference_fasta), ch_bt2indices])
//                         return [meta, fastas, listmaps, listids]
//                     }
//                 }
//             } else {
//                 illuminaPresent = ch_samples
//                       .filter { it[0].platform == "ILLUMINA" }
//                       .count()
//                       .map { it > 0 }
//                 ch_bt2_indices = Channel.empty()
//                 // Then, depending on the presence of 'ILLUMINA' samples, optionally run BOWTIE2_BUILD
//                 if (illuminaPresent) {
//                     // If there are 'ILLUMINA' samples, we prepare and run BOWTIE2_BUILD
//                     println("ILLUMINA samples found, performing BOWTIE2_BUILD: Local.")
//                     ch_reference_fasta
//                         .map { fasta ->
//                             def basen = fasta.baseName
//                             return [ [id: basen], fasta ]
//                         }
//                         .set { fastaForBowtieBuild }

//                     BOWTIE2_BUILD_LOCAL(fastaForBowtieBuild)


//                     fastaForBowtieBuild
//                         .join(BOWTIE2_BUILD_LOCAL.out.index, by: 0) // Join by the first element ('id')
//                         .set { fastaWithIndexChannel }
//                     ch_mapped_assemblies.combine(fastaWithIndexChannel.map{
//                         meta, fasta, index -> {
//                             return [fasta, index]
//                         }
//                     }).map{
//                         meta, fastas, listmaps, listids, singlefasta, fastaWithIndex -> {
//                             fastas.add([singlefasta, fastaWithIndex]) // Add FASTA and Index to 'fastas'
//                             return [meta, fastas, listmaps, listids]
//                         }
//                     }.set { ch_mapped_assemblies }
//                 } else {
//                     println("No ILLUMINA samples found, skipping BOWTIE2_BUILD: Local.")
//                 }
//             }
//         }

//         ch_mapped_assemblies = ch_mapped_assemblies.combine(
//             MAP_LOCAL_ASSEMBLY_TO_FASTA.out.map.map {  meta, mapfile ->  return mapfile  }
//         ).combine(
//             MAP_LOCAL_ASSEMBLY_TO_FASTA.out.accessions.map {  meta, gcfids ->  return gcfids  }
//         )

//         ch_mapped_assemblies.map { meta, fastas, listmaps, listids, mapfile, gcfids -> {
//                 listmaps.add(mapfile)
//                 listids.add(gcfids)
//                 return [meta, fastas, listmaps, listids]
//         }
//         }.set{ ch_mapped_assemblies }
//     }
//     if (!params.skip_kraken2 && !params.skip_realignment){
//         DOWNLOAD_ASSEMBLY(
//             ch_reports_to_download.map {
//                 meta, report ->  return [ meta, report ]
//             },
//             ch_assembly_txt
//         )

//         // get all where meta.platform == ILLUMINA and run BOWTIE2_BUILD on the fasta
//         // get all where meta.platform == OXFORD and run MINIMAP2_BUILD on the fasta
//         DOWNLOAD_ASSEMBLY.out.fasta.branch{
//                 longreads: it[0].platform =~ 'OXFORD'
//                 shortreads: it[0].platform =~ 'ILLUMINA'
//         }.set { ch_platform_split }

//         BOWTIE2_BUILD_DWNLD(
//             ch_platform_split.shortreads
//         )
//         ch_platform_split.shortreads.join(BOWTIE2_BUILD_DWNLD.out.index)
//         .map{meta, fasta, index -> [meta, [fasta, index]] }.set { merged_shortreads_index }

//         ch_platform_split.longreads.map{meta, fasta -> [meta, [fasta]] }.set { merged_longreads_only }

//         ch_fullset = merged_shortreads_index.mix(merged_longreads_only)



//         ch_mapped_assemblies.join(ch_fullset)
//             .join(DOWNLOAD_ASSEMBLY.out.gcfids)
//             .join(DOWNLOAD_ASSEMBLY.out.mapfile).set{ ch_mapped_assemblies }


//         ch_mapped_assemblies.map { meta, fastas, listmaps, listids, fasta, gcfids, mapfile -> {
//                 listmaps.add(mapfile)
//                 listids.add(gcfids)
//                 fastas.add(fasta)
//                 return [meta, fastas, listmaps, listids ]
//         }
//         }.set{ ch_mapped_assemblies }
//     }

//     COMBINE_MAPFILES(
//         ch_mapped_assemblies.map { meta, fastas, listmaps, listids ->  return [ meta, listmaps, listids ] }
//     )
//     ch_mapped_assemblies = ch_mapped_assemblies.join(COMBINE_MAPFILES.out.mergefiles)
//         .map {
//             meta, fastas, listmaps, listids, mergedmap, mergedids -> {
//                 return [ meta, fastas, mergedmap, mergedids ]
//             }
//         }

//     if (params.get_features){
//         FEATURES_DOWNLOAD(
//             ch_mapped_assemblies.map { meta, fastas, listmaps, listids  ->  return [ meta, listids ] },
//             ch_assembly_txt
//         )

//         FEATURES_TO_BED(
//             FEATURES_DOWNLOAD.out.features
//         )

//         ch_bedfiles = FEATURES_TO_BED.out.bed
//     }

    emit:
        versions = ch_versions
        ch_feature_bedfiles = ch_bedfiles
        ch_preppedfiles = ch_mapped_assemblies
}
