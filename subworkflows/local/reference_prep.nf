/* groovylint-disable Indentation */
include { MAKE_FILE } from '../../modules/local/make_file'
include { DOWNLOAD_ASSEMBLY } from '../../modules/local/download_assembly'
include { MAP_LOCAL_ASSEMBLY_TO_FASTA } from '../../modules/local/map_assembly_to_fasta'
include { MAP_TAXID_ASSEMBLY } from '../../modules/local/map_taxid_assembly'
include { FEATURES_DOWNLOAD } from '../../modules/local/download_features'
include { FEATURES_TO_BED } from '../../modules/local/convert_features_to_bed'
include { COMBINE_MAPFILES } from '../../modules/local/combine_mapfiles'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_LOCAL } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_DWNLD } from '../../modules/nf-core/bowtie2/build/main'
include { HISAT2_BUILD as HISAT2_BUILD_LOCAL } from '../../modules/nf-core/hisat2/build/main'
include { HISAT2_BUILD as HISAT2_BUILD_DWNLD } from '../../modules/nf-core/hisat2/build/main'

workflow  REFERENCE_PREP {
    take:
    ch_samples
    ch_reference_fasta
    ch_assembly_txt
    ch_pathogens_file

    main:
    ch_versions = Channel.empty()
    ch_accessions = Channel.empty()
    ch_prepfiles = Channel.empty()

    ch_features = ch_samples.map{ meta, report -> {
            return [ meta,  []]
        }
    }

    ch_cds_to_taxids = ch_samples.map{ meta, report -> {
            return [ meta,  []]
        }
    }
    ch_fastas = ch_samples.map{ meta, report -> {
            return [ meta,  []]
        }
    }

    ch_cds = ch_samples.map{ meta, report -> {
            return [ meta,  []]
        }
    }

    ch_bedfiles = ch_samples.map { meta, report ->
        return [meta, []]
    }

    ch_mapped_assemblies = ch_samples.map{ meta, report -> {
            return [meta, [], [], []  ]
        }
    }

    ch_reports_to_download = ch_samples.map{ meta, report -> {
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


    if (!params.skip_realignment) {
        if (params.reference_fasta || params.get_pathogens) {
            MAP_LOCAL_ASSEMBLY_TO_FASTA(
                ch_reference_fasta
                .map { fasta ->
                    // Now fasta is an individual file, not a list
                    def basename = fasta.baseName
                    return [ [id: basename],  fasta]
                },
                ch_assembly_txt,
                ch_pathogens_file
            )
            // add ch_reference_fasta to all ch_fastas
            // add ch_reference_fasta to all ch_fastas
            ch_fastas = ch_fastas.combine(ch_reference_fasta.collect()).map{
                meta, fastas, fasta -> {
                    fastas.add(fasta)
                    return [meta, fastas]
                }
            }

            if (params.use_bt2) {
                if (params.bt2_indices) {
                    // If bt2_indices parameter is provided, create a channel from the provided path
                    ch_bt2_indices = params.bt2_indices ? Channel.from(params.bt2_indices.split(" ").collect { it  }) : Channel.empty()

                    fastaWithIndexChannel = ch_reference_fasta.merge(ch_bt2_indices)
                    ch_mapped_assemblies.combine(fastaWithIndexChannel.collect({
                        return [it[0], it[1]]
                    }, flat:false).toList()).map {
                        meta, fastas, listmaps, listids, fastaWithIndex ->
                            fastas.addAll([fastaWithIndex])
                            return [meta, fastas, listmaps, listids]
                    }.set { ch_mapped_assemblies }
                } else {
                    // If bt2_indices parameter is not provided, build the indices locally
                    // ////////////////////////////////////////////////////////////////////////////////

                    fastaForBowtieBuild = ch_reference_fasta
                        .map { fasta ->
                            def basen = fasta.baseName
                            return [ [id: basen], fasta ]
                        }

                    BOWTIE2_BUILD_LOCAL(fastaForBowtieBuild)

                    fastaWithIndexChannel = fastaForBowtieBuild
                        .join(BOWTIE2_BUILD_LOCAL.out.index, by: 0) // Join by the first element ('id')

                    ch_mapped_assemblies = ch_mapped_assemblies
                        .combine(fastaWithIndexChannel.collect({ return [it[1], it[2]] }, flat: false).toList())
                        .map { meta, fastas, listmaps, listids, fastaWithIndex ->
                            fastas.addAll([fastaWithIndex])
                            return [meta, fastas, listmaps, listids]
                        }

                    //////////////////////////////////////////////////////////////////////////////////
                }
            } else {
                // Case when `use_bt2` is false, just add the fastas directly
                ch_mapped_assemblies.combine(ch_reference_fasta.collect().toList()).map {
                    meta, fastas, listmaps, listids, fasta ->
                        fasta.each{ f -> {
                                return fastas.add([f])
                            }
                        }
                        return [meta, fastas, listmaps, listids]
                }.set { ch_mapped_assemblies }
            }
            // Collect the maps from `MAP_LOCAL_ASSEMBLY_TO_FASTA.out.map`
            MAP_LOCAL_ASSEMBLY_TO_FASTA.out.map
            .collect { meta, gcfmaps -> return gcfmaps }
            .set { merged_map }

            MAP_LOCAL_ASSEMBLY_TO_FASTA.out.accessions
            .collect { meta, accessions -> return accessions }
            .set { merged_map_ids }

            // Combine the maps into the third list (listmaps) of the ch_mapped_assemblies structure
            ch_mapped_assemblies.combine(merged_map.toList())
            .map{
                meta, fastas, listmaps, listids, map -> {
                    listmaps.addAll(map)
                    return [meta, fastas, listmaps, listids]
                }
            }.combine(
                merged_map_ids.toList()
            ).map{
                meta, fastas, listmaps, listids, ids -> {
                    listids.addAll(ids)
                    return [meta, fastas, listmaps, listids]
                }
            }.set{ ch_mapped_assemblies }
        }
    }
    // get the size of the ch_reports_to_download
    // if the size is greater than 0, then download the reports
    if ((!params.skip_refpull ) && (!params.skip_realignment && !(params.skip_kraken2 && (!params.organisms && !params.organisms_file)) ) ){
        DOWNLOAD_ASSEMBLY(
            ch_reports_to_download.map {
                meta, report ->  return [ meta, report ]
            },
            ch_assembly_txt
        )

        if (params.use_bt2) {
            BOWTIE2_BUILD_DWNLD(
                DOWNLOAD_ASSEMBLY.out.fasta
            )
            DOWNLOAD_ASSEMBLY.out.fasta.join(BOWTIE2_BUILD_DWNLD.out.index)
            .map{ meta, fasta, index -> [meta, [fasta, index]] }.set { merged_index }
        } else {
            DOWNLOAD_ASSEMBLY.out.fasta.map{meta, fasta -> [meta, [fasta]] }.set { merged_index }
        }
        // merge all the fasta outputs from the DOWNLOAD_ASSEMBLY process into ch_fastas on id meta
        ch_fastas = ch_fastas.join(DOWNLOAD_ASSEMBLY.out.fasta)
            .map { meta, fastas, fasta -> {
                fastas.add(fasta)
                return [meta, fastas]
            }
        }


        ch_mapped_assemblies.join(merged_index)
            .join(DOWNLOAD_ASSEMBLY.out.gcfids)
            .join(DOWNLOAD_ASSEMBLY.out.mapfile).set{ ch_mapped_assemblies }

        ch_mapped_assemblies.map { meta, fastas, listmaps, listids, fasta, gcfids, mapfile -> {
                listmaps.add(mapfile)
                listids.add(gcfids)
                fastas.add([fasta])
                return [meta, fastas, listmaps, listids ]
        }
        }.set{ ch_mapped_assemblies }
    }

    COMBINE_MAPFILES(
        ch_mapped_assemblies.map { meta, fastas, listmaps, listids ->  return [ meta, listmaps, listids ] }
    )

    ch_mapped_assemblies = ch_mapped_assemblies.join(COMBINE_MAPFILES.out.mergefiles)
        .map {
            meta, fastas, listmaps, listids, mergedmap, mergedids -> {
                return [ meta, fastas, mergedmap, mergedids ]
            }
        }

    MAP_TAXID_ASSEMBLY(
        ch_mapped_assemblies.map{meta, fastas, mergedmap, mergedids -> return [meta, mergedmap] },
        ch_assembly_txt
    )
    if ((params.use_diamond ) && ( !params.skip_refpull && ( !params.skip_realignment && !params.skip_features ) ) ) {
        try {
            // Attempt to use the FEATURES_DOWNLOAD process
            FEATURES_DOWNLOAD(
                ch_mapped_assemblies.map { meta, fastas, listmaps, listids ->
                    return [meta, listids]
                },
                ch_assembly_txt,
                true
            )
            ch_cds = FEATURES_DOWNLOAD.out.proteins
            ch_cds_to_taxids = FEATURES_DOWNLOAD.out.mapfile
            ch_versions = ch_versions.mix(FEATURES_DOWNLOAD.out.versions)
        /* groovylint-disable-next-line CatchException */
        } catch (Exception e) {
            // On failure, fallback to an alternative channel
            ch_cds = ch_samples.map { meta, report ->
                return [meta, []]
            }
            ch_cds_to_taxids = ch_samples.map { meta, report ->
                return [meta, []]
            }
        }
        try {
            FEATURES_TO_BED(
                FEATURES_DOWNLOAD.out.features
            )
            ch_features = FEATURES_DOWNLOAD.out.features
            ch_versions = ch_versions.mix(FEATURES_TO_BED.out.versions)
            ch_bedfiles = FEATURES_TO_BED.out.bed
        } catch (Exception e) {
            ch_bedfiles = ch_samples.map { meta, report ->
                return [meta, []]
            }
        }
    }
    ch_mapped_assemblies = MAP_TAXID_ASSEMBLY.out.taxidmerged.join(
        ch_mapped_assemblies.map{meta, fastas, mergedmap, mergedids -> return [meta, fastas, mergedids] }
    ).map{ meta, mergedmap, fastas, mergedids -> {
            return [meta, fastas, mergedmap, mergedids]
    }
    }
    emit:
        versions = ch_versions
        ch_bedfiles = ch_bedfiles
        ch_preppedfiles = ch_mapped_assemblies
        ch_reference_cds = ch_cds
        ch_cds_to_taxids = ch_cds_to_taxids
        features = ch_features
        fastas = ch_fastas
}
