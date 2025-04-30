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
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_INDEX as MINIMAP2_INDEX_LOCAL } from '../../modules/nf-core/minimap2/index/main'

workflow  REFERENCE_PREP {
    take:
    ch_samples
    ch_reference_fasta
    ch_assembly_txt
    ch_pathogens_file

    main:
    ch_versions = Channel.empty()
    ch_accessions = Channel.empty()
    ch_features = Channel.empty()
    ch_prepfiles = Channel.empty()
    ch_fastas = Channel.empty()

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
        ch_versions = ch_versions.mix(MAKE_FILE.out.versions)
        ch_reports_to_download = ch_reports_to_download.combine(
            ch_organisms
        ).map{
            meta, report, organisms -> {
                report.add(organisms)
                return [meta, report]
            }
        }
    }
    if (!params.skip_realignment || params.force_pull) {
        // get the size of the ch_reports_to_download
        // if the size is greater than 0, then download the reports
        if ((!params.skip_kraken2 && params.force_pull)|| (!params.skip_refpull) && (!params.skip_realignment && !(params.skip_kraken2 && (!params.organisms && !params.organisms_file)) ) ){
            DOWNLOAD_ASSEMBLY(
                ch_reports_to_download.map {
                    meta, report ->  return [ meta, report ]
                },
                ch_assembly_txt
            )
            ch_versions = ch_versions.mix(DOWNLOAD_ASSEMBLY.out.versions)

            DOWNLOAD_ASSEMBLY.out.fasta.map{meta, fasta -> [meta, [fasta]] }.set { merged_index }

            if (params.use_bt2){
                BOWTIE2_BUILD(
                    merged_index
                )
                // join BOWTIE2_BUILD.out.index with mergeD_index
                ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
                merged_index.join(BOWTIE2_BUILD.out.index).set{ merged_index }

            } else if (!params.use_hisat2){
                MINIMAP2_INDEX(
                    merged_index
                )
                ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
                merged_index.join(MINIMAP2_INDEX.out.index).set{ merged_index }
            } else {
                merged_index = DOWNLOAD_ASSEMBLY.out.fasta.map{meta, fasta -> [meta, [fasta], null] }
            }
            merged_index = merged_index.map{
                meta, fastas, indices -> {
                    def ff = fastas
                    if (indices){
                        ff.add(indices)
                    }
                    return [meta, ff]
                }
            }

            ch_versions = ch_versions.mix(DOWNLOAD_ASSEMBLY.out.versions)


            ch_mapped_assemblies
                .join(merged_index)
                .join(DOWNLOAD_ASSEMBLY.out.mapfile)
                .join(DOWNLOAD_ASSEMBLY.out.gcfids)
                .set{ ch_mapped_assemblies }


            ch_mapped_assemblies.map { meta, fastas, listmaps, listids, fastadwnl, mapfile, gcfids -> {
                    listmaps.add(mapfile)
                    listids.add(gcfids)
                    fastas.add(fastadwnl)
                    return [meta, fastas, listmaps, listids ]
            }
            }.set{ ch_mapped_assemblies }

        }
        if (params.reference_fasta || params.get_pathogens) {
            dff = ch_reference_fasta
                .map { fasta ->
                    // Now fasta is an individual file, not a list
                    def basename = fasta.baseName
                    return [ [id: basename],  fasta]
                }.combine(ch_pathogens_file).combine(ch_assembly_txt)

            MAP_LOCAL_ASSEMBLY_TO_FASTA(
                dff
            )
            // add ch_reference_fasta to all ch_fastas
            ch_local_fastas = ch_reference_fasta
            // Initialize an empty channel for seen indices

            // Create the seen indices channel based on parameters
            if ( params.use_bt2 && params.bt2_indices ) {
                println "bt2 indices being used"
                ch_seen_indices = params.bt2_indices.split(" ")
            } else if ( params.mmi ) {
                println "mmi indices being used"
                ch_seen_indices = params.mmi.split(" ")
            } else {
                // no indices provided; create an empty list to avoid problems later
                ch_seen_indices = []
            }
            // for each ch_seen_indices make it a file
            ch_seen_indices = ch_seen_indices.collect { f ->
                return file(f)
            }
            ch_reference_fasta.map{m->[[m]]}.collect().map{
                fasta -> {
                    def idx=0;
                    def size_of_indices = ch_seen_indices.size()
                    def fastas=[]
                    fasta.each{ f->
                        // println "fasta: ${f}, ${size_of_indices}"
                        if (size_of_indices > idx){
                            // if the size of the indices is less than the index, then add the fasta to the fastas list
                            f.add(ch_seen_indices[idx])
                        }
                        fastas.add(f)
                        idx++
                    }
                    return fastas
                }
            }.flatMap { innerList -> innerList }.set{ ch_fastas_tmp }

            ch_fastas_tmp.map{
                fasta -> {
                    def basename = fasta.first().baseName
                    return [ [id: basename],  fasta]
                }
            }.set{ ch_fastas_tmp }

            ch_fastas_tmp.branch{
                bt2_index: params.use_bt2 && it[1].size() <= 1
                mmap2_index: !params.use_hisat2 && it[1].size() <= 1
                other: true
            }
            .set{ ch_brahc}

            MINIMAP2_INDEX_LOCAL(
                ch_brahc.mmap2_index
            )
            BOWTIE2_BUILD_LOCAL(
                ch_brahc.bt2_index
            )

            ch_brahc.bt2_index
                .join(BOWTIE2_BUILD_LOCAL.out.index)
                .map{
                    meta, file, index -> {
                        return [meta, [file.first(), index]]
                    }
                }
                .set{ indx }

            ch_brahc.mmap2_index
                .join(MINIMAP2_INDEX_LOCAL.out.index)
                .map{
                    meta, file, index -> {
                        return [meta, [file.first(), index]]
                    }
                }
                .set{ indx2 }


            ch_indexes = indx2.mix(indx).mix(ch_brahc.other).map{
                _, fasta -> { return fasta}
            }.collect(flat:false).map{fasta->return [fasta]}

            // // Collect the maps from `MAP_LOCAL_ASSEMBLY_TO_FASTA.out.map`
            MAP_LOCAL_ASSEMBLY_TO_FASTA.out.map
                .collect { meta, gcfmaps -> return gcfmaps }
                .set { merged_map }
            ch_versions = ch_versions.mix(MAP_LOCAL_ASSEMBLY_TO_FASTA.out.versions)

            MAP_LOCAL_ASSEMBLY_TO_FASTA.out.accessions
                .collect { meta, accessions -> return accessions }
                .set { merged_map_ids }
            ch_versions = ch_versions.mix(MAP_LOCAL_ASSEMBLY_TO_FASTA.out.versions)

            // combine merged_map_ids, merged_map and ch_indexes together
            ch_indexes.combine(
                merged_map.map{gcfmaps->return [gcfmaps]}
            ).combine(
                merged_map_ids.map{gcfids->return [gcfids]}
            ).set{ ch_indexes }

            ch_mapped_assemblies.combine(ch_indexes).map{
                meta, fastas, listmaps, listids, fasta_local, maps_local, gcfids_local -> {
                    fastas.addAll(fasta_local)
                    listmaps.addAll(maps_local)
                    listids.addAll(gcfids_local)
                    return [meta, fastas, listmaps, listids ]
                }
            }.set{ ch_mapped_assemblies }



        }
    }


    COMBINE_MAPFILES(
        ch_mapped_assemblies.map { meta, fastas, listmaps, listids ->  return [ meta, listmaps, listids ] }
    )
    ch_fastas = ch_mapped_assemblies.map{
        meta, fastas, listmaps, listids -> {
            def fasta_list = []
            fastas.each{
                f ->{
                    // if a list add f else add first element
                    if (f instanceof List){
                        fasta_list.add(f.first())
                    } else {
                        fasta_list.add(f)
                    }
                }
            }
            return [meta, fasta_list]
        }
    }
    ch_versions = ch_versions.mix(COMBINE_MAPFILES.out.versions)

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
    ch_versions = ch_versions.mix(MAP_TAXID_ASSEMBLY.out.versions)

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
        versions = ch_versions
        features = ch_features
}
