//
// Check input samplesheet and get read channels
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

include { MINIMAP2_ALIGN as FILTER_MINIMAP2 } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX as FILTER_MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
include { KRAKEN2_KRAKEN2 as FILTER_KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main'
include { REMOVE_HOSTREADS } from '../../modules/local/remove_unaligned'
include { SAMTOOLS_INDEX as FILTERED_SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS as FILTERED_STATS } from '../../modules/nf-core/samtools/stats/main'
include { CHECK_GZIPPED_READS } from '../../modules/local/check_reads_exist'
include { DOWNLOAD_DB as FILTER_DB_DOWNLOAD } from '../../modules/local/download_db'
include { FETCH_HOST_REFS } from '../../modules/local/fetch_host_refs'
include { FETCH_HOST_REFS as DEACON_FETCH_HOST_REFS } from '../../modules/local/fetch_host_refs'
include { DEACON_INDEX_FETCH } from '../../modules/local/deacon_index_fetch'
include { DEACON_INDEX_BUILD } from '../../modules/local/deacon_index_build'
include { DEACON_FILTER } from '../../modules/local/deacon_filter'

// Deacon prebuilt indexes that `deacon index fetch` knows by name.
def deaconPrebuilt() { return ['panhuman-1', 'panmouse-1'] }

// Classify a --deacon_index / hosts.config `deacon_index` value: an .idx path
// (local or remote) is used as-is, anything else is a name to fetch.
def deaconIndexSource(String v) {
    if (v.endsWith('.idx') || v.contains('/')) {
        return ['path', v]
    }
    if (!deaconPrebuilt().contains(v)) {
        log.warn "Deacon index '${v}' is not a known prebuilt (${deaconPrebuilt().join(', ')}); trying 'deacon index fetch ${v}' anyway."
    }
    return ['fetch', v]
}

// Cache key for a Deacon index built from user FASTAs: name + a short hash of
// path/size/mtime, so a changed file under the same name is rebuilt.
def hostFastaKey(List paths) {
    def files = paths.collect { p -> file(p, checkIfExists: true) }
    def base  = files.size() == 1
        ? files[0].name.replaceAll(/\.(fa|fna|fasta|fas)(\.gz)?$/, '').replaceAll(/\.gz$/, '')
        : 'custom_host'
    def sig = files.collect { f -> "${f.toUriString()}:${f.size()}:${f.lastModified()}" }.join('|')
    return "${base}_${sig.md5().take(8)}".replaceAll(/[^A-Za-z0-9._-]/, '_')
}

workflow HOST_REMOVAL {
    take:
        ch_reads
        genome

    main:
        // // Remove the human reads first
        ch_bt2_index = Channel.empty()
        ch_filt_illumina = Channel.empty()
        ch_filt_oxfo = Channel.empty()
        ch_filtered_stats = Channel.empty()
        ch_host_removal_stats = Channel.empty()
        ch_reference_fasta = Channel.empty()
        ch_filter_db = Channel.empty()
        supported_filter_dbs = [
            'human': [
                'url': 'https://zenodo.org/records/8339700/files/k2_Human_20230629.tar.gz?download=1',
                'checksum': '9d388703b1fa7c2e269bb63acf1043dbec7bb62da0a57c4fb1c41d8ab7f9c953',
                'size': '5G'
            ],
            'viral': [
                'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz',
                'checksum': 'adf5deba8a62f995609592aa86e2f7aac7e49162e995e132a765b96edb456f99',
                'size': '553M'
            ],

        ]
        // minimap2 route (default). With --use_deacon the Deacon route below is
        // used instead and resolves its own index. Resolve the de-hosting
        // reference, in priority order:
        //   1. --remove_reference_file : a local FASTA the user supplied
        //   2. --genome <key> with a `fasta` path (the iGenomes entries)
        //   3. --genome <key> with `accessions` (the named host targets in
        //      conf/hosts.config, e.g. human / mosquito-any / tick-any) - the
        //      genomes are pulled from NCBI once and cached by FETCH_HOST_REFS.
        // Everything downstream consumes `ch_host_fasta`, a value channel, so the
        // fetched and the pre-existing routes behave identically.
        def genome_entry = (params.genome && params.genomes) ? params.genomes[params.genome] : null
        ch_host_fasta = Channel.empty()
        def run_reference_removal = false

        def ref_list = params.remove_reference_file
            ? params.remove_reference_file.toString().split(',').collect { it.trim() }.findAll { it }
            : []
        def run_deacon = false

        if (params.use_deacon) {
            // ── Deacon (minimizer-based) host depletion ────────────────────────
            // Index resolution, first match wins:
            //   1. --deacon_index <path.idx>          a local/remote prebuilt index
            //   2. --deacon_index <name>              a Deacon prebuilt (panhuman-1,
            //                                         panmouse-1), auto-fetched
            //   3. --remove_reference_file *.idx      reuse the existing flag for an index
            //   4. --remove_reference_file a.fa[,b.fa] build an index from local FASTA(s)
            //   5. --genome <key> with `deacon_index` (hosts.config, e.g. human ->
            //      panhuman-1), unless --deacon_prefer_prebuilt false
            //   6. --genome <key> with `fasta` (iGenomes)      -> build
            //   7. --genome <key> with `accessions`            -> fetch from NCBI, build
            //   8. nothing given                               -> panhuman-1
            // Fetched/built indexes are cached in --deacon_index_dir.
            run_deacon = true
            def mode = null
            def idx_val = null
            def entry_idx = genome_entry?.deacon_index ? genome_entry.deacon_index.toString() : null

            if (params.deacon_index) {
                def r_ = deaconIndexSource(params.deacon_index.toString()); mode = r_[0]; idx_val = r_[1]
            } else if (ref_list && ref_list.every { it.endsWith('.idx') }) {
                if (ref_list.size() > 1) {
                    error "--use_deacon accepts a single prebuilt index; got ${ref_list.size()} .idx files in --remove_reference_file. Build one combined index from the FASTAs instead."
                }
                def r_ = ['path', ref_list[0]]; mode = r_[0]; idx_val = r_[1]
            } else if (ref_list) {
                if (ref_list.any { it.endsWith('.idx') }) {
                    error "--remove_reference_file mixes .idx and FASTA files; pass either one Deacon index or FASTA file(s)."
                }
                def r_ = ['build_local', ref_list]; mode = r_[0]; idx_val = r_[1]
            } else if (entry_idx && params.deacon_prefer_prebuilt) {
                def r_ = deaconIndexSource(entry_idx); mode = r_[0]; idx_val = r_[1]
            } else if (genome_entry && genome_entry.fasta) {
                def r_ = ['build_local', [genome_entry.fasta.toString()]]; mode = r_[0]; idx_val = r_[1]
            } else if (genome_entry && genome_entry.accessions) {
                def r_ = ['build_host', genome_entry.accessions]; mode = r_[0]; idx_val = r_[1]
            } else {
                log.info "--use_deacon with no index, --remove_reference_file or --genome: defaulting to the prebuilt 'panhuman-1' index."
                def r_ = ['fetch', 'panhuman-1']; mode = r_[0]; idx_val = r_[1]
            }

            ch_deacon_index = Channel.empty()
            if (mode == 'path') {
                log.info "Deacon: using index ${idx_val}"
                ch_deacon_index = Channel.value(file(idx_val, checkIfExists: true))
            } else if (mode == 'fetch') {
                log.info "Deacon: fetching prebuilt index '${idx_val}' (≈3-4 GB, cached after the first run)"
                DEACON_INDEX_FETCH(Channel.of(idx_val))
                ch_deacon_index = DEACON_INDEX_FETCH.out.index.first()
            } else if (mode == 'build_local') {
                def key = (genome_entry && !ref_list) ? params.genome.toString() : hostFastaKey(idx_val)
                log.info "Deacon: building index '${key}' from ${idx_val.join(', ')}"
                DEACON_INDEX_BUILD(
                    Channel.of([ key, idx_val.collect { file(it, checkIfExists: true) } ])
                )
                ch_deacon_index = DEACON_INDEX_BUILD.out.index.first()
            } else if (mode == 'build_host') {
                log.info "Deacon: building index for host target '${params.genome}' from ${idx_val} (fetched from NCBI)"
                DEACON_FETCH_HOST_REFS(
                    Channel.of([ params.genome, idx_val, file(params.assembly ?: params.assembly_summary_refseq ?: "$projectDir/assets/NO_FILE") ])
                )
                DEACON_INDEX_BUILD(
                    DEACON_FETCH_HOST_REFS.out.fasta.map { target, fasta -> [ target, [ fasta ] ] }
                )
                ch_deacon_index = DEACON_INDEX_BUILD.out.index.first()
            }

            DEACON_FILTER(ch_reads, ch_deacon_index)
            ch_filtered_reads = DEACON_FILTER.out.reads
            ch_host_removal_stats = DEACON_FILTER.out.stats
                .filter{ !it[0].insilico }.collect{it[1]}.ifEmpty([])
        } else if (params.remove_reference_file){
            if (ref_list.any { it.endsWith('.idx') }) {
                error "--remove_reference_file points at a Deacon index (.idx); add --use_deacon, or give a FASTA for minimap2."
            }
            ch_host_fasta = Channel.value(file(params.remove_reference_file, checkIfExists: true))
            run_reference_removal = true
        } else if (genome_entry && genome_entry.fasta) {
            ch_host_fasta = Channel.value(file(genome_entry.fasta))
            run_reference_removal = true
        } else if (genome_entry && genome_entry.accessions) {
            def host_cache = params.host_reference_dir ?: "${params.outdir}/host_references"
            def host_label = genome_entry.description ?: params.genome
            println "Host target '${params.genome}' (${host_label}) will be de-hosted against ${genome_entry.accessions}; genomes cached in ${host_cache}"
            FETCH_HOST_REFS(
                Channel.of([ params.genome, genome_entry.accessions, file(params.assembly ?: params.assembly_summary_refseq ?: "$projectDir/assets/NO_FILE") ])
            )
            ch_host_fasta = FETCH_HOST_REFS.out.fasta.map { target, fasta -> fasta }.first()
            run_reference_removal = true
        }

        if (run_reference_removal){
            // Run minimap2 module on all LONGREAD platforms reads and the same on ILLUMINA reads
            // if ch_aligned_for_filter.shorteads is not empty
            // Run minimap2 on all for host removal - as host removal outperforms bowtie2 for host false negative rate https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9040843/

            FILTER_MINIMAP2(
                ch_reads.combine(ch_host_fasta),
                true,
                true,
                true
            )

            ch_bam_hosts = FILTER_MINIMAP2.out.bam

            // Join the BAM with the original (fastp-trimmed) reads so that
            // REMOVE_HOSTREADS can use QNAME-based FASTQ filtering.  This is
            // required for paired-end samples whose R1/R2 files may be
            // desynchronised (Casava 1.8+ headers, orphan reads after QC) and
            // where positional-pair BAM flags (-f 12) are therefore unreliable.
            REMOVE_HOSTREADS(
                ch_bam_hosts.join(ch_reads)
            )
            ch_filtered_reads = REMOVE_HOSTREADS.out.reads
            // Simulated (meta.insilico) datasets are excluded from the MultiQC feed.
            ch_host_removal_stats = REMOVE_HOSTREADS.out.stats
                .filter{ !it[0].insilico }.collect{it[1]}.ifEmpty([])

            // Continue processing the final reads
            FILTERED_SAMTOOLS_INDEX(
                ch_bam_hosts
            )

            ch_bai_files = ch_bam_hosts.join(FILTERED_SAMTOOLS_INDEX.out.bai)
            FILTERED_STATS(
                ch_bai_files,
                ch_host_fasta.map { fasta -> [ [], fasta ] }
            )
            // Simulated (meta.insilico) datasets are excluded from the MultiQC feed.
            ch_filtered_stats = FILTERED_STATS.out.stats
                .filter{ !it[0].insilico }.collect{it[1]}.ifEmpty([])
        } else if (params.filter_kraken2){
            if (supported_filter_dbs.containsKey(params.filter_kraken2)) {
                println "Kraken db ${params.filter_kraken2} will be downloaded if it cannot be found. This requires ${supported_filter_dbs[params.filter_kraken2]['size']} of space."
                FILTER_DB_DOWNLOAD(
                    params.filter_kraken2,
                    supported_filter_dbs[params.filter_kraken2]['url'],
                    supported_filter_dbs[params.filter_kraken2]['checksum']
                )
                /* groovylint-disable-next-line UnnecessaryGetter */
                ch_db = FILTER_DB_DOWNLOAD.out.k2d.map { file -> file.getParent() }
                ch_filter_db = ch_db
            } else {
                println "Kraken local filter db ${params.filter_kraken2} will be used."
                ch_filter_db = file(params.filter_kraken2)
            }
            FILTER_KRAKEN2 (
                ch_reads,
                ch_filter_db,
                true,
                false,
            )
            ch_reads = FILTER_KRAKEN2.out.unclassified_reads_fastq
        }
        // Shared by the minimap2 and Deacon routes: check the filtered output and
        // fall back to the original reads when de-hosting left nothing (or
        // removed the outputs because every read was host).
        if (run_reference_removal || run_deacon){
            CHECK_GZIPPED_READS(ch_filtered_reads, 4)
            ch_valid_reads = CHECK_GZIPPED_READS.out.check_result
            ch_orig_reads = ch_valid_reads.filter({
                it[1].name == 'emptyfile.txt'
            }).join(ch_reads).map({
                meta, result, reads -> return [meta, reads]
            })
            ch_filtered_reads = ch_valid_reads.filter({
                it[1].name == 'minimum_reads_check.txt'
            }).join(ch_filtered_reads).map({
                meta, result, reads -> return [meta, reads]
            })
            ch_reads = ch_orig_reads.mix(ch_filtered_reads)
        }

    emit:
        unclassified_reads = ch_reads
        stats_filtered = ch_filtered_stats
        host_removal_stats = ch_host_removal_stats
}
