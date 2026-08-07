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

include { SAMPLESHEET_CHECK   } from '../../modules/local/samplesheet_check'
include { GENERATE_SAMPLESHEET } from '../../modules/local/samplesheet_generate'
include { SRA_RESOLVE          } from '../../modules/local/sra_resolve'
include { SRA_FETCH_ENA        } from '../../modules/local/sra_fetch_ena'
include { SRA_FETCH_SRATOOLS   } from '../../modules/local/sra_fetch_sratools'

workflow INPUT_CHECK {


    main:
    ch_meta = file("$projectDir/assets/NO_FILE_meta_csv")
    versions = Channel.empty()


    if (params.fastq_1 || params.bam) {
        // Create a synthetic samplesheet from fastq/bam params.
        // params.fastq_1 may be a local path OR an SRA/ENA accession (SRR…, PRJNA…);
        // accessions are detected below and downloaded before anything else runs.
        GENERATE_SAMPLESHEET(
            [
                sampleName: params.sample ?: (is_accession(params.fastq_1) ? params.fastq_1.trim() : 'sample'),
                platform  : params.platform ?: '',
                fastq_1   : params.bam ? null : params.fastq_1,
                fastq_2   : params.bam ? null : params.fastq_2,
                bam       : params.bam,
                seq_summary: params.seq_summary,
                trim      : params.bam ? 'false' : params.trim,
                type      : params.type
            ]
        )

        ch_rows = GENERATE_SAMPLESHEET.out.csv.splitCsv(header: true, sep: ',')
        // GENERATE_SAMPLESHEET has no meta output; ch_meta stays as the
        // NO_FILE sentinel already set above when using --fastq_1 directly.
        versions = GENERATE_SAMPLESHEET.out.versions
    } else if (params.input) {
        // Use the provided samplesheet
        ch_rows = SAMPLESHEET_CHECK(file(params.input))
            .csv
            .splitCsv(header: true, sep: ',')

        SAMPLESHEET_CHECK.out.meta
            .set { ch_meta }
        versions = SAMPLESHEET_CHECK.out.versions
    }
    else {
        error "ERROR: Must specify one of --input, --fastq_1 or --bam"
    }

    // ── Split rows into local files vs SRA/ENA accessions ────────────────────
    // Rows whose fastq_1 looks like an accession (SRR/ERR/DRR, SRX, SRS/SAMN,
    // SRP/PRJNA…) are routed through the download path; everything else keeps
    // the existing local-path behaviour untouched.
    ch_rows
        .branch { row ->
            sra: is_sra_row(row)
            local: true
        }
        .set { ch_split }

    ch_local_reads = ch_split.local.map { create_fastq_channel(it) }

    // ── Resolve accessions → concrete runs ───────────────────────────────────
    // Every accession in the run is batched into ONE request file so the ENA /
    // NCBI lookup happens once, not once per sample. The file is headerless;
    // resolve_sra.py detects and skips a header if one is ever added.
    ch_requests = ch_split.sra
        .map { row -> "${row.sample},${row.fastq_1.trim()}" }
        .collectFile(name: 'sra_requests.csv', newLine: true, sort: true)

    SRA_RESOLVE(ch_requests)

    // Lookup of accession → original samplesheet row, so per-sample columns
    // (trim, type, controls, run metadata…) survive the round-trip.
    ch_lookup = ch_split.sra
        .map { row -> [ (row.fastq_1.trim()): row ] }
        .reduce([:]) { acc, entry -> acc + entry }

    ch_fetch_rows = SRA_RESOLVE.out.manifest
        .splitCsv(header: true, sep: ',')
        .combine(ch_lookup)
        .map { manifest_row, lookup -> merge_sra_row(manifest_row, lookup) }

    ch_fetch_rows
        .branch { row ->
            ena: row.sra_source == 'ena'
            sratools: true
        }
        .set { ch_fetch }

    SRA_FETCH_ENA(
        ch_fetch.ena.map { row ->
            [ row.sra_run, row.single_end, row.sra_url_1, row.sra_url_2, row.sra_md5_1, row.sra_md5_2 ]
        }
    )
    SRA_FETCH_SRATOOLS(
        ch_fetch.sratools.map { row -> [ row.sra_run, row.single_end ] }
    )

    // SRA processes only run when accessions were supplied, so these channels are
    // frequently empty — collect().ifEmpty([]) keeps the mix safe in that case.
    versions = versions
        .mix(SRA_RESOLVE.out.versions.collect().ifEmpty([]))
        .mix(SRA_FETCH_ENA.out.versions.collect().ifEmpty([]))
        .mix(SRA_FETCH_SRATOOLS.out.versions.collect().ifEmpty([]))

    // Rejoin the downloaded FASTQs with their metadata row. The join key is the
    // run accession, which resolve_sra.py guarantees is unique per manifest.
    ch_downloaded = SRA_FETCH_ENA.out.reads.mix(SRA_FETCH_SRATOOLS.out.reads)

    ch_sra_reads = ch_fetch_rows
        .map { row -> [ row.sra_run, row ] }
        .join(ch_downloaded)
        .map { run, row, files -> create_fastq_channel(row, files instanceof List ? files : [ files ]) }

    ch_local_reads
        .mix(ch_sra_reads)
        .toList()
        .flatMap { resolve_control_types(it) }
        .set { reads }

    emit:
        reads                                     // channel: [ val(meta), [ reads ] ]
        ch_meta
        versions
}

// ── SRA/ENA accession detection ─────────────────────────────────────────────
// Route a row down the download path if EITHER check_samplesheet.py flagged it
// (is_sra column) or the fastq_1 value matches an accession pattern here. The
// two checks are equivalent today; keeping both means a future divergence
// degrades to "downloaded anyway" rather than "treated as a missing file".
// The --fastq_1 param path has no is_sra column, so the regex is what fires there.
def is_sra_row(row) {
    if (row.containsKey('is_sra') && row.is_sra != null &&
        row.is_sra.toString().toLowerCase() == 'true') {
        return true
    }
    return is_accession(row.fastq_1)
}

// Delegates to WorkflowTaxtriage.isAccession (lib/WorkflowTaxtriage.groovy) so
// this subworkflow and the pre-flight check in workflows/taxtriage.nf can never
// disagree about what counts as an accession.
def is_accession(value) {
    return WorkflowTaxtriage.isAccession(value)
}

// Build the samplesheet row for one resolved SRA run by layering the manifest
// (run accession, platform, paired-ness, URLs) on top of the user's original
// row, so samplesheet columns the user set still win where they were provided.
def merge_sra_row(LinkedHashMap manifest_row, Map lookup) {
    def base = lookup[manifest_row.request] ?: [:]
    def row  = new LinkedHashMap(base)

    row.sample = manifest_row.sample
    // Honour an explicit platform column; otherwise use the instrument platform
    // reported by ENA/NCBI (ILLUMINA / OXFORD / PACBIO).
    def declared = base.platform?.toString()?.trim()
    row.platform = declared ? declared : (manifest_row.platform ?: 'ILLUMINA')

    row.single_end       = manifest_row.single_end
    row.directory        = 'false'
    row.is_fasta         = 'false'
    row.needscompressing = null

    row.sra_run    = manifest_row.run_accession
    row.sra_source = manifest_row.source
    row.sra_url_1  = manifest_row.fastq_1_url
    row.sra_url_2  = manifest_row.fastq_2_url
    row.sra_md5_1  = manifest_row.fastq_1_md5
    row.sra_md5_2  = manifest_row.fastq_2_md5

    // fastq_1/fastq_2 are filled in from the staged files after download.
    row.fastq_1 = null
    row.fastq_2 = null

    return row
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
//
// `resolved` is the list of already-staged Path objects for SRA-downloaded
// samples. When it is supplied the path/existence checks below are skipped —
// the files came straight out of a fetch process, so they exist by
// construction — and single/paired-ness is taken from what was actually
// downloaded rather than from the samplesheet.
def create_fastq_channel(LinkedHashMap row, List resolved = null) {
    // create meta map
    def meta = [:]
    // if fastq_2 is not a column then set it as null for all rows


    meta.id         = row.sample
    meta.platform = row.platform ? row.platform : 'ILLUMINA'
    // capitalize the platform
    meta.platform = meta.platform.toUpperCase()

    // ── Pre-aligned (BAM/CRAM) input ────────────────────────────────────────
    // A populated `bam` column means the sample is already aligned: every
    // upfront read step (compression, trimming, QC, host removal, classifier,
    // reference download) is bypassed and the file feeds coverage stats +
    // match_paths.py directly.  fastq_1/fastq_2 are ignored for these rows.
    // SRA-downloaded rows never carry a bam column, so the two paths are
    // mutually exclusive.
    meta.bam = (row.containsKey('bam') && row.bam) ? row.bam.trim() : null
    meta.is_bam = meta.bam ? true : false
    if (meta.is_bam && !file(meta.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> alignment (BAM/CRAM) file does not exist!\n${meta.bam}"
    }

    if (resolved != null) {
        // Sort so <run>_1 precedes <run>_2; a lone file sorts to position 0.
        def files = resolved.findAll { it != null }.sort { it.getName() }
        if (!files) {
            exit 1, "ERROR: no FASTQ files were downloaded for accession ${row.sra_run} (sample ${row.sample})"
        }
        meta.fastq_1 = files[0].toString()
        meta.fastq_2 = files.size() > 1 ? files[1].toString() : null
        meta.single_end = files.size() > 1 ? false : true
        meta.sra_run = row.sra_run
        meta.sra_source = row.sra_source
    } else {
        meta.fastq_1 = row.fastq_1
        // Check if 'fastq_2' exists in 'row'
        if (row.containsKey('fastq_2')) {
            meta.fastq_2 = row.fastq_2
        } else {
            meta.fastq_2 = null
        }
        // if meta.fastq_2 it is not single end, set meta.single_end as true else meta.single_end is false
        meta.single_end = (row.fastq_2 && !meta.is_bam) ? false : true
    }
    meta.needscompressing = row.needscompressing ? row.needscompressing : null

    // if meta.needscompressing is null or false AND the filename ends with .fastq or .fq then set to true
    // if (!meta.needscompressing && (meta.fastq_1.endsWith('.fastq') || meta.fastq_1.endsWith('.fq'))) {
    //     meta.needscompressing = true
    // }
    meta.aligner  = row.aligner ? row.aligner : 'minimap2'
    // if meta.aligner is not minimap2, hisat2, or bowtie2 then exit and send error
    if (meta.aligner != 'minimap2' && meta.aligner != 'hisat2' && meta.aligner != 'bowtie2') {
        exit 1, "ERROR: Please check input samplesheet -> aligner is not specified as minimap2, hisat2, or bowtie2 \n${meta.sample}"
    }
    // fastq_1 may be a ';'-delimited list of paths (multi-file input).
    // Check every path individually.  Skipped entirely for pre-aligned samples
    // and for SRA downloads, whose files come straight out of a fetch process.
    if (resolved == null && !meta.is_bam) {
        if (!meta.fastq_1) {
            exit 1, "ERROR: Please check input samplesheet -> a sample must have either fastq_1 or bam populated!\n${meta.id}"
        }
        meta.fastq_1.split(';').each { rawPath ->
            def p = rawPath.trim()
            if (p && !file(p).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> sequence file does not exist!\n${p}"
            }
        }
        if (!meta.single_end && meta.fastq_2 && !file(meta.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 sequence file (FASTQ) does not exist!\n${meta.fastq_2}"
        }
    }

    if (meta.is_bam) {
        meta.trim = false
    } else if (row.trim && row.trim.toString().toLowerCase() == "true"){
        meta.trim = true
    } else if (!row.trim  || (row.trim && row.trim.toString().toLowerCase() == "false")){
        meta.trim = false
    }
    meta.type = row.type
    meta.directory = (resolved == null && row.directory) ?  row.directory.toBoolean() : null
    meta.sequencing_summary = row.sequencing_summary ? file(row.sequencing_summary) : null
    // FASTA input flag: set by check_samplesheet.py when fastq_1 is a FASTA file.
    // Samples with is_fasta=true bypass all QC/trimming steps but still run
    // through host removal, Kraken2, and alignment.
    meta.is_fasta = (resolved == null && row.containsKey('is_fasta') && row.is_fasta != null && row.is_fasta.toString().toLowerCase() == 'true') ? true : false
    // Optional minimap2 preset override (e.g. map-ont, map-pb, map-hifi, lr:hq,
    // sr, splice, splice:hq, asm5, ava-pb, ava-ont …).  When set, overrides the
    // platform-derived default in the MINIMAP2_ALIGN module for BOTH host removal
    // and the main alignment step.  Leave blank to keep the existing behaviour.
    meta.minimap2_preset = (row.containsKey('minimap2_preset') && row.minimap2_preset) ? row.minimap2_preset.trim() : null

    // Control sample columns
    meta.control = (row.control && row.control.toString().toUpperCase() == "TRUE") ? true : false
    // Normalize spaces → underscores to match how check_samplesheet.py transforms sample names
    meta.negative = (row.containsKey('negative') && row.negative) ? row.negative.trim().replaceAll(/\s+/, '_') : null
    meta.positive = (row.containsKey('positive') && row.positive) ? row.positive.trim().replaceAll(/\s+/, '_') : null
    meta.control_type = null  // resolved later in resolve_control_types()

    // ── Run-level metadata (optional samplesheet columns) ────────────────────
    meta.run_id          = (row.containsKey('run_id')          && row.run_id)          ? row.run_id.trim()          : null
    meta.latitude        = (row.containsKey('latitude')        && row.latitude)        ? row.latitude.trim()        : null
    meta.longitude       = (row.containsKey('longitude')       && row.longitude)       ? row.longitude.trim()       : null
    meta.depth           = (row.containsKey('depth')           && row.depth)           ? row.depth.trim()           : null
    meta.salinity        = (row.containsKey('salinity')        && row.salinity)        ? row.salinity.trim()        : null
    meta.collection_time = (row.containsKey('collection_time') && row.collection_time) ? row.collection_time.trim() : null
    meta.location        = (row.containsKey('location')        && row.location)        ? row.location.trim()        : null

    // ── Build the [meta, reads] tuple ────────────────────────────────────────
    def fastq_meta = []

    if (resolved != null) {
        // SRA path: files are already staged Paths, in _1/_2 order.
        def files = resolved.findAll { it != null }.sort { it.getName() }
        fastq_meta = [ meta, files ]
        return fastq_meta
    }

    // fastq_1 may be a ';'-delimited list of file paths for multi-file inputs
    // (e.g. multiple FASTA assemblies fed to a single minimap2 splice run).
    if (meta.is_bam) {
        // Pre-aligned: the "reads" slot carries the alignment file so the tuple
        // shape stays [meta, [file]] for every consumer.
        meta.directory = false
        meta.is_fasta = false
        meta.needscompressing = null
        fastq_meta = [ meta, [ file(meta.bam) ] ]
    } else if (meta.directory) {
        if (meta.platform == 'OXFORD' || meta.platform == 'PACBIO') {
            fastq_meta = [ meta, [ file(meta.fastq_1) ] ]
        } else {
            exit 1, "ERROR: Please check input samplesheet -> the platform is not specified as OXFORD or PACBIO \n${meta.sample}"
        }
    } else {
        // Resolve all paths listed in fastq_1 (single path or ';'-separated list)
        def read1_files = meta.fastq_1.split(';').collect { it.trim() }.findAll { it }.collect { p ->
            if (!file(p).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> sequence file does not exist!\n${p}"
            }
            file(p)
        }

        if (meta.single_end || read1_files.size() > 1) {
            // Single-end OR multi-file: pass all read1 files as the reads list
            fastq_meta = [ meta, read1_files ]
        } else {
            // Paired-end single file: combine read1 + read2
            if (meta.fastq_2 && !file(meta.fastq_2).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> Read 2 sequence file (FASTQ) does not exist!\n${meta.fastq_2}"
            }
            fastq_meta = [ meta, [ read1_files[0], file(meta.fastq_2) ] ]
        }
    }

    return fastq_meta
}

// Second pass: infer which samples are controls by checking whether their
// sample name appears in any other sample's negative or positive column.
// Also honours an explicit control column if present in the samplesheet.
def resolve_control_types(List all_samples) {
    // Collect all sample names referenced as negative or positive controls.
    // Normalize spaces → underscores to match check_samplesheet.py which
    // converts spaces to underscores in the sample name column but not in
    // the negative/positive reference columns.
    def neg_names = all_samples.collect { it[0].negative }
        .findAll { it }
        .collect { it.replaceAll(/\s+/, '_') }
        .toSet()
    def pos_names = all_samples.collect { it[0].positive }
        .findAll { it }
        .collect { it.replaceAll(/\s+/, '_') }
        .toSet()

    all_samples.each { sample ->
        def name = sample[0].id
        // If this sample's name is referenced as a negative control by another sample
        if (name in neg_names) {
            sample[0].control = true
            sample[0].control_type = 'negative'
        }
        // If this sample's name is referenced as a positive control by another sample
        else if (name in pos_names) {
            sample[0].control = true
            sample[0].control_type = 'positive'
        }
        // Also honour explicit control column (if it was TRUE but not referenced)
        // — this sample is a control but we can't determine type from references
        // so leave control_type as null; it will still run in the controls branch
        // but without --control_type flag
    }
    println "Final sample metadata after resolving control types:"
    all_samples.each { sample ->
        println "\t${sample[0].id}: control=${sample[0].control}, control_type=${sample[0].control_type}, negative=${sample[0].negative}, positive=${sample[0].positive}"
    }
    return all_samples
}
