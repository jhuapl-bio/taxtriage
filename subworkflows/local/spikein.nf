//
// SPIKEIN: fixed-background, varying-organism-load series.
//
// The BACKGROUND subworkflow answers "how deep must I sequence?" by subsampling a
// real FASTQ. SPIKEIN answers the other limit-of-detection question — "how much of
// this organism must be present?" — by holding the background at full depth and
// mixing in a defined number of reads from each spike organism's own reference.
//
// Flow:
//   spike-in sheet (csv/tsv/xlsx: accession, count, replicates[, level])
//     -> PARSE_SPIKEIN        normalise + group rows into levels
//     -> FETCH_SPIKEIN_REFS   one reference FASTA per accession
//     -> SPIKEIN_POOL         one simulated read pool per accession (ISS / NanoSim)
//     -> SPIKE_INTO_BACKGROUND  exact counts drawn per (level, replicate),
//                               concatenated onto the background
//
// Datasets are named <background>_background_ss_<mode>_c<level>_r<rep> — the same
// grammar the dilution series uses — so they flow through the existing insilico
// injection path and land in the In-Silico report tab. The manifest records
// kind=spikein so the report knows c<level> is a spike amount, not a depth.
//

include { PARSE_SPIKEIN          } from '../../modules/local/parse_spikein'
include { FETCH_SPIKEIN_REFS     } from '../../modules/local/fetch_spikein_refs'
include { SPIKEIN_POOL_ISS       } from '../../modules/local/spikein_pool'
include { SPIKEIN_POOL_NANOSIM   } from '../../modules/local/spikein_pool'
include { SPIKE_INTO_BACKGROUND  } from '../../modules/local/spike_into_background'


workflow SPIKEIN {
    take:
    ch_background        // tuple(meta, reads) — cleaned reads of each background source
    ch_assembly_summary  // path — assembly_summary_refseq.txt (or NO_FILE)

    main:
    ch_versions = Channel.empty()

    ch_sheet = Channel.value(file(params.spikein_sheet, checkIfExists: true))

    // ── 1. Normalise the sheet ──────────────────────────────────────────────
    PARSE_SPIKEIN(ch_sheet, params.sim_series_replicates ?: 1)
    ch_versions = ch_versions.mix(PARSE_SPIKEIN.out.versions)
    ch_spikein_tsv = PARSE_SPIKEIN.out.tsv

    // ── 2. One reference per distinct accession ─────────────────────────────
    ch_accessions = PARSE_SPIKEIN.out.accessions
        .splitText()
        .map { it.trim() }
        .filter { it }

    FETCH_SPIKEIN_REFS(ch_accessions.combine(ch_assembly_summary))
    ch_versions = ch_versions.mix(FETCH_SPIKEIN_REFS.out.versions.first())

    // ── 3. One simulated pool per accession ─────────────────────────────────
    // The pool must cover the largest single request for that accession, times a
    // factor so replicates draw different reads rather than the same set.
    ch_max_count = ch_spikein_tsv
        .splitCsv(header: true, sep: '\t')
        .map { row -> [row.accession, row.count as long] }
        .groupTuple()
        .map { acc, counts -> [acc, counts.max()] }

    def pool_factor = (params.spikein_pool_factor ?: 3) as int
    // Illumina pools are drawn as PAIRS, and `iss -n` counts single reads, so a
    // paired pool needs twice the records to yield the requested number of pairs.
    def bg_paired = !(params.background_reads2 == null && params.background_platform?.toString()?.toUpperCase() in ['OXFORD', 'PACBIO'])
    def sim_platform = (params.background_platform ?: (params.background_reads2 ? 'ILLUMINA' : 'OXFORD')).toString().toUpperCase()
    def is_ont = sim_platform in ['OXFORD', 'NANOPORE', 'ONT']

    ch_pool_input = FETCH_SPIKEIN_REFS.out.reference
        .join(ch_max_count)
        .map { acc, ref, maxc ->
            // `iss -n` counts single reads while a paired pool is drawn as PAIRS,
            // so a paired pool needs twice the records to yield the requested pairs.
            def n = Math.max(1000L, (maxc as long) * pool_factor * (is_ont ? 1 : 2))
            [acc, ref, n]
        }

    ch_nanosim_training = params.nanosim_training
        ? Channel.value(file(params.nanosim_training, checkIfExists: true))
        : Channel.value(file("$projectDir/assets/NO_FILE"))

    if (is_ont && !params.nanosim_training) {
        error "--spikein_sheet with an ONT background needs --nanosim_training so NanoSim can " +
              "simulate the spike reads. Provide it, or run the spike-in against an Illumina background."
    }

    // One simulator per platform, each in the same public biocontainer the
    // pipeline's existing ISS / NanoSim modules use.
    if (is_ont) {
        SPIKEIN_POOL_NANOSIM(ch_pool_input, ch_nanosim_training)
        ch_pools    = SPIKEIN_POOL_NANOSIM.out.pool
        ch_versions = ch_versions.mix(SPIKEIN_POOL_NANOSIM.out.versions.first())
    } else {
        SPIKEIN_POOL_ISS(ch_pool_input, params.iss_model ?: 'miseq')
        ch_pools    = SPIKEIN_POOL_ISS.out.pool
        ch_versions = ch_versions.mix(SPIKEIN_POOL_ISS.out.versions.first())
    }

    // Every pool must be staged into the single mixing task.
    ch_all_pools = ch_pools
        .map { acc, files -> files }
        .flatten()
        .collect()

    // ── 4. Mix into each background ─────────────────────────────────────────
    // Tag insilico-style so the datasets reuse the existing injection path; the
    // 'background' token in the id is what puts them in the In-Silico tab.
    ch_tagged_bg = ch_background.map { meta, reads ->
        def m = meta.collectEntries { k, v -> [k, v] }
        m.parent_id = meta.id
        m.id        = "${meta.id}_background"
        m.insilico  = true
        m.background = true
        m.spikein   = true
        m.control   = false
        m.trim      = false
        [m, reads]
    }

    def mode = params.sim_subsample_mode ?: 'randomized'
    def seed = params.spikein_seed ?: (params.sim_subsample_seed ?: 42)

    log.info "SPIKEIN: sheet=${params.spikein_sheet} mode=${mode} seed=${seed} platform=${sim_platform}"

    SPIKE_INTO_BACKGROUND(
        ch_tagged_bg,
        ch_spikein_tsv,
        ch_all_pools,
        mode,
        seed
    )
    ch_versions = ch_versions.mix(SPIKE_INTO_BACKGROUND.out.versions)

    // ── 5. Regroup the flat dataset FASTQs into one sample per dataset ──────
    ch_spikein_reads = SPIKE_INTO_BACKGROUND.out.reads
        .flatMap { meta, files ->
            def fl = (files instanceof List) ? files : [files]
            fl.collect { f ->
                def n = f.getName()
                def dsid = n.replaceAll(/\.spikein(_R[12])?\.fastq\.gz$/, '')
                tuple(dsid, meta, f)
            }
        }
        .groupTuple(by: 0)
        .map { dsid, metas, fastqs ->
            def m = metas[0].collectEntries { k, v -> [k, v] }
            m.id = dsid
            def cm = (dsid =~ /_c(\d+)_r\d+$/)
            if (cm.find()) {
                m.read_count = cm.group(1) as Integer
            }
            m.subsample      = true
            m.spikein        = true
            m.subsample_mode = params.sim_subsample_mode
            def sorted = fastqs.sort { it.getName() }
            [m, sorted.size() == 1 ? sorted[0] : sorted]
        }

    emit:
    spikein_reads = ch_spikein_reads                     // tuple(meta, reads) — new samples for ALIGNMENT
    manifests     = SPIKE_INTO_BACKGROUND.out.manifest   // *_spikein_manifest.tsv
    versions      = ch_versions
}
