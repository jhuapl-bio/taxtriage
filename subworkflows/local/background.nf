//
// BACKGROUND: dilution-series subsampling of a real (natural) FASTQ.
//
// Mirrors the subsampling stage of INSILICO, but the master pool is a
// user-provided FASTQ (single or paired) rather than synthetic ISS/NanoSim
// reads. The full-depth background sample is classified + reference-prepped once
// upstream (it enters the pipeline as a normal sample, tagged background_source);
// here we take its cleaned reads and subsample them into the same spike-in /
// dilution series. Each derived dataset is tagged insilico-style so it reuses the
// existing injection path (clone parent refs -> ALIGNMENT_PER_SAMPLE_INSILICO)
// and is picked up by the In-Silico suite report tab (platform token 'background').
//

include { SUBSAMPLE_INSILICO } from '../../modules/local/subsample_insilico'


// Resolve the read-count series from an explicit list or the start/step/n generator.
def bgResolveSeriesCounts() {
    if (params.sim_series_counts) {
        return params.sim_series_counts.toString()
            .split(/[,\s]+/).findAll { it }.collect { it as long }
    }
    if (params.sim_series_start != null && params.sim_series_step != null && params.sim_series_n) {
        def start = params.sim_series_start as long
        def step  = params.sim_series_step as long
        def n     = params.sim_series_n as int
        return (0..<n).collect { start + it * step }
    }
    error "background_reads is set but no series was defined. Set --sim_series_counts " +
          "(e.g. '100,500,1000') OR --sim_series_start/--sim_series_step/--sim_series_n."
}


// Regroup the flat per-dataset FASTQs into one (meta, reads) tuple per dataset.
def bgRegroup(ch_subsample_reads) {
    ch_subsample_reads
        .flatMap { meta, files ->
            def fl = (files instanceof List) ? files : [files]
            fl.collect { f ->
                def n = f.getName()
                def dsid = n.replaceAll(/\.subsample(_R[12])?\.fastq\.gz$/, '')
                tuple(dsid, meta, f)
            }
        }
        .groupTuple(by: 0)
        .map { dsid, metas, fastqs ->
            def pmeta = metas[0]
            def m = pmeta.collectEntries { k, v -> [k, v] }
            m.id = dsid
            def cm = (dsid =~ /_c(\d+)_r\d+$/)
            if (cm.find()) {
                m.read_count = cm.group(1) as Integer
            }
            m.subsample      = true
            m.background     = true
            m.subsample_mode = params.sim_subsample_mode
            def sorted = fastqs.sort { it.getName() }
            [m, sorted.size() == 1 ? sorted[0] : sorted]
        }
}


workflow BACKGROUND {
    take:
    ch_master   // tuple(meta, reads) — the classified background sample's cleaned reads
                //   (meta carries background_source=true, id=<background_name>)

    main:
    ch_versions  = Channel.empty()
    ch_manifests = Channel.empty()

    // Tag the master insilico-style so its subsamples reuse the insilico injection
    // path. Platform token 'background' distinguishes it in the id / report tab.
    ch_tagged_master = ch_master.map { meta, reads ->
        def m = meta.collectEntries { k, v -> [k, v] }
        m.parent_id  = meta.id
        m.id         = "${meta.id}_background"
        m.insilico   = true
        m.background = true
        m.control    = false
        m.trim       = false
        [m, reads]
    }

    def counts = bgResolveSeriesCounts().join(',')
    def reps   = params.sim_series_replicates ?: 1
    def mode   = params.sim_subsample_mode ?: 'randomized'
    def seed   = params.sim_subsample_seed ?: 42

    log.info "BACKGROUND subsampling: mode=${mode} counts=[${counts}] replicates=${reps} seed=${seed}"

    SUBSAMPLE_INSILICO(
        ch_tagged_master,
        mode,
        counts,
        reps,
        seed
    )
    ch_versions  = ch_versions.mix(SUBSAMPLE_INSILICO.out.versions)
    ch_manifests = SUBSAMPLE_INSILICO.out.manifest

    ch_background_reads = bgRegroup(SUBSAMPLE_INSILICO.out.reads)

    emit:
    background_reads = ch_background_reads   // tuple(meta, reads) — new samples for ALIGNMENT
    manifests        = ch_manifests          // *_subsample_manifest.tsv
    versions         = ch_versions
}
