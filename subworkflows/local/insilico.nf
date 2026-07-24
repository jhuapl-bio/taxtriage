//
// INSILICO: Optional in-silico read simulation from classification results
//
// Generates simulated FASTQ reads (Illumina via InSilicoSeq, ONT via NanoSim)
// and emits them as new samples with insilico meta so they can flow through the
// normal ALIGNMENT → REPORT pipeline.  In REPORT, these insilico samples are
// split out, processed through ALIGNMENT_PER_SAMPLE, and their resulting JSONs
// are used as insilico controls for the real (non-control) samples.
//
// Flow:
//   Kraken2 top_report (or custom abundance) + merged_taxid + reference FASTA
//     -> MAKE_SIMULATED_SAMPLES (abundance.tsv + reference.fasta per sample)
//     -> INSILICOSEQ_SIMULATE  (paired R1/R2 fastq.gz)   [if generate_iss]
//     -> PREPARE_NANOSIM_INPUTS + NANOSIM_SIMULATE (ONT fastq.gz) [if generate_nanosim]
//   Simulated reads are tagged with insilico meta and emitted as new samples.
//

include { MAKE_SIMULATED_SAMPLES } from '../../modules/local/make_simulated_samples'
include { INSILICOSEQ_SIMULATE   } from '../../modules/local/insilicoseq'
include { PREPARE_NANOSIM_INPUTS } from '../../modules/local/prepare_nanosim_inputs'
include { NANOSIM_SIMULATE       } from '../../modules/local/nanosim'
include { SUBSAMPLE_INSILICO     } from '../../modules/local/subsample_insilico'


// ── Resolve the read-count series from either an explicit list or a generator ──
def resolveSeriesCounts() {
    if (params.sim_series_counts) {
        return params.sim_series_counts.toString()
            .split(/[,\s]+/)
            .findAll { it }
            .collect { it as long }
    }
    if (params.sim_series_start != null && params.sim_series_step != null && params.sim_series_n) {
        def start = params.sim_series_start as long
        def step  = params.sim_series_step as long
        def n     = params.sim_series_n as int
        return (0..<n).collect { start + it * step }
    }
    error "sim_subsample is enabled but no series was defined. Set --sim_series_counts " +
          "(e.g. '100,500,1000') OR the generator trio --sim_series_start/--sim_series_step/--sim_series_n."
}


// ── Regroup the flat list of per-dataset FASTQs emitted by SUBSAMPLE_INSILICO ──
// into one (meta, reads) tuple per derived dataset, re-tagged as its own sample.
def regroupSubsampled(ch_subsample_reads) {
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
            // parent_id is preserved from the master insilico meta so reference-prep
            // cloning downstream still points at the real originating sample.
            def cm = (dsid =~ /_c(\d+)_r\d+$/)
            if (cm.find()) {
                m.read_count = cm.group(1) as Integer
            }
            m.subsample       = true
            m.subsample_mode  = params.sim_subsample_mode
            def sorted = fastqs.sort { it.getName() }
            [m, sorted.size() == 1 ? sorted[0] : sorted]
        }
}


workflow INSILICO {
    take:
    ch_tops              // tuple(meta, top_report.tsv)       from CLASSIFIER/TOP_HITS
    ch_merged_taxid      // tuple(meta, merged_taxid.tsv)     from REFERENCE_PREP
    ch_fastas            // tuple(meta, [fasta1, fasta2, ...]) from REFERENCE_PREP

    main:
    ch_versions = Channel.empty()
    ch_insilico_reads = Channel.empty()
    ch_subsample_manifests = Channel.empty()

    ch_empty_file = file("$projectDir/assets/NO_FILE")

    // ── Resolve abundance source ────────────────────────────────────────
    if (params.sim_abundance) {
        ch_abundance_input = file(params.sim_abundance, checkIfExists: true)
    } else {
        ch_abundance_input = ch_empty_file
    }

    // ── Join inputs by sample meta ──────────────────────────────────────
    ch_sim_input = ch_tops
        .join(ch_merged_taxid)
        .join(ch_fastas)
    // Now: tuple(meta, top_report, merged_taxid, [fasta1, fasta2, ...])

    // ── Step 1: Generate abundance files + per-sample reference FASTAs ──
    MAKE_SIMULATED_SAMPLES(
        ch_sim_input,
        params.sim_nsamples    ?: 1,
        params.sim_ranks       ?: 'S S1 S2 S3',
        params.sim_minreads    ?: 3,
        params.sim_exclude_taxids ?: '9606',
        params.sim_include_taxids ?: '',
        ch_abundance_input
    )
    ch_versions = ch_versions.mix(MAKE_SIMULATED_SAMPLES.out.versions)

    // Flatten per-sample outputs so each sample runs independently
    ch_samples = MAKE_SIMULATED_SAMPLES.out.samples.transpose()
    // Now: tuple(meta, abundance.tsv, reference.fasta) per sample

    // ── Step 2a: InSilicoSeq (Illumina paired-end) ──────────────────────
    if (params.generate_iss) {
        ch_iss_input = ch_samples.map { meta, abu, ref ->
            [meta, ref, abu]
        }

        def iss_nreads = params.sim_nreads ?: 100000

        INSILICOSEQ_SIMULATE(
            ch_iss_input,
            iss_nreads,
            params.iss_model ?: 'miseq'
        )
        ch_versions = ch_versions.mix(INSILICOSEQ_SIMULATE.out.versions)

        // Tag ISS reads with insilico meta so they can flow as new samples
        ch_iss_tagged = INSILICOSEQ_SIMULATE.out.reads.map { meta, reads ->
            def insilico_meta = meta.collectEntries { k, v -> [k, v] }
            insilico_meta.parent_id = meta.id
            insilico_meta.id = "${meta.id}_insilico_iss"
            insilico_meta.insilico = true
            insilico_meta.control = false
            insilico_meta.platform = 'ILLUMINA'
            insilico_meta.single_end = false
            insilico_meta.trim = false
            insilico_meta.read_count = iss_nreads
            [insilico_meta, reads]
        }
        ch_insilico_reads = ch_insilico_reads.mix(ch_iss_tagged)
    }

    // ── Step 2b: NanoSim (ONT long reads) ───────────────────────────────
    if (params.generate_nanosim) {
        if (!params.nanosim_training) {
            log.warn "WARNING: --generate_nanosim is set but --nanosim_training is not provided. NanoSim simulation will be skipped."
        } else {
            def iss_nreads = params.sim_nreads ?: 100000
            def ont_divisor = params.sim_ont_divisor ?: 40
            def ont_nreads  = Math.max(1, (int)(iss_nreads / ont_divisor))

            ch_nanosim_prep_input = ch_samples.combine(ch_merged_taxid, by: 0)

            PREPARE_NANOSIM_INPUTS(
                ch_nanosim_prep_input,
                ont_nreads
            )
            ch_versions = ch_versions.mix(PREPARE_NANOSIM_INPUTS.out.versions)
            ch_genomes = PREPARE_NANOSIM_INPUTS.out.genomes
            ch_nanosim_training = Channel.value(file(params.nanosim_training, checkIfExists: true))
            ch_nanosim_input = PREPARE_NANOSIM_INPUTS.out.nanosim_inputs
                .combine(ch_nanosim_training).combine(ch_genomes, by: 0)

            NANOSIM_SIMULATE(
                ch_nanosim_input
            )
            ch_versions = ch_versions.mix(NANOSIM_SIMULATE.out.versions)

            // Tag NanoSim reads with insilico meta
            ch_nanosim_tagged = NANOSIM_SIMULATE.out.reads.map { meta, reads ->
                def insilico_meta = meta.collectEntries { k, v -> [k, v] }
                insilico_meta.parent_id = meta.id
                insilico_meta.id = "${meta.id}_insilico_nanosim"
                insilico_meta.insilico = true
                insilico_meta.control = false
                insilico_meta.platform = 'OXFORD'
                insilico_meta.single_end = true
                insilico_meta.trim = false
                insilico_meta.read_count = ont_nreads
                [insilico_meta, reads]
            }
            ch_insilico_reads = ch_insilico_reads.mix(ch_nanosim_tagged)
        }
    }

    // ── Step 3: Subsample master datasets into spike-in / dilution series ────
    // Up to here, ch_insilico_reads holds the *master* synthetic datasets (one
    // per platform per sample). When --sim_subsample is set we replace each
    // master with a series of derived datasets (nested "consistent" spike-ins,
    // independent "randomized" spike-ins, or a replicated dilution series). Each
    // derived dataset is recorded as a list of read indices into its master so
    // the FASTQs can be deleted and reconstructed later.
    if (params.sim_subsample) {
        def counts = resolveSeriesCounts().join(',')
        def reps   = params.sim_series_replicates ?: 1
        def mode   = params.sim_subsample_mode ?: 'randomized'
        def seed   = params.sim_subsample_seed ?: 42

        log.info "INSILICO subsampling: mode=${mode} counts=[${counts}] replicates=${reps} seed=${seed}"

        SUBSAMPLE_INSILICO(
            ch_insilico_reads,   // master datasets (tagged with insilico meta)
            mode,
            counts,
            reps,
            seed
        )
        ch_versions = ch_versions.mix(SUBSAMPLE_INSILICO.out.versions)

        // Only the derived subsamples flow downstream; the master is the source pool.
        ch_insilico_reads = regroupSubsampled(SUBSAMPLE_INSILICO.out.reads)

        // Manifests (per parent×platform) carry authoritative counts for the report suite tab.
        ch_subsample_manifests = SUBSAMPLE_INSILICO.out.manifest
    }

    emit:
    insilico_reads = ch_insilico_reads          // tuple(insilico_meta, reads) — new samples for ALIGNMENT
    subsample_manifests = ch_subsample_manifests // *_subsample_manifest.tsv — authoritative counts for report
    versions       = ch_versions
}
