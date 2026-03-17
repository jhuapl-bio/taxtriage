//
// INSILICO: Optional in-silico read simulation from classification results
//
// Generates simulated FASTQ reads (Illumina via InSilicoSeq, ONT via NanoSim),
// aligns them back to the reference, and produces control JSONs that can be
// compared against real sample results in the report.
//
// Flow:
//   Kraken2 top_report (or custom abundance) + merged_taxid + reference FASTA
//     -> MAKE_SIMULATED_SAMPLES (abundance.tsv + reference.fasta per sample)
//       -> INSILICOSEQ_SIMULATE  (paired R1/R2 fastq.gz)   [if generate_iss]
//       -> PREPARE_NANOSIM_INPUTS + NANOSIM_SIMULATE (ONT fastq.gz) [if generate_nanosim]
//     -> ALIGN_INSILICO_READS (align simulated reads back -> JSON control)
//

include { MAKE_SIMULATED_SAMPLES } from '../../modules/local/make_simulated_samples'
include { INSILICOSEQ_SIMULATE   } from '../../modules/local/insilicoseq'
include { PREPARE_NANOSIM_INPUTS } from '../../modules/local/prepare_nanosim_inputs'
include { NANOSIM_SIMULATE       } from '../../modules/local/nanosim'
include { ALIGN_INSILICO_READS   } from '../../modules/local/align_insilico'


workflow INSILICO {
    take:
    ch_tops              // tuple(meta, top_report.tsv)       from CLASSIFIER/TOP_HITS
    ch_merged_taxid      // tuple(meta, merged_taxid.tsv)     from REFERENCE_PREP
    ch_fastas            // tuple(meta, [fasta1, fasta2, ...]) from REFERENCE_PREP

    main:
    ch_versions = Channel.empty()
    ch_iss_reads = Channel.empty()
    ch_nanosim_reads = Channel.empty()
    ch_insilico_json = Channel.empty()

    ch_empty_file = file("$projectDir/assets/NO_FILE")
    ch_pathogens = Channel.fromPath("$projectDir/assets/pathogen_sheet.csv", checkIfExists: true)
    if (params.pathogens) {
        ch_pathogens = Channel.fromPath(params.pathogens, checkIfExists: true)
    }

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
        ch_iss_reads = INSILICOSEQ_SIMULATE.out.reads
        ch_versions = ch_versions.mix(INSILICOSEQ_SIMULATE.out.versions)

        // ── Step 3a: Align ISS reads back to reference → insilico JSON ──
        // Combine: ISS reads + shared_reference + merged_taxid
        // ISS reads: (meta, [R1.fq.gz, R2.fq.gz])
        // shared_reference: (meta, reference_sequences.fasta) from MAKE_SIMULATED_SAMPLES
        // merged_taxid: (meta, merged_taxid.tsv)
        ch_iss_align_input = INSILICOSEQ_SIMULATE.out.reads
            .join(MAKE_SIMULATED_SAMPLES.out.shared_reference)
            .combine(ch_merged_taxid.map { meta, taxid -> taxid }.first())
        // Now: (meta, [reads], reference.fasta, merged_taxid.tsv)

        // Restructure for ALIGN_INSILICO_READS input
        ch_iss_align = ch_iss_align_input.map { meta, reads, ref, mapping ->
            [meta, reads, ref, mapping]
        }

        ALIGN_INSILICO_READS(
            ch_iss_align,
            ch_pathogens.first(),
            params.assembly ? file(params.assembly) : [],
            params.minmapq ?: 0,
            params.taxdump ? file(params.taxdump) : file("$projectDir/assets/NO_FILE")
        )
        ch_insilico_json = ALIGN_INSILICO_READS.out.json
        ch_versions = ch_versions.mix(ALIGN_INSILICO_READS.out.versions)
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

            ch_nanosim_training = Channel.value(file(params.nanosim_training, checkIfExists: true))

            ch_nanosim_input = PREPARE_NANOSIM_INPUTS.out.nanosim_inputs
                .combine(ch_nanosim_training)

            NANOSIM_SIMULATE(
                ch_nanosim_input
            )
            ch_nanosim_reads = NANOSIM_SIMULATE.out.reads
            ch_versions = ch_versions.mix(NANOSIM_SIMULATE.out.versions)

            // TODO: Add NanoSim alignment step similar to ISS when needed
        }
    }

    emit:
    iss_reads      = ch_iss_reads       // tuple(meta, [R1.fastq.gz, R2.fastq.gz])
    nanosim_reads  = ch_nanosim_reads   // tuple(meta, [ont_reads.fastq.gz])
    insilico_json  = ch_insilico_json   // tuple(meta, insilico.json) for report comparison
    versions       = ch_versions
}
