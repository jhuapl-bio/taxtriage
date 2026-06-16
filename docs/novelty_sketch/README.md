# Novelty detection branch — sketch

Reference-free / open-set detection of organisms with **no near neighbor in the DB**, so
taxtriage can (a) recover a **genus-or-higher** call for divergent organisms and (b) emit a
per-sample **novelty signal** when the DB explains almost nothing in the sample.

The current path is closed-set: Kraken2 top hits → pull reference FASTA → nucleotide
alignment → confidence scoring. Anything with no reference falls through silently. This
branch catches exactly that residual.

## Idea in one line

Take what falls through (reads not aligned to any reference, plus de novo contigs), run a
**translated-search LCA** (protein space stays alignable across genus/family/order), and
fold the result into the existing confidence scoring as a novelty score + candidate calls.

## Files and where they go

| Sketch file | Destination in repo |
|---|---|
| `extract_unmapped.nf` | `modules/local/extract_unmapped.nf` |
| `mmseqs_taxonomy.nf` | `modules/local/mmseqs_taxonomy.nf` |
| `novelty_score.nf` | `modules/local/novelty_score.nf` |
| `novelty_score.py` | `bin/novelty_score.py` (chmod +x) |
| `novelty.nf` | `subworkflows/local/novelty.nf` |

(The `../../modules/...` include paths in `novelty.nf` assume it lives in
`subworkflows/local/`; same for the relative bin call.)

## How to run it

```
nextflow run main.nf <your usual args> \
    --detect_novelty \
    --novelty_db /path/to/mmseqs_seqTaxDB
```

`--detect_novelty` is the master switch; `--novelty_db` points at a prebuilt mmseqs
seqTaxDB (build once, cache — see DB note). Everything else has defaults.

## The score

Transparent, weighted, z-scored — sits next to your custom scoring instead of replacing it:

```
novelty = w_dark * z(dark_fraction)
        + w_rank * z(highrank_only_fraction)
        + w_idnt * z(lowident_tail_mass)
```

- **dark_fraction** — reads explained by *nothing* (not K2-classified, not ref-aligned, not
  protein-assigned). A spike vs. baseline = "the DB doesn't have anything like this."
- **highrank_only_fraction** — reads the translated search places only at genus/family/
  order, never species. This is the direct "we see it at the genus level" signal.
- **lowident_tail_mass** — share of best hits below ~50% aa identity: divergent-but-
  homologous content (a candidate new genus/family rather than pure dark matter).

`z()` is taken against negative/no-template **controls** when available, otherwise the
in-run across-sample distribution (`--baseline auto` with `--run-summaries`). Flag at
`novelty >= 2.0` (≈ 2σ out) by default. Weights/threshold are `task.ext.args` / params, so
you can tune per platform without code changes.

Outputs per sample:
- `*.novelty.summary.tsv` — one row, joins to your mqc/confidence tables by sample id.
- `*.novelty.candidates.tsv` — one row per genus+ candidate taxon with read support.

## Wiring into `workflows/taxtriage.nf`

The subworkflow now does the unmapped-read extraction and the `meta` read-accounting itself
(via `EXTRACT_UNMAPPED`, which surfaces `total_reads`/`ref_aligned`/`k2_classified` as `env`
outputs and folds them onto `meta`). So the caller just hands it the standard channels you
already have in scope after ALIGNMENT/CLASSIFIER/ASSEMBLY:

```groovy
include { NOVELTY } from '../subworkflows/local/novelty'

if (params.detect_novelty) {
    NOVELTY(
        ALIGNMENT.out.bams,                       // [meta, bam, csi]  reference-merged per sample
        ch_kraken2_report,                        // [meta, kreport]   from CLASSIFIER
        ch_denovo,                                // [meta, contigs]   ASSEMBLY.out.ch_denovo_assembly
        file(params.novelty_db, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(NOVELTY.out.versions)
}
```

Place this just after the `ASSEMBLY(...)` block (so `ch_denovo` exists). If `--detect_novelty`
is set without `--use_denovo`, `ch_denovo` may be empty — the subworkflow's `remainder: true`
join handles that and falls back to unmapped-reads-only.

**Join into the report.** Add `NOVELTY.out.summary` to the `input_alignment_files` join with
`remainder: true` + a `NO_FILE` placeholder, exactly like `ch_annotate_report_tsv` does today,
so the novelty columns land in the final per-sample report and MultiQC:

```groovy
ch_novelty_tsv = ALIGNMENT.out.bams.map { meta, bam, csi ->
    [meta, file("$projectDir/assets/NO_FILE")]
}
if (params.detect_novelty) {
    ch_novelty_tsv = ch_novelty_tsv
        .join(NOVELTY.out.summary, remainder: true)
        .map { meta, ph, real -> [meta, real ?: ph] }
}
// ...then add `.join(ch_novelty_tsv)` to the input_alignment_files chain.
```

## Params to add to `nextflow.config` / `nextflow_schema.json`

```groovy
params.detect_novelty   = false   // master switch
params.novelty_db       = null    // path to prebuilt mmseqs seqTaxDB (UniRef50+taxonomy)
params.novelty_flag_z   = 2.0     // flag threshold
params.novelty_weights  = '0.5,0.3,0.2'  // w_dark,w_rank,w_idnt
params.novelty_idnt_cut = 50.0    // aa %% identity tail cutoff
```

## Two-pass baseline (optional, better)

For a true across-sample z-score, run scoring twice: pass 1 emits per-sample summaries,
`collect()` them, then re-run `NOVELTY_SCORE` with `--run-summaries` = the concatenated
file (mark controls with `is_control=1` so they anchor the baseline). The sketch wires the
single-pass version; the `run_summaries` input is already plumbed for the 2-pass swap.

## Build order (matches the earlier recommendation)

1. `MMSEQS_TAXONOMY` on unmapped reads + `dark_fraction`/`highrank_only` scoring — answers
   both "give me a genus call" and "flag the weird sample" with minimal change.
2. Add `lowident_tail_mass` once `convertalis`/diamond pident is flowing.
3. Pool in de novo contigs (already emitted by ASSEMBLY) for more signal on flagged samples.
4. Later: RdRp palmprint (palmscan) / phylogenetic placement (EPA-ng) for *characterizing*
   the flagged novelty — a separate module that consumes `candidates.tsv`.

## DB note

`MMSEQS_TAXONOMY` needs a prebuilt mmseqs **seqTaxDB**. UniRef50 + NCBI taxonomy is a good
default (sensitive, manageable size); build once and cache, point `params.novelty_db` at it.
Kaiju (nr_euk or RVDB) is a drop-in alternative — swap the module, see the commented block
at the bottom of `mmseqs_taxonomy.nf`. For Kaiju, use match-length/read-length as the
divergence proxy instead of pident in `novelty_score.py`.
```
