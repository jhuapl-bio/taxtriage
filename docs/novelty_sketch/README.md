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

| Sketch file            | Destination in repo                  |
| ---------------------- | ------------------------------------ |
| `mmseqs_downloaddb.nf` | `modules/local/mmseqs_downloaddb.nf` |
| `extract_unmapped.nf`  | `modules/local/extract_unmapped.nf`  |
| `mmseqs_taxonomy.nf`   | `modules/local/mmseqs_taxonomy.nf`   |
| `novelty_score.nf`     | `modules/local/novelty_score.nf`     |
| `novelty_score.py`     | `bin/novelty_score.py` (chmod +x)    |
| `novelty.nf`           | `subworkflows/local/novelty.nf`      |

(The `../../modules/...` include paths in `novelty.nf` assume it lives in
`subworkflows/local/`; same for the relative bin call.)

## How to run it

```
nextflow run main.nf <your usual args> \
    --detect_novelty \
    --novelty_db UniProtKB        # name -> auto-download, OR a local path -> used as-is
```

`--detect_novelty` is the master switch. `--novelty_db` is **local-first, like `--db`**:

- if it's an existing path/dir → used directly as a prebuilt seqTaxDB (no download);
- otherwise it's treated as an `mmseqs databases` name and downloaded once, then cached.

Default is `UniProtKB`. Everything else has defaults.

## The DB step (mirrors Kraken2/Centrifuge `--db`)

`MMSEQS_DOWNLOADDB` wraps `mmseqs databases <Name> <out>/seqTaxDB tmp`, exactly the form from
the MMseqs2 wiki (`mmseqs databases UniProtKB/Swiss-Prot outpath/swissprot tmp`). It always
writes the db under the fixed prefix `seqTaxDB` inside a `mmseqs_db/` directory, so everything
downstream references `<dir>/seqTaxDB` no matter which source db was chosen.

**Caching.** The module uses `storeDir "${params.novelty_db_cache}/<db_name>"`. `storeDir`
persists the finished db to that folder and **skips the process whenever the files already
exist** — so the multi-GB download happens once and is reused by `-resume` (within a run) _and_
by every later run (across runs). The per-db-name subfolder means switching `--novelty_db`
won't clobber a previously downloaded db.

Supported names (taxonomy-enabled, usable as seqTaxDB): `UniProtKB` (default),
`UniProtKB/Swiss-Prot`, `UniProtKB/TrEMBL`, `UniRef100`, `UniRef90`, `UniRef50`, `NR`, `GTDB`,
`SILVA`, `Kalamari`. (Profile/nucleotide-only entries like Pfam/PDB/NT are not seqTaxDBs and
won't work for `mmseqs taxonomy`.)

Resolution snippet for `workflows/taxtriage.nf` (drop in near the Kraken2 `--db` block):

```groovy
include { MMSEQS_DOWNLOADDB } from '../modules/local/mmseqs_downloaddb'

ch_novelty_db = Channel.empty()
if (params.detect_novelty) {
    if (params.novelty_db && file(params.novelty_db).exists()) {
        println "Using local mmseqs seqTaxDB at ${params.novelty_db}"
        ch_novelty_db = Channel.value(file(params.novelty_db, checkIfExists: true))
    } else {
        def dbname = params.novelty_db ?: 'UniProtKB'
        println "mmseqs seqTaxDB '${dbname}' not found locally; downloading via 'mmseqs databases' " +
                "(cached at ${params.novelty_db_cache})."
        MMSEQS_DOWNLOADDB(dbname)
        ch_novelty_db = MMSEQS_DOWNLOADDB.out.db.first()   // value channel, reused per sample
    }
}
```

> Local-db note: a local `--novelty_db` should be a **directory** containing files prefixed
> `seqTaxDB` (or set `--novelty_db_prefix` to your prefix). This matches how `--db` points at a
> directory for Kraken2.

## The score

Transparent, weighted, z-scored — sits next to your custom scoring instead of replacing it:

```
novelty = w_dark * z(dark_fraction)
        + w_rank * z(highrank_only_fraction)
        + w_idnt * z(lowident_tail_mass)
```

- **dark_fraction** — reads explained by _nothing_ (not K2-classified, not ref-aligned, not
  protein-assigned). A spike vs. baseline = "the DB doesn't have anything like this."
- **highrank_only_fraction** — reads the translated search places only at genus/family/
  order, never species. The direct "we see it at the genus level" signal.
- **lowident_tail_mass** — share of best hits below ~50% aa identity: divergent-but-
  homologous content (a candidate new genus/family rather than pure dark matter).

`z()` is taken against negative/no-template **controls** when available, otherwise the
in-run across-sample distribution (`--baseline auto` with `--run-summaries`). Flag at
`novelty >= 2.0` (≈ 2σ out) by default.

Outputs per sample:

- `*.novelty.summary.tsv` — one row, joins to your mqc/confidence tables by sample id.
- `*.novelty.candidates.tsv` — one row per genus+ candidate taxon with read support.

## Wiring into `workflows/taxtriage.nf`

The subworkflow does the unmapped-read extraction and the `meta` read-accounting itself (via
`EXTRACT_UNMAPPED`, which surfaces `total_reads`/`ref_aligned`/`k2_classified` as `env` outputs
and folds them onto `meta`). So the caller just hands it the standard channels already in scope
after ALIGNMENT/CLASSIFIER/ASSEMBLY:

```groovy
include { NOVELTY } from '../subworkflows/local/novelty'

if (params.detect_novelty) {
    NOVELTY(
        ALIGNMENT.out.bams,        // [meta, bam, csi]  reference-merged per sample
        ch_kraken2_report,         // [meta, kreport]   from CLASSIFIER
        ch_denovo,                 // [meta, contigs]   ASSEMBLY.out.ch_denovo_assembly
        ch_novelty_db              // resolved above
    )
    ch_versions = ch_versions.mix(NOVELTY.out.versions)
}
```

Place this just after the `ASSEMBLY(...)` block (so `ch_denovo` exists). If `--detect_novelty`
is set without `--use_denovo`, `ch_denovo` may be empty — the `remainder: true` join handles
that and falls back to unmapped-reads-only.

**Join into the report.** Add `NOVELTY.out.summary` to the `input_alignment_files` join with
`remainder: true` + a `NO_FILE` placeholder, exactly like `ch_annotate_report_tsv` does today:

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

## Params to add (`nextflow.config` / `nextflow_schema.json`)

```groovy
params.detect_novelty   = false              // master switch
params.novelty_db       = 'UniProtKB'        // mmseqs db NAME (auto-download) or local path
params.novelty_db_cache = "${projectDir}/dbs/mmseqs"  // storeDir cache (persists across runs)
params.novelty_db_prefix= 'seqTaxDB'         // base name inside the db dir
params.novelty_flag_z   = 2.0                // flag threshold
params.novelty_weights  = '0.5,0.3,0.2'      // w_dark,w_rank,w_idnt
params.novelty_idnt_cut = 50.0               // aa %% identity tail cutoff
```

## Two-pass baseline (optional, better)

For a true across-sample z-score, run scoring twice: pass 1 emits per-sample summaries,
`collect()` them, then re-run `NOVELTY_SCORE` with `--run-summaries` = the concatenated file
(mark controls with `is_control=1` so they anchor the baseline). The sketch wires the
single-pass version; the `run_summaries` input is already plumbed for the 2-pass swap.

## Build order

1. `MMSEQS_TAXONOMY` on unmapped reads + `dark_fraction`/`highrank_only` scoring — answers
   both "give me a genus call" and "flag the weird sample" with minimal change.
2. Add `lowident_tail_mass` once `convertalis` pident is flowing.
3. Pool in de novo contigs (already emitted by ASSEMBLY) for more signal on flagged samples.
4. Later: RdRp palmprint (palmscan) / phylogenetic placement (EPA-ng) for _characterizing_
   the flagged novelty — a separate module that consumes `candidates.tsv`.

## DB note

`mmseqs taxonomy` needs a taxonomy-enabled seqTaxDB. `UniProtKB` (default) is broad;
`UniRef50` is smaller/faster; `GTDB` is bacteria/archaea-focused. Kaiju (nr_euk / RVDB) is a
drop-in alternative — swap the module, see the commented block at the bottom of
`mmseqs_taxonomy.nf`. For Kaiju, use match-length/read-length as the divergence proxy instead
of pident in `novelty_score.py`.
