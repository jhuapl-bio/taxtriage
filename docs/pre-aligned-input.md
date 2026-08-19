# Pre-Aligned Input (BAM / CRAM / SAM)

TaxTriage can start from reads you have already aligned yourself, instead of from raw FASTQ. Point it at a BAM, CRAM or SAM file and the pipeline skips everything upstream of alignment and goes straight to coverage statistics and [TASS scoring](tass-scoring.md).

This is the right entry point when you already have an alignment you trust, when you want to re-score an old run under different TASS weights without re-aligning, or when your alignment came from a tool or reference set TaxTriage would not have chosen on its own.

> The one thing to get right is the reference. TaxTriage needs to know what your BAM was aligned *against*, and a BAM header does not carry that information in full. See [What TaxTriage needs the reference for](#what-taxtriage-needs-the-reference-for) below.

---

## Two ways to supply an alignment

### Single sample: `--bam`

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  --bam my_sample.bam \
  --reference_fasta refs.fasta \
  --type stool \
  --sample my_sample \
  --outdir ./results \
  -profile local,docker \
  -resume
```

### Samplesheet: the `bam` column

Add a `bam` column and put the alignment path there instead of a path in `fastq_1`:

```
sample,platform,bam,type
prealigned_a,ILLUMINA,/data/prealigned_a.bam,stool
prealigned_b,OXFORD,/data/prealigned_b.cram,nasal
```

A samplesheet must contain **either** a `fastq_1` column or a `bam` column. If both are present, a row that fills in `bam` is treated as pre-aligned and its `fastq_1` value is ignored, so you can mix the two input styles in one sheet:

```
sample,platform,fastq_1,bam,type
raw_sample,ILLUMINA,/data/raw_R1.fastq.gz,,stool
prealigned,ILLUMINA,,/data/prealigned.bam,stool
```

In a mixed sheet the FASTQ rows run the full pipeline as normal and only the pre-aligned rows take the short path. All the usual extra columns (`type`, `positive`/`negative`, metadata columns) work the same on both.

### Accepted formats

| Extension | Handling |
|---|---|
| `.bam` | Used as-is |
| `.sam` | Compressed to BAM first |
| `.cram` | Decoded to BAM first. Needs a reference, see [CRAM](#cram) |

Your file does **not** need to be sorted or indexed beforehand. `PREPARE_BAM` reads the `@HD SO:` header tag and coordinate-sorts the file only if it is not already in coordinate order, then writes a `.csi` index. Both are required by `samtools coverage` and by the read-count indexing in `match_paths.py`.

---

## What gets skipped

Pre-aligned samples bypass every step that operates on raw reads or on reference selection, because the alignment has already settled both questions:

- read QC and quality trimming (`trim` is forced to `false` on these rows)
- host read removal
- Kraken2 / Bracken classification
- reference download from NCBI

When **every** sample in the run is pre-aligned, TaxTriage reports a `BAM-only run detected` message and additionally forces `--skip_kraken2` and `--skip_refpull` on. Any of the following that you enabled are switched off with a warning, because each one needs raw reads or de novo contigs that a BAM cannot provide:

`--use_denovo`, `--use_diamond`, `--annotate`, `--microbert`, `--novelty`, `--generate_iss`, `--generate_nanosim`, `--reference_assembly`, `--get_variants`

Two consequences worth knowing:

- **No Kraken2 database is required** for a BAM-only run. The usual "if `--skip_kraken2` is false you must provide `--db`" check does not apply, and neither does the rule that skipping Kraken2 obliges you to supply a reference or organism list. Pre-aligned samples are exempt because their references are already fixed by the BAM.
- **No VF/AMR annotation, novelty detection or de novo assembly** is available for these samples. If you need those, start from FASTQ.

In a mixed FASTQ/BAM run none of the global flags are forced off. Only the pre-aligned samples bypass those stages; your FASTQ samples still get classification, annotation and everything else.

---

## What TaxTriage needs the reference for

A BAM header lists reference **names and lengths**, but not their **sequence**. `match_paths.py` can read the lengths straight off the header, which is enough for coverage, breadth and read counts. Two other things need actual bases:

1. **The accession → taxid map** (`-m`), without which detections cannot be tied to organisms.
2. **Sourmash sketching, the shared-window conflict report and the ANI matrix** (`-f`), which drive the minhash component and conflict-based read removal (see [TASS Scoring §3](tass-scoring.md#3-stage-2-conflict-region-detection-conflict_regionspy)).

So supply `--reference_fasta` pointing at the reference your BAM was aligned against whenever you have it. It is fuzzy-matched against the NCBI assembly summary to build the taxid map, exactly as a local reference is on the FASTQ path.

### When you don't supply a reference

If a pre-aligned sample has neither `--reference_fasta` nor `--get_pathogens`, TaxTriage reconstructs the reference sequence by calling a consensus off the alignment itself (`BAM_CONSENSUS`, using `samtools consensus`), then treats that FASTA as though you had supplied it. You will see:

```
NOTE: pre-aligned input without --reference_fasta -> reconstructing reference
sequence from the alignment (samtools consensus).
```

The consensus is run with `-a`, so every reference stays at full length with uncovered positions written as `N`. That matters because the shared-window sketching slides fixed windows along these sequences and would otherwise be comparing mismatched coordinates. References whose consensus comes back almost entirely `N`, meaning nothing in the BAM header that no read actually supports, are dropped rather than passed on as empty sketches that would drag every ANI comparison toward zero.

> ⚠️ **A consensus is not as good as the real reference, and it biases TASS in a specific direction.** It only covers positions that have aligned reads, so absent regions look like absent sequence rather than unsequenced sequence. More importantly, a multi-mapping read contributes to *every* reference it was placed on, which makes related organisms look more similar to each other than they are. Since the ANI and shared-window numbers drive conflict-based read removal, that inflated similarity makes removal **more aggressive** than it should be. Pass `--reference_fasta` whenever the true reference is available.

### Controlling it

| Value | Behaviour |
|---|---|
| `--bam_consensus` unset (default) | Auto: consensus is derived only when a BAM sample has no `--reference_fasta` and no `--get_pathogens` |
| `--bam_consensus true` | Always derive a consensus, even alongside a supplied reference |
| `--bam_consensus false` | Never derive one |

Setting `--bam_consensus false` with no reference available is a valid choice, but it leaves the run with no bases at all, so sourmash/ANI comparison and conflict-based read removal are both disabled. TaxTriage warns you and suggests rebalancing:

```
WARNING: pre-aligned input without --reference_fasta and --bam_consensus false:
no reference sequence is available, so sourmash/ANI comparison and conflict-based
read removal are disabled. Set --minhash_weight 0 to rebalance the TASS weights.
```

Take that advice seriously. Minhash carries one of the three largest default weights (0.31, see [TASS Scoring §9.2](tass-scoring.md#92-default-weights)), so leaving it weighted while its input is unavailable depresses every score in the run.

### Consensus tuning

| Parameter | Default | Description |
|---|---|---|
| `--consensus_min_depth <int>` | `1` | `samtools consensus -d`: minimum read depth to call a base |
| `--consensus_min_bases <int>` | `500` | Drop consensus records with fewer than this many non-`N` bases |
| `--consensus_min_mapq <int>` | unset | `samtools consensus --min-MQ` |
| `--consensus_mode <simple\|bayesian>` | `simple` | `samtools consensus -m` |

---

## CRAM

CRAM is decoded to BAM before anything else runs, and decoding needs the reference the file was compressed against. TaxTriage uses `--reference_fasta` for this by default. If the CRAM was compressed against a different reference, point `--cram_reference` at that one:

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  --bam my_sample.cram \
  --cram_reference /refs/the_one_the_cram_was_made_with.fasta \
  --reference_fasta refs.fasta \
  --type stool \
  --outdir ./results \
  -profile local,docker
```

---

## MAPQ filtering

A pre-aligned file is taken **as given**. The `--minmapq` default of 5 applies to alignments TaxTriage produces itself and is deliberately not applied to yours, on the assumption that you already filtered to taste.

If you do want a filter applied, set `--bam_minmapq <int>` and it is passed to `samtools view -q` during preparation.

---

## Parameter summary

| Parameter | Description |
|---|---|
| `--bam <path>` | Pre-aligned `.bam`, `.cram` or `.sam`. Single-sample shortcut; the samplesheet equivalent is a `bam` column |
| `--reference_fasta <path>` | The reference the alignment was made against. Strongly recommended |
| `--bam_minmapq <int>` | Optional MAPQ filter on pre-aligned input. Unset means take the file as given |
| `--cram_reference <path>` | Reference used to decode CRAM, when different from `--reference_fasta` |
| `--bam_consensus <bool>` | Recover reference sequence from the alignment. Auto when no reference is supplied |
| `--consensus_min_depth <int>` | Minimum depth to call a consensus base (default `1`) |
| `--consensus_min_bases <int>` | Minimum non-`N` bases to keep a consensus record (default `500`) |
| `--consensus_min_mapq <int>` | Minimum MAPQ for consensus calling |
| `--consensus_mode <str>` | `simple` (default) or `bayesian` |

---

## Next Steps

- [Samplesheet](samplesheet.md): the `bam` column alongside every other input option
- [CLI Parameters](cli-parameters.md#pre-aligned-input-bam--cram--sam): full parameter reference
- [TASS Scoring](tass-scoring.md): what happens once coverage statistics exist
- [Troubleshooting](troubleshooting.md): common failures
