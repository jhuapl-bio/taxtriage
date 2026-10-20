# Region Masking: FCS Contamination & Targeted TASS

TaxTriage is built for untargeted, agnostic detection, so by default every organism is scored over the **whole length** of its reference genome: breadth, depth, Gini evenness and every other TASS input use all of it.

Sometimes only part of a genome should count:

- **Contaminated references.** NCBI's [Foreign Contamination Screen (FCS)](https://github.com/ncbi/fcs) finds stretches of assemblies that come from another organism - a human fragment inside a bacterial assembly, a cloning vector, a misassembled contig from a different genus. Reads that land there are evidence for the _contaminant_, not the organism the assembly is named after.
- **Targeted questions.** You may only care about a toxin locus, a set of marker genes, an amplicon panel, or a lineage-defining region, and want TASS computed from those areas **only**.

Region masking handles both. It is off unless you ask for it, and a run without it is byte-for-byte unchanged.

---

## Two modes

| Mode                  | You give                                                                                                | What is scored                                            |
| --------------------- | ------------------------------------------------------------------------------------------------------- | --------------------------------------------------------- |
| `exclude` _(default)_ | Regions to **ignore**: NCBI FCS spans (`--fcs`), adaptor/vector reports, your own BED (`--mask_bed`)    | Everything except those regions                           |
| `target`              | Regions to **score**: your BED (`--mask_bed`) or genes picked from annotation (`--mask_feature_select`) | Only the targets, minus any FCS contaminant spans in them |

In target mode, references with no target at all are dropped from scoring (`--mask_untargeted drop`, the default) so TASS comes **only** from the targeted areas. Use `--mask_untargeted keep` to score untargeted references over their whole length as usual.

---

## How it works

Masking is applied once, while `match_paths.py` reads the BAM, and every downstream metric inherits it.

1. **Each reference gets a set of _kept_ intervals:** the complement of the exclude set, or the targets themselves (minus exclusions).
2. **Every read is clipped to the kept intervals.** A read with at least `--mask_read_overlap` (default 0.5) of its aligned span in masked space is dropped entirely, before any metric (including the high-MAPQ fraction) sees it. A kept read contributes only its unmasked bases.
3. **Coordinate compression.** The kept intervals are laid end to end into a new coordinate space whose length is the kept bp, and the surviving read pieces are translated into it. The masked bases simply do not exist in that genome.

```
reference      |=====[########]=========[####]=========|   # = masked
kept space     |=====|=========|=========|              length = kept bp
```

Because the genome the scorer sees has no masked bases in it, nothing else needs to know about masking:

- **Breadth** = covered kept bases ÷ kept bp. A masked region can neither raise breadth (contaminant pile-ups) nor lower it (a region nobody expected reads in).
- **Mean depth, RPKM, abundance confidence** use the kept length.
- **Gini evenness** is computed over the compressed genome, so a pile-up in a masked contaminant span cannot make coverage look uneven, and a deliberately uncovered region cannot look like a gap.
- **Species / genus aggregation** sums the kept lengths of member references, exactly as it sums full lengths in an unmasked run.

The original BAM is never modified: the mask is applied to the reads as they are scored. If you also want the masked alignment as a file - to open in IGV, or to feed another tool - add `--mask_write_bam`. It writes `masking/bam/<sample>.masked.bam` (+ `.bai`) using exactly the same rule as scoring: every alignment with at least `--mask_read_overlap` of its span masked is removed (primary, secondary and supplementary), unmapped reads pass through, and a kept alignment that partly overlaps the mask carries a `ZM:i` tag with the number of masked bases. The BAM header records the mask in an `@PG` line (`ID:taxtriage-mask`).

### Optional penalty

Ignoring masked reads is usually what you want. If reads piling into masked space should also count **against** an organism - typically, an assembly whose "signal" is mostly its contaminant contig - set `--mask_penalty_power`:

$$
\text{TASS}_{\text{masked}} = \text{TASS} \times (1 - f_{\text{masked}})^{p}
$$

where $f_{\text{masked}}$ is the fraction of the organism's primary reads that fell in masked space. `p = 0` (default) turns this off; `p = 1` scales TASS linearly. In target mode the same factor rewards organisms whose reads concentrate on the targets. The factor is applied at strain, species and genus level and recorded as `mask_penalty` on every organism.

### What is _not_ masked

- **Conflict / shared-region read removal** and the **MinHash** comparison run on the full alignment before scoring. They decide which reads belong to which organism; masking then decides which of those reads count.
- For that reason the species/genus **union-breadth and LCA read-count overrides** (which come from that unmasked step) are switched off in masked runs, and species/genus breadth is pooled from the masked strain values instead. The log says so: `[mask] species/genus union-breadth & LCA read-count overrides disabled`.
- The **Alignment Trends** depth profile is computed in the compressed (scored) coordinates. The report's **positional breadth histogram** is drawn in real genome coordinates with the masked spans overlaid (see [Visualization](#visualization)).

---

## NCBI FCS

FCS-GX screens assemblies against a large cross-kingdom database and assigns every contaminant span an **action**:

| Action    | Meaning                                | Masked by default        |
| --------- | -------------------------------------- | ------------------------ |
| `EXCLUDE` | The whole sequence is contaminant      | yes                      |
| `TRIM`    | Contaminant at a sequence end          | yes                      |
| `FIX`     | Contaminant inside a sequence          | yes                      |
| `REVIEW`  | Possible contaminant, lower confidence | add with `--fcs_actions` |
| `INFO`    | Informational (e.g. organelle)         | add with `--fcs_actions` |

NCBI publishes the results for every GenBank and RefSeq assembly it has screened, under [`genomes/ASSEMBLY_REPORTS`](https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/):

| File                         | Content                                                                               |
| ---------------------------- | ------------------------------------------------------------------------------------- |
| `fcs_details_genbank.txt.gz` | One row per contaminant span in GenBank assemblies (≈53 M rows, 4.6 GB uncompressed)  |
| `fcs_details_refseq.txt.gz`  | The same for RefSeq assemblies                                                        |
| `fcs_summary_genbank.txt.gz` | One row per GenBank assembly: contaminated length, fraction, action counts, divisions |
| `fcs_summary_refseq.txt.gz`  | The same for RefSeq assemblies                                                        |

Details columns: `seq_id start_pos end_pos seq_len action contam_type coverage contam_details assembly_accession` (1-based, inclusive - converted to BED automatically).

### Using it

```bash
# Let TaxTriage download NCBI's reports once (cached in --fcs_db)
nextflow run jhuapl-bio/taxtriage ... --fcs --fcs_db /shared/taxtriage/fcs_db

# Or point at local copies (plain or .gz, comma-separated)
nextflow run jhuapl-bio/taxtriage ... \
    --fcs_details /data/fcs_details_genbank.txt,/data/fcs_details_refseq.txt.gz \
    --fcs_summary /data/fcs_summary_genbank.txt,/data/fcs_summary_refseq.txt
```

The mask is built **once per run** (`BUILD_MASK`), restricted to the references that actually appear in the samples' BAM headers, so the multi-GB details file is streamed a single time - about 90 s for the full GenBank file on a laptop. The matching rows are published as `masking/fcs_subset.txt`; pass that back as `--fcs_details` on later runs against the same references to skip the scan entirely.

### Accession matching

FCS rows use GenBank sequence accessions (`CP012345.1`, `JAAXYZ010000083.1`); TaxTriage references are usually RefSeq (`NZ_CP012345.1`). Matching is prefix- and version-tolerant in both directions, so a GenBank FCS row masks the corresponding `NZ_` RefSeq copy. RefSeq-only accessions (`NC_…`) are matched by `fcs_details_refseq`.

### RefSeq references and the assembly flag

RefSeq assemblies are already screened, and many contaminated sequences are removed or trimmed before they reach RefSeq, so a run built on RefSeq references often has **few or no** per-sequence FCS spans to mask. That is expected, not a failure.

The per-assembly summary is still useful. With `--fcs_summary` (or `--fcs`), every organism whose reference assembly NCBI found contaminated gets an `fcs_assembly` annotation - assembly, contaminated bp, contaminated fraction, number of contaminant sequences - and `masking/fcs_assemblies.tsv` lists every assembly in the run. This changes no score; it tells you which detections rest on an assembly NCBI has flagged.

### Narrowing what is masked

```bash
--fcs_actions EXCLUDE,TRIM,FIX,REVIEW   # stricter: include low-confidence calls
--fcs_actions all                        # every FCS row
--fcs_divisions 'prok|synt'              # only prokaryotic and synthetic contaminants
--fcs_min_len 200                        # ignore spans shorter than 200 bp
--mask_pad 50                            # widen every span by 50 bp each side
```

### Running FCS yourself

For references NCBI has not screened (custom assemblies, MAGs, new genomes), run FCS directly and hand TaxTriage the reports. See the [FCS wiki](https://github.com/ncbi/fcs/wiki) for installation; the full GX database is several hundred GiB and is best run on a large-memory host.

```bash
# FCS-GX: contamination from other organisms
python3 fcs.py screen genome --fasta my_assembly.fna --out-dir gx_out/ \
    --gx-db "$GXDB_LOC/gxdb" --tax-id 1280
# FCS-adaptor: adaptor and vector
./run_fcsadaptor.sh --fasta-input my_assembly.fna --output-dir adaptor_out/ --prok \
    --container-engine singularity --image fcs-adaptor.sif

nextflow run jhuapl-bio/taxtriage ... \
    --fcs_gx_report gx_out/my_assembly.fna.1280.fcs_gx_report.txt \
    --fcs_adaptor_report adaptor_out/fcs_adaptor_report.txt
```

Per-assembly `*_fcs_report.txt` files from an assembly's NCBI FTP folder use the same format and are accepted by `--fcs_gx_report` too.

---

## Custom masks: BED or masked FASTA

Any mask you make yourself works, in either of two forms:

- **BED** (`--mask_bed`): `accession  start  end  [label]`, 0-based, half-open, `.gz` accepted. In exclude mode these regions are ignored; in target mode they are the only regions scored.
- **Masked FASTA** (`--mask_fasta`): a copy of the reference FASTA with the unwanted bases masked - what FCS-GX `clean`, `dustmasker`, RepeatMasker or `bedtools maskfasta` produce. Masked bases are never scored, in either mode. `--mask_fasta_case hard` (default) reads runs of `N`, `soft` reads lowercase bases, `both` reads both.

```bash
nextflow run jhuapl-bio/taxtriage ... --mask_fasta refs.masked.fasta.gz
nextflow run jhuapl-bio/taxtriage ... --mask_fasta refs.softmasked.fa --mask_fasta_case soft
```

Give TaxTriage the masked FASTA **as a mask, not as the alignment reference**. Aligning to an N-masked reference would also stop reads from landing in those regions, and they would then pile onto other organisms instead of being visibly set aside.

Reference FASTAs already contain native `N`s (scaffold gaps, single ambiguity bases). Those are masked too. That is harmless because no reads align to them, but `--mask_fasta_min_run 10` skips isolated `N`s and keeps the mask summary tidy.

### Making a test mask

`bin/make_dummy_masked_fasta.py` writes a masked copy of a reference FASTA, which is a quick way to see what masking does to a run before using a real mask. Give it either a region table or ask for random regions:

```bash
# Fixed regions: seq_id, start, end, kind (hard = N, soft = lowercase), note; 0-based, half-open, end may be "end"
bin/make_dummy_masked_fasta.py -i refs.fasta -r regions.tsv -o refs.masked.fasta.gz

# 1-3 random regions per accession listed in a TaxTriage mapping file (reproducible with --seed)
bin/make_dummy_masked_fasta.py -i refs.fasta --mapping sample.dwnld.gcfmapping.tsv \
    --random 3 --seed 42 -o refs.masked.fasta.gz

nextflow run jhuapl-bio/taxtriage ... --mask_fasta refs.masked.fasta.gz --mask_compare --mask_write_bam
```

The regions actually applied are written next to the output (`refs.masked.fasta.regions.tsv`). `--max_frac` (default 0.05) caps how much of each sequence random regions cover and `--soft_prob` (default 0.2) sets the share of lowercase regions.

What to expect: organisms lose roughly the reads that fell in their masked spans, breadth is recomputed over the remaining bases (it can rise slightly when the masked span held little coverage), and an organism whose reads were concentrated in a masked span can drop below the detection minimum and disappear. In `all.odr.html` the header shows **⚠ Region masking: …**, the Table tab gains the Masked Reads / Masked bp / Mask Penalty columns, and the breadth tracks shade the masked spans; `masking/all.masking.html` lists any detection that disappeared as `removed`.

## Targeted TASS

```bash
# Your own regions (0-based, half-open BED; .gz ok)
nextflow run jhuapl-bio/taxtriage ... --mask_mode target --mask_bed stec_loci.bed

# Genes picked by name from NCBI feature tables or GFF3
nextflow run jhuapl-bio/taxtriage ... --mask_mode target \
    --mask_feature_table GCF_000008865.2_ASM886v2_feature_table.txt.gz \
    --mask_feature_select '^(stx1a|stx1b|stx2a|stx2b|eae)$'

# Targets, minus anything FCS flags inside them
nextflow run jhuapl-bio/taxtriage ... --mask_mode target --mask_bed panel.bed --fcs
```

`--mask_feature_select` is a case-insensitive regex tested against `gene`/`symbol`, `locus_tag`, `Name`, `ID` and `product`. `--mask_feature_types` (default `CDS,gene`) limits which feature types are considered. If no target matches any aligned reference the run stops with a clear error rather than silently scoring nothing.

BED column 1 must name the reference accession. A quick way to get the accessions in your run: `samtools view -H sample.bam | grep '^@SQ'`.

---

## Parameters

| Parameter               | Default            | Description                                                                                               |
| ----------------------- | ------------------ | --------------------------------------------------------------------------------------------------------- |
| `--fcs`                 | `false`            | Mask NCBI FCS contaminant spans. Downloads NCBI's reports unless `--fcs_details` is given.                |
| `--fcs_details`         | –                  | Comma-separated local `fcs_details_*.txt[.gz]` (or a previous `fcs_subset.txt`). Implies `--fcs`.         |
| `--fcs_summary`         | –                  | Comma-separated `fcs_summary_*.txt[.gz]` for the per-assembly table and organism flag.                    |
| `--fcs_gx_report`       | –                  | Comma-separated FCS-GX reports from your own runs or an assembly's FTP folder.                            |
| `--fcs_adaptor_report`  | –                  | Comma-separated FCS-adaptor reports.                                                                      |
| `--fcs_actions`         | `EXCLUDE,TRIM,FIX` | FCS actions to mask; `all` for every row.                                                                 |
| `--fcs_divisions`       | –                  | Regex on the contaminant division (`prok`, `anml:primates`, `synt`, …).                                   |
| `--fcs_min_len`         | –                  | Ignore FCS spans shorter than this (bp).                                                                  |
| `--fcs_db`              | `<outdir>/fcs_db`  | Where downloaded FCS reports are cached.                                                                  |
| `--mask_bed`            | –                  | Comma-separated BED(s): regions to ignore (exclude) or to score (target).                                 |
| `--mask_fasta`          | –                  | Comma-separated masked FASTA(s) of the references; masked bases are never scored.                         |
| `--mask_fasta_case`     | `hard`             | `hard` = runs of N, `soft` = lowercase, `both`.                                                           |
| `--mask_fasta_min_run`  | `1`                | Ignore masked runs shorter than N bp (e.g. `10` to skip isolated ambiguity Ns).                           |
| `--mask_mode`           | `exclude`          | `exclude` or `target`.                                                                                    |
| `--mask_feature_table`  | –                  | Target mode: NCBI feature tables to select genes from.                                                    |
| `--mask_gff`            | –                  | Target mode: GFF3 files to select genes from.                                                             |
| `--mask_feature_select` | –                  | Target mode: regex on gene / locus_tag / product / Name / ID.                                             |
| `--mask_feature_types`  | `CDS,gene`         | Feature types considered for targets.                                                                     |
| `--mask_untargeted`     | `drop`             | Target mode: drop references with no target from scoring, or `keep` them whole.                           |
| `--mask_read_overlap`   | `0.5`              | Drop a read when this fraction of it is masked. `1.0` = only fully masked reads.                          |
| `--mask_penalty_power`  | `0`                | `TASS × (1 − masked read fraction)^p`. 0 = no penalty.                                                    |
| `--mask_pad`            | `0`                | Pad every interval by N bp.                                                                               |
| `--mask_min_region`     | `0`                | Ignore intervals shorter than N bp.                                                                       |
| `--mask_report`         | `true`             | Write `masking/all.masking.html` when masking is on.                                                      |
| `--mask_compare`        | `false`            | Also score every sample **unmasked** so the masking report shows TASS before/after. Doubles scoring time. |
| `--mask_write_bam`      | `false`            | Also write `masking/bam/<sample>.masked.bam` (+ `.bai`) with the masked alignments removed.               |
| `--mask_report_top`     | `24`               | Coverage tracks drawn in the masking report.                                                              |

---

## Outputs

Everything lands in `<outdir>/masking/`:

| File                                    | Content                                                                                                                                                                            |
| --------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mask.exclude.bed`                      | Regions never scored (FCS spans, plus `--mask_bed` in exclude mode). Column 4 labels the source, e.g. `FCS:FIX:prok:b-proteobacteria:…`. Loads straight into IGV next to the BAMs. |
| `mask.target.bed`                       | Target mode: the targeted regions before FCS subtraction.                                                                                                                          |
| `mask.*.summary.tsv`                    | Per reference: length, region bp and %, interval count, labels, source files.                                                                                                      |
| `mask.*.build.json`                     | Build statistics: rows scanned, rows kept per action, assemblies hit.                                                                                                              |
| `fcs_assemblies.tsv`                    | Per assembly in the run: NCBI FCS source, contaminated fraction/bp/sequences, action counts, divisions, bp masked in this run. `not_screened` = NCBI has no FCS result for it.     |
| `fcs_subset.txt`                        | The FCS rows that matched this run's references - reuse as `--fcs_details`.                                                                                                        |
| `all.masking.html` / `all.masking.tsv`  | The masking report (below).                                                                                                                                                        |
| `unmasked/<sample>.unmasked.paths.json` | With `--mask_compare`: the unmasked baseline scores.                                                                                                                               |
| `bam/<sample>.masked.bam` (+ `.bai`)    | With `--mask_write_bam`: the sample's BAM minus every alignment the mask removed.                                                                                                  |

In every per-sample JSON (`alignment/<sample>.paths.json`):

- `metadata.masking` - `applied` (true when the mask changed what was scored for this sample), mode, BEDs, read overlap, references masked, masked / kept bp, reads considered / removed, penalty power, FCS assembly counts, and `masked_bam` with its kept / dropped alignment counts when `--mask_write_bam` is on. `null` when masking is off.
- Per organism (strain, species, genus): `length` is the **scored** length, plus `length_full`, `masked_bp`, `mask_bp_fraction`, `mask_reads_total`, `mask_reads_masked`, `mask_read_fraction`, `mask_penalty`, and `fcs_assembly` when available.
- Per contig: `length_full`, `masked_bp`, `reads_masked`.
- Per strain: `mask_regions` (the masked or targeted intervals, original coordinates) and `breadth_histogram.mask_b64` (per-bin masked %).

---

## Visualization

**Interactive report (`report/all.odr.html`).** When any sample was masked:

- the header line ends with **⚠ Region masking: …** (FCS contaminant regions excluded / masked regions excluded / TASS from targeted regions only, and how many samples). Hovering it lists each sample's mode, references and bp masked, reads removed, penalty and masked BAM. The PDF export cover carries the same line;
- the Table tab gains **Masked Reads**, **Masked Reads %**, **Full Genome Length (bp)**, **Masked bp**, **Masked bp %** and **Mask Penalty** columns (plus **FCS Assembly** / **FCS Contaminated bp** when an NCBI FCS summary was used). **Genome Length (bp)** and **Breadth %** are over the scored bases. Rows from unmasked samples are blank in these columns;
- Histograms → Breadth of Coverage draws the positional track in real genome coordinates with the masked spans shaded grey (exclude) or the targets shaded amber (target); hovering a bin shows how much of it was masked or targeted.

Reports from unmasked runs look exactly as before.

**Masking report (`masking/all.masking.html`).** A standalone page with:

- per-sample read accounting (references masked, bp masked / kept, reads removed);
- a sortable, filterable organism table - TASS (and, with `--mask_compare`, the unmasked TASS and the change), reads kept / masked, scored and masked bp, breadth, penalty, NCBI FCS assembly flag; detections that disappear after masking are listed as `removed`;
- coverage tracks for the most-affected organisms with the masked or targeted spans drawn in;
- the NCBI FCS assembly table and the per-reference mask contents.

**IGV.** `masking/mask.exclude.bed` (or `mask.target.bed`) loads as a track next to `alignment/*.bam`; with `--mask_write_bam`, load `masking/bam/<sample>.masked.bam` beside the original to see exactly which alignments were removed.

---

## Running the pieces locally

All three scripts live in `bin/` and run outside Nextflow. `build_mask.py` and `mask_report.py` use only the Python standard library (`pysam` is needed only to read references from a BAM).

```bash
# 1. Build a mask for the references in a BAM (or FASTA / .fai / list of accessions)
bin/build_mask.py \
    --fcs_details fcs_details_genbank.txt fcs_details_refseq.txt.gz \
    --fcs_summary fcs_summary_refseq.txt fcs_summary_genbank.txt \
    --mapping sample.combined.gcfmap.tsv \
    --refs sample.dwnld.references.bam \
    --out fcs_mask.bed --summary_out fcs_mask.summary.tsv \
    --assembly_out fcs_assemblies.tsv --subset_out fcs_subset.txt

# 2. Score with it (any existing match_paths.py command line, plus):
bin/match_paths.py ... --mask_bed fcs_mask.bed --fcs_assembly_table fcs_assemblies.tsv
#    targeted:
bin/match_paths.py ... --mask_bed targets.bed --mask_mode target --mask_exclude_bed fcs_mask.bed

# 3. Visualize (baseline = the same samples scored without --mask_bed)
bin/mask_report.py --jsons masked/*.paths.json --baseline unmasked/*.paths.json \
    --assembly_table fcs_assemblies.tsv --mask_summary fcs_mask.summary.tsv --prefix my.masking
```

`build_mask.py --help` lists every source format (FCS details, FCS-GX, FCS-adaptor, NCBI feature tables, GFF3, BED). With no `--refs`, every FCS row is written - millions of intervals for the GenBank file - so always restrict to your references.

---

## Troubleshooting

| Symptom                                                                    | Cause / fix                                                                                                                                                                 |
| -------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `[mask] WARNING: no BAM reference matched any BED sequence name`           | BED column 1 does not use the BAM's accessions. Compare with `samtools view -H`. Version and `NZ_` differences are handled; other renames (e.g. `chr1`) are not.            |
| FCS on, but `mask.exclude.bed` is absent and nothing is masked             | None of the run's references has FCS spans - common for RefSeq. Check `fcs_assemblies.tsv` for assembly-level flags, and use `fcs_details_refseq` for `NC_` references.     |
| `--mask_mode target, but no target interval matched any aligned reference` | The BED / `--mask_feature_select` names nothing that reads aligned to. Check accessions and the regex.                                                                      |
| An organism disappears in target mode                                      | It has no target, so `--mask_untargeted drop` removed it from scoring. Use `keep` to score it whole.                                                                        |
| TASS barely changes with FCS on                                            | Few reads hit the masked spans. Set `--mask_penalty_power 1` if reads in contaminant spans should count against the organism, and `--mask_compare` to see the before/after. |
