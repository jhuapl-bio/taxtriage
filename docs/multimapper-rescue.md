# Ambiguous-Read (Multimapper) Rescue

A read that aligns beautifully to a reference can still be assigned **MAPQ 0**. MAPQ does not measure how good an alignment is — it measures how confident the aligner is that it picked the _right_ reference. When two or more genomes in the reference set sit at 97–98% ANI, a read from any conserved region matches both equally well, the aligner cannot choose, and it reports MAPQ 0 by definition.

That creates a problem for a hard `--minmapq` cut. The reads being discarded are not bad reads; they are reads about which the _reference set_ is ambiguous. On a real sample this can remove three quarters of the true signal:

```
487 reads containing the target accession in the FASTQ
~600 alignments to it post-minimap2 (some false-positive)
120 alignments surviving --minmapq 1
```

The missing ~370 are near-perfect alignments sitting at MAPQ 0 because two or three neighbours in the database are 98% identical.

---

## What rescue does

`--minmapq` remains a hard filter for uniquely-placed reads. What changes is **what gets phred-scored** for a read the aligner marked as a tie. Instead of MAPQ, the alignment is scored on its own observed divergence, on the same phred scale, and held to the same threshold:

$$
\text{aln}_{\text{phred}} = -10 \cdot \log_{10}\!\left(\frac{\text{NM}}{\text{aligned query bases}}\right)
$$

`NM` is the edit distance (mismatches + inserted + deleted bases) and the aligned query length excludes soft- and hard-clipped bases, so the ratio is the observed read-to-reference divergence over the aligned block. A perfect match is capped at Q60. When `NM` is absent, minimap2's gap-compressed divergence tag `de` is used instead; when neither exists the read is **not** rescued.

| Observed divergence | `aln_phred` |
| ------------------- | ----------- |
| 0% (perfect)        | 60 (capped) |
| 1%                  | 20          |
| 2%                  | ~17         |
| 5%                  | ~13         |
| 10%                 | 10          |
| 30%                 | ~5          |

So at the default `--minmapq 5`, an ambiguous read aligning at 2% divergence scores ~Q17 and survives; an ambiguous read aligning at 40% divergence scores ~Q4 and is dropped exactly as it was before.

## What rescue does not do

- **It does not rescue uniquely-placed reads.** A read at MAPQ 1–4 has one home and the aligner still doubts it — that is real evidence of a poor alignment, and it stays filtered. Only alignments at or below `--rescue_max_mapq` (default `0`, the true ties) are eligible.
- **It does not multiply counts across the ANI cluster.** Only the primary record of each read is ever counted; secondary and supplementary alignments are skipped. A read matching three genomes remains **one** observation, not three.
- **It does not invent information.** At 98% ANI, a short read from a conserved region genuinely does not contain enough signal to name its source genome. No BAM setting recovers that. Rescue keeps the read as evidence for the **cluster**; deciding whether one specific accession is present is what the separate metrics below are for.

---

## The gate, in order

Each check is evaluated only for reads that already failed the MAPQ cut, using C-level pysam attributes plus one or two tag lookups, so the cost on a multi-million-read BAM is negligible and the BAM is streamed exactly **once**.

1. `MAPQ <= --rescue_max_mapq` **and** the low MAPQ looks ambiguous (MAPQ 0, or a bwa `XA` tag).
2. Aligned (unclipped) fraction of the read ≥ `--rescue_min_aln_frac`. Guards against a conserved fragment anchoring a mostly-clipped read.
3. Aligned query length ≥ `--rescue_min_aln_len`.
4. Properly paired, if `--rescue_require_proper_pair` is set (paired-end only; single-end reads are unaffected).
5. `aln_phred` ≥ `--rescue_min_aln_phred`, plus the optional absolute ceiling `--rescue_max_nm_rate`.

The aligned-fraction default follows `--platform`, so Illumina paired-end, single-end Illumina, and ONT all work without tuning:

| Platform | Default `rescue_min_aln_frac` |
| -------- | ----------------------------- |
| Illumina | 0.95                          |
| ONT      | 0.80                          |
| PacBio   | 0.85                          |

---

## Reading the result

Every reference in `*.paths.json` now reports the two kinds of evidence separately, so a claim can be qualified rather than flattened into one number:

| Field                    | Meaning                                                                  |
| ------------------------ | ------------------------------------------------------------------------ |
| `numreads`               | All kept reads (unique + rescued)                                        |
| `numreads_unique`        | Reads with MAPQ ≥ `--minmapq` — evidence favouring **this** accession    |
| `numreads_rescued`       | Strong ambiguous reads — evidence for the **ANI cluster**                |
| `rescued_fraction`       | `numreads_rescued / numreads`                                            |
| `mean_rescued_aln_phred` | Mean alignment phred of the rescued reads                                |
| `highmapq_fraction`      | Fraction of mapped reads at MAPQ ≥ `--minmapq` (scales breadth and Gini) |

Rescued reads deliberately do **not** count toward `highmapq_fraction`, so cluster-level evidence cannot masquerade as accession-specific evidence in the [TASS score](tass-scoring.md). `--rescue_counts_as_highmapq` flips that if you want breadth and Gini computed over the rescued set too.

A defensible reading of a high-ANI cluster:

| Evidence                                  | Interpretation                                      |
| ----------------------------------------- | --------------------------------------------------- |
| High `numreads_unique`                    | Supports this specific accession                    |
| High `rescued_fraction`, low unique count | Supports the species complex, not one genome        |
| Broad, even coverage (low Gini)           | More convincing than a pile-up over conserved genes |
| Target-specific SNPs / unique regions     | The strongest accession-level evidence available    |

A reference whose reads are almost entirely rescued is a cluster-level call. Report it at the species or species-complex level and use the unique reads — and their genome-wide distribution — to decide whether that one accession is itself supported.

---

## Parameters

| Parameter                        | Default | Description                                                                                     |
| -------------------------------- | ------- | ----------------------------------------------------------------------------------------------- |
| `--rescue_multimapped`           | `true`  | Enable the rescue. Use `--rescue_multimapped false` for a hard `--minmapq` cut on every read.   |
| `--rescue_min_aln_phred <float>` | `null`  | Phred bar for the alignment itself. Null uses `--minmapq`, so one threshold governs both cases. |
| `--rescue_max_mapq <int>`        | `0`     | Only alignments at or below this MAPQ are eligible. `0` = true ties only.                       |
| `--rescue_min_aln_frac <float>`  | `null`  | Minimum unclipped fraction of the read. Null uses the platform preset above.                    |
| `--rescue_max_nm_rate <float>`   | `null`  | Optional extra ceiling on NM per aligned base (e.g. `0.02`).                                    |
| `--rescue_min_aln_len <int>`     | `50`    | Minimum aligned query length in bp.                                                             |
| `--rescue_require_proper_pair`   | `false` | Paired-end only: restrict rescues to reads flagged `0x2`.                                       |
| `--rescue_counts_as_highmapq`    | `false` | Count rescued reads toward `highmapq_fraction` (breadth/Gini scaling).                          |

### Examples

Default behaviour — one threshold, applied to MAPQ for unique reads and to alignment phred for ties:

```bash
nextflow run jhuapl-bio/taxtriage --input samplesheet.csv --minmapq 5 ...
```

Stricter: only ambiguous reads at ≤1% divergence count, and they must be properly paired:

```bash
nextflow run jhuapl-bio/taxtriage --input samplesheet.csv \
    --minmapq 5 --rescue_min_aln_phred 20 --rescue_require_proper_pair ...
```

Disable entirely (pre-existing hard-cut behaviour):

```bash
nextflow run jhuapl-bio/taxtriage --input samplesheet.csv --rescue_multimapped false ...
```

### Checking one accession by hand

The `samtools` equivalent of the default gate, for a quick sanity check against a known ground-truth accession:

```bash
samtools view -c -F 0x904 \
  -e 'mapq>=5 || (mapq==0 && (qlen-sclen)>=0.95*qlen && [NM]<=0.02*(qlen-sclen))' \
  aligned_sorted.bam NZ_CP076232
```

`-F 0x904` drops unmapped, secondary and supplementary records, so each read is counted once.

---

## See also

- [TASS Scoring](tass-scoring.md) — how `highmapq_fraction` scales breadth and Gini
- [Detection Rescue](detection-rescue.md) — a different mechanism: re-surfacing below-cutoff organisms in the report
- [CLI Parameters](cli-parameters.md)
