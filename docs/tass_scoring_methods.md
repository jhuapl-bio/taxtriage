# TASS Scoring Methodology

## Complete Pipeline: From Raw Sequencing Data to a Final TASS Score

---

## 1. Pipeline Overview

The TASS (Threat Agnostic Sentinel Surveillance) pipeline takes raw sequencing alignment data and produces a confidence score between 0 and 1 for each detected organism. A higher score means we're more confident the organism is genuinely present in the sample.

The pipeline scripts that contribute to the TASS scoring mentioned below are available [here](../bin):

1. **Input parsing & per-reference stats** — Pull alignment metrics from BAM, depth/bedgraph, and FASTA files (`match_paths.py`)
2. **Conflict resolution & signature comparison** — Find and handle cases where reads map ambiguously to multiple organisms (`conflict_regions.py`)
3. **Score computation** — Calculate individual features and combine them into a final TASS score (`optimize_weights.py`)

### Input Files

| File Type                                    | What it's for                                                                                      |
| -------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| **BAM** (`.bam`)                             | Aligned reads against the reference database; source of read counts, mapping quality, and coverage |
| **BEDGRAPH / Depth** (`.bedgraph`, `.depth`) | Per-position coverage depth for each reference                                                     |
| **FASTA** (`.fasta`, `.fa`)                  | Reference genome sequences; used to build k-mer signatures for conflict detection                  |
| **Match file** (`.tsv`)                      | Maps accession IDs to taxonomy (taxid, organism name, description)                                 |
| **Kraken2 report** (optional)                | Independent taxonomic classification for cross-checking                                            |
| **HMP distributions** (optional)             | Human Microbiome Project abundance data for body-site context                                      |

---

## 2. Stage 1: Input Parsing and Per-Reference Stats

### 2.1 BAM Processing (`match_paths.py`)

The BAM file is read using `pysam`. For each reference, we extract:

- **`numreads`** — Number of primary, non-supplementary alignments
- **`meanmapq`** — Average mapping quality across all aligned reads
- **`highmapq_fraction`** — Fraction of reads with MAPQ ≥ threshold (e.g., MAPQ ≥ 5)
- **`covered_regions`** — List of `(start, end, depth)` tuples from the bedgraph

### 2.2 Coverage from Depth/Bedgraph

The bedgraph file is parsed into a table with columns `(chrom, start, end, depth)`. For each reference:

```
covered_bases = sum of (end - start + 1) for every region where depth > 0

coverage = covered_bases / reference_length
```

### 2.3 Abundance Metrics (RPM / RPKM)

Computed in `compute_scores_per()`. These normalize read counts by sequencing depth and genome length:

```
RPM  = numreads / (total_reads / 1,000,000)

RPKM = numreads / ((total_reads / 1,000,000) × (reference_length / 1,000))
```

### 2.4 MAPQ Normalization

From `utils.py::normalize_mapq()`. This rescales mapping quality to a 0–1 range:

```
mapq_score = (mapq - min_mapq) / (max_mapq - min_mapq)

where max_mapq = 60, min_mapq = 0
```

The result is clipped to [0, 1] to handle any out-of-range values.

---

## 3. Stage 2: Conflict Region Detection (`conflict_regions.py`)

### 3.1 What This Stage Does

Closely related organisms (e.g., _E. coli_ vs _Shigella_) share large regions of their genomes. This means reads from one organism can legitimately align to the other — a problem called **cross-mapping**. Left unchecked, cross-mapping inflates the apparent confidence for organisms that aren't really there.

This stage finds those shared regions and removes a proportional share of ambiguous reads before scoring.

### 3.2 Region Merging

Raw bedgraph intervals are merged into larger contiguous regions using `merge_bedgraph_regions()`. The default strategy (`"jump"`) groups adjacent intervals when:

1. The depth change between consecutive intervals is ≤ a threshold (set automatically from the 97.5th percentile of depth)
2. The gap between intervals is ≤ 10% of the reference length (configurable)
3. The merged group doesn't exceed `max_group_size` intervals

The result is a smaller set of meaningful regions per reference, each described by `(chrom, start, end, mean_depth)`.

### 3.3 Signature Generation

For each merged region, a **sourmash MinHash signature** is created from the reference FASTA sequence:

```python
mh = MinHash(n=0, ksize=kmer_size, scaled=scaled)   # default ksize=31, scaled=100
mh.add_sequence(fasta.fetch(chrom, start, end))
sig = SourmashSignature(mh, name="ref:start-end")
```

- If FASTA files are available, the sequence is pulled directly from them.
- Regions shorter than `kmer_size` are skipped (can't build a valid signature).

### 3.4 SBT Index and Search

Signatures are organized into a **Sequence Bloom Tree (SBT)** index, then each signature is searched against it:

```python
for search_result in sbt.search(sig_query, threshold=min_threshold):
    jaccard = search_result.score    # Jaccard similarity between two regions
```

Any pair of regions from _different_ references with Jaccard similarity ≥ `min_threshold` (default 0.2) is flagged as a potential conflict.

### 3.5 Conflict Group Construction

`build_conflict_groups()` builds a graph where:

- **Nodes** are genomic regions `(ref, start, end)`
- **Edges** connect regions from different references that are similar enough (Jaccard ≥ threshold)

Connected components (found via BFS) define **conflict groups** — clusters of regions across multiple organisms that share significant k-mer content.

### 3.6 Proportional Read Removal

For each conflict group, `finalize_proportional_removal()` handles the cleanup:

- Reads overlapping any region in the conflict group are collected
- A proportional subset is randomly removed (mode: `'random'`) to reduce cross-mapping inflation
- Removed read IDs are logged per reference

### 3.7 Comparison Metrics (Δ Breadth)

After read removal, `compare_metrics()` calculates how much each reference changed:

```
Δ All     = total_reads - pass_filtered_reads          (number of reads removed)
Δ All%    = 100 × (reads_removed / total_reads)        (% of reads removed)
Δ⁻¹ Breadth = breadth_new / breadth_old               (how much breadth remained after removal)
```

**`Δ⁻¹ Breadth`** and **`Δ All%`** are passed directly into the minhash score (see Section 7).

### 3.8 Shared-Window Mode (Optional)

When `compare_to_reference_windows=True`, an alternative approach is used:

1. Fixed-size windows (default: 2000 bp, step: 500 bp) are tiled across each reference FASTA
2. Each window is sketched into a sourmash signature
3. Windows are compared across all FASTAs to find shared regions
4. A scoring heuristic determines which reads to drop per conflict

### 3.9 ANI Matrix Generation

`generate_ani_matrix()` builds an organism-level Average Nucleotide Identity (ANI) matrix:

1. Composite signatures are built per taxid from all its accessions
2. Pairwise ANI is estimated from MinHash containment between organism signatures
3. The result is written to `organism_ani_matrix.csv` as a symmetric matrix

---

## 4. Stage 3: Score Computation (`optimize_weights.py`)

All component scores are computed in `compute_scores_per()` and combined in `compute_tass_score()`.

---

## 5. Component 1: Breadth Log Score (Genome Coverage)

### 5.1 What It Measures

The breadth score captures how much of the reference genome has at least one read aligned to it. If reads are spread across the genome, that's a good sign the organism is genuinely present. If reads pile up in just one spot, it's more likely cross-mapping or contamination.

This fraction is passed through a sigmoid curve to produce a 0–1 confidence value.

### 5.2 Sigmoid Function

Defined in `breadth_score_sigmoid()`:

$$
\text{breadth\_sigmoid}(c) = \frac{1}{1 + e^{-s \cdot (c - m)}}
$$

Where:

- $c$ = coverage fraction (covered_bases / reference_length), range [0, 1]
- $m$ = midpoint — the coverage level at which the sigmoid returns 0.5 (default: 0.01, i.e., 1%)
- $s$ = steepness — controls how quickly the score jumps from 0 to 1 (default: 12,000)

**Overflow protection:** If $s \cdot (c - m) \geq 50$, the function returns 1.0. If $\leq -50$, it returns 0.0.

### 5.3 MAPQ Scaling

Low mapping quality (MAPQ ≈ 0) suggests reads are aligning ambiguously — they may not truly belong to this organism. The breadth sigmoid is scaled down when most reads have low MAPQ:

$$
\text{mapq\_scale} = (\text{highmapq\_fraction})^{p}
$$

Where $p$ = `mapq_breadth_power` (default: 2.0).

### 5.4 Final Breadth Log Score

$$
\boxed{\text{breadth\_log\_score} = \text{breadth\_sigmoid}(c, m, s) \times (\text{highmapq\_fraction})^{p}}
$$

**Example:** Coverage = 15%, midpoint = 0.01, steepness = 12,000, 90% of reads have high MAPQ:

```
sigmoid(0.15) ≈ 1.0
mapq_scale = 0.9² = 0.81
breadth_log_score = 1.0 × 0.81 = 0.81
```

### 5.5 Tunable Parameters by Sample Type

| Sample Type   | Midpoint | Steepness | Why                                                |
| ------------- | -------- | --------- | -------------------------------------------------- |
| Standard      | 0.01     | 12,000    | 1% coverage = moderate confidence                  |
| Sterile/Blood | 0.0005   | 50,000    | Even 0.05% coverage is meaningful in sterile sites |

---

## 6. Component 2: Gini Coefficient (Coverage Uniformity)

### 6.1 What It Measures

The Gini score measures how evenly reads are spread across the genome. If an organism is truly present, reads should appear across many different genomic positions — not just cluster in one conserved region.

The pipeline inverts the classical Gini coefficient so that **uniform coverage → high score** and **clumped coverage → low score**.

### 6.2 Step-by-Step Computation

The full computation lives in `getGiniCoeff()` and has five sub-steps.

#### Step 1: Build a Transformed Coverage Histogram

Each region's depth is log-transformed to reduce the effect of extreme outliers:

$$
d'_i = \log_{10}(1 + d_i)
$$

A histogram maps each transformed depth value to the number of bases at that depth. Uncovered bases are recorded at $\log_{10}(1 + 0) = 0$.

#### Step 2: Compute Raw Gini from the Lorenz Curve

The Gini coefficient summarizes inequality in coverage:

1. Sort coverage values in ascending order
2. Build cumulative Lorenz points: $(x_i, y_i)$ where $x_i$ = cumulative fraction of bases, $y_i$ = cumulative fraction of total coverage
3. Compute the area under the Lorenz curve via trapezoidal integration:

$$
A_L = \sum_{i} \frac{(x_{i+1} - x_i)(y_i + y_{i+1})}{2}
$$

$$
G_{\text{raw}} = 1 - 2 A_L
$$

$G_{\text{raw}} \in [0, 1]$: 0 = perfectly uniform, 1 = maximally unequal.

#### Step 3: Invert and Scale the Raw Gini

Since low Gini (uniform coverage) is what we want, we invert it:

$$
G_{\text{transformed}} = \text{clamp}\!\Big(\alpha \cdot \sqrt{1 - G_{\text{raw}}},\; 0,\; 1\Big)
$$

Where $\alpha = 1.8$ (default). Examples:

- $G_{\text{raw}} = 0$ (perfectly uniform) → $G_{\text{transformed}} = \min(1.8, 1.0) = 1.0$
- $G_{\text{raw}} = 0.9$ (very clumped) → $G_{\text{transformed}} = 1.8 \times \sqrt{0.1} \approx 0.57$

#### Step 4: Length-Based Scaling Factor

Larger genomes need more reads to achieve uniform coverage. A log-scale bonus rewards organisms with bigger genomes:

$$
S_{\text{length}} = 1 + R \cdot \log_{10}\!\Big(\max\!\big(\frac{\min(L, L_{\max})}{L_{\text{base}}},\; 1\big)\Big)
$$

Where:

- $L$ = genome length
- $L_{\text{base}}$ = baseline length (default: 500,000 bp)
- $L_{\max}$ = length cap (default: 10⁹ bp)
- $R$ = reward factor (default: 2)

#### Step 5: Positional Dispersion Factor

This measures how spread out the covered regions are across the genome. Even if depth is uniform, we want reads scattered in different places, not all in one block.

$$
\bar{m} = \frac{1}{n}\sum_{i=1}^{n} \frac{s_i + e_i}{2}
\qquad
\sigma^2 = \frac{1}{n}\sum_{i=1}^{n}(m_i - \bar{m})^2
\qquad
D = \sqrt{\frac{\sigma^2}{L^2 / 12}}
$$

$L^2/12$ is the maximum variance of a uniform distribution over $[0, L]$. $D \in [0, 1]$: higher means more spatially spread out.

#### Step 6: Final Gini Score

$$
\boxed{\text{gini\_coefficient} = \min\!\Big(1.0,\;\; G_{\text{transformed}} \times S_{\text{length}} \times (1 + \beta \cdot D)\Big)}
$$

Where $\beta$ = dispersion weight (default: 0.5).

**Example:** Moderate uniformity ($G_{\text{raw}} = 0.3$), genome length 5 Mbp, well-spread regions ($D = 0.7$):

```
G_transformed  = 1.8 × √(0.7) ≈ 1.0 (clamped to 1.0)
S_length       = 1 + 2 × log10(5,000,000 / 500,000) = 1 + 2 × 1 = 3.0
Dispersion     = 1 + 0.5 × 0.7 = 1.35
Final          = min(1.0, 1.0 × 3.0 × 1.35) = 1.0 (clamped to 1.0)
```

---

## 7. Component 3: Minhash Reduction (Signature-Based Uniqueness)

### 7.1 What It Measures

The minhash component asks: _how much did this organism's signal change after we removed cross-mapped reads?_

An organism that loses most of its reads during conflict resolution probably wasn't genuinely present — those reads belonged to something else. An organism that retains its reads is more likely real.

### 7.2 Raw Minhash Score

The conflict detection pipeline produces per-reference metrics at the species level. Two values are key:

- **`Δ⁻¹ Breadth`** ($B_r$): Ratio of post-removal breadth to pre-removal breadth. 1.0 = no breadth lost; 0.5 = half the breadth was from cross-mapped reads.
- **`Δ All%`** ($\Delta\%$): Percentage of reads removed during conflict resolution.

The raw minhash score penalizes organisms that lost a lot of reads:

$$
\text{penalty} = \frac{1}{1 + e^{-k \cdot (\Delta\% - x_0)}}
\qquad k = 0.90,\; x_0 = -10.0
$$

$$
\text{minhash\_score} = \min(1.0,\;\; B_r \times \text{penalty})
$$

**Fallback (no comparison data available):**

$$
\text{minhash\_score}_{\text{fallback}} = \text{rpm\_confidence} \times 0.5
$$

Where `rpm_confidence` is a sigmoid over the fraction of total reads:

$$
\text{rpm\_confidence} = \frac{1}{1 + e^{-50000 \cdot (f_{\text{reads}} - 0.0001)}}
$$

### 7.3 Confidence Gating

Before the minhash score is used, it's gated by a coverage confidence factor. This prevents organisms with very little coverage from getting inflated minhash scores:

$$
\text{conf} = w_b \cdot \text{breadth\_sigmoid}(c) + w_g \cdot G_{\text{score}}
$$

By default $w_b = 1.0$ and $w_g = 0.0$, so gating is based purely on breadth.

$$
\boxed{\text{minhash\_reduction} = \text{clamp}(\text{conf},\; 0,\; 1)}
$$

> **Note:** The value stored as `minhash_reduction` is the _confidence_ value, not the raw minhash score itself. This means the minhash component represents how confident we are that the organism's alignment signal is real — primarily driven by whether the organism has meaningful genome breadth.

---

## 8. Additional Scoring Components

### 8.1 HMP Percentile (Body-Site Abundance Context)

This compares the observed abundance to what the Human Microbiome Project found for this organism at the relevant body site. If an organism is at typical abundance for that site, it gets a neutral score. Unusually high or low abundance is flagged.

$$
f_{\text{obs}} = \frac{\text{numreads}}{\text{total\_reads}}
\qquad
z = \frac{f_{\text{obs}} - \mu_{\text{HMP}}}{\sigma_{\text{HMP}}}
\qquad
P_{\text{base}} = \Phi(z)
$$

The base percentile is then adjusted based on how far from normal the observation is:

| Z-score range          | Adjusted percentile               | What it means                     |
| ---------------------- | --------------------------------- | --------------------------------- |
| $z < -1.0$             | $P_{\text{base}}^2$               | Penalize — below typical          |
| $-1.0 \leq z \leq 1.0$ | $P_{\text{base}}$                 | Neutral — near expected           |
| $1.0 < z \leq 2.0$     | $1 - (1 - P_{\text{base}})^{1.5}$ | Boost — above typical             |
| $z > 2.0$              | $1 - (1 - P_{\text{base}})^{2}$   | Strong boost — well above typical |

### 8.2 Plasmid Score

Within each species, strains are compared by their plasmid coverage. Computed in `match_paths.py`.

**Absolute quality** (does this plasmid have real coverage?):

$$
Q = \sqrt{B_{\text{plasmid}} \times G_{\text{plasmid}}}
$$

Where $B_{\text{plasmid}}$ = `breadth_score_sigmoid(plasmid_coverage)` and $G_{\text{plasmid}}$ = `getGiniCoeff(plasmid_regions)`. The geometric mean requires both breadth and uniformity to be decent — one alone isn't enough.

**Relative disparity** (how does this strain compare to others in the same species?):

$$
D_{\text{rel}} = \begin{cases}
\min\!\big(1.0,\; 0.7 \cdot \frac{c_{\text{strain}}}{c_{\max}} + 0.3 \cdot \frac{r_{\text{strain}}}{r_{\max}}\big) & \text{if multiple strains have plasmids} \\
1.0 & \text{if this is the only strain with a plasmid}
\end{cases}
$$

**Final plasmid score:**

$$
\text{plasmid\_score} = \min(1.0,\;\; Q \times D_{\text{rel}})
$$

Strains with no plasmid accessions receive `plasmid_score = 0`.

### 8.3 Low-Abundance Confidence

A sigmoid in log₁₀-RPM space that rewards organisms with meaningful RPM even when absolute read counts are low. Particularly useful for sterile or blood-site samples.

$$
\text{RPM} = \frac{\text{numreads}}{\text{total\_reads}} \times 10^6
$$

$$
\text{abundance\_confidence} = \frac{1}{1 + e^{-s \cdot (\log_{10}(\text{RPM}) - \log_{10}(m))}}
$$

Where $s$ = steepness (default: 2.0) and $m$ = RPM midpoint (default: 5.0).

### 8.4 K2 Disparity Score

Compares the pipeline's species-level classification against Kraken2's independent classification. If Kraken2 assigns reads to a different genus, this reduces confidence. **Disabled by default** (`k2_disparity_score_weight = 0`).

### 8.5 DIAMOND Identity

Average amino acid identity from a DIAMOND BLASTX protein alignment, weighted by CDS count. Provides protein-level confirmation of the taxonomic assignment. **Disabled by default** (`diamond_identity = 0`).

---

## 9. Final TASS Score Formula

### 9.1 The Core Weighted Sum

Computed in `compute_tass_score()`:

$$
\boxed{
\text{TASS}_{\text{base}} = \sum_{i} w_i \cdot x_i
}
$$

Expanded:

$$
\text{TASS}_{\text{base}} = w_b \cdot \text{breadth\_log\_score}
\;+\; w_m \cdot \text{minhash\_reduction}
\;+\; w_g \cdot \text{gini\_coefficient}
\;+\; w_h \cdot \text{hmp\_percentile}
$$

$$
\;+\; w_d \cdot \text{disparity}
\;+\; w_q \cdot \text{mapq\_score}
\;+\; w_{k2} \cdot \text{k2\_disparity\_score}
\;+\; w_{di} \cdot \text{diamond\_identity}
\;+\; w_{ac} \cdot \text{abundance\_confidence}
$$

### 9.2 Default Weights

| Component              | Weight   | CLI Flag                      |
| ---------------------- | -------- | ----------------------------- |
| `breadth_log_score`    | **0.26** | `--breadth_weight`            |
| `minhash_reduction`    | **0.29** | `--minhash_weight`            |
| `gini_coefficient`     | **0.45** | `--gini_weight`               |
| `hmp_percentile`       | 0.00     | `--hmp_weight`                |
| `disparity`            | 0.00     | `--disparity_weight`          |
| `mapq_score`           | 0.00     | `--mapq_score`                |
| `k2_disparity_score`   | 0.00     | `--k2_disparity_score_weight` |
| `diamond_identity`     | 0.00     | `--diamond_identity`          |
| `plasmid_bonus_weight` | 0.19     | `--plasmid_bonus_weight`      |

**Why don't the weights add up to 1.0?** If the three primary weights (breadth=0.40, minhash=0.55, gini=0.15) sum to 1.0 when normalized, and plasmid adds another 0.19. This is intentional — the final score is clamped to [0, 1] at the end, so overshooting is fine. It allows components to reinforce each other when the support for an organism presence is good/high.

### 9.3 Additive Modifiers

**Plasmid bonus** (applied after the base sum):

$$
\text{TASS}_1 = \text{TASS}_{\text{base}} + w_{\text{plasmid}} \times \text{plasmid\_score}
$$

Where $w_{\text{plasmid}}$ = `plasmid_bonus_weight` (default: 0.0 — disabled).

### 9.4 Multiplicative Abundance Gate (Optional)

When `--abundance_gate` is enabled, the full score is multiplied by the abundance confidence sigmoid. This collapses scores for organisms with very low read fractions:

$$
\text{TASS}_2 = \text{TASS}_1 \times \text{abundance\_confidence}
$$

**Note:** This is exclusive with the additive `plasmid_bonus_weight` approach.

### 9.5 Power Transform (Optional)

When `score_power` ≠ 1.0, a power transform reshapes the score distribution:

$$
\text{TASS}_3 = (\text{clamp}(\text{TASS}_2, 0, 1))^{p}
$$

Where $p$ = `score_power`. Some examples:

- $p = 0.5$: score 0.09 → 0.30, score 0.95 → 0.97 (lifts low scores)
- $p = 0.3$: score 0.09 → 0.52, score 0.95 → 0.98 (aggressive lift)
- $p = 1.0$: no change (default)

### 9.6 Final Clamping

$$
\boxed{\text{TASS}_{\text{final}} = \text{clamp}(\text{TASS}_3,\; 0,\; 1)}
$$

---

## 10. Full Pipeline Sequence Diagram

```
BAM file ──────────────────┐
                           ▼
                    ┌──────────────┐
BEDGRAPH/depth ───▶│ match_paths  │──▶ per-reference stats:
                   │   .py        │    numreads, meanmapq, coverage,
Match file ───────▶│              │    covered_regions, highmapq_fraction
                   └──────┬───────┘
                          │
          ┌───────────────┼───────────────┐
          ▼               ▼               ▼
   ┌─────────────┐ ┌────────────┐  ┌────────────────┐
   │ conflict_   │ │ optimize_  │  │ map_taxid.py   │
   │ regions.py  │ │ weights.py │  │ body_site_     │
   │             │ │            │  │ normalization   │
   │ • merge     │ │ • gini     │  │ • taxonomy     │
   │   bedgraph  │ │ • breadth  │  │   lookup       │
   │ • sketch    │ │ • minhash  │  │ • pathogen     │
   │   signatures│ │ • HMP      │  │   classification│
   │ • SBT search│ │ • plasmid  │  └────────────────┘
   │ • conflict  │ │ • TASS     │
   │   groups    │ │   scoring  │
   │ • read      │ └─────┬──────┘
   │   removal   │       │
   └──────┬──────┘       │
          │              │
          ▼              ▼
   comparison_df    TASS scores
   (Δ Breadth,     per organism
    Δ All%)              │
          │              ▼
          └────────▶ Final TSV Report
```

---

### Extra

#### Example command line for running TASS calculations:

To run `match_paths.py` you will need 4 files minimum along with a directory of taxonomy ([taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz)). Make sure you're in the `bin` folder and adjust paths as needed. The required files are:

1. BAM file with reads aligned to the reference database
2. FASTA file with reference sequences (optional but improves conflict detection)
3. Bedgraph file with coverage depth per position for each reference
4. Match file mapping accessions to taxids and organism names

5. Alignment files can be generated with standard tools like [Minimap2](https://github.com/lh3/minimap2) or [Bowtie2](https://github.com/BenLangmead/bowtie2), using the same reference database that TASS will use for scoring. Example command can be:

```bash
minimap2 \
     -ax sr \
    --split-prefix Miseq_Run_A.Miseq_Run_A.dwnld.references.prefix \
    Miseq_Run_A.dwnld.references.fasta \
    Miseq_Run_A_1.fastp.fastq.gz Miseq_Run_A_2.fastp.fastq.gz \
    -L \
    -a | samtools sort | samtools view -b -h -o Miseq_Run_A.Miseq_Run_A.dwnld.references.bam
```

2. Ensure (for speed reasons) that the reference FASTA file is indexed with `samtools faidx` and that the BAM file is indexed with `samtools index` before running the TASS pipeline. You can run `samtools faidx` and `samtools index` with:

```bash
samtools faidx Miseq_Run_A.dwnld.references.fasta
samtools index Miseq_Run_A.Miseq_Run_A.dwnld.references.bam
```

3. The bedgraph can be generated with `bedtools genomecov` or similar. Example command:

```bash
bedtools genomecov -ibam Miseq_Run_A.Miseq_Run_A.dwnld.references.bam -bg > Miseq_Run_A.bedgraph
```

4. The match file should be a TSV with columns: `accession`, `taxid`, `organism_name`, `description`. This can be generated from the reference FASTA headers and a taxonomy lookup. Example:

| Acc         | Mapped_Value |
| ----------- | ------------ |
| NC_000913.3 | 511145       |
| NC_000964.3 | 224308       |
| NC_001781.1 | 11250        |
| NC_001803.1 | 208893       |

You can also provide custom names with the `--namecol` option if you don't want the taxdump names to be used by default. You'll need to have the column added to the match file and specify the column name like `--namecol 2`. Example:
| Acc | Mapped_Value | Custom_Name |
|---|---|---|
| NC_000913.3 | 511145 | Escherichia coli str. K-12 substr. MG1655 |
| NC_000964.3 | 224308 | Bacillus subtilis subsp. subtilis str. 168 |
| NC_001781.1 | 11250 | Human respiratory syncytial virus |
| NC_001803.1 | 208893 | Pseudomonas aeruginosa PAO1
|

Finally, the TASS scoring can be run with:

```bash
match_paths.py \
       -i Miseq_Run_A.Miseq_Run_A.dwnld.references.bam \
       -o Miseq_Run_A.paths.json \
       -b Miseq_Run_A.bedgraph \
       -s Miseq_Run_A \
       -t nasal \
       --output_dir search_results \
       -f Miseq_Run_A.dwnld.references.fasta \
       -p pathogen_sheet.csv \
       -m Miseq_Run_A.merged.taxid.tsv \
       --k2 Miseq_Run_A.filtered.report \
       --taxid_removal_stats \
       --rank genus \
       --minmapq 5 \
       --compare_references \
       --taxdump taxdump \
       --platform ILLUMINA \
       --thresholds_json sampletype_best_thresholds.json
```

:warning: The above command is an example and may require adjustments based on your specific file paths, sample types, and desired parameters.

The output will be in JSON format, which can be generated downstream as ODR and TSV with:

```
create_report.py \
  -i Miseq_Run_A.paths.json \
  -u Miseq_Run_A.organisms.report.txt \
  -o Miseq_Run_A.organisms.report.pdf \
  --show_potentials \
  --rank genus
```

## 11. Aggregation Levels

Scores are computed at multiple taxonomic levels:

1. **Accession level** — Raw per-reference statistics
2. **Strain level** (key) — Accessions grouped by strain; plasmid accessions are tagged and scored separately
3. **Species level** (subkey) — Aggregated across strains; minhash comparison operates at this level
4. **Genus level** (toplevelkey) — Highest aggregation; used for HMP lookups and final reporting

At each level, metrics are either summed (`numreads`, `covered_bases`), averaged (MAPQ), or re-computed from scratch (`gini`, `breadth` from merged regions).

---

## 12. Weight Optimization (`ground_truth.py`)

When ground-truth labels are available (e.g., simulated data where read names encode their true origin), the pipeline can automatically optimize scoring weights using:

1. **Differential Evolution** — broad search across the full weight space
2. **SLSQP** — local refinement starting from the best solution found above

The objective is to maximize the score gap between true positives and false positives, while keeping weights within valid ranges. Optimized weights are written back into the pipeline arguments for the final scoring pass.

---

## 13. Summary of Key Equations

### Breadth Log Score

$$
\text{breadth\_log\_score} = \frac{1}{1 + e^{-s(c-m)}} \times (\text{highmapq\_fraction})^p
$$

### Gini Coefficient Score

$$
\text{gini\_coefficient} = \min\!\Big(1,\; \alpha\sqrt{1-G} \;\times\; \big(1 + R\log_{10}\tfrac{L}{L_b}\big) \;\times\; (1+\beta D)\Big)
$$

### Minhash Reduction

$$
\text{minhash\_reduction} = \text{clamp}\!\big(w_b \cdot \text{breadth\_sigmoid}(c) + w_g \cdot G_{\text{score}},\; 0,\; 1\big)
$$

### TASS Score

$$
\text{TASS} = \text{clamp}\!\bigg(\Big[\sum_i w_i x_i + w_{\text{plasm}} \cdot P\Big] \times \text{gate}^{[ab\_gate]}\bigg)^{p}
$$

Where $\text{gate}^{[ab\_gate]}$ = plasmid_weight if the gate is enabled, or 1.0 if disabled. $p$ = score_power (1.0 if disabled).
