# TASS Scoring Methodology

## Complete Pipeline: From Raw Sequencing Data to a Final TASS Score

---

See [notations](#notations) for info on the mathematical symbols used in the formulas below.

## 1. Pipeline Overview

The TASS (Threat Agnostic Sentinel Surveillance) pipeline takes raw sequencing alignment data and produces a confidence score between 0 and 1 for each detected organism. A higher score means we're more confident the organism is genuinely present in the sample.

The pipeline scripts that contribute to the TASS scoring mentioned below are available [here](../bin):

1. **Input parsing & per-reference stats** — Pull alignment metrics from BAM, depth/bedgraph, and FASTA files (`match_paths.py`)
2. **Conflict resolution & signature comparison** — Find and handle cases where reads map ambiguously to multiple organisms (`conflict_regions.py`)
3. **Score computation** — Calculate individual features and combine them into a final TASS score (`optimize_weights.py`)

### Input Files

| File Type                                    | What it's for                                                                                                                                                     |
| -------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **BAM** (`.bam`)                             | Aligned reads against the reference database; source of read counts, mapping quality, and coverage                                                                |
| **BEDGRAPH / Depth** (`.bedgraph`, `.depth`) | Per-position coverage depth for each reference. **Optional** when `--compare_references` is active (coverage is derived from the BAM index instead of a bedgraph) |
| **FASTA** (`.fasta`, `.fa`)                  | Reference genome sequences; used to build k-mer signatures for conflict detection. Required when `--compare_references` is active                                 |
| **Match file** (`.tsv`)                      | Maps accession IDs to taxonomy (taxid, organism name, description)                                                                                                |
| **Kraken2 report** (optional)                | Independent taxonomic classification for cross-checking                                                                                                           |
| **HMP distributions** (optional)             | Human Microbiome Project abundance data for body-site context                                                                                                     |

---

## 2. Stage 1: Input Parsing and Per-Reference Stats

### 2.1 BAM Processing (`match_paths.py`)

The BAM file is read using `pysam`. For each reference, we extract:

- **`numreads`** — Number of primary, non-supplementary alignments
- **`meanmapq`** — Average mapping quality across all aligned reads
- **`highmapq_fraction`** — Fraction of reads with MAPQ ≥ threshold (e.g., MAPQ ≥ 5)
- **`covered_regions`** — List of `(start, end, depth)` from the bedgraph

\*_Note:_ When `--compare_references` is active, coverage is estimated from BAM index statistics instead of a bedgraph file (see Section 3.8). Additionally, we will use the term _contig_ interchangeably with reference accessions be they scaffolds, plasmids, contigs, or full chromosomes.

### 2.2 Coverage from Depth/Bedgraph

The bedgraph file is parsed into a table with columns `(chrom, start, end, depth)`. For each reference:

```
covered_bases = sum of (end - start + 1) for every region where depth > 0

coverage = covered_bases / reference_length
```

> **Bedgraph is optional when `--compare_references` is active.** In that mode the pipeline tiles fixed-size windows across the reference FASTAs directly and compares them via k-mer sketches (Section 3.8). Coverage for scoring is then estimated from the BAM index statistics rather than from an external bedgraph file. If no bedgraph is supplied without `--compare_references`, `determine_conflicts()` raises a `ValueError` immediately so the user is not left with a empty results file.

### 2.3 Abundance Metrics (RPM / RPKM)

Computed in `compute_scores_per()`. Raw read counts aren't directly comparable between samples because some samples have more total reads than others, and some organisms have bigger genomes (which naturally attract more reads). RPM and RPKM normalize for these differences:

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

Closely related organisms (e.g., _E. coli_ vs _Shigella_) share large regions of their genomes. This means reads from one organism can legitimately align to the other — a problem called **cross-mapping**. Cross-mapping inflates the apparent confidence for organisms that aren't really there.

This step finds those shared regions and removes a proportional share of ambiguous reads before scoring.

### 3.2 Region Merging

Raw bedgraph intervals are merged into larger contiguous regions using `merge_bedgraph_regions()`. The default strategy (`"jump"`) groups adjacent intervals when:

1. The depth change between consecutive intervals is ≤ a threshold (set automatically from the 97.5th percentile of depth)
2. The gap between intervals is ≤ 10% of the reference length (editable via `--gap_allowance`)
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

`build_conflict_groups()` builds a network (graph) of similar regions:

- **Nodes** are genomic regions `(ref, start, end)` — each node represents one stretch of one organism's genome.
- **Edges** connect regions from different organisms that look too similar (Jaccard similarity ≥ threshold). **Jaccard similarity** is a simple measure of overlap: it counts how many k-mers two regions share divided by the total unique k-mers between them. A value of 1.0 means identical; 0.0 means nothing in common.

Connected components (found via BFS — breadth-first search) define **conflict groups** clusters of regions across multiple organisms that share significant k-mer content and therefore could cause reads to be ambiguously assigned.

### 3.6 Best-Alignment Read Removal

For each read with multi-reference alignments, `build_removed_ids_best_alignment()` keeps the alignment with the highest composite score and schedules the rest for removal. The score combines four terms:

```math
\text{score} = \text{MAPQ}
             + w_{AS} \cdot AS
             - w_{pen} \cdot N_{\text{alt}}
             + w_{ani} \cdot \Omega_{\text{ani}}
```

Where:

- **MAPQ** — raw mapping quality reported by the aligner (0–60)
- $w_{AS}$ — alignment-score weight (default `0.0`; AS tie-breaks only when enabled)
- $N_{\text{alt}}$ — number of competing references that share at least one high-ANI window with this contig; $w_{pen} = 50.0$ by default
- $\Omega_{\text{ani}}$ — **ANI dominance term** (see below); $w_{ani} = 10.0$ by default

#### ANI Dominance Term $\Omega_{\text{ani}}$

Before the BAM is scanned, `_count_reads_per_ref()` builds a fast index of the pre-removal read count for every reference from the BAM index statistics. For each competing reference that shares a high-ANI window with the current contig, we measure how much more (or less) support our contig has:

```math
\Omega_{\text{ani}} = \sum_{\text{alt} \in A} J_{\max}(\text{contig},\, \text{alt}) \cdot \log_2\!\left(\frac{N_{\text{contig}}}{N_{\text{alt}}}\right)
```

Where:

- $A$ — the set of competing contigs that share at least one window with this contig
- $J_{\max}(\text{contig}, \text{alt})$ — the maximum window-level Jaccard similarity observed between the contig and that competitor (from the shared-window index). Jaccard 1.0 = nearly identical sequence; 0.1 = only slightly similar
- $N_{\text{contig}}$, $N_{\text{alt}}$ — pre-removal read counts for the respective references (minimum 1 to avoid log(0))
- $\log_2(N_{\text{contig}} / N_{\text{alt}})$ — the log-ratio of read counts. This is positive when our reference has more reads than the competitor (reward) and negative when it has fewer (penalty)

**Worked example — O104:H4 vs Shigella:** Suppose O104:H4 has 80,000 reads and a Shigella contig has 5,000, and they share a window with $J_{\max} = 0.95$:

```
ani_dominance  = 0.95 × log₂(80000 / 5000)
               = 0.95 × log₂(16)
               = 0.95 × 4.0 = 3.8

score boost (w_ani=10) → +38 points for O104:H4 alignment
score penalty          → −38 points for Shigella alignment
```

A single high-ANI window with a 16× read-count advantage adds ≈ 38 points, decisively overcoming any MAPQ noise between the two references. Reads that align equally well to both will therefore be assigned to O104:H4.

#### Parameters

The actual values passed to `build_removed_ids_best_alignment()` differ by mode:

| Parameter          | `--compare_references` mode | Standard (bedgraph) mode | Effect                                                         |
| ------------------ | --------------------------- | ------------------------ | -------------------------------------------------------------- |
| `penalize_weight`  | 0                           | 23.0                     | Penalty per competing contig in $N_{\text{alt}}$               |
| `as_weight`        | 1.0                         | 0.0                      | Weight on the aligner AS tag                                   |
| `ani_boost_weight` | 10.0                        | _(not used)_             | Scales $\Omega_{\text{ani}}$; raise to make dominance stronger |
| `min_alt_count`    | 2                           | 1                        | Minimum competing alignments before removal applies            |

In `--compare_references` mode the ANI dominance term ($\Omega_{\text{ani}}$) does all the heavy lifting, so `penalize_weight` is set to 0. In the standard bedgraph path, `as_weight` is disabled (0.0) and the competing-contig penalty is 23.0.

### 3.7 Comparison Metrics (Δ Breadth)

After read removal, `compare_metrics()` calculates how much each reference changed. The results are written to the conflict-comparison xlsx/csv. Each row represents one reference and includes:

**Core removal metrics** (always present):

| Column        | Formula                                            | Meaning                                          |
| ------------- | -------------------------------------------------- | ------------------------------------------------ |
| `Δ All`       | $N_{\text{total}} - N_{\text{pass}}$               | Number of reads removed                          |
| `Δ All%`      | $100 \times N_{\text{removed}} / N_{\text{total}}$ | Percentage of reads removed                      |
| `Δ⁻¹ Breadth` | $B_{\text{post}} / B_{\text{pre}}$                 | Fraction of genome breadth that survived removal |

**Abundance metrics** (pre- and post-removal, so downstream scoring can see the delta):

| Column      | Formula                                                                                    |
| ----------- | ------------------------------------------------------------------------------------------ |
| `RPKM Pre`  | $N_{\text{pre}} \;\big/\; \big(\tfrac{R_{\text{total}}}{10^6} \times \tfrac{L}{10^3}\big)$ |
| `RPKM Post` | Same formula using $N_{\text{post}}$                                                       |
| `RPM Pre`   | $N_{\text{pre}} \;\big/\; \big(R_{\text{total}} / 10^6\big)$                               |
| `RPM Post`  | Same formula using $N_{\text{post}}$                                                       |

Where $N_{\text{pre}}$ and $N_{\text{post}}$ are read counts before and after removal, $R_{\text{total}}$ is the total reads in the BAM, and $L$ is the reference length in bp.

**Shared-window ANI summary** (populated when `--compare_references` is active; derived from `compute_shared_window_stats()`):

| Column             | Meaning                                                                                         |
| ------------------ | ----------------------------------------------------------------------------------------------- |
| `Shared Windows`   | Number of tiled windows on this reference that overlap at least one window on another reference |
| `Conflicting Refs` | Count of distinct other references sharing at least one window                                  |
| `Mean Shared ANI`  | Average Jaccard similarity across all shared windows (proxy for ANI; 0→1)                       |
| `Max Shared ANI`   | Highest single-window Jaccard seen for this reference — flags the most similar competitor       |
| `Shared BP`        | Approximate total base-pairs covered by shared windows (windows may overlap)                    |

`Δ⁻¹ Breadth` and `Δ All%` are passed directly into the minhash score (see Section 7). The ANI summary columns are informational in the xlsx but their per-window Jaccard values feed the ANI dominance term $\Omega_{\text{ani}}$ during read removal (Section 3.6).

### 3.8 Shared-Window Mode (`--compare_references`)

When `--compare_references` is passed to `match_paths.py`, an alternative — and more thorough — conflict detection approach is used. Rather than sketching merged bedgraph regions, the pipeline tiles fixed-size windows across every reference FASTA and compares them directly to find shared sequence.

#### How it works

1. **Tiling** — Each reference FASTA is divided into windows of `--window_size` bp (default: 900,000 bp) with a step of `--step_size` bp (default: 900,000 bp; set to `window_size/2` for 50% overlap). Windows with more than 5% N-bases are skipped.
2. **Sketching** — Each window is converted to a MinHash signature (k-mer size 31, scaled 200 by default).
3. **Pairwise search** — A `ProcessPoolExecutor` is used to compare all query windows against all target windows from other FASTAs. Pairs with Jaccard ≥ `jaccard_threshold` (default 0.10) are recorded as shared.

#### Performance design

`report_shared_windows_across_fastas()` is optimised to minimise inter-process communication overhead at scale:

- **Worker initializer** — Target window hash-sets are sent to each worker process _once_ via the pool initializer (`_init_search_worker`). They are deserialized and stored in a module-level dict (`_SEARCH_WORKER_TARGETS`), keyed by FASTA label. Subsequent job submissions carry only the query chunk — no per-job target pickling.
- **Frozenset Jaccard** — Similarity is computed with Python's native `frozenset &` intersection (C-level), not via `MinHash` objects. This avoids Python-object construction overhead on every comparison: $J = |A \cap B| / |A \cup B|$ where $A$ and $B$ are `frozenset`s of hash values. An early-exit is applied when $|A \cap B| = 0$.
- **Intra-FASTA skip** — Same-FASTA comparisons are skipped in $O(1)$ by checking the FASTA label key before entering the inner loop.
- **Job granularity** — One job per query chunk (chunk size = ceil(N_queries / n_workers)), compared to the earlier approach of N_tiles × chunk jobs. This reduces scheduler overhead for large databases.

#### Parameters

| CLI flag               | Default | Effect                                                        |
| ---------------------- | ------- | ------------------------------------------------------------- |
| `--window_size`        | 900,000 | Window width in bp                                            |
| `--step_size`          | 900,000 | Step between window starts; `window_size/2` gives 50% overlap |
| `--min_threshold`      | 0.4     | Jaccard threshold for the region-based (non-window) path      |
| `--compare_references` | off     | Enable this mode                                              |

> **Note:** The bedgraph (`-b`) is **not required** when `--compare_references` is active. The pipeline infers coverage from BAM index statistics instead. The FASTA files (`-f`) **are** required because windows are built from reference sequence.

### 3.9 ANI Matrix Generation

`generate_ani_matrix()` builds an organism-level Average Nucleotide Identity (ANI) matrix:

1. Composite signatures are built per taxid from all its accessions
2. Pairwise ANI is estimated from MinHash containment between organism signatures
3. The result is written to `organism_ani_matrix.csv` as a symmetric matrix

---

### 3.10 Pre-Removal Read Count Index (`_count_reads_per_ref`)

Before any reads are removed, `_count_reads_per_ref()` builds a dict of `{reference: read_count}` for every contig in the BAM. It uses the BAM index statistics (`pysam.AlignmentFile.get_index_statistics()`) which is an O(N_references) operation — no reads are iterated. If the index is unavailable, it falls back to a full single-pass scan.

This index is what enables the $N_{\text{contig}} / N_{\text{alt}}$ ratio inside the ANI dominance term $\Omega_{\text{ani}}$ (Section 3.6).

### 3.11 Shared-Window Conflict Summary (`compute_shared_window_stats`)

After the shared-window search, `compute_shared_window_stats()` aggregates the per-window results into a per-reference summary that is stored in the comparison xlsx (Section 3.7). For each reference it computes:

| Field                 | Computation                                                                      |
| --------------------- | -------------------------------------------------------------------------------- |
| `n_shared_windows`    | `len(windows)` — count of tiled windows on this ref that have at least one match |
| `n_conflicting_refs`  | `len({w.alt_contig for w in windows})` — distinct competitors                    |
| `mean_window_jaccard` | `statistics.mean([w.jaccard for w in windows])`                                  |
| `max_window_jaccard`  | `max(w.jaccard for w in windows)`                                                |
| `shared_bp`           | `sum(w.end - w.start for w in windows)` (windows may overlap; approximate)       |

`max_window_jaccard` is the value used as $J_{\max}$ in the ANI dominance formula, so it directly affects which organism "wins" when reads are ambiguous (Section 3.6).

## 4. Stage 3: Score Computation (`optimize_weights.py`)

All component scores are computed in `compute_scores_per()` and combined in `compute_tass_score()`.

---

## 5. Component 1: Breadth Log Score (Genome Coverage)

### 5.1 What It Measures

The breadth score captures how much of the reference genome has at least one read aligned to it. If reads are spread across the genome, that's a good sign the organism is genuinely present. If reads pile up in just one spot, it's more likely cross-mapping or contamination.

This fraction is passed through a sigmoid curve to produce a 0–1 confidence value.

### 5.2 Sigmoid Function

A **sigmoid** is an S-shaped curve that smoothly converts any input value into a number between 0 and 1. Think of it like a dimmer switch: instead of flipping abruptly from "off" to "on," it gradually ramps up. We use it here because we don't want the score to jump from 0 to 1 the instant coverage crosses a threshold — we want a smooth transition.

Defined in `breadth_score_sigmoid()`:

```math
\text{breadth\_sigmoid}(c) = \frac{1}{1 + e^{-s \cdot (c - m)}}
```

Where:

- $c$ = the coverage fraction (covered_bases / reference_length). This is just "what percentage of the genome has at least one read on it?" It ranges from 0 (nothing covered) to 1 (everything covered).
- $m$ = the midpoint — the coverage level where the sigmoid outputs 0.5, meaning "we're 50/50 on whether this is real." Default is 0.01 (i.e., 1% of the genome covered).
- $s$ = the steepness — controls how quickly the curve transitions from low confidence to high confidence around the midpoint. A higher number makes the jump sharper. Default is 12,000.
- $e$ = Euler's number (~2.718), just a mathematical constant that gives the curve its smooth S-shape.

**Overflow protection:** If $s \cdot (c - m) \geq 50$, the function returns 1.0. If $\leq -50$, it returns 0.0.

### 5.3 MAPQ Scaling

Low mapping quality (MAPQ ≈ 0) suggests reads are aligning ambiguously — they could plausibly belong to multiple organisms, so we can't be sure they belong to this one. When most reads have low MAPQ, we scale the breadth score down to reflect that uncertainty:

```math
\text{mapq\_scale} = (\text{highmapq\_fraction})^{p}
```

Where $p$ (the exponent) = `mapq_breadth_power` (default: 2.0). Squaring the high-MAPQ fraction means that if, say, only 70% of reads are high-quality, the penalty is 0.7² = 0.49 — a meaningful reduction. This makes the score more conservative when read quality is mixed.

### 5.4 Final Breadth Log Score

```math
\boxed{\text{breadth\_log\_score} = \text{breadth\_sigmoid}(c, m, s) \times (\text{highmapq\_fraction})^{p}}
```

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

Inspiration for Gini is derived from a paper that used gini to determine gene expression information from [GeneGini](https://doi.org/10.1016/j.cels.2018.01.003) (DOI: 10.1016/j.cels.2018.01.003).

### 6.1 What It Measures

The Gini score measures how evenly reads are spread across the genome. If an organism is truly present, reads should appear across many different genomic positions — not just cluster in one conserved region.

The pipeline inverts the classical Gini coefficient so that **uniform coverage → high score** and **clumped coverage → low score**.

### 6.2 Step-by-Step Computation

The full computation lives in `getGiniCoeff()` and has five sub-steps.

#### Step 1: Build a Transformed Coverage Histogram

Each region's depth is log-transformed to reduce the effect of extreme outliers. Imagine one region has a depth of 10,000 while most are around 10 — without the log transform, that one outlier would dominate the calculation. The log squashes those extremes:

```math
d'_i = \log_{10}(1 + d_i)
```

Here, $d_i$ is the raw depth at position $i$, and $d'_i$ is the transformed (compressed) version. The "+1" inside the log prevents taking log of zero (which is undefined). Uncovered bases get $\log_{10}(1 + 0) = 0$.

A histogram is then built that maps each transformed depth value to the number of bases at that depth.

#### Step 2: Compute Raw Gini from the Lorenz Curve

The Gini coefficient is borrowed from economics (where it measures wealth inequality) and repurposed here to measure how unevenly reads are distributed across the genome. The idea: if coverage is perfectly even everywhere, Gini is 0. If all coverage is piled into one tiny spot, Gini approaches 1.

To calculate it:

1. Sort coverage values from lowest to highest
2. Build a **Lorenz curve** — a line that shows "what fraction of total coverage is accounted for by the bottom X% of genome positions?" In a perfectly equal world, the bottom 50% of positions would hold 50% of coverage. In practice, low-coverage positions hold very little.
3. Measure the area under the Lorenz curve ($A_L$) using a standard trapezoid method:

```math
A_L = \sum_{i} \frac{(x_{i+1} - x_i)(y_i + y_{i+1})}{2}
```

Here $(x_i, y_i)$ are the Lorenz curve points, where $x_i$ is the cumulative fraction of bases and $y_i$ is the cumulative fraction of total coverage at that point.

```math
G_{\text{raw}} = 1 - 2 A_L
```

$G_{\text{raw}}$ falls in the range [0, 1]: 0 means perfectly uniform coverage, 1 means maximally unequal (all reads in one spot).

#### Step 3: Invert and Scale the Raw Gini

Since a low raw Gini means uniform coverage (which is a _good_ sign), we flip and scale it so that uniform coverage gives a _high_ score:

```math
G_{\text{transformed}} = \text{clamp}\!\Big(\alpha \cdot \sqrt{1 - G_{\text{raw}}},\; 0,\; 1\Big)
```

Where $\alpha$ (alpha) is a scaling multiplier set to 1.8 by default. Setting it above 1.0 means that even moderately uniform distributions can reach the maximum score of 1.0. The square root softens the transformation so that small improvements in uniformity still get meaningful credit. The `clamp(…, 0, 1)` at the end ensures the result never goes below 0 or above 1. Examples:

- $G_{\text{raw}} = 0$ (perfectly uniform) → $G_{\text{transformed}} = \min(1.8, 1.0) = 1.0$
- $G_{\text{raw}} = 0.9$ (very clumped) → $G_{\text{transformed}} = 1.8 \times \sqrt{0.1} \approx 0.57$

#### Step 4: Length-Based Scaling Factor

Larger genomes need more reads to achieve uniform coverage, so it's harder to get a good Gini score for a big genome. To compensate, we give a bonus that scales with genome size. We use a log scale so the bonus grows gradually rather than exploding for very large genomes:

```math
S_{\text{length}} = 1 + R \cdot \log_{10}\!\Big(\max\!\big(\frac{\min(L, L_{\max})}{L_{\text{base}}},\; 1\big)\Big)
```

Where:

- $L$ = the genome length (in base pairs)
- $L_{\text{base}}$ = a baseline genome length (default: 500,000 bp). Genomes shorter than this get no bonus.
- $L_{\max}$ = a length cap (default: 10⁹ bp) to prevent absurdly large genomes from getting too much of a boost.
- $R$ = the reward factor (default: 2), which controls how generous the bonus is. Higher R = more reward for larger genomes.

#### Step 5: Positional Dispersion Factor

This measures how spread out the covered regions are physically along the genome. Even if the depth of coverage is uniform, we want reads scattered across many different locations and not all clumped in one contiguous block. Think of it like this: if you're checking whether someone actually read a whole book, you'd be more convinced if they could talk about many passages from every chapter with general detail instead of having a highly detailed summary of just one chapter.

```math
\bar{m} = \frac{1}{n}\sum_{i=1}^{n} \frac{s_i + e_i}{2}
\qquad
\sigma^2 = \frac{1}{n}\sum_{i=1}^{n}(m_i - \bar{m})^2
\qquad
D = \sqrt{\frac{\sigma^2}{L^2 / 12}}
```

Breaking this down:

- First, we find the **midpoint** of each covered region: $(s_i + e_i) / 2$, where $s_i$ is the start position and $e_i$ is the end position of region $i$.
- $\bar{m}$ (m-bar) is the **average midpoint** across all $n$ regions. It tells us where the "center of mass" of coverage is.
- $\sigma^2$ (sigma squared) is the **variance** of those midpoints — how spread out they are from the average. If all regions are near the same spot, variance is low. If they're scattered across the genome, variance is high.
- **$L^2/12$** is the key normalization term. It represents the **maximum possible variance** you'd see if regions were spread out as evenly as possible (uniformly distributed) across a genome of length $L$. This comes from a well-known statistical fact: a uniform distribution over any range [0, L] has a variance of exactly L²/12. By dividing the actual variance by this, we get a ratio that tells us "how spread out are our regions compared to the best-case scenario?" The denominator has been shown with some slight benchmarking to be sufficient for most genomes but could use some additional tuning in the future.
- $D$ (the Dispersion factor) is the square root of that ratio, which brings it back to the same scale as the original positions. $D$ falls in the range [0, 1]: a value near 0 means everything is clustered together, while a value near 1 means the regions are spread out almost as well as theoretically possible.

#### Step 6: Final Gini Score

```math
\boxed{\text{gini\_coefficient} = \min\!\Big(1.0,\;\; G_{\text{transformed}} \times S_{\text{length}} \times (1 + \beta \cdot D)\Big)}
```

Where $\beta$ (beta) is the **dispersion weight** (default: 0.5). It controls how much the positional spread of reads (the $D$ factor from Step 5) influences the final Gini score. At 0.5, if regions are perfectly spread out ($D = 1.0$), the score gets a 50% boost (multiplied by 1.5). If $\beta$ were 0, dispersion would be ignored entirely. If $\beta$ were 1.0, perfect dispersion would double the score. The default of 0.5 strikes a balance: spatial spread matters, but it's not the dominant factor.

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

- **`Δ⁻¹ Breadth`** ($B_r$): The ratio of post-removal breadth to pre-removal breadth. Think of it as "what fraction of genome coverage survived after we removed ambiguous reads?" A value of 1.0 means no breadth was lost (good — the signal was all unique). A value of 0.5 means half the coverage came from reads that could've belonged to a different organism (concerning).
- **`Δ All%`** ($\Delta\%$): The percentage of total reads that were removed during conflict resolution. The Δ (delta) here means "change" — specifically how many reads changed (got removed).

The raw minhash score blends two retention signals — breadth retention and read retention — to differentiate true positives from false positives, particularly within high-ANI conflict groups (e.g. _E. coli_ vs _Shigella_) where $B_r$ alone is nearly identical for every member:

- **`Δ⁻¹ Breadth`** ($B_r$): fraction of genome breadth surviving conflict removal (0–1).
- **Read retention** ($R_r$): fraction of reads surviving conflict removal, i.e. `Pass Filtered Reads / Total Reads`.

```math
\text{minhash\_score} = \min\!\big(1.0,\;\; 0.30 \cdot B_r + 0.70 \cdot R_r\big)
```

Read retention is weighted more heavily (70%) because it differentiates conflict-group members much more cleanly: the dominant true-positive organism retains most of its reads during removal, while false-positive relatives lose the majority of theirs.

**Example — O104:H4 vs Shigella in the same conflict group:**

```
O104:H4  (TP):  Br=0.12, read_retention=0.85  →  0.30×0.12 + 0.70×0.85 = 0.631
Shigella (FP):  Br=0.12, read_retention=0.30  →  0.30×0.12 + 0.70×0.30 = 0.246
```

Both organisms have the same $B_r$ (they share similar genome breadth reduction in the high-ANI group), but their read-retention scores clearly separate them.

**Fallback (no comparison data available):**

When conflict comparison data isn't available, we fall back to a simpler estimate based on how many reads the organism has relative to the total:

```math
\text{minhash\_score}_{\text{fallback}} = \text{rpm\_confidence} \times 0.5
```

Where `rpm_confidence` is a sigmoid over the fraction of total reads ($f_{\text{reads}}$). The very steep sigmoid (steepness = 50,000) acts almost like a switch: if the organism has more than a tiny fraction of reads (above 0.01%), confidence jumps to ~1.0; below that, it drops to ~0:

```math
\text{rpm\_confidence} = \frac{1}{1 + e^{-50000 \cdot (f_{\text{reads}} - 0.0001)}}
```

### 7.3 Confidence Gating

Before the minhash score is used, it's "gated" by a coverage confidence factor. Think of a gate as a checkpoint: if an organism doesn't have enough genome coverage, we don't trust its minhash score regardless of how good it looks. This prevents organisms with very sparse coverage from getting artificially inflated scores:

```math
\text{conf} = w_b \cdot \text{breadth\_sigmoid}(c) + w_g \cdot G_{\text{score}}
```

Here $w_b$ and $w_g$ are weights that control how much the breadth sigmoid vs. the Gini score contribute to the gate. By default $w_b = 1.0$ and $w_g = 0.0$, so gating is based purely on breadth (genome coverage).

```math
\boxed{\text{minhash\_reduction} = \text{clamp}(\text{conf},\; 0,\; 1)}
```

> **Note:** The value stored as `minhash_reduction` is the _confidence_ value, not the raw minhash score itself. This means the minhash component represents how confident we are that the organism's alignment signal is real — primarily driven by whether the organism has meaningful genome breadth.

### 7.4 Low-Read Penalty with Coverage Bypass

Organisms with very few reads relative to the total are more likely to be false positives, so their `minhash_confidence` is scaled down by a coverage-aware penalty. Rather than applying the penalty at full strength regardless of coverage evidence, it is **faded out when the organism has solid genome breadth** — preventing a well-covered dominant organism from being dragged down by the global read-fraction sigmoid:

```math
w_{\text{pen}} = 0.35 \times \bigl(1 - \min\!\bigl(1.0,\; \text{cov\_conf} + 0.20 \times \text{mapq\_norm}\bigr)\bigr)
```

where $\text{mapq\_norm} = \min(1.0, \text{meanmapq} / 60)$. The additional MapQ bypass term (`0.20 × mapq_norm`) fades the penalty further when reads are high-quality — a high-MAPQ organism with few reads (e.g. a rare pathogen in a sterile-site sample) should not be penalised as harshly as a low-MAPQ scatter hit.

```math
\text{minhash\_confidence} \mathrel{*}= \big(1 - w_{\text{pen}}\big) + \big(w_{\text{pen}} \times \text{rpm\_weight}\big)
```

Where `rpm_weight` is a steep sigmoid over the fraction of total reads:

```math
\text{rpm\_weight} = \frac{1}{1 + e^{-50000 \cdot (f_{\text{reads}} - 0.0001)}}
```

and `cov_conf = breadth_sigmoid(coverage)`.

**Effect at different coverage levels (with mapq_norm = 0, i.e. low MAPQ):**

| `cov_conf`          | $w_{\text{pen}}$ | Penalty at `rpm_weight=0.1` | Meaning                              |
| ------------------- | ---------------- | --------------------------- | ------------------------------------ |
| 0.95 (100% breadth) | 0.017            | ×0.998 (≈ none)             | Dominant organism — penalty bypassed |
| 0.50 (moderate)     | 0.175            | ×0.843 (−16%)               | Moderate dampening                   |
| 0.10 (sparse)       | 0.315            | ×0.717 (−28%)               | Full penalty applies                 |

When `mapq_norm = 1.0` (MAPQ 60), the effective coverage bypass is raised by 0.20, reducing $w_{\text{pen}}$ further at every coverage level.

This ensures that a high-confidence organism like O104:H4 — which spans the full genome — is never unfairly penalised simply because its read fraction is computed against total sample reads that include many other organisms.

### 7.5 RPM-Normalised Minhash Reduction (Multi-Member Groups)

Within each parent group (organisms sharing the same `toplevelkey` — i.e. the same species or taxid group), minhash scores are further adjusted based on each organism's share of total RPM within that group. This is the final normalisation step and runs _after_ the confidence gate and read-penalty.

#### Solo organisms (only member in their group)

A single organism in its group gets an **exclusivity boost**: the gap between its current `minhash_reduction` and the coverage-gate ceiling is partially filled, proportional to `rpm_confidence_weight`:

```math
\text{gap} = \max(0,\; \text{minhash\_confidence} - \text{raw\_minhash})
\qquad
\text{boost} = \text{gap} \times 0.9 \times \text{rpm\_weight}
\qquad
\text{adjusted} = \min(\text{minhash\_confidence},\; \text{raw\_minhash} + \text{boost})
```

The `minhash_confidence` ceiling ensures a low-coverage solo organism can't inflate itself beyond what its coverage evidence supports.

#### Multi-member groups — dominant organism

The member with the highest RPM in the group is **boosted toward its coverage-gate ceiling** rather than scaled down. This is the key change for high-ANI groups like O104:H4 vs Shigella:

```math
\text{gap} = \max(0,\; \text{minhash\_confidence} - \text{raw\_minhash})
\qquad
\text{boost} = \text{gap} \times \text{rpm\_share}
\qquad
\text{adjusted} = \min(\text{minhash\_confidence},\; \text{raw\_minhash} + \text{boost})
```

Where $\text{rpm\_share} = \text{rpm}_{\text{this}} / \sum_{\text{group}} \text{rpm}$.

**Example:** O104:H4 has 60% of group RPM, $\text{raw\_minhash} = 0.90$, $\text{minhash\_confidence} = 0.95$:

```
gap     = 0.95 − 0.90 = 0.05
boost   = 0.05 × 0.60 = 0.03
adjusted = min(0.95, 0.90 + 0.03) = 0.93   ← was 0.72 before this change
```

#### Multi-member groups — losing organisms

All non-dominant members are penalised using a two-part multiplicative scale:

**Part 1 — RPM-share floor:**

```math
\text{base\_scale} = \underbrace{0.10}_{\text{floor}} + 0.90 \times \text{rpm\_share}
```

The floor of 0.10 (reduced from the previous 0.30) prevents partial false positives from retaining an inflated minhash score. With a 10% floor, Shigella at 8% RPM share gets a base scale of ≈ 0.17.

**Part 2 — Coverage-ratio penalty:**

```math
r_{\text{cov}} = \min\!\Big(1.0,\; \frac{c_{\text{this}}}{c_{\text{dominant}}}\Big)
\qquad
\text{cov\_scale} = 0.25 + 0.75 \times r_{\text{cov}}
\qquad
\text{scale} = \text{base\_scale} \times \text{cov\_scale}
```

If the dominant organism has 20% genome coverage and this organism only has 4%, the coverage ratio is 0.20, giving `cov_scale = 0.25 + 0.75×0.20 = 0.40`. This extra penalty captures the case where a low-coverage organism has a similar RPM share simply because its reads all piled into one conserved region.

**Full worked example — Shigella vs O104:H4:**

```
rpm_share  = 0.21   cov = 8%   top_cov = 20%

base_scale = 0.10 + 0.90 × 0.21 = 0.289
r_cov      = 8 / 20 = 0.40
cov_scale  = 0.25 + 0.75 × 0.40 = 0.55
scale      = 0.289 × 0.55 = 0.159

adjusted   = 0.90 × 0.159 = 0.143   ← was 0.72× before
```

The combined effect is that in a high-ANI conflict group, the organism with the most reads and broadest coverage is strongly promoted, while relatives (Shigella, other E. coli strains) have their minhash contribution cut to a small fraction of its naive value.

---

## 8. Additional Scoring Components

### 8.1 HMP Percentile (Body-Site Abundance Context)

This compares the observed abundance to what the Human Microbiome Project found for this organism at the relevant body site. If an organism is at typical abundance for that site, it gets a neutral score. Unusually high or low abundance is flagged.

```math
f_{\text{obs}} = \frac{\text{numreads}}{\text{total\_reads}}
\qquad
z = \frac{f_{\text{obs}} - \mu_{\text{HMP}}}{\sigma_{\text{HMP}}}
\qquad
P_{\text{base}} = \Phi(z)
```

Breaking this down:

- $f_{\text{obs}}$ is the **observed fraction** of reads belonging to this organism — simply its read count divided by total reads.
- $\mu_{\text{HMP}}$ (mu) is the **average abundance** for this organism at this body site, according to the Human Microbiome Project data.
- $\sigma_{\text{HMP}}$ (sigma) is the **standard deviation** of that abundance — how much variation is normal.
- $z$ is the **z-score**: how many standard deviations away from normal our observation is. A z-score of 0 means "exactly average," +2 means "well above average," and -2 means "well below average."
- $\Phi(z)$ (Phi of z) converts the z-score into a **percentile** on a bell curve. For example, a z-score of 0 gives the 50th percentile; a z-score of +2 gives the 97.5th percentile.

The base percentile is then adjusted based on how far from normal it is:

| Z-score range          | Adjusted percentile               | What it means                     |
| ---------------------- | --------------------------------- | --------------------------------- |
| $z < -1.0$             | $P_{\text{base}}^2$               | Penalize — below typical          |
| $-1.0 \leq z \leq 1.0$ | $P_{\text{base}}$                 | Neutral — near expected           |
| $1.0 < z \leq 2.0$     | $1 - (1 - P_{\text{base}})^{1.5}$ | Boost — above typical             |
| $z > 2.0$              | $1 - (1 - P_{\text{base}})^{2}$   | Strong boost — well above typical |

### 8.2 Plasmid Score

Within each species, strains are compared by their plasmid coverage. Calc. in `match_paths.py`.

**Absolute quality** (does this plasmid have good coverage?):

```math
Q = \sqrt{B_{\text{plasmid}} \times G_{\text{plasmid}}}
```

Where $B_{\text{plasmid}}$ is the breadth sigmoid score for the plasmid's coverage and $G_{\text{plasmid}}$ is the Gini coefficient for the plasmid's depth distribution. The square root of this product is called a **geometric mean** — it's a way of averaging two values that requires _both_ to be decent. The reason we use the geometric mean is because it accounts for means in extreme value scenarios i.e. there is a more likely scenario that a small plasmid will either have no coverage or super good coverage. If either breadth or uniformity is near zero, the geometric mean is poor. You can't fake a good plasmid score with good breadth of coverage but terrible uniformity, or vice versa.

**Relative disparity** (how does this strain compare to others in the same species?):

```math
D_{\text{rel}} = \begin{cases}
\min\!\big(1.0,\; 0.7 \cdot \frac{c_{\text{strain}}}{c_{\max}} + 0.3 \cdot \frac{r_{\text{strain}}}{r_{\max}}\big) & \text{if multiple strains have plasmids} \\
1.0 & \text{if this is the only strain with a plasmid}
\end{cases}
```

In other terms: when there are multiple strains of the same species (matched from their shared species level taxid), we compare each strain's plasmid coverage ($c_{\text{strain}}$) and read count ($r_{\text{strain}}$) against the best strain in the group ($c_{\max}$ and $r_{\max}$). Coverage gets 70% of the weight and read count gets 30%. If this is the only strain with a plasmid, there's nothing to compare against, so it gets a 1.0.

**Final plasmid score:**

```math
\text{plasmid\_score} = \min(1.0,\;\; Q \times D_{\text{rel}})
```

Strains with no plasmid accessions receive `plasmid_score = 0`. This is not a penalty so that we don't strains/species that lack plasmids in the database/references.

### 8.3 Low-Abundance Confidence

This is a sigmoid curve applied in log₁₀-RPM space that gives credit to organisms that have a meaningful number of reads per million, even when the absolute read count is small. It's especially useful for sterile-site or blood samples where you don't expect many reads from any organism, but even a few reads can be clinically significant.

```math
\text{RPM} = \frac{\text{numreads}}{\text{total\_reads}} \times 10^6
```

```math
\text{abundance\_confidence} = \frac{1}{1 + e^{-s \cdot (\log_{10}(\text{RPM}) - \log_{10}(m))}}
```

Where $s$ = steepness (default: 2.0) and $m$ = RPM midpoint (default: 5.0).

### 8.4 K2 Disparity Score

Compares the pipeline's species-level classification against Kraken2's (or Centrifuge if enabled instead) classification. If Kraken2/Centrifuge assigns reads to a different genus, this reduces confidence. **Disabled by default** (`k2_disparity_score_weight = 0`).

### 8.5 DIAMOND Identity

Average amino acid identity from a DIAMOND BLASTX protein alignment, weighted by CDS count. Provides protein-level confirmation of the taxonomic assignment. **Disabled by default** (`diamond_identity = 0`). We have not tested this since v1.0 so consider it less reliable of a tool to implement. The associated modules are still functional, however.

---

## 9. Final TASS Score Formula

### 9.1 The Core Weighted Sum

Computed in `compute_tass_score()`. The final score is simply a **weighted sum** — each component score ($x_i$) is multiplied by its relative weight ($w_i$), and then they're all added together and clamped between 0 and 1. Components with higher weights have more influence on the final result:

```math
\boxed{
\text{TASS}_{\text{base}} = \sum_{i} w_i \cdot x_i
}
```

Expanded:

```math
\text{TASS}_{\text{base}} = w_b \cdot \text{breadth\_log\_score}
\;+\; w_m \cdot \text{minhash\_reduction}
\;+\; w_g \cdot \text{gini\_coefficient}
\;+\; w_h \cdot \text{hmp\_percentile}
```

```math
\;+\; w_d \cdot \text{disparity}
\;+\; w_q \cdot \text{mapq\_score}
\;+\; w_{k2} \cdot \text{k2\_disparity\_score}
\;+\; w_{di} \cdot \text{diamond\_identity}
\;+\; w_{ac} \cdot \text{abundance\_confidence}
```

### 9.2 Default Weights

| Component              | Weight   | CLI Flag                      | Type     |
| ---------------------- | -------- | ----------------------------- | -------- |
| `breadth_log_score`    | **0.27** | `--breadth_weight`            | included |
| `minhash_reduction`    | **0.31** | `--minhash_weight`            | included |
| `gini_coefficient`     | **0.42** | `--gini_weight`               | included |
| `hmp_percentile`       | 0.00     | `--hmp_weight`                | included |
| `disparity`            | 0.00     | `--disparity_weight`          | included |
| `mapq_score`           | 0.00     | `--mapq_score`                | included |
| `k2_disparity_score`   | 0.00     | `--k2_disparity_score_weight` | included |
| `diamond_identity`     | 0.00     | `--diamond_identity`          | included |
| `plasmid_bonus_weight` | 0.00     | `--plasmid_bonus_weight`      | bonus    |

**What if all weights don't add up to 1.0?** The three active primary weights (breadth=0.27, minhash=0.31, gini=0.42) sum to 1.0, which is the intended normalization. The plasmid bonus and abundance confidence weights sit outside this normalized pool and are additive on top of the base sum. This is intentional — the final score is clamped to $[0, 1]$ at the end, so overshooting is fine. It allows components to reinforce each other when the support for an organism is strong.

### 9.3 Additive Modifiers

**Plasmid bonus** (applied after the base sum):

```math
\text{TASS}_1 = \text{TASS}_{\text{base}} + w_{\text{plasmid}} \times \text{plasmid\_score}
```

Where $w_{\text{plasmid}}$ = `plasmid_bonus_weight` (default: 0.0 — disabled).

This can also work with additional additive modifiers to be added in the future.

### 9.4 Multiplicative Abundance Gate (Optional)

When `--abundance_gate` (not recommended for most samples) is enabled, the full score is multiplied by the abundance confidence sigmoid. This collapses scores for organisms with very low read fractions:

```math
\text{TASS}_2 = \text{TASS}_1 \times \text{abundance\_confidence}
```

**Note:** This is exclusive with the additive `plasmid_bonus_weight` approach.

### 9.5 Power Transform (Optional)

When `score_power` isn't set to the default of 1.0, a power transform reshapes the score distribution. This is useful for adjusting how "generous" or "strict" the final scores are. A power less than 1 rises low scores (making the system more lax), while a power greater than 1 would push low scores even lower (more conservative):

```math
\text{TASS}_3 = (\text{clamp}(\text{TASS}_2, 0, 1))^{p}
```

Where $p$ = `score_power`. Some examples:

- $p = 0.5$: score 0.09 → 0.30, score 0.95 → 0.97 (rises low scores, making borderline organisms more visible)
- $p = 0.3$: score 0.09 → 0.52, score 0.95 → 0.98 (dramatically rises low scores)
- $p = 1.0$: no change (default — the score passes through as-is)

### 9.6 Final Clamping

```math
\boxed{\text{TASS}_{\text{final}} = \text{clamp}(\text{TASS}_3,\; 0,\; 1)}
```

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

> **Note:** The bedgraph is **not required** when using `--compare_references`. In that mode the pipeline tiles windows across the FASTA files directly and derives coverage from the BAM index. You can omit `-b` entirely from the command below if you are using `--compare_references`.

4. The match file should be a TSV with columns: `accession`, `taxid`, `organism_name`, `description`. This can be generated from the reference FASTA headers and a taxonomy lookup. Example:

| Acc         | Mapped_Value |
| ----------- | ------------ |
| NC_000913.3 | 511145       |
| NC_000964.3 | 224308       |
| NC_001781.1 | 11250        |
| NC_001803.1 | 208893       |

You can also provide custom names with the `--namecol` option if you don't want the taxdump names to be used by default. You'll need to have the column index (0-index) added to the match file and specify the column name like `--namecol 2`. Example:

| Acc         | Mapped_Value | Custom_Name                                |
| ----------- | ------------ | ------------------------------------------ |
| NC_000913.3 | 511145       | Escherichia coli str. K-12 substr. MG1655  |
| NC_000964.3 | 224308       | Bacillus subtilis subsp. subtilis str. 168 |
| NC_001781.1 | 11250        | Human respiratory syncytial virus          |
| NC_001803.1 | 208893       | Pseudomonas aeruginosa PAO1                |

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

:warning: The above command is an example and may require adjustments based on your specific file paths, sample types, and desired parameters. Not all of these args are required as well, so you can omit those that don't apply to your use case. For instance, if you don't have Kraken2/Centrifuge results for comparison, you can leave out the `--k2` argument.

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
4. **Genus level** (toplevelkey) — Highest aggregation, used for HMP lookups and final reporting. You can adjust this for higher level lookups such as order or family.

At each level, metrics are either summed (`numreads`, `covered_bases`), averaged (MAPQ), or re-computed from scratch (`gini`, `breadth` from merged regions). The toplevelkey can also be adjusted with the `--rank` parameter in TaxTriage.

---

## 12. Weight Optimization (`ground_truth.py`)

When ground-truth labels are available (e.g., simulated data where read names encode their true origin), the pipeline can automatically optimize scoring weights using:

1. [**Differential Evolution**](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html) — broad search across the full weight space
2. [**SLSQP**](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html) — local refinement starting from the best solution found above

The objective is to maximize the score gap between true positives and false positives, while keeping weights within valid ranges. Optimized weights are written back into the pipeline arguments for the final scoring pass.

See info in the [insilico testing docs](insilico_simulation.md) for how reads were generated to establish a ground truth of hits post-alignment. When optimizing, the reads are mapped to ground truth by reading the fastq read names, which are formatted like this:

`@CP000253.1_621_4/1`

Where the first part is the reference accession is the organism name. From these headers, we can determine which reads belong to which organisms and use that as a basis for evaluating how well different weight combinations separate true positives from false positives. Keep it in mind when establishing your own benchmarks. See the arguments within `match_paths.py` for how to enable weight optimization and adjust the optimization parameters. You can start with just applying `--optimize` provided your BAM file reads have names formatted like the above example. All of the read ids should have a match the the reference FASTA file. Also make sure the accession matches the references in the `-m` argument (match file) so that the taxid lookup works correctly. For example:

| Acc         | Mapped_Value |
| ----------- | ------------ |
| NC_000913.3 | 511145       |
| CP000253.1  | 93601        |
| ......      | ......       |

You can see that CP000253.1 matches the read id when parsing the read info out i.e. `_621_4/1`. We recommend using InsilicoSeq to generate simulate reads.

Lastly, outside of the scope of this document, we've collected various clinical datasets with spike-in positives and have been used for benchmarking to create the "optimized" weights for common clinical sample types. This is found [here](../assets/sampletype_best_thresholds.json) in the repo. The specified sample type is matched based on the platform $[ILLUMINA, NANOPORE]$ and body site $[nasal, oral, skin, gut, blood, sterile]$. If you have a dataset with known positives that you'd like to use for optimizing weights, you can add it to this JSON file and re-run the pipeline with `--thresholds_json` pointing to your updated file.

---

## 13. Summary of Key Equations

### Breadth Log Score

```math
\text{breadth\_log\_score} = \frac{1}{1 + e^{-s(c-m)}} \times (\text{highmapq\_fraction})^p
```

### Gini Coefficient Score

```math
\text{gini\_coefficient} = \min\!\Big(1,\; \alpha\sqrt{1-G} \;\times\; \big(1 + R\log_{10}\tfrac{L}{L_b}\big) \;\times\; (1+\beta D)\Big)
```

### Minhash Score (raw)

```math
\text{minhash\_score} = \min\!\big(1.0,\;\; 0.30 \cdot B_r + 0.70 \cdot R_r\big)
```

### Minhash Reduction (confidence-gated)

```math
\text{minhash\_reduction} = \text{clamp}\!\big(w_b \cdot \text{breadth\_sigmoid}(c) + w_g \cdot G_{\text{score}},\; 0,\; 1\big)
```

further adjusted by the low-read penalty (Section 7.4) using a MAPQ-aware coverage bypass.

### TASS Score

```math
\text{TASS} = \text{clamp}\!\bigg(\Big[\sum_i w_i x_i + w_{\text{plasm}} \cdot P\Big] \times \text{gate}^{[ab\_gate]}\bigg)^{p}
```

In words: add up all the weighted component scores, toss in the plasmid bonus, optionally multiply by the abundance gate (if enabled; otherwise it's just ×1.0 which does nothing), raise to the power $p$ (which is 1.0 by default, meaning no change), and clamp the result to [0, 1].

#### Notations

Common math symbols you'll see in the doc:

- **∈** — means "is in the range of" or "falls between." For example, `D ∈ [0, 1]` just means D can be any value from 0 to 1.
- **Greek letters** — these are just variable names, like nicknames for values:
  - **α** (alpha) — a scaling multiplier
  - **β** (beta) — a weighting factor that controls how much influence something has
  - **σ²** (sigma squared) — the variance, which tells you how spread out values are from the average
  - **μ** (mu) — the average (mean) of a set of values
  - **Φ** (Phi) — the standard normal cumulative distribution function (basically: "what percentile does this value fall at on a bell curve?")
  - **Δ** (delta) — means "change in" or "difference"
- **Σ** (capital sigma) — means "add up all the values." For example, `Σ wᵢ x xᵢ` means "multiply each weight by its matching score, then add them all together."
- **log₁₀** — the base-10 log. It answers "10 raised to what power gives me this number?" For example, log₁₀(1000) = 3 because 10³ = 1000. We use logs to compress big ranges of values into more manageable ones.
- **e** — Euler's number (~2.718), a mathematical constant used in sigmoid ("S-shaped") curves.
- **clamp(value, 0, 1)** — if the value goes below 0, force it to 0; if it goes above 1, force it to 1. Think of it as a buffer that keep a number from going too high.
- **min(a, b)** — whichever value is smaller.

---
