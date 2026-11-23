# Adaptive TASS Weighting (opt-in)

By default TaxTriage scores every organism in a sample with **one** weight vector, chosen per sample type and platform from `assets/sampletype_best_thresholds.json` (or set on the command line). This **TASS (default)** is always computed. `--adaptive_tass` adds a second, **adaptive TASS** per organism: its weights are shifted according to the situation the organism is in **within that sample**, and ambiguous calls with weak unique support are marked down.

In the pipeline outputs (`tass_score`, cutoffs, `.odr.txt` / PDF reports) TASS (default) stays the score; the adaptive score is written beside it. In the HTML report you can switch the whole report to the adaptive score (see below).

It is **off by default**. With it off, the pipeline's command lines, its `-resume` cache and every output are unchanged.

```bash
nextflow run jhuapl-bio/taxtriage ... --adaptive_tass --enable_matrix
```

`--enable_matrix` is optional but needed for the `ani` term (the main reason to use this), see [Terms](#terms).

---

## Why

A fixed vector is fitted to "the average sample of this type". It cannot tell the difference between:

- an organism that is alone in a clean sample, and the same numbers sitting next to a near-identical relative (Vaccinia vs Cowpox, _E. coli_ vs _Shigella_), where broad coverage is expected for both and only **unique** (MinHash) support separates them;
- a dominant organism and a trace call in a complex community, where read count and breadth say less;
- a single organism and a species/genus roll-up built from several children.

The adaptive layer moves weight towards the evidence that is diagnostic in each of those situations.

---

## Terms

Each organism gets a context from the other calls at the same rank (strain, species, genus) in the same sample:

| Term            | Context                                                                                                             | Weight shift                                     | Also                      |
| --------------- | ------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------ | ------------------------- |
| `complexity`    | Sample richness × evenness (Shannon, at species level). Active above 0.45.                                          | + MinHash, + Gini                                | advisory cutoff +         |
| `ani`           | Organisms of a **different species** with ANI ≥ `--adaptive_tass_high_ani` (default `--ani_threshold`) are present. | + MinHash, − breadth, + Gini                     | ambiguity penalty         |
| `sibling`       | Other species of the same genus are present (their share of the genus reads).                                       | + MinHash, − breadth (strain level, share > 25%) | smaller ambiguity penalty |
| `low_abundance` | Under 3% of the reads at its rank.                                                                                  | + MinHash, + Gini                                |                           |
| `rollup`        | A species/genus built from ≥ 2 children.                                                                            | + breadth, + Gini, − MinHash                     | advisory cutoff −         |

Weights are shifted, clipped at 0 and re-normalised so the weighted core (breadth, MinHash, Gini, disparity, HMP, MAPQ) keeps its original total. Plasmid bonus, abundance confidence (only reduced, never raised) and the score power are left as they are. The organism is then scored with the normal TASS formula and multiplied by

$$
\text{penalty} = 1 - p_{\max}\cdot s\cdot \text{ambiguity}\cdot(1 - \text{minhash\_reduction})
$$

where $\text{ambiguity} = 0.7\cdot\text{competitor share}\cdot\text{ANI severity} + 0.3\cdot\text{sibling share}$, $s$ is `--adaptive_tass_strength` and $p_{\max}$ is `--adaptive_tass_max_penalty` (0.35). A call with strong unique support keeps its score even next to a close relative; a weakly supported one drops.

**Strains of the same species never count against each other.** Which strain is present is an _assignment_ question, already reported (without changing TASS) by the strain-assignment step ([TASS Scoring §11](tass-scoring.md#11-aggregation-levels)). Counting them here would mark down a real species whose reads are split over several reference assemblies.

With the pipeline's grouping (`--rank genus`) the species is each organism's `subkey` and the sibling group its genus. If `match_paths.py` is run with `--rank species`, the top-level group is used as the species and the `sibling` term is inactive; the choice is recorded as `species_field` / `genus_field` in `metadata.adaptive_tass`.

The shifts and penalty are ported from the `adaptive_tass_weighting` prototype (`adaptive_weighting.py`). They are heuristics, not fitted values; see [Calibration](#calibration).

---

## Parameters

| Parameter                     | Default      | Meaning                                                                                                                                                                                                         |
| ----------------------------- | ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--adaptive_tass`             | `false`      | Also compute the adaptive score.                                                                                                                                                                                |
| `--adaptive_tass_strength`    | `1.0`        | Scale of every shift and the penalty, 0–1. `0` makes the adaptive score equal the standard score.                                                                                                               |
| `--adaptive_tass_terms`       | `all`        | Comma-separated subset of `complexity,ani,sibling,low_abundance,rollup`, or `all` / `none`.                                                                                                                     |
| `--adaptive_tass_base`        | `sampletype` | Starting weights for the adaptive score from the thresholds JSON: `sampletype` (same `type\|platform` entry as the standard score), `platform` (`all\|<platform>`), `global` (`all\|all`, as in the prototype). |
| `--adaptive_tass_high_ani`    | `null`       | ANI fraction for competitors; `null` = `--ani_threshold`. Needs `--enable_matrix`.                                                                                                                              |
| `--adaptive_tass_max_penalty` | `0.35`       | Largest multiplicative mark-down at strength 1 (max 0.9).                                                                                                                                                       |

Out-of-range values stop the run immediately with a message (they are never silently clamped).

`--adaptive_tass_base` only changes where the adaptive score starts. The standard score always uses the usual sample-type weights and cutoffs; explicit weight flags apply to both.

---

## Outputs

Nothing is removed or renamed, and every existing value is what a run without `--adaptive_tass` writes. With `--adaptive_tass` the per-sample JSON gains:

- **`tass_score`** — unchanged: the standard score (this is what the report filters on).
- **`tass_score_adaptive`** — the adaptive score, on the same 0–100 scale, with the same region-masking penalty and its own monotone roll-up floor (a species / genus is never below its best child).
- **`adaptive_tass`** on every organism — the audit trail: `standard` (the standard raw score), `static` (the raw score under the adaptive base weights, the same as `standard` unless `--adaptive_tass_base` is set) and `adaptive` raw scores, `penalty`, the `weights` used, `weight_shift`, `complexity`, `relative_abundance`, `max_ani`, `competitors` (up to five, with ANI and reads), `competitor_pressure`, `sibling_pressure`, `ambiguity`, `n_children`, `cutoff_shift` and plain-language `reasons`. `status` is `applied`, `static` (expected fallback) or `error`.
- **`metadata.adaptive_tass`** — settings, `primary_score: "standard"`, sample complexity, whether ANI and MinHash comparison were available, per-rank counts against the standard score (changed / raised / lowered, mean and max change) and errors.

### In the HTML report

When any sample was run with `--adaptive_tass`, the report shows an **Adaptive TASS** section in the right-hand filter panel:

- **TASS used in report**: `TASS (default)` or `Adaptive TASS`. This applies to the **whole report**: every tab, table, chart, filter, cutoff check, species / genus roll-up and data export uses the chosen score. With `Adaptive TASS`, each organism that has an adaptive score uses it; organisms (or samples) without one keep TASS (default). The cutoffs do not change (the suggested cutoff shift is never applied). An **Adaptive TASS** pill shows in the banner while it is in use (click it for the comparison). Saved sessions remember the choice; the records inside a session file are always stored with TASS (default).
- **Heatmap: show change (adaptive − default)**: colours the heatmap (Value = TASS Score) by the difference. A masking before / after heatmap view, when present, turns this off and is always computed on TASS (default).
- **Show adaptive TASS beside TASS in tables** (or **Show TASS (default) beside TASS** while adaptive is in use): a small chip with the other score after each TASS in the Summary table (`A▴/A▾` = adaptive, `D▴/D▾` = default), green / red when that score would flip the call at the current cutoff.
- **Hover**: any TASS in the Summary or Table tab, or a heatmap cell, shows both scores (marking the one in use), the ambiguity mark-down and the reasons.
- A one-line summary: how many strains score higher / lower adaptively, and how many calls adaptive TASS adds or drops at the current cutoffs.
- **Adaptive vs TASS (default)** opens a comparison: what the adaptive terms do, per-sample settings and change counts, a scatter of TASS (default) (x) against adaptive (y) with the cutoff, and a sortable table (pass / fail under each score, ambiguity mark-down, suggested cutoff shift, reasons) with CSV export.

The Table tab also gets **TASS (adaptive)**, **TASS Change (adaptive)**, **Adaptive Penalty**, **Adaptive Cutoff Shift** and **Adaptive Reasons** columns. Samples run without the flag have none of these values, and the section is hidden when no sample has them.

To inspect one sample from the command line:

```bash
python bin/adaptive_tass.py summarize results/alignment/<sample>.paths.json --level all
```

### Cutoffs are not moved

The prototype also raises or lowers the pass/fail cutoff per organism. That shift is computed and reported (`cutoff_shift`, in TASS points) but **not applied**: the report keeps the sample's usual cutoff. It is there to be calibrated first.

---

## Safety rails

- Off by default; when off, `match_paths.py` never imports the module and the task script is byte-identical.
- On or off, the standard `tass_score` of every organism is identical (checked on `Miseq_Run_A`: 34/34 organisms, every other output field unchanged).
- `strength 0` and `terms none` make the adaptive score equal the standard score (tested).
- Adapted weights are checked on every organism: finite, non-negative, core total unchanged, non-core weights untouched.
- If the MinHash reference comparison did not run, the whole sample keeps its static weights (`metadata.adaptive_tass.reason`), because every shift and the penalty depend on it.
- Any error for one organism keeps that organism's static score; errors are counted in `metadata.adaptive_tass.errors` and logged.
- With `--optimize`, the weights are fitted on static scores; a warning says the adaptive layer is applied on top.

---

## Measured effect

Replays of the bundled Kraken2 test samples (`Nanopore_Run_A`, stool / ONT and `Miseq_Run_A`, nasal / Illumina) with the pipeline's own `ALIGNMENT_PER_SAMPLE` inputs:

- Without `--adaptive_tass`: output identical to the unmodified code.
- `--adaptive_tass --adaptive_tass_strength 0`: every TASS identical; only the audit fields are added.
- `--adaptive_tass`: the standard TASS equals the static run's TASS for 33/33 and 34/34 organisms. Without an ANI matrix the changes are small (complexity / low-abundance / roll-up terms only): most organisms move < 1 point, the largest −5.7 (a minor _E. coli_ O157:H7 strain with weak MinHash support, 70.7 → 65.0, which crosses the stool/ONT strain cutoff of 69.6).
- With a high-ANI pair injected into the ANI matrix (_E. coli_ ↔ _Salmonella_ at 0.975, a stand-in for _E. coli_ / _Shigella_), the weakly unique members of the pair drop 8–17 points and the reasons name the competitor.

On the prototype's own scenarios, scored with TaxTriage's real formula and the `all|all` weights: the Cowpox decoy next to Vaccinia goes from 78.6 (passes 75) to 62.0 while Vaccinia stays at 88; the _Shigella_ decoy in a complex stool sample drops the most while _E. coli_ and _Salmonella_ are unchanged.

---

## Calibration

The magnitudes are heuristics. Before relying on adaptive scores (or applying the cutoff shift), fit them on labelled data — spike-ins, simulations or adjudicated calls — the same way the sample-type weights were fitted. Every input and output needed for that is in the per-organism `adaptive_tass` block.

Tests: `python3 bin/test_adaptive_tass.py`.
