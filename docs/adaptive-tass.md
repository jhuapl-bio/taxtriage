# Adaptive TASS Weighting (opt-in)

By default TaxTriage scores every organism in a sample with **one** weight vector, chosen per sample type and platform from `assets/sampletype_best_thresholds.json` (or set on the command line). `--adaptive_tass` adds a per-organism layer on top: each organism's weights are shifted according to the situation it is in **within that sample**, and ambiguous calls with weak unique support are marked down.

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

| Parameter                     | Default      | Meaning                                                                                                                                                                                                                  |
| ----------------------------- | ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--adaptive_tass`             | `false`      | Turn the layer on.                                                                                                                                                                                                       |
| `--adaptive_tass_strength`    | `1.0`        | Scale of every shift and the penalty, 0–1. `0` reproduces the static scores exactly.                                                                                                                                     |
| `--adaptive_tass_terms`       | `all`        | Comma-separated subset of `complexity,ani,sibling,low_abundance,rollup`, or `all` / `none`.                                                                                                                              |
| `--adaptive_tass_base`        | `sampletype` | Starting weights **and cutoffs** from the thresholds JSON: `sampletype` (the usual `type\|platform` entry), `platform` (`all\|<platform>`), `global` (`all\|all`, as in the prototype). Explicit weight flags still win. |
| `--adaptive_tass_high_ani`    | `null`       | ANI fraction for competitors; `null` = `--ani_threshold`. Needs `--enable_matrix`.                                                                                                                                       |
| `--adaptive_tass_max_penalty` | `0.35`       | Largest multiplicative mark-down at strength 1 (max 0.9).                                                                                                                                                                |

Out-of-range values stop the run immediately with a message (they are never silently clamped).

With `--adaptive_tass_base platform|global`, `tass_score_static` is the static score **with that base**, not the sample-type default; run once without `--adaptive_tass` for the default static scores.

---

## Outputs

Nothing is removed or renamed. With `--adaptive_tass` the per-sample JSON gains:

- **`tass_score`** — the adaptive score (this is what the report filters on).
- **`tass_score_static`** — exactly what a run without `--adaptive_tass` reports for that organism (the monotone roll-up floor is mirrored).
- **`adaptive_tass`** on every organism — the audit trail: `static` / `adaptive` raw scores, `penalty`, the `weights` used, `weight_shift`, `complexity`, `relative_abundance`, `max_ani`, `competitors` (up to five, with ANI and reads), `competitor_pressure`, `sibling_pressure`, `ambiguity`, `n_children`, `cutoff_shift` and plain-language `reasons`. `status` is `applied`, `static` (expected fallback) or `error`.
- **`metadata.adaptive_tass`** — settings, sample complexity, whether ANI and MinHash comparison were available, per-rank counts (changed / raised / lowered, mean and max change) and errors.

The interactive report adds **TASS (static)**, **TASS Change (adaptive)**, **Adaptive Penalty**, **Adaptive Cutoff Shift** and **Adaptive Reasons** to the Table tab.

To inspect one sample from the command line:

```bash
python bin/adaptive_tass.py summarize results/alignment/<sample>.paths.json --level all
```

### Cutoffs are not moved

The prototype also raises or lowers the pass/fail cutoff per organism. That shift is computed and reported (`cutoff_shift`, in TASS points) but **not applied**: the report keeps the sample's usual cutoff. It is there to be calibrated first.

---

## Safety rails

- Off by default; when off, `match_paths.py` never imports the module and the task script is byte-identical.
- `strength 0` and `terms none` give exactly the static scores (tested).
- Adapted weights are checked on every organism: finite, non-negative, core total unchanged, non-core weights untouched.
- If the MinHash reference comparison did not run, the whole sample keeps its static weights (`metadata.adaptive_tass.reason`), because every shift and the penalty depend on it.
- Any error for one organism keeps that organism's static score; errors are counted in `metadata.adaptive_tass.errors` and logged.
- With `--optimize`, the weights are fitted on static scores; a warning says the adaptive layer is applied on top.

---

## Measured effect

Replays of the bundled Kraken2 test samples (`Nanopore_Run_A`, stool / ONT and `Miseq_Run_A`, nasal / Illumina) with the pipeline's own `ALIGNMENT_PER_SAMPLE` inputs:

- Without `--adaptive_tass`: output identical to the unmodified code.
- `--adaptive_tass --adaptive_tass_strength 0`: every TASS identical; only the audit fields are added.
- `--adaptive_tass`: `tass_score_static` equals the static run's TASS for 33/33 and 34/34 organisms. Without an ANI matrix the changes are small (complexity / low-abundance / roll-up terms only): most organisms move < 1 point, the largest −5.7 (a minor _E. coli_ O157:H7 strain with weak MinHash support, 70.7 → 65.0, which crosses the stool/ONT strain cutoff of 69.6).
- With a high-ANI pair injected into the ANI matrix (_E. coli_ ↔ _Salmonella_ at 0.975, a stand-in for _E. coli_ / _Shigella_), the weakly unique members of the pair drop 8–17 points and the reasons name the competitor.

On the prototype's own scenarios, scored with TaxTriage's real formula and the `all|all` weights: the Cowpox decoy next to Vaccinia goes from 78.6 (passes 75) to 62.0 while Vaccinia stays at 88; the _Shigella_ decoy in a complex stool sample drops the most while _E. coli_ and _Salmonella_ are unchanged.

---

## Calibration

The magnitudes are heuristics. Before relying on adaptive scores (or applying the cutoff shift), fit them on labelled data — spike-ins, simulations or adjudicated calls — the same way the sample-type weights were fitted. Every input and output needed for that is in the per-organism `adaptive_tass` block.

Tests: `python3 bin/test_adaptive_tass.py`.
