# Detection Rescue & Report Surfacing

TaxTriage hides low-confidence detections by default: an organism whose [TASS score](tass-scoring.md) falls below the cutoff is dropped from the tables, charts, and heatmap so the report stays readable. But "below the cutoff" is not the same as "absent" — a sub-threshold organism may still have reference reads aligned, carry virulence/AMR genes, or be placed by the reference-free [novelty](novelty-detection.md) branch.

**Rescue** is the set of mechanisms that bring those suppressed-but-supported detections back into view *without* changing any score. Every rescue is opt-in or clearly badged, lives only in the detection tables (the Summary tab's Detections view and the full Table tab), and never alters the KPIs, charts, or heatmap, which stay restricted to passing detections.

This page covers the three rescue toggles, the differentiating UI each uses, and the related surfacing updates to the interactive report. For the upstream novelty workflow itself (backends, databases, params) see [Novelty Detection](novelty-detection.md).

---

## The three rescue toggles

All three live in the report sidebar under the filter controls. They are independent and can be combined; rows surfaced by more than one are de-duplicated.

| Toggle (sidebar label) | Default | What it surfaces |
|---|---|---|
| **Roll up threshold** (`filter-rollup`) | off | A strain kept visible when its parent **species or genus aggregate** passes the cutoff, even though the strain's own score does not. |
| **Show below-cutoff organisms with VF/AMR hits** (`filter-below-vfamr`) | off | A below-cutoff organism whose **genus has a Virulence-Factor or AMR gene hit** in the same sample. Requires `--annotate`. |
| **Show sub-threshold & novelty-supported organisms** (`filter-novelty-sub`) | off | Below-cutoff organisms that still have **reference alignments**, and no-alignment organisms with a **genus/species novelty hit** in the same sample. Shown whenever novelty data is present. |

### 1. Threshold roll-up

A high-confidence call at the genus or species level can vanish simply because no single strain underneath it cleared the cutoff. With roll-up on, a strain is kept visible when its parent species *or* genus aggregate passes — so a strong genus call is never erased by strain-level dilution. Rescued rows are tagged with a coloured **↑ rollup** badge naming the level that rescued them (`↑ Species` / `↑ Genus`), and an amber bar in the charts.

### 2. Below-cutoff VF/AMR

When protein annotation is enabled (`--annotate`) and one or more VF or AMR genes were detected for an organism's genus **in the same sample**, this re-surfaces that organism as a faded row even though its score is below the cutoff — flagging potential pathogen signal the score alone would suppress. These rows carry a pink **↓ below cutoff** badge and a red rail.

### 3. Sub-threshold & novelty-supported organisms

This rescue handles the case the other two miss: a sample that *does* have passing detections but also contains rows with sub-threshold or alignment-free signal. Consider a sample with three organisms — one above cutoff with alignments, one below cutoff with alignments, and one with no alignments at all that the novelty classifier nonetheless places at the genus/species level. The first is already shown; this toggle surfaces the other two:

| Surfaced row | Condition | Badge | Rail |
|---|---|---|---|
| **Sub-threshold, aligned** (`__belowCutoffAligned`) | Below the cutoff, but `# Reads Aligned > 0`. | <span>↓ sub-threshold (aligned)</span> (blue) | blue |
| **Novelty, no alignment** (`__noveltyNoAlign`) | `# Reads Aligned = 0`, but the active novelty backend placed this organism's **genus or species** in the same sample. | ✦ novelty, no alignment (orange) | orange |

The novelty match reuses the same per-row logic as the Novelty column, so the badge and the column's signal always agree. Rows already shown by the VF/AMR toggle are skipped to avoid duplicates, and all the usual sidebar filters (text, category, High-Consequence, molecular type, view level) apply. The badge tooltip explains exactly why each row is visible (its TASS, aligned-read count, or the novelty rank that backs it).

> Like every rescue, these rows appear **only** in the detection tables. The heatmap, TASS chart, coverage plots, and KPI cards stay restricted to detections that pass the cutoff.

---

## Annotation of unaligned samples

VF/AMR annotation (`--annotate`) runs on the **de novo contigs** of every sample, producing `annotate/<sample>.annotate_report.tsv`. Previously these hits only reached the comparison report when they could be attached to an *aligned* organism, so a sample with **no reference alignment** (e.g. a shallow dilution that assembled contigs but mapped nothing) lost its annotation entirely — the Novelty tab showed signal, but the Summary tab's Annotation Summary and the VF/AMR views were blank for that sample.

The standalone `annotate_report.tsv` files are now plumbed directly into the comparison report:

- The workflow collects every per-sample `annotate_report.tsv` and passes it to `CREATE_COMPARISON_REPORT` (`subworkflows/local/report.nf` → `modules/local/create_comparison_report.nf`).
- `bin/make_report.py` gains `--annotate_reports`. For any sample **not already covered** by the merged annotation XLSX, it builds supplemental `per_gene_hits` / `amr_genes` / `genus_summary` rows from that sample's TSV, routing AMR vs. virulence by property, and lets the existing pathogen-by-taxid stamping resolve the canonical genus.

The result: unaligned samples now show their VF/AMR annotation in the Summary tab, the VF/AMR tab, and the per-row Annotation column, consistent with the Novelty tab. This is automatic whenever `--annotate` is set — there is no extra user flag. (`--annotate_reports` is an internal `make_report.py` argument the module supplies; samples already represented by the merged XLSX are never double-counted.)

---

## Related report surfacing updates

These are interactive-report changes that make rescued and edge-case data easier to read. See [Interactive Report](interactive-report.md) for the full tour.

### Right-panel shows every sample

The sample list in the right panel is now built from the authoritative per-sample metadata, not just from samples that produced detections. Samples with no hits, no alignments, or everything filtered out (e.g. shallow dilutions) now appear in the list so they can be coloured, hidden, renamed, or inspected like any other.

### Low-reads / no-alignment warning icon

Each sample row that has few reads or no alignments shows an amber ⚠ icon. Hovering it explains that the sample's plots and scores may be empty or unreliable, so a missing bar or flat heatmap row reads as "expected for this depth" rather than a bug.

### Cross-sample charts: detail on hover

The cross-sample views moved their numbers off the plot and onto hover to cut clutter:

- **Feature Compare matrix** — cells are colour swatches only; hovering shows the exact value, the metric description, and the per-sample breakdown.
- **Co-occurrence matrix** — cell tooltips now also list the shared / detected sample names.
- **Genera Comparison Across Samples** — a **Show values** checkbox toggles numeric labels back on. Per-segment numbers are drawn only where they fit (so they never overlap), the per-genus total sits once at the bar end, and the chart scrolls horizontally if a bar runs wide. With values off, the full breakdown is on hover.

---

*See also: [Novelty Detection](novelty-detection.md) · [TASS Scoring](tass-scoring.md) · [Interactive Report](interactive-report.md) · [CLI Parameters](cli-parameters.md)*
