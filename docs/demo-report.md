---
title: Demo Report
hide:
  - toc
---

<style>.md-content__inner{max-width:none}</style>

# Interactive Report Demo

A live, fully interactive TaxTriage report built from the example dataset — the
same artifact the pipeline writes to `multiqc/heatmap.html` at the end of a run.
Nothing is uploaded and nothing is computed server-side; the whole report is a
single self-contained HTML file carrying its own data.

[Open in a new tab :material-open-in-new:](demo/index.html){ .md-button .md-button--primary target=\_blank }
[Read the report guide](interactive-report.md){ .md-button }

<!-- Raw HTML is not path-rewritten by MkDocs the way Markdown links are, and
     pages are served at directory URLs, so this must step up a level. -->
<div class="demo-frame">
  <iframe src="../demo/index.html" title="TaxTriage interactive report demo" loading="lazy"></iframe>
</div>

## What you're looking at

| Tab                   | Shows                                                                                       |
| --------------------- | ------------------------------------------------------------------------------------------- |
| **Summary**           | Per-sample overview with the top-scoring organisms and QC flags                             |
| **Heatmap**           | Organism × sample matrix; switch the cell value between TASS, abundance, coverage and depth |
| **TASS**              | Score decomposition — how each component contributed to the final confidence                |
| **Coverage**          | Per-reference coverage tracks with the conflict regions marked                              |
| **Sunburst**          | Taxonomic composition, drillable from kingdom down to species                               |
| **Proteins / VF-AMR** | Virulence factor and antimicrobial resistance gene hits                                     |
| **Explore**           | Free-form filtering across every metric in the run                                          |

Full documentation of each panel lives in [Interactive Report](interactive-report.md),
and the scoring behind the heatmap values is described in [TASS Scoring](tass-scoring.md).

## Try the sample QC rules

The demo opens with two **sample QC rules** already active, so the feature is
visible rather than hidden behind a dialog. They flag three of the six samples —
one for shallow read depth, three for having fewer than four distinct organisms
above TASS 99 — and you can see the effect as an outlined column in the Heatmap,
a tinted group header in the Table and Summary, and a badge in Metadata &
Mapping. Hover any marker for the exact reason.

Open **Sample QC / Flags → Filters & flags…** in the right-hand panel to edit
them: change a threshold and watch the flagged list update as you type, add a
rule on any metadata column or detection metric, or switch a rule from _flag it_
to _flag & hide it_ to drop those samples from every view. **Remove all rules**
clears them entirely.

!!! tip "These thresholds are illustrative"

    They are tuned to trip on the example dataset so the markers are visible.
    In a real run the defaults come from the pipeline's `--report_flag_*`
    parameters — see [Report Sample-QC Flags](cli-parameters.md#report-sample-qc-flags).

!!! note "This is example data"

    The demo is built from the pipeline's test dataset so the layout and
    interactions are representative, but the organisms and scores are not from a
    real clinical sample. Run the pipeline on your own data to generate a report
    like this — see [Quick Start](quick-start.md).
