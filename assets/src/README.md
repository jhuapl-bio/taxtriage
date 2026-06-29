# Interactive report — source layout

`assets/heatmap.html` is a **thin shell** (~3.4k lines): the `<head>`, the HTML
body markup (banner, tabs, panels, overlays), and `<link>` / `<script src>`
references to the external CSS and JS in this directory. It is directly
openable in a browser (with `assets/heatmap_boot.js` present for demo data).

The previous single ~30k-line file lives on only as the *generated*
self-contained output; it is no longer committed.

## What's here

```
src/css/
  geo.css        <style id="geo-host-extra-css"> — geo/leaflet host styles
  main.css       the main stylesheet
  loading.css    <style id="tt-loading-css"> — loading overlay
  compare.css    comparison-view styles
src/js/
  early.js       small <head> script — reads window.HEATMAP_BOOT
  00_preamble.js … 36_guided_help_tour.js
                 the main bundle, one file per author section (heatmap, tass,
                 sunburst, coverage, proteins/VF-AMR, table, histogram, explore,
                 summary, cross-sample, compare-coverages, longitudinal, map,
                 file upload, session state, help tour, …)
```

The `src/js/*.js` files load in filename order as separate classic `<script>`
tags. Classic scripts share one global scope, so they call into each other
freely — they are split for readability, not isolated as ES modules. Don't add
duplicate top-level `const`/`let`/`class` names across files (that throws across
separate scripts), and keep the numeric ordering correct.

## How the one-file outputs are produced

`bin/report_template.py:inline_template()` folds every `src/` reference back
inline (each `<link>` → `<style>`, each `<script src>` → `<script>`) to yield a
single portable file. Both builders call it:

| Builder | Output | Used by |
| --- | --- | --- |
| `bin/make_report.py` | `all.comparison.report.html` | nextflow pipeline |
| `scripts/inline_boot_json.py` | `_site/index.html` | GitHub Pages deploy |

In the pipeline, `assets/src/` is staged next to `heatmap.html` in the report
task's workdir (see `modules/local/create_comparison_report.nf` and
`subworkflows/local/report.nf`); the inliner also falls back to the template's
real path, so it resolves the parts whether they're staged or symlinked.

## Editing

Edit the shell (`assets/heatmap.html`) and/or the files here directly, then:

```bash
# preview the assembled self-contained file
python bin/report_template.py -o /tmp/report.html

# rebuild the GitHub Pages page
python scripts/inline_boot_json.py
```

`scripts/split_report_template.py` is the one-shot bootstrap that produced this
layout from the old monolith; you normally won't need it again.
