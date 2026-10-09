# Alignment Trends

Alignment Trends answers one question across many samples: **which stretches of a
reference keep coming up empty, thin, or piled up?** A gap that recurs in most
samples that cover a reference well usually means the samples do not carry that
part of the assembly: a deletion, a divergent or novel locus, or an assembly
artefact. A pile-up that recurs points at repeats, rRNA operons, mobile
elements, plasmid copy number or loci that attract contaminating reads.

It is available in two places that give the same numbers for the same settings:

| Where              | How                                                     | Output                                                                                          |
| ------------------ | ------------------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| Interactive report | **Trends → Alignment Trends** sub-tab                   | _Live analysis_ (adjustable) and _Pipeline results_ (as the run computed them)                  |
| Pipeline           | on by default (`--alignment_trends false` turns it off) | `<outdir>/alignment_trends/` tables, JSON, XLSX, PNGs, and the report's _Pipeline results_ view |
| Standalone         | `bin/alignment_trends.py -i *.paths.json -o prefix`     | same as the pipeline                                                                            |

## How it works

1. **Depth profile.** For every reference (strain key) `match_paths.py` writes a
   `depth_profile` into the per-sample JSON: windowed mean depth and breadth for
   _every_ accession of the reference, including ones with no reads, each in its
   own coordinates. The window is the smallest `100 × 2^k` bp that gives at most
   `--depth_profile_windows` (default 400) windows over the reference, so all
   samples of a reference share a grid. The format is documented in
   `bin/depth_profile.py`; it adds well under 10% to a JSON.
2. **Normalise.** A sample's window depth is divided by that sample's mean depth
   on the reference, so a 2× library and a 200× library compare directly.
3. **Classify each window** per sample: **zero** (no read touches it), **low**
   (below `low_frac` × the sample mean, default 0.2), **high** (above `high_frac` ×
   the mean, default 3), or normal. Absolute cutoffs (`--trend_low_abs`,
   `--trend_high_abs`) override the relative ones.
4. **Only count samples that could have seen the window.** A sample counts toward
   a window when it expects at least `min_reads_per_window` reads there (its
   reads × window length ÷ reference length, default 20). With fewer, an empty
   window is Poisson noise rather than a gap. This also keeps tiny contigs from
   generating spurious "gaps".
5. **Frequencies and regions.** Per window: the fraction of counted samples that
   are zero, low-or-zero, and high. Consecutive windows (within one contig) where
   a fraction is at least `min_freq` (default 50%) merge into a recurrent region.

## Report sub-tab (Trends → Alignment Trends)

The sub-tab has two views.

**Pipeline results** (shown whenever the pipeline step ran, i.e. unless `--alignment_trends false`):
`make_report.py` embeds `all.alignment_trends.json` (`--alignment_trends_json`), and
this view reports the comparison exactly as the pipeline computed it — the settings
used, totals, the per-reference summary, every recurrent region (filter by type and
reference, download TSV) and the per-sample rows. Clicking a reference or region
opens it in the live view with the pipeline's settings applied and filters off.

**Live analysis** recomputes the same thing in the browser from the per-sample
depth profiles, so every cutoff can be changed. It opens with more permissive
defaults than the pipeline so shallow samples still show up — **Window × 4**,
**min reads 3**, **min reads / window 1** (pipeline: window × 1, 10 reads, 20 reads
per window); the other cutoffs match. Opening a reference from _Pipeline
results_ switches the controls to the pipeline's settings.

- **Reference** picker, ordered by how many samples count; the _All references_
  table at the bottom summarises every reference and opens one on click.
- **Heatmap** — one row per sample, log2(depth ÷ sample mean); black = zero,
  faded = the sample does not count there, `*` = below the expected-reads cutoff.
- **Frequency track** — fraction of counted samples zero / low-or-zero / high,
  recurrent regions shaded, grey where too few samples count. Drag to zoom.
- **Recurrent regions** table — click a row to zoom to it; download as TSV.
- **Respect filters** (on by default) limits the analysis to detections that pass
  the report's active filters and TASS cutoff; hidden samples are always excluded.

## Pipeline parameters

| Parameter                                | Default   | Meaning                                                              |
| ---------------------------------------- | --------- | -------------------------------------------------------------------- |
| `--depth_profile_windows`                | 400       | windows per reference in the JSON; 0 disables profiles (and the tab) |
| `--alignment_trends`                     | true      | run the cross-sample analysis; `false` turns it off                  |
| `--trend_min_samples`                    | 2         | counted samples needed per reference / window                        |
| `--trend_min_reads`                      | 10        | reads a sample needs on a reference to count                         |
| `--trend_min_reads_per_window`           | 20        | expected reads per window for a sample to count there                |
| `--trend_low_frac` / `--trend_high_frac` | 0.2 / 3.0 | relative low / high cutoffs                                          |
| `--trend_low_abs` / `--trend_high_abs`   | –         | absolute cutoffs (override the relative ones)                        |
| `--trend_min_freq`                       | 0.5       | recurrence fraction for a region                                     |
| `--trend_min_region_windows`             | 1         | minimum windows per region                                           |
| `--trend_plots`                          | 10        | PNGs for the top N references                                        |
| `--trend_matrix`                         | false     | also write the per-sample × window matrix                            |

The analysis uses the real (non-control, non-simulated) samples.

### Cost

The step reads only the per-sample JSONs, so it adds little to a run. On a
9-sample test set (43 MB of JSON, including a 293 Mb / 26,577-scaffold draft
assembly) it took about 5 s and 265 MB of RAM with 10 plots, and wrote about
2.8 MB of tables (the windows TSV and xlsx are most of that; the JSON embedded in
the report is about 160 KB). The depth profiles add well under 10% to a
per-sample JSON and no measurable time to `match_paths.py`. Draft assemblies with
more than 500 scaffolds that got no reads list those scaffolds only by count and
total bp, so a fragmented reference does not bloat the JSON; the tab reports
them as "no reads in any sample".

## Outputs (`<outdir>/alignment_trends/`)

| File                               | Contents                                                                                                    |
| ---------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `all.alignment_trends.summary.tsv` | one row per reference: samples counted, window, % of the reference in recurrent zero / low / high regions   |
| `all.alignment_trends.regions.tsv` | recurrent regions: type, contig, start, end, mean / max frequency, affected samples, `whole_contig`         |
| `all.alignment_trends.windows.tsv` | per window: counted samples, zero / low / high counts and frequencies, normalised depth summary             |
| `all.alignment_trends.samples.tsv` | per sample × reference: reads, mean depth, breadth, expected reads per window, counted, % windows per class |
| `all.alignment_trends.json`        | all of the above (minus the matrix)                                                                         |
| `all.alignment_trends.xlsx`        | summary / regions / samples / windows sheets (when openpyxl is available)                                   |
| `all.alignment_trends.plots/`      | heatmap + frequency PNG per top reference (when matplotlib is available)                                    |

## Running it on existing results

Any set of TaxTriage per-sample JSONs (or combined `all.odr.json` files) made
with a version that writes depth profiles can be analysed together, including
samples from different runs:

```bash
bin/alignment_trends.py -i run1/alignment/*.paths.json run2/report/all.odr.json \
    -o cohort.alignment_trends --plots 20 --min_reads_per_window 10
```

Dropping those JSONs onto any TaxTriage report also fills the Alignment Trends tab.
