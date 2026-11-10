# Illumina Index Hopping

On patterned Illumina flow cells (NovaSeq, HiSeq X/4000, NextSeq 2000, …) a small fraction of reads swap an index with another library in the same lane. A sample can therefore pick up a few reads of an organism that is abundant in a neighbouring library. TaxTriage can read the demultiplexer's index-hopping table and **flag** every detection that this kind of bleed-through could explain.

It is a flag only: TASS scores, thresholds and calls are unchanged. Flagged detections are highlighted in the reports so a reviewer can judge them, and they can be hidden in the HTML report with one click.

---

## Quick start

```bash
nextflow run jhuapl-bio/taxtriage -r main \
  --input samplesheet.csv --outdir results -profile docker \
  --index_hopping Reports/Index_Hopping_Counts.csv
```

The table can be **CSV, TSV or XLSX** (`.gz` is fine for CSV/TSV), for example bcl-convert's `Reports/Index_Hopping_Counts.csv` or the same table exported from BaseSpace or Excel:

```text
Lane  SampleID             Sample_Project     index       index2      # Reads   % of Hopped Reads  % of All Reads
1     IDR2600053924-01-01  MTB-PHOENIX-11955  TGCCGGTCAG  TTGTATCAGG  1861265                      0.003531
1     IDR2600055130-01-01  MTB-PHOENIX-11955  CACTCAATTC  CACAACTTAA  1287738                      0.002443
1     IDR2600055131-01-01  MTB-PHOENIX-11955  TCTCACACGC  TCACTGTCCG  2116690                      0.004016
1     IDR2600055132-01-01  MTB-PHOENIX-11955  TCAATGGAGA  GTTGTAGTGC  3024655                      0.005739
1                          (hopped)           TGCCGGTCAG  CACAACTTAA  412       0.18               0.0000008
```

Column names are matched case-insensitively and title rows above the header are skipped (common in Excel exports). Only `# Reads` and either `index` or `SampleID` are required; `Lane`, `index2` and the rest are used when present. Thousands separators (`1,861,265`) are fine. Legacy `.xls` is not supported; save it as `.xlsx`, `.csv` or `.tsv`.

### Try it

`examples/index_hopping_demo.tsv` is a made-up table for the example run (`-profile test`, i.e. `examples/Samplesheet.csv`). It puts the three MiSeq samples in one lane with their read totals and no hopped rows, so the default 0.1 % rate applies. `Nanopore_Run_A` is left out (ONT has no index pairs) and is reported as not matched.

```bash
nextflow run jhuapl-bio/taxtriage -r main -profile test,docker --outdir test_output \
  --index_hopping examples/index_hopping_demo.tsv
# or on an existing run:
bin/index_hopping.py -i test_output/alignment/*_Miseq.paths.json test_output/alignment/Miseq_Run_A.paths.json \
  --hopping examples/index_hopping_demo.tsv -o annotated
```

Or, without rerunning anything, open an existing `report/all.odr.html` and drop the table on its upload zone (see [In the HTML report](#in-the-html-report)).

Ten of the negative control's low-level hits (e.g. 10 _P. aeruginosa_ reads where ~17 could hop in from the 16,568 in `Miseq_Run_A`) are flagged. Hits no other sample carries in bulk, such as its 43 _S. cerevisiae_ reads, are not.

---

## Which samples are compared

Each TaxTriage sample is matched to a `SampleID` in the table:

1. the samplesheet `index_sample_id` column, when set;
2. otherwise the `sample` name, ignoring case and punctuation (`IDR2600053924-01-01` = `idr2600053924_01_01`), and also with a bcl-convert FASTQ suffix removed (`IDR2600053924-01-01_S1_L001` → `IDR2600053924-01-01`).

Only samples in **the same lane of the same table** are compared with each other. Samples from different flow cells (different tables) never are. Negative and positive controls are included, since a negative control is the classic place for hopped reads to show up. Samples that are not in a table get `metadata.index_hopping.applied = false` and nothing else changes for them.

### Several runs in one TaxTriage run

Give each row its own table with the samplesheet `index_hopping` column. It overrides `--index_hopping` for that row; several tables for one row are `;`-separated:

```csv
sample,platform,fastq_1,fastq_2,type,index_hopping,index_sample_id
S1,ILLUMINA,s1_R1.fq.gz,s1_R2.fq.gz,sputum,runA/Index_Hopping_Counts.csv,
S2,ILLUMINA,s2_R1.fq.gz,s2_R2.fq.gz,sputum,runA/Index_Hopping_Counts.csv,
S3,ILLUMINA,s3_R1.fq.gz,s3_R2.fq.gz,sputum,runB/hopping.xlsx,IDR2600055131-01-01
NEG,ILLUMINA,neg_R1.fq.gz,neg_R2.fq.gz,sputum,runB/hopping.xlsx,IDR2600055132-01-01
```

---

## How a detection is judged

For a detection (strain, species or genus) with **N** reads in sample **R**, every other sample **D** of the lane that has _more_ reads of the same taxon is a possible donor:

```text
expected  = Σ over donors D:  reads_D(taxon) × rate(D → R)
p         = P(X ≥ N),  X ~ Poisson(factor × expected)
flagged   = expected > 0  and  p ≥ alpha
```

So a detection is flagged when hopping **cannot be ruled out**. With the defaults (`--index_hopping_factor 2`, `--index_hopping_alpha 0.01`) and a donor with 50,000 reads at a 0.1 % rate (50 reads expected), 40 reads in the recipient are flagged and 5,000 are not.

The hop rate is taken from the table, most specific first:

| Rate source | When                                                                                               | Rate                                                                                                                      |
| ----------- | -------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `pair`      | The table lists hopped index combinations of D and R (i7 of one, i5 of the other)                  | (hopped reads of the two combinations + 1) ÷ (reads of D + reads of R)                                                    |
| `lane`      | The lane has hopped rows, but none for this pair                                                   | hopped ÷ (hopped + correctly paired reads) in the lane: the chance a read hops anywhere, an upper bound for one recipient |
| `default`   | The table has no hopped rows (only each sample's own pair, like the example above), or no `index2` | `--index_hopping_rate` (0.1 %)                                                                                            |

When two samples share several lanes, the largest of their lane rates is used. A row is a sample's _own_ (correctly paired) count when it has a `SampleID` and is that sample's pair with the most reads in the lane; every other row (blank `SampleID`, or another pair) is a hopped combination.

!!! note "Unique dual indexes"
With unique dual indexes (UDIs) a single hop produces an invalid index pair, which the demultiplexer discards. A read only lands in another sample if both indexes change, which is far rarer. The single-hop rates above are therefore an upper bound for UDI libraries, and suit combinatorial dual indexing. Because this is a flag, not a filter, erring on the side of flagging is deliberate; raise `--index_hopping_alpha` or lower `--index_hopping_factor` for fewer flags.

Read counts are the aligned reads TaxTriage reports for each sample, so host removal and QC affect donor and recipient alike. Taxa are matched at each level separately (strain to strain, species to species, genus to genus), so the species and genus rows also catch hops whose reads were assigned to a different strain reference in the recipient.

---

## Parameters

| Parameter                | Default | Description                                                                        |
| ------------------------ | ------- | ---------------------------------------------------------------------------------- |
| `--index_hopping`        | –       | Comma-separated index-hopping table(s) (CSV / TSV / XLSX) applied to every sample. |
| `--index_hopping_rate`   | `0.001` | Hop rate when a table has no hopped rows (0.1 %).                                  |
| `--index_hopping_factor` | `2`     | Multiplier on the expected hopped reads, for uncertainty in the rate.              |
| `--index_hopping_alpha`  | `0.01`  | Flag when P(reads ≥ observed \| hopping) ≥ alpha.                                  |

| Samplesheet column | Description                                                           |
| ------------------ | --------------------------------------------------------------------- |
| `index_hopping`    | This row's table(s); overrides `--index_hopping`. `;`-separated.      |
| `index_sample_id`  | This row's `SampleID` in the table, when it is not the `sample` name. |

The step runs only when `--index_hopping` is set or a samplesheet row has an `index_hopping` value; otherwise the pipeline is unchanged.

---

## Where the flags show up

- **HTML report (`report/all.odr.html`).** A purple _Index hop N · m/n_ badge in the banner (N flagged strain-level detections, m of n samples in the table) opens a per-sample summary; see [In the HTML report](#in-the-html-report). Table columns _Index Hop Suspect_, _Index Hop Expected Reads_, _Index Hop Explained %_, _Index Hop p_, _Index Hop Source_ (top donor), _Index Hop Source Reads_, _Index Hop Rate %_ and _Index Hop Rate Source_. A default Organism QC rule (`Index Hop Suspect == true → flag`) highlights the suspects; switch the Organism QC view to _Hide flagged_ to drop them from every tab, or _Only flagged_ to review them.
- **Per-sample PDF (`report/<sample>.odr.pdf`).** Flagged rows carry an orange note, e.g. `[index hop? ~44 of 90 reads expected from Sample_7]`. The `.odr.txt` sheet is unchanged.
- **JSON.** Every organism node gets an `index_hop` block (`suspect`, `observed_reads`, `expected_reads`, `p_value`, `explained_fraction`, `top_donor`, `top_donor_reads`, `top_donor_rate`, `rate_source`, and up to five `donors`), and each sample's `metadata.index_hopping` records the table, the matched `SampleID`, lanes, lane rates, donors and counts. These are in `report/all.odr.json` and `index_hopping/json/<sample>.paths.json`. `alignment/<sample>.paths.json` stays as scored.
- **`index_hopping/index_hopping.tsv`.** One row per flagged detection (sample, level, taxon, observed vs expected reads, p, top donor, rate and its source, TASS, call).
- **`index_hopping/index_hopping.json`.** The per-sample summary.

### Hover details

- **Icons.** In the HTML detection tables (Summary _Detections_ and the Table tab), a flagged row gets a small orange <i class="fas fa-shuffle"></i> and a possible source a faint <i class="fas fa-share"></i> after the organism name. The Detections legend explains both when hopping data is present.
- **Flagged-row card (<i class="fas fa-shuffle"></i>).** Reads seen vs. expected from hopping (% explained), the main source sample with its reads and TASS, the hop rate used and where it came from, the probability that hopping alone gives this many reads, and what to check next.
- **Source card (<i class="fas fa-share"></i>).** Each sample the reads may have hopped into, with the reads there, the share expected from this sample and the p-value there.
- **Column help.** Hovering an _Index Hop …_ column header in the Table tab explains the column.
- **PDF.** The Color Key explains the `[index hop?]` and `[possible hop source]` notes, only when the sample has them.

### The source side

The organism the reads likely hopped **from** gets a light, informational note: no flag, no highlight, TASS unchanged. Example: _P. aeruginosa_ in `Miseq_Run_A` when the negative control's few _P. aeruginosa_ reads are flagged. A donor is noted when it is expected to have sent at least 0.5 reads into a flagged detection in another sample.

- **HTML report.** A faint arrow (<i class="fas fa-share"></i>) after the organism name in the detection tables; its tooltip names the sample(s) and the expected reads. Two columns: _Index Hop Source Of_ (e.g. `Negative_Control_Miseq (~16.6)`) and _Index Hop Reads Out_. The badge panel counts these rows per sample under _Source of_.
- **Per-sample PDF.** A small grey note on the row, e.g. `[possible hop source → Negative_Control_Miseq (~18.8)]`.
- **JSON.** `index_hop_source` on the organism node (`n_recipients`, `expected_out`, `recipients` with `sample`, `expected_reads`, `observed_reads`) and `metadata.index_hopping.n_source` per level.

## In the HTML report

The report can also run the check itself, so a table can be applied to (or tried on) any report without rerunning the pipeline:

- **Load a table.** Drop the hopping table (CSV, TSV or XLSX) on the report's _Upload Data_ zone, or use _Load hopping table…_ in the badge's panel. A dropped file is recognised by its header (`# Reads` plus an `index` column), so it can be dropped together with `.paths.json` files: the detections are merged first and the hops computed over all of them. `.gz` tables have to be unzipped first in the browser.
- **Same result as the pipeline.** The browser uses the same parser, rates and test as `bin/index_hopping.py`, and matches taxa the same way (by level and taxonomic ID), so the columns agree to the last decimal. For a report the pipeline already annotated, it reuses each sample's matched `SampleID`, including samplesheet `index_sample_id` mappings.
- **Badge panel.** One row per sample: its `SampleID` in the table, lane(s), the hop rate used (`per pair`, `lane x %` or `default x %`), number of donor samples, and flagged rows at strain / species / genus level. It also lists samples not in the table, and how many flagged rows pass the current TASS cutoff and filters. A flagged row below the cutoff is still flagged, just not shown, which is why _Only flagged_ can come up empty.
- **Settings.** _Rate %_, _Factor_ and _Alpha_ (same meaning as the pipeline parameters) and _Recompute_ re-run the check on the loaded table. _Only flagged_ / _Hide flagged_ switch the Organism QC view, and _Clear table_ drops the browser result. On a pipeline-annotated report this returns to the pipeline's values.
- **Only loaded samples are compared.** `all.odr.html` includes the negative and positive controls (hidden or shown with the _Show controls_ switch, which never removes them from the check), so a table loaded in the browser compares them like the pipeline does. Reports from before this change left the controls out; drop `alignment/<control>.paths.json` on those to include them. Loading or clearing more `.paths.json` files recomputes automatically.

## Standalone use

The script runs on existing outputs too:

```bash
bin/index_hopping.py -i results/alignment/*.paths.json \
    --hopping Index_Hopping_Counts.csv \
    -o annotated --summary index_hopping.tsv --report index_hopping.json
bin/make_report.py -i annotated/*.paths.json -t assets/heatmap.html -o all.odr.html
```

Without `--sample-map` every table applies to every sample.
