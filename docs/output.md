# Output

This page describes the files and directories produced by TaxTriage. All paths are relative to `--outdir` (your specified output folder).

---

## Most Important Files

These are the outputs you should review first after a successful run:

| File                   | Location                              | Description                                                                                  |
| ---------------------- | ------------------------------------- | -------------------------------------------------------------------------------------------- |
| Combined ODR PDF       | `report/all.odr.pdf`                  | All-sample Organism Discovery Report with confidence-ranked pathogen table                   |
| Per-sample ODR PDF     | `report/<sample>.odr.pdf`             | Single-sample ODR                                                                            |
| Interactive Comparison | `report/all.odr.html`                 | Multi-sample interactive comparison report. See [Interactive Report](interactive-report.md). |
| MultiQC Report         | `report/multiqc_report.html`          | QC stats and alignment metrics across all samples                                            |
| Krona Plot             | `report/combined_krona_kreports.html` | Interactive radial abundance visualization from Kraken2                                      |
| Microbial Sheet        | `report/<sample\|all>.odr.txt`        | Full tabular data backing the ODR PDF                                                        |
| Annotated Workbook     | `report/<sample\|all>.odr.xlsx`       | Excel workbook of the ODR data, including VF/AMR annotation and a Metadata sheet             |

---

## Directory Structure

```
<OUTDIR>/
├── report/
│   ├── all.odr.pdf                       # Combined Organism Discovery Report
│   ├── <sample>.odr.pdf                  # Per-sample ODR
│   ├── all.odr.html                      # Multi-sample interactive comparison report
│   ├── all.odr.txt                       # Combined microbial sheet
│   ├── all.odr.xlsx                      # Combined annotated workbook (incl. VF/AMR + Metadata)
│   ├── <sample>.odr.txt                  # Per-sample microbial sheet
│   ├── <sample>.odr.xlsx                 # Per-sample annotated workbook
│   ├── multiqc_report.html               # MultiQC report
│   └── combined_krona_kreports.html      # Krona plot
│
├── alignment/
│   ├── <sample>.<sample>.dwnld.references.bam      # Post-alignment BAM (+ .csi index)
│   ├── <sample>.histo.txt                          # Coverage histogram
│   ├── <sample>.paths.json                         # Per-sample data backing the interactive report
│   └── <sample>_removal_stats_by_taxid.xlsx        # Conflict / LCA read-removal stats per taxid
│
├── samtools/
│   └── <sample>.txt                      # Per-reference coverage
│
├── bcftools/                             # (--reference_assembly only)
│   ├── <sample>.<taxid>.vcf.gz          # Variant calls
│   └── <sample>.consensus.fa            # Consensus assembly
│
├── fastqc/                              # (Illumina)
│   ├── *_fastqc.html
│   └── *_fastqc.zip
│
├── nanoplot/                            # (ONT)
│
├── mergedkrakenreport/
│   └── krakenreport.merged_mqc.tsv      # Top hits per sample from Kraken2
│
├── top/
│   ├── <sample>.top_report.tsv          # Per-sample top hits table
│   ├── <sample>.toptaxids.txt           # Top taxids selected for alignment
│   └── <sample>.topnames.txt            # Top organism names
│
├── download/
│   └── <sample>.dwnld.references.fasta  # Downloaded reference sequences
│
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    ├── pipeline_dag.svg
    └── software_versions.yml
```

---

## Output File Details

### FastQC (`fastqc/`)

Standard FastQC outputs for Illumina samples:

- `*_fastqc.html` — Interactive quality metrics report
- `*_fastqc.zip` — Archived report with raw data

> Note: FastQC in the MultiQC report shows **untrimmed** reads and may contain adapter sequences.

**Sequence count distribution:**
![FastQC sequence counts](images/mqc_fastqc_counts.png)

**Mean quality scores:**
![FastQC mean quality](images/mqc_fastqc_quality.png)

**Adapter content:**
![FastQC adapter content](images/mqc_fastqc_adapter.png)

### MultiQC (`report/multiqc_report.html`)

Aggregated report across all samples, including:

- FastQC / NanoPlot summaries
- Alignment statistics from samtools
- Kraken2 classification summary
- Software version traceability

### Organism Discovery Report PDF

The main deliverable. Example report:

<img src="../images/pathogens.report.example.png" width="60%">

Each table row is one detected organism with the following key columns:

| Column                   | Description                                                                                                     |
| ------------------------ | --------------------------------------------------------------------------------------------------------------- |
| **Detected Organism**    | Organism name                                                                                                   |
| **Specimen ID**          | NCBI Taxid                                                                                                      |
| **Sample Type**          | body site (blood, stool, etc.)                                                                                  |
| **% Reads**              | Percentage of total reads assigned                                                                              |
| **# Reads Aligned**      | Reads mapped above minimum MAPQ                                                                                 |
| **% Aligned Reads**      | Fraction of dataset aligned                                                                                     |
| **Coverage**             | Fraction of reference genome covered                                                                            |
| **HHS Percentile**       | Abundance rank vs. HMP healthy human dataset                                                                    |
| **IsAnnotated**          | Present in pathogen sheet                                                                                       |
| **Microbial Category**   | Unknown / Commensal / Primary / Opportunistic / Potential — see [Microbial Categories](microbial-categories.md) |
| **High Consequence**     | Always shown in PDF regardless of TASS cutoff                                                                   |
| **Gini Coefficient**     | Coverage distribution inequality (0=uneven, 1=uniform)                                                          |
| **Mean BaseQ**           | Average base quality across aligned reads                                                                       |
| **Mean MapQ**            | Average mapping quality of aligned reads                                                                        |
| **Mean Coverage**        | Mean coverage across all contigs/chromosomes                                                                    |
| **Mean Depth**           | Mean depth across all positions                                                                                 |
| **Minhash Score**        | False-positive penalty score                                                                                    |
| **Breadth Weight Score** | Log-based breadth of coverage component                                                                         |
| **TASS Score**           | Final combined confidence score (0–1)                                                                           |
| **K2 Reads**             | Number of Kraken2-assigned reads                                                                                |
| **Siblings Score**       | Read abundance vs. related organisms in same genus                                                              |

See [TASS Scoring](tass-scoring.md) for full definitions of each metric.

### Interactive Comparison Report (`report/all.odr.html`)

A self-contained, browser-based report that compares every sample in the run side by side, with a TASS heatmap, summary table, coverage/sunburst/explore views, a per-sample-type TASS cutoff slider, species/genus roll-up views, and a built-in Export-to-PDF button. No server is required — the file can be emailed or hosted as-is. See the dedicated [Interactive Report](interactive-report.md) page for a full walkthrough.

### Microbial Sheet (`report/<sample>.odr.txt`)

Tabular plain-text version of the ODR data. Contains all organisms including commensals and potentials regardless of `--show_commensals` / `--show_potentials` flags.

### Annotated Workbook (`report/<sample>.odr.xlsx`)

Excel version of the same data with VF/AMR protein-annotation columns (when annotation is enabled) and an appended `Metadata` sheet describing the run.

### BAM Files (`alignment/<sample>.<sample>.dwnld.references.bam`)

Post-alignment BAM files (with a `.csi` index) for all top-hit organisms. Can be viewed in IGV or any BAM-compatible viewer.

### Per-sample JSON (`alignment/<sample>.paths.json`)

Structured per-sample results that back the interactive comparison report. These are the preferred input to `make_report.py`.

### Read-removal Stats (`alignment/<sample>_removal_stats_by_taxid.xlsx`)

Per-taxid breakdown of how many reads were removed during conflict resolution, including the LCA-aware (species/genus) accounting used when scoring closely related strains. See [TASS Scoring](tass-scoring.md#11-aggregation-levels).

### Coverage File (`samtools/<sample>.txt`)

Per-reference coverage summary produced by `samtools coverage` for each top-hit organism.

### Variant Files (`bcftools/`) — Optional

Available when `--reference_assembly` is enabled:

- VCF files with variant calls per taxid
- Consensus FASTA assembled from the reference + variants

---

## Pathogen Sheet

The pipeline uses a curated pathogen annotation sheet with ~1,600 taxa at [`assets/pathogen_sheet.csv`](https://github.com/jhuapl-bio/taxtriage/blob/main/assets/pathogen_sheet.csv).

### Updating the Pathogen Sheet

You can modify the sheet or provide your own with `--pathogens`. Required columns:

1. `name` — Organism name (doesn't need to match NCBI exactly)
2. `taxid` — NCBI taxonomy ID
3. `general_classification` — `Primary`, `Opportunistic`, `Potential`, or `Commensal` (see [Microbial Categories](microbial-categories.md) for definitions, site-aware resolution, and interpretation)
4. `high_consequence` — `TRUE`/`FALSE` — always shown in PDF regardless of TASS score

Optional but recommended:

- `pathogenic_sites` — Comma-separated body sites where this organism is pathogenic
- `commensal_sites` — Sites where this organism is commensal (overrides general classification for those sites)
- `assembly_accession` — Curated `GCF_*` / `GCA_*` assembly to pin for this organism. When set, it is downloaded directly (accession-first) instead of the pipeline picking an assembly by taxid; see [Assembly Selection Order](cli-parameters.md#assembly-selection-order). Auto-generated during the database build and kept as the last column; leave blank to let TaxTriage choose.

To request new organisms be added to the default sheet, open a GitHub issue.

---

## Working Directory (`work/`)

Nextflow stores all intermediate process files in `work/`. These are **not** final outputs but are used when resuming a pipeline with `-resume`.

To free disk space after a successful run:

```bash
rm -rf work/
```

> ⚠️ Deleting `work/` prevents future `-resume` from reusing cached results.
