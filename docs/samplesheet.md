# Samplesheet

TaxTriage uses a comma-separated samplesheet (`.csv`) to specify the input samples. This page describes the format, required and optional columns, metadata support, and sample types.

---

## Basic Format

The samplesheet must have a header row and at least the mandatory columns:

```
sample,platform,fastq_1,fastq_2,sequencing_summary,trim,type
NB03,OXFORD,examples/data/fastq_demux/NB03,,,FALSE,nasal
BC05_flu,OXFORD,examples/data/BC05.fastq.gz,,,FALSE,nasal
longreads,OXFORD,examples/data/nanosim_metagenome.fastq.gz,,,FALSE,gut
shortreads,ILLUMINA,examples/data/iss_reads_R1.fastq.gz,examples/data/iss_reads_R2.fastq.gz,,TRUE,blood
```

> ⚠️ **Sample names must be unique per row.** Spaces are automatically converted to underscores.

---

## Column Reference

| Column               | Required                            | Description                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| -------------------- | ----------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `sample`             | **MANDATORY**                       | Unique sample name. If the same name appears on multiple rows, those reads are concatenated before processing. Spaces are converted to `_`.                                                                                                                                                                                                                                                                                                            |
| `platform`           | **MANDATORY**                       | Sequencing platform: `ILLUMINA`, `OXFORD`, or `PACBIO`. Defaults to `ILLUMINA` if blank. Other platforms may work as a proxy but are not officially supported.                                                                                                                                                                                                                                                                                         |
| `bam`                | OPTIONAL                            | Path to a pre-aligned `.bam`, `.cram` or `.sam` file. When set, the row skips QC, trimming, host removal, classification and reference download, and `fastq_1` is ignored. A sheet must have either this column or `fastq_1`. See [Pre-Aligned Input](#pre-aligned-input-bam-column).                                                                                                                                                                  |
| `fastq_1`            | **MANDATORY** (unless `bam` is set) | Path to the primary input file. Accepts gzipped or uncompressed FASTQ (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`) **or** FASTA (`.fa`, `.fasta`, `.fna`, `.fa.gz`, `.fasta.gz`, `.fna.gz`). May also be a directory of FASTQ files (ONT demux mode), **or an SRA/ENA accession** to download automatically. See [SRA / ENA Accessions](#sra--ena-accessions), [FASTA Input](#fasta-input) and [Multi-File Input](#multi-file-input) for extended options. |
| `fastq_2`            | OPTIONAL                            | Path to the second FASTQ file for paired-end Illumina reads. Must be gzipped. Ignored when `fastq_1` is FASTA, contains multiple files, or is an accession.                                                                                                                                                                                                                                                                                            |
| `sequencing_summary` | OPTIONAL                            | Path to a Nanopore sequencing summary file — enables pycoQC plots.                                                                                                                                                                                                                                                                                                                                                                                     |
| `trim`               | OPTIONAL                            | `TRUE` or `FALSE` — whether to run adapter trimming on this sample.                                                                                                                                                                                                                                                                                                                                                                                    |
| `type`               | OPTIONAL (Recommended)              | Sample body site (e.g., `blood`, `stool`). Used for HMP-based filtering and pathogen annotation. See [Sample Types](#sample-types) below.                                                                                                                                                                                                                                                                                                              |
| `minimap2_preset`    | OPTIONAL                            | Override the minimap2 alignment preset for this sample. See [minimap2 Presets](#minimap2-presets) below.                                                                                                                                                                                                                                                                                                                                               |
| `positive`           | OPTIONAL                            | Name of the positive control sample (matching `sample` column) to compare against.                                                                                                                                                                                                                                                                                                                                                                     |
| `negative`           | OPTIONAL                            | Name of the negative control sample to compare against.                                                                                                                                                                                                                                                                                                                                                                                                |

---

## SRA / ENA Accessions

Instead of a file path, `fastq_1` may hold a public sequence archive accession. TaxTriage resolves it, downloads the reads, and feeds them into the pipeline exactly as if you had supplied local FASTQs. `fastq_2` is not used — whether a run is paired-end is read from the archive, not guessed.

```
sample,platform,fastq_1,fastq_2,sequencing_summary,trim,type
covid_run,,SRR13191702,,,FALSE,nasal
ont_survey,,ERR1160846,,,FALSE,stool
whole_project,,PRJNA681875,,,FALSE,stool
```

### Accepted accession types

| Level           | Prefixes                                       | Result                                              |
| --------------- | ---------------------------------------------- | --------------------------------------------------- |
| Run             | `SRR`, `ERR`, `DRR`                            | One sample, keeping the name in the `sample` column |
| Experiment      | `SRX`, `ERX`, `DRX`                            | One sample per child run                            |
| Sample          | `SRS`, `ERS`, `DRS`, `SAMN`, `SAMEA`, `SAMD`   | One sample per child run                            |
| Study / project | `SRP`, `ERP`, `DRP`, `PRJNA`, `PRJEB`, `PRJDB` | One sample per child run                            |

When an accession expands to more than one run, each resulting sample is named `<sample>_<run_accession>` (e.g. `whole_project_SRR13191701`) so names stay unique and traceable. Every other column on that row — `trim`, `type`, `positive`/`negative`, `minimap2_preset`, metadata columns — is inherited by all of its runs.

### Paired-end detection

Nothing needs to be declared. The pipeline reads ENA's file listing for the run: two `_1`/`_2` FASTQs means paired-end, one file means single-end. This is checked against the actual files rather than the submitter's declared `library_layout`, so a run that claims to be `PAIRED` but only ships one FASTQ is correctly treated as single-end. On the sra-tools path the same conclusion comes from what `fasterq-dump --split-3` writes out.

### Platform

Leave `platform` blank and the instrument platform reported by the archive is used — `ILLUMINA`, `OXFORD` (Oxford Nanopore) or `PACBIO` — which in turn drives the default minimap2 preset. Filling in the `platform` column overrides that.

### How the download works

1. **ENA first.** ENA mirrors most of SRA as ready-made `fastq.gz`, already split into `_1`/`_2`. These are downloaded directly and verified against ENA's md5 — no `.sra` archive, no conversion step.
2. **sra-tools fallback.** If ENA has no files for the run (very recent submissions, unmirrored runs), the accession is expanded via NCBI eutils and fetched with `prefetch` + `fasterq-dump`. Force this path for everything with `--sra_force_sratools`.

### Caching

Downloads are written to `<outdir>/sra_downloads/<run_accession>/` and are **skipped entirely** on any later run where they already exist — both with `-resume` and on a fresh run. Set `--sra_cache_dir /shared/sra` to place the cache on shared storage so an accession is only ever downloaded once across projects or users.

> ℹ️ A file that happens to be _named_ like an accession (e.g. `SRR13191702.fastq.gz`) is still treated as a local path. Only bare accessions — no dots, no slashes — trigger a download.

> ⚠️ **Controls and multi-run accessions.** The `positive`/`negative` columns match on final sample names. If a control is itself given as a project/experiment accession, its samples get the `_<run_accession>` suffix, so point the control columns at the suffixed name (or use a run accession for controls so the name stays exactly what you wrote).

---

## Pre-Aligned Input (`bam` column)

If you have already aligned your reads, add a `bam` column and give the path to a `.bam`, `.cram` or `.sam` file instead of a path in `fastq_1`. Those samples skip QC, trimming, host removal, classification and reference download, going straight to coverage statistics and TASS scoring.

```
sample,platform,bam,type
prealigned_a,ILLUMINA,/data/prealigned_a.bam,stool
prealigned_b,OXFORD,/data/prealigned_b.cram,nasal
```

A samplesheet must have **either** a `fastq_1` column or a `bam` column. When both are present, a row that fills in `bam` is pre-aligned and its `fastq_1` value is ignored, so the two input styles can be mixed in one sheet:

```
sample,platform,fastq_1,bam,type
raw_sample,ILLUMINA,/data/raw_R1.fastq.gz,,stool
prealigned,ILLUMINA,,/data/prealigned.bam,stool
```

The FASTQ rows run the full pipeline; only the pre-aligned rows take the short path. `trim` is forced to `false` on pre-aligned rows, and the file does not need to be sorted or indexed beforehand.

Supply `--reference_fasta` (the reference the alignment was made against) whenever you have it. Without it TaxTriage reconstructs the reference from the alignment, which works but biases conflict-based read removal. See [Pre-Aligned Input](pre-aligned-input.md) for the full guide, including CRAM handling and the consensus parameters.

---

## FASTA Input

In addition to FASTQ files, `fastq_1` accepts FASTA-format sequence files. This is useful when you already have assembled or consensus sequences and want to skip read-level quality control.

**Supported FASTA extensions:** `.fa`, `.fasta`, `.fna`, `.fa.gz`, `.fasta.gz`, `.fna.gz`

### What changes for FASTA samples

FASTA inputs **skip** these steps automatically (no flag needed):

- Adapter trimming (TrimGalore / Porechop)
- Read compression (pigz)
- Read counting
- FASTP quality filtering
- FastQC / NanoPlot quality plots

FASTA inputs **still run** through:

- Host removal (minimap2 handles FASTA queries natively)
- Kraken2 classification (if not disabled)
- Reference alignment (minimap2 / Bowtie2 / HISAT2)
- All downstream scoring and reporting

> ℹ️ If host removal is configured (`--remove_reference_file` or `--genome`), unaligned FASTA sequences are converted to FASTQ internally by `samtools fastq` before continuing downstream. Quality scores in the converted file are placeholder values and are intentionally excluded from FastQC / NanoPlot.

### Example

```
sample,platform,fastq_1,type,minimap2_preset
ont_assembly,OXFORD,results/assembly.fa,respiratory,map-ont
hifi_consensus,PACBIO,results/hifi_consensus.fasta,blood,map-hifi
illumina_contigs,ILLUMINA,results/spades_contigs.fa,gut,asm5
rna_iso,OXFORD,results/iso_seq.fa,lung,splice:hq
```

---

## Multi-File Input

Multiple input files for a single sample can be listed in `fastq_1` as a **semicolon-separated** (`;`) list. This is commonly used when a sample spans multiple sequencing runs, flow cells, or batches, or when passing multiple query FASTA files to minimap2 (e.g. RNA-seq splice presets).

### Rules

- All files in the list must be the **same type** — all FASTA or all FASTQ. Mixing raises a validation error.
- Multi-file inputs are always treated as **single-end** (`fastq_2` is ignored if present).
- Each path is validated individually at startup.

### FASTA example (multiple assemblies as one sample)

```
sample,platform,fastq_1,minimap2_preset
iso_rna,OXFORD,run1.fa;run2.fa;run3.fa,splice
ont_multi,OXFORD,flowcell1.fa;flowcell2.fa;flowcell3.fa,map-ont
assembly_comparison,ILLUMINA,asm_chr1.fa;asm_chr2.fa,asm5
```

### FASTQ example (multiple runs merged into one sample)

```
sample,platform,fastq_1
hifi_combined,PACBIO,pass1.fastq.gz;pass2.fastq.gz;pass3.fastq.gz
ont_combined,OXFORD,run_a.fastq.gz;run_b.fastq.gz
```

> ℹ️ For multi-file FASTQ, adapter trimming (TrimGalore / FASTP) expects 1–2 files and may not behave as intended. Set `trim=FALSE` for multi-file FASTQ samples if trimming is not needed.

---

## minimap2 Presets

The optional `minimap2_preset` column sets the minimap2 alignment preset **per sample**, overriding the platform-derived default. This applies to both host removal and the main alignment step.

When this column is absent or blank, the preset is selected automatically from `platform`:

| Platform                  | Default preset |
| ------------------------- | -------------- |
| `ILLUMINA`                | `sr`           |
| `PACBIO`                  | `map-hifi`     |
| `OXFORD` (and all others) | `map-ont`      |

### Available presets

| Preset            | Best for                                                |
| ----------------- | ------------------------------------------------------- |
| `map-ont`         | Oxford Nanopore genomic reads                           |
| `map-pb`          | PacBio CLR genomic reads                                |
| `map-hifi`        | PacBio HiFi / CCS reads (v2.19+)                        |
| `lr:hq`           | Nanopore Q20+ high-quality reads (v2.27+)               |
| `sr`              | Short-read paired-end genomic reads                     |
| `splice`          | Spliced long reads (strand unknown)                     |
| `splice -uf -k14` | Noisy Nanopore direct RNA-seq (pass as `task.ext.args`) |
| `splice:hq`       | PacBio Kinnex / Iso-seq (RNA-seq)                       |
| `splice:sr`       | Short-read RNA-seq (v2.29+)                             |
| `asm5`            | Intra-species assembly-to-assembly alignment            |
| `asm10`           | Cross-species assembly alignment                        |
| `asm20`           | Divergent assembly alignment                            |
| `ava-pb`          | PacBio read overlap                                     |
| `ava-ont`         | Nanopore read overlap                                   |

> ℹ️ Presets that require additional flags beyond the preset name (e.g. `splice -uf -k14`) should pass the extra flags via `task.ext.args` in your Nextflow config rather than in the `minimap2_preset` column.

### Example

```
sample,platform,fastq_1,type,minimap2_preset
direct_rna,OXFORD,direct_rna.fa,lung,splice
iso_seq,PACBIO,iso_seq.fq.gz,lung,splice:hq
rna_short,ILLUMINA,r1.fq.gz;r2.fq.gz,gut,splice:sr
ont_q20,OXFORD,q20_reads.fq.gz,blood,lr:hq
asm_vs_ref,ILLUMINA,assembly.fa,sterile,asm5
```

---

## Using Controls

You can specify positive and negative lab controls per sample:

```
sample,platform,fastq_1,fastq_2,sequencing_summary,trim,type,negative,positive
shortreads,ILLUMINA,examples/data/iss_reads_R1.fastq.gz,examples/data/iss_reads_R2.fastq.gz,,TRUE,blood,Negative Control Miseq,Positive Control Miseq
Positive Control Miseq,ILLUMINA,examples/data/controls/positive/pos_R1.fastq.gz,examples/data/controls/positive/pos_R2.fastq.gz,,FALSE,nasal,,
Negative Control Miseq,ILLUMINA,examples/data/controls/negative/neg_R1.fastq.gz,examples/data/controls/negative/neg_R2.fastq.gz,,FALSE,nasal,,
```

> ❗ A positive or negative control entry whose `sample` value doesn't appear in any other row is ignored — no comparison will occur for that row.

---

## Sample Types

The `type` column influences downstream filtering and annotation. Supported values:

| Value     | HMP Filtering | Notes                       |
| --------- | ------------- | --------------------------- |
| `blood`   | No            |                             |
| `csf`     | No            | Cerebrospinal fluid         |
| `sterile` | No            |                             |
| `nasal`   | No            |                             |
| `oral`    | Yes           | Matches HMP oral dataset    |
| `vaginal` | Yes           | Matches HMP vaginal dataset |
| `gut`     | Yes           | Matches HMP gut dataset     |
| `stool`   | Yes           | Matches HMP stool dataset   |
| `skin`    | Yes           | Matches HMP skin dataset    |
| `abscess` | No            |                             |
| `lung`    | No            |                             |
| `ear`     | No            |                             |
| `sputum`  | No            |                             |
| `urine`   | No            |                             |

> If a type is not in this list or is left blank, the HMP distribution filtering step is skipped for that sample. The type still affects pathogen annotation labels in the final PDF.

---

## Metadata Support

You can enrich reports with metadata in two ways:

### Option 1: Extra Columns in the Samplesheet

Add any additional columns beyond the standard ones. They will appear in the report and can be used for filtering. Example:

```
sample,platform,fastq_1,type,altitude,collection_date
Sample_A,ILLUMINA,reads_R1.fastq.gz,blood,1200,2024-03-12
```

**Grouping DNA + RNA libraries into one specimen.** Give both libraries the same `specimen` value; the interactive report can then merge them (see the **Cross-Sample Organism Analysis** tab's _Merge_ toggle) so prevalence counts the specimen once:

```
sample,platform,fastq_1,type,specimen
Patient01_DNA,ILLUMINA,p01_dna_R1.fastq.gz,nasal,Patient01
Patient01_RNA,ILLUMINA,p01_rna_R1.fastq.gz,nasal,Patient01
Patient02_DNA,ILLUMINA,p02_dna_R1.fastq.gz,nasal,Patient02
```

### Option 2: Separate Metadata CSV (`--meta`)

Provide a separate CSV file with at minimum a `sample` column matching your samplesheet:

```bash
--meta /path/to/metadata.csv
```

> ❗ Values in `--meta` take **priority** over samplesheet columns if there are conflicts.

### Recognized Metadata Columns

Beyond `sample`, several column names are **recognized** by the interactive report and power dedicated views (Map, Longitudinal, Geographic Comparison, Host & Disease). Names are matched exactly (snake_case), so spelling matters. See [Interactive Report → Run Metadata tab](interactive-report.md#run-metadata-tab-and-metadata-driven-views) for what each one drives.

| Column                                                                                                      | Used in Multi-Run Report | Description                                                                                                                                                                                                                                                                                                                                                                         |
| ----------------------------------------------------------------------------------------------------------- | ------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`                                                                                                    | Yes                      | Must match samplesheet sample names                                                                                                                                                                                                                                                                                                                                                 |
| `run_id`                                                                                                    | Yes                      | Groups samples by run in the ODR                                                                                                                                                                                                                                                                                                                                                    |
| `specimen`                                                                                                  | Yes                      | Groups multiple samples (e.g. a DNA and an RNA library from the same swab) into one **specimen**. Samples sharing a `specimen` value can be merged in the interactive report's **Cross-Sample Organism Analysis** tab so prevalence and per-organism TASS / coverage / reads are computed per specimen rather than per sample. Leave blank for samples that are their own specimen. |
| `latitude`                                                                                                  | Yes                      | Sample collection latitude; enables the **Map** tab and Geographic Comparison                                                                                                                                                                                                                                                                                                       |
| `longitude`                                                                                                 | Yes                      | Sample collection longitude (pair with `latitude`)                                                                                                                                                                                                                                                                                                                                  |
| `collection_time`                                                                                           | Yes                      | Drives **Longitudinal Analysis**. Format: `YYYY-MM-DD HH:MM:SS` (e.g., `2024-03-12 09:45:00`) or `M/D/YYYY`                                                                                                                                                                                                                                                                         |
| `sample_origin_country`                                                                                     | Yes                      | Country for the **Geographic Comparison** choropleth                                                                                                                                                                                                                                                                                                                                |
| `sample_origin_state_province_territory`                                                                    | Yes                      | State / province / territory for the **Geographic Comparison** choropleth (default level)                                                                                                                                                                                                                                                                                           |
| `host_scientific_name`                                                                                      | Yes                      | Host species; a grouping in **Host & Disease**                                                                                                                                                                                                                                                                                                                                      |
| `host_disease`                                                                                              | Yes                      | Host disease / symptoms; default **Host & Disease** grouping and the **Symptom × Organism** matrix. Multi-value: comma-separated symptoms are split                                                                                                                                                                                                                                 |
| `environmental_site`                                                                                        | Yes                      | Environmental sampling site; a grouping in **Host & Disease**                                                                                                                                                                                                                                                                                                                       |
| `depth`                                                                                                     | Yes                      | Numeric (m); shown in the metadata table                                                                                                                                                                                                                                                                                                                                            |
| `salinity`                                                                                                  | Yes                      | Numeric (PSU); shown in the metadata table                                                                                                                                                                                                                                                                                                                                          |
| `location`                                                                                                  | Yes                      | Free-text location; shown in the metadata table                                                                                                                                                                                                                                                                                                                                     |
| `sequencing_instrument`, `sequencing_platform`, `library_preparation_kit`, `sequencing_protocol_primer_set` | Yes                      | Sequencing metadata; rendered with friendly labels                                                                                                                                                                                                                                                                                                                                  |
| Any other column                                                                                            | Yes                      | Included in report for filtering/display                                                                                                                                                                                                                                                                                                                                            |

An example ODR showing metadata in action: [https://jhuapl-bio.github.io/taxtriage/](https://jhuapl-bio.github.io/taxtriage/)

---

## Running Without a Samplesheet (Single Sample)

For a quick single-sample run, skip the CSV and provide FASTQ files directly:

```bash
nextflow run main.nf \
  -profile test,docker -resume \
  --fastq_1 examples/data/iss_reads_R1.fastq.gz \
  --fastq_2 examples/data/iss_reads_R2.fastq.gz \
  --type blood \
  --sample my_sample
```

The pipeline will create a temporary samplesheet internally.

`--fastq_1` also accepts an accession, in which case the reads are downloaded first:

```bash
nextflow run main.nf \
  -profile docker -resume \
  --fastq_1 SRR13191702 \
  --type nasal \
  --outdir results
```

Paired-end layout and platform come from the archive, so neither `--fastq_2` nor `--platform` is needed. If `--sample` is omitted the accession itself is used as the sample name. A project accession here fans out into one sample per run, same as in a samplesheet.

| Parameter    | Required | Description                                                          |
| ------------ | -------- | -------------------------------------------------------------------- |
| `--fastq_1`  | Yes      | Path to primary input file (FASTQ or FASTA), or an SRA/ENA accession |
| `--fastq_2`  | No       | Path to second FASTQ (paired-end). Ignored for accessions            |
| `--sample`   | No       | Sample name (defaults to FASTQ basename)                             |
| `--type`     | No       | Sample type (e.g., `blood`)                                          |
| `--platform` | No       | Platform (defaults to `ILLUMINA`)                                    |

---

## Example Samplesheet

A full example is provided in the repository at [`examples/Samplesheet.csv`](https://github.com/jhuapl-bio/taxtriage/blob/main/examples/Samplesheet.csv).

---

## Next Steps

- [Running the Pipeline](running-the-pipeline.md) — how to pass your samplesheet to TaxTriage
- [CLI Parameters](cli-parameters.md) — full parameter reference
