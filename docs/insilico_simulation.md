# In-Silico Simulation Pipeline

This page describes the in-silico read simulation workflow in TaxTriage, from specification of params through simulated read generation, alignment, and final comparison metrics in the ODR.

## Overview

The in silico simulation pipeline makes synthetic (simulated) sequencing reads from the organisms detected in each samples' Kraken2 classification. The reads are treated as new samples that flow through the standard alignment pipeline (minimap2, bowtie2, or hisat2). Their alignment results are then compared against the non-control samples' results to compute precision, recall, F1, and accuracy.

Two simulators are supported and can be run either separately or together:

- [**InSilicoSeq (ISS)**](https://github.com/HadrienG/InSilicoSeq) — Illumina paired-end reads with realistic error profiles
- [**NanoSim**](https://github.com/bcgsc/NanoSim) — Oxford Nanopore long reads with configurable error models

When both are enabled, each produces a separate insilico sample per non control sample, and the report makes one metrics table per simulator.

## Parameters

### Simulation Control

| Parameter | Default | Type | Description |
|-----------|---------|-------------|-------------|
| `--generate_iss` | `false` | Boolean | Enable Illumina read simulation via InSilicoSeq |
| `--generate_nanosim` | `false` | Boolean | Enable ONT read simulation via NanoSim |
| `--sim_nreads` | `100000` | Integer | Total number of reads to simulate per sample |
| `--sim_nsamples` | `1` | Integer | Number of simulated sample variants to generate per real sample |
| `--sim_ranks` | `'S S1 S2 S3'` | String | Kraken2 taxonomic rank codes to include from the top report |
| `--sim_minreads` | `3` | Integer | Minimum `clade_fragments_covered` threshold for organism inclusion |
| `--sim_exclude_taxids` | `'9606'` | String | Comma-separated taxids to exclude (default excludes human) |
| `--sim_include_taxids` | `''` | String | Comma-separated taxids to force-include regardless of filters |
| `--sim_abundance` | none | Path | Optional custom abundance TSV file (bypasses Kraken2 report parsing) |
| `--sim_random_abundance` | `false` | Boolean | Use Dirichlet random sampling for abundances instead of observed proportions |

### ISS-Specific

| Parameter | Default | Type | Description |
|-----------|---------|-------------|-------------|
| `--iss_model` | `'miseq'` | String | ISS error profile model. Options: `miseq`, `hiseq`, `novaseq` |

### NanoSim-Specific

| Parameter | Default | Type | Description |
|-----------|---------|-------------|-------------|
| `--nanosim_training` | required | Path | Path to NanoSim training model directory. Required when `--generate_nanosim` is set |
| `--nanosim_base` | required, default is "training" | String | Base filename for training model files within the training directory. Check the basename of the files in the model folder when specifying |
| `--sim_ont_divisor` | `40` | Integer | ONT read count divisor: `ont_nreads = sim_nreads / sim_ont_divisor` |

## Pipeline Architecture

```
Real Sample
  ├── Kraken2 --> top_report.tsv
  ├── Reference Prep --> merged_taxid.tsv + reference FASTAs
  │
MAKE_SIMULATED_SAMPLES
  ├── Discovers ALL organisms from FASTA + merged_taxid (not just top 3 from Kraken2)
  ├── Builds per-accession abundance profile
  └── Outputs: abundance.tsv + reference.fasta per simulated sample
        │
        ├── INSILICOSEQ_SIMULATE (if --generate_iss)
        │      Produces: {sample}_insilico_iss_R1.fastq.gz, _R2.fastq.gz
        │      Tagged: platform=ILLUMINA, single_end=false
        │
        └── PREPARE_NANOSIM_INPUTS --> NANOSIM_SIMULATE (if --generate_nanosim)
               Produces: {sample}_insilico_nanosim*.fastq.gz
               Tagged: platform=OXFORD, single_end=true
                 │
        Injected into ch_reads as NEW SAMPLES
                 │
        ALIGNMENT (standard pipeline — same containers as real samples)
                 │
        REPORT (3-way branch: control / insilico / non-control)
          ├── Insilico --> ALIGNMENT_PER_SAMPLE_INSILICO --> insilico JSONs
          └── Non-control --> ALIGNMENT_PER_SAMPLE (receives insilico JSONs as --insilico_controls)
                              │
                       match_paths.py computes fold-changes & missing organisms
                              │
                       create_report.py renders per-type metrics tables
```

## Step-by-Step Details

### Step 1: Organism Discovery and Abundance Profile

The `make_simulated_samples.py` script takes 3 inputs per sample:

1. **top_report.tsv** — Kraken2 classification report with columns including `taxid`, `rank`, `clade_fragments_covered`, `abundance`, and `number_fragments_assigned`
2. **merged_taxid.tsv** — Reference prep mapping file with columns: `Acc`, `Assembly`, `Organism_Name`, `Description`, `Mapped_Value` (taxid)
3. **Reference FASTA files** — The downloaded reference sequences or those provided with the `--reference_fasta` param.

#### Organism Selection

The script uses a 2 step method to find all organisms to simulate:

**Source 1: Kraken2 top_report** — Organisms are included if they match the specified ranks (default: `S, S1, S2, S3`), have `clade_fragments_covered >= --sim_minreads`, and their taxid is not in the exclusion list.

**Source 2: FASTA-first discovery** — The `discover_organisms_from_fasta()` function cross-references the accessions present in the reference FASTA with the merged_taxid.tsv mapping. This captures organisms that have downloaded reference sequences but may not appear in the Kraken2 report at the expected rank level (e.g., strain level accessions when the report only has species level entries). Newly discovered organisms receive a default abundance of `max(1.0, --sim_minreads)`.

The two sources are merged, ensuring every organism with reference sequences gets simulated.

#### Abundance File Format
The output `abundance.tsv` is a file at the accession level:

```
NC_003310.1	0.15373765867418904
NZ_AP023069.1	0.061780265963376
NZ_AP023070.1	0.784482075362476
```

Values are abundances that sum to 1.0 (100%). When `--sim_random_abundance` is set, abundances are sampled from a standard distribution. Otherwise, they are the same as to the observed `number_fragments_assigned` numbers from the Kraken2 report.

### Step 2a: InSilicoSeq (Illumina Simulation)
ISS generates paired end Illumina reads using kernel density estimation (KDE) error models trained on real sequencing data.

**Container:** `biocontainers/insilicoseq:2.0.1`

**Command:**
```bash
iss generate \
    --genomes reference.fasta \
    --model miseq \
    --output {sample_id}.iss \
    --mode kde \
    --abundance_file abundance.tsv \
    -n 100000 \
    --cpus {cpus}
```

**Output:** `{sample_id}.iss_R1.fastq.gz` and `{sample_id}.iss_R2.fastq.gz`

The simulated reads are tagged with metadata:

| Field | Value |
|-------|-------|
| `meta.id` | `{parent_sample_id}_insilico_iss` |
| `meta.parent_id` | `{parent_sample_id}` |
| `meta.insilico` | `true` |
| `meta.control` | `false` |
| `meta.platform` | `ILLUMINA` |
| `meta.single_end` | `false` |

### Step 2b: NanoSim (ONT Simulation)

Nanosim requires a preparation step that converts the accessions abundance file into organism level values.

#### Preparation (prepare_nanosim_inputs.py)

Converts the ISS style abundance file into NanoSim's metagenome format:

**genome_list.tsv** — Maps organism names to individual FASTA files:
```
Escherichia_coli    genomes/Escherichia_coli.fasta
Zika_virus          genomes/Zika_virus.fasta
```

**size_file.tsv** — Header line with total read count, then organism level abundances as percentages:
```
Size    2500
Escherichia_coli    45.62
Zika_virus          54.38
```

The ONT read count is computed as `sim_nreads / sim_ont_divisor` (default: 100000 / 40 = 2500 reads).

#### Simulation

**Container:** `biocontainers/nanosim:3.2.3`

**Command:**
```bash
simulator.py metagenome \
    -gl genome_list.tsv \
    -a size_file.tsv \
    -c {training_dir}/{model_base} \
    -o {sample_id}.nanosim \
    --fastq \
    --perfect \
    -t {cpus}
```

**Output:** `{sample_id}.nanosim_aligned_reads.fastq.gz`

The `--perfect` flag generates reads without sequencing errors, useful for baseline detection benchmarking.

Tagged metadata:

| Field | Value |
|-------|-------|
| `meta.id` | `{parent_sample_id}_insilico_nanosim` |
| `meta.parent_id` | `{parent_sample_id}` |
| `meta.insilico` | `true` |
| `meta.control` | `false` |
| `meta.platform` | `OXFORD` |
| `meta.single_end` | `true` |

### Step 3: Integration into Main Alignment Pipeline

The pipeline creates the necessary supporting data for each insilico sample by cloning the parent sample's reference prep data:

- **Reference FASTA files** — Same as sample
- **Mapping file** — Same as sample (merged_taxid.tsv)
- **Kraken2 report** — Placeholder NO_FILE (insilico samples skip classification)
- **Assembly analysis** — Placeholder NO_FILE2

This means insilico reads align against the **full reference set** (not just the organisms used to generate them), which is needed for detecting false positives.

### Step 4: Report Stage — 3-Way Branch

In the reporting step, all results from alignment are split into branches:

```groovy
alignments.branch {
    control: it[0].control == true      // Lab controls (negative/positive)
    insilico: it[0].insilico == true    // In-silico simulated samples
    noncontrol: true                     // Real samples
}
```

**Branch processing:**

1. **Control samples** → `ALIGNMENT_PER_SAMPLE_CONTROLS` — Processed first, outputs collected as control JSONs
2. **Insilico samples** → `ALIGNMENT_PER_SAMPLE_INSILICO` — Processed independently, outputs collected as insilico JSONs keyed by `parent_id`
3. **Non-control samples** → `ALIGNMENT_PER_SAMPLE` — Receives both lab control JSONs and insilico JSONs as inputs

When a parent sample has multiple insilico children (one ISS, one NanoSim), all their JSONs are collected into a list and passed together via `--insilico_controls`.

### Step 5: Comparison Annotation (match_paths.py)

For each non control sample, `match_paths.py` receives the insilico JSON(s) and performs:

#### Combined Annotation

All insilico JSONs are loaded into a single point. For each organism in the sample, the script does:

- **`insilico_tass`** — The insilico sample's TASS score for this organism
- **`insilico_reads`** — The insilico sample's read count for this organism
- **`tass_fold_over_insilico`** — Ratio: sample TASS / insilico TASS
- **`reads_fold_over_insilico`** — Ratio: sample reads / insilico reads

#### Per-Simulator-Type Annotation

The script classifies insilico files by simulator method based on filename patterns (`_insilico_iss` vs `_insilico_nanosim`). For each type, it creates separate comparison data stored under:

- `insilico_comparison_iss`
- `insilico_comparison_nanosim`

#### Missing Organism Detection

2 categories are identified:

**Missing from sample (False Negatives):** Organisms present in the insilico simulation but absent from the real sample. Detected at species/subkey level via `find_missing_positive_controls()`. Each missing organism is enriched with its microbial category (Primary, Opportunistic, Potential, Commensal) from the pathogens list.

**Sample-only (False Positives):** Organisms detected in the real sample but are absent from the insilico simulation

### Step 6: Report Rendering (create_report.py)

#### Metrics Table

TP/FP/FN/TN per microbial category:

| Classification | Condition |
|---------------|-----------|
| **True Positive (TP)** | Organism is in the simulation AND its TASS score >= confidence threshold |
| **False Positive (FP)** | Organism is NOT in the simulation AND its TASS score >= confidence threshold |
| **False Negative (FN)** | Organism is in the simulation BUT its TASS score < threshold, OR it is entirely missing from the sample |
| **True Negative (TN)** | Organism is NOT in the simulation AND its TASS score < threshold |

Derived metrics per category and total:

- **Precision** = TP / (TP + FP)
- **Recall** = TP / (TP + FN)
- **F1** = 2 * Precision * Recall / (Precision + Recall)
- **Accuracy** = (TP + TN) / (TP + FP + FN + TN)

When multiple simulator types are present, one table is made for each type with headers like "InSilicoSeq (Illumina) Metrics" and "NanoSim (ONT) Metrics".

#### Missing Organisms Detail Table

Below each metrics table, a red-themed detail table lists each organism that was present in the simulation but missing from the sample:

| Column | Description |
|--------|-------------|
| Organism | Species/strain name |
| Taxid | Taxonomic ID |
| Category | Microbial category (Primary, Opportunistic, etc.) |
| InSilico TASS | TASS score in the simulation (0-100) |
| InSilico Reads | Read count in the simulation |
| Status | "Missing from sample (FN)" |

#### Sample-Only Organisms

After the missing organisms table, a text section lists organisms detected in the real sample but absent from the simulation. These count as FP in the metrics above.

#### Ctrl Column in the Main Organism Table

The bar plot in the Ctrl column shows change comparisons with symbols indicating the control type:

| Symbol | Color | Meaning |
|--------|-------|---------|
| **−** (minus) | Dark gray or red | Change vs. negative control |
| **+** (plus) | Green | Change vs. positive control |
| **∞** (infinity loop) | Teal | Change vs. in silico control |

The fold values show `#x` for TASS fold-change and `#x rd` for read count fold-change. When an organism is in the sample but not in the simulation, the teal-ish row displays `∞ not in sim`.

## Example Usage

### ISS Only (Illumina)

```bash
nextflow run jhuapl-bio/taxtriage \
    --generate_iss \
    --sim_nreads 100000 \
    --iss_model miseq \
    --input samplesheet.csv \
    --db /path/to/kraken2_db \
    --outdir results
```

### NanoSim Only (ONT)

```bash
nextflow run jhuapl-bio/taxtriage \
    --generate_nanosim \
    --nanosim_training /path/to/training_model \
    --sim_nreads 100000 \
    --sim_ont_divisor 40 \
    --input samplesheet.csv \
    --db /path/to/kraken2_db \
    --outdir results
```

### Both Simulators

```bash
nextflow run jhuapl-bio/taxtriage \
    --generate_iss \
    --generate_nanosim \
    --nanosim_training /path/to/training_model \
    --sim_nreads 100000 \
    --iss_model miseq \
    --sim_ont_divisor 40 \
    --input samplesheet.csv \
    --db /path/to/kraken2_db \
    --outdir results
```

When both are enabled, each real sample produces two insilico children (e.g., `sample1_insilico_iss` and `sample1_insilico_nanosim`), and the report renders separate metrics tables for each.

### Custom Abundance Profile

```bash
nextflow run jhuapl-bio/taxtriage \
    --generate_iss \
    --sim_abundance /path/to/custom_abundance.tsv \
    --input samplesheet.csv \
    --db /path/to/kraken2_db \
    --outdir results
```

The custom abundance file should be a two-column TSV: `taxid<TAB>abundance`.

## Interpreting Results

### High Precision, Low Recall

The pipeline identifies what it detects, but misses some organisms that were simulated. Check the "Missing In-Silico Organisms" table for which organisms were not recovered. Likely causes: bad coverage, too short sequences to compare against, or the organism(s) are below the TASS threshold.

### Low Precision, High Recall

The pipeline detects most simulated organisms but also reports organisms that were not in the simulation. Check the "Sample-Only Organisms" section. These false positives may indicate cross-mapping between similar reference sequences, contamination in the reference database, or organisms that align well but weren't part of the simulation input.

### Per-Type Differences

If ISS and NanoSim produce different metrics, this reflects the impact of read length and error profile on detection accuracy. Long ONT reads may resolve repeat regions better but have lower total coverage; short Illumina reads provide higher coverage but may cross-map between similar organisms.

## Output JSON Fields

The `match_paths.py` output JSON for non-control samples includes:

```json
{
  "metadata": {
    "insilico_controls_used": ["sample1_insilico_iss.paths.json", "sample1_insilico_nanosim.paths.json"],
    "insilico_simulator_types": ["iss", "nanosim"],
    "missing_insilico_controls": [
      {
        "name": "Zika virus",
        "id": "64320",
        "microbial_category": "Primary",
        "pos_tass_score": 0.87,
        "pos_numreads": 5000,
        "control_source": "insilico",
        "missing_control": true
      }
    ],
    "missing_insilico_by_type": {
      "iss": [...],
      "nanosim": [...]
    }
  },
  "organisms": [
    {
      "toplevelkey": "64320",
      "insilico_comparison": {
        "insilico_tass": 0.92,
        "insilico_reads": 5000,
        "tass_fold_over_insilico": 0.0,
        "reads_fold_over_insilico": 0.0,
        "missing_from_insilico": false
      },
      "insilico_comparison_iss": { ... },
      "insilico_comparison_nanosim": { ... }
    }
  ]
}
```
