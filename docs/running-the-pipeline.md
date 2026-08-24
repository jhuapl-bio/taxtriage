# Running the Pipeline

This page covers the core execution patterns for TaxTriage — from the standard run command to profiles, offline mode, cloud execution, and background operation.

---

## Standard Run

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  --input samplesheet.csv \
  --outdir ./results \
  -profile docker \
  -resume
```

> ⚠️ For laptops and workstations, start with `-profile local,docker` to apply conservative resource defaults before tuning parameters manually.

Nextflow creates the following in your working directory:

```
work/               # Intermediate working files per process
<OUTDIR>/           # Final results (set via --outdir)
.nextflow_log       # Execution log
```

---

## Running Straight from an SRA / ENA Accession

You do not need to download reads yourself. Wherever TaxTriage accepts a FASTQ path it also accepts a public archive accession, and it will resolve, fetch and feed the reads in as if you had supplied local files.

For a single sample, pass the accession to `--fastq_1` instead of a path:

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  --fastq_1 ERR1160846 \
  --type stool \
  --sample ont_survey \
  --outdir ./results \
  -profile local,docker \
  -resume
```

That is the same run as this row in a samplesheet, so use whichever fits your workflow:

```
sample,platform,fastq_1,fastq_2,sequencing_summary,trim,type
ont_survey,,ERR1160846,,,FALSE,stool
```

Note what you do **not** have to supply. `--fastq_2` is unnecessary because paired-end layout is read from the archive's actual file listing rather than guessed, and `--platform` is unnecessary because the instrument platform comes from the archive too, which in turn selects the right minimap2 preset (`ILLUMINA`, `OXFORD` for Oxford Nanopore, or `PACBIO`). Fill either one in only if you want to override what the archive reports. If you omit `--sample`, the accession itself becomes the sample name.

### One accession, many samples

Run accessions (`SRR`, `ERR`, `DRR`) give you one sample. Every other accession level expands into one pipeline sample per child run, which is the main thing to be aware of before launching:

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  --fastq_1 PRJNA681875 \
  --type stool \
  --sample whole_project \
  --outdir ./results \
  -profile local,docker \
  -resume
```

This does not produce one sample called `whole_project`. It produces one sample per run in the project, each named `<sample>_<run_accession>` (`whole_project_SRR13191701`, `whole_project_SRR13191702`, and so on) so the names stay unique and you can trace any row back to its source run. Everything else you passed, such as `--type stool`, is inherited by all of them. A project can hold hundreds of runs, so check how many you are about to pull before committing the compute.

The same expansion applies to experiment (`SRX`) and sample (`SRS`, `SAMN`) accessions. See [Samplesheet → Accepted accession types](samplesheet.md#accepted-accession-types) for the full table.

### Downloads are cached

Reads land in `<outdir>/sra_downloads/<run_accession>/` and are skipped on any later run where they already exist, on a fresh run as well as with `-resume`. Point `--sra_cache_dir` at shared storage so an accession is only ever downloaded once across projects or users:

```bash
  --sra_cache_dir /shared/sra
```

For the retrieval mechanics (ENA fast path, sra-tools fallback) see [Samplesheet → SRA / ENA Accessions](samplesheet.md#sra--ena-accessions); for `--sra_force_sratools`, `--sra_max_size` and `--ncbi_api_key` see [CLI Parameters → SRA / ENA Download](cli-parameters.md#sra--ena-download).

---

## Starting from an Alignment You Already Have

If your reads are already aligned, hand TaxTriage the BAM instead of FASTQ and it skips straight to coverage statistics and scoring:

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  --bam my_sample.bam \
  --reference_fasta refs.fasta \
  --type stool \
  --sample my_sample \
  --outdir ./results \
  -profile local,docker \
  -resume
```

`.cram` and `.sam` work too, and the file does not need to be sorted or indexed first. The samplesheet equivalent is a `bam` column, which can be mixed with ordinary `fastq_1` rows in the same sheet.

Pre-aligned samples skip read QC, trimming, host removal, classification and reference download, so no Kraken2 database is needed for a BAM-only run. De novo assembly, VF/AMR annotation and novelty detection are unavailable for them, since all three need raw reads.

Pass `--reference_fasta` pointing at whatever the alignment was made against whenever you can. A BAM header carries reference names and lengths but no bases, and without a reference TaxTriage has to reconstruct one from the alignment, which skews conflict-based read removal. See [Pre-Aligned Input](pre-aligned-input.md) for the full picture.

---

## Profiles

Profiles control resource allocation and container runtime. Multiple profiles can be chained with commas (order matters — later profiles override earlier ones):

```bash
-profile test,docker
```

### Container Profiles

| Profile       | Description                        |
| ------------- | ---------------------------------- |
| `docker`      | Use Docker for all containers      |
| `singularity` | Use Singularity for all containers |

### Execution Profiles

| Profile      | Description                                              |
| ------------ | -------------------------------------------------------- |
| `local`      | Reduced resource limits for laptops/workstations         |
| `test`       | Minimal test dataset (pulls from GitHub)                 |
| `test_viral` | Test dataset with the viral Kraken2 database             |
| `mce`        | Uses the pathogen FASTA sheet for alignment (no Kraken2) |

---

## Key Execution Flags

| Flag      | Description                                                        |
| --------- | ------------------------------------------------------------------ |
| `-resume` | Resume from the last successful step (uses Nextflow's cache)       |
| `-latest` | Pull the latest commit from the specified branch                   |
| `-r main` | Use the `main` branch (or `stable`, or a version tag like `1.3.1`) |
| `-bg`     | Run Nextflow in the background, detached from the terminal         |

---

## Offline / No Internet Mode

For air-gapped environments, provide all remote resources as local files:

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  --input examples/Samplesheet.csv \
  --db "./k2_viral" \
  -r main -latest \
  --skip_kraken2 \
  --outdir tmp \
  --reference_fasta ./refer.fasta \
  --assembly examples/assembly_summary_refseq.txt \
  -profile local,docker \
  -resume --demux
```

Three files are required locally:

1. A Kraken2 database directory (or use `--skip_kraken2` with a local FASTA)
2. A reference FASTA file (`--reference_fasta`)
3. An NCBI assembly summary text file (`--assembly`)

> ⚠️ Nextflow may still attempt anonymous connections to AWS S3 by default. Configure AWS credentials or use `--assembly` to bypass NCBI lookups.

---

## Using a Local Database

To use a Kraken2 database already on your filesystem (e.g., a decompressed `k2_viral` directory):

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  --input examples/Samplesheet.csv \
  --db "./k2_viral" \
  -r main -latest \
  --outdir output_viral_local \
  -profile local,docker \
  -resume
```

For a full list of downloadable databases, see [CLI Parameters → Databases](cli-parameters.md#taxonomic-classification-and-databases).

---

## Resource Limits

Override memory and CPU defaults:

```bash
--max_memory 10GB    # Max RAM for the pipeline (primarily used by Kraken2)
--max_cpus 3         # Max CPU cores
--low_memory         # Read Kraken2 DB from disk instead of loading into RAM (slower)
```

---

## Removing Host Reads

To remove human reads before classification:

```bash
--remove_taxids "9606"          # Remove human reads via Kraken2 taxid
--remove_reference_file host.fasta  # Remove reads aligned to a host FASTA
--genome GRCh37                 # Auto-download an iGenomes host reference
```

---

## Using a Backup Reference FASTA

Ensure specific organisms are always aligned, regardless of Kraken2 results:

**Option A — Local FASTA file** (header must follow `>accession description` format):

```fasta
>NC_003663.2 Cowpox virus, complete genome
```

```bash
--reference_fasta ./my_references.fasta
```

**Option B — NCBI taxids** (requires internet):

```bash
--organisms "10243 2331"
```

---

## Reproducibility

Pin a specific pipeline version with `-r`:

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  -r 1.3.1 \
  --input samplesheet.csv \
  --outdir results \
  -profile docker
```

Version numbers are logged in all reports. Find releases at the [GitHub releases page](https://github.com/jhuapl-bio/taxtriage/releases).

---

## Running in the Background

Use Nextflow's built-in `-bg` flag to detach the process from your terminal:

```bash
nextflow run ... -bg
```

Alternatively, use `screen` or `tmux`:

```bash
screen -S taxtriage
# Run your nextflow command inside the screen session
# Detach: Ctrl+A, D
```

---

## Custom Configuration

### Overriding Container Versions

To use a different tool version for a specific process, create a custom config file:

```groovy
// custom.config
process {
  withName: 'KRAKEN2_KRAKEN2' {
    container = 'docker.io/biocontainers/kraken2:2.1.3--pl5321h9f5acd7_0'
  }
}
```

Pass it with:

```bash
-c custom.config
```

### Using Institutional nf-core Configs

Many HPC clusters have pre-configured profiles at [nf-core/configs](https://github.com/nf-core/configs). These are automatically loaded at runtime. Use:

```bash
-profile <institute>
```

---

## Next Steps

- [CLI Parameters](cli-parameters.md) — complete parameter reference
- [Cloud & Seqera](cloud-and-seqera.md) — AWS/Seqera Tower execution
- [Output](output.md) — understanding your results
