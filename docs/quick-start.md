# Quick Start

This guide walks you through a test run of TaxTriage to verify your installation and get familiar with the pipeline output.

---

## Prerequisites

Make sure you have installed:

- [Nextflow](installation.md#1-install-nextflow) (`>= 24.04.0`)
- [Docker](installation.md#a-docker-recommended) or [Singularity](installation.md#b-singularity-hpc)

---

## Test Run (Recommended First Step)

This command pulls test data from GitHub and runs the full pipeline. It takes approximately **10–15 minutes**.

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  -r main -latest \
  -profile test,docker \
  -resume
```

Replace `docker` with `singularity` if running on an HPC:

> ```bash
> -profile test,singularity
> ```

---

## Running on Local FASTQ Files

Once you've confirmed the test run works, run the pipeline on your own data:

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  --input samplesheet.csv \
  --outdir ./results \
  -profile local,docker \
  -resume
```

> **Tip:** Always use `-profile local,docker` for laptops and workstations before tuning advanced parameters.

---

## Running on a Public SRA Accession

No local reads needed — put a bare accession in `fastq_1` and TaxTriage downloads
it for you. `SRR13191702` is a small nasal-swab run, which makes it a good first
real dataset.

Create `srr_samplesheet.csv`:

```
sample,platform,fastq_1,fastq_2,sequencing_summary,trim,type
covid_run,,SRR13191702,,,FALSE,nasal
```

Then run it:

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  --input srr_samplesheet.csv \
  --outdir ./results_srr \
  --db "viral" --download_db \
  -profile local,docker \
  -resume
```

Leave `platform` blank and the instrument platform is taken from the archive.
Paired-end is detected from the files themselves, so nothing needs declaring.
Downloads land in `<outdir>/sra_downloads/<run_accession>/` and are skipped on
later runs — see [Samplesheet](samplesheet.md#sra--ena-accessions) for accession
types and caching.

---

## Starting from a BAM You Already Have

If reads are already aligned, skip straight to the scoring stages by supplying
the alignment and the reference it was aligned against:

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  --bam my_sample.bam \
  --reference_fasta refs.fasta \
  --sample my_sample \
  --type stool \
  --outdir ./results_bam \
  -profile local,docker \
  -resume
```

For more than one alignment, use a samplesheet with a `bam` column instead of
`fastq_1`:

```
sample,platform,bam,type
prealigned_a,ILLUMINA,/data/prealigned_a.bam,stool
prealigned_b,OXFORD,/data/prealigned_b.cram,nasal
```

CRAM works the same way. A sheet needs **either** `fastq_1` or `bam`; if both
columns exist, a row that fills in `bam` is treated as pre-aligned. See
[Pre-Aligned Input](pre-aligned-input.md) for mixing both styles and for what
the pipeline skips.

---

## Quick Example with Auto-Downloaded Viral Database

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  --outdir tmp_viral \
  --input examples/Samplesheet.csv \
  -r main -latest \
  --db "viral" --download_db \
  -profile local,docker \
  -resume
```

---

## Running from a Local Clone

If you'd like to make local changes or run without re-pulling from GitHub every time:

```bash
# Clone the repository
git clone https://github.com/jhuapl-bio/taxtriage.git
cd taxtriage

# Run from local main.nf
nextflow run ./main.nf -profile test,docker -resume
```

> ⚠️ If you get an error about uncommitted changes when pulling from the remote URL, run:
>
> ```bash
> nextflow drop -f https://github.com/jhuapl-bio/taxtriage
> ```

---

## Offline / Air-Gapped Mode

If no internet is available, provide local copies of the 4 required files/folders:

| File                       | Parameter                                  | Source                                                                                           |
| -------------------------- | ------------------------------------------ | ------------------------------------------------------------------------------------------------ |
| Kraken2 database directory | `--db ./k2_viral`                          | [K2 Indices](https://benlangmead.github.io/aws-indexes/k2)                                       |
| Reference FASTA            | `--reference_fasta ./refer.fasta`          |                                                                                                  |
| NCBI assembly summary      | `--assembly ./assembly_summary_refseq.txt` | [NCBI Assembly Summary](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt) |
| Taxdump                    | `--taxdump ./taxdump`                      | [NCBI Taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)                         |

Example command:

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  --input examples/Samplesheet.csv \
  --db "k2_viral" -r main -latest \
  --skip_kraken2 \
  --outdir tmp \
  --reference_fasta ./refer.fasta \
  --assembly examples/assembly_summary_refseq.txt \
  -profile local,docker \
  -resume --demux
```

> ⚠️ Using `--skip_kraken2` with a local FASTA **only performs alignment** — metagenomics classification is skipped.

---

## Available Databases for `--download_db`

| Name          | Size   | Notes                          |
| ------------- | ------ | ------------------------------ |
| `viral`       | 553 MB | Good for quick viral screening |
| `standard8`   | 7.5 GB | Standard 8 GB Kraken2          |
| `flukraken2`  | 180 MB | Influenza-focused              |
| `minikraken2` | 7.5 GB | Standard (alternate mirror)    |
| `test`        | 112 MB | Minimal test database          |

For the full database list, see [CLI Parameters → Databases](cli-parameters.md#taxonomic-classification-and-databases).

---

## Next Steps

- Prepare your input data: [Samplesheet](samplesheet.md)
- Explore all run options: [Running the Pipeline](running-the-pipeline.md)
- See what outputs to expect: [Output](output.md)
