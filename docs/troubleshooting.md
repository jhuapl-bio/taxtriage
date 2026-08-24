# Troubleshooting & FAQ

This page covers common questions, error messages, and debugging strategies for TaxTriage.

---

## General Tips

### Reading the Nextflow Log

Red text and `ERROR` messages in the Nextflow log do **not** always mean the pipeline failed — many steps emit warnings as part of normal operation. Look for these indicators of an actual failure:

```
ERROR ~ Error executing process > 'NFCORE_TAXTRIAGE:TAXTRIAGE:<MODULE_NAME>'
```

The key line to read is `Command error:`, which contains the actual error message:

```
Command error:
  [CRITICAL] At least the first FASTQ file is required. On line 2.
```

### Why Does My Pipeline Run Slowly?

Several steps are inherently slow:

- **MultiQC** (final step) — especially with many samples
- **NanoPlot** — can be very slow for ONT data; disable with `--skip_plots` if not needed
- **Alignment** — more top hits = more reference FASTAs = longer alignment
- **De novo / reference assembly** — avoid unless needed (`--use_denovo`, `--reference_assembly`)
- **Kraken2 loading** — use `--low_memory` if RAM is limited, but this is much slower

---

## Input Issues

### What is a Samplesheet?

A samplesheet is a `.csv` file with one row per sample, specifying the FASTQ file paths, platform, and optional metadata. See the full [Samplesheet](samplesheet.md) reference page.

### My pipeline fails at SAMPLESHEET_CHECK

```
Command error:
  [CRITICAL] At least the first FASTQ file is required. On line 2.
```

The `fastq_1` column in your samplesheet is missing a value on that row. Each row must have a valid path to a FASTQ file (or a directory if using `--demux`). Verify your samplesheet and re-run with `-resume`.

### What is the difference between `--remove_taxids` and `--genome`/`--remove_reference_file`?

- `--remove_taxids "9606"` — removes reads **classified by Kraken2** as human (post-classification). Does not catch reads that Kraken2 missed.
- `--remove_reference_file` / `--genome` — aligns all reads against a host reference and removes those that map (pre-classification). More comprehensive but requires a host FASTA or iGenomes key.

> Neither option guarantees complete host removal. Review data privacy requirements before sharing processed data.

---

## Common Module Errors

### FASTP

**Symptom:** FASTP fails with exit status `255`

```
Command error:
  [WARNING] No reads passed filter...
```

**Causes and fixes:**

- Read quality is lower than `--minq` threshold — lower `--minq` or use `--skip_fastp`
- FASTQ file path is wrong or the file is empty — check your samplesheet paths

### KRAKEN2_KRAKEN2

**Memory error:**

```
Command error:
  Kraken2 DB could not be loaded into memory
```

- The Kraken2 database is larger than available RAM.
- Solution A: `--low_memory` — reads the DB from disk (much slower)
- Solution B: `--max_memory 13GB` — explicitly set memory limit (Nextflow 24.x–25.x; on 26+ set `process.resourceLimits` in a config instead — see [CLI Parameters](cli-parameters.md#workflow-control-and-execution))

**Invalid or corrupt database:**

Verify the database directory contains exactly these three files:

- `hash.k2d`
- `opts.k2d`
- `taxo.k2d`

If any are missing or the directory is empty, re-download the database.

### Bowtie2 / Minimap2 (Alignment)

**Empty reference FASTA:**

```
download/<sample>.dwnld.references.fasta  (empty or missing)
```

- If using Kraken2 (default), check internet connectivity — NCBI FASTA downloads use `curl`
- Check `top/<sample>.top_report.tsv` — it must have at least one line of data
- If the file exists but is empty: you may be behind a firewall; use `--reference_fasta` with a local FASTA instead

**Minimap2: `no SQ lines in header` / `Parse error`:**

```
[WARNING] For a multi-part index, no @SQ lines will be outputted. Please use --split-prefix.
[E::sam_parse1] no SQ lines present in the header
```

This usually means the reference FASTA was loaded incorrectly. The pipeline automatically retries with increased memory up to 3 times. If it persists:

- Check that your reference FASTA file is not corrupt
- Reduce the number of top hits with `--top_hits_count` or `--top_per_taxa`
- Try `--no_split_prefix` to disable split-prefix mode

**Out of memory during alignment (e.g. exit code `143`):**

Minimap2 loads the reference index into RAM. The pipeline now backs off automatically: on each retry it shrinks the `-I` index batch (8G → 4G → 2G → 1G) and `-K` query batch (500M → 200M → 100M → 50M) while stepping the memory request up (24 → 40 → 56 → 72 GB), retrying up to 3 times. Most OOM kills resolve on a later attempt with no action needed.

If it still fails:

- Set `--mmap2_I` / `--mmap2_K` to smaller values to force low-RAM mode from attempt 1 (e.g. `--mmap2_I 1G --mmap2_K 50M`)
- Reduce `--top_hits_count` or `--top_per_taxa` to limit the number of references
- Increase Docker/Singularity memory limits
- Split the run into smaller batches of samples

### SAMTOOLS_HIST_COVERAGE

Make sure a `.bam` file exists in the `alignment/` directory for each sample. An empty or missing BAM indicates alignment failed upstream.

### SpAdes / FLYE (De Novo Assembly)

De novo assembly is memory-intensive. If it fails:

1. Remove `--reference_fasta` if a large FASTA was provided
2. Lower `--top_hits_count` or `--top_per_taxa` to reduce input reads
3. Increase the memory limit with `--max_memory`

### SAMTOOLS_SORT

Usually CPU-intensive, but can fail if RAM is insufficient. Ensure you have sufficient memory set with `--max_memory`, and that the input BAM is not empty.

### SAMTOOLS_VIEW

Typically fails only due to insufficient disk space or timeout. Check available storage on your system.

---

## Assemblies Not Downloading

**Symptom:** The pipeline stalls or fails when trying to pull assemblies from NCBI.

**Possible causes:**

- Firewall or proxy blocking outbound `curl` requests
- NCBI FTP temporarily unavailable
- Missing SSL/TLS certificates in your environment

**Solutions:**

- Use `--assembly /path/to/local_assembly_summary.txt` with a downloaded copy of the NCBI assembly summary
- Use `--reference_fasta` with manually downloaded FASTA files
- Work with your IT team to allow outbound HTTPS requests to NCBI FTP

---

## Uncommitted Changes Warning

If you see:

```
Error: Uncommitted changes in the repository
```

Run:

```bash
nextflow drop -f https://github.com/jhuapl-bio/taxtriage
```

Then re-run the pipeline. This only applies when running from the remote URL with `-latest` or `-r main`.

---

## NanoPlot Failing on Empty Data

If NanoPlot fails with an empty directory or empty FASTQ, that sample had no reads pass QC. Check the fastp log for that sample and consider lowering `--minq` or disabling fastp with `--skip_fastp`.

---

## A Species Scores 80+ but All Its Strains Score Below 10

This gap is itself the diagnosis: strain-level scores are computed directly from covered bases ÷ genome length, while species and genus levels apply roll-up coverage overrides. When the two disagree that sharply, the override is the suspect, not the alignment.

**Confirm it** by comparing two columns of the same species row in the interactive report:

| Column      | Suspicious pattern    |
| ----------- | --------------------- |
| `Coverage`  | high — 80–100         |
| `Breadth %` | near zero — 0.00–0.05 |

Those two describe the same organism and should not disagree. `Breadth %` is honest (covered bases ÷ full genome length); `Coverage` is the roll-up override. If `Coverage` is high while `Breadth %` is ~0, the override has latched onto a small accession.

**Why it happens.** The species-level coverage override takes the maximum breadth fraction across member accessions, and an accession is one BAM reference — one _contig_. For chromosome-level assemblies one accession is the genome and this is correct. For scaffold-level draft assemblies (many eukaryotic references, and any `GCA_` assembly still in thousands of pieces) a single read covering one ~900 bp scaffold yields a 0.94 breadth fraction that is then applied to the whole genome. That value drives the breadth term _and_ opens the minhash gate — the two heaviest default weights — so TASS lands near 90 for an organism that is absent.

**What to do.** Current defaults (`--rep_breadth_min_frac 0.01`, `--rep_breadth_min_len 50000`) already require an accession to be a meaningful fraction of the genome before it can set the group maximum. If you still see the pattern, the offending accession is clearing those bars — raise either value and re-run the scoring step. Check the `[LCA] representative-breadth eligibility:` line in the `match_paths.py` log to see how many accessions were excluded.

Legitimate low-breadth detections do exist — a genuine sterile-site pathogen at very low titre. Distinguish them by checking whether reads are spread across _many_ large accessions (real) or concentrated on a few tiny ones (artefact), and by whether the same organism scores consistently across replicate samples.

See [TASS Scoring 11.3](tass-scoring.md#113-size-eligibility-for-the-representative-coverage-maximum) for the full mechanism.

---

## A Sample Is Missing From the Interactive Report

The report can hide a sample for three different reasons, and they undo in
different places. Check them in this order:

| Cause                                                   | How you can tell                                                                                                            | Fix                                                                                                                                        |
| ------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| A **sample QC rule** with a _flag & hide_ action        | The **Sample QC / Flags** panel in the right sidebar reports a flagged / hidden count; open the dialog to see which samples | Open **Filters & flags…**, then either untick **Hide all flagged samples**, change the rule's action back to _flag it_, or delete the rule |
| Hidden **by hand**, via the eye icon on its sidebar row | Its row in the Samples list is dimmed with a crossed-out eye                                                                | Click the eye icon again, or use **Hide All / Show All**                                                                                   |
| The **sidebar search** auto-hide                        | A query is present in the Samples search box and **Hide filtered-out samples** is ticked                                    | Clear the search box, or untick that option                                                                                                |

None of these delete anything — the sample's data is still in the file, and the
Metadata & Mapping table keeps listing it as part of the run inventory. Only the
red **trash icon** on a sidebar row actually removes a sample from the loaded
dataset, and even that leaves the file on disk untouched.

If the sample is missing from the **Samples list itself**, it never reached the
report: check that its `<sample>.paths.json` was produced by the run, and see
[Interactive Report](interactive-report.md) for what a sample with no
alignments looks like.

If the flags were not something you configured, they came from the pipeline —
see [Report Sample-QC Flags](cli-parameters.md#report-sample-qc-flags) for the
`--report_flag_*` parameters that seed them.

---

## Debugging a Specific Step

When a module fails, you can inspect and rerun its exact command:

1. Find the working directory hash from the Nextflow log output:

   ```
   [28/5412a9] process > NFCORE_TAXTRIAGE:TAXTRIAGE:KRAKEN2_KRAKEN2
   ```

2. Navigate to the directory:

   ```bash
   cd work/28/
   # Tab-complete the full path:
   cd work/28/5412a9<TAB>
   ```

3. Inspect the command and its output:

   ```bash
   cat .command.sh   # The exact command that was run
   cat .command.out  # Stdout
   cat .command.err  # Stderr
   ```

4. Rerun the command:
   ```bash
   bash .command.sh
   ```

> Non-global scripts (Python/bash tools) are usually at `../../../bin/<script.py>`. Prepend that path if running from the work directory directly.

> ⚠️ Some commands require input files from the original run. Ensure those files are still present before rerunning.

---

## AWS / Nextflow Credential Warnings

```
WARN: Nextflow will attempt to connect to AWS as an anonymous user
```

This is a warning, not a failure, for most runs. However, if you use private S3 buckets, configure your AWS credentials:

```bash
export AWS_ACCESS_KEY_ID=...
export AWS_SECRET_ACCESS_KEY=...
export AWS_DEFAULT_REGION=...
```

Or use an IAM role if running on EC2.

---

## Getting Help

- **Slack:** [#taxtriage on nf-core Slack](https://nfcore.slack.com/channels/taxtriage) ([join here](https://nf-co.re/join/slack))
- **GitHub Issues:** [github.com/jhuapl-bio/taxtriage/issues](https://github.com/jhuapl-bio/taxtriage/issues)
- **Seqera/AWS:** Contact brian.merritt@jhuapl.edu for access to the JHU/APL Seqera instance
