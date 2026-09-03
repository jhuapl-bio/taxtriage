# CLI Parameters

Complete reference for all TaxTriage parameters, organized by category. Parameters prefixed with `--` are pipeline parameters (double hyphen); Nextflow-native flags use a single hyphen (e.g., `-resume`, `-profile`).

---

## Input and Sample Information

| Parameter                             | Description                                                                                                                                                                                                                                                                              |
| ------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--input <samplesheet.csv>`           | CSV samplesheet. See [Samplesheet](samplesheet.md) for column reference.                                                                                                                                                                                                                 |
| `--outdir <path>`                     | Output directory (will be created if it doesn't exist).                                                                                                                                                                                                                                  |
| `--fastq_1 <path\|accession>`         | Path to the primary FASTQ file, **or** an SRA/ENA accession (`SRR`/`ERR`/`DRR`, `SRX`, `SRS`/`SAMN`, `SRP`/`PRJNA`) to download automatically. Required if `--input` is not provided (single-sample mode). See [Samplesheet → SRA / ENA Accessions](samplesheet.md#sra--ena-accessions). |
| `--fastq_2 <path>`                    | Path to the second FASTQ file for paired-end Illumina. Optional. Ignored when `fastq_1` is an accession — paired-end layout is read from the archive.                                                                                                                                    |
| `--platform [ILLUMINA,OXFORD,PACBIO]` | Sequencing platform. Defaults to `ILLUMINA` if not specified.                                                                                                                                                                                                                            |
| `--type <str>`                        | Sample body site (e.g., `blood`, `stool`). See [Samplesheet → Sample Types](samplesheet.md#sample-types).                                                                                                                                                                                |
| `--sample <name>`                     | Sample name for single-sample runs. Defaults to `"sample"`.                                                                                                                                                                                                                              |
| `--sequencing_summary <path>`         | Path to a Nanopore sequencing summary file for pycoQC plotting.                                                                                                                                                                                                                          |
| `--demux`                             | Required when the samplesheet `fastq_1` column contains a **directory** rather than individual files.                                                                                                                                                                                    |
| `--meta <path>`                       | Path to a separate metadata CSV file. See [Samplesheet → Metadata](samplesheet.md#metadata-support).                                                                                                                                                                                     |

### SRA / ENA Download

Applies whenever `--fastq_1` or a samplesheet `fastq_1` cell holds an accession.

| Parameter                | Description                                                                                                                                                                                  |
| ------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--sra_cache_dir <path>` | Persistent cache for downloaded FASTQs. Defaults to `<outdir>/sra_downloads`. Point at shared storage to reuse downloads across runs and users — a cached run is skipped, not re-downloaded. |
| `--sra_force_sratools`   | Skip the ENA fast path and always fetch with `prefetch` + `fasterq-dump`. Slower; use when ENA's mirror is stale or unreachable.                                                             |
| `--sra_max_size <size>`  | Value passed to `prefetch --max-size` (e.g. `50G`). Defaults to `u` (unlimited).                                                                                                             |
| `--ncbi_api_key <key>`   | Optional NCBI eutils API key. Raises the rate limit when expanding accessions ENA doesn't know about.                                                                                        |

### Pre-Aligned Input (BAM / CRAM / SAM)

Start from reads you have already aligned. Pre-aligned samples skip QC, trimming, host removal, classification and reference download. See [Pre-Aligned Input](pre-aligned-input.md) for the full guide.

| Parameter                             | Description                                                                                                                                                                                                                                    |
| ------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--bam <path>`                        | Pre-aligned `.bam`, `.cram` or `.sam` file. Single-sample shortcut; the samplesheet equivalent is a `bam` column. Sorting and indexing are handled for you.                                                                                    |
| `--bam_minmapq <int>`                 | Optional MAPQ filter applied to pre-aligned input. Unset (the default) takes the alignment as given; `--minmapq` is not applied to it.                                                                                                         |
| `--cram_reference <path>`             | Reference FASTA used to decode CRAM input, when it differs from `--reference_fasta`.                                                                                                                                                           |
| `--bam_consensus <bool>`              | Recover reference sequence from the alignment with `samtools consensus` so sourmash/ANI have bases to work with. Auto-enabled when a BAM sample has no `--reference_fasta`; set `false` to run without the minhash/conflict component instead. |
| `--consensus_min_depth <int>`         | `samtools consensus -d`: minimum read depth to call a base. Default `1`.                                                                                                                                                                       |
| `--consensus_min_bases <int>`         | Drop consensus records with fewer than this many non-`N` bases. Default `500`.                                                                                                                                                                 |
| `--consensus_min_mapq <int>`          | `samtools consensus --min-MQ`. Unset by default.                                                                                                                                                                                               |
| `--consensus_mode <simple\|bayesian>` | `samtools consensus -m`. Default `simple`.                                                                                                                                                                                                     |

> Supply `--reference_fasta` (the reference the alignment was made against) whenever you have it. A consensus derived from the alignment overstates similarity between related organisms and makes conflict-based read removal more aggressive.

---

## Workflow Control and Execution

| Parameter                                  | Description                                                                                                                                             |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-resume`                                  | Resume from cached intermediate results.                                                                                                                |
| `-latest`                                  | Pull the latest commit from the specified branch at runtime.                                                                                            |
| `-r [main, stable, 1.x.x]`                 | Specify the GitHub branch or version tag.                                                                                                               |
| `-profile [local,test,docker,singularity]` | One or more execution profiles, comma-separated.                                                                                                        |
| `--max_memory <N>GB`                       | **Legacy (Nextflow 24.x–25.x).** Maximum RAM allocated to the pipeline. Kraken2 requires at least as much RAM as the database size. See the note below. |
| `--max_cpus <N>`                           | **Legacy (Nextflow 24.x–25.x).** Maximum number of CPU cores. See the note below.                                                                       |
| `--low_memory`                             | Load the Kraken2 database from disk via I/O rather than into RAM (much slower but lower memory).                                                        |

!!! warning "`--max_cpus` / `--max_memory` / `--max_time` are legacy knobs"

    These three parameters are the pre-`resourceLimits` way of capping resources, and they are **deprecated from Nextflow 26 onward** — treat them as applying to Nextflow **24.x and 25.x** only. Nextflow 24.04 introduced the native [`process.resourceLimits`](https://www.nextflow.io/docs/latest/reference/process.html#resourcelimits) directive, and the nf-core template dropped the `max_*` parameters in v3.0.0; the old `check_max()` helper they depended on is gone from this pipeline entirely, because function definitions are rejected by the strict config parser that became the default in Nextflow 26.04.

    On Nextflow 26 and later, set the ceiling natively instead — in a profile or a `-c` config:

    ```groovy
    process {
        resourceLimits = [ cpus: 16, memory: '64.GB', time: '24.h' ]
    }
    ```

    For now the parameters still take effect on every supported version, because `conf/base.config` feeds them straight into `resourceLimits`. They will not be maintained past the 25.x line, so new configs and new documentation should use `resourceLimits` directly.

> **Air-gapped / no-internet runs:** The pipeline pins `nf-tower@1.17.3` in `nextflow.config` so Nextflow uses the locally cached plugin without fetching remote metadata. If you hit a `Conversion = '4'` error on a machine that has never run the pipeline before, the plugin cache is empty — run once with internet access to populate it, or use `NXF_OFFLINE=false` temporarily. You can also pass Nextflow's native `-offline` flag (single hyphen) to fully suppress all network calls once plugins are cached:
>
> ```bash
> nextflow run . -offline [other params]
> ```

---

## Preprocessing and Read Filtering

| Parameter                  | Description                                                                                                |
| -------------------------- | ---------------------------------------------------------------------------------------------------------- |
| `--trim`                   | Enable adapter trimming (`"true"` or `"false"`). Defaults to `false`.                                      |
| `--skip_fastp`             | Skip quality filtering with fastp entirely.                                                                |
| `--minq <N>`               | Minimum read quality score. Default: `7` for ONT, `20` for Illumina. Disabled when `--skip_fastp` is set.  |
| `--subsample <N>`          | Subsample to N reads per sample before classification. Useful for very large datasets.                     |
| `--downsample`             | Use BBNorm (bbmap) to reduce redundant reads. Recommended for deep sequencing (e.g., >40 GB/sample).       |
| `--decompress_pre_megahit` | Decompress reads before MEGAHIT assembly. Required for some HPC environments with I/O errors at this step. |

---

## Taxonomic Classification and Databases

| Parameter                              | Description                                                                                                                                                                     |
| -------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--db <path_or_name>`                  | Kraken2 database. Provide a local path to a directory containing `.k2d` files, or a preset name with `--download_db`.                                                           |
| `--download_db`                        | Auto-download the preset database named in `--db` to `--outdir`.                                                                                                                |
| `--download_taxdump`                   | Download NCBI taxonomy `nodes.dmp` and `names.dmp`. Adds genus information to the PDF report.                                                                                   |
| `--taxdump <dir>`                      | Use a local NCBI taxonomy dump directory.                                                                                                                                       |
| `--k2_confidence <float>`              | Minimum Kraken2 confidence score to classify a read. See [Kraken2 confidence docs](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#confidence-scoring). |
| `--skip_kraken2`                       | Skip Kraken2 classification. Must provide `--reference_fasta` or `--organisms`.                                                                                                 |
| `--centrifuge`                         | Use Centrifuge instead of Kraken2 (specify Centrifuge DB with `--db`).                                                                                                          |
| `--organisms "<taxid1 taxid2>"`        | Force-include specific organisms by taxid or name, separated by spaces.                                                                                                         |
| `--organisms_file <file>`              | A file listing organisms (taxids or names) to include.                                                                                                                          |
| `--fuzzy`                              | Match organism names by fuzzy string matching rather than exact taxid.                                                                                                          |
| `--top_per_taxa "<taxid:amount:rank>"` | Minimum species count per parent taxon. Example: `"10239:20:S 2:20:S"` = ≥20 viral and ≥20 bacterial species. Rank codes: `G`, `S`, `P`, `F`, `O`, `C`.                         |
| `--top_hits_count <N>`                 | Minimum top hits per rank level. If `--top_per_taxa` produces a larger set, that takes priority.                                                                                |

### Supported Downloadable Databases

| Name                | Size   | URL                                                                                                       |
| ------------------- | ------ | --------------------------------------------------------------------------------------------------------- |
| `standard8`         | 7.5 GB | [Download](https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz)                   |
| `viral`             | 553 MB | [Download](https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20230605.tar.gz)                           |
| `flukraken2`        | 180 MB | [Download](https://media.githubusercontent.com/media/jhuapl-bio/mytax/master/databases/flukraken2.tar.gz) |
| `test`              | 112 MB | [Download](https://github.com/jhuapl-bio/datasets/raw/main/databases/kraken2/test_metagenome.tar.gz)      |
| `minikraken2`       | 7.5 GB | [Download](https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v2_8GB_201904.tgz)                      |
| `pluspf`            | 77 GB  | [Download](https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240112.tar.gz)                          |
| `pluspf8`           | 7.5 GB | [Download](https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08gb_20240112.tar.gz)                     |
| `eupath`            | 11 GB  | [Download](https://genome-idx.s3.amazonaws.com/kraken/k2_eupathdb48_20230407.tar.gz)                      |
| `centrifuge refseq` | 7 GB   | [Download](https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed%2Bh%2Bv.tar.gz)                    |

Additional Kraken2 databases: [benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2)

---

## Reference Selection and Assemblies

| Parameter                           | Description                                                                                                                                                                                                                                   |
| ----------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--reference_fasta <path>`          | Local FASTA file, space-separated list of FASTA files, or directory. Header format must be `>accession description` (space between accession and organism name).                                                                              |
| `--recursive_reference`             | Search `--reference_fasta` directory recursively.                                                                                                                                                                                             |
| `--get_variants`                    | Run mpileup for variant calling. Auto-enabled with `--reference_assembly`.                                                                                                                                                                    |
| `--assembly <path_or_url>`          | Local path to a pre-downloaded NCBI **RefSeq** assembly summary file (`assembly_summary_refseq.txt`). Skips the live download; required for air-gapped runs.                                                                                  |
| `--assembly_summary_genbank <path>` | Local path to a pre-downloaded NCBI **GenBank** assembly summary file (`assembly_summary_genbank.txt`). Only meaningful when `--assembly` is also provided — supplies the GenBank half of the lookup table without downloading it at runtime. |
| `--enable_genbank`                  | Also download the NCBI GenBank assembly summary at runtime (in addition to RefSeq). Extends the taxid resolution lookup to GenBank-only accessions (GCA entries not mirrored in RefSeq). Disabled by default; adds ~500 MB to startup I/O.    |
| `--custom_accession_map <path>`     | TSV file with `accession` and `taxid` columns for manually mapping non-standard accessions not found in any assembly summary or NCBI database. Applied after the assembly summary lookup and before the NCBI E-utilities backup.              |
| `--use_denovo`                      | Enable de novo assembly (MEGAHIT for Illumina, FLYE for ONT). Disabled by default.                                                                                                                                                            |
| `--use_megahit_longreads`           | Use MEGAHIT instead of FLYE for ONT de novo assembly.                                                                                                                                                                                         |
| `--reference_assembly`              | Enable reference-based assembly with variant calling. Disabled by default.                                                                                                                                                                    |
| `--generate_iss`                    | Generate InSilicoSeq (Illumina) simulated reads for benchmarking.                                                                                                                                                                             |
| `--generate_nanosim`                | Generate NanoSim (ONT) simulated reads for benchmarking.                                                                                                                                                                                      |
| `--nanosim_training <path>`         | Required with `--generate_nanosim`. Path to the NanoSim model directory.                                                                                                                                                                      |
| `--thresholds_json <path>`          | JSON file with optimal TASS weights per sample type and platform. Defaults to `assets/sampletype_best_thresholds.json`.                                                                                                                       |

> ⚠️ The `--reference_fasta` FASTA header **must** follow NCBI's format: `>NC_003663.2 Cowpox virus, complete genome`. The accession is followed by a space and then the organism description.

### Taxid Resolution Order

For every genome pulled into the pipeline, TaxTriage needs to resolve the sequence accession to a numeric NCBI taxid. It works through the following steps in order, stopping as soon as a match is found:

1. **RefSeq assembly summary** — The primary source. The GCF/GCA assembly accession is looked up against the `#assembly_accession` column in `assembly_summary_refseq.txt` (downloaded live, or supplied via `--assembly`). This resolves the vast majority of standard reference genomes.

2. **GenBank assembly summary** — If `--enable_genbank` is set (or `--assembly_summary_genbank` is provided alongside `--assembly`), the same lookup is repeated against `assembly_summary_genbank.txt`. GenBank contains a much larger set of GCA-prefixed assemblies that are not mirrored in RefSeq; enabling this catches those genomes at the cost of a larger download and slightly longer startup.

3. **Custom accession map (`--custom_accession_map`)** — A user-supplied TSV with `accession` and `taxid` columns. Applied to any accession that still has no match after both assembly summary lookups. This is the right tool for local reference FASTAs that use non-standard or private accessions (e.g. `CUSTOM_SEQ_001`), or for nuccore accessions (like `OR833055.1`) that exist in NCBI's sequence databases but have no corresponding assembly summary entry. The `accession` column must match the nuccore accession in the FASTA header (not the GCF/GCA assembly ID).

   Example TSV:

   ```
   accession	taxid
   OR833055.1	2697049
   CUSTOM_SEQ_001	12345
   ```

4. **NCBI E-utilities fallback** — Always active when internet access is available. For any accession still unresolved after steps 1–3, TaxTriage queries NCBI's `esearch` + `esummary` E-utilities APIs by nuccore accession to retrieve the taxid directly. Providing `--email` (your registered NCBI email) increases the rate limit from 3 to 10 requests/second. Failures are non-fatal — the taxid is left blank and a warning is printed rather than aborting the run. On air-gapped systems this step silently produces no results; use `--custom_accession_map` as a substitute.

> **Air-gapped tip:** Use `--assembly` + `--assembly_summary_genbank` (if needed) to supply pre-downloaded summary files, and `--custom_accession_map` for any accessions not covered by those summaries. This fully replaces all network-dependent taxid resolution.

### Assembly Selection Order

The [Taxid Resolution Order](#taxid-resolution-order) above runs **accession → taxid**. This section is the reverse direction — **taxid → assembly** — i.e. deciding _which_ genome to download for each top-hit organism (`download_fastas.py`, used by the `DOWNLOAD_ASSEMBLY` module). TaxTriage picks an assembly in this order, stopping at the first that resolves:

1. **Curated accession (accession-first)** — If the organism's row in the pathogen sheet (`assets/pathogen_sheet.csv`, or your `--pathogens` sheet) has a non-empty `assembly_accession`, that exact assembly (`GCF_*` / `GCA_*`) is used, provided it exists in the supplied assembly summary. This lets you **pin a specific, validated genome** per organism and bypass the heuristic ranking below. Matching is keyed by both taxid and organism name, so it works in taxid mode and in `--fuzzy` (name) mode.

2. **Local taxid match with priority tiers** — Otherwise, every assembly whose taxid matches is ranked and the best is chosen, highest priority first:

   | Priority | `assembly_level` / category             |
   | -------- | --------------------------------------- |
   | 0 (best) | `representative genome`                 |
   | 1        | `reference genome`                      |
   | 2        | `Complete Genome`                       |
   | 3        | anything else (scaffold / contig level) |

   Ties — and equal-tier candidates — are broken by **source order**: the RefSeq summary (`assembly_summary_refseq.txt`) is consulted before GenBank (`assembly_summary_genbank.txt`), so a RefSeq assembly is preferred over an equivalent GenBank one. (`--enable_genbank` only _adds_ GenBank as a lower-priority source; it never overrides a RefSeq match.)

3. **NCBI E-utilities fallback (`--ncbi_fallback`)** — For any organism still unresolved after the local summary lookup, `download_fastas.py` can query NCBI's `esearch` + `esummary` APIs — by accession first, then by taxid — to locate an assembly. This is used by the offline database build; provide `--email` to raise the NCBI rate limit. It is left off inside the pipeline containers by default so runs stay offline-safe; the local taxid match (step 2) is the in-pipeline fallback.

> The `assembly_accession` column is populated automatically when the curated database is (re)built (`add_assembly_column.py`), and is appended as the **last** column so the sheet's other columns keep their positions. Rows with no downloaded assembly are left blank and simply fall through to step 2.

---

## Alignment and Mapping

| Parameter                         | Description                                                                                                                                                                               |
| --------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--minmapq <N>`                   | Minimum mapping quality (MAPQ) for aligned reads. Default: `5`.                                                                                                                           |
| `--min_mapq_host <N>`             | Minimum MAPQ for host read removal. Only active with `--remove_reference_file`.                                                                                                           |
| `--mmap2_window <N>`              | Minimap2 window size. Default: `10`.                                                                                                                                                      |
| `--mmap2_fraction_filter <float>` | Minimap2 fraction filter for repetitive minimizers. Default: `0.0009`.                                                                                                                    |
| `--mmap2_I <size>`                | Minimap2 `-I` index batch size (e.g. `4G`). Overrides the automatic per-attempt default. Controls how much of the reference is loaded into RAM at once — smaller = less RAM, more passes. |
| `--mmap2_K <size>`                | Minimap2 `-K` query/minibatch size (e.g. `100M`). Overrides the automatic per-attempt default. Smaller = less RAM.                                                                        |

> **Automatic RAM back-off.** When `--mmap2_I` / `--mmap2_K` are not set, minimap2 scales both down on each retry so out-of-memory failures (e.g. exit code `143`) succeed on a later attempt without manual tuning:
>
> | Attempt | `-I` (index) | `-K` (query) | Memory request |
> | ------- | ------------ | ------------ | -------------- |
> | 1       | 8G           | 500M         | 24 GB          |
> | 2       | 4G           | 200M         | 40 GB          |
> | 3       | 2G           | 100M         | 56 GB          |
> | 4       | 1G           | 50M          | 72 GB          |
>
> `MINIMAP2_ALIGN` / `FILTER_MINIMAP2` retry up to 3 times (4 attempts total). User-supplied `--mmap2_I` / `--mmap2_K` always take priority over these defaults.
> | `--skip_realignment` | Skip the realignment step. No alignment metrics report will be generated. |
> | `--use_bt2` | Use Bowtie2 instead of minimap2 for Illumina reads. |
> | `--use_hisat2` | Use Hisat2 instead of minimap2 for Illumina reads. |
> | `--no_split_prefix` | Disable minimap2 `split-prefix` mode. Reduces RAM savings but may improve speed. |
> | `--loose` | Less sensitive false-positive removal — compares signatures across aligned queries only rather than whole genomes. Useful for samples with many variants. Slower and more I/O intensive. |

---

## Host Removal

| Parameter                          | Description                                                                                                |
| ---------------------------------- | ---------------------------------------------------------------------------------------------------------- |
| `--genome <key>`                   | Auto-download an iGenomes reference for host removal (e.g., `GRCh37`, `BDGP6`, `GRCz10`).                  |
| `--remove_reference_file <path>`   | FASTA file — reads aligned to these accessions are removed.                                                |
| `--include_singletons_hostremoval` | Retain singleton reads during paired-end host removal. Default: `FALSE`. Ignored for MEGAHIT/diamond runs. |

---

## Pathogen Discovery and Reporting Filters

> Primary pathogens are always shown; the `--show_*` flags only reveal the lower-severity categories. See [Microbial Categories](microbial-categories.md) for what each category means and how it is resolved per body site.

| Parameter                  | Description                                                                                                                              |
| -------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `--pathogens <path>`       | Custom pathogen annotation sheet (same format as `assets/pathogen_sheet.csv`).                                                           |
| `--get_pathogens`          | Use the pathogen FASTA file for alignment in addition to (or instead of if `--skip_kraken2` is used) Kraken2 top hits. Default: `FALSE`. |
| `--remove_taxids "<ids>"`  | Space-separated list of taxids to remove from Kraken2 output. Use `"9606"` to remove human.                                              |
| `--compress_species`       | Collapse all ODR annotations to species level. Default: `FALSE`.                                                                         |
| `--show_commensals`        | Show commensal organisms in the final report. Default: hidden.                                                                           |
| `--show_potentials`        | Show potential pathogens in the final report. Default: hidden.                                                                           |
| `--show_opportunistics`    | Show opportunistic pathogens in the final report. Default: hidden.                                                                       |
| `--show_unidentified`      | Show unannotated organisms in the final report. Default: hidden.                                                                         |
| `--remove_commensal`       | Drop commensal taxids from the candidate list before alignment. Ignored for sterile sample types. Default: `FALSE`.                      |
| `--add_irregular_top_hits` | Add organisms with irregular abundance distributions (z-score outliers) to top hits. Default: `false`.                                   |
| `--distributions <path>`   | Custom `.tsv.gz` or `.tsv` abundance distributions file (replaces HMP defaults).                                                         |
| `--min_conf <float>`       | Minimum TASS confidence score to appear in the PDF report.                                                                               |
| `--negative <name>`        | Override all non-control samples to use this negative control name (from samplesheet).                                                   |
| `--positive <name>`        | Override all non-control samples to use this positive control name.                                                                      |

---

## TASS Scoring and Threshold Parameters

| Parameter                          | Description                                                                                                                                                                                                                                                                                                                         |
| ---------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--gini_weight <float>`            | Weight for the Gini coefficient (coverage distribution inequality).                                                                                                                                                                                                                                                                 |
| `--breadth_weight <float>`         | Weight for log-based breadth of coverage score.                                                                                                                                                                                                                                                                                     |
| `--minhash_weight <float>`         | Weight for minhash false-positive filtering score.                                                                                                                                                                                                                                                                                  |
| `--mapq_weight <float>`            | Weight for MAPQ normalization score. Disabled by default.                                                                                                                                                                                                                                                                           |
| `--hmp_weight <float>`             | Weight for HMP abundance percentile. Recommend disabling for unknown sample types.                                                                                                                                                                                                                                                  |
| `--disparity_score_weight <float>` | Bonus weight for more prevalent organisms. Unused by default.                                                                                                                                                                                                                                                                       |
| `--dispersion_factor <float>`      | Dispersion factor in TASS score calculation.                                                                                                                                                                                                                                                                                        |
| `--reward_factor <float>`          | Reward factor in TASS score calculation.                                                                                                                                                                                                                                                                                            |
| `--disable_auto_weights`           | Disable clinically derived best cutoffs. Uses manual or default values instead.                                                                                                                                                                                                                                                     |
| `--ani_threshold <float>`          | ANI threshold for strain/species similarity. Default: `0.95` (95%).                                                                                                                                                                                                                                                                 |
| `--rep_breadth_min_frac <float>`   | Minimum share of a group's total reference length an accession must represent before it may set that group's representative coverage (the per-group **maximum** of covered-bp ÷ reference length). Default: `0.01` (1%). Set to `0` to disable.                                                                                     |
| `--rep_breadth_min_len <int>`      | Absolute minimum accession length (bp) before an accession may set its group's representative coverage. Applied together with `--rep_breadth_min_frac` — an accession qualifies if it clears **either** bar. Default: `50000`. Set to `0` to disable.                                                                               |
| `--no_read_removal`                | Disable conflict read removal for scoring. Conflict regions and removal stats are still computed and written to disk, but no reads are removed — coverage, Gini, `numreads` and TASS all derive from the complete original alignment set. Default: `false`. See [TASS Scoring 3.6](tass-scoring.md#cli-flags-that-control-removal). |
| `--taxid_removal_stats`            | Write the taxid-aggregated `removal_stats_by_taxid.xlsx` alongside the standard `removal_stats.xlsx`, grouping accessions by taxid. Default: `true`; set `--taxid_removal_stats false` to skip the extra report.                                                                                                                    |

See [TASS Scoring](tass-scoring.md) for a full explanation of how each metric contributes.

!!! note "Some scoring knobs are `match_paths.py`-only"

    `--no_read_removal` and `--taxid_removal_stats` above are pipeline parameters and work with `nextflow run`. The remaining `match_paths.py` scoring knobs — including `--skip_read_removal_scoring`, `--dominance_protect_ratio`, `--exclude_descriptions`, the benchmarking reports (`--report_metrics`, `--report_confusion_xlsx`) and the control-sample flags — are **not** plumbed through, so passing them to `nextflow run` has no effect. They apply when you invoke `match_paths.py` directly, the usual reason being an on/off scoring comparison or a benchmarking run. See [TASS Scoring 3.6 → CLI flags that control removal](tass-scoring.md#cli-flags-that-control-removal) and [3.12](tass-scoring.md#312-reference-exclusion-removal-reporting-and-control-samples).

> **On the two `rep_breadth_*` parameters.** These guard against scaffold-level draft assemblies, where the reference is thousands of ~1 kb contigs rather than a few chromosomes. Without them, a single read covering one small scaffold sets the coverage for the entire species — inflating both the breadth and minhash terms at once. The defaults comfortably admit any chromosome-level reference; raise them if your panel contains large scaffolds that still act as noise magnets, or set both to `0` to restore the previous behaviour. See [TASS Scoring 11.3](tass-scoring.md#113-size-eligibility-for-the-representative-coverage-maximum).

---

## Novelty Detection (Reference-Free / Open-Set)

The novelty branch classifies the **de novo assembly** (all post-QC reads, aligned or not) with a pluggable backend against a broad database, surfacing organisms the reference panel never aligned. By default it predicts genes with Pyrodigal first and classifies those (`--disable_gene` classifies the whole contigs). It is additive and never alters TASS scores. See [Novelty Detection](novelty-detection.md) for the full walkthrough and the per-backend database table.

| Parameter                                  | Default                     | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| ------------------------------------------ | --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--novelty <backend>`                      | `false`                     | Master switch **and** backend selector: `mmseqs2` (translated LCA), `kaiju` (protein translated search), or `bracken` (kraken2 + bracken abundance). `false` disables the branch.                                                                                                                                                                                                                                                                                   |
| `--disable_gene`                           | `false`                     | Skip the default Pyrodigal ORF-prediction step and classify the whole de novo contigs instead of predicted genes. By default (flag off) genes are predicted first, so a near-complete genome contributes several LCA rows; with `--disable_gene` each contig is one query.                                                                                                                                                                                          |
| `--novelty_db <name\|alias\|path>`         | `Kalamari`                  | Database for the chosen backend. mmseqs2: seqTaxDB name or dir. kaiju: alias (`test`, `viruses`/`viral`, `nr`, `nr_euk`, `refseq`, `refseq_nr`, `refseq_ref`, `progenomes`, `fungi`, `plasmids`, `rvdb`), `*.tgz`, URL, or dir. bracken: alias (`standard[_8gb\|_16gb]`, `viral`, `minusb`, `pluspf[_8gb\|_16gb]`, `pluspfp[_8gb\|_16gb]`, `core_nt`, `gtdb`, `eupathdb`, `ncbi_reference`; `k2_` prefix optional), `*.tar.gz`, URL, or dir. Aliases auto-download. |
| `--novelty_db_cache <path>`                | `${projectDir}/dbs/mmseqs`  | `storeDir` cache for the mmseqs2 download.                                                                                                                                                                                                                                                                                                                                                                                                                          |
| `--novelty_db_prefix <str>`                | `seqTaxDB`                  | Base name of the mmseqs2 database files inside the db directory.                                                                                                                                                                                                                                                                                                                                                                                                    |
| `--novelty_kaiju_db_cache <path>`          | `${projectDir}/dbs/kaiju`   | `storeDir` cache for the kaiju download.                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `--novelty_kraken2_db_cache <path>`        | `${projectDir}/dbs/kraken2` | `storeDir` cache for the bracken/kraken2 download.                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `--novelty_kaiju_db_baseurl <url>`         | kaiju S3 bucket             | Index bucket kaiju aliases are fetched from.                                                                                                                                                                                                                                                                                                                                                                                                                        |
| `--novelty_kraken2_db_baseurl <url>`       | genome-idx S3 bucket        | Index bucket kraken2/bracken aliases are fetched from.                                                                                                                                                                                                                                                                                                                                                                                                              |
| `--novelty_kaiju_mode <greedy\|mem>`       | `greedy`                    | Kaiju run mode. `greedy` (`-a greedy -e 5`) is more sensitive on divergent contigs; `mem` (`-a mem`) uses exact matches only and avoids long conserved contigs collapsing to a genus LCA. kaiju backend only.                                                                                                                                                                                                                                                       |
| `--novelty_weights <w_dark,w_rank,w_idnt>` | `0.5,0.3,0.2`               | Weights for the three novelty-score signals: dark fraction, high-rank-only fraction, low-identity tail.                                                                                                                                                                                                                                                                                                                                                             |
| `--novelty_flag_z <float>`                 | `2.0`                       | Score threshold (≈ z-score) above which a sample is flagged as novel.                                                                                                                                                                                                                                                                                                                                                                                               |
| `--novelty_idnt_cut <float>`               | `50.0`                      | Amino-acid % identity below which a best hit counts toward the low-identity tail mass (mmseqs2 pident).                                                                                                                                                                                                                                                                                                                                                             |
| `--novelty_min_reads <int>`                | `2`                         | Candidate floor: minimum hits (contigs for kaiju/mmseqs2, reads for bracken) to report a taxon. `1` surfaces single-contig hits and disables the depth-scaled cutoff.                                                                                                                                                                                                                                                                                               |
| `--novelty_bracken_readlen <int>`          | `100`                       | bracken `-r` read length (kmer distribution). bracken backend only.                                                                                                                                                                                                                                                                                                                                                                                                 |
| `--novelty_bracken_level <str>`            | `S`                         | bracken `-l` rank level (S = species). bracken backend only.                                                                                                                                                                                                                                                                                                                                                                                                        |

**Candidate-gate.** Whether a genus+ taxon is reported as a candidate is controlled by `novelty_score.py`, whose effective cutoff is `max(--novelty_min_reads, ceil(min_cand_frac × placeable_hits))`. The floor is the top-level `--novelty_min_reads` (default `2`). The depth-scaling fraction `min_cand_frac` (default `0.002`) is a script default, not a pipeline flag; the `NOVELTY_SCORE` module forces it to `0` when `--novelty_min_reads ≤ 1`. For very shallow residuals, set `--novelty_min_reads 1`.

> `--novelty` requires de novo contigs, so pair it with `--use_denovo`. It works well alongside a sparse `--reference_fasta`/`--organisms` panel, where the closed-set table is intentionally limited and the novelty branch recovers higher-rank genera. Report-side, the [Detection Rescue](detection-rescue.md) toggles can re-surface sub-threshold and novelty-supported detections without changing any score.

---

## Annotation and Functional Analysis

| Parameter                    | Description                                                                                                        |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `--annotate`                 | Enable VF/AMR annotation (also triggers de novo assembly). Default: `false`.                                       |
| `--annotate_proteins <path>` | FAA file of proteins for annotation. Default: `assets/bvbrc_specialty_genes_with_sequences_taxids_and_sites.faa`.  |
| `--annotate_meta <path>`     | Metadata TSV for annotation proteins. Default: `assets/bvbrc_specialty_genes_with_sequences_taxids_and_sites.tsv`. |
| `--pident <float>`           | Minimum percent identity for annotation hits. Default: `75`.                                                       |

> VF/AMR annotation runs on each sample's de novo contigs. Samples that assembled contigs but aligned **no** reference now carry their annotation into the comparison report's Summary and VF/AMR views (previously only the Novelty tab showed signal for them). This is automatic when `--annotate` is set — see [Detection Rescue → Annotation of unaligned samples](detection-rescue.md#annotation-of-unaligned-samples).

---

## MicroBERT and AI/ML Predictions

| Parameter                | Description                                                                      |
| ------------------------ | -------------------------------------------------------------------------------- |
| `--microbert <path>`     | Path to a MicroBERT model directory to enable AI/ML-based taxonomic predictions. |
| `--microbert_maxlen <N>` | Maximum sequence chunk length for MicroBERT. Default: `2000` (2 kb).             |

See [Pipeline Modules → MicroBERT](pipeline-modules.md#microbert-module-optional) for setup details.

When `--microbert` is supplied, the pipeline clusters each sample's aligned reads (`MMSEQS_EASYCLUSTER`), runs the model over the cluster representatives (`MICROBERT_PREDICT`), and produces a per-accession confidence profile (`MICROBERT_PARSE` → `*.microbert.accession.tsv`). The model's per-sequence top-1 probability is aggregated **by the reference accession each read aligned to** (read-support-weighted), so the score attaches to the _detected organism_ rather than to the model's own predicted-taxid namespace. `match_paths.py` then joins this profile to each organism via its reference accession(s) and stores it as `mmbert` / `mmbert_model`.

This surfaces two columns in the reports:

| Column                    | Meaning                                                                                               |
| ------------------------- | ----------------------------------------------------------------------------------------------------- |
| **MicrobeRT Probability** | Model confidence (0–100 %) averaged over the reads aligned to that organism's reference(s).           |
| **MicrobeRT Model**       | Name of the MicroBERT model directory used (e.g. `VHDB-nucleotide-transformer-v2-50m-multi-species`). |

They appear in the [Interactive Report](interactive-report.md) **Table** tab (visible by default) and the `mmbert%` column of the Organism Discovery Report PDF.

> **Model domain matters.** The probability reflects what the _chosen model_ was trained to recognise. A viral model (e.g. VHDB) scores every organism's reads on a "virus-likeness" basis, so a bacterium may receive a moderate score and a virus a high one — interpret the value relative to the model, not as an absolute identity confidence.

---

## Output, Reporting, and Visualization

| Parameter                      | Description                                                                                                                                                                                                                                                                                                                           |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--skip_plots`                 | Skip all plots.                                                                                                                                                                                                                                                                                                                       |
| `--skip_multiqc`               | Skip MultiQC report generation.                                                                                                                                                                                                                                                                                                       |
| `--skip_matrix`                | Skip ANI similarity matrix (faster compute).                                                                                                                                                                                                                                                                                          |
| `--integrate_strains_table`    | Merge strain rows into the same table row as the parent species. Default: `false` (separate table).                                                                                                                                                                                                                                   |
| `--sort_alphabetical`          | Sort PDF tables alphabetically by organism name instead of High Consequence then TASS score.                                                                                                                                                                                                                                          |
| `--no_subkey`                  | Do not split organisms into species and strain sub-tables.                                                                                                                                                                                                                                                                            |
| `--offline_report`             | Download the interactive report's CDN libraries (D3, xlsx, jsPDF, Leaflet, Font Awesome) at build time and embed them inline, so `all.odr.html` opens without internet. Requires network on the **pipeline** host. Default: `false`.                                                                                                  |
| `--offline_report_files <dir>` | Directory of local copies of those libraries (and their fonts/marker images) to embed inline — a fully offline build with **no network**. Takes precedence over `--offline_report`. Prepare it with `python scripts/fetch_offline_report_libs.py`. See [Interactive Report → Offline Reports](interactive-report.md#offline-reports). |

---

## Report Sample-QC Flags

Whole-**sample** criteria that seed the interactive report's **Sample QC / Flags** rule set. These are report _defaults only_: every sample still runs through the full workflow and appears in every output. A sample that matches is marked in the Heatmap, Table, Metadata & Mapping and Summary tabs and, with `--report_flag_action hide`, also removed from those views (reversibly — the report keeps a one-click toggle). Users can edit, add to or delete every rule live in the report.

!!! note "What the counts count"

    The detection-based criteria (`--report_flag_min_organisms`, `--report_flag_min_detections`) are computed from the **whole dataset** in the report, not from whatever the report is currently displaying: an organism suppressed by the TASS slider, a name filter or a hidden column still counts, as long as it clears the criterion's own TASS cutoff. Only two things are left out — species and genus roll-up rows (counting them would inflate every distinct-organism figure two- or three-fold), and anything on `--report_flag_exclude_taxids`, which defaults to human. Detections the *pipeline* never wrote are of course not there to count.

    These sample flags are also entirely separate from the report's per-organism **follow-up list** (the star on a detection): one marks samples, the other marks organism–sample pairs, and neither reads the other's state.

| Parameter                              | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| -------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--report_flag_min_reads <n>`          | Flag any sample with fewer than `n` total reads.                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `--report_flag_min_aligned_reads <n>`  | Flag any sample with fewer than `n` aligned reads. Set to `1` to catch samples that produced no alignments at all.                                                                                                                                                                                                                                                                                                                                                          |
| `--report_flag_min_organisms <n>`      | Flag any sample with fewer than `n` **distinct** organisms scoring at or above `--report_flag_organism_tass`.                                                                                                                                                                                                                                                                                                                                                               |
| `--report_flag_organism_tass <tass>`   | TASS cutoff (0–100) used by `--report_flag_min_organisms`. Defaults to `--min_conf` when set, else `75`.                                                                                                                                                                                                                                                                                                                                                                    |
| `--report_flag_min_detections <n>`     | Flag any sample with fewer than `n` detections passing their own threshold.                                                                                                                                                                                                                                                                                                                                                                                                 |
| `--report_flag_metadata <spec>`        | Metadata criteria, `;`-separated. Each clause is `field:op:value` (or the shorthand `field=value`, meaning equals). Operators: `==` `!=` `contains` `!contains` `regex` `empty` `!empty` `<` `<=` `>` `>=`. The field is looked up in the run metadata first, then in the sample metadata.                                                                                                                                                                                  |
| `--report_flag_logic <any\|all>`       | Flag when **any** criterion matches (default) or only when **all** of them do.                                                                                                                                                                                                                                                                                                                                                                                              |
| `--report_flag_action <flag\|hide>`    | What happens to a matching sample: `flag` marks it but keeps it visible (default); `hide` also removes it from every chart and table.                                                                                                                                                                                                                                                                                                                                       |
| `--report_flag_missing`                | Treat a missing or blank value as a match. Off by default, so a criterion can never flag a sample purely because the field was never populated.                                                                                                                                                                                                                                                                                                                             |
| `--report_flag_view <all\|hide\|only>` | Which samples the report **opens** on. `all` (default) shows everything, with a rule's own `hide` action still applying; `hide` hides every flagged sample; `only` inverts it — the flagged samples are shown and everything that passed is hidden. Switchable at any time from the dropdown beside **Filters…** in the report's Sample QC panel (All samples / Hide flagged / Only flagged), and never destructive: nothing is dropped from the data, only from the views. |
| `--report_flag_exclude_taxids <ids>`   | Taxids or organism names that never count toward a detection figure — the distinct-organism counts, the detection counts and any aggregated column. Defaults to `9606` (human), so `--report_flag_min_organisms` means the same thing on a dehosted run and on one where human is still in the table. Comma- or space-separated; taxids match exactly, names match anywhere in the organism / species / genus text. Pass `''` to count host like any other organism.        |
| `--report_flag_rules <file.json>`      | A JSON file holding a full rule list. Replaces every other `--report_flag_*` parameter.                                                                                                                                                                                                                                                                                                                                                                                     |

Example — flag anything shallow, uninformative, or missing its collection site:

```bash
nextflow run . -profile docker \
    --input samplesheet.csv \
    --report_flag_min_reads 500000 \
    --report_flag_min_organisms 1 --report_flag_organism_tass 75 \
    --report_flag_metadata "environmental_site:empty:;host_disease:contains:influenza"
```

### Classifier vs alignment

Three of the `derived` fields compare what Kraken2/Centrifuge assigned (the `K2 Reads` column) against what the aligner corroborated (`# Reads Aligned`) — the gap between the two is the usual signature of a database artefact or a conserved-region pile-up, and neither column shows it alone:

| Field                      | Meaning                                                                                                                                                                                                                                             |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `k2_reads_sum`             | Total classifier reads across the sample's detections. The counterpart to `reads_aligned_sum`.                                                                                                                                                      |
| `aligned_to_k2_ratio`      | `reads_aligned_sum ÷ k2_reads_sum`. Well below 1 means the classifier is carrying the sample on its own. **Missing** (never `0`) when the sample has no classifier reads at all, so a run without a classifier column cannot trip a ratio rule.     |
| `unsupported_k2_organisms` | How many distinct organisms the classifier called loudly and the aligner did not back — at least `k2min` classifier reads (default 50) but fewer than 5% of them aligned. Per-organism, so one bad call in an otherwise healthy sample still shows. |

`unsupported_k2_organisms` is usually the more actionable of the three: the sample-level ratio tells you something is off, this tells you how many organisms are responsible. Both honour `--report_flag_exclude_taxids`, which matters here — host is often the loudest classifier call in the sample.

```bash
nextflow run . -profile docker \
    --input samplesheet.csv \
    --report_flag_rules classifier_qc.json     # unsupported_k2_organisms >= 1
```

### Rules file format

`--report_flag_rules` takes either a bare list of rule objects or an object with `logic`, `missing_fails`, `view`, `exclude_taxids` and `rules` (the two report-wide settings may also be given here instead of on the command line):

```json
{
  "logic": "any",
  "missing_fails": false,
  "view": "all",
  "exclude_taxids": ["9606"],
  "rules": [
    { "source": "meta", "field": "total_reads", "op": "<", "value": 500000, "action": "hide" },
    { "source": "derived", "field": "unique_taxids_above_tass", "op": "<", "value": 1, "tass": 75, "action": "flag" },
    { "source": "runmeta", "field": "sample_origin_country", "op": "empty", "value": "", "action": "flag" },
    { "source": "data", "field": "Coverage", "agg": "max", "op": "<", "value": 10, "action": "flag" },
    { "source": "derived", "field": "unsupported_k2_organisms", "op": ">=", "value": 1, "k2min": 50, "action": "flag" }
  ]
}
```

`source` is one of:

| Source    | Reads                                                                                                                                                                                                                                                                                       |
| --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `meta`    | A per-sample pipeline metric: `total_reads`, `aligned_reads`, `total_organism_reads`, `num_keys`, `num_subkeys`, `num_species_groups`, `platform`, `sample_type`, `control_type`, …                                                                                                         |
| `derived` | A count computed from the detections: `unique_taxids_above_tass` (with `tass`), `unique_taxids`, `passing_detections`, `detections`, `high_consequence`, `unique_genera`, `max_tass`, `reads_aligned_sum`, `k2_reads_sum`, `aligned_to_k2_ratio`, `unsupported_k2_organisms` (with `k2min`) |
| `runmeta` | Any metadata column (run metadata first, sample metadata as a fallback)                                                                                                                                                                                                                     |
| `data`    | Any numeric detection column, rolled up per sample by `agg` (`max`, `min`, `mean`, `sum`, `count`)                                                                                                                                                                                          |

A malformed rule is reported on stderr and skipped rather than failing the run. See [Interactive Report → Sample QC flags](interactive-report.md#sample-qc-flags).

---

## Deprecated Parameters

These parameters are no longer actively supported:

| Parameter        | Description                                                            |
| ---------------- | ---------------------------------------------------------------------- |
| `--metaphlan`    | Use MetaPhlAn. Requires a `.pkl` and `.btl2` Bowtie2 index path.       |
| `--use_diamond`  | Use DIAMOND on features pulled from NCBI post de novo assembly.        |
| `--get_features` | Pull CDS features from each assembly and compute per-feature coverage. |

---

## Nextflow Core Arguments

These are Nextflow-native flags (single hyphen):

| Flag               | Description                                          |
| ------------------ | ---------------------------------------------------- |
| `-profile`         | Comma-separated list of execution profiles.          |
| `-resume`          | Resume from cached results.                          |
| `-r <branch/tag>`  | Branch or version tag.                               |
| `-latest`          | Autopull the latest commit for the specified branch. |
| `-bg`              | Run in background, detached from terminal.           |
| `-c <config>`      | Supply a custom Nextflow config file.                |
| `-work-dir <path>` | Override working directory (default: `./work`).      |
