# Profiles

A profile is a bundle of pre-set parameters and resource ceilings, selected with
Nextflow's `-profile` flag. TaxTriage ships three kinds, and a real run normally
chains one of each:

```bash
nextflow run https://github.com/jhuapl-bio/taxtriage \
  -r stable -latest \
  -profile desktop,docker \
  --input samplesheet.csv --outdir results
```

| Kind          | Examples                                                            | What it sets                                                           |
| ------------- | ------------------------------------------------------------------- | ---------------------------------------------------------------------- |
| **Workload**  | `test`, `low`, `viral`, `desktop`, `deep`, `local`, `mce`, `seqera` | Resource ceilings, database choice, reporting and sensitivity defaults |
| **Container** | `docker`, `singularity`, `apptainer`, `podman`, `conda`, `wave`     | Which runtime executes each process                                    |
| **Executor**  | `slurm`, `sge`, `cloud`                                             | Where jobs are submitted                                               |

!!! note "Order matters"

    Profiles are applied left to right and later ones win. `-profile deep,docker`
    and `-profile docker,deep` are not the same command. Put the workload profile
    first and the container profile last, and remember that anything you pass on
    the CLI (`--db`, `--max_memory`, …) overrides all of them.

---

## Workload profiles at a glance

| Profile                 | CPUs / RAM / time ceiling   | Default `--db`              | Built for                                                                                                 |
| ----------------------- | --------------------------- | --------------------------- | --------------------------------------------------------------------------------------------------------- |
| `test`                  | 2 / 7 GB / 10 h             | `test` (112 MB)             | CI and install verification. Pulls a tiny bundled samplesheet.                                            |
| `test_viral`            | 2 / 3 GB / 10 h             | `viral` (~550 MB)           | Same, but exercising the viral database path.                                                             |
| `low`                   | 4 / 8 GB / 20 h             | `standard8` (~7.5 GB)       | Laptops and small VMs. `low_memory = true`, so Kraken2 reads the DB from disk instead of RAM.             |
| `viral`                 | 4 / 8 GB / 20 h             | `viral` (~550 MB)           | Virus-only triage on modest hardware. `min_reads_align = 1`, `sensitive = true`.                          |
| `desktop`               | 6 / 24 GB / 12 h            | `pluspfp16` (~16 GB)        | Mid-tier workstation doing routine bacterial + fungal + protozoan triage. Annotation on, de novo off.     |
| `local`                 | 4 / 100 GB / 24 h           | `standard` (~80 GB)         | A workstation or HPC node with real memory behind it. Annotation on, full standard DB in RAM.             |
| `deep`                  | 32 / 256 GB / 72 h          | `standard` (~80 GB)         | Multi-GB FASTQ (deep WGS, high-output NovaSeq, PromethION). Turns on read reduction see below.            |
| `mce`                   | 16 CPUs / 24 h (no RAM cap) | `standard` + pathogen FASTA | MCE-style runs that align against the curated pathogen reference set. `min_conf = 0.4`, `strains = true`. |
| `seqera` / `mce_seqera` | 32 / 500 GB / 32 h          | S3-hosted prebuilt DB       | Seqera Platform / AWS Batch. Reads databases and the host reference from S3 rather than downloading.      |

The ceilings above are `resourceLimits` a cap, not a request. A process asking
for 4 CPUs still gets 4 on the `deep` profile; the cap only stops a process from
asking for more than the machine has.

---

## Choosing one

| If your situation is…                               | Start with          | Why                                                                                                                                  |
| --------------------------------------------------- | ------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| "I just installed this and want to see it run"      | `test,docker`       | Downloads ~112 MB, finishes in 10–15 minutes, needs no samplesheet of your own.                                                      |
| Laptop, 8–16 GB RAM, a handful of samples           | `low,docker`        | The 8 GB capped database plus `--low_memory` is the only combination that reliably fits. Expect reduced species-level resolution.    |
| Outbreak/clinical question that is viral only       | `viral,docker`      | Skips the cost of a bacterial database entirely and drops the alignment floor to 1 read so low-titre virus is not filtered out.      |
| Office workstation, 32 GB RAM, routine surveillance | `desktop,docker`    | Widest taxonomic coverage that still fits in a capped 16 GB index, with annotation enabled.                                          |
| HPC node or big workstation, 128 GB+ RAM            | `local,singularity` | Full standard database resident in RAM the accuracy baseline.                                                                        |
| Multi-GB FASTQ per sample                           | `deep,singularity`  | The only profile that normalises read depth before classification, plus much wider memory/time ladders and minimap2 index splitting. |
| Aligning against the curated pathogen sheet         | `mce,docker`        | Adds `--reference_fasta` for the pathogen set and strain-level output.                                                               |
| Seqera Platform / AWS Batch                         | `seqera`            | Points at S3-hosted databases so nothing is re-downloaded per run.                                                                   |
| SLURM cluster                                       | add `slurm`         | e.g. `-profile deep,singularity,slurm`.                                                                                              |

---

## What `deep` does differently

`deep` is the only profile that reduces read volume before classification, because
at multi-GB scale the classifier and aligner are otherwise intractable:

1. **bbnorm** (`--downsample`, on) k-mer coverage normalisation, **Illumina only**.
   Runs after trimming, before host removal. Flattens the high-abundance background
   while leaving rare organisms intact. `bbnorm_mindepth` is set to `2` rather than
   bbnorm's own default of `6`, which would delete exactly the low-abundance
   organisms triage is looking for.
2. **seqtk** (`--subsample`, 100 M reads) a blunt random cap, all platforms. This
   is the only reduction ONT/PacBio gets, since bbnorm's k-mer table is meaningless
   on noisy long reads.

!!! warning "bbnorm destroys relative abundance"

    If your interpretation depends on read-proportion abundance rather than
    presence/absence plus coverage, run `--downsample false` and rely on
    `--subsample` alone.

`deep` also loads `conf/deep_process.config`, which sets per-process memory ladders
that retry upward after an OOM (Kraken2 120 → 240 GB, bbnorm 96 → 192 GB, minimap2
64 → 160 GB). Those directives live in a separate file on purpose `withName:`
blocks written inside a profile are silently overwritten by `conf/modules.config`.

---

## Resource and database limits read before you launch

!!! danger "The database, not the pipeline, decides whether your run fits"

    Kraken2 memory-maps its index. Peak RAM for classification is roughly the
    on-disk size of the database, not a function of how many reads you have.
    A machine with 16 GB of RAM cannot hold the ~80 GB `standard` database no
    matter what `--max_memory` says, and the failure mode is an exit code `137`
    (OOM kill) well into the run after you have already paid for trimming and
    host removal.

**Database footprints.** Sizes below are approximate and grow with every upstream
release; check [benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2)
for the current figures, and see the
[downloadable database table](cli-parameters.md#supported-downloadable-databases).

| Database    | On disk ≈ RAM to load       | Coverage                                             |
| ----------- | --------------------------- | ---------------------------------------------------- |
| `viral`     | ~550 MB                     | Viruses only                                         |
| `standard8` | ~7.5 GB (hard-capped build) | Bacteria/archaea/viral/human, aggressively minimised |
| `pluspfp16` | ~16 GB (hard-capped build)  | Adds protozoa, fungi, plant                          |
| `standard`  | ~80 GB and rising           | Full RefSeq bacteria/archaea/viral/human             |
| `pluspfp`   | 100 GB+                     | Full standard plus protozoa, fungi, plant            |
| `core_nt`   | 100 GB+                     | Broadest, effectively HPC-only                       |

The capped `_8gb` / `_16gb` builds are not smaller downloads of the same data —
they are lossy, minimised indexes. They trade sensitivity and species-level
resolution for fitting on the machine. A species missing from a `standard8` run is
not evidence of absence.

**Disk, not just memory.** Plan for roughly **3× the database tarball** during
setup (download + extraction + the extracted copy), plus the Nextflow `work/`
directory, which on a deep run routinely exceeds the size of the input FASTQ
several times over. `cleanup` is deliberately left `false` so `-resume` works;
run `nextflow clean -f` once you are done inspecting a completed run.

**Downloads happen once, if you let them.** `--download_db` caches into a
per-name `storeDir`, so the multi-GB pull is skipped on later runs. Pointing
`--outdir` somewhere new each time, or running on ephemeral cloud storage,
re-downloads everything. On a shared or offline system, download once and pass
`--db /path/to/db` instead.

**Other things that will bite you:**

- **`--low_memory` is the escape hatch, and it is slow.** It streams the index
  from disk instead of loading it, which can turn minutes into hours but on
  a fast NVMe it is usually better than not running at all. `low` sets it; `viral`,
  `desktop`, `local` and `deep` do not.
- **Minimap2 has its own ceiling.** The reference index is separate from the
  Kraken2 database, and a high `--top_hits_count` (the `deep` profile pulls 50)
  can push it past a single index budget. `deep` sets `split_prefix = true` with
  `--mmap2_I 8G` for this reason. The pipeline also backs memory off automatically
  across retries. See [Troubleshooting](troubleshooting.md).
- **Ceilings are caps, not reservations.** A profile does not check that your
  machine actually has what it claims. `-profile deep` on a 32 GB box will happily
  submit a job requesting 120 GB and get it killed. Match the profile to the
  hardware, or override with `--max_cpus` / `--max_memory` / `--max_time`.
- **On a cluster, the profile is not the whole story.** Your scheduler's per-job
  memory and wall-time limits apply on top. A 72-hour `deep` run on a partition
  with a 24-hour limit dies at hour 24.
- **Test profiles are not starting points for real work.** `test`, `test_viral`
  and (despite the name in its own header) `low` and `viral` set
  `input` to a bundled example samplesheet. Always pass your own `--input`, and
  drop `test` entirely once you are past verification.

---

## Overriding a profile

Anything a profile sets can be overridden on the command line the CLI wins:

```bash
# desktop ceilings, but a bigger database and more headroom
nextflow run https://github.com/jhuapl-bio/taxtriage \
  -r stable -latest \
  -profile desktop,docker \
  --db standard --max_memory 96.GB \
  --input samplesheet.csv --outdir results
```

```bash
# deep, but preserving relative abundance
nextflow run https://github.com/jhuapl-bio/taxtriage \
  -r stable -latest \
  -profile deep,docker \
  --downsample false --subsample 50000000 \
  --input samplesheet.csv --outdir results
```

For anything you will reuse, put the overrides in your own `-c custom.config`
rather than in the command. See
[Running the Pipeline](running-the-pipeline.md#profiles) and
[CLI Parameters](cli-parameters.md).
