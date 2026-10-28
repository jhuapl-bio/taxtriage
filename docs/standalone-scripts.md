# Standalone Scripts

Not everything in TaxTriage has to run inside Nextflow. These scripts are useful on their own, for example to build references before a run, prepare test data, or post-process results. They live in two places:

| Location                                                                | What's there                                                                                                                                                    |
| ----------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`scripts/`](https://github.com/jhuapl-bio/taxtriage/tree/main/scripts) | Side tools the workflow never calls: reference builders (below) and developer helpers.                                                                          |
| [`bin/`](https://github.com/jhuapl-bio/taxtriage/tree/main/bin)         | Scripts the workflow runs in its processes, which also work from the command line. They stay in `bin/` because Nextflow puts that folder on each task's `PATH`. |

Run them from a TaxTriage checkout (`git clone https://github.com/jhuapl-bio/taxtriage`): several `bin/` scripts import their neighbours (`fasta_utils.py`, `determine_priority_assembly.py`, ...), so they need to stay next to them. All are Python ≥ 3.8; extra packages are listed per script, and every script prints its options with `--help`.

---

## Reference building

### `scripts/pull_contigs.py` - a FASTA from a list of contigs

TaxTriage pulls whole assemblies for taxids or names (`--organisms`, `--organisms_file`, `--get_pathogens`); it has no input for a list of individual contigs. To align against exactly a set of nuccore accessions, build the reference up front and pass the two outputs to the pipeline:

```bash
# contigs.txt: one accession per line, versioned or not ('#' lines and a header row are ignored)
python3 scripts/pull_contigs.py -i contigs.txt -o refs/ \
    --taxdump /path/to/taxdump -e you@org.org     # names.dmp + nodes.dmp the pipeline uses

nextflow run jhuapl-bio/taxtriage ... \
    --reference_fasta refs/contigs.fasta \
    --custom_accession_map refs/contigs.accession_map.tsv
```

`--custom_accession_map` is authoritative (see [Taxid Resolution Order](cli-parameters.md#taxid-resolution-order)), so every contig gets exactly the taxid in the map rather than one guessed from its FASTA header.

**Where the data comes from**

- **Sequences:** `--source-fasta` first (e.g. the prebuilt `pathogens.fasta`), then NCBI `efetch`, batched by total length. A plain FASTA with an up-to-date `.fai` is read by random access. A gzipped FASTA, or one whose `.fai` is older than the FASTA or doesn't match it, is streamed once instead.
- **Metadata:** `--source-map` first (a `pathogens.mapping.taxid.tsv`, a 4-column gcfmapping, or a custom accession map), then NCBI `esummary`. Contigs whose record names no assembly get a nuccore → assembly `elink`, in small batches (`--elink-batch-size`, default 20) that are split automatically when NCBI drops one; `--no-elink` skips it. `--assembly` (assembly summaries) fills in missing taxids and adds assembly-level columns.
- **Taxonomy:** with `--taxdump`, merged taxids are remapped, deleted or unknown taxids are flagged, and names and species-level taxids come from that dump. `--map-rank species` writes the species taxid to the map instead of the record's own (often strain-level) taxid. `--no-map-names` leaves the `name` column out so the pipeline names organisms from its own taxonomy.
- **Offline:** `--offline --source-fasta pathogens.fasta --source-map pathogens.mapping.taxid.tsv` never contacts NCBI.

**Resume.** Each batch of sequences is fsync'd to the FASTA and checkpointed in `<outdir>/.<prefix>.state/`, and NCBI metadata and assembly lookups are cached there too. Re-running the same command continues where it stopped: a half-written batch is truncated away, failed NCBI batches are retried, accessions added to the list are fetched, and removed ones are dropped from the FASTA. `--restart` starts over.

**Outputs** (prefix defaults to the input file name):

| File                                    | Contents                                                                                                                                                                                                                       |
| --------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `<prefix>.fasta`, `.fasta.fai`          | sequences (80-column, `>ACC.V description`) → `--reference_fasta`                                                                                                                                                              |
| `<prefix>.accession_map.tsv`            | `accession`, `taxid`, `name` → `--custom_accession_map`                                                                                                                                                                        |
| `<prefix>.mapping.taxid.tsv`            | the pipeline's map format: `Acc`, `Assembly`, `Organism_Name`, `Description`, `Mapped_Value`                                                                                                                                   |
| `<prefix>.metadata.tsv`                 | requested id, delivered id, assembly and its source, taxid and status (`ok` / `merged` / `deleted` / `not_in_taxdump`), original taxid, rank, species, length, sequence source, notes (`version_changed`, `taxid_merged`, ...) |
| `<prefix>.assemblies.txt`               | unique assembly accessions                                                                                                                                                                                                     |
| `<prefix>.missing.tsv`                  | each requested id that was not delivered or has no taxid, with the reason                                                                                                                                                      |
| `<prefix>.summary.json`, `<prefix>.log` | run counts; timestamped log, appended across resumes                                                                                                                                                                           |

`--strict` exits with status 1 if anything is missing or lacks a taxid. NCBI requests are rate-limited (faster with `--api-key` or `$NCBI_API_KEY`) and retried with backoff. Only the Python standard library is needed.

### `scripts/pathogen_db/build_pathogen_db.sh` - rebuild the pathogen reference

Builds the reference that `--get_pathogens` downloads (`pathogens.fasta.gz`) and its accession → taxid mapping from the [pathogen sheet](pathogen-sheet.md): one assembly per sheet taxid.

```bash
scripts/pathogen_db/build_pathogen_db.sh -o pathogen_build/ -e you@org.org -t 8
```

| Step                                                            | Tool                                  | Output                                                                 |
| --------------------------------------------------------------- | ------------------------------------- | ---------------------------------------------------------------------- |
| 1. assembly summary (downloaded if `-a` is not given)           | curl / wget                           | `assembly_summary_refseq.txt`                                          |
| 2. taxids + download input (`taxid⇥assembly_accession`)         | -                                     | `pathogens.taxids.txt`, `pathogens.download_input.tsv`                 |
| 3. download: sheet accession first, then taxid match, then NCBI | `bin/download_fastas.py`              | `pathogens.fasta`, `pathogens.gcfmapping.tsv`, `pathogens.missing.txt` |
| 4. index + compress                                             | `samtools faidx`, bgzip / pigz / gzip | `pathogens.fasta.fai`, `pathogens.fasta.gz`                            |
| 5. accession → taxid                                            | `bin/append_taxid.py`                 | `pathogens.mapping.taxid.tsv`                                          |
| 6. sheet name vs. downloaded organism                           | `validate_names.py`                   | `pathogens.mapping.taxid.validated.tsv`                                |
| 7. sheet + `assembly_accession`; next run's download input      | `add_assembly_column.py`              | `pathogen_sheet.with_assembly.csv`                                     |

Options: `-s` sheet (default `assets/pathogen_sheet.csv`), `-a` assembly summary, `-e` NCBI email (enables the NCBI fallback and taxid backup), `-t` threads, `-p` prefix, `--fresh-input` (rebuild the download input from the sheet instead of reusing the last run's), `--skip-download` (reuse `pathogens.fasta` and only rebuild the maps), `--no-compress`. Everything is logged to `pathogens.build.log`. Needs biopython and pandas; samtools and bgzip/pigz are optional.

The `.fai` is always rebuilt. An index left over from an older FASTA silently points at the wrong bytes for any tool that trusts it.

The two helpers also run on their own:

- `validate_names.py --mapsheet M --pathogens_sheet S -o OUT [--mismatches MM]` checks, for each sheet row, whether its taxid downloaded an assembly and whether the downloaded organism name matches the sheet name.
- `add_assembly_column.py --pathogens_sheet S --mapsheet M -o OUT [--download_input DI] [--position end]` adds an `assembly_accession` column to the sheet. When a taxid didn't download this run, its existing accession is kept.

---

## Workflow scripts that also run standalone (`bin/`)

### References and taxonomy

| Script                     | Does                                                                                                            | Needs     |
| -------------------------- | --------------------------------------------------------------------------------------------------------------- | --------- |
| `download_fastas.py`       | taxids / names / assembly accessions → one FASTA plus a gcfmapping, from an assembly summary with NCBI fallback | biopython |
| `fuzzy_match_assembly.py`  | local FASTA headers → assembly accession, by organism name against assembly summaries                           | biopython |
| `append_taxid.py`          | gcfmapping → `mapping.taxid.tsv` (`--custom-map`, `--ncbi-backup`)                                              | pandas    |
| `ncbi_taxid.py`            | `ncbi_taxid.py ACC [ACC ...]` → taxid via E-utilities                                                           | -         |
| `resolve_fasta_taxids.py`  | taxid of a local FASTA from its header accessions                                                               | -         |
| `extract_fasta_records.py` | keep only named records of a FASTA                                                                              | -         |
| `build_mask.py`            | region-mask BED from NCBI FCS reports, annotation or BEDs ([Region Masking](masking.md))                        | pysam     |

### Getting reads

| Script               | Does                                                  | Needs |
| -------------------- | ----------------------------------------------------- | ----- |
| `resolve_sra.py`     | SRA/ENA accessions → a downloadable FASTQ manifest    | -     |
| `fetch_ena_fastq.py` | download ENA FASTQs with md5 verification and retries | -     |

### Simulation and spike-ins ([In-Silico Simulation](in-silico.md))

| Script                          | Does                                                          | Needs |
| ------------------------------- | ------------------------------------------------------------- | ----- |
| `make_simulated_samples.py`     | abundance files + reference FASTAs for ISS / NanoSim          | pysam |
| `prepare_nanosim_inputs.py`     | NanoSim metagenome inputs from ISS-style abundances           | pysam |
| `parse_spikein_sheet.py`        | normalise a spike-in sheet                                    | -     |
| `subsample_insilico_reads.py`   | master synthetic FASTQ → spike-in / dilution-series datasets  | -     |
| `spike_into_background.py`      | draw exact read counts from simulated pools into a background | -     |
| `reconstruct_insilico_reads.py` | rebuild a subsampled FASTQ from a master FASTQ + read index   | -     |

### Results and reports

| Script                              | Does                                                                                                                        | Needs                |
| ----------------------------------- | --------------------------------------------------------------------------------------------------------------------------- | -------------------- |
| `export_data.py`                    | multi-tab export of report data (XLSX, CSV, wide or stacked tables), with filters such as `--min-tass` and `--passing-only` | openpyxl (XLSX)      |
| `combine_samples_json.py`           | merge per-sample `.paths.json` into `all.samples.json`                                                                      | pandas               |
| `make_report.py`                    | multi-run comparison HTML report                                                                                            | pandas               |
| `alignment_trends.py`               | cross-sample depth trends ([Alignment Trends](alignment-trends.md))                                                         | matplotlib, openpyxl |
| `mask_compare.py`, `mask_report.py` | TASS before vs. after masking; visual summary of a masked run                                                               | -                    |
| `get_top_hits.py`                   | top-hit organisms from a Kraken2 report                                                                                     | pandas               |
| `make_kt_files.py`                  | Kraken report → Krona input                                                                                                 | -                    |

---

## Developer scripts

The rest of `scripts/` supports development rather than analysis: `test_*.js` / `test_*.py` exercise the interactive report and metadata CLI, and `migrate_wiki.py`, `write_docs_ref.py`, `docs_release_stub.py` and `fetch_offline_report_libs.py` maintain this documentation site and the offline report bundle.
