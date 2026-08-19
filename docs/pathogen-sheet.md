---
title: Pathogen Sheet
hide:
  - toc
---

<style>.md-content__inner{max-width:none}</style>

# Pathogen Annotation Sheet

Every organism TaxTriage can flag, with the annotations that drive
[Microbial Categories](microbial-categories.md) and the
[TASS confidence score](tass-scoring.md).

This table reads
[`assets/pathogen_sheet.csv`](https://github.com/jhuapl-bio/taxtriage/blob/main/assets/pathogen_sheet.csv)
straight from the `main` branch of the pipeline repository **when the page
loads** — it is not bundled into this site. Merge a change to that CSV and it
shows up here on the next refresh, with no redeploy. Columns the sheet does not
currently carry are left blank.

!!! tip "Using this table"

    Search matches organism names, synonyms, tax IDs and full lineage. Facets
    combine with AND across groups and OR within a group, and each group's counts
    reflect the *other* active filters. Click any row for the full record
    including references. **Export CSV** downloads exactly what you have filtered
    to — handy for building a custom `--pathogens` sheet.

<div class="pt-app" id="pathogen-app">
  <div class="pt-empty">Loading the pathogen sheet…</div>
</div>

## Column reference

| Column                   | Meaning                                                                                                                  |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------ |
| `name`                   | Organism name, matched against the taxonomic classification output                                                       |
| `taxid`                  | NCBI Taxonomy identifier                                                                                                 |
| `general_classification` | `primary`, `opportunistic`, `potential`, or `commensal` — see [Microbial Categories](microbial-categories.md)            |
| `status`                 | `established` where the pathogen–site association is well documented, `putative` where it is suggested but not confirmed |
| `high_consequence`       | Flags organisms warranting immediate escalation regardless of abundance                                                  |
| `pathogenic_sites`       | Body sites where detection is considered clinically meaningful                                                           |
| `commensal_sites`        | Body sites where the organism is expected flora, which down-weights its score                                            |
| `alternative_names`      | Synonyms and former names, also searched during matching                                                                 |
| `pathology`              | Free-text disease association                                                                                            |
| `host_organism`          | Host the annotation applies to                                                                                           |
| `kingdom` … `genus`      | Lineage, used for genus-level rollups                                                                                    |
| `mol_type`               | `dna` or `rna`, which determines the alignment path                                                                      |
| `reference`              | Literature supporting the annotation                                                                                     |
| `assembly_accession`     | Default reference assembly downloaded for alignment                                                                      |

## Contributing an organism

The sheet is a plain CSV — additions are pull requests against the pipeline repo.
See [Contributing](contributing.md#adding-organisms-to-the-pathogen-sheet) for
the required fields and review expectations.

To run with your own sheet instead of the bundled one:

```bash
nextflow run jhuapl-bio/taxtriage \
    --input samplesheet.csv \
    --pathogens /path/to/my_pathogen_sheet.csv \
    -profile docker
```
