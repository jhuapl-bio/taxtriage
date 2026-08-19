# TaxTriage Documentation

[![DOI](https://img.shields.io/badge/doi-10.1093/bioinformatics/btag119-blue)](https://doi.org/10.1093/bioinformatics/btag119)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)

[![Docker](https://img.shields.io/badge/install%20run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![Singularity](https://img.shields.io/badge/install%20run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## See it first

[Interactive report demo](demo-report.md){ .md-button .md-button--primary }
[Example ODR (PDF)](images/pathogens.report.example.pdf){ .md-button }
[Pathogen sheet](pathogen-sheet.md){ .md-button }

The **demo** is a live report built from the example dataset — the same artifact
the pipeline writes at the end of a run. The **ODR** is the static PDF
deliverable. The **pathogen sheet** is every organism TaxTriage can flag,
searchable and filterable.

## About TaxTriage

TaxTriage is a flexible, containerized bioinformatics pipeline designed to identify pathogens within complex samples (e.g., respiratory swabs, lesion swabs, whole blood) using untargeted DNA or RNA sequencing data. It supports short-read (Illumina) and long-read (ONT, PacBio) platforms, incorporating quality control, organism classification, read mapping, and a unified confidence metric for all identified organisms.

The final analysis output is an **Organism Discovery Report (ODR)** — 2 files (a PDF and an interactive HTML) with summaries of all intermediate data supporting pathogen identification. TaxTriage is designed for broad deployment and early-stage outbreak investigations and is **not** intended as a standalone diagnostic capability.

---

## Quick Navigation

| Section                                         | Description                                                                        |
| ----------------------------------------------- | ---------------------------------------------------------------------------------- |
| [Installation](installation.md)                 | Install Nextflow, Docker, and Singularity                                          |
| [Quick Start](quick-start.md)                   | Run your first test pipeline in minutes                                            |
| [Samplesheet](samplesheet.md)                   | Format and prepare your input samplesheet                                          |
| [Running the Pipeline](running-the-pipeline.md) | Commands, profiles, and execution modes                                            |
| [CLI Parameters](cli-parameters.md)             | Full reference for all pipeline parameters                                         |
| [Pipeline Modules](pipeline-modules.md)         | Step-by-step breakdown of the workflow                                             |
| [Output](output.md)                             | Understanding output files and directories                                         |
| [TASS Scoring](tass-scoring.md)                 | How the confidence score is calculated                                             |
| [Microbial Categories](microbial-categories.md) | Primary, Opportunistic, Potential, Commensal — what they mean and how to read them |
| [In-Silico Simulation](in-silico.md)            | Simulated read validation                                                          |
| [Cloud & Seqera](cloud-and-seqera.md)           | Running on AWS with Nextflow Tower / Seqera                                        |
| [Geneious Plugin](geneious-plugin.md)           | Geneious Prime integration                                                         |
| [Troubleshooting](troubleshooting.md)           | FAQ and common errors                                                              |
| [Contributing](contributing.md)                 | How to contribute to TaxTriage                                                     |
| [Citations](citations.md)                       | Tools and publications to cite                                                     |

---

## Pipeline Overview

TaxTriage ingests raw FASTQ data and processes it through the following major stages:

1. **Read QC** — FastQC / NanoPlot / pycoQC
2. **Trimming** — Trimgalore (Illumina) / Porechop (ONT)
3. **Host Removal** — Minimap2 against host reference
4. **Metagenomics Classification** — Kraken2 (+ Krona plots)
5. **Top Hits Assignment** — Selects organisms for alignment
6. **Reference Preparation** — Downloads assemblies from NCBI
7. **Alignment** — BWA-MEM2 (Illumina) / Minimap2 (ONT/PacBio)
8. **Stats** — Coverage, depth, MAPQ via samtools
9. **TASS Scoring** — Confidence metric calculation
10. **Report Generation** — MultiQC + Organism Discovery Report PDF

![TaxTriage Schematic](images/taxtriage_schematics.png)

---

## Collaborative Efforts

### UW Madison — Geneious Prime Plugin

Dave O'Connor's Laboratory at the University of Wisconsin-Madison developed a custom Geneious Prime plugin to run TaxTriage analyses directly from within Geneious. See the [Geneious Plugin](geneious-plugin.md) page for details.

---

## Citation

If you use TaxTriage, please cite:

> Brian B. Merritt & Jeremy D. Ratcliff et. al, TaxTriage: an open-source metagenomic sequencing data analysis pipeline enabling putative pathogen detection. doi: [10.1093/bioinformatics/btag119/8571885](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btag119/8571885)

See the [Citations](citations.md) page for the full list of tools to cite.

---

## Copyright

Copyright 2022–2026 The Johns Hopkins University Applied Physics Laboratory LLC. See [LICENSE](https://github.com/jhuapl-bio/taxtriage/blob/main/LICENSE) for details.
