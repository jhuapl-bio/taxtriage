[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/taxtriage/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularityCE](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23taxtriage-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/taxtriage)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

![](assets/taxtriage_schematics.png)

#### Alerts

:warning: Git History was recently updated, causing a conflict in updating already cloned repos when running the test profile or called `-latest -r main/stable`. As a result you must run `nextflow drop https://github.com/jhuapl-bio/taxtriage` first. This only applies to pipelines run by calling the remote repo and the previously mentioned parameters

**TaxTriage** is a bioinformatics best-practice analysis pipeline for APHL pipeline for triage classification reports.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/taxtriage/results).

Tax Triage is designed as a pipeline for the purpose of giving an initial triage of taxonomic classifications, using Kraken2 database(s), that can then be ingested into a CLIA-style report format. It is under active development, but in the current state it is capable of running a set number of samples end-to-end using a user-created samplesheet in `.csv` format. The output format is a `HTML` which is highly interactive and distributable. This pipeline uses the `nextflow` ecosystem and is also available as a module in [Basestack](https://github.com/jhuapl-bio/Basestack). Currently, Basestack is undergoing improvements to allow easier usage of nextflow pipelines (includes TaxTriage) that is scheduled for release in early August.

Efforts are underway to provide full support of this pipeline on [nf-core](nf-core.re) to provide a seamless deployment methodology. The pipeline also requires installation of [Docker](https://docker.com) or [Singularity](https://docs.sylabs.io/) (*CE ONLY*) for the individual modules within it. Because these modules are separate from the source code of TaxTriage, we recommend following the examples outlined in the [usage details](docs/usage.md) first to automatically run the pipeline and install all dependencies while also giving you some example outputs and a better feel for how the pipeline operates.

[See Here for full usage details](docs/usage.md)

[See Here for troubleshooting & FAQ](docs/troubleshooting.md)

## Installation

TaxTriage requires 2 primary installs for it to work

1. Nextflow
2. Singularity or Docker (recommended)

### 1. Nextflow 

Follow instructions [here](https://nf-co.re/docs/usage/installation) or run these commands in your WSL2, Native Linux, or Mac environment


```
# Make sure that Java v11+ is installed:
java -version
 
# Install Nextflow
curl -fsSL get.nextflow.io | bash

```

Note, this command requires sudo to move to your home path. If you are on an HPC, make sure that nextflow is in your $PATH if not globally available
Place it in your `$PATH`

```
# Add Nextflow binary to your user's PATH:
mv nextflow ~/bin/
```

If installing globally, requiring sudo, type: 

```
sudo mv nextflow /usr/local/bin
```

When complete, verify installation with `nextflow -v` to see the version

### 2. Containerization Approach Install

Choose *A* (Recommended - Docker) or *B*. If on a HPC, talk with your IT to get B. Singularity setup. You do NOT need to install both software tools.

#### A. Docker 

Follow these steps for your OS [here](https://docs.docker.com/engine/install/) - IF on WSL2 (Windows), choose Docker Desktop for Windows and it should be available automatically in your WSL environment

#### B. Singularity

[Install Instructions](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)


## Quick Start

Make sure you have either Docker or Singularity installed, as well as Nextflow

### Test Run

This will pull the test data and run the pipeline. It should take ~10-15 minutes.

```
nextflow run https://github.com/jhuapl-bio/taxtriage -r main -latest -profile test,docker --outdir test_output -resume
```

❗If you want singularity instead, make sure to specify that in the profile instead of docker like: `test,singularity`

### Cloud

Follow the steps [here](docs/usage.md#aws-with-nextflow-tower)

### Local Data

1. Run the command:

```
nextflow drop -f https://github.com/jhuapl-bio/taxtriage
nextflow pull https://github.com/jhuapl-bio/taxtriage
cp -r ~/.nextflow/assets/jhuapl-bio/taxtriage ~/taxtriage
cd ~/taxtriage
```

:warning: If you get an error on uncommitted changes please run the `nextflow drop -f https://github.com/jhuapl-bio/taxtriage`

2. Running Kraken2 and FASTQC report with the k2_viral db

### Running it with the local config (for laptops/workstations) with limited RAM

```
nextflow run https://github.com/jhuapl-bio/taxtriage  \
  --outdir tmp_viral \
  -resume \
  --input examples/Samplesheet.csv \
  --taxtab "default" -r main -latest \
  --db "viral" --download-db \
  -profile local,docker
```

### Running it by overriding some parameters from the local config

```
nextflow run https://github.com/jhuapl-bio/taxtriage \
   --input examples/Samplesheet.csv -r main -latest \
   --db viral --download_db --skip_assembly \
   --outdir tmp --max_memory 10GB --max_cpus 3   \
   -profile docker  -resume --demux --remove_taxids "'9606'"
```

:warning: Please see the contents of test or local config to figure out what the defaults are for those profiles

Remember, if you are doing a single taxid, wrap it with '' inside the "" quote

#### Using a db that is on your local filesystem

```

nextflow run https://github.com/jhuapl-bio/taxtriage \
   --input examples/Samplesheet.csv \
   --db "k2_viral" -r main -latest \
   --outdir tmp_viral  \
   --profile local,docker \
   -resume
```

### Running it without internet availability

This will use a local assembly text and reference fasta, assuming the reference FASTA is called `refer.fasta`

You will need 3 files locally on your system

1. assembly
2. reference_fasta
3. db

```

nextflow run https://github.com/jhuapl-bio/taxtriage \
   --input examples/Samplesheet.csv \
   --db "k2_viral" -r main -latest \
   --outdir tmp --reference_fasta ./refer.fasta \
   -profile local,docker \
   -resume \
   --demux \
   --assembly examples/assembly_summary_refseq.txt

```

#### Using a Custom Taxonomy

```

nextflow run https://github.com/jhuapl-bio/taxtriage \
   --input examples/Samplesheet_flu.csv \
   --db viral --download_db -r main -latest \
   --outdir tmp_viral \
   -profile local,docker  \
   --assembly data/databases/flukraken2/library/influenza-fixed.fna --assembly_file_type kraken2 \
   -resume

```

#### Running on local nf files (test config)

:warning: Make sure you're in the `jhuaplbio/taxtriage` repo first!

```
nextflow run ./main.nf \
   -profile test,docker \
   --outdir tmp_local \
   --input examples/Samplesheet.csv \
   --db ~/Desktop/mytax/test_metagenome \
   -resume
```

If you want to download the databases from scratch, you can see them here
Make sure to Download these databases to your `Desktop` or wherever you are the most comfortable. Remember the location and specify the `--db` parameter as the absolute path. For example `~/Desktop/flukraken2`. Also, remove the `--download-db` parameter

- [standard-8](https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz)
- [viral](https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20230605.tar.gz)
- [flukraken2](https://media.githubusercontent.com/media/jhuapl-bio/mytax/master/databases/flukraken2.tar.gz)

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)).

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run https://github.com/jhuapl-bio/taxtriage -profile test,docker --outdir ./outdir
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker` or `singularity` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.

4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```console
           nextflow run https://github.com/jhuapl-bio/taxtriage -r main -latest --outdir test_output -profile <local,docker/singularity>
   ```

## Modules

0. Subsample (OPTIONAL)
1. Guppyplex (Oxford Nanopore Only)
2. QC Plotting part 1 (pycoQC – Oxford Nanopore)
3. Trimming (Trimgalore – Illumina, Porechop – Oxford Nanopore)
4. Filtering ( Kraken2 – Illlumina, Oxford Nanopore)
5. QC Plotting part 2 (FastQC – Illlumina, Nanoplot – Oxford Nanopore)
6. Classification ( Kraken2 – Illumina, Oxford Nanopore, Krona Plots)
7. Alignment for Stats ( BWAMEM2 – Illumina, Minimap2 – Oxford Nanopore)

- :warning:Currently, the only realignment is going to be based on a taxid call. For example, if there will not be a complete realignment of "order" despite there being multiple species all within that order. For the most part, this is limited to more specific ranks like species, strain, subspecies etc.

8. Report Generation ( MultiQC – Illumina, Oxford Nanopore)


## Credits

TaxTriage was originally written by Brian Merritt, MS Bioinformatics.

We thank the following people for their extensive assistance in the development of this pipeline:


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  taxtriage for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

## Copyright

##############################################################################################
Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
All rights reserved.
Permission is hereby granted, free of charge, to any person obtaining a copy of this
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify,
merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
OR OTHER DEALINGS IN THE SOFTWARE.

## Acknowledgements

##############################################################################################

This software tool was supported by the Cooperative Agreement Number NU60OE000104, funded by the Centers
for Disease Control and Prevention through the Association of Public Health Laboratories. Its contents
are solely the responsibility of the authors and do not necessarily represent the official views of the Centers
for Disease Control and Prevention, the Department of Health and Human Services, or the Association of Public Health
Laboratories.
