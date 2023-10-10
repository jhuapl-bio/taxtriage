# nf-core/taxtriage: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/taxtriage/usage](https://nf-co.re/taxtriage/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```console
--input '[path to samplesheet file]'
```

### Multiple Samples 

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,platform,fastq_1,fastq_2,sequencing_summary,trim
NB03,OXFORD,examples/data/fastq_demux/NB03,,,FALSE
BC05_flu,OXFORD,examples/data/BC05.fastq.gz,,,FALSE
longreads,OXFORD,examples/data/nanosim_metagenome.fastq.gz,,,FALSE
shortreads,ILLUMINA,examples/data/iss_reads_R1.fastq.gz,examples/data/iss_reads_R2.fastq.gz,,TRUE
```




### Multiple Samples AND Platforms 

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:



```console
sample,platform,fastq_1,fastq_2,sequencing_summary,trim
NB03,OXFORD,examples/data/fastq_demux/NB03,,,FALSE
BC05_flu,OXFORD,examples/data/BC05.fastq.gz,,,FALSE
longreads,OXFORD,examples/data/nanosim_metagenome.fastq.gz,,,FALSE
shortreads,ILLUMINA,examples/data/iss_reads_R1.fastq.gz,examples/data/iss_reads_R2.fastq.gz,,TRUE
```

### Samplesheet Information

| Column     | Description                                                                                                                                                                            |
| ---------  | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`   | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`  | OPTIONAL (if not using 'from'). Full path to FastQ file for Illumina short reads 1 OR OXFORD reads. File MUST be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                             |
| `fastq_2`  | OPTIONAL. Full path to FastQ file for Illumina short reads 2. File MUST be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `platform` | Platform used, [ILLUMINA, OXFORD]                                                            |
| `trim` | TRUE/FALSE, do you want to run trimming on the sample?                                       |
| `sequencing_summary` | OPTIONAL. If detected, output plots based on the the sequencing summary file for that sample | 

An [example samplesheet](../examples/Samplesheet.csv) has been provided with the pipeline alongside some demo data.

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/taxtriage --input samplesheet.csv --outdir ./$OUTDIR -profile docker
```

:warning: please be aware that for local instances, due to limited resources, we have defined several config profiles. Please try `-profile test,docker` for a trial run of dummy data and then use `-profile local,docker` first on your own data before manually adjusting parameters

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work                # Directory containing the nextflow working files
<OUTIDR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### CLI Parameters Possible and Explained

| Parameter     | Description                                                                                                                                                                            |
| ---------  | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--input <samplesheet.csv_path>`   | CSV file containing samplesheet information for each sample. See [Samplesheet section](#samplesheet-information)  |
| `--outdir <path_to_output_to_be_made>`  | Output folder.                                             |
| `--skip_assembly`  | Skip the assembly process. Currently this is recommmended as spades/flye are not functioning properly |
| `--skip_fastp`  | TRUE/FALSE,  do not filter with fastp |
| `--skip_plots`  | TRUE/FALSE,  do not make any plots  |
| `--remove_taxids  <numbers[]>`  | "taxidA taxidB..." a list of one or more taxids to remove from the kraken report prior to downstream analysis. Use "'9606'" for human reads |
| `--k2_confidence  <numbers[]>`  | Minimum confidence to classify a read using KRAKEN2 only. See [here](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#confidence-scoring) for more information. |
| `--skip_realignment`  | TRUE/FALSE, Skip realignment step. You will not get a metrics report as a result                    |            
| `--minq <number>`     | What minimum quality would you want in your samples. Disabled if you run --skip_fastp. Default is 7 for Oxford Nanopore, 20 for Illumina         |
| `--trim`     | Remove adapters from data prior to qc filtering. Trimgalore for Illumina, Porechop for Illumina (SLOW)  |
| `--subsample <number>`     | Take a subsample of n reads from each sample. Useful if your data size is very large and you want a quick triage analysis      |
| `--reference_fasta <filepath>`     | Location of a reference fasta to run alignment on RATHER than downloading references from NCBI. Useful if you know what you're looking for or have no internet connection to pull said references     |
| `--db <path_to_kraken2_database>` | Database to be used. IF `--low_memory` is called it will read the database from the fileystem. If not called, it will load it all into memory first so ensure that the memory available (limited as well by `--max_memory` is enough to hold the database). If using with --download-db, choose from download options {minikraken2, flukraken2} instead of using a path. [See here for a full list](#supported-default-databases)  |
| `--download_db` | Download the preset database indicated in `--db` to `--outdir` |
| `--max_memory <number>GB` | Max RAM you want to dedicate to the pipeline. Most of it will be used during Kraken2 steps so ensure you have more memory than the size of the `--db` called |
| `-latest` | Use the latest revision (-r). If you want to autopull the latest commit for a branch from https://github.com/jhuapl-bio/taxtriage, specify this. Used in Basestack and the Cloud (default toggle) |
| `--low_memory <number>` | If you don't have enough memory to load the kraken2 db, call this to read it from the filesystem for each read. THIS IS MUCH SLOWER THAN THE DEFAULT METHOD |
| `--max_cpus <number>` | Max CPUs you want to dedicate to the pipeline | 
| `--top_per_taxa ,taxid:amount:rank]` | One or more 3 element values for the minimum taxa at a rank. Example "10239:20:S 2:20:S" is minimum 20 virus species (10239 is Virus) AND (separated by space for a new definition) 20 Bacteria (2) species. Possible Rank codes are G,S,P,F,O,C | 
| `--demux` | If your Samplesheet contains a folder (rather than 1-2 fastq files), you MUST call this flag | 
| `--top_hits_count` | Minimum taxa to filter out per rank level. If top_per_taxa specified, whatever is the larger of the 2 values for the results takes priority. | 
| `-resume` | Resume the run from where it left off. IF not called the pipeline will restart from the Samplesheet check each time | 
| `-r [main, stable, etc.]` | Specify the branch/revision name to use if pulling from github (not local main.nf file) | 
| `-profile [local,test,test_viral,docker,singularity]` | Default profile, 2 tests, Docker, or Singularity for execution reasons | 

### Supported Default databases

                                                                                               
| Database Name | Location | Size | 
| ------  | ------ | --- |
| standard8  | [Download](https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz) | 7.5G |
| viral  | [Download](https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20230605.tar.gz) | 553M |
| flukraken2  | [Download](https://media.githubusercontent.com/media/jhuapl-bio/mytax/master/databases/flukraken2.tar.gz) | 180M |
| test    |    [Download](https://github.com/jhuapl-bio/datasets/raw/main/databases/kraken2/test_metagenome.tar.gz)   |    112M   |
| minikraken2  | [Download](https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v2_8GB_201904.tgz) | 7.5G |

## AWS with Nextflow Tower

While we support Taxtriage in both Basestack and native and local CLI deployment, you can also import and run code from Nextflow Tower. This process can be convoluted in order to get it applicable for cloud environments but, once fully setup, becauses very easy to reproduce at low cost on AWS. Please be aware of sensitivity of data when working in the cloud

1. First, you must create a [Nextflow Tower](https://cloud.tower.nf/) account
    - Further documentation for the following steps can be found [here](https://help.tower.nf/22.2/)
2. Request a user account to brian.merritt@jhuapl.edu
    - Information for accessing the S3 Buckets and creating a credential-specific Compute environment will be included in the email
3. At this point follow all steps for setting up AWS in the following link. [View Steps Here](images/Cloud_AWS_NFTower_only/Cloud_AWS_NFTower_only.pdf)


#### Working with your own data NF Tower

When all necessary items have been setup, you'll need to load your data you want to run into an S3 bucket. Be aware that all filesystem refrences to files/directories need to be relative to the S3 bucket. However, with things like the Samplesheet, the paths can be relative only within the S3 bucket so you don't need to use things like `S3://bucketname/path/to/file `

<img src="images/s3bucket.png" alt="drawing" height="300"/>
<img src="images/s3bucket2.png" alt="drawing" height="300"/>

In the above example, I've made an s3 bucket with necessary permissions for nextflow-tower to run. I have a directory with the same test-data you can find in the `examples` directory in the root level of this bit of code [here](../examples)


Once the AWS system is setup, let's head back to Nextflow Tower. On the left, you can select a drop-down and open up your own personal launchpad (and other features)

<img src="images/cloud_showcase_nftower.png" alt="drawing" width="40%" height="20%"/>

<img src="images/addpipeline.png"  width="20%" height="30%">


<img src="images/taxtriagelaunchpad1.png"  width="50%" height="50%">


MAKE SURE that the compute environment matches the one you set up when you set your credentials in AWS with NF tower

If you expand the pipeline parameters, you can mimic what I've written for my example with your own paths for the S3 bucket and example data. Note that these are going to be identical to the parameters available at [here](#cli-parameters-possible-and-explained)


<img src="images/taxtriagelaunchpad2.png"  width="50%" height="50%">


Next, Run the pipeline from the Launchpad. Be aware that all inputs to the S3 bucket are tailored for my example. You own inputs will vary including things like `low_memory` or database name used

[cloud_run_2.webm](https://user-images.githubusercontent.com/50592701/192596313-7e30f285-dc1d-4c62-99d2-5791a5d8c0e9.webm)

[cloud_run_3.webm](https://user-images.githubusercontent.com/50592701/192596272-46007980-cc07-46c3-978f-e1846adbfffb.webm)

[cloud_run_1.webm](https://user-images.githubusercontent.com/50592701/192596324-57162d50-2738-4b7f-ba8c-2fa473ca1433.webm)

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/taxtriage releases page](https://github.com/nf-core/taxtriage/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)


### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

## Custom configuration


### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
