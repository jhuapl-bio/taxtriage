## Troubleshooting and FAQ
WORK IN PROGRESS

### 1. Input

#### What is a Samplesheet?

A Samplesheet is a .csv (comman delimited) file that contains the necessary locations of all of your input data. It consists of a set of columns defined at this [location] (https://github.com/jhuapl-bio/taxtriage/blob/main/docs/usage.md#samplesheet-input)
This is required and denoted by the parameter `--input` before it if you're running from the command line like `nextflow run https://github.com/jhuapl-bio/taxtriage --input <name of samplesheet...`

Each Row corresponds to a single sample (the first column in the example samplesheet [here](https://github.com/jhuapl-bio/taxtriage/blob/main/examples/Samplesheet_cli.csv).
You can use this example samplesheet as a reference in making your own samplesheet. Each row MUST have a sample name as it is used to track all downstream modules/steps that are used

See this [table](https://github.com/jhuapl-bio/taxtriage/blob/main/docs/usage.md#samplesheet-information) for a description of all columns

One important thing to note is that the `fastq_1` column MUST have a value if using a FILE that is a compressed fastq (`.gz` extension), `from` is directory of fastq files that may or may not be compressed with `.gz`. A paired end read set for Illumina must have a file in `fastq_1` and another in `fastq_2` for each row.  
IF you're using `from` you must have a directory of fastq file(s) in the path

### 2. Nextflow Tower
### 3. AWS (Compute env and S3)
### 4. Installation
### 5. Basestack deployment
### 6. Nextflow Logs of TaxTriage

ERROR - RED Text is all over the place, did my pipeline fail?

IF you're running TaxTriage in the Command line or viewing the Execution logs in Nextflow Tower runs, you may experience an error like: 

```
#### ERROR ~ Error executing process > 'NFCORE_TAXTRIAGE:TAXTRIAGE:INPUT_CHECK:SAMPLESHEET_CHECK (Samplesheet_cli.csv)'

Caused by:
  Process `NFCORE_TAXTRIAGE:TAXTRIAGE:INPUT_CHECK:SAMPLESHEET_CHECK (Samplesheet_cli.csv)` terminated with an error exit status (1)

Command executed:

  check_samplesheet.py \
      Samplesheet_cli.csv \
      samplesheet.valid.csv
  
  cat <<-END_VERSIONS > versions.yml
  "NFCORE_TAXTRIAGE:TAXTRIAGE:INPUT_CHECK:SAMPLESHEET_CHECK":
      python: $(python --version | sed 's/Python //g')
  END_VERSIONS

Command exit status:
  1

Command output:
  (empty)

Command error:
  [CRITICAL] At least the first FASTQ file is required. On line 2.

Work dir:
  /Users/merribb1/Documents/Projects/APHL/taxtriage/work/b3/cb238415397372a2caa35c06c77c51

Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

 -- Check '.nextflow.log' file for details
```

The title of the module, in this example, is the INPUT_CHECK:SAMPLESHEET_CHECK where we are importing the information from a samplesheet. This name may be different depending on the step in the pipeline
It is important to note that ```Command error: [CRITICAL] 
``` line indicates what went wrong. In this case, the second line in the samplesheet is incorrect as it needs a FASTQ File `At least the first FASTQ file is required. On line 2.`

I can simply edit that column and rerun with the `-resume` flag on the CLI or by clicking the dots in the upper right of the failed job and hitting "resume" once I've made my changes
to the samplesheet on S3 on AWS.


#### ERROR - my pipeline fails at NanoPlot?

Check to make sure that the directory/file being used is not empty. If no data is present NanoPlot will fail


#### ERROR - I cant get the Assemblies to download

There is a good chance that you are behind a firewall or require specifical certificates to be setup, please try to adjust things to work in your organization's internet as this step uses `curl` to pull data


#### What is the difference between remove_taxids and filter when getting rid of host reads?

`filter` is when you want to use a human reads database (or any other database) that automatically removes any and all classified reads. `remove_taxids` is a list of taxids, separated by a space, that will be removed post-kraken2 output with your `--db` database. 
Please be aware that these options will not catch all organisms of interest, so please be aware of data privacy concerns before passing read data off that have gone through these steps


   
