## Troubleshooting and FAQ

WORK IN PROGRESS

### Input

#### What is a Samplesheet?

A Samplesheet is a .csv (comman delimited) file that contains the necessary locations of all of your input data. It consists of a set of columns defined at this [location] (https://github.com/jhuapl-bio/taxtriage/blob/main/docs/usage.md#samplesheet-input)
This is required and denoted by the parameter `--input` before it if you're running from the command line like `nextflow run https://github.com/jhuapl-bio/taxtriage --input <name of samplesheet...`

Each Row corresponds to a single sample (the first column in the example samplesheet [here](https://github.com/jhuapl-bio/taxtriage/blob/main/examples/Samplesheet_cli.csv).
You can use this example samplesheet as a reference in making your own samplesheet. Each row MUST have a sample name as it is used to track all downstream modules/steps that are used

See this [table](https://github.com/jhuapl-bio/taxtriage/blob/main/docs/usage.md#samplesheet-information) for a description of all columns

One important thing to note is that the `fastq_1` column MUST have a value if using a FILE that is a compressed fastq (`.gz` extension), `from` is directory of fastq files that may or may not be compressed with `.gz`. A paired end read set for Illumina must have a file in `fastq_1` and another in `fastq_2` for each row.
IF you're using `from` you must have a directory of fastq file(s) in the path

### Nextflow Logs of TaxTriage

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
It is important to note that `Command error: [CRITICAL]
` line indicates what went wrong. In this case, the second line in the samplesheet is incorrect as it needs a FASTQ File `At least the first FASTQ file is required. On line 2.`

I can simply edit that column and rerun with the `-resume` flag on the CLI or by clicking the dots in the upper right of the failed job and hitting "resume" once I've made my changes
to the samplesheet on S3 on AWS.

#### ERROR - my pipeline fails at NanoPlot?

Check to make sure that the directory/file being used is not empty. If no data is present NanoPlot will fail

#### ERROR - I cant get the Assemblies to download

There is a good chance that you are behind a firewall or require specifical certificates to be setup, please try to adjust things to work in your organization's internet as this step uses `curl` to pull data

#### What is the difference between remove_taxids and filter when getting rid of host reads?

`filter` is when you want to use a human reads database (or any other database) that automatically removes any and all classified reads. `remove_taxids` is a list of taxids, separated by a space, that will be removed post-kraken2 output with your `--db` database.
Please be aware that these options will not catch all organisms of interest, so please be aware of data privacy concerns before passing read data off that have gone through these steps

### Why is my run taking so long?

This could be a variety of reasons. Specific modules are prone to take much longer than others, namely the MultiQC (final step), consensus/assembly, or any of the alignment steps. Additionally, plotting for Oxford Nanopore can take a large amount of time when using NanoPlot. Consider disabling many of these features if you simply want the baseline kraken results with --skip\_{name_of_step(s)}.

### Common errors seen

For the most part, many of the steps do not have a mandatory success error code attached to them. That is, a sample can be incomplete or have errors at a certain step, but the pipeline can continue. A few pipelines **do** require successful completion, however.

Here is a list of steps that are required to complete successfully as well as probably solutions to your problem.

### SAMPLESHEET_CHECK

### FASTP

**Description**: This is the QC step that removes low quality reads from your samples.

**Possible Problems & Fixes**

FASTP halts with a red error code:

```
ERROR ~ Error executing process > 'NFCORE_TAXTRIAGE:TAXTRIAGE:FASTP (shortreads)'

Caused by:
  Process `NFCORE_TAXTRIAGE:TAXTRIAGE:FASTP (shortreads)` terminated with an error exit status (255)

Command executed:

  [ ! -f  shortreads_1.fastq.gz ] && ln -sf iss_reads_R1.fastq.gz shortreads_1.fastq.gz
  [ ! -f  shortreads_2.fastq.gz ] && ln -sf iss_reads_R2.fastq.gz shortreads_2.fastq.gz
  fastp \
      --in1 shortreads_1.fastq.gz \
      --in2 shortreads_2.fastq.gz \
      --out1 shortreads_1.fastp.fastq.gz \
      --out2 shortreads_2.fastp.fastq.gz \
      --json shortreads.fastp.json \
      --html shortreads.fastp.html \
       \
       \
       \
      --thread 2 \
      --detect_adapter_for_pe \
        -q 500  \
      2> shortreads.fastp.log

  cat <<-END_VERSIONS > versions.yml
  "NFCORE_TAXTRIAGE:TAXTRIAGE:FASTP":
      fastp: $(fastp --version 2>&1 | sed -e "s/fastp //g")
  END_VERSIONS

Command exit status:
  255

Command output:
  (empty)

Work dir:
  /Users/merribb1/Documents/Projects/APHL/taxtriage/work/b5/6e19a5de67b703d5f44cbe1aed94fa

Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

 -- Check '.nextflow.log' file for details
```

The above error is the result of either:
a. The quality of the reads is lower than the default. For ONT, it is 7 and for Illumina it is 20. - Solution: Update the `--minq` parameter to be lower or specify `--skip_fastp`
b. The reads themselves are empty. Make sure you are specifying the correct path to the files in your samplesheet.

### KRAKEN2_KRAKEN2

**Description**: This the metagenomics step "agnostically" search for all organisms in your samples files in a quick manner. It is tied into the "top hits" and pathogens discovery parameters to downsamples the downstream alignments that are needed.

**Possible Problems & Fixes**

- Memory error: Oftentime, Kraken2 requires loading the **entire** db into RAM by default. Ensure that your RAM limits are greater than the size of the db directory. You can check with `du -sh $DB` where $DB is the path of the database.
  a. Set `--low_memory` to read the DB in by I/O.
  b. Set `--max_memory` if your memory is lower than the default required (~36 GB). For example: `--max_memory 13GB`
- Invalild or corrupt DB:
  a. Check that the $DB directory is not corrupt. Usually, this has 3 `k2d` files called:
  1. `hash.k2d`
  2. `opts.k2d`
  3. `taxo.k2d`

### Bowtie2 (Index and Align) / Minimap2 (Align)

- Incorrect or Unknown Ref. FASTA file
  a. If using Kraken2 (default, no `--skip_kraken2` specified), ensure that you have Internet connectivity to pull the FASTA files. You can check the `download` folder and look for the `download/<samplename>.output.references.fasta`. If it is empty or not-present, then the sample failed to get the taxids.
  b. Check that you have 1 or more tophits. You can see that in `tops/<samplename>.top_report.tsv` that there is more than one line
  c. Out of Memory - This can be altered by raising the Memory in `Docker` (if using) or by limiting the number of top hits if your computer can't handle all of the FASTA file information.
  - Minimap2: Loads the FASTA into RAM as an index so ensure that the FASTA **File** is under the limits on your RAM. The pipeline will re-attempt with more tries if the RAM is too low by raising it up to 3 consecutive times.

```
Command executed:

  minimap2 \
        -ax map-ont  \
      -t 2 -I 0G \
      longreads.dwnld.references.fasta \
      longreads.classified.fastq.gz \
       \
      -L \
      -a | samtools sort | samtools view  -q 5  -@ 2 -b -h -o longreads.bam  -q 5


  cat <<-END_VERSIONS > versions.yml
  "NFCORE_TAXTRIAGE:TAXTRIAGE:ALIGNMENT:MINIMAP2_ALIGN":
      minimap2: $(minimap2 --version 2>&1)
  END_VERSIONS

Command exit status:
  1

Command output:
  (empty)

Command error:
  [M::mm_idx_gen::0.065*1.19] collected minimizers
  [M::mm_idx_gen::0.102*1.09] sorted minimizers
  [WARNING] For a multi-part index, no @SQ lines will be outputted. Please use --split-prefix.
  [M::main::0.105*1.07] loaded/built the index for 1 target sequence(s)
  [M::mm_mapopt_update::0.108*1.06] mid_occ = 10
  [M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 1
  [M::mm_idx_stat::0.110*1.06] distinct minimizers: 87097 (99.73% are singletons); average occurrences: 1.005; average spacing: 5.352; total length: 468348
  [E::sam_parse1] no SQ lines present in the header
  [W::sam_read1_sam] Parse error at line 4
  samtools sort: truncated file. Aborting
  [main_samview] fail to read the header from "-".
```

- Bowtie2: Most of the is from the `build` portion of the pipeline.

### SAMTOOLS_HIST_COVERAGE

This step is very low level. Make sure that you have a `bam` file in the `minimap2/bowtie2` directories for each sample.

### SpAdes / Flye

This step is not recommended but is useful for generating either contigs or complete chromosomes from your data. The biggest concern for these is memory allocation required. The best options for this are to

1. Remove the `--reference_fasta` if specified and your FASTA is large
2. Remove or lower the `--top_hits` or `--top_per_taxa` parameters to try to assembly with less reads

### SAMTOOLS_SORT

This step is usually CPU intensive rather than with RAM. However, you still need a sufficient amount of RAM so ensure that (above descriptions) that you've set memory limits as needed.

Or, ensure that you have a `bam` file that is **not** empty!

### SAMTOOLS_VIEW

The only concern with this is either storage capacity (HDD/SSD) or timeouts. This is a low-compute-requiring command to run.
