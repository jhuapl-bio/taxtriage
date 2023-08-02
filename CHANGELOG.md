# nf-core/taxtriage: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of nf-core/taxtriage, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

check_samplesheet.py would compare all the suffixes after the initial "." which would result in an error if samples have discrepancies after the initial ".". Solutoin was to compare only the last 2 suffixes that always should be ".fastq" and ".gz".

### `Dependencies`

### `Deprecated`
