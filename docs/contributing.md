# Contributing

Thank you for your interest in contributing to TaxTriage! We manage tasks and bugs through GitHub issues and welcome contributions of all kinds — bug reports, feature requests, documentation improvements, and code.

---

## Getting Help First

Before contributing, check the [Troubleshooting](troubleshooting.md) page. For questions about usage or development, reach out on:

- **Slack:** [#taxtriage on nf-core Slack](https://nfcore.slack.com/channels/taxtriage) ([join here](https://nf-co.re/join/slack))
- **GitHub Issues:** [github.com/jhuapl-bio/taxtriage/issues](https://github.com/jhuapl-bio/taxtriage/issues)

---

## Contribution Workflow

### 1. Check for Existing Issues

Search [existing issues](https://github.com/jhuapl-bio/taxtriage/issues) to avoid duplicating work. If no issue exists, create one so others know you're working on it.

### 2. Fork and Branch

[Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [jhuapl-bio/taxtriage repository](https://github.com/jhuapl-bio/taxtriage) to your GitHub account, then create a feature branch from `dev`.

### 3. Make Changes

Follow the [pipeline conventions](#pipeline-contribution-conventions) below. For new parameters, update the JSON schema:

```bash
nf-core schema build
```

This requires [nf-core tools](https://github.com/nf-core/tools) >= 1.10.

### 4. Submit a Pull Request

Open a pull request against the **`dev` branch**. Wait for CI tests to pass and code review to complete before merging.

---

## Automated Tests

Pull requests trigger two sets of automated tests via GitHub Actions:

### Lint Tests

`nf-core lint` checks the pipeline against nf-core guidelines. Run locally with:

```bash
nf-core lint <pipeline-directory>
```

Fix any warnings or failures listed before submitting your PR.

### Pipeline Tests

The pipeline is tested end-to-end on a minimal test dataset. Tests run against both the **latest** and **minimum required** versions of Nextflow. If tests fail, review the error messages in the GitHub Actions log.

---

## Pipeline Contribution Conventions

- All new processes should use containers from [Biocontainers](https://biocontainers.pro/) where possible
- New parameters should be documented in `nextflow_schema.json` (via `nf-core schema build`)
- Follow DSL2 conventions — one container per process
- Add test coverage for any new feature using the minimal test dataset

---

## Patch Releases (Bug Fixes)

In the rare event a release contains a bug:

1. Create a `patch` branch from `upstream/master` on your fork
2. Fix the bug and bump the patch version (`X.Y.Z+1`)
3. Open a PR directly against `master` from the `patch` branch

---

## Adding Organisms to the Pathogen Sheet

The curated pathogen sheet (`assets/pathogen_sheet.csv`) can be extended without code changes:

1. Add new rows to the CSV with at minimum: `name`, `taxid`, `general_classification`, `high_consequence`
2. Optionally add `pathogenic_sites` and `commensal_sites` columns for site-specific annotation. Prefer these over changing `general_classification` — they make the annotation conditional on body site instead of global. See [Microbial Categories](microbial-categories.md#8-curating-categories) for curation guidance.
3. Optionally set `assembly_accession` (`GCF_*` / `GCA_*`) to pin a specific validated genome for that organism — it is used **accession-first** at download time, falling back to taxid-based selection when blank (see [Assembly Selection Order](cli-parameters.md#assembly-selection-order)). This column is appended last and regenerated automatically when the database is rebuilt, so hand-edits to it may be overwritten.
4. Open a GitHub issue to request additions to the default sheet

---

## Resources

- [GitHub collaboration guide](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests)
- [Git learning resources](https://try.github.io/)
- [nf-core developer guidelines](https://nf-co.re/developers/guidelines)
- [nf-core/configs documentation](https://github.com/nf-core/configs#documentation)
