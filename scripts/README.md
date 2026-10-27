# scripts/

Scripts that run **outside** the Nextflow workflow. Full usage: [Standalone Scripts](https://jhuapl-bio.github.io/taxtriage/standalone-scripts/) (docs) / [wiki](https://github.com/jhuapl-bio/taxtriage/wiki/Standalone-Scripts).

## Reference building

| Script | Purpose |
| --- | --- |
| `pull_contigs.py` | List of contig accessions → `--reference_fasta` + `--custom_accession_map` (local FASTA / NCBI, taxdump-normalised taxids, resumable). Standard library only. |
| `pathogen_db/build_pathogen_db.sh` | Rebuild `pathogens.fasta(.gz)` + `pathogens.mapping.taxid.tsv` from `assets/pathogen_sheet.csv` (the reference `--get_pathogens` downloads). |
| `pathogen_db/validate_names.py` | Sheet name vs. downloaded organism name per taxid. |
| `pathogen_db/add_assembly_column.py` | Add `assembly_accession` to the pathogen sheet; write the next build's download input. |

```bash
python3 scripts/pull_contigs.py -i contigs.txt -o refs/ --taxdump taxdump/ -e you@org.org
scripts/pathogen_db/build_pathogen_db.sh -o pathogen_build/ -e you@org.org -t 8
```

Many workflow scripts in `../bin/` also work standalone (`download_fastas.py`, `append_taxid.py`, `resolve_sra.py`, `export_data.py`, `spike_into_background.py`, ...); they are indexed on the same docs page.

## Developer helpers

`test_*.js`, `test_*.py`, `bench_*.js`, `verify_*.js`: interactive-report and metadata tests.
`migrate_wiki.py`, `write_docs_ref.py`, `docs_release_stub.py`, `fetch_offline_report_libs.py`, `split_report_template.py`, `inline_boot_json.py`, `make_heatmap_boot.py`, `extract_heatmap_boot.py`, `reconstruct_paths_json.py`: docs-site and report-bundle maintenance.
