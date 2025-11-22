#!/bin/bash

echo "Testing corrected Nextflow command..."

cd /var/folders/pw/r9x37tvs2q5d6c044x2fhp200000gn/T/taxtriage_20250930_171048

export NXF_ANSI_LOG=false
echo "Set NXF_ANSI_LOG=false"

"/Users/dho/Library/Application Support/Geneious/plugins/com.jhuapl.taxtriage.geneious.TaxTriagePlugin/bin/nextflow" run https://github.com/jhuapl-bio/taxtriage \
  -r main \
  -profile docker \
  --input /var/folders/pw/r9x37tvs2q5d6c044x2fhp200000gn/T/taxtriage_20250930_171048/config/samplesheet.csv \
  --outdir /var/folders/pw/r9x37tvs2q5d6c044x2fhp200000gn/T/taxtriage_20250930_171048/output \
  -work-dir /var/folders/pw/r9x37tvs2q5d6c044x2fhp200000gn/T/taxtriage_20250930_171048/work \
  --platform ILLUMINA \
  --db /Users/dho/.geneious/plugin-data/taxtriage/databases/viral \
  --download_db false \
  --max_cpus 14 \
  --max_memory 8.GB