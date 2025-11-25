#!/bin/bash

# Test the EXACT command the plugin is running with the extracted Nextflow

echo "Testing extracted bundled Nextflow with exact command..."
echo ""

NEXTFLOW="/Users/dho/Library/Application Support/Geneious/plugins/com.jhuapl.taxtriage.geneious.TaxTriagePlugin/bin/nextflow"
WORKSPACE="/var/folders/pw/r9x37tvs2q5d6c044x2fhp200000gn/T/taxtriage_20250930_171048"

cd "$WORKSPACE"

echo "Testing command with -ansi-log false flag..."
echo ""

# Run EXACTLY what the plugin is running
"$NEXTFLOW" -ansi-log false run https://github.com/jhuapl-bio/taxtriage \
  -r main \
  -profile docker \
  --input "$WORKSPACE/config/samplesheet.csv" \
  --outdir "$WORKSPACE/output" \
  -work-dir "$WORKSPACE/work" \
  --platform ILLUMINA \
  --db /Users/dho/.geneious/plugin-data/taxtriage/databases/viral \
  --download_db false \
  --max_cpus 14 \
  --max_memory 8.GB