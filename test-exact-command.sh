#!/bin/bash

# Test the exact command the plugin is trying to run

set -e

echo "Testing exact Nextflow command from plugin..."
echo ""

NEXTFLOW="/Users/dho/Library/Application Support/Geneious/plugins/com.jhuapl.taxtriage.geneious.TaxTriagePlugin/bin/nextflow"
WORKSPACE="/var/folders/pw/r9x37tvs2q5d6c044x2fhp200000gn/T/taxtriage_20250930_165738"

echo "Nextflow binary: $NEXTFLOW"
echo "Workspace: $WORKSPACE"
echo ""

cd "$WORKSPACE"

echo "Running command..."
"$NEXTFLOW" run https://github.com/jhuapl-bio/taxtriage \
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