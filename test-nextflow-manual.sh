#!/bin/bash

# Manual test of Nextflow workflow using plugin-generated config files

set -e

WORKSPACE="/private/var/folders/pw/r9x37tvs2q5d6c044x2fhp200000gn/T/taxtriage_20250930_163456"

echo "=========================================="
echo "Manual Nextflow Workflow Test"
echo "=========================================="
echo ""
echo "Workspace: $WORKSPACE"
echo ""

# Check that workspace exists
if [ ! -d "$WORKSPACE" ]; then
    echo "ERROR: Workspace directory not found!"
    echo "Please update WORKSPACE variable to point to latest taxtriage_* directory"
    exit 1
fi

# Check config files exist
echo "Checking config files..."
if [ ! -f "$WORKSPACE/config/samplesheet.csv" ]; then
    echo "ERROR: samplesheet.csv not found!"
    exit 1
fi
echo "  ✓ samplesheet.csv"

if [ ! -f "$WORKSPACE/config/params.json" ]; then
    echo "ERROR: params.json not found!"
    exit 1
fi
echo "  ✓ params.json"

# Check database exists
DB_PATH="/Users/dho/.geneious/plugin-data/taxtriage/databases/viral"
echo ""
echo "Checking database..."
if [ -d "$DB_PATH" ]; then
    echo "  ✓ Database found at: $DB_PATH"
    ls -lh "$DB_PATH" | head -5
else
    echo "  ✗ Database NOT found at: $DB_PATH"
    echo "  Will attempt to download during workflow..."
fi

echo ""
echo "=========================================="
echo "Starting Nextflow Workflow"
echo "=========================================="
echo ""

cd "$WORKSPACE"

# Run the Nextflow command exactly as the plugin would
nextflow run https://github.com/jhuapl-bio/taxtriage \
  -r main \
  -profile docker \
  --input "$WORKSPACE/config/samplesheet.csv" \
  --outdir "$WORKSPACE/output" \
  -work-dir "$WORKSPACE/work" \
  --platform ILLUMINA \
  --db "$DB_PATH" \
  --download_db false \
  --max_cpus 14 \
  --max_memory 8.GB \
  2>&1 | tee /tmp/nextflow-manual-run.log

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Workflow completed with exit code: $EXIT_CODE"
echo "=========================================="
echo ""
echo "Output directory: $WORKSPACE/output"
echo "Work directory: $WORKSPACE/work"
echo "Full log: /tmp/nextflow-manual-run.log"
echo ""

if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ SUCCESS!"
    echo ""
    echo "Output files:"
    ls -lh "$WORKSPACE/output" 2>/dev/null || echo "  (empty)"
else
    echo "✗ FAILED"
    echo ""
    echo "Check the log for errors: /tmp/nextflow-manual-run.log"
fi
