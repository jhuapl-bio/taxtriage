#!/bin/bash

# Test script to verify bundled Nextflow can run the TaxTriage workflow
# Uses files from the latest plugin working directory

set -e

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Testing Bundled Nextflow with TaxTriage${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# Extract bundled Nextflow from plugin
PLUGIN_FILE="/Users/dho/.geneious/plugins/TaxTriage.gplugin"
EXTRACT_DIR="/tmp/test-bundled-nextflow"
NEXTFLOW_BIN="$EXTRACT_DIR/com.jhuapl.taxtriage.geneious.TaxTriagePlugin/bin/nextflow"

echo "Step 1: Extracting bundled Nextflow from plugin..."
rm -rf "$EXTRACT_DIR"
mkdir -p "$EXTRACT_DIR"
cd "$EXTRACT_DIR"
unzip -q "$PLUGIN_FILE"

if [ ! -f "$NEXTFLOW_BIN" ]; then
    echo -e "${RED}ERROR: Nextflow binary not found in plugin!${NC}"
    exit 1
fi

chmod +x "$NEXTFLOW_BIN"
echo -e "${GREEN}✓ Bundled Nextflow extracted${NC}"
echo ""

# Test Nextflow version
echo "Step 2: Testing bundled Nextflow binary..."
"$NEXTFLOW_BIN" -version
echo ""

# Use the latest working directory from plugin
WORKSPACE="/var/folders/pw/r9x37tvs2q5d6c044x2fhp200000gn/T/taxtriage_20250930_164417"

if [ ! -d "$WORKSPACE" ]; then
    echo -e "${YELLOW}WARNING: Latest workspace not found${NC}"
    echo "Looking for most recent taxtriage workspace..."
    WORKSPACE=$(ls -td /var/folders/pw/r9x37tvs2q5d6c044x2fhp200000gn/T/taxtriage_* 2>/dev/null | head -1)
    if [ -z "$WORKSPACE" ]; then
        echo -e "${RED}ERROR: No taxtriage workspace found!${NC}"
        echo "Please run the Geneious plugin first to generate workspace files"
        exit 1
    fi
fi

echo "Step 3: Using workspace: $WORKSPACE"
echo ""

# Verify required files exist
echo "Step 4: Verifying workspace files..."
if [ ! -f "$WORKSPACE/config/samplesheet.csv" ]; then
    echo -e "${RED}ERROR: samplesheet.csv not found!${NC}"
    exit 1
fi
echo "  ✓ samplesheet.csv"

# Check database
DB_PATH="/Users/dho/.geneious/plugin-data/taxtriage/databases/viral"
if [ -d "$DB_PATH" ]; then
    echo "  ✓ Database found at: $DB_PATH"
else
    echo -e "${YELLOW}  ⚠ Database not found, will attempt to download${NC}"
fi
echo ""

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Running TaxTriage Workflow${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Workspace: $WORKSPACE"
echo "Nextflow: $NEXTFLOW_BIN"
echo ""

cd "$WORKSPACE"

# Run the workflow using bundled Nextflow
"$NEXTFLOW_BIN" run https://github.com/jhuapl-bio/taxtriage \
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
  2>&1 | tee /tmp/bundled-nextflow-workflow.log

EXIT_CODE=$?

echo ""
echo -e "${GREEN}========================================${NC}"
echo "Workflow exit code: $EXIT_CODE"
echo -e "${GREEN}========================================${NC}"
echo ""

if [ $EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}✓ SUCCESS!${NC}"
    echo ""
    echo "Output directory:"
    ls -lh "$WORKSPACE/output" 2>/dev/null || echo "  (empty)"
else
    echo -e "${RED}✗ FAILED${NC}"
    echo ""
    echo "Check the log for errors:"
    echo "  /tmp/bundled-nextflow-workflow.log"
fi

echo ""
echo "Full log saved to: /tmp/bundled-nextflow-workflow.log"
