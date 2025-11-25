#!/bin/bash

# Simulate exactly what the plugin does to run Nextflow
# This will help us identify the issue

set -e

echo "=========================================="
echo "Simulating Plugin Nextflow Execution"
echo "=========================================="
echo ""

# Create a new test workspace
WORKSPACE="/tmp/test-plugin-nextflow-$(date +%Y%m%d_%H%M%S)"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

echo "Created test workspace: $WORKSPACE"
mkdir -p config input output preprocessed/bbtools_work work

# Copy files from a recent plugin run
RECENT_RUN="/var/folders/pw/r9x37tvs2q5d6c044x2fhp200000gn/T/taxtriage_20250930_174123"

echo "Copying config files from recent run..."
cp "$RECENT_RUN/config/samplesheet.csv" config/
cp "$RECENT_RUN/config/params.json" config/
cp "$RECENT_RUN/config/nextflow.config" config/ 2>/dev/null || echo "No nextflow.config to copy"

echo "Copying preprocessed files..."
cp -r "$RECENT_RUN/preprocessed/bbtools_work/"* preprocessed/bbtools_work/

echo ""
echo "Files in workspace:"
ls -R "$WORKSPACE"
echo ""

# Use the same Nextflow binary the plugin uses
NEXTFLOW="/Users/dho/Library/Application Support/Geneious/plugins/com.jhuapl.taxtriage.geneious.TaxTriagePlugin/bin/nextflow"

echo "=========================================="
echo "Test 1: Check Nextflow binary"
echo "=========================================="
"$NEXTFLOW" -version
echo ""

# Set up environment EXACTLY as the plugin does
export JAVA_HOME=/usr
export JAVA_CMD=java
export NXF_HOME="$HOME/.nextflow"
export NXF_ANSI_LOG=false

echo "=========================================="
echo "Test 2: Run Nextflow with exact plugin command"
echo "=========================================="
echo "Working directory: $WORKSPACE"
echo "Environment:"
echo "  JAVA_HOME=$JAVA_HOME"
echo "  NXF_HOME=$NXF_HOME"
echo "  NXF_ANSI_LOG=$NXF_ANSI_LOG"
echo ""

# Build the exact command the plugin uses
CMD=(
    "$NEXTFLOW"
    "run"
    "https://github.com/jhuapl-bio/taxtriage"
    "-r" "main"
    "-profile" "docker"
    "--input" "$WORKSPACE/config/samplesheet.csv"
    "--outdir" "$WORKSPACE/output"
    "-work-dir" "$WORKSPACE/work"
    "--platform" "ILLUMINA"
    "--db" "/Users/dho/.geneious/plugin-data/taxtriage/databases/viral"
    "--download_db" "false"
    "--max_cpus" "14"
    "--max_memory" "8.GB"
)

echo "Command: ${CMD[@]}"
echo ""
echo "Starting Nextflow..."
echo ""

# Run and capture output
"${CMD[@]}" 2>&1 | tee "$WORKSPACE/nextflow.log"

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Results"
echo "=========================================="
echo "Exit code: $EXIT_CODE"
echo ""
echo "Work directory contents:"
ls -la "$WORKSPACE/work" 2>/dev/null || echo "(empty)"
echo ""
echo "Output directory contents:"
ls -la "$WORKSPACE/output" 2>/dev/null || echo "(empty)"
echo ""
echo "Full log saved to: $WORKSPACE/nextflow.log"
echo "Workspace: $WORKSPACE"
