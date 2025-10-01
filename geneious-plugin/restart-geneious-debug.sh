#!/bin/bash

# Script to deploy TaxTriage plugin and restart Geneious with debug output
# This allows viewing all console output including System.out.println statements

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}TaxTriage Plugin Debug Deployment${NC}"
echo -e "${GREEN}========================================${NC}"

# Configuration
PLUGIN_DIR="/Users/dho/Documents/taxtriage/geneious-plugin"
BUILD_DIR="$PLUGIN_DIR/build"
PLUGIN_FILE="$BUILD_DIR/TaxTriage.gplugin"
GENEIOUS_PLUGINS_DIR="$HOME/.geneious/plugins"
GENEIOUS_APP="/Applications/Geneious Prime.app"

# Check if build exists
if [ ! -f "$PLUGIN_FILE" ]; then
    echo -e "${YELLOW}Plugin not built yet. Building now...${NC}"
    cd "$PLUGIN_DIR"
    ant clean distribute
fi

echo -e "${GREEN}Step 1: Stopping Geneious if running...${NC}"
if pgrep -f "Geneious" > /dev/null; then
    echo "  Geneious is running, sending quit command..."
    osascript -e 'tell application "Geneious Prime" to quit'

    # Wait for Geneious to fully quit (max 10 seconds)
    for i in {1..10}; do
        if ! pgrep -f "Geneious" > /dev/null; then
            echo "  Geneious quit successfully"
            break
        fi
        echo "  Waiting for Geneious to quit... ($i/10)"
        sleep 1
    done

    # Force quit if still running
    if pgrep -f "Geneious" > /dev/null; then
        echo -e "${YELLOW}  Force quitting Geneious...${NC}"
        pkill -9 -f "Geneious"
        sleep 2
    fi
else
    echo "  Geneious is not running"
fi

echo -e "${GREEN}Step 2: Cleaning old plugin installations...${NC}"
# Remove old plugin files
rm -rf "$GENEIOUS_PLUGINS_DIR/com.jhuapl.taxtriage.geneious.TaxTriagePlugin"
rm -f "$GENEIOUS_PLUGINS_DIR/TaxTriage.gplugin"
rm -f "$GENEIOUS_PLUGINS_DIR"/*.gplugin 2>/dev/null || true
echo "  Cleaned old installations"

echo -e "${GREEN}Step 3: Installing new plugin...${NC}"
# Ensure plugins directory exists
mkdir -p "$GENEIOUS_PLUGINS_DIR"

# Copy the plugin
cp "$PLUGIN_FILE" "$GENEIOUS_PLUGINS_DIR/"
echo "  Copied: $PLUGIN_FILE"
echo "  To: $GENEIOUS_PLUGINS_DIR/"

# Verify installation
if [ -f "$GENEIOUS_PLUGINS_DIR/TaxTriage.gplugin" ]; then
    echo -e "${GREEN}  ✓ Plugin installed successfully${NC}"
    echo "  Plugin size: $(ls -lh "$GENEIOUS_PLUGINS_DIR/TaxTriage.gplugin" | awk '{print $5}')"
else
    echo -e "${RED}  ✗ Plugin installation failed${NC}"
    exit 1
fi

echo ""
echo -e "${YELLOW}========================================${NC}"
echo -e "${YELLOW}Starting Geneious from shell...${NC}"
echo -e "${YELLOW}All debug output will appear here${NC}"
echo -e "${YELLOW}========================================${NC}"
echo ""
echo "Watch for these debug messages:"
echo "  - 'TAXTRIAGE PLUGIN EXECUTION STARTED'"
echo "  - 'PREPROCESSING CONFIGURATION:'"
echo "  - 'WORKFLOW CONFIGURATION COMPLETE'"
echo "  - 'About to call executeWorkflow'"
echo "  - 'Calling executeWorkflow method...'"
echo ""
echo -e "${GREEN}Starting Geneious Prime...${NC}"
echo ""

# Start Geneious from command line to capture all output
# Using open -W to wait and capture output
"$GENEIOUS_APP/Contents/MacOS/JavaApplicationStub" 2>&1 | tee /tmp/geneious-debug.log

# Note: The script will block here until Geneious is closed
