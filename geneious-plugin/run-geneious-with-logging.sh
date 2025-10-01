#!/bin/bash

# Alternative script: Runs Geneious with output logged to a file
# You can tail the log file in VS Code while Geneious runs

set -e

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Configuration
PLUGIN_DIR="/Users/dho/Documents/taxtriage/geneious-plugin"
BUILD_DIR="$PLUGIN_DIR/build"
PLUGIN_FILE="$BUILD_DIR/TaxTriage.gplugin"
GENEIOUS_PLUGINS_DIR="$HOME/.geneious/plugins"
GENEIOUS_APP="/Applications/Geneious Prime.app"
LOG_FILE="/tmp/geneious-taxtriage-debug.log"

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}TaxTriage Plugin Debug Deployment${NC}"
echo -e "${GREEN}(Logging to file)${NC}"
echo -e "${GREEN}========================================${NC}"

# Check if build exists
if [ ! -f "$PLUGIN_FILE" ]; then
    echo -e "${YELLOW}Building plugin...${NC}"
    cd "$PLUGIN_DIR"
    ant clean distribute
fi

# Stop Geneious
echo -e "${GREEN}Stopping Geneious...${NC}"
if pgrep -f "Geneious" > /dev/null; then
    osascript -e 'tell application "Geneious Prime" to quit' 2>/dev/null || true
    sleep 3
    pkill -9 -f "Geneious" 2>/dev/null || true
    sleep 1
fi

# Clean and install plugin
echo -e "${GREEN}Installing plugin...${NC}"
rm -rf "$GENEIOUS_PLUGINS_DIR/com.jhuapl.taxtriage.geneious.TaxTriagePlugin"
rm -f "$GENEIOUS_PLUGINS_DIR"/*.gplugin 2>/dev/null || true
mkdir -p "$GENEIOUS_PLUGINS_DIR"
cp "$PLUGIN_FILE" "$GENEIOUS_PLUGINS_DIR/"

echo -e "${GREEN}✓ Plugin installed${NC}"
echo ""
echo -e "${YELLOW}========================================${NC}"
echo -e "${YELLOW}Starting Geneious with logging${NC}"
echo -e "${YELLOW}========================================${NC}"
echo ""
echo "Log file: ${LOG_FILE}"
echo ""
echo "To view log in real-time, run in another terminal:"
echo -e "${GREEN}  tail -f ${LOG_FILE}${NC}"
echo ""
echo "Or open in VS Code:"
echo -e "${GREEN}  code ${LOG_FILE}${NC}"
echo ""

# Clear previous log
> "$LOG_FILE"

# Start Geneious and log output
echo "Starting Geneious Prime..."
"$GENEIOUS_APP/Contents/MacOS/JavaApplicationStub" > "$LOG_FILE" 2>&1 &

GENEIOUS_PID=$!
echo "Geneious started with PID: $GENEIOUS_PID"
echo ""
echo "Showing last 30 lines of log (updates every 2 seconds):"
echo "Press Ctrl+C to stop watching (Geneious will keep running)"
echo ""

# Tail the log file
tail -f "$LOG_FILE"
