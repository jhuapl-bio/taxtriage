#!/bin/bash

# TaxTriage Geneious Plugin Test Harness
# This script helps with iterative debugging of the plugin

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLUGIN_DIR="$SCRIPT_DIR"
BUILD_DIR="$PLUGIN_DIR/build"
PLUGIN_FILE="$BUILD_DIR/TaxTriage.gplugin"
GENEIOUS_PLUGINS_DIR="$HOME/Library/Application Support/Geneious/plugins"
TEST_DATA_DIR="/private/var/folders/pw/r9x37tvs2q5d6c044x2fhp200000gn/T"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== TaxTriage Geneious Plugin Test Harness ===${NC}"

# Function to build the plugin
build_plugin() {
    echo -e "${YELLOW}Building plugin...${NC}"
    cd "$PLUGIN_DIR"
    ant clean distribute
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Build successful${NC}"
        return 0
    else
        echo -e "${RED}✗ Build failed${NC}"
        return 1
    fi
}

# Function to install the plugin
install_plugin() {
    echo -e "${YELLOW}Installing plugin to Geneious...${NC}"

    # Create plugins directory if it doesn't exist
    mkdir -p "$GENEIOUS_PLUGINS_DIR"

    # Remove old version if exists
    rm -rf "$GENEIOUS_PLUGINS_DIR/com.jhuapl.taxtriage.geneious.TaxTriagePlugin"
    rm -f "$GENEIOUS_PLUGINS_DIR/TaxTriage.gplugin"

    # Copy new plugin
    cp "$PLUGIN_FILE" "$GENEIOUS_PLUGINS_DIR/"

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Plugin installed to: $GENEIOUS_PLUGINS_DIR${NC}"
        return 0
    else
        echo -e "${RED}✗ Installation failed${NC}"
        return 1
    fi
}

# Function to find recent TaxTriage output directories
find_test_data() {
    echo -e "${YELLOW}Looking for recent TaxTriage output directories...${NC}"

    # Find directories matching the pattern
    find "$TEST_DATA_DIR" -type d -name "taxtriage_*" -maxdepth 2 -mtime -7 2>/dev/null | head -5
}

# Function to test import with a specific directory
test_import() {
    local test_dir="$1"

    if [ -z "$test_dir" ]; then
        echo -e "${RED}Please provide a test directory${NC}"
        return 1
    fi

    echo -e "${YELLOW}Testing import with: $test_dir${NC}"

    # Check what's in the directory
    echo "Contents:"
    ls -la "$test_dir" | head -10

    # Check for BAM files
    echo -e "\nBAM files:"
    find "$test_dir" -name "*.bam" -type f 2>/dev/null

    # Check for reference files
    echo -e "\nReference files:"
    find "$test_dir" -name "*.fasta" -o -name "*.fa" -o -name "*.gb" -o -name "*.gbk" 2>/dev/null

    # Check for output directories
    echo -e "\nOutput structure:"
    find "$test_dir/output" -type d 2>/dev/null | head -10
}

# Function to run the plugin programmatically (if we can)
run_plugin_test() {
    echo -e "${YELLOW}Attempting to run plugin test...${NC}"

    # This would need to be implemented based on how the plugin can be triggered
    # For now, we'll just provide instructions
    echo "To test the plugin:"
    echo "1. Open Geneious Prime"
    echo "2. Go to Tools -> Plugins"
    echo "3. Check that TaxTriage plugin is installed"
    echo "4. Run the TaxTriage workflow"
    echo "5. Check console output (View -> Show Console)"
}

# Function to check plugin logs
check_logs() {
    echo -e "${YELLOW}Checking Geneious logs...${NC}"

    # Look for Geneious log files
    LOG_DIR="$HOME/Library/Logs/Geneious"
    if [ -d "$LOG_DIR" ]; then
        echo "Recent log entries:"
        tail -n 50 "$LOG_DIR"/*.log 2>/dev/null | grep -i "taxtriage\|bam\|reference"
    else
        echo "Log directory not found at: $LOG_DIR"
    fi
}

# Main menu
show_menu() {
    echo -e "\n${GREEN}Choose an option:${NC}"
    echo "1) Build plugin"
    echo "2) Build and install plugin"
    echo "3) Find test data"
    echo "4) Test import with specific directory"
    echo "5) Check logs"
    echo "6) Full cycle (build, install, find data)"
    echo "q) Quit"

    read -p "Enter choice: " choice

    case $choice in
        1)
            build_plugin
            ;;
        2)
            build_plugin && install_plugin
            ;;
        3)
            find_test_data
            ;;
        4)
            read -p "Enter test directory path: " test_path
            test_import "$test_path"
            ;;
        5)
            check_logs
            ;;
        6)
            build_plugin && install_plugin && find_test_data
            ;;
        q)
            echo -e "${GREEN}Goodbye!${NC}"
            exit 0
            ;;
        *)
            echo -e "${RED}Invalid choice${NC}"
            ;;
    esac
}

# Main loop
while true; do
    show_menu
done