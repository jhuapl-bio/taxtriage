#!/bin/bash

echo "Starting final comprehensive test of TaxTriage Geneious plugin..."
echo ""

# Make the test script executable
chmod +x /Users/dho/Documents/taxtriage/final_test_execution.py

# Run the comprehensive test
python3 /Users/dho/Documents/taxtriage/final_test_execution.py

# Capture the exit code
EXIT_CODE=$?

echo ""
echo "Test execution completed with exit code: $EXIT_CODE"

# Create a summary file
SUMMARY_FILE="/Users/dho/Documents/taxtriage/test_summary.txt"
echo "TaxTriage Plugin Test Summary" > "$SUMMARY_FILE"
echo "=============================" >> "$SUMMARY_FILE"
echo "Date: $(date)" >> "$SUMMARY_FILE"
echo "Exit Code: $EXIT_CODE" >> "$SUMMARY_FILE"

case $EXIT_CODE in
    0)
        echo "Status: SUCCESS - Plugin ready for installation" >> "$SUMMARY_FILE"
        ;;
    1)
        echo "Status: PARTIAL SUCCESS - Minor issues detected" >> "$SUMMARY_FILE"
        ;;
    2)
        echo "Status: FAILURE - Critical issues need resolution" >> "$SUMMARY_FILE"
        ;;
    *)
        echo "Status: UNKNOWN - Unexpected exit code" >> "$SUMMARY_FILE"
        ;;
esac

echo ""
echo "Summary saved to: $SUMMARY_FILE"