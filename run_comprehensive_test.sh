#!/bin/bash

# Make the comprehensive test executable and run it
chmod +x /Users/dho/Documents/taxtriage/comprehensive_test.py

echo "Executing comprehensive TaxTriage plugin test suite..."
echo ""

# Run the Python test script
python3 /Users/dho/Documents/taxtriage/comprehensive_test.py

echo ""
echo "Test execution completed."

# Check if detailed report was created
REPORT_FILE="/Users/dho/Documents/taxtriage/detailed_test_report.json"
if [ -f "$REPORT_FILE" ]; then
    echo ""
    echo "Detailed test report available at: $REPORT_FILE"
fi