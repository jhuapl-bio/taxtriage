#!/bin/bash

# Make the test runner executable and run it
# chmod +x /Users/dho/Documents/taxtriage/test_runner.sh

# Execute the test runner
bash /Users/dho/Documents/taxtriage/test_runner.sh

# Display the test report
echo ""
echo "=========================================="
echo "          TEST REPORT SUMMARY"
echo "=========================================="
cat /Users/dho/Documents/taxtriage/test_report.txt
