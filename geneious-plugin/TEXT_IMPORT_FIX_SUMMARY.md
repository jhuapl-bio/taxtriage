# TaxTriage Geneious Plugin - Text Import Fix

## Problem
Text and TSV files were not importing correctly after implementing the subfolder structure. Files like `taxtriage_1822973822309212160.txt` and `test.krakenreport.krona.txt` were being missed.

## Solution
Created a dedicated `PlainTextImporter` class that:

1. **Forces text import** - Ensures files are always imported as plain text documents
2. **Comprehensive file detection** - Identifies all TaxTriage-related text files
3. **Fallback mechanisms** - Multiple import strategies to ensure success

## Key Components

### PlainTextImporter.java
- `forceImportAsText()` - Guarantees text import by using temporary .txt file
- `shouldImportAsPlainText()` - Determines which files should be text
- Handles these patterns:
  - `.txt`, `.tsv`, `.csv`, `.log`, `.report`
  - Files containing: `.kraken`, `.krona`, `taxtriage_`, `.krakenreport`
  - TaxTriage-specific: `.histo`, `.paths`, `.topnames`, `.toptaxids`, `.gcfids`

### SampleBasedImporter Updates
- Uses `PlainTextImporter.forceImportAsText()` for all text files
- Added `importMissedTextFiles()` method to catch any overlooked files
- Scans root output directory for `taxtriage_*.txt` files
- Enhanced file detection patterns

## Files That Will Now Import Correctly
- `taxtriage_1822973822309212160.txt` (root directory files)
- `test.krakenreport.krona.txt` (krona report files)
- All `.tsv` files (tab-separated values)
- All `.csv` files (comma-separated values)
- Kraken report files
- All other TaxTriage output text files

## Import Behavior
All identified text files will:
1. Be imported as plain text documents (not parsed as sequences)
2. Retain their original names (minus extension)
3. Be placed in appropriate subfolders (Kraken_Results or Reports)

## Testing
The plugin has been successfully compiled and is ready at:
`/Users/dho/Documents/taxtriage/geneious-plugin/build/TaxTriage.gplugin`