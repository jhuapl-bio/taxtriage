# GenBank Field Fix for Reference Matching

## Problem Identified
The BAM files use the full VERSION identifier (e.g., "NC_045512.2") as the reference name, but the downloaded GenBank files had:
- LOCUS: NC_045512 (missing the .2 suffix)
- ACCESSION: NC_045512 (missing the .2 suffix)
- VERSION: NC_045512.2 (correct)

This mismatch prevented Geneious from properly matching references during BAM import.

## Solution Implemented
Created `GenBankFileFixer.java` that automatically fixes GenBank files to ensure all three fields (LOCUS, ACCESSION, VERSION) contain the same identifier.

### Before Fix:
```
LOCUS       NC_045512              29903 bp
ACCESSION   NC_045512
VERSION     NC_045512.2
```

### After Fix:
```
LOCUS       NC_045512.2     29903 bp
ACCESSION   NC_045512.2
VERSION     NC_045512.2
```

## Integration Points

1. **During Download** - NCBIReferenceDownloader automatically fixes GenBank files immediately after downloading from NCBI

2. **Before Import** - SampleBasedImporter checks and fixes any existing GenBank files in the minimap2 and genbank_downloads folders before importing

3. **Validation** - The fixer includes a validation method to verify all fields match

## How It Works

The GenBankFileFixer:
1. Reads the GenBank file and finds the VERSION field value
2. Updates the LOCUS line to use the VERSION value while preserving formatting
3. Updates the ACCESSION line to use the VERSION value
4. Validates that all three fields now match

## Benefits

- **Automatic Reference Matching** - BAM files can now find their references because the names match exactly
- **No Manual Intervention** - Files are automatically fixed during the import process
- **Backwards Compatible** - Works with existing GenBank files and newly downloaded ones

## Testing

Tested on actual GenBank file NC_045512.2.gb:
- Successfully changed LOCUS from "NC_045512" to "NC_045512.2"
- Successfully changed ACCESSION from "NC_045512" to "NC_045512.2"
- Validation confirms all fields now match

## Plugin Status

The updated plugin with GenBank field fixing is ready at:
`/Users/dho/Documents/taxtriage/geneious-plugin/build/TaxTriage.gplugin`

This fix, combined with co-locating references in the minimap2 folder, should resolve all BAM import reference matching issues.