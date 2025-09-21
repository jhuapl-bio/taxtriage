# Co-Located References Update

## Problem
BAM files were failing to find reference sequences during import because the references and BAM files were in different directories, making it harder for Geneious to match them automatically.

## Solution
Modified the plugin to download GenBank reference files directly into the `minimap2` folder where the BAM files are located. This co-location strategy improves the chances of successful reference matching during BAM import.

## Key Changes

### 1. Automatic GenBank Download
- Added `downloadGenBankReferencesForBamFiles()` method to SampleBasedImporter
- Automatically extracts reference accessions from BAM files
- Downloads missing GenBank files from NCBI
- Saves them in the same `minimap2` folder as the BAM files

### 2. Import Order
The plugin now follows this sequence:
1. **Download GenBank references** to minimap2 folder (co-located with BAMs)
2. **Import text files** using DirectTextImporter
3. **Import GenBank references** from minimap2 folder to References folder
4. **Import BAM files** with auto-reference detection

### 3. Reference Search Locations
The plugin now checks for GenBank references in multiple locations:
- `minimap2/` folder (co-located with BAM files) - PRIMARY
- `genbank_downloads/` folder (legacy location)
- `download/` folder (workflow output)

## Benefits

1. **Better Reference Matching** - Co-locating references with BAM files improves Geneious's ability to automatically detect and use the correct reference sequences

2. **Automatic Download** - The plugin automatically downloads any missing GenBank references from NCBI based on the accessions found in BAM headers

3. **Complete Annotations** - GenBank files include full annotations (genes, CDS, etc.) unlike the workflow-generated FASTA files

## Console Output
The plugin now provides detailed console output showing:
- Which references are needed by BAM files
- Which GenBank files are being downloaded
- Where files are being imported from
- Success/failure status of each import

## Testing
The updated plugin is ready at:
`/Users/dho/Documents/taxtriage/geneious-plugin/build/TaxTriage.gplugin`

When you run the plugin, it will:
1. Automatically detect references needed by BAM files
2. Download missing GenBank files to the minimap2 folder
3. Import all files with proper folder organization
4. Use co-located references for better BAM import success