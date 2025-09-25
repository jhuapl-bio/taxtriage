# BAM Import Code Extract for Geneious Plugin

## Problem Description
BAM files are not automatically finding their reference sequences during programmatic import, even though:
1. GenBank reference files are successfully imported to the References folder
2. Manual drag-and-drop of BAM files works correctly and finds references
3. The GenBank files have correct field matching (LOCUS/ACCESSION/VERSION all match, e.g., NC_045512.2)
4. BAM and GenBank files are co-located in the minimap2 folder

When importing programmatically, Geneious shows a dialog saying it can't find reference sequences. When the user manually selects the References folder in the dialog, the import works correctly.

## Current Import Flow in SampleBasedImporter.java

### Step 1: Import GenBank References (WORKING)
```java
// From SampleBasedImporter.java lines 354-420
// Check for GenBank downloads in minimap2 folder (co-located with BAM files)
Path minimap2Dir = outputDir.resolve("minimap2");
if (Files.exists(minimap2Dir)) {
    // Debug output shows GenBank files are found
    System.out.println("  === GenBank Import Debug ===");
    System.out.println("  Sample name: " + sample);
    System.out.println("  Minimap2 directory: " + minimap2Dir);

    // Fix GenBank files to ensure LOCUS/ACCESSION/VERSION match
    int fixedCount = GenBankFileFixer.fixGenBankFilesInDirectory(minimap2Dir);

    // Find GenBank files (using special method since they're named by accession, not sample)
    List<File> gbFiles = findGenBankReferenceFiles(Arrays.asList(minimap2Dir));

    System.out.println("  GenBank reference files found: " + gbFiles.size());

    for (File file : gbFiles) {
        // Validate and fix if needed
        if (!GenBankFileFixer.validateGenBankFile(file)) {
            boolean fixed = GenBankFileFixer.fixGenBankFile(file);
        }

        System.out.println("  Importing GenBank from minimap2 to References: " + file.getName());

        // Import the documents
        List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(file, progressListener);

        if (docs != null && !docs.isEmpty()) {
            WritableDatabaseService referencesFolder = folders.get("References");

            // Add each document to the References folder
            for (AnnotatedPluginDocument doc : docs) {
                AnnotatedPluginDocument copiedDoc = referencesFolder.addDocumentCopy(doc, ProgressListener.EMPTY);
                if (copiedDoc != null) {
                    allImported.add(copiedDoc);
                    System.out.println("      - Added to References: " + copiedDoc.getName());
                }
            }
        }
    }
}
```

### Step 2: Collect Reference Documents (lines 506-528)
```java
// Collect all imported reference documents for BAM linking
List<AnnotatedPluginDocument> importedReferences = new ArrayList<>();
if (folders.containsKey("References")) {
    try {
        WritableDatabaseService referencesFolder = folders.get("References");
        // Filter to get only reference sequences (GenBank documents)
        for (AnnotatedPluginDocument doc : allImported) {
            if (doc.getDocumentClass() != null &&
                (doc.getDocumentClass().getName().contains("NucleotideSequence") ||
                 doc.getDocumentClass().getName().contains("DefaultNucleotideSequence"))) {
                importedReferences.add(doc);
            }
        }
        System.out.println("  Collected " + importedReferences.size() + " reference sequences for BAM import");
        for (AnnotatedPluginDocument ref : importedReferences) {
            System.out.println("    - Reference: " + ref.getName());
        }
    } catch (Exception e) {
        logger.warning("Could not collect reference documents: " + e.getMessage());
    }
}
```

### Step 3: Import BAM Files (FAILING - lines 551-594)
```java
// Current implementation that still shows reference dialog
for (File bamFile : bamFiles) {
    System.out.println("  Importing BAM file: " + bamFile.getName());
    try {
        // Import BAM directly to the database folder containing references
        // This SHOULD ensure Geneious can find references automatically
        WritableDatabaseService targetFolder;

        if (!importedReferences.isEmpty() && folders.containsKey("References")) {
            // Import to References folder where the reference sequences are
            targetFolder = folders.get("References");
            System.out.println("    Importing BAM directly to References folder with " +
                              importedReferences.size() + " reference(s)");

            // Import BAM file directly to the folder
            List<AnnotatedPluginDocument> docs = PluginUtilities.importDocumentsToDatabase(
                bamFile, targetFolder, progressListener);

            if (docs != null && !docs.isEmpty()) {
                allImported.addAll(docs);
                System.out.println("    Successfully imported BAM with " + docs.size() +
                                  " alignment(s) to References folder");
                for (AnnotatedPluginDocument doc : docs) {
                    System.out.println("      - " + doc.getName());
                }
            }
        }
    } catch (Exception e) {
        logger.log(Level.WARNING, "Failed to import BAM file: " + bamFile.getName(), e);
    }
}
```

## Previous Attempts That Also Failed

### Attempt 1: Pass References to BamFileImporter
```java
// This was tried but still showed the dialog
if (!importedReferences.isEmpty()) {
    docs = BamFileImporter.importBamFile(bamFile, progressListener, importedReferences);
}
```

### Attempt 2: Use BamImportHelper with References
```java
// In BamFileImporter.java
if (referenceDocuments != null && !referenceDocuments.isEmpty()) {
    List<AnnotatedPluginDocument> result = BamImportHelper.importBamWithReferences(
        bamFile, referenceDocuments, progressListener);
}
```

### Attempt 3: Import with Options
```java
// In BamImportHelper.java - tried passing reference URNs
Map<String, String> options = new HashMap<>();
options.put("referenceSequences", refUrns.toString());
options.put("referenceURNs", refUrns.toString());
options.put("findReferences", "true");
options.put("autoFindReferences", "true");

List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(
    bamFile, options, progressListener);
```

## Key Information

### BAM File Headers
The BAM files contain reference names like:
- `NC_045512.2` (full version identifier)
- `NC_049394.1`

### GenBank Files After Fixing
After GenBankFileFixer processes them:
- LOCUS: `NC_045512.2` (was `NC_045512`)
- ACCESSION: `NC_045512.2` (was `NC_045512`)
- VERSION: `NC_045512.2` (unchanged)

### What Works
1. Manual drag-and-drop of BAM files when references are already in Geneious
2. Selecting the References folder when prompted by the dialog
3. GenBank import to References folder
4. Text file import

### What Doesn't Work
Automatic reference discovery during programmatic BAM import using:
- `PluginUtilities.importDocumentsToDatabase(bamFile, referencesFolder, progressListener)`
- `PluginUtilities.importDocuments(bamFile, progressListener)` with reference documents passed
- `BamImportHelper.importBamWithReferences()` with reference documents

## Questions for the Developer

1. **Is there a specific import option or parameter** that tells Geneious where to look for reference sequences during BAM import?

2. **Does `importDocumentsToDatabase()` support any options** to specify reference location or to suppress the reference selection dialog?

3. **Is there a way to pre-register or index references** so they're found automatically during BAM import?

4. **What is the difference between manual drag-and-drop and programmatic import** that causes references to be found in one case but not the other?

5. **Are there any timing issues** where references need to be fully indexed/committed before BAM import?

6. **Is there an alternative API method** specifically for importing alignment files with known references?

## Console Output When Dialog Appears
```
  Importing BAM file: test.viral.bam
    Importing BAM directly to References folder with 3 reference(s)
[Dialog appears asking to select reference location]
[User selects References folder]
    Successfully imported BAM with 1 alignment(s) to References folder
      - test.viral on NC_045512.2
```

## Environment
- Geneious Prime (version not specified in code)
- Plugin uses Geneious Public API
- Java 11 compilation target

## Files Involved
1. `/Users/dho/Documents/taxtriage/geneious-plugin/src/com/jhuapl/taxtriage/geneious/importers/SampleBasedImporter.java`
2. `/Users/dho/Documents/taxtriage/geneious-plugin/src/com/jhuapl/taxtriage/geneious/importers/BamFileImporter.java`
3. `/Users/dho/Documents/taxtriage/geneious-plugin/src/com/jhuapl/taxtriage/geneious/importers/BamImportHelper.java`
4. `/Users/dho/Documents/taxtriage/geneious-plugin/src/com/jhuapl/taxtriage/geneious/importers/GenBankFileFixer.java`

## What We Need
A way to programmatically import BAM files that automatically finds and links to reference sequences already in the Geneious database, without showing the reference selection dialog.