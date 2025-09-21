# TaxTriage Geneious Plugin Import Analysis and Solution

## Problem Analysis

### Root Causes of BAM Import Failure

1. **Double Import of References**
   - References were being imported TWICE: once as GenBank and once from FASTA
   - This created confusion and potential naming conflicts
   - Located in `TaxTriageResultImporter.java` lines 187-205 (GenBank + FASTA fallback)

2. **Database Indexing Timing Issue**
   - References were imported and immediately followed by BAM import
   - Geneious database indexing is asynchronous
   - BAM importer couldn't find references that weren't yet searchable
   - Even with `ForceIndexingHelper.forceIndexing()`, timing was insufficient

3. **Manual vs Programmatic Import Difference**
   - Manual import works because user imports references first, waits, then imports BAM
   - Natural delay allows database indexing to complete
   - Programmatic import happens too quickly in sequence

## Solution: Three-Phase Import Architecture

### New Architecture Overview

The solution restructures the import into three completely separate operations:

1. **Phase 1: References Only** - Import ONLY GenBank files (no FASTA)
2. **Phase 2: Text Files** - Import text reports after references are committed
3. **Phase 3: BAM Files** - Import BAM files after references are fully indexed

### Implementation Details

#### 1. New Import Phases Enum

```java
public enum ImportPhase {
    REFERENCES_ONLY,  // Phase 1: Import only GenBank reference sequences
    TEXT_ONLY,        // Phase 2: Import only text files
    BAM_ONLY,         // Phase 3: Import only BAM files
    ALL               // Import all phases sequentially
}
```

#### 2. Phase-Specific Import Methods

- `importReferencesOnly()` - Downloads and imports only GenBank files from NCBI
- `importTextOnly()` - Imports only text report files
- `importBamOnly()` - Imports only BAM files, finding existing references

#### 3. Database Commit Management

- `waitForDatabaseCommit()` - Ensures references are fully committed and indexed
- Extended timeout from 10 seconds to 30 seconds
- Additional 2-second buffer for database operations

#### 4. Enhanced Reference Resolution

- `findAvailableReferences()` - Locates existing references in database for BAM import
- Better logging to diagnose reference matching issues
- Eliminates FASTA fallback to prevent double-import

## Key Code Changes

### File: `TaxTriageResultImporter.java`

1. **Added ImportPhase enum** (lines 41-47)
2. **Added phase-specific import methods** (lines 1137-1360)
3. **Enhanced database commit waiting** (lines 1271-1294)
4. **Deprecated FASTA fallback** to prevent double-import

### File: `PhaseBasedImportWorkflow.java` (NEW)

Demonstrates proper usage of the new phase-based import:
- Complete workflow with proper timing
- Individual phase import for troubleshooting
- Usage examples and documentation

## Usage Examples

### Option 1: Complete Phase-Based Workflow (Recommended)

```java
TaxTriageResultImporter importer = new TaxTriageResultImporter();

// Import all phases with proper timing
List<AnnotatedPluginDocument> docs = importer.importResults(
    outputDir, progressListener, ImportPhase.ALL);
```

### Option 2: Individual Phases (For Troubleshooting)

```java
TaxTriageResultImporter importer = new TaxTriageResultImporter();

// Phase 1: References only
List<AnnotatedPluginDocument> refs = importer.importResults(
    outputDir, progressListener, ImportPhase.REFERENCES_ONLY);

// Wait or do other work...
Thread.sleep(10000);

// Phase 2: Text files
List<AnnotatedPluginDocument> texts = importer.importResults(
    outputDir, progressListener, ImportPhase.TEXT_ONLY);

// Phase 3: BAM files (should now find references)
List<AnnotatedPluginDocument> bams = importer.importResults(
    outputDir, progressListener, ImportPhase.BAM_ONLY);
```

### Option 3: Using the Workflow Helper

```java
// Automatically handles timing and phase separation
List<AnnotatedPluginDocument> allDocs = PhaseBasedImportWorkflow
    .executePhaseBasedImport(outputDir, progressListener);
```

## Why This Fixes the BAM Import Issue

1. **Eliminates Double-Import**
   - Only GenBank files are imported in Phase 1
   - No FASTA fallback prevents duplicate references
   - Clean reference namespace in database

2. **Proper Database Timing**
   - References are committed and indexed before BAM import
   - 30-second indexing timeout plus 2-second buffer
   - 5-second delays between phases in workflow

3. **Enhanced Reference Resolution**
   - BAM import now searches existing database for references
   - Better logging for troubleshooting name mismatches
   - References are guaranteed to be available and searchable

4. **Backward Compatibility**
   - Original `importResults()` method still works (uses `ImportPhase.ALL`)
   - Existing code doesn't need to change
   - New phase-based approach is opt-in

## Testing and Validation

### To Test the Fix:

1. **Use Phase-Based Import:**
   ```java
   PhaseBasedImportWorkflow.executePhaseBasedImport(outputDir, progressListener);
   ```

2. **Check Logs for:**
   - "Phase 1 Complete: Persisted X reference sequences to Geneious database"
   - "Database indexing completed successfully"
   - "Phase 3 Complete: Successfully imported X BAM documents"

3. **Manual Verification:**
   - References should appear in Geneious after Phase 1
   - BAM files should import with proper reference alignment after Phase 3
   - No duplicate reference entries

### Expected Improvements:

- ✅ BAM files should now import programmatically
- ✅ References should be found automatically
- ✅ No more double-import of reference sequences
- ✅ Better error reporting and troubleshooting

## Troubleshooting

If BAM import still fails after Phase 3:

1. **Check Reference Names:**
   - Look for exact name matches between BAM headers and imported references
   - Check logs for "Available ref:" entries

2. **Verify Database State:**
   - Ensure references have URNs assigned
   - Check that indexing completed successfully

3. **Manual Import Test:**
   - Try manually importing a single BAM file through Geneious UI
   - If manual works but programmatic doesn't, increase wait times

4. **Fallback Options:**
   - Use individual phase imports with longer delays
   - Import references manually, then run Phase 3 only

This solution addresses the core timing and double-import issues that prevented BAM files from finding their reference sequences programmatically.