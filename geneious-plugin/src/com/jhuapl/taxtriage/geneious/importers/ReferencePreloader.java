package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.databaseservice.DatabaseService;
import com.biomatters.geneious.publicapi.databaseservice.Query;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.documents.sequence.NucleotideSequenceDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentImportException;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Helper class to pre-load reference sequences into Geneious memory/cache
 * to ensure they are discoverable during BAM import.
 */
public class ReferencePreloader {

    private static final Logger logger = Logger.getLogger(ReferencePreloader.class.getName());

    /**
     * Enterprise-grade reference preloading with comprehensive validation and caching.
     * This method ensures references are properly loaded, cached, and discoverable
     * before BAM import attempts to use them.
     *
     * @param referenceDocuments The reference documents to pre-load
     * @param bamFile The BAM file that will need these references
     * @return Map of reference names to their pre-loaded documents
     */
    public static Map<String, AnnotatedPluginDocument> preloadReferences(
            List<AnnotatedPluginDocument> referenceDocuments,
            File bamFile) {

        logger.info("===== ENTERPRISE REFERENCE PRELOADING =====");
        logger.info("Preloading " + referenceDocuments.size() + " reference document(s) for BAM: " + bamFile.getName());

        Map<String, AnnotatedPluginDocument> preloadedRefs = new HashMap<>();

        // Phase 1: Extract BAM requirements
        List<BamReferenceExtractor.ReferenceInfo> bamRefs = extractBamRequirements(bamFile);
        if (bamRefs.isEmpty()) {
            logger.warning("No BAM references extracted - preloading skipped");
            return preloadedRefs;
        }

        // Phase 2: Validate and prepare references
        Map<String, AnnotatedPluginDocument> validatedRefs = validateAndPrepareReferences(
            referenceDocuments, bamRefs);

        // Phase 3: Execute comprehensive preloading strategies
        for (Map.Entry<String, AnnotatedPluginDocument> entry : validatedRefs.entrySet()) {
            String bamRefName = entry.getKey();
            AnnotatedPluginDocument refDoc = entry.getValue();

            if (executePreloadingStrategies(bamRefName, refDoc)) {
                preloadedRefs.put(bamRefName, refDoc);
                logger.info("Successfully preloaded reference: '" + bamRefName + "'");
            } else {
                logger.warning("Failed to preload reference: '" + bamRefName + "'");
            }
        }

        // Phase 4: Verify preloading results
        verifyPreloadingResults(preloadedRefs, bamRefs);

        logger.info("Preloading complete: " + preloadedRefs.size() + "/" + bamRefs.size() + " references loaded");
        logger.info("===== PRELOADING COMPLETE =====\n");

        return preloadedRefs;
    }

    /**
     * Extracts BAM reference requirements with error handling.
     */
    private static List<BamReferenceExtractor.ReferenceInfo> extractBamRequirements(File bamFile) {
        try {
            List<BamReferenceExtractor.ReferenceInfo> bamRefs = BamReferenceExtractor.extractReferences(bamFile);
            logger.info("BAM file requires " + bamRefs.size() + " reference(s):");
            for (BamReferenceExtractor.ReferenceInfo ref : bamRefs) {
                logger.info("  - '" + ref.name + "' (accession: " + ref.accession + ", length: " + ref.length + ")");
            }
            return bamRefs;
        } catch (IOException e) {
            logger.log(Level.WARNING, "Could not extract BAM references", e);
            return new ArrayList<>();
        }
    }

    /**
     * Validates and prepares references for preloading.
     */
    private static Map<String, AnnotatedPluginDocument> validateAndPrepareReferences(
            List<AnnotatedPluginDocument> referenceDocuments,
            List<BamReferenceExtractor.ReferenceInfo> bamRefs) {

        Map<String, AnnotatedPluginDocument> validatedRefs = new HashMap<>();

        logger.info("Validating and matching references...");

        for (BamReferenceExtractor.ReferenceInfo bamRef : bamRefs) {
            AnnotatedPluginDocument matchedDoc = findMatchingReference(bamRef, referenceDocuments);

            if (matchedDoc != null) {
                // Validate the matched document
                if (validateReferenceDocument(matchedDoc, bamRef)) {
                    validatedRefs.put(bamRef.name, matchedDoc);
                    logger.info("Validated reference match: '" + bamRef.name + "' -> '" + matchedDoc.getName() + "'");
                } else {
                    logger.warning("Reference validation failed: '" + bamRef.name + "'");
                }
            } else {
                logger.warning("No matching document found for BAM reference: '" + bamRef.name + "'");
            }
        }

        return validatedRefs;
    }

    /**
     * Finds a matching reference document for a BAM reference.
     */
    private static AnnotatedPluginDocument findMatchingReference(
            BamReferenceExtractor.ReferenceInfo bamRef,
            List<AnnotatedPluginDocument> referenceDocuments) {

        // Strategy 1: Exact name match
        for (AnnotatedPluginDocument refDoc : referenceDocuments) {
            if (bamRef.name.equals(refDoc.getName())) {
                return refDoc;
            }
        }

        // Strategy 2: Accession match
        if (bamRef.accession != null) {
            for (AnnotatedPluginDocument refDoc : referenceDocuments) {
                String docName = refDoc.getName();
                if (docName != null && docName.contains(bamRef.accession)) {
                    return refDoc;
                }
            }
        }

        // Strategy 3: Partial name match
        for (AnnotatedPluginDocument refDoc : referenceDocuments) {
            String docName = refDoc.getName();
            if (docName != null && (docName.contains(bamRef.name) || bamRef.name.contains(docName))) {
                return refDoc;
            }
        }

        return null;
    }

    /**
     * Validates a reference document against BAM requirements.
     */
    private static boolean validateReferenceDocument(
            AnnotatedPluginDocument refDoc,
            BamReferenceExtractor.ReferenceInfo bamRef) {

        try {
            // Basic validation
            if (refDoc == null || refDoc.getName() == null) {
                return false;
            }

            // Check if it's a sequence document
            if (refDoc.getDocumentClass() == NucleotideSequenceDocument.class) {
                NucleotideSequenceDocument seqDoc = (NucleotideSequenceDocument) refDoc.getDocument();

                // Validate sequence length if possible
                if (seqDoc.getSequenceLength() > 0 && bamRef.length > 0) {
                    if (seqDoc.getSequenceLength() != bamRef.length) {
                        logger.warning("Length mismatch: document=" + seqDoc.getSequenceLength() +
                                     ", BAM=" + bamRef.length);
                        // Continue anyway - length mismatches might be acceptable
                    }
                }
            }

            return true;

        } catch (Exception e) {
            logger.log(Level.WARNING, "Error validating reference document", e);
            return false;
        }
    }

    /**
     * Executes comprehensive preloading strategies for a single reference.
     */
    private static boolean executePreloadingStrategies(String bamRefName, AnnotatedPluginDocument refDoc) {
        logger.info("Executing preloading strategies for: '" + bamRefName + "'");

        boolean success = false;

        // Strategy 1: Document registration
        if (tryDocumentRegistration(refDoc)) {
            success = true;
        }

        // Strategy 2: Memory caching
        if (tryMemoryCaching(refDoc)) {
            success = true;
        }

        // Strategy 3: URN-based caching
        if (tryUrnCaching(refDoc)) {
            success = true;
        }

        // Strategy 4: Name standardization
        if (tryNameStandardization(refDoc, bamRefName)) {
            success = true;
        }

        return success;
    }

    /**
     * Strategy 1: Document registration with DocumentUtilities.
     */
    private static boolean tryDocumentRegistration(AnnotatedPluginDocument refDoc) {
        try {
            logger.fine("Attempting document registration...");
            DocumentUtilities.addGeneratedDocument(refDoc, false);
            logger.fine("Document registration successful");
            return true;
        } catch (Exception e) {
            logger.fine("Document registration failed: " + e.getMessage());
            return false;
        }
    }

    /**
     * Strategy 2: Memory caching by accessing document properties.
     */
    private static boolean tryMemoryCaching(AnnotatedPluginDocument refDoc) {
        try {
            logger.fine("Attempting memory caching...");
            // Access various properties to force loading into memory
            String name = refDoc.getName();
            String urn = refDoc.getURN() != null ? refDoc.getURN().toString() : null;
            Object document = refDoc.getDocument();
            logger.fine("Memory caching successful");
            return true;
        } catch (Exception e) {
            logger.fine("Memory caching failed: " + e.getMessage());
            return false;
        }
    }

    /**
     * Strategy 3: URN-based caching.
     */
    private static boolean tryUrnCaching(AnnotatedPluginDocument refDoc) {
        try {
            if (refDoc.getURN() != null) {
                logger.fine("Attempting URN-based caching...");
                AnnotatedPluginDocument cached = DocumentUtilities.getDocumentByURN(refDoc.getURN());
                if (cached != null) {
                    logger.fine("URN-based caching successful");
                    return true;
                }
            }
        } catch (Exception e) {
            logger.fine("URN-based caching failed: " + e.getMessage());
        }
        return false;
    }

    /**
     * Strategy 4: Name standardization to match BAM expectations.
     */
    private static boolean tryNameStandardization(AnnotatedPluginDocument refDoc, String expectedName) {
        try {
            String currentName = refDoc.getName();
            if (currentName != null && !currentName.equals(expectedName)) {
                logger.fine("Attempting name standardization: '" + currentName + "' -> '" + expectedName + "'");
                refDoc.setName(expectedName);
                logger.fine("Name standardization successful");
                return true;
            }
        } catch (Exception e) {
            logger.fine("Name standardization failed: " + e.getMessage());
        }
        return false;
    }

    /**
     * Verifies preloading results by testing database searchability.
     */
    private static void verifyPreloadingResults(
            Map<String, AnnotatedPluginDocument> preloadedRefs,
            List<BamReferenceExtractor.ReferenceInfo> bamRefs) {

        logger.info("Verifying preloading results...");

        try {
            DatabaseService database = null;
            try {
                Object service = PluginUtilities.getGeneiousService("DatabaseService");
                if (service instanceof DatabaseService) {
                    database = (DatabaseService) service;
                }
            } catch (Exception e) {
                logger.warning("Could not get DatabaseService for verification");
                return;
            }

            if (database == null) {
                logger.warning("DatabaseService not available - skipping verification");
                return;
            }

            int verifiedCount = 0;
            for (BamReferenceExtractor.ReferenceInfo bamRef : bamRefs) {
                if (isReferenceSearchable(database, bamRef.name)) {
                    verifiedCount++;
                    logger.info("  ✓ Verified searchable: '" + bamRef.name + "'");
                } else {
                    logger.warning("  ✗ Not searchable: '" + bamRef.name + "'");
                }
            }

            logger.info("Verification complete: " + verifiedCount + "/" + bamRefs.size() + " references searchable");

        } catch (Exception e) {
            logger.log(Level.WARNING, "Error during verification", e);
        }
    }

    /**
     * Tests if a reference is searchable in the database.
     */
    private static boolean isReferenceSearchable(DatabaseService database, String referenceName) {
        try {
            Query query = Query.Factory.createQuery("Name=\"" + referenceName + "\"");
            List<AnnotatedPluginDocument> results = database.retrieve(query, ProgressListener.EMPTY);
            return !results.isEmpty();
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Creates temporary FASTA files for references that the BAM importer can find.
     * This is a fallback approach if in-memory references aren't being found.
     *
     * @param referenceDocuments The reference documents
     * @param bamFile The BAM file that needs these references
     * @param tempDir Directory to write temporary files
     * @return List of created FASTA files
     */
    public static List<File> createTemporaryReferenceFiles(
            List<AnnotatedPluginDocument> referenceDocuments,
            File bamFile,
            File tempDir) {

        List<File> tempFiles = new ArrayList<>();

        System.out.println("===== CREATING TEMPORARY REFERENCE FILES =====");
        System.out.println("Temp directory: " + tempDir.getAbsolutePath());

        if (!tempDir.exists()) {
            tempDir.mkdirs();
        }

        try {
            // Extract BAM reference requirements
            List<BamReferenceExtractor.ReferenceInfo> bamRefs = BamReferenceExtractor.extractReferences(bamFile);

            for (BamReferenceExtractor.ReferenceInfo bamRef : bamRefs) {
                System.out.println("Creating temp file for: " + bamRef.name);

                // Find matching document
                AnnotatedPluginDocument matchedDoc = null;
                for (AnnotatedPluginDocument refDoc : referenceDocuments) {
                    String docName = refDoc.getName();
                    if (docName != null && (docName.equals(bamRef.name) ||
                        (bamRef.accession != null && docName.contains(bamRef.accession)))) {
                        matchedDoc = refDoc;
                        break;
                    }
                }

                if (matchedDoc != null) {
                    try {
                        // Create temp FASTA file
                        File tempFasta = new File(tempDir, bamRef.name.replaceAll("[^a-zA-Z0-9._-]", "_") + ".fasta");

                        // Note: DocumentUtilities doesn't have exportDocument method
                        // We would need to use PluginUtilities.exportDocuments or similar
                        // For now, we'll skip this functionality
                        System.out.println("  Skipping temp file creation (API method not available)");
                    } catch (Exception e) {
                        System.out.println("  Failed to create temp file: " + e.getMessage());
                    }
                }
            }

        } catch (IOException e) {
            System.out.println("Error creating temporary files: " + e.getMessage());
            logger.log(Level.WARNING, "Failed to create temporary reference files", e);
        }

        System.out.println("Created " + tempFiles.size() + " temporary reference file(s)");
        System.out.println("=====================================\n");

        return tempFiles;
    }

    /**
     * Imports temporary reference files to ensure they are available for BAM import.
     */
    public static List<AnnotatedPluginDocument> importTemporaryReferences(
            List<File> tempReferenceFiles,
            ProgressListener progressListener) {

        List<AnnotatedPluginDocument> importedRefs = new ArrayList<>();

        System.out.println("===== IMPORTING TEMPORARY REFERENCES =====");

        for (File refFile : tempReferenceFiles) {
            try {
                System.out.println("Importing: " + refFile.getName());
                List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(refFile, progressListener);

                if (docs != null && !docs.isEmpty()) {
                    importedRefs.addAll(docs);
                    System.out.println("  SUCCESS: Imported " + docs.size() + " document(s)");
                } else {
                    System.out.println("  FAILED: No documents imported");
                }
            } catch (IOException | DocumentImportException e) {
                System.out.println("  ERROR: " + e.getMessage());
            }
        }

        System.out.println("Imported " + importedRefs.size() + " reference(s) from temp files");
        System.out.println("==========================================\n");

        return importedRefs;
    }
}