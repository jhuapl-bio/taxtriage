package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentImportException;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import com.biomatters.geneious.publicapi.plugin.DocumentFileImporter;
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
 * Utility class for importing BAM files into Geneious programmatically.
 * This class provides multiple strategies for importing BAM files,
 * leveraging Geneious's built-in SAM/BAM import capabilities.
 */
public class BamFileImporter {

    private static final Logger logger = Logger.getLogger(BamFileImporter.class.getName());

    // Known IDs for SAM/BAM importers in different Geneious versions
    private static final String[] KNOWN_SAM_IMPORTER_IDS = {
        "com.biomatters.plugins.fileimportexport.samimporter.SamImporterPlugin",
        "com.biomatters.plugins.samImporter",
        "SAM/BAM Importer"
    };

    /**
     * Imports a BAM file using the best available method.
     *
     * @param bamFile The BAM file to import
     * @param progressListener Progress listener for the import operation (can be null or ProgressListener.EMPTY)
     * @return List of imported documents, or empty list if import fails
     */
    public static List<AnnotatedPluginDocument> importBamFile(File bamFile, ProgressListener progressListener) {
        return importBamFile(bamFile, progressListener, null);
    }

    /**
     * Imports a BAM file using the best available method, with optional reference documents.
     *
     * @param bamFile The BAM file to import
     * @param progressListener Progress listener for the import operation (can be null or ProgressListener.EMPTY)
     * @param referenceDocuments Optional list of reference sequence documents (can be null)
     * @return List of imported documents, or empty list if import fails
     */
    public static List<AnnotatedPluginDocument> importBamFile(File bamFile, ProgressListener progressListener,
                                                              List<AnnotatedPluginDocument> referenceDocuments) {
        if (!bamFile.exists()) {
            logger.warning("BAM file does not exist: " + bamFile.getAbsolutePath());
            return new ArrayList<>();
        }

        // Use EMPTY listener if null provided
        if (progressListener == null) {
            progressListener = ProgressListener.EMPTY;
        }

        // Log if reference documents are available
        if (referenceDocuments != null && !referenceDocuments.isEmpty()) {
            logger.info("Have " + referenceDocuments.size() + " reference documents available for BAM import");
            for (AnnotatedPluginDocument ref : referenceDocuments) {
                logger.info("  Reference: " + ref.getName());
            }
        }

        // First, ensure the BAM file has an index
        boolean hasIndex = BamIndexer.ensureIndexExists(bamFile);
        if (!hasIndex) {
            logger.warning("Could not ensure BAM index exists for: " + bamFile.getName());
            // Continue anyway - some versions of Geneious may not require an index
        }

        // Try different import strategies in order of preference

        // If we have reference documents, use BamImportHelper's reference-aware import
        if (referenceDocuments != null && !referenceDocuments.isEmpty()) {
            logger.info("Using reference-aware BAM import with " + referenceDocuments.size() + " reference(s)");
            List<AnnotatedPluginDocument> result = BamImportHelper.importBamWithReferences(
                bamFile, referenceDocuments, progressListener);
            if (!result.isEmpty()) {
                logger.info("Successfully imported BAM with reference linking");
                return result;
            }
            logger.warning("Reference-aware import failed, falling back to standard import");
        }

        // Strategy 1: Use PluginUtilities.importDocuments (simplest and most reliable)
        List<AnnotatedPluginDocument> result = importUsingPluginUtilities(bamFile, progressListener);
        if (!result.isEmpty()) {
            return result;
        }

        // Strategy 2: Try to get a specific SAM/BAM importer
        result = importUsingSpecificImporter(bamFile, progressListener);
        if (!result.isEmpty()) {
            return result;
        }

        // Strategy 3: Let Geneious auto-detect the importer
        result = importUsingAutoDetection(bamFile, progressListener);
        if (!result.isEmpty()) {
            return result;
        }

        logger.warning("All BAM import strategies failed for: " + bamFile.getName());
        return new ArrayList<>();
    }

    /**
     * Import using the basic PluginUtilities.importDocuments method.
     */
    private static List<AnnotatedPluginDocument> importUsingPluginUtilities(
            File bamFile, ProgressListener progressListener) {
        try {
            System.out.println("========== BAM IMPORT ATTEMPT ==========");
            System.out.println("Method: PluginUtilities.importDocuments");
            System.out.println("BAM file: " + bamFile.getAbsolutePath());
            System.out.println("File exists: " + bamFile.exists());
            System.out.println("File size: " + bamFile.length() + " bytes");
            System.out.println("Has index (.bai): " + new File(bamFile.getAbsolutePath() + ".bai").exists());

            logger.info("========== BAM IMPORT ATTEMPT ==========");
            logger.info("Method: PluginUtilities.importDocuments");
            logger.info("BAM file: " + bamFile.getAbsolutePath());
            logger.info("File exists: " + bamFile.exists());
            logger.info("File size: " + bamFile.length() + " bytes");
            logger.info("Has index (.bai): " + new File(bamFile.getAbsolutePath() + ".bai").exists());

            // Log BAM header info
            try {
                List<BamReferenceExtractor.ReferenceInfo> bamRefs = BamReferenceExtractor.extractReferences(bamFile);
                System.out.println("BAM contains " + bamRefs.size() + " reference(s):");
                logger.info("BAM contains " + bamRefs.size() + " reference(s):");
                for (BamReferenceExtractor.ReferenceInfo ref : bamRefs) {
                    System.out.println("  - Reference: '" + ref.name + "' (length: " + ref.length + ", accession: " + ref.accession + ")");
                    logger.info("  - Reference: '" + ref.name + "' (length: " + ref.length + ", accession: " + ref.accession + ")");
                }
            } catch (Exception e) {
                System.out.println("Could not extract BAM references for logging: " + e.getMessage());
                logger.warning("Could not extract BAM references for logging: " + e.getMessage());
            }

            // Use the simplest import method
            System.out.println("Calling PluginUtilities.importDocuments...");
            logger.info("Calling PluginUtilities.importDocuments...");
            List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(bamFile, progressListener);

            if (docs != null && !docs.isEmpty()) {
                System.out.println("SUCCESS: Imported " + docs.size() + " document(s) from BAM file");
                logger.info("SUCCESS: Imported " + docs.size() + " document(s) from BAM file");
                for (AnnotatedPluginDocument doc : docs) {
                    System.out.println("  - Imported: " + doc.getName() + " (Type: " + doc.getDocumentClass() + ")");
                    logger.info("  - Imported: " + doc.getName() + " (Type: " + doc.getDocumentClass() + ")");
                }
                return docs;
            } else {
                System.out.println("FAILED: PluginUtilities.importDocuments returned null or empty");
                logger.warning("FAILED: PluginUtilities.importDocuments returned null or empty");
            }

        } catch (IOException e) {
            logger.log(Level.WARNING, "IOException during BAM import: " + e.getMessage(), e);
        } catch (DocumentImportException e) {
            logger.log(Level.WARNING, "DocumentImportException during BAM import: " + e.getMessage(), e);
        } catch (Exception e) {
            logger.log(Level.WARNING, "Unexpected error during BAM import: " + e.getMessage(), e);
        }

        logger.info("========== END BAM IMPORT ATTEMPT ==========");
        return new ArrayList<>();
    }

    /**
     * Import using a specific SAM/BAM importer if available.
     */
    private static List<AnnotatedPluginDocument> importUsingSpecificImporter(
            File bamFile, ProgressListener progressListener) {

        // Try to find a SAM/BAM importer by known IDs
        for (String importerId : KNOWN_SAM_IMPORTER_IDS) {
            try {
                logger.info("Looking for SAM/BAM importer with ID: " + importerId);
                DocumentFileImporter importer = PluginUtilities.getDocumentFileImporter(importerId);

                if (importer != null) {
                    logger.info("Found SAM/BAM importer, attempting import");

                    // Create import options if needed
                    Map<String, String> options = new HashMap<>();
                    // Add any specific import options here if needed

                    List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(
                        bamFile, options, progressListener);

                    if (docs != null && !docs.isEmpty()) {
                        logger.info("Successfully imported " + docs.size() + " document(s) using specific importer");
                        return docs;
                    }
                }
            } catch (Exception e) {
                logger.log(Level.FINE, "Could not use importer ID: " + importerId, e);
            }
        }

        return new ArrayList<>();
    }

    /**
     * Import using auto-detection to find the appropriate importer.
     */
    private static List<AnnotatedPluginDocument> importUsingAutoDetection(
            File bamFile, ProgressListener progressListener) {
        try {
            logger.info("Attempting BAM import using auto-detection for: " + bamFile.getName());

            // Let Geneious auto-detect the appropriate importer
            DocumentFileImporter importer = PluginUtilities.getDocumentFileImporterForSingleFile(bamFile);

            if (importer != null) {
                logger.info("Auto-detected importer: " + importer.getClass().getName());

                // Use the auto-detected importer
                List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(bamFile, progressListener);

                if (docs != null && !docs.isEmpty()) {
                    logger.info("Successfully imported " + docs.size() + " document(s) using auto-detected importer");
                    return docs;
                }
            } else {
                logger.warning("No suitable importer auto-detected for BAM file");
            }

        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to import BAM using auto-detection", e);
        }

        return new ArrayList<>();
    }

    /**
     * Checks if BAM import is supported in the current Geneious installation.
     *
     * @return true if BAM import appears to be supported
     */
    public static boolean isBamImportSupported() {
        // Check if we can find a SAM/BAM importer
        for (String importerId : KNOWN_SAM_IMPORTER_IDS) {
            try {
                DocumentFileImporter importer = PluginUtilities.getDocumentFileImporter(importerId);
                if (importer != null) {
                    return true;
                }
            } catch (Exception e) {
                // Continue checking other IDs
            }
        }

        // Check if a test BAM file would be recognized
        try {
            File testFile = new File("test.bam");
            DocumentFileImporter importer = PluginUtilities.getDocumentFileImporterForSingleFile(testFile);
            if (importer != null) {
                return true;
            }
        } catch (Exception e) {
            // Ignore
        }

        logger.warning("BAM import support not detected in current Geneious installation");
        return false;
    }

    /**
     * Gets information about available BAM import capabilities.
     *
     * @return A string describing the BAM import capabilities
     */
    public static String getBamImportInfo() {
        StringBuilder info = new StringBuilder();
        info.append("BAM Import Capabilities:\n");

        // Check for samtools
        if (BamIndexer.isSamtoolsAvailable()) {
            info.append("- Samtools available for indexing\n");
        } else {
            info.append("- Samtools not found (indexing may fail)\n");
        }

        // Check for BAM importer
        boolean foundImporter = false;
        for (String importerId : KNOWN_SAM_IMPORTER_IDS) {
            try {
                DocumentFileImporter importer = PluginUtilities.getDocumentFileImporter(importerId);
                if (importer != null) {
                    info.append("- Found SAM/BAM importer: ").append(importerId).append("\n");
                    foundImporter = true;
                    break;
                }
            } catch (Exception e) {
                // Continue
            }
        }

        if (!foundImporter) {
            info.append("- No specific SAM/BAM importer found\n");
        }

        return info.toString();
    }
}