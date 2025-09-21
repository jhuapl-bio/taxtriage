package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.URN;
import com.biomatters.geneious.publicapi.plugin.DocumentImportException;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * BAM importer using Geneious API to allow automatic reference detection.
 * Based on the provided example code.
 */
public class GeneiousBamImporter {

    private static final Logger logger = Logger.getLogger(GeneiousBamImporter.class.getName());

    /**
     * Imports a BAM file and allows Geneious to automatically detect reference sequences
     * that are already present in the local database.
     *
     * @param bamFile The BAM file to be imported
     * @param targetFolder The destination folder in Geneious
     * @param progressListener Progress tracking
     * @return List of imported documents
     */
    public static List<AnnotatedPluginDocument> importBamWithAutoReferenceDetection(
            File bamFile,
            WritableDatabaseService targetFolder,
            ProgressListener progressListener) {

        List<AnnotatedPluginDocument> importedDocs = new ArrayList<>();

        logger.info("Importing BAM file with auto-reference detection: " + bamFile.getName());
        System.out.println("Importing BAM file: " + bamFile.getAbsolutePath());

        try {
            // Ensure BAM index exists
            BamIndexer.ensureIndexExists(bamFile);

            // Import the BAM file using the standard Geneious importer
            // This will trigger the search for reference sequences within the database
            if (targetFolder != null) {
                // Import to specific folder
                importedDocs = PluginUtilities.importDocumentsToDatabase(bamFile, targetFolder, progressListener);
            } else {
                // Import to current location
                importedDocs = PluginUtilities.importDocuments(bamFile, progressListener);
            }

            if (!importedDocs.isEmpty()) {
                logger.info("Successfully imported BAM file: " + bamFile.getName() +
                           " (" + importedDocs.size() + " document(s))");
                System.out.println("  Successfully imported BAM with " + importedDocs.size() + " document(s)");

                for (AnnotatedPluginDocument doc : importedDocs) {
                    logger.info("  Imported: " + doc.getName() + " (Type: " + doc.getDocumentClass() + ")");
                }
            } else {
                logger.warning("BAM import returned no documents: " + bamFile.getName());
                System.out.println("  Warning: No documents imported from BAM file");
            }

        } catch (IOException | DocumentImportException e) {
            logger.log(Level.WARNING, "Failed to import BAM file: " + bamFile.getName(), e);
            System.out.println("  Error importing BAM: " + e.getMessage());
        }

        return importedDocs;
    }

    /**
     * Imports BAM file after ensuring references are committed to database.
     *
     * @param bamFile The BAM file
     * @param referenceDocuments Reference documents that should already be in database
     * @param targetFolder Target folder for import
     * @param progressListener Progress tracking
     * @return List of imported documents
     */
    public static List<AnnotatedPluginDocument> importBamWithReferences(
            File bamFile,
            List<AnnotatedPluginDocument> referenceDocuments,
            WritableDatabaseService targetFolder,
            ProgressListener progressListener) {

        logger.info("Importing BAM with " + referenceDocuments.size() + " known references");

        // Log reference information
        for (AnnotatedPluginDocument ref : referenceDocuments) {
            logger.info("  Reference available: " + ref.getName() +
                       " (URN: " + ref.getURN() + ")");
        }

        // Add a delay to ensure database indexing is complete
        try {
            logger.info("Waiting for database indexing...");
            Thread.sleep(3000); // 3 second delay
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }

        // Import using auto-detection (references should be found)
        return importBamWithAutoReferenceDetection(bamFile, targetFolder, progressListener);
    }
}