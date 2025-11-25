package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.documents.PluginDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentFileImporter;
import com.biomatters.geneious.publicapi.plugin.DocumentOperationException;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import jebl.util.ProgressListener;

import java.io.File;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Proper BAM importer that uses the DocumentFileImporter API directly
 * to import BAM files with their reference sequences.
 *
 * This implementation follows the guidance from the Geneious developer
 * to pass both BAM and reference files to the importer's
 * importDocumentsFromMultipleFilesReturningUnimported() method.
 */
public class ProperBamImporter {

    private static final Logger logger = Logger.getLogger(ProperBamImporter.class.getName());

    /**
     * Import a BAM file along with its reference sequences (GenBank files)
     * into the specified database folder.
     *
     * @param bamFile The BAM file to import
     * @param referenceFiles List of GenBank reference files (can be empty)
     * @param destination The database folder to import into
     * @param progressListener Progress tracking
     * @return List of imported documents (alignments)
     * @throws DocumentOperationException if import fails
     */
    public static List<AnnotatedPluginDocument> importBamWithReferences(
            File bamFile,
            List<File> referenceFiles,
            WritableDatabaseService destination,
            ProgressListener progressListener) throws DocumentOperationException {

        if (bamFile == null || !bamFile.exists()) {
            throw new DocumentOperationException("BAM file does not exist: " + bamFile);
        }

        if (destination == null) {
            throw new DocumentOperationException("Destination folder is null");
        }

        System.out.println("\n=== ProperBamImporter: Starting BAM import ===");
        System.out.println("  BAM file: " + bamFile.getName());
        System.out.println("  Reference files: " + referenceFiles.size());
        for (File ref : referenceFiles) {
            System.out.println("    - " + ref.getName());
        }

        // Get the BAM importer
        DocumentFileImporter importer = PluginUtilities.getDocumentFileImporters().stream()
            .filter(imp ->
                Arrays.stream(imp.getPermissibleExtensions()).anyMatch(extension ->
                    extension.equalsIgnoreCase(".bam"))
            ).findFirst()
            .orElseThrow(() -> new DocumentOperationException("No BAM importer available"));

        System.out.println("  Found BAM importer: " + importer.getClass().getSimpleName());

        // Prepare the list of files to import
        // The BAM importer expects both the BAM file and reference files together
        List<File> files = new ArrayList<>();

        // Add reference files first (important for the importer to find them)
        files.addAll(referenceFiles);

        // Add the BAM file
        files.add(bamFile);

        System.out.println("  Total files to import together: " + files.size());

        // Create an ImportCallback that stores the results
        List<AnnotatedPluginDocument> imported = new ArrayList<>();
        DocumentFileImporter.ImportCallback callback = new DocumentFileImporter.ImportCallback() {
            @Override
            public AnnotatedPluginDocument addDocument(PluginDocument document) {
                AnnotatedPluginDocument apd = DocumentUtilities.createAnnotatedPluginDocument(document);
                imported.add(apd);
                System.out.println("    Imported document: " + apd.getName());
                return apd;
            }

            @Override
            public AnnotatedPluginDocument addDocument(AnnotatedPluginDocument annotatedDocument) {
                imported.add(annotatedDocument);
                System.out.println("    Imported document: " + annotatedDocument.getName());
                return annotatedDocument;
            }
        };

        try {
            // Call the importer's method directly with all files
            // This is the key - passing both BAM and reference files together
            System.out.println("  Calling importDocumentsFromMultipleFilesReturningUnimported...");

            // According to the developer's example:
            // importer.importDocumentsFromMultipleFilesReturningUnimported(
            //     importer.getOptions(files, ProgressListener.EMPTY),
            //     files,
            //     callback,
            //     ProgressListener.EMPTY)

            List<File> unused = importer.importDocumentsFromMultipleFilesReturningUnimported(
                importer.getOptions(files, progressListener != null ? progressListener : ProgressListener.EMPTY),
                files,
                callback,
                progressListener != null ? progressListener : ProgressListener.EMPTY
            );

            if (!unused.isEmpty()) {
                System.out.println("  Warning: " + unused.size() + " file(s) were not imported:");
                for (File f : unused) {
                    System.out.println("    - " + f.getName());
                }
            }

            // Now add the imported documents to the destination folder
            System.out.println("  Adding " + imported.size() + " document(s) to destination folder...");
            List<AnnotatedPluginDocument> finalDocuments = new ArrayList<>();

            for (AnnotatedPluginDocument toAdd : imported) {
                try {
                    AnnotatedPluginDocument added = destination.addDocumentCopy(
                        toAdd,
                        progressListener != null ? progressListener : ProgressListener.EMPTY
                    );
                    if (added != null) {
                        finalDocuments.add(added);
                        System.out.println("    Added to folder: " + added.getName());
                    }
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Failed to add document to folder: " + toAdd.getName(), e);
                    System.out.println("    Failed to add: " + toAdd.getName() + " - " + e.getMessage());
                }
            }

            System.out.println("  Successfully imported " + finalDocuments.size() + " document(s)");
            System.out.println("=== ProperBamImporter: Import complete ===\n");

            return finalDocuments;

        } catch (Exception e) {
            logger.log(Level.SEVERE, "BAM import failed", e);
            throw new DocumentOperationException("Failed to import BAM file: " + e.getMessage(), e);
        }
    }

    /**
     * Import a BAM file without explicit reference files.
     * The importer will try to find references in the database.
     *
     * @param bamFile The BAM file to import
     * @param destination The database folder to import into
     * @param progressListener Progress tracking
     * @return List of imported documents (alignments)
     * @throws DocumentOperationException if import fails
     */
    public static List<AnnotatedPluginDocument> importBam(
            File bamFile,
            WritableDatabaseService destination,
            ProgressListener progressListener) throws DocumentOperationException {

        return importBamWithReferences(bamFile, new ArrayList<>(), destination, progressListener);
    }
}