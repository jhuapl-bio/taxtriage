package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentImportException;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Helper class to import documents into organized folder structure within Geneious.
 * Creates subfolders for different types of data and imports documents to appropriate locations.
 */
public class FolderBasedImporter {

    private static final Logger logger = Logger.getLogger(FolderBasedImporter.class.getName());

    private WritableDatabaseService rootService;
    private Map<String, WritableDatabaseService> folderServices;
    private String runName;

    /**
     * Creates a new folder-based importer.
     * @param runName Name for the main folder (e.g., "TaxTriage_Run_20240916")
     */
    public FolderBasedImporter(String runName) {
        this.runName = runName;
        this.folderServices = new HashMap<>();
        initializeRootService();
    }

    /**
     * Initializes the root database service.
     */
    private void initializeRootService() {
        try {
            // Get the local database service
            Object service = PluginUtilities.getGeneiousService(PluginUtilities.LOCAL_DATABASE_SERVICE_UNIQUE_ID);
            if (service instanceof WritableDatabaseService) {
                this.rootService = (WritableDatabaseService) service;
                logger.info("Got local database service");
            } else {
                // Fallback to generic WritableDatabaseService
                service = PluginUtilities.getGeneiousService("WritableDatabaseService");
                if (service instanceof WritableDatabaseService) {
                    this.rootService = (WritableDatabaseService) service;
                    logger.info("Got generic writable database service");
                }
            }
        } catch (Exception e) {
            logger.log(Level.WARNING, "Could not get database service", e);
        }
    }

    /**
     * Creates the folder structure for TaxTriage results.
     * @param progressListener Progress tracking
     * @return true if folder structure was created successfully
     */
    public boolean createFolderStructure(ProgressListener progressListener) {
        if (rootService == null) {
            logger.warning("No database service available - cannot create folder structure");
            return false;
        }

        try {
            logger.info("Creating folder structure for: " + runName);
            System.out.println("===== CREATING FOLDER STRUCTURE =====");
            System.out.println("Main folder: " + runName);

            // Create main folder for this TaxTriage run
            WritableDatabaseService mainFolder = rootService.createChildFolder(runName);
            if (mainFolder == null) {
                // Folder might already exist, try to get it
                mainFolder = rootService.getChildService(runName);
            }

            if (mainFolder == null) {
                logger.warning("Could not create or access main folder: " + runName);
                return false;
            }

            folderServices.put("main", mainFolder);

            // Create subfolders for different data types
            String[] subfolders = {
                "References",      // For GenBank reference sequences
                "Text_Reports",    // For text/CSV reports
                "Alignments",      // For BAM alignments
                "Kraken_Results",  // For Kraken classification results
                "FASTQ_Files",     // For FASTQ read files
                "Consensus",       // For consensus sequences
                "VCF_Files"        // For variant call files
            };

            for (String folderName : subfolders) {
                try {
                    WritableDatabaseService subfolder = mainFolder.createChildFolder(folderName);
                    if (subfolder == null) {
                        // Folder might already exist
                        subfolder = mainFolder.getChildService(folderName);
                    }

                    if (subfolder != null) {
                        folderServices.put(folderName, subfolder);
                        System.out.println("  Created subfolder: " + folderName);
                    } else {
                        System.out.println("  Could not create subfolder: " + folderName);
                    }
                } catch (Exception e) {
                    logger.warning("Error creating subfolder " + folderName + ": " + e.getMessage());
                }
            }

            System.out.println("=====================================\n");
            return true;

        } catch (Exception e) {
            logger.log(Level.SEVERE, "Failed to create folder structure", e);
            return false;
        }
    }

    /**
     * Imports documents to a specific subfolder.
     * @param file File to import
     * @param folderName Name of the subfolder (e.g., "References", "Alignments")
     * @param progressListener Progress tracking
     * @return List of imported documents
     */
    public List<AnnotatedPluginDocument> importToFolder(File file, String folderName,
                                                         ProgressListener progressListener)
                                                         throws IOException, DocumentImportException {

        WritableDatabaseService targetFolder = folderServices.get(folderName);

        if (targetFolder == null) {
            logger.warning("Target folder not found: " + folderName + ". Using main folder.");
            targetFolder = folderServices.get("main");
        }

        if (targetFolder == null) {
            logger.warning("No folder available. Falling back to standard import.");
            return PluginUtilities.importDocuments(file, progressListener);
        }

        logger.info("Importing " + file.getName() + " to folder: " + folderName);
        System.out.println("Importing to " + folderName + ": " + file.getName());

        try {
            // Import documents to the specific database folder
            List<AnnotatedPluginDocument> imported = PluginUtilities.importDocumentsToDatabase(
                file, targetFolder, progressListener);

            if (imported != null && !imported.isEmpty()) {
                logger.info("Successfully imported " + imported.size() + " document(s) to " + folderName);
                return imported;
            } else {
                logger.warning("No documents imported from: " + file.getName());
                return new ArrayList<>();
            }

        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to import to folder: " + folderName, e);
            throw new DocumentImportException("Failed to import " + file.getName() + " to " + folderName, e);
        }
    }

    /**
     * Imports multiple files to a specific subfolder.
     * @param files Files to import
     * @param folderName Name of the subfolder
     * @param progressListener Progress tracking
     * @return Combined list of all imported documents
     */
    public List<AnnotatedPluginDocument> importMultipleToFolder(List<File> files, String folderName,
                                                                ProgressListener progressListener)
                                                                throws IOException, DocumentImportException {

        List<AnnotatedPluginDocument> allImported = new ArrayList<>();

        int fileCount = 0;
        for (File file : files) {
            if (progressListener != null) {
                progressListener.setMessage("Importing " + file.getName() + " to " + folderName);
                progressListener.setProgress((double) fileCount / files.size());
            }

            List<AnnotatedPluginDocument> imported = importToFolder(file, folderName, progressListener);
            allImported.addAll(imported);
            fileCount++;
        }

        return allImported;
    }

    /**
     * Gets the reference documents that were imported to the References folder.
     * @return List of reference documents
     */
    public List<AnnotatedPluginDocument> getReferenceDocuments() {
        WritableDatabaseService refFolder = folderServices.get("References");
        if (refFolder == null) {
            return new ArrayList<>();
        }

        try {
            // Get all documents from the References folder
            // Note: This might require using the retrieve method with an empty query
            // The exact method depends on the Geneious API version
            return new ArrayList<>(); // Placeholder - needs implementation based on API
        } catch (Exception e) {
            logger.log(Level.WARNING, "Could not retrieve reference documents", e);
            return new ArrayList<>();
        }
    }

    /**
     * Imports BAM files with their associated references.
     * References should already be imported to the References folder.
     * @param bamFile BAM file to import
     * @param progressListener Progress tracking
     * @return List of imported BAM documents
     */
    public List<AnnotatedPluginDocument> importBamWithReferences(File bamFile,
                                                                 ProgressListener progressListener)
                                                                 throws IOException, DocumentImportException {

        // First ensure references are available
        List<AnnotatedPluginDocument> references = getReferenceDocuments();

        logger.info("Importing BAM file: " + bamFile.getName());
        logger.info("Available references: " + references.size());

        // Import the BAM to the Alignments folder
        return importToFolder(bamFile, "Alignments", progressListener);
    }

    /**
     * Gets the main folder service.
     */
    public WritableDatabaseService getMainFolder() {
        return folderServices.get("main");
    }

    /**
     * Gets a specific subfolder service.
     */
    public WritableDatabaseService getSubfolder(String folderName) {
        return folderServices.get(folderName);
    }

    /**
     * Checks if the folder structure has been created.
     */
    public boolean isFolderStructureReady() {
        return !folderServices.isEmpty() && folderServices.get("main") != null;
    }
}