package com.jhuapl.taxtriage.geneious.importer;

import com.biomatters.geneious.publicapi.databaseservice.DatabaseServiceException;
import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.implementations.TextDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentOperationException;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Builds the folder structure in Geneious to mirror TaxTriage output directories.
 *
 * This class creates nested folders in Geneious that match the directory structure
 * of TaxTriage output, organizing imported documents in a logical hierarchy.
 * It handles special TaxTriage folders and maintains path mappings for efficient
 * document placement.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class FolderStructureBuilder {

    private static final Logger logger = Logger.getLogger(FolderStructureBuilder.class.getName());

    /** Map of directory paths to Geneious folder documents */
    private final Map<String, AnnotatedPluginDocument> folderCache;

    /** The root folder for TaxTriage imports */
    private AnnotatedPluginDocument rootFolder;

    /** The database service for creating folders */
    private WritableDatabaseService databaseService;

    /** Known TaxTriage output directories and their descriptions */
    private static final Map<String, String> TAXTRIAGE_FOLDERS = new HashMap<>();

    static {
        TAXTRIAGE_FOLDERS.put("fastp", "Quality control and preprocessing results");
        TAXTRIAGE_FOLDERS.put("kraken2", "Taxonomic classification results");
        TAXTRIAGE_FOLDERS.put("bracken", "Abundance estimation results");
        TAXTRIAGE_FOLDERS.put("count", "Read counting and statistics");
        TAXTRIAGE_FOLDERS.put("pipeline_info", "Workflow execution information");
        TAXTRIAGE_FOLDERS.put("removetaxidsclassification", "Filtered classification results");
        TAXTRIAGE_FOLDERS.put("multiqc", "Quality control reports");
        TAXTRIAGE_FOLDERS.put("assembly", "Assembly results");
        TAXTRIAGE_FOLDERS.put("annotation", "Genome annotation results");
        TAXTRIAGE_FOLDERS.put("variants", "Variant calling results");
        TAXTRIAGE_FOLDERS.put("reports", "Analysis reports and summaries");
    }

    /**
     * Creates a new folder structure builder.
     *
     * @param databaseService the database service for creating folders
     */
    public FolderStructureBuilder(WritableDatabaseService databaseService) {
        this.databaseService = databaseService;
        this.folderCache = new HashMap<>();
    }

    /**
     * Creates or gets the root folder for TaxTriage imports.
     *
     * @param workflowId the workflow ID for naming the root folder
     * @return the root folder document
     * @throws DatabaseServiceException if folder creation fails
     */
    public AnnotatedPluginDocument createRootFolder(String workflowId) throws DatabaseServiceException {
        if (rootFolder != null) {
            return rootFolder;
        }

        String folderName = "TaxTriage_Results_" + workflowId;
        String description = "TaxTriage workflow results imported on " + new Date();

        rootFolder = createFolderDocument(folderName, description);
        folderCache.put("", rootFolder);

        logger.info("Created root folder: " + folderName);
        return rootFolder;
    }

    /**
     * Creates or gets a folder for the given path relative to the output directory.
     *
     * @param outputDirectory the TaxTriage output directory
     * @param targetFile the file that will be placed in this folder
     * @return the folder document for the file's directory
     * @throws DatabaseServiceException if folder creation fails
     */
    public AnnotatedPluginDocument getOrCreateFolder(File outputDirectory, File targetFile)
            throws DatabaseServiceException {

        if (rootFolder == null) {
            throw new IllegalStateException("Root folder not created. Call createRootFolder() first.");
        }

        // Calculate relative path from output directory to file's parent
        Path outputPath = outputDirectory.toPath().toAbsolutePath();
        Path filePath = targetFile.toPath().toAbsolutePath();
        Path relativePath = outputPath.relativize(filePath.getParent());

        String pathKey = relativePath.toString();

        // Return root folder for files directly in output directory
        if (pathKey.equals(".") || pathKey.isEmpty()) {
            return rootFolder;
        }

        // Check cache first
        if (folderCache.containsKey(pathKey)) {
            return folderCache.get(pathKey);
        }

        // Create folder hierarchy
        return createFolderHierarchy(relativePath);
    }

    /**
     * Creates the complete folder hierarchy for the given path.
     *
     * @param relativePath the relative path from the output directory
     * @return the deepest folder in the hierarchy
     * @throws DatabaseServiceException if folder creation fails
     */
    private AnnotatedPluginDocument createFolderHierarchy(Path relativePath) throws DatabaseServiceException {
        AnnotatedPluginDocument parentFolder = rootFolder;
        StringBuilder currentPath = new StringBuilder();

        for (Path component : relativePath) {
            if (currentPath.length() > 0) {
                currentPath.append(File.separator);
            }
            currentPath.append(component.toString());

            String pathKey = currentPath.toString();

            // Check if folder already exists
            if (folderCache.containsKey(pathKey)) {
                parentFolder = folderCache.get(pathKey);
            } else {
                // Create new folder
                String folderName = component.toString();
                String description = getFolderDescription(folderName);

                AnnotatedPluginDocument newFolder = createFolderDocument(folderName, description);
                folderCache.put(pathKey, newFolder);

                logger.info("Created folder: " + pathKey);
                parentFolder = newFolder;
            }
        }

        return parentFolder;
    }

    /**
     * Creates a folder document with the given name and description.
     *
     * @param name the folder name
     * @param description the folder description
     * @return the created folder document
     */
    private AnnotatedPluginDocument createFolderDocument(String name, String description) {
        // Create a simple document to represent the folder
        TextDocument folderDoc;
        try {
            folderDoc = new TextDocument(name, description, TextDocument.Format.Plain);
        } catch (DocumentOperationException e) {
            throw new RuntimeException("Failed to create folder document", e);
        }

        return DocumentUtilities.createAnnotatedPluginDocument(folderDoc);
    }

    /**
     * Gets a description for a folder based on its name.
     *
     * @param folderName the folder name
     * @return a descriptive string
     */
    private String getFolderDescription(String folderName) {
        // Check for known TaxTriage folders
        if (TAXTRIAGE_FOLDERS.containsKey(folderName.toLowerCase())) {
            return TAXTRIAGE_FOLDERS.get(folderName.toLowerCase());
        }

        // Check for sample folders (usually start with sample name)
        if (folderName.matches("^[A-Za-z0-9_-]+$") && !folderName.toLowerCase().contains("report")) {
            return "Sample-specific analysis results for " + folderName;
        }

        // Default description
        return "Analysis results in " + folderName + " directory";
    }

    /**
     * Gets the total number of folders created.
     *
     * @return the folder count
     */
    public int getFolderCount() {
        return folderCache.size();
    }

    /**
     * Gets all created folder paths.
     *
     * @return array of folder path strings
     */
    public String[] getCreatedFolderPaths() {
        return folderCache.keySet().toArray(new String[0]);
    }

    /**
     * Clears the folder cache and resets the builder.
     */
    public void reset() {
        folderCache.clear();
        rootFolder = null;
    }

    /**
     * Gets the cached folder for a given path.
     *
     * @param relativePath the relative path from output directory
     * @return the folder document, or null if not cached
     */
    public AnnotatedPluginDocument getCachedFolder(String relativePath) {
        return folderCache.get(relativePath);
    }

    /**
     * Checks if a folder exists for the given path.
     *
     * @param relativePath the relative path to check
     * @return true if the folder exists in cache
     */
    public boolean hasCachedFolder(String relativePath) {
        return folderCache.containsKey(relativePath);
    }

    /**
     * Creates a summary of the folder structure created.
     *
     * @return formatted summary string
     */
    public String getFolderStructureSummary() {
        StringBuilder summary = new StringBuilder();
        summary.append("Created Folder Structure:\n");

        if (rootFolder != null) {
            summary.append("- ").append(rootFolder.getName()).append(" (root)\n");
        }

        // Sort paths for better readability
        folderCache.entrySet().stream()
                .filter(entry -> !entry.getKey().isEmpty()) // Skip root folder
                .sorted(Map.Entry.comparingByKey())
                .forEach(entry -> {
                    String path = entry.getKey();
                    AnnotatedPluginDocument folder = entry.getValue();
                    int depth = path.split(File.separator.equals("\\") ? "\\\\" : File.separator).length;
                    String indent = "  ".repeat(depth);
                    summary.append("-").append(indent).append(folder.getName()).append("\n");
                });

        return summary.toString();
    }

    /**
     * Validates that the folder structure is consistent.
     *
     * @return true if the structure is valid
     */
    public boolean validateStructure() {
        if (rootFolder == null) {
            logger.warning("No root folder created");
            return false;
        }

        // Check that all cached folders are valid
        for (Map.Entry<String, AnnotatedPluginDocument> entry : folderCache.entrySet()) {
            if (entry.getValue() == null) {
                logger.warning("Null folder document for path: " + entry.getKey());
                return false;
            }
        }

        logger.info("Folder structure validation passed: " + folderCache.size() + " folders");
        return true;
    }
}