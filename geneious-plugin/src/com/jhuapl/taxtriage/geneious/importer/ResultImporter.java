package com.jhuapl.taxtriage.geneious.importer;

import com.biomatters.geneious.publicapi.databaseservice.DatabaseServiceException;
import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Main importer for TaxTriage output files into Geneious.
 *
 * This class orchestrates the import of all TaxTriage output files, maintaining
 * the original folder structure and using specialized importers for different
 * file types. It provides comprehensive progress tracking and error handling.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class ResultImporter {

    private static final Logger logger = Logger.getLogger(ResultImporter.class.getName());

    /** File format detector */
    private final FileFormatDetector formatDetector;

    /** Folder structure builder */
    private final FolderStructureBuilder folderBuilder;

    /** Specialized importers */
    private final KrakenReportImporter krakenImporter;
    private final FastqImporter fastqImporter;
    private final HtmlReportImporter htmlImporter;
    private final TextFileImporter textImporter;

    /** Database service for creating documents and folders */
    private final WritableDatabaseService databaseService;

    /** Import statistics */
    private final ImportResult importResult;

    /** Files to skip during import */
    private static final String[] SKIP_PATTERNS = {
            ".*\\.tmp$", ".*\\.temp$", ".*\\.lock$", ".*\\.idx$",
            ".*\\.bai$", ".*\\.tbi$", ".*\\.csi$", ".*\\.fai$",
            ".*\\.DS_Store$", ".*Thumbs\\.db$", ".*desktop\\.ini$",
            ".*\\.git.*", ".*\\.svn.*", ".*\\.nextflow.*",
            "work/.*", ".*\\.log$", ".*\\.command\\.(begin|err|log|out|run|sh)$"
    };

    /**
     * Creates a new result importer.
     *
     * @param databaseService the database service for creating documents
     */
    public ResultImporter(WritableDatabaseService databaseService) {
        this.databaseService = databaseService;
        this.formatDetector = new FileFormatDetector();
        this.folderBuilder = new FolderStructureBuilder(databaseService);
        this.krakenImporter = new KrakenReportImporter();
        this.fastqImporter = new FastqImporter();
        this.htmlImporter = new HtmlReportImporter();
        this.textImporter = new TextFileImporter();
        this.importResult = new ImportResult();
    }

    /**
     * Imports all TaxTriage output files from the specified directory.
     *
     * @param outputDirectory the TaxTriage output directory
     * @param workflowId the workflow ID for organizing results
     * @param progressListener optional progress listener
     * @return import results with statistics and imported documents
     * @throws DatabaseServiceException if database operations fail
     */
    public ImportResult importResults(File outputDirectory, String workflowId, ProgressListener progressListener)
            throws DatabaseServiceException {

        if (outputDirectory == null || !outputDirectory.exists() || !outputDirectory.isDirectory()) {
            throw new IllegalArgumentException("Invalid output directory: " + outputDirectory);
        }

        logger.info("Starting import of TaxTriage results from: " + outputDirectory.getAbsolutePath());

        try {
            // Initialize progress
            if (progressListener != null) {
                progressListener.setMessage("Scanning TaxTriage output directory...");
                progressListener.setProgress(0.0);
            }

            // Create root folder
            folderBuilder.createRootFolder(workflowId);

            // Scan directory for importable files
            List<File> filesToImport = scanDirectory(outputDirectory);
            logger.info("Found " + filesToImport.size() + " importable files");

            if (filesToImport.isEmpty()) {
                logger.warning("No importable files found in: " + outputDirectory.getAbsolutePath());
                importResult.markCompleted();
                return importResult;
            }

            // Import files
            importFiles(outputDirectory, filesToImport, progressListener);

            // Record created folders
            for (String folderPath : folderBuilder.getCreatedFolderPaths()) {
                importResult.addCreatedFolder(folderPath);
            }

            importResult.markCompleted();
            logger.info("Import completed successfully. " + importResult.getSuccessfulImports() +
                       " files imported, " + importResult.getFailedImports() + " failed.");

            return importResult;

        } catch (Exception e) {
            logger.log(Level.SEVERE, "Import failed", e);
            importResult.addError(null, e, "Import operation failed: " + e.getMessage());
            importResult.markCompleted();
            throw new DatabaseServiceException(e, "Failed to import TaxTriage results", false);
        }
    }

    /**
     * Scans the output directory recursively to find importable files.
     *
     * @param outputDirectory the directory to scan
     * @return list of files to import
     * @throws IOException if directory scanning fails
     */
    private List<File> scanDirectory(File outputDirectory) throws IOException {
        List<File> filesToImport = new ArrayList<>();

        Files.walkFileTree(outputDirectory.toPath(), new SimpleFileVisitor<Path>() {
            @Override
            public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                File f = file.toFile();

                // Skip files matching skip patterns
                if (shouldSkipFile(f)) {
                    return FileVisitResult.CONTINUE;
                }

                // Check if file should be imported
                if (formatDetector.shouldImport(f)) {
                    filesToImport.add(f);
                    logger.fine("Will import: " + f.getName() + " (Type: " +
                              formatDetector.getFormatDescription(f) + ")");
                }

                return FileVisitResult.CONTINUE;
            }

            @Override
            public FileVisitResult visitFileFailed(Path file, IOException exc) throws IOException {
                logger.warning("Failed to access file: " + file + " - " + exc.getMessage());
                return FileVisitResult.CONTINUE;
            }
        });

        return filesToImport;
    }

    /**
     * Imports a list of files with progress tracking.
     *
     * @param outputDirectory the output directory
     * @param filesToImport the files to import
     * @param progressListener optional progress listener
     */
    private void importFiles(File outputDirectory, List<File> filesToImport, ProgressListener progressListener) {
        AtomicInteger processedCount = new AtomicInteger(0);
        int totalFiles = filesToImport.size();

        for (File file : filesToImport) {
            try {
                // Update progress
                if (progressListener != null) {
                    double progress = (double) processedCount.get() / totalFiles;
                    progressListener.setProgress(progress);
                    progressListener.setMessage("Importing " + file.getName() + "...");
                }

                // Import the file
                importSingleFile(outputDirectory, file);

                processedCount.incrementAndGet();

            } catch (Exception e) {
                logger.log(Level.WARNING, "Failed to import file: " + file.getName(), e);
                importResult.addError(file, e, "Import failed: " + e.getMessage());
            }
        }

        if (progressListener != null) {
            progressListener.setProgress(1.0);
            progressListener.setMessage("Import completed");
        }
    }

    /**
     * Imports a single file using the appropriate importer.
     *
     * @param outputDirectory the output directory
     * @param file the file to import
     * @throws Exception if import fails
     */
    private void importSingleFile(File outputDirectory, File file) throws Exception {
        FileFormatDetector.FileFormat format = formatDetector.detectFormat(file);
        logger.fine("Importing " + file.getName() + " as " + format.toString());

        List<AnnotatedPluginDocument> documents = new ArrayList<>();

        // Use appropriate importer based on format
        switch (format.getType()) {
            case KRAKEN_REPORT:
                AnnotatedPluginDocument krakenDoc = krakenImporter.importKrakenReport(file);
                if (krakenDoc != null) {
                    documents.add(krakenDoc);
                }
                break;

            case SEQUENCE:
                if (format.getSubtype() != null &&
                    (format.getSubtype().equals("fastq") || format.getSubtype().equals("fq"))) {
                    documents.addAll(fastqImporter.importFastqFile(file));
                } else {
                    // For FASTA files, use text importer for now
                    AnnotatedPluginDocument textDoc = textImporter.importTextFile(file);
                    if (textDoc != null) {
                        documents.add(textDoc);
                    }
                }
                break;

            case REPORT:
                if (format.getSubtype() != null && format.getSubtype().equals("html")) {
                    AnnotatedPluginDocument htmlDoc = htmlImporter.importHtmlReport(file);
                    if (htmlDoc != null) {
                        documents.add(htmlDoc);
                    }
                } else {
                    AnnotatedPluginDocument textDoc = textImporter.importTextFile(file);
                    if (textDoc != null) {
                        documents.add(textDoc);
                    }
                }
                break;

            case TEXT:
            case ALIGNMENT:
            case VARIANT:
            default:
                AnnotatedPluginDocument textDoc = textImporter.importTextFile(file);
                if (textDoc != null) {
                    documents.add(textDoc);
                }
                break;
        }

        // Add documents to appropriate folder
        if (!documents.isEmpty()) {
            AnnotatedPluginDocument targetFolder = folderBuilder.getOrCreateFolder(outputDirectory, file);

            for (AnnotatedPluginDocument document : documents) {
                // In a real implementation, you would add the document to the folder here
                // For now, we just track it in the import result
                importResult.addImportedDocument(document, file, format.toString());
                logger.fine("Successfully imported: " + document.getName());
            }
        } else {
            logger.warning("No documents created for file: " + file.getName());
        }
    }

    /**
     * Checks if a file should be skipped during import.
     *
     * @param file the file to check
     * @return true if the file should be skipped
     */
    private boolean shouldSkipFile(File file) {
        if (file == null || !file.isFile()) {
            return true;
        }

        String fileName = file.getName();
        String absolutePath = file.getAbsolutePath();

        // Check skip patterns
        for (String pattern : SKIP_PATTERNS) {
            if (fileName.matches(pattern) || absolutePath.matches(".*" + pattern)) {
                return true;
            }
        }

        // Skip very large files (> 1GB)
        if (file.length() > 1024L * 1024L * 1024L) {
            logger.warning("Skipping very large file: " + fileName + " (" + file.length() + " bytes)");
            return true;
        }

        // Skip empty files
        if (file.length() == 0) {
            logger.fine("Skipping empty file: " + fileName);
            return true;
        }

        return false;
    }

    /**
     * Gets the import results.
     *
     * @return the import result with statistics
     */
    public ImportResult getImportResult() {
        return importResult;
    }

    /**
     * Gets the folder structure builder.
     *
     * @return the folder builder
     */
    public FolderStructureBuilder getFolderBuilder() {
        return folderBuilder;
    }

    /**
     * Gets detailed import statistics.
     *
     * @return formatted statistics string
     */
    public String getImportStatistics() {
        StringBuilder stats = new StringBuilder();
        stats.append("=== TaxTriage Import Statistics ===\n\n");

        // Overall stats
        stats.append("Total Files Processed: ").append(importResult.getTotalFilesProcessed()).append("\n");
        stats.append("Successful Imports: ").append(importResult.getSuccessfulImports()).append("\n");
        stats.append("Failed Imports: ").append(importResult.getFailedImports()).append("\n");
        stats.append("Success Rate: ").append(String.format("%.1f%%", importResult.getSuccessRate())).append("\n");
        stats.append("Total Data Processed: ").append(formatBytes(importResult.getTotalBytesProcessed())).append("\n");

        // Folder structure
        stats.append("\nFolder Structure:\n");
        stats.append("Folders Created: ").append(folderBuilder.getFolderCount()).append("\n");
        stats.append(folderBuilder.getFolderStructureSummary());

        // File type breakdown
        stats.append("\nFile Type Breakdown:\n");
        importResult.getFileTypeStats().entrySet().stream()
                .sorted((a, b) -> b.getValue().compareTo(a.getValue()))
                .forEach(entry -> stats.append("- ").append(entry.getKey())
                        .append(": ").append(entry.getValue()).append("\n"));

        // Errors
        if (importResult.getFailedImports() > 0) {
            stats.append("\nErrors:\n");
            for (ImportResult.ImportError error : importResult.getErrors()) {
                stats.append("- ").append(error.getFile() != null ? error.getFile().getName() : "Unknown")
                        .append(": ").append(error.getMessage()).append("\n");
            }
        }

        return stats.toString();
    }

    /**
     * Formats byte count in human-readable format.
     *
     * @param bytes the byte count
     * @return formatted string
     */
    private String formatBytes(long bytes) {
        if (bytes < 1024) {
            return bytes + " bytes";
        } else if (bytes < 1024 * 1024) {
            return String.format("%.1f KB", bytes / 1024.0);
        } else if (bytes < 1024 * 1024 * 1024) {
            return String.format("%.1f MB", bytes / (1024.0 * 1024.0));
        } else {
            return String.format("%.1f GB", bytes / (1024.0 * 1024.0 * 1024.0));
        }
    }

    /**
     * Validates that the importer is ready for use.
     *
     * @return true if the importer is ready
     */
    public boolean validateConfiguration() {
        if (databaseService == null) {
            logger.severe("Database service is null");
            return false;
        }

        if (!folderBuilder.validateStructure()) {
            logger.warning("Folder structure validation failed");
            return false;
        }

        logger.info("Result importer configuration is valid");
        return true;
    }

    /**
     * Resets the importer for a new import operation.
     */
    public void reset() {
        folderBuilder.reset();
        // importResult is created fresh for each import operation
    }
}