package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentImportException;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Direct text importer that creates text documents without relying on file extension detection.
 * Uses temp file approach to guarantee text import.
 */
public class DirectTextImporter {

    private static final Logger logger = Logger.getLogger(DirectTextImporter.class.getName());

    /**
     * Directly imports a text document by converting to .txt first.
     */
    public static List<AnnotatedPluginDocument> directImportText(File file, WritableDatabaseService targetFolder,
                                                                 ProgressListener progressListener)
            throws IOException, DocumentImportException {

        if (!file.exists()) {
            throw new DocumentImportException("File does not exist: " + file.getAbsolutePath());
        }

        logger.info("Direct text import: " + file.getName());

        // Read file content
        byte[] content = Files.readAllBytes(file.toPath());

        // Get the base name without extension
        String originalName = file.getName();
        String baseName = originalName;
        if (baseName.contains(".")) {
            baseName = baseName.substring(0, baseName.lastIndexOf('.'));
        }

        // Create temp file with .txt extension to force text import
        File tempFile = File.createTempFile("taxtriage_text_", ".txt");
        tempFile.deleteOnExit();

        try {
            // Write content to temp file
            Files.write(tempFile.toPath(), content);

            // Import the temp file
            List<AnnotatedPluginDocument> docs;
            if (targetFolder != null) {
                // Import directly to the target folder
                docs = PluginUtilities.importDocumentsToDatabase(tempFile, targetFolder, progressListener);
            } else {
                docs = PluginUtilities.importDocuments(tempFile, progressListener);
            }

            // Fix the name to match original
            if (docs != null && !docs.isEmpty()) {
                for (AnnotatedPluginDocument doc : docs) {
                    try {
                        doc.setName(baseName);
                    } catch (Exception e) {
                        logger.warning("Could not rename document: " + e.getMessage());
                    }
                }
                logger.info("Successfully imported text file: " + originalName + " (" + docs.size() + " documents)");
            } else {
                logger.warning("No documents imported from: " + originalName);
                docs = new ArrayList<>();
            }

            return docs;

        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to import text file: " + file.getName(), e);
            throw new DocumentImportException("Failed to import " + file.getName() + ": " + e.getMessage());

        } finally {
            // Clean up temp file
            try {
                tempFile.delete();
            } catch (Exception e) {
                // Ignore cleanup errors
            }
        }
    }

    /**
     * Batch import text files directly to a folder.
     */
    public static List<AnnotatedPluginDocument> batchImportTextFiles(List<File> files,
                                                                     WritableDatabaseService targetFolder,
                                                                     ProgressListener progressListener) {
        List<AnnotatedPluginDocument> allDocs = new ArrayList<>();

        for (int i = 0; i < files.size(); i++) {
            File file = files.get(i);

            if (progressListener != null && progressListener.isCanceled()) {
                break;
            }

            if (progressListener != null) {
                progressListener.setProgress((double) i / files.size());
                progressListener.setMessage("Importing text: " + file.getName());
            }

            try {
                List<AnnotatedPluginDocument> docs = directImportText(file, targetFolder, progressListener);
                allDocs.addAll(docs);

            } catch (Exception e) {
                logger.log(Level.WARNING, "Failed to import text file: " + file.getName(), e);
                System.out.println("      Error importing " + file.getName() + ": " + e.getMessage());
            }
        }

        if (progressListener != null) {
            progressListener.setProgress(1.0);
        }

        return allDocs;
    }

    /**
     * Test if a file should be imported as text.
     */
    public static boolean isTextFile(File file) {
        if (!file.isFile()) {
            return false;
        }

        String name = file.getName().toLowerCase();

        // Comprehensive list of text file patterns
        return name.endsWith(".txt") ||
               name.endsWith(".tsv") ||
               name.endsWith(".csv") ||
               name.endsWith(".log") ||
               name.endsWith(".report") ||
               name.contains(".kraken") ||
               name.contains(".krona") ||
               name.contains("kraken") ||
               name.contains("krona") ||
               name.startsWith("taxtriage_") ||
               name.contains("taxtriage") ||
               name.endsWith(".histo") ||
               name.endsWith(".paths") ||
               name.endsWith(".topnames") ||
               name.endsWith(".toptaxids") ||
               name.endsWith(".gcfids") ||
               name.endsWith(".combined") ||
               name.contains("report") ||
               name.contains("count");
    }

    /**
     * Alternative import method using multiple attempts.
     */
    public static List<AnnotatedPluginDocument> robustTextImport(File file,
                                                                 WritableDatabaseService targetFolder,
                                                                 ProgressListener progressListener) {
        List<AnnotatedPluginDocument> docs = null;

        // Try 1: Direct import with .txt extension
        try {
            docs = directImportText(file, targetFolder, progressListener);
            if (docs != null && !docs.isEmpty()) {
                return docs;
            }
        } catch (Exception e) {
            logger.info("First import attempt failed, trying alternative: " + e.getMessage());
        }

        // Try 2: Import to null folder then copy
        try {
            docs = directImportText(file, null, progressListener);
            if (docs != null && !docs.isEmpty() && targetFolder != null) {
                for (AnnotatedPluginDocument doc : docs) {
                    targetFolder.addDocumentCopy(doc, progressListener);
                }
                return docs;
            }
        } catch (Exception e) {
            logger.warning("Alternative import also failed: " + e.getMessage());
        }

        return new ArrayList<>();
    }
}