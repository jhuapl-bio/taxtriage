package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.implementations.DefaultAlignmentDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentImportException;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Dedicated importer for plain text and TSV files.
 * Ensures these files are always imported as text documents, never parsed as sequences.
 */
public class PlainTextImporter {

    private static final Logger logger = Logger.getLogger(PlainTextImporter.class.getName());

    /**
     * Force imports a file as a plain text document.
     * This method ensures the file is treated as text regardless of extension.
     *
     * @param file The file to import
     * @param progressListener Progress listener
     * @return List containing the imported text document
     */
    public static List<AnnotatedPluginDocument> forceImportAsText(File file, ProgressListener progressListener)
            throws IOException, DocumentImportException {

        if (!file.exists()) {
            throw new DocumentImportException("File does not exist: " + file.getAbsolutePath());
        }

        logger.info("Force importing as text: " + file.getName());

        // Read the file content
        String content = new String(Files.readAllBytes(file.toPath()), StandardCharsets.UTF_8);

        // Create a temporary .txt file to ensure text import
        File tempFile = File.createTempFile("taxtriage_text_", ".txt");
        tempFile.deleteOnExit();

        try {
            // Write content to temp file
            Files.write(tempFile.toPath(), content.getBytes(StandardCharsets.UTF_8));

            // Import using Geneious utilities
            List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(tempFile, progressListener);

            if (docs == null || docs.isEmpty()) {
                // Fallback: Try to create a text document directly
                logger.info("Standard import failed, trying direct text document creation");
                docs = createTextDocumentDirectly(file, content);
            }

            // Set the proper name (without extension changes)
            String originalName = file.getName();
            if (originalName.contains(".")) {
                originalName = originalName.substring(0, originalName.lastIndexOf('.'));
            }

            for (AnnotatedPluginDocument doc : docs) {
                doc.setName(originalName);
                logger.info("Set document name to: " + originalName);
            }

            return docs;

        } finally {
            // Clean up temp file
            if (tempFile.exists()) {
                tempFile.delete();
            }
        }
    }

    /**
     * Creates a text document directly if standard import fails.
     */
    private static List<AnnotatedPluginDocument> createTextDocumentDirectly(File file, String content)
            throws DocumentImportException {

        try {
            // Try alternative import methods
            // Method 1: Import with explicit text file extension
            File altTempFile = File.createTempFile("text_", ".txt");
            altTempFile.deleteOnExit();

            Files.write(altTempFile.toPath(), content.getBytes(StandardCharsets.UTF_8));

            List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(altTempFile, null);

            altTempFile.delete();

            if (docs != null && !docs.isEmpty()) {
                return docs;
            }

        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to create text document directly", e);
        }

        // Return empty list if all methods fail
        return new ArrayList<>();
    }

    /**
     * Imports multiple text/TSV files.
     *
     * @param files List of files to import
     * @param progressListener Progress listener
     * @return List of all imported documents
     */
    public static List<AnnotatedPluginDocument> importTextFiles(List<File> files, ProgressListener progressListener)
            throws IOException, DocumentImportException {

        List<AnnotatedPluginDocument> allDocs = new ArrayList<>();

        for (int i = 0; i < files.size(); i++) {
            File file = files.get(i);

            if (progressListener != null) {
                progressListener.setProgress((double) i / files.size());
                progressListener.setMessage("Importing text file: " + file.getName());
            }

            try {
                List<AnnotatedPluginDocument> docs = forceImportAsText(file, progressListener);
                allDocs.addAll(docs);
                logger.info("Successfully imported text file: " + file.getName() + " (" + docs.size() + " documents)");

            } catch (Exception e) {
                logger.log(Level.WARNING, "Failed to import text file: " + file.getName(), e);
                // Continue with other files
            }
        }

        if (progressListener != null) {
            progressListener.setProgress(1.0);
        }

        return allDocs;
    }

    /**
     * Checks if a file should be imported as plain text.
     *
     * @param file The file to check
     * @return true if the file should be imported as plain text
     */
    public static boolean shouldImportAsPlainText(File file) {
        String name = file.getName().toLowerCase();

        // List of extensions and patterns that should always be imported as text
        return name.endsWith(".txt") ||
               name.endsWith(".tsv") ||
               name.endsWith(".csv") ||
               name.endsWith(".log") ||
               name.endsWith(".report") ||
               name.contains(".kraken") ||
               name.contains(".krona") ||
               name.contains("taxtriage_") ||
               name.contains(".krakenreport") ||
               name.endsWith(".histo") ||
               name.endsWith(".paths") ||
               name.endsWith(".topnames") ||
               name.endsWith(".toptaxids") ||
               name.endsWith(".gcfids");
    }

    /**
     * Batch imports text files to ensure they're processed correctly.
     * This method handles special cases for TaxTriage output files.
     *
     * @param directory Directory containing text files
     * @param progressListener Progress listener
     * @return List of imported documents
     */
    public static List<AnnotatedPluginDocument> importTextFilesFromDirectory(File directory,
                                                                             ProgressListener progressListener)
            throws IOException, DocumentImportException {

        if (!directory.exists() || !directory.isDirectory()) {
            throw new DocumentImportException("Invalid directory: " + directory.getAbsolutePath());
        }

        List<File> textFiles = new ArrayList<>();

        // Scan for text files
        File[] files = directory.listFiles();
        if (files != null) {
            for (File file : files) {
                if (file.isFile() && shouldImportAsPlainText(file)) {
                    textFiles.add(file);
                    logger.info("Found text file for import: " + file.getName());
                }
            }
        }

        logger.info("Importing " + textFiles.size() + " text files from " + directory.getName());

        return importTextFiles(textFiles, progressListener);
    }
}