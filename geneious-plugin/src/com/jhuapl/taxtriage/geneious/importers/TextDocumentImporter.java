package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentImportException;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Imports text and CSV files as text documents, preventing CSV parser from trying to parse as sequences.
 */
public class TextDocumentImporter {

    private static final Logger logger = Logger.getLogger(TextDocumentImporter.class.getName());

    /**
     * Imports a text, CSV, or TSV file as a text document.
     * This prevents the CSV/TSV parser from trying to interpret it as sequences.
     *
     * @param file The text, CSV, or TSV file to import
     * @param progressListener Progress listener
     * @return List containing the imported text document
     */
    public static List<AnnotatedPluginDocument> importAsTextDocument(File file, ProgressListener progressListener)
            throws IOException, DocumentImportException {

        if (!file.exists()) {
            throw new DocumentImportException("File does not exist: " + file.getAbsolutePath());
        }

        // For CSV/TSV files, rename to .txt to prevent CSV/TSV parser
        File importFile = file;
        boolean isTempFile = false;
        String lowerName = file.getName().toLowerCase();

        if (lowerName.endsWith(".csv") || lowerName.endsWith(".tsv") ||
            lowerName.endsWith(".krona") || lowerName.endsWith(".kraken")) {
            // Create a temp .txt file with the content
            String content = new String(Files.readAllBytes(file.toPath()));
            importFile = File.createTempFile("taxtriage_", ".txt");
            importFile.deleteOnExit();
            Files.write(importFile.toPath(), content.getBytes());
            isTempFile = true;
            logger.info("Converting " + file.getName() + " to text document");
        }

        // Import the file
        List<AnnotatedPluginDocument> documents = PluginUtilities.importDocuments(importFile, progressListener);

        // Set proper name (without extension)
        String name = file.getName();
        if (name.contains(".")) {
            name = name.substring(0, name.lastIndexOf('.'));
        }

        for (AnnotatedPluginDocument doc : documents) {
            doc.setName(name);
        }

        // Clean up temp file
        if (isTempFile) {
            importFile.delete();
        }

        logger.info("Imported as text document: " + file.getName() + " (" + documents.size() + " document(s))");
        return documents;
    }

    /**
     * Determines if a file should be imported as a text document.
     * @param file The file to check
     * @return true if the file should be imported as text
     */
    public static boolean shouldImportAsText(File file) {
        String name = file.getName().toLowerCase();

        // Import these file types as text documents
        return name.endsWith(".txt") ||
               name.endsWith(".csv") ||
               name.endsWith(".tsv") ||
               name.endsWith(".report") ||
               name.endsWith(".krona") ||
               name.endsWith(".kraken") ||
               name.endsWith(".log");
    }
}