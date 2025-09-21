package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentImportException;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Imports BAM files and their references together in a single operation.
 * This approach mimics how Geneious might handle minimap2 results internally.
 */
public class CombinedImporter {

    private static final Logger logger = Logger.getLogger(CombinedImporter.class.getName());

    /**
     * Imports BAM file and its references by creating a temporary directory
     * with both files together, then importing the directory.
     */
    public static List<AnnotatedPluginDocument> importBamWithReferencesAtomic(
            File bamFile,
            List<File> referenceFiles,
            ProgressListener progressListener) {

        List<AnnotatedPluginDocument> importedDocs = new ArrayList<>();

        System.out.println("===== ATOMIC BAM+REFERENCE IMPORT =====");
        System.out.println("BAM file: " + bamFile.getName());
        System.out.println("Reference files: " + referenceFiles.size());

        Path tempDir = null;
        try {
            // Create temporary directory
            tempDir = Files.createTempDirectory("geneious_bam_import_");
            System.out.println("Created temp directory: " + tempDir);

            // Copy reference files to temp directory FIRST
            for (File refFile : referenceFiles) {
                Path destPath = tempDir.resolve(refFile.getName());
                Files.copy(refFile.toPath(), destPath, StandardCopyOption.REPLACE_EXISTING);
                System.out.println("  Copied reference: " + refFile.getName());
            }

            // Copy BAM file and its index to temp directory
            Path bamDest = tempDir.resolve(bamFile.getName());
            Files.copy(bamFile.toPath(), bamDest, StandardCopyOption.REPLACE_EXISTING);
            System.out.println("  Copied BAM: " + bamFile.getName());

            // Copy BAM index if it exists
            File baiFile = new File(bamFile.getAbsolutePath() + ".bai");
            if (baiFile.exists()) {
                Path baiDest = tempDir.resolve(baiFile.getName());
                Files.copy(baiFile.toPath(), baiDest, StandardCopyOption.REPLACE_EXISTING);
                System.out.println("  Copied BAM index: " + baiFile.getName());
            }

            // Import all files from the directory together
            System.out.println("Importing directory contents atomically...");
            File[] filesToImport = tempDir.toFile().listFiles((dir, name) ->
                !name.endsWith(".bai"));  // Don't directly import .bai files

            if (filesToImport != null) {
                for (File file : filesToImport) {
                    try {
                        System.out.println("  Importing: " + file.getName());
                        List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(
                            file, progressListener);

                        if (docs != null && !docs.isEmpty()) {
                            importedDocs.addAll(docs);
                            System.out.println("    Imported " + docs.size() + " document(s)");
                        }
                    } catch (Exception e) {
                        System.out.println("    Failed to import: " + e.getMessage());
                    }
                }
            }

            System.out.println("Total imported: " + importedDocs.size() + " document(s)");

        } catch (IOException e) {
            System.out.println("Error during atomic import: " + e.getMessage());
            logger.log(Level.WARNING, "Failed atomic import", e);
        } finally {
            // Clean up temp directory
            if (tempDir != null) {
                try {
                    Files.walk(tempDir)
                        .sorted((a, b) -> -a.compareTo(b))
                        .forEach(path -> {
                            try {
                                Files.delete(path);
                            } catch (IOException e) {
                                // Ignore
                            }
                        });
                    System.out.println("Cleaned up temp directory");
                } catch (IOException e) {
                    logger.log(Level.WARNING, "Failed to clean up temp directory", e);
                }
            }
        }

        System.out.println("======================================\n");
        return importedDocs;
    }

    /**
     * Alternative approach: Import references first, wait, then import BAM.
     * But ensures they're imported in the same "batch" or session.
     */
    public static List<AnnotatedPluginDocument> importInSequence(
            File bamFile,
            List<File> referenceFiles,
            ProgressListener progressListener) {

        List<AnnotatedPluginDocument> allDocs = new ArrayList<>();

        System.out.println("===== SEQUENTIAL IMPORT (SAME SESSION) =====");

        try {
            // First, import all reference files
            System.out.println("Step 1: Importing " + referenceFiles.size() + " reference file(s)");
            List<AnnotatedPluginDocument> refDocs = new ArrayList<>();

            for (File refFile : referenceFiles) {
                System.out.println("  Importing reference: " + refFile.getName());
                List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(
                    refFile, progressListener);

                if (docs != null && !docs.isEmpty()) {
                    refDocs.addAll(docs);
                    System.out.println("    SUCCESS: Imported " + docs.size() + " document(s)");
                    for (AnnotatedPluginDocument doc : docs) {
                        System.out.println("      - " + doc.getName());
                    }
                }
            }

            allDocs.addAll(refDocs);
            System.out.println("Total references imported: " + refDocs.size());

            // Short wait to ensure references are indexed
            System.out.println("\nStep 2: Waiting 2 seconds for indexing...");
            Thread.sleep(2000);

            // Now import the BAM file
            System.out.println("\nStep 3: Importing BAM file: " + bamFile.getName());
            List<AnnotatedPluginDocument> bamDocs = PluginUtilities.importDocuments(
                bamFile, progressListener);

            if (bamDocs != null && !bamDocs.isEmpty()) {
                allDocs.addAll(bamDocs);
                System.out.println("  SUCCESS: Imported " + bamDocs.size() + " BAM document(s)");
            } else {
                System.out.println("  FAILED: BAM import returned no documents");
            }

        } catch (Exception e) {
            System.out.println("Error during sequential import: " + e.getMessage());
            e.printStackTrace();
            logger.log(Level.WARNING, "Failed sequential import", e);
        }

        System.out.println("\nTotal documents imported: " + allDocs.size());
        System.out.println("==========================================\n");

        return allDocs;
    }
}