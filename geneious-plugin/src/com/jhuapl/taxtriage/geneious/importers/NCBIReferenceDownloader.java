package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.databaseservice.DatabaseServiceException;
import com.biomatters.geneious.publicapi.databaseservice.Query;
import com.biomatters.geneious.publicapi.databaseservice.RetrieveCallback;
import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentImportException;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import jebl.util.ProgressListener;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Downloads reference sequences from NCBI in GenBank format.
 * This provides full annotation information for the reference sequences.
 */
public class NCBIReferenceDownloader {

    private static final Logger logger = Logger.getLogger(NCBIReferenceDownloader.class.getName());

    // NCBI E-utilities base URLs
    private static final String EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
    private static final String ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";

    /**
     * Downloads GenBank files for the given accessions from NCBI.
     *
     * @param accessions List of NCBI accession numbers
     * @param outputDir Directory to save the GenBank files
     * @param progressListener Progress listener
     * @return Map of accession to downloaded file path
     */
    public static Map<String, File> downloadGenBankFiles(List<String> accessions,
                                                         Path outputDir,
                                                         ProgressListener progressListener) {
        Map<String, File> downloadedFiles = new HashMap<>();

        System.out.println("\n=== NCBIReferenceDownloader: Starting GenBank downloads ===");
        System.out.println("  Accessions to download: " + accessions.size());
        for (String acc : accessions) {
            System.out.println("    - " + acc);
        }
        System.out.println("  Output directory: " + outputDir);

        if (accessions == null || accessions.isEmpty()) {
            logger.warning("No accessions provided for download");
            System.out.println("  ERROR: No accessions provided for download");
            return downloadedFiles;
        }

        try {
            Files.createDirectories(outputDir);
            System.out.println("  Output directory created/verified: " + outputDir);
        } catch (IOException e) {
            logger.log(Level.SEVERE, "Could not create output directory", e);
            System.out.println("  ERROR: Could not create output directory: " + e.getMessage());
            return downloadedFiles;
        }

        int count = 0;
        for (String accession : accessions) {
            if (progressListener != null && progressListener.isCanceled()) {
                break;
            }

            if (progressListener != null) {
                progressListener.setProgress((double) count / accessions.size());
                progressListener.setMessage("Downloading " + accession + " from NCBI...");
            }

            System.out.println("\n  Processing accession " + (count + 1) + "/" + accessions.size() + ": " + accession);

            // Check if file already exists
            File existingFile = outputDir.resolve(accession + ".gb").toFile();
            File gbFile = null;

            if (existingFile.exists() && existingFile.length() > 100) {
                System.out.println("    ✓ GenBank file already exists: " + existingFile.getName() + " (" + existingFile.length() + " bytes)");
                System.out.println("    Skipping download from NCBI");
                logger.info("Using existing GenBank file for " + accession);
                gbFile = existingFile;
            } else {
                System.out.println("    Downloading from NCBI...");
                try {
                    gbFile = downloadSingleGenBankFile(accession, outputDir);
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Failed to download GenBank file for " + accession, e);
                    System.out.println("    ✗ Download failed: " + e.getMessage());
                    e.printStackTrace(System.out);
                    count++;
                    continue;
                }
            }

            try {
                if (gbFile != null && gbFile.exists()) {
                    // Fix the GenBank file to ensure LOCUS and ACCESSION match VERSION
                    logger.info("Fixing GenBank file fields for " + accession);
                    System.out.println("    Fixing GenBank file fields...");
                    boolean fixed = GenBankFileFixer.fixGenBankFile(gbFile);
                    if (fixed) {
                        logger.info("Successfully fixed GenBank file for " + accession);
                        System.out.println("    ✓ GenBank file fixed successfully");
                    } else {
                        logger.warning("Could not fix GenBank file for " + accession + " - may cause import issues");
                        System.out.println("    ⚠ Warning: Could not fix GenBank file - may cause import issues");
                    }

                    downloadedFiles.put(accession, gbFile);
                    logger.info("GenBank file ready for " + accession);
                } else {
                    System.out.println("    ✗ File not available");
                }
            } catch (Exception e) {
                logger.log(Level.WARNING, "Failed to process GenBank file for " + accession, e);
                System.out.println("    ✗ Processing failed: " + e.getMessage());
                e.printStackTrace(System.out);
            }

            count++;

            // Add a small delay between requests to be respectful to NCBI servers
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                break;
            }
        }

        if (progressListener != null) {
            progressListener.setProgress(1.0);
        }

        logger.info("Downloaded " + downloadedFiles.size() + " of " + accessions.size() + " GenBank files");
        System.out.println("\n=== Download Summary ===");
        System.out.println("  Successfully downloaded: " + downloadedFiles.size() + " of " + accessions.size() + " GenBank files");
        if (downloadedFiles.size() < accessions.size()) {
            System.out.println("  ⚠ Warning: Some downloads failed. Check errors above.");
        }
        System.out.println("=== NCBIReferenceDownloader: Complete ===\n");
        return downloadedFiles;
    }

    /**
     * Downloads a single GenBank file from NCBI.
     */
    private static File downloadSingleGenBankFile(String accession, Path outputDir) throws IOException {
        // Build the E-fetch URL for GenBank format
        String urlString = String.format("%s?db=nuccore&id=%s&rettype=gbwithparts&retmode=text",
                                        EFETCH_URL, accession);

        System.out.println("    Connecting to NCBI: " + urlString);
        logger.fine("Downloading from: " + urlString);

        URL url = new URL(urlString);
        HttpURLConnection conn = (HttpURLConnection) url.openConnection();
        conn.setRequestMethod("GET");
        conn.setConnectTimeout(30000);
        conn.setReadTimeout(30000);

        System.out.println("    Sending HTTP request...");
        int responseCode = conn.getResponseCode();
        System.out.println("    HTTP response code: " + responseCode);

        if (responseCode != HttpURLConnection.HTTP_OK) {
            String errorMsg = "HTTP error " + responseCode + " downloading " + accession;
            System.out.println("    ✗ " + errorMsg);
            throw new IOException(errorMsg);
        }

        // Save to file
        File outputFile = outputDir.resolve(accession + ".gb").toFile();
        System.out.println("    Saving to: " + outputFile.getAbsolutePath());

        try (InputStream is = conn.getInputStream();
             BufferedInputStream bis = new BufferedInputStream(is);
             FileOutputStream fos = new FileOutputStream(outputFile);
             BufferedOutputStream bos = new BufferedOutputStream(fos)) {

            byte[] buffer = new byte[8192];
            int bytesRead;
            long totalBytes = 0;

            while ((bytesRead = bis.read(buffer)) != -1) {
                bos.write(buffer, 0, bytesRead);
                totalBytes += bytesRead;
            }

            System.out.println("    Downloaded " + totalBytes + " bytes");
            logger.fine("Downloaded " + totalBytes + " bytes for " + accession);
        }

        conn.disconnect();

        // Verify the file contains GenBank data
        if (outputFile.length() < 100) {
            logger.warning("Downloaded file for " + accession + " seems too small");
            System.out.println("    ✗ Warning: Downloaded file seems too small (" + outputFile.length() + " bytes)");
            outputFile.delete();
            return null;
        }

        return outputFile;
    }

    /**
     * Imports GenBank files using sample-based importer for organized storage.
     *
     * @param genBankFiles Map of accession to GenBank file
     * @param referenceNameMap Map of accession to desired reference name from BAM
     * @param progressListener Progress listener
     * @param sampleImporter Sample-based importer for organized storage
     * @return List of imported documents with corrected names
     */
    public static List<AnnotatedPluginDocument> importGenBankFiles(Map<String, File> genBankFiles,
                                                                   Map<String, String> referenceNameMap,
                                                                   ProgressListener progressListener,
                                                                   SampleBasedImporter sampleImporter) {
        // For now, just use the legacy method
        // The sample importer will handle the actual file organization
        return importGenBankFiles(genBankFiles, referenceNameMap, progressListener);
    }

    // REMOVED: Method using deleted FolderBasedImporter class
    // Use the SampleBasedImporter version instead

    /**
     * Imports GenBank files into Geneious and renames them to match BAM reference names.
     *
     * @param genBankFiles Map of accession to GenBank file
     * @param referenceNameMap Map of accession to desired reference name from BAM
     * @param progressListener Progress listener
     * @return List of imported documents with corrected names
     */
    public static List<AnnotatedPluginDocument> importGenBankFiles(Map<String, File> genBankFiles,
                                                                   Map<String, String> referenceNameMap,
                                                                   ProgressListener progressListener) {
        List<AnnotatedPluginDocument> importedDocs = new ArrayList<>();

        logger.info("===== IMPORTING GENBANK FILES =====");
        logger.info("Files to import: " + genBankFiles.size());
        logger.info("Reference name map has " + referenceNameMap.size() + " entries:");
        for (Map.Entry<String, String> mapEntry : referenceNameMap.entrySet()) {
            logger.info("  Map: '" + mapEntry.getKey() + "' -> '" + mapEntry.getValue() + "'");
        }

        for (Map.Entry<String, File> entry : genBankFiles.entrySet()) {
            String accession = entry.getKey();
            File gbFile = entry.getValue();

            logger.info("Processing GenBank file for accession: " + accession);

            if (!gbFile.exists()) {
                logger.warning("GenBank file does not exist: " + gbFile);
                continue;
            }

            try {
                // Import the GenBank file
                logger.info("Importing file: " + gbFile.getName() + " (" + gbFile.length() + " bytes)");
                List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(gbFile, progressListener);
                logger.info("Imported " + docs.size() + " document(s) from " + gbFile.getName());

                // Rename documents to match BAM reference names if needed
                String targetName = referenceNameMap != null ? referenceNameMap.get(accession) : null;
                logger.info("Target name for accession '" + accession + "': " +
                           (targetName != null ? "'" + targetName + "'" : "null (no mapping)"));

                if (targetName != null && !targetName.trim().isEmpty()) {
                    for (AnnotatedPluginDocument doc : docs) {
                        if (doc == null) {
                            logger.warning("Encountered null document in import list");
                            continue;
                        }
                        String originalName = doc.getName();
                        if (originalName == null) {
                            logger.warning("Document has null name, skipping rename");
                            continue;
                        }
                        try {
                            // Set the name to match what's in the BAM file
                            doc.setName(targetName);
                            logger.info("RENAMED: '" + originalName + "' -> '" + targetName + "'");
                        } catch (Exception e) {
                            logger.log(Level.WARNING, "Could not rename document from '" + originalName +
                                      "' to '" + targetName + "'", e);
                        }
                    }
                } else {
                    for (AnnotatedPluginDocument doc : docs) {
                        logger.info("NOT RENAMED (no mapping): '" + doc.getName() + "'");
                    }
                }

                importedDocs.addAll(docs);

            } catch (IOException | DocumentImportException e) {
                logger.log(Level.WARNING, "Failed to import GenBank file: " + gbFile, e);
            }
        }

        logger.info("===== GENBANK IMPORT COMPLETE =====");
        logger.info("Total imported documents: " + importedDocs.size());

        return importedDocs;
    }

    /**
     * Uses Geneious's built-in NCBI search capability to retrieve sequences.
     * This is an alternative approach that leverages Geneious's existing NCBI integration.
     */
    public static List<AnnotatedPluginDocument> searchAndRetrieveFromNCBI(List<String> accessions,
                                                                          ProgressListener progressListener) {
        List<AnnotatedPluginDocument> retrievedDocs = new ArrayList<>();

        try {
            // Try to get the NCBI database service from Geneious
            // This would typically be available if NCBI search is enabled in Geneious
            logger.info("Attempting to use Geneious NCBI search for " + accessions.size() + " accessions");

            // Note: The actual implementation would depend on the specific Geneious API
            // for accessing the NCBI database service, which may vary by version

            for (String accession : accessions) {
                if (progressListener != null && progressListener.isCanceled()) {
                    break;
                }

                if (progressListener != null) {
                    progressListener.setMessage("Searching NCBI for " + accession + "...");
                }

                // Create a query for the specific accession
                Query query = Query.Factory.createQuery(accession);

                // Note: This is a placeholder - the actual implementation would use
                // the appropriate Geneious API to search and retrieve from NCBI
                logger.info("Would search NCBI for: " + accession);
            }

        } catch (Exception e) {
            logger.log(Level.WARNING, "Could not use Geneious NCBI search", e);
        }

        return retrievedDocs;
    }

    /**
     * Creates a mapping between accessions and the full reference names from BAM.
     */
    public static Map<String, String> createReferenceNameMap(List<BamReferenceExtractor.ReferenceInfo> references) {
        Map<String, String> map = new HashMap<>();

        for (BamReferenceExtractor.ReferenceInfo ref : references) {
            if (ref.accession != null) {
                map.put(ref.accession, ref.name);
            }
        }

        return map;
    }

    /**
     * Extracts unique accessions from reference information.
     */
    public static List<String> extractAccessions(Collection<BamReferenceExtractor.ReferenceInfo> references) {
        Set<String> accessions = new LinkedHashSet<>();

        for (BamReferenceExtractor.ReferenceInfo ref : references) {
            if (ref.accession != null && !ref.accession.isEmpty()) {
                accessions.add(ref.accession);
            }
        }

        return new ArrayList<>(accessions);
    }
}