package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.databaseservice.DatabaseServiceException;
import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.plugin.DocumentImportException;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import com.jhuapl.taxtriage.geneious.documents.TaxTriageResultDocument;
import com.jhuapl.taxtriage.geneious.utils.FileTypeUtil;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * Enterprise-grade importer for organizing TaxTriage analysis results into hierarchical folder structures.
 *
 * <p>This class provides intelligent organization of complex TaxTriage workflow outputs into
 * a logical, sample-based folder hierarchy within Geneious. It handles multiple data types
 * including taxonomic reports, sequence alignments, reference sequences, and quality control
 * metrics, ensuring that results are easily navigable and accessible for downstream analysis.</p>
 *
 * <h3>Folder Organization Strategy:</h3>
 * <ul>
 *   <li><strong>Main Folder:</strong> Named by run timestamp for unique identification</li>
 *   <li><strong>Sample Folders:</strong> Individual folders for each detected sample</li>
 *   <li><strong>Type-specific Subfolders:</strong> Organized by data type and processing status</li>
 * </ul>
 *
 * <h3>Supported Data Types:</h3>
 * <ul>
 *   <li><strong>Alignment Results:</strong> BAM files with separate folders for original and deduplicated reads</li>
 *   <li><strong>Reference Sequences:</strong> GenBank files downloaded from NCBI with field validation</li>
 *   <li><strong>Taxonomic Reports:</strong> Kraken2 reports, Krona visualizations, and summary statistics</li>
 *   <li><strong>Quality Control:</strong> MultiQC reports and processing metrics</li>
 * </ul>
 *
 * <h3>Smart Sample Detection:</h3>
 * <ul>
 *   <li><strong>File Pattern Recognition:</strong> Automatic sample name extraction from filenames</li>
 *   <li><strong>Cross-directory Scanning:</strong> Sample detection across multiple output directories</li>
 *   <li><strong>Deduplication Awareness:</strong> Intelligent handling of original vs. deduplicated files</li>
 * </ul>
 *
 * <h3>Reference Management:</h3>
 * <ul>
 *   <li><strong>GenBank Integration:</strong> Automatic download of reference sequences from NCBI</li>
 *   <li><strong>Field Validation:</strong> GenBank file validation and automatic field correction</li>
 *   <li><strong>Co-location Strategy:</strong> References stored alongside corresponding BAM files</li>
 * </ul>
 *
 * <h3>Import Process:</h3>
 * <ol>
 *   <li><strong>Discovery:</strong> Scan output directory to identify samples and data types</li>
 *   <li><strong>Structure Creation:</strong> Build hierarchical folder structure in Geneious</li>
 *   <li><strong>Reference Processing:</strong> Download, validate, and import reference sequences</li>
 *   <li><strong>Text File Import:</strong> Process all text-based reports and metrics</li>
 *   <li><strong>Alignment Import:</strong> Import BAM files with reference linking</li>
 *   <li><strong>Quality Assurance:</strong> Validation and error reporting</li>
 * </ol>
 *
 * <h3>Performance Optimizations:</h3>
 * <ul>
 *   <li><strong>Efficient Scanning:</strong> Single-pass directory traversal with pattern matching</li>
 *   <li><strong>Deduplication:</strong> Prevents duplicate imports of identical files</li>
 *   <li><strong>Batch Processing:</strong> Optimized bulk operations for large datasets</li>
 *   <li><strong>Progress Tracking:</strong> Real-time feedback for long-running operations</li>
 * </ul>
 *
 * <h3>Error Recovery:</h3>
 * <ul>
 *   <li><strong>Graceful Failures:</strong> Individual file failures don't abort entire import</li>
 *   <li><strong>GenBank Fallback:</strong> Multiple strategies for reference sequence acquisition</li>
 *   <li><strong>Validation Retry:</strong> Automatic retry with field correction for malformed files</li>
 * </ul>
 *
 * @author TaxTriage Development Team
 * @version 2.0
 * @since 1.0
 */
public class SampleBasedImporter {

    private static final Logger logger = Logger.getLogger(SampleBasedImporter.class.getName());

    /**
     * Checks if a file is a text file using FileTypeUtil.
     */
    private static boolean isTextFile(File file) {
        return FileTypeUtil.isTextFile(file);
    }

    private WritableDatabaseService rootService;
    private WritableDatabaseService mainFolder;
    private Map<String, Map<String, WritableDatabaseService>> sampleFolders; // sample -> (subfolder -> service)
    private String runName;
    private Path outputDir;
    private WritableDatabaseService targetDatabaseService; // The selected folder to import into

    /**
     * Creates a new sample-based importer.
     * @param runName Name for the main folder (e.g., "TaxTriage_output_20240916")
     * @param outputDir Path to TaxTriage output directory
     */
    public SampleBasedImporter(String runName, Path outputDir) {
        this.runName = runName;
        this.outputDir = outputDir;
        this.sampleFolders = new HashMap<>();
        initializeRootService();
    }

    /**
     * Creates a new sample-based importer with a specific target database.
     * @param runName Name for the main folder (e.g., "TaxTriage_output_20240916")
     * @param outputDir Path to TaxTriage output directory
     * @param targetDatabase The database service to import into (e.g., currently selected folder)
     */
    public SampleBasedImporter(String runName, Path outputDir, WritableDatabaseService targetDatabase) {
        this.runName = runName;
        this.outputDir = outputDir;
        this.sampleFolders = new HashMap<>();
        this.targetDatabaseService = targetDatabase;
        if (targetDatabase != null) {
            this.rootService = targetDatabase;
        } else {
            initializeRootService();
        }
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
     * Discovers samples from the TaxTriage output directory.
     * @return Set of sample names
     */
    public Set<String> discoverSamples() {
        Set<String> samples = new HashSet<>();

        try {
            // Look for sample-specific files in various directories
            Path[] checkDirs = {
                outputDir.resolve("kraken2"),
                outputDir.resolve("minimap2"),
                outputDir.resolve("alignment"),
                outputDir.resolve("consensus"),
                outputDir.resolve("download")
            };

            for (Path dir : checkDirs) {
                if (Files.exists(dir)) {
                    Files.list(dir)
                        .filter(Files::isRegularFile)
                        .map(p -> p.getFileName().toString())
                        .forEach(filename -> {
                            // Extract sample name from filename
                            // Patterns: sample.kraken2.report.txt, sample.bam, sample.consensus.fasta
                            String sample = extractSampleName(filename);
                            if (sample != null && !sample.isEmpty()) {
                                samples.add(sample);
                            }
                        });
                }
            }
        } catch (IOException e) {
            logger.log(Level.WARNING, "Error discovering samples", e);
        }

        if (samples.isEmpty()) {
            // Default to a single sample called "results"
            samples.add("results");
        }

        logger.info("Discovered samples: " + samples);
        return samples;
    }

    /**
     * Extracts sample name from a filename.
     */
    private String extractSampleName(String filename) {
        // Remove common suffixes to get sample name
        String[] suffixes = {
            ".kraken2.report.txt",
            ".krakenreport.krona.txt",
            ".consensus.fasta",
            ".dwnld.references.fasta",
            ".dwnld.references.bam",
            ".bam.bai",
            ".bam",
            ".combined.gcfids.txt",
            ".topnames.txt",
            ".toptaxids.txt",
            ".histo.txt",
            ".paths.txt"
        };

        for (String suffix : suffixes) {
            if (filename.endsWith(suffix)) {
                String sample = filename.substring(0, filename.length() - suffix.length());
                // Handle double sample names like "test.test.dwnld.references.bam"
                if (sample.contains(".")) {
                    String[] parts = sample.split("\\.");
                    if (parts.length >= 2 && parts[0].equals(parts[1])) {
                        return parts[0]; // Return single sample name
                    }
                }
                return sample;
            }
        }
        return null;
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

            // Create main folder for this TaxTriage run
            try {
                mainFolder = rootService.createChildFolder(runName);
            } catch (DatabaseServiceException e) {
                logger.warning("Could not create main folder: " + runName + " - " + e.getMessage());
            }

            if (mainFolder == null) {
                // Folder might already exist, try to get it
                mainFolder = rootService.getChildService(runName);
            }

            if (mainFolder == null) {
                logger.warning("Could not create or access main folder: " + runName);
                return false;
            }

            // Discover samples
            Set<String> samples = discoverSamples();

            // Create folder structure for each sample
            for (String sample : samples) {
                createSampleFolders(sample, progressListener);
            }

            logger.info("Folder structure creation completed");
            return true;

        } catch (Exception e) {
            logger.log(Level.SEVERE, "Failed to create folder structure", e);
            return false;
        }
    }

    /**
     * Creates folders for a specific sample based on available data.
     * ONLY creates per-sample folders: Alignments-original, Alignments-deduplicated
     * Shared files (Reports, References, Kraken_Results) are handled at root level
     */
    private void createSampleFolders(String sample, ProgressListener progressListener) {
        logger.info("Creating folders for sample: " + sample);

        Map<String, WritableDatabaseService> folders = new HashMap<>();

        // Create sample folder
        WritableDatabaseService sampleFolder = null;
        try {
            sampleFolder = mainFolder.createChildFolder(sample);
        } catch (DatabaseServiceException e) {
            logger.warning("Could not create sample folder: " + sample + " - " + e.getMessage());
        }

        if (sampleFolder == null) {
            sampleFolder = mainFolder.getChildService(sample);
        }

        if (sampleFolder == null) {
            logger.warning("Could not create sample folder: " + sample);
            return;
        }

        // Define subfolder types and their corresponding output directories
        // ONLY PER-SAMPLE FOLDERS: Alignments-original, Alignments-deduplicated
        Map<String, List<Path>> subfolderMapping = new HashMap<>();

        // Check if deduplication was performed (look for .dedup.bam files)
        boolean hasDeduplicatedFiles = hasDeduplicatedBamFiles(sample,
            Arrays.asList(outputDir.resolve("minimap2"), outputDir.resolve("alignment")));

        if (hasDeduplicatedFiles) {
            logger.info("Deduplication detected - creating separate alignment folders");
            // Create separate folders for original and deduplicated alignments
            subfolderMapping.put("Alignments-original", Arrays.asList(
                outputDir.resolve("minimap2"),
                outputDir.resolve("alignment")
            ));

            subfolderMapping.put("Alignments-deduplicated", Arrays.asList(
                outputDir.resolve("minimap2"),
                outputDir.resolve("alignment")
            ));
        } else {
            logger.info("No deduplication detected - creating single Alignments folder");
            // Just create single Alignments folder if no deduplication
            subfolderMapping.put("Alignments", Arrays.asList(
                outputDir.resolve("minimap2"),
                outputDir.resolve("alignment")
            ));
        }

        // Only create subfolders that have content
        for (Map.Entry<String, List<Path>> entry : subfolderMapping.entrySet()) {
            String folderName = entry.getKey();
            List<Path> sourceDirs = entry.getValue();

            boolean shouldCreateFolder = false;

            // Special handling for alignment folders
            if (folderName.equals("Alignments-original")) {
                // Check for non-deduplicated BAM files
                shouldCreateFolder = hasOriginalBamFiles(sample, sourceDirs);
                logger.fine("Checking for original BAM files: " + shouldCreateFolder);
            } else if (folderName.equals("Alignments-deduplicated")) {
                // Check for deduplicated BAM files
                shouldCreateFolder = hasDeduplicatedBamFiles(sample, sourceDirs);
                logger.fine("Checking for deduplicated BAM files: " + shouldCreateFolder);
            } else if (folderName.equals("Alignments")) {
                // Check for any BAM files
                shouldCreateFolder = hasOriginalBamFiles(sample, sourceDirs) || hasDeduplicatedBamFiles(sample, sourceDirs);
                logger.fine("Checking for any BAM files: " + shouldCreateFolder);
            }

            if (shouldCreateFolder) {
                WritableDatabaseService subfolder = null;
                try {
                    subfolder = sampleFolder.createChildFolder(folderName);
                } catch (DatabaseServiceException e) {
                    logger.warning("Could not create subfolder: " + folderName + " - " + e.getMessage());
                }

                if (subfolder == null) {
                    subfolder = sampleFolder.getChildService(folderName);
                }

                if (subfolder != null) {
                    folders.put(folderName, subfolder);
                    logger.info("Created subfolder: " + sample + "/" + folderName);
                }
            }
        }

        sampleFolders.put(sample, folders);
    }

    /**
     * Checks if there are original (non-deduplicated) BAM files for the sample.
     */
    private boolean hasOriginalBamFiles(String sample, List<Path> dirs) {
        logger.fine("Checking for original BAM files for sample: " + sample);
        for (Path dir : dirs) {
            if (Files.exists(dir)) {
                logger.fine("Checking directory: " + dir);
                try {
                    List<Path> originalFiles = Files.list(dir)
                        .filter(Files::isRegularFile)
                        .filter(p -> {
                            String filename = p.getFileName().toString();
                            // Original BAM files: contain sample name, end with .bam, but NOT .dedup.bam
                            boolean isOriginal = filename.contains(sample) &&
                                                 filename.endsWith(".bam") &&
                                                 !filename.endsWith(".dedup.bam");
                            if (isOriginal) {
                                logger.fine("Found original file: " + filename);
                            }
                            return isOriginal;
                        })
                        .collect(Collectors.toList());

                    if (!originalFiles.isEmpty()) {
                        logger.info("Found " + originalFiles.size() + " original BAM file(s)");
                        return true;
                    }
                } catch (IOException e) {
                    logger.warning("Error checking for original BAM files: " + e.getMessage());
                }
            }
        }
        logger.fine("No original BAM files found");
        return false;
    }

    /**
     * Checks if there are deduplicated BAM files for the sample.
     */
    private boolean hasDeduplicatedBamFiles(String sample, List<Path> dirs) {
        logger.fine("Checking for deduplicated BAM files for sample: " + sample);
        for (Path dir : dirs) {
            if (Files.exists(dir)) {
                logger.fine("Checking directory: " + dir);
                try {
                    List<Path> dedupFiles = Files.list(dir)
                        .filter(Files::isRegularFile)
                        .filter(p -> {
                            String filename = p.getFileName().toString();
                            boolean isDedup = filename.contains(sample) && filename.endsWith(".dedup.bam");
                            if (isDedup) {
                                logger.fine("Found deduplicated file: " + filename);
                            }
                            return isDedup;
                        })
                        .collect(java.util.stream.Collectors.toList());

                    if (!dedupFiles.isEmpty()) {
                        logger.info("Found " + dedupFiles.size() + " deduplicated BAM file(s)");
                        return true;
                    }
                } catch (IOException e) {
                    logger.warning("Error checking for deduplicated BAM files: " + e.getMessage());
                }
            }
        }
        System.out.println("    No deduplicated BAM files found");
        return false;
    }

    /**
     * Checks if there is content for a sample in the given directories.
     */
    private boolean hasContentForSample(String sample, List<Path> dirs) {
        for (Path dir : dirs) {
            if (Files.exists(dir)) {
                try {
                    boolean hasContent = Files.list(dir)
                        .filter(Files::isRegularFile)
                        .anyMatch(p -> {
                            String filename = p.getFileName().toString();
                            return filename.contains(sample) ||
                                   filename.equals("assembly_summary_refseq.txt") ||
                                   filename.startsWith("complete.") ||
                                   filename.startsWith("multiqc");
                        });
                    if (hasContent) return true;
                } catch (IOException e) {
                    logger.warning("Error checking directory: " + dir);
                }
            }
        }
        return false;
    }

    /**
     * Imports files for a specific sample to appropriate subfolders.
     */
    public List<AnnotatedPluginDocument> importSampleFiles(String sample, ProgressListener progressListener)
            throws IOException, DocumentImportException {

        List<AnnotatedPluginDocument> allImported = new ArrayList<>();
        Map<String, WritableDatabaseService> folders = sampleFolders.get(sample);

        if (folders == null || folders.isEmpty()) {
            logger.warning("No folders created for sample: " + sample);
            return allImported;
        }

        System.out.println("\n===== IMPORTING FILES FOR SAMPLE: " + sample + " =====");

        // First, check if we need to download GenBank references for BAM files
        downloadGenBankReferencesForBamFiles(sample, progressListener);

        // Import references (GenBank files from minimap2 folder where BAM files are)
        if (folders.containsKey("References")) {
            // Check for GenBank downloads in minimap2 folder (co-located with BAM files)
            Path minimap2Dir = outputDir.resolve("minimap2");
            if (Files.exists(minimap2Dir)) {
                // Debug: Show all files in minimap2 directory
                System.out.println("  === GenBank Import Debug ===");
                System.out.println("  Sample name: " + sample);
                System.out.println("  Minimap2 directory: " + minimap2Dir);
                try {
                    List<Path> allFiles = Files.list(minimap2Dir).collect(Collectors.toList());
                    System.out.println("  Total files in minimap2: " + allFiles.size());
                    for (Path p : allFiles) {
                        String filename = p.getFileName().toString();
                        if (filename.endsWith(".gb") || filename.endsWith(".gbk") || filename.endsWith(".genbank")) {
                            System.out.println("    GenBank file found: " + filename);
                        }
                    }
                } catch (IOException e) {
                    System.out.println("  Error listing directory: " + e.getMessage());
                }

                // First, fix any GenBank files to ensure proper field matching
                System.out.println("  Checking and fixing GenBank files in minimap2 folder...");
                int fixedCount = GenBankFileFixer.fixGenBankFilesInDirectory(minimap2Dir);
                if (fixedCount > 0) {
                    System.out.println("  Fixed " + fixedCount + " GenBank file(s) for proper field matching");
                }

                // Use special method for GenBank files since they're named by accession, not sample
                List<File> gbFiles = findGenBankReferenceFiles(Arrays.asList(minimap2Dir));

                System.out.println("  GenBank reference files found: " + gbFiles.size());
                for (File f : gbFiles) {
                    System.out.println("    - " + f.getName());
                }

                for (File file : gbFiles) {
                    try {
                        // Validate and fix the GenBank file if needed
                        System.out.println("  Processing GenBank file: " + file.getName());
                        if (!GenBankFileFixer.validateGenBankFile(file)) {
                            System.out.println("    File has mismatched fields, attempting to fix...");
                            boolean fixed = GenBankFileFixer.fixGenBankFile(file);
                            if (fixed) {
                                System.out.println("    Successfully fixed field mismatches");
                                // Re-validate after fixing
                                if (!GenBankFileFixer.validateGenBankFile(file)) {
                                    System.out.println("    WARNING: File still has issues after fix attempt");
                                }
                            } else {
                                System.out.println("    WARNING: Could not fix field mismatches");
                                // Continue anyway, import might still work
                            }
                        } else {
                            System.out.println("    File validation passed");
                        }

                        System.out.println("  Importing GenBank from minimap2 to References: " + file.getName());
                        // First import the documents
                        List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(file, progressListener);

                        if (docs != null && !docs.isEmpty()) {
                            WritableDatabaseService referencesFolder = folders.get("References");
                            System.out.println("    Imported " + docs.size() + " document(s) from file");

                            // Add each document to the References folder
                            for (AnnotatedPluginDocument doc : docs) {
                                try {
                                    AnnotatedPluginDocument copiedDoc = referencesFolder.addDocumentCopy(doc, ProgressListener.EMPTY);
                                    if (copiedDoc != null) {
                                        allImported.add(copiedDoc);
                                        System.out.println("      - Added to References: " + copiedDoc.getName());
                                    }
                                } catch (Exception e) {
                                    System.out.println("      - Failed to add to References: " + doc.getName() + " - " + e.getMessage());
                                    logger.log(Level.WARNING, "Failed to add document to References folder", e);
                                }
                            }
                        } else {
                            System.out.println("    Warning: No documents imported from " + file.getName());
                        }
                    } catch (Exception e) {
                        logger.log(Level.WARNING, "Failed to import GenBank file: " + file.getName(), e);
                        System.out.println("    Error: " + e.getMessage());
                    }
                }
                logger.info("Imported " + gbFiles.size() + " GenBank reference(s) for sample: " + sample);
            }

            // Also check legacy genbank_downloads folder if it exists
            Path genBankDir = outputDir.resolve("genbank_downloads");
            if (Files.exists(genBankDir)) {
                // Fix GenBank files in this directory too
                System.out.println("  Checking and fixing GenBank files in genbank_downloads folder...");
                int fixedCount = GenBankFileFixer.fixGenBankFilesInDirectory(genBankDir);
                if (fixedCount > 0) {
                    System.out.println("  Fixed " + fixedCount + " GenBank file(s)");
                }

                // Use special method for GenBank files since they're named by accession, not sample
                List<File> gbFiles = findGenBankReferenceFiles(Arrays.asList(genBankDir));

                for (File file : gbFiles) {
                    try {
                        System.out.println("  Importing GenBank from genbank_downloads to References: " + file.getName());
                        // First import the documents
                        List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(file, progressListener);

                        if (docs != null && !docs.isEmpty()) {
                            WritableDatabaseService referencesFolder = folders.get("References");
                            // Add each document to the References folder
                            for (AnnotatedPluginDocument doc : docs) {
                                try {
                                    AnnotatedPluginDocument copiedDoc = referencesFolder.addDocumentCopy(doc, ProgressListener.EMPTY);
                                    if (copiedDoc != null) {
                                        allImported.add(copiedDoc);
                                        System.out.println("    - Added to References: " + copiedDoc.getName());
                                    }
                                } catch (Exception ex) {
                                    logger.log(Level.WARNING, "Failed to add document to References folder: " + doc.getName(), ex);
                                }
                            }
                        }
                    } catch (Exception e) {
                        logger.log(Level.WARNING, "Failed to import GenBank file: " + file.getName(), e);
                    }
                }
            }

            // Skip workflow-generated FASTA files (test.dwnld.references.fasta)
            // as these don't have annotations
        }

        // Import ALL text files first using direct import
        importAllTextFiles(sample, folders, allImported, progressListener);

        // Import HTML reports (non-text)
        if (folders.containsKey("Reports")) {
            List<File> htmlFiles = findFilesForSample(sample,
                Arrays.asList(outputDir.resolve("report"), outputDir.resolve("pipeline_info")),
                Arrays.asList(".html"));

            for (File file : htmlFiles) {
                try {
                    List<AnnotatedPluginDocument> docs = importToFolder(file, folders.get("Reports"), progressListener);
                    allImported.addAll(docs);
                    System.out.println("  Imported HTML: " + file.getName());
                } catch (Exception e) {
                    logger.warning("Failed to import HTML file: " + file.getName());
                }
            }
        }

        // Collect all imported reference documents for BAM linking
        List<AnnotatedPluginDocument> importedReferences = new ArrayList<>();
        if (folders.containsKey("References")) {
            try {
                WritableDatabaseService referencesFolder = folders.get("References");
                // Note: We track the references we imported earlier in allImported
                // Filter to get only reference sequences (GenBank documents)
                for (AnnotatedPluginDocument doc : allImported) {
                    // Check if this is a sequence document that could be a reference
                    if (doc.getDocumentClass() != null &&
                        (doc.getDocumentClass().getName().contains("NucleotideSequence") ||
                         doc.getDocumentClass().getName().contains("DefaultNucleotideSequence"))) {
                        importedReferences.add(doc);
                    }
                }
                System.out.println("  Collected " + importedReferences.size() + " reference sequences for BAM import");
                for (AnnotatedPluginDocument ref : importedReferences) {
                    System.out.println("    - Reference: " + ref.getName());
                }
            } catch (Exception e) {
                logger.warning("Could not collect reference documents: " + e.getMessage());
            }
        }

        // Import BAM files using proper DocumentFileImporter API
        // Check if we have separate folders for original and deduplicated alignments
        boolean hasDeduplicationFolders = folders.containsKey("Alignments-original") &&
                                          folders.containsKey("Alignments-deduplicated");

        if (folders.containsKey("Alignments") || hasDeduplicationFolders || folders.containsKey("References")) {
            List<File> bamFiles = findFilesForSample(sample,
                Arrays.asList(outputDir.resolve("minimap2"), outputDir.resolve("alignment")),
                Arrays.asList(".bam"));

            // Find GenBank reference files co-located with BAM files
            Path minimap2Dir = outputDir.resolve("minimap2");
            List<File> colocatedReferences = new ArrayList<>();

            if (Files.exists(minimap2Dir)) {
                try {
                    Files.list(minimap2Dir)
                        .filter(p -> {
                            String name = p.getFileName().toString().toLowerCase();
                            return name.endsWith(".gb") || name.endsWith(".gbk") || name.endsWith(".genbank");
                        })
                        .map(Path::toFile)
                        .forEach(colocatedReferences::add);

                    if (!colocatedReferences.isEmpty()) {
                        System.out.println("  Found " + colocatedReferences.size() + " GenBank reference(s) co-located with BAM files");
                    }
                } catch (IOException e) {
                    logger.warning("Error scanning for co-located references: " + e.getMessage());
                }
            }

            // Separate BAM files into original and deduplicated
            List<File> originalBamFiles = new ArrayList<>();
            List<File> dedupBamFiles = new ArrayList<>();

            System.out.println("\n  Categorizing " + bamFiles.size() + " BAM file(s):");
            for (File bamFile : bamFiles) {
                String filename = bamFile.getName();
                if (filename.endsWith(".dedup.bam")) {
                    dedupBamFiles.add(bamFile);
                    System.out.println("    - " + filename + " → deduplicated");
                } else if (filename.endsWith(".bam") && !filename.contains(".dedup")) {
                    originalBamFiles.add(bamFile);
                    System.out.println("    - " + filename + " → original");
                }
            }

            System.out.println("\n  Folder structure: hasDeduplicationFolders=" + hasDeduplicationFolders);
            System.out.println("    Alignments-original folder exists: " + folders.containsKey("Alignments-original"));
            System.out.println("    Alignments-deduplicated folder exists: " + folders.containsKey("Alignments-deduplicated"));
            System.out.println("    Alignments folder exists: " + folders.containsKey("Alignments"));

            // Import original BAM files
            if (!originalBamFiles.isEmpty() && hasDeduplicationFolders) {
                WritableDatabaseService originalFolder = folders.get("Alignments-original");
                if (originalFolder != null) {
                    System.out.println("\n  Importing " + originalBamFiles.size() + " original BAM file(s) to Alignments-original folder:");
                    for (File bamFile : originalBamFiles) {
                        importBamFile(bamFile, originalFolder, colocatedReferences, progressListener, allImported);
                    }
                } else {
                    System.out.println("\n  ERROR: Alignments-original folder not found, cannot import original BAM files");
                }
            }

            // Import deduplicated BAM files
            if (!dedupBamFiles.isEmpty() && hasDeduplicationFolders) {
                WritableDatabaseService dedupFolder = folders.get("Alignments-deduplicated");
                if (dedupFolder != null) {
                    System.out.println("\n  Importing " + dedupBamFiles.size() + " deduplicated BAM file(s) to Alignments-deduplicated folder:");
                    for (File bamFile : dedupBamFiles) {
                        importBamFile(bamFile, dedupFolder, colocatedReferences, progressListener, allImported);
                    }
                } else {
                    System.out.println("\n  ERROR: Alignments-deduplicated folder not found, cannot import deduplicated BAM files");
                }
            }

            // If no deduplication folders, import all BAM files to single Alignments folder
            if (!hasDeduplicationFolders) {
                WritableDatabaseService targetFolder = folders.get("Alignments");
                if (targetFolder == null) {
                    targetFolder = folders.get("References");
                }

                if (targetFolder == null) {
                    System.out.println("    Error: No suitable folder for BAM import (need Alignments or References folder)");
                } else {
                    System.out.println("\n  Importing BAM files to " +
                                     (folders.get("Alignments") == targetFolder ? "Alignments" : "References") + " folder:");
                    for (File bamFile : bamFiles) {
                        importBamFile(bamFile, targetFolder, colocatedReferences, progressListener, allImported);
                    }
                }
            }
        }

        System.out.println("Imported " + allImported.size() + " documents for sample: " + sample);
        System.out.println("=====================================\n");

        return allImported;
    }

    /**
     * Helper method to import a single BAM file
     */
    private void importBamFile(File bamFile, WritableDatabaseService targetFolder,
                               List<File> colocatedReferences, ProgressListener progressListener,
                               List<AnnotatedPluginDocument> allImported) {
        System.out.println("      Processing: " + bamFile.getName());

        try {
            // Use the proper BAM importer with references
            List<AnnotatedPluginDocument> bamDocs;

            if (!colocatedReferences.isEmpty()) {
                System.out.println("        Using ProperBamImporter with " + colocatedReferences.size() + " reference file(s)");
                bamDocs = ProperBamImporter.importBamWithReferences(
                    bamFile,
                    colocatedReferences,
                    targetFolder,
                    progressListener
                );
            } else {
                System.out.println("        Using ProperBamImporter without explicit reference files");
                System.out.println("        (Importer will look for references in the database)");
                bamDocs = ProperBamImporter.importBam(
                    bamFile,
                    targetFolder,
                    progressListener
                );
            }

            if (bamDocs != null && !bamDocs.isEmpty()) {
                allImported.addAll(bamDocs);
                System.out.println("        Successfully imported " + bamDocs.size() + " document(s) from BAM");
            } else {
                System.out.println("        Warning: No documents imported from BAM file");
            }

        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to import BAM file: " + bamFile.getName(), e);
            System.out.println("        Error importing BAM: " + e.getMessage());
            e.printStackTrace();
        }
    }

    /**
     * Imports ALL text files using direct import method.
     * This ensures all text files are properly imported without parsing.
     */
    private void importAllTextFiles(String sample, Map<String, WritableDatabaseService> folders,
                                    List<AnnotatedPluginDocument> allImported,
                                    ProgressListener progressListener) {

        System.out.println("\n  === Scanning for all text files ===");

        // Collect all text files from all directories
        List<File> allTextFiles = new ArrayList<>();
        Set<String> processedFiles = new HashSet<>();

        try {
            // Scan root output directory
            if (Files.exists(outputDir)) {
                Files.walk(outputDir, 1) // Only immediate files
                    .filter(Files::isRegularFile)
                    .map(Path::toFile)
                    .filter(file -> isTextFile(file))
                    .forEach(file -> {
                        if (!processedFiles.contains(file.getName())) {
                            allTextFiles.add(file);
                            processedFiles.add(file.getName());
                            System.out.println("    Found text file: " + file.getName());
                        }
                    });
            }

            // Scan all subdirectories
            Path[] scanDirs = {
                outputDir.resolve("kraken2"),
                outputDir.resolve("kreport"),
                outputDir.resolve("mergedsubspecies"),
                outputDir.resolve("mergedkrakenreport"),
                outputDir.resolve("report"),
                outputDir.resolve("top"),
                outputDir.resolve("combine"),
                outputDir.resolve("count"),
                outputDir.resolve("alignment"),
                outputDir.resolve("get"),
                outputDir.resolve("pipeline_info")
            };

            for (Path dir : scanDirs) {
                if (Files.exists(dir)) {
                    Files.list(dir)
                        .filter(Files::isRegularFile)
                        .map(Path::toFile)
                        .filter(file -> isTextFile(file))
                        .forEach(file -> {
                            if (!processedFiles.contains(file.getName())) {
                                allTextFiles.add(file);
                                processedFiles.add(file.getName());
                                System.out.println("    Found text file in " + dir.getFileName() + ": " + file.getName());
                            }
                        });
                }
            }

        } catch (IOException e) {
            logger.log(Level.WARNING, "Error scanning for text files", e);
        }

        System.out.println("  Total text files found: " + allTextFiles.size());
        System.out.println("\n  === Importing text files ===");

        // Now import all text files
        for (File textFile : allTextFiles) {
            String filename = textFile.getName();

            // Check if this file is relevant to our sample or is a general file
            if (filename.contains(sample) ||
                filename.startsWith("taxtriage_") ||
                filename.contains("kraken") ||
                filename.contains("krona") ||
                filename.contains("report") ||
                filename.contains("count") ||
                filename.contains("complete") ||
                filename.contains("multiqc")) {

                // Determine target folder
                WritableDatabaseService targetFolder = null;
                String folderName = "";

                if (filename.contains("kraken") || filename.contains("krona")) {
                    targetFolder = folders.get("Kraken_Results");
                    folderName = "Kraken_Results";
                } else {
                    targetFolder = folders.get("Reports");
                    folderName = "Reports";
                }

                if (targetFolder != null) {
                    try {
                        System.out.println("  Importing to " + folderName + ": " + filename);

                        // Read text file content using Files.readAllBytes
                        try {
                            byte[] fileBytes = Files.readAllBytes(textFile.toPath());
                            String content = new String(fileBytes, "UTF-8");

                            // Determine file format and result type
                            String fileFormat = determineFileFormat(filename);
                            String resultType = determineResultType(filename);

                            // Create TaxTriageResultDocument with the content
                            TaxTriageResultDocument resultDoc = new TaxTriageResultDocument(
                                filename, content, resultType, textFile.getAbsolutePath(), fileFormat);

                            // Convert to AnnotatedPluginDocument using DocumentUtilities
                            AnnotatedPluginDocument annotatedDoc = DocumentUtilities.createAnnotatedPluginDocument(resultDoc);

                            // Add to target folder using addDocumentCopy
                            AnnotatedPluginDocument copiedDoc = targetFolder.addDocumentCopy(annotatedDoc, ProgressListener.EMPTY);

                            if (copiedDoc != null) {
                                allImported.add(copiedDoc);
                                System.out.println("      Successfully imported text file: " + filename + " (" + fileFormat + ")");
                            } else {
                                System.out.println("      Warning: Failed to copy document to folder: " + filename);
                            }

                        } catch (IOException ioException) {
                            logger.log(Level.WARNING, "Failed to read text file: " + filename, ioException);
                            System.out.println("      Error reading file: " + ioException.getMessage());
                        } catch (Exception docException) {
                            logger.log(Level.WARNING, "Failed to create document from text file: " + filename, docException);
                            System.out.println("      Error creating document: " + docException.getMessage());
                        }

                    } catch (Exception e) {
                        logger.warning("Failed to import text file: " + filename + " - " + e.getMessage());
                        System.out.println("    ERROR: Failed to import " + filename + ": " + e.getMessage());
                    }
                }
            }
        }

        System.out.println("  === Text file import complete ===");
        System.out.println();
    }

    /**
     * Finds files for a specific sample in given directories.
     */
    private List<File> findFilesForSample(String sample, List<Path> dirs, List<String> extensions) {
        List<File> files = new ArrayList<>();

        for (Path dir : dirs) {
            if (Files.exists(dir)) {
                try {
                    Files.list(dir)
                        .filter(Files::isRegularFile)
                        .filter(p -> {
                            String name = p.getFileName().toString();
                            // Check if file is for this sample or a general output file
                            boolean isSampleFile = name.contains(sample) ||
                                                  name.startsWith("taxtriage_") || // Include taxtriage_ prefixed files
                                                  name.contains(".krakenreport") || // Include krakenreport files
                                                  name.contains(".krona") || // Include krona files
                                                  name.startsWith("complete.") ||
                                                  name.startsWith("multiqc") ||
                                                  name.equals("assembly_summary_refseq.txt") ||
                                                  name.equals("count.txt");
                            // Check extension
                            boolean hasValidExt = extensions.stream().anyMatch(name::endsWith);
                            return isSampleFile && hasValidExt;
                        })
                        .map(Path::toFile)
                        .forEach(files::add);
                } catch (IOException e) {
                    logger.warning("Error reading directory: " + dir);
                }
            }
        }

        return files;
    }

    /**
     * Finds GenBank reference files in directories.
     * These files are named by accession number, not sample name.
     */
    private List<File> findGenBankReferenceFiles(List<Path> dirs) {
        List<File> files = new ArrayList<>();

        for (Path dir : dirs) {
            if (Files.exists(dir)) {
                try {
                    Files.list(dir)
                        .filter(Files::isRegularFile)
                        .filter(p -> {
                            String name = p.getFileName().toString().toLowerCase();
                            // GenBank files have specific extensions
                            return name.endsWith(".gb") || name.endsWith(".gbk") || name.endsWith(".genbank");
                        })
                        .map(Path::toFile)
                        .forEach(files::add);
                } catch (IOException e) {
                    logger.warning("Error reading directory for GenBank files: " + dir);
                }
            }
        }

        return files;
    }

    /**
     * Imports a file to a specific folder.
     */
    private List<AnnotatedPluginDocument> importToFolder(File file, WritableDatabaseService folder,
                                                         ProgressListener progressListener)
                                                         throws IOException, DocumentImportException {

        System.out.println("  Importing: " + file.getName());

        try {
            List<AnnotatedPluginDocument> docs = PluginUtilities.importDocumentsToDatabase(
                file, folder, progressListener);

            if (docs != null && !docs.isEmpty()) {
                logger.info("Successfully imported " + docs.size() + " document(s) from " + file.getName());
            }

            return docs != null ? docs : new ArrayList<>();

        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to import file: " + file.getName(), e);
            // Don't fail completely, just skip this file
            return new ArrayList<>();
        }
    }

    /**
     * Gets all discovered samples.
     */
    public Set<String> getSamples() {
        return new HashSet<>(sampleFolders.keySet());
    }

    /**
     * Gets the root database service for importing shared files.
     * @return WritableDatabaseService for the root database
     */
    public WritableDatabaseService getRootDatabaseService() {
        if (rootService == null) {
            initializeRootService();
        }
        return rootService;
    }

    /**
     * Gets the main folder (TaxTriage run folder) for importing shared files.
     * @return WritableDatabaseService for the main folder
     */
    public WritableDatabaseService getMainFolder() {
        return mainFolder;
    }

    /**
     * Downloads GenBank references for BAM files and saves them in the minimap2 folder.
     * This co-locates references with BAM files for better import success.
     */
    private void downloadGenBankReferencesForBamFiles(String sample, ProgressListener progressListener) {
        Path minimap2Dir = outputDir.resolve("minimap2");
        if (!Files.exists(minimap2Dir)) {
            return;
        }

        System.out.println("\n  === Checking for BAM references to download ===");

        // Find BAM files for this sample
        List<File> bamFiles = findFilesForSample(sample,
            Arrays.asList(minimap2Dir),
            Arrays.asList(".bam"));

        if (bamFiles.isEmpty()) {
            System.out.println("  No BAM files found, skipping GenBank download");
            return;
        }

        // Extract references from BAM files
        Set<String> allAccessions = new HashSet<>();
        for (File bamFile : bamFiles) {
            try {
                List<BamReferenceExtractor.ReferenceInfo> refs =
                    BamReferenceExtractor.extractReferences(bamFile);
                for (BamReferenceExtractor.ReferenceInfo ref : refs) {
                    if (ref.accession != null && !ref.accession.isEmpty()) {
                        allAccessions.add(ref.accession);
                        System.out.println("  Found reference in BAM: " + ref.name + " (" + ref.accession + ")");
                    }
                }
            } catch (Exception e) {
                logger.warning("Could not extract references from BAM: " + bamFile.getName());
            }
        }

        if (allAccessions.isEmpty()) {
            System.out.println("  No accessions found in BAM files");
            return;
        }

        // Check which GenBank files already exist
        List<String> toDownload = new ArrayList<>();
        for (String accession : allAccessions) {
            File gbFile = minimap2Dir.resolve(accession + ".gb").toFile();
            if (!gbFile.exists()) {
                toDownload.add(accession);
            } else {
                System.out.println("  GenBank file already exists: " + gbFile.getName());
            }
        }

        if (toDownload.isEmpty()) {
            System.out.println("  All GenBank references already downloaded");
            return;
        }

        // Download missing GenBank files to minimap2 folder
        System.out.println("\n  Downloading " + toDownload.size() + " GenBank reference(s) to minimap2 folder...");
        try {
            Map<String, File> downloaded = NCBIReferenceDownloader.downloadGenBankFiles(
                toDownload, minimap2Dir, progressListener);

            for (Map.Entry<String, File> entry : downloaded.entrySet()) {
                System.out.println("    Downloaded: " + entry.getValue().getName());
            }

            if (downloaded.size() < toDownload.size()) {
                System.out.println("  Warning: Only " + downloaded.size() + " of " +
                                  toDownload.size() + " files downloaded successfully");
            }
        } catch (Exception e) {
            logger.log(Level.WARNING, "Error downloading GenBank files", e);
            System.out.println("  Error downloading GenBank files: " + e.getMessage());
        }

        System.out.println("  === GenBank download complete ===");
        System.out.println();
    }

    /**
     * Determines the file format based on filename extension.
     * @param filename The name of the file
     * @return File format string (e.g., "TXT", "TSV", "CSV", "HTML", "JSON")
     */
    private String determineFileFormat(String filename) {
        String lowerName = filename.toLowerCase();

        if (lowerName.endsWith(".tsv")) {
            return "TSV";
        } else if (lowerName.endsWith(".csv")) {
            return "CSV";
        } else if (lowerName.endsWith(".html") || lowerName.endsWith(".htm")) {
            return "HTML";
        } else if (lowerName.endsWith(".json")) {
            return "JSON";
        } else if (lowerName.endsWith(".log")) {
            return "LOG";
        } else {
            return "TXT"; // Default to TXT for .txt and other text files
        }
    }

    /**
     * Determines the result type based on filename patterns.
     * @param filename The name of the file
     * @return Result type string describing the content
     */
    private String determineResultType(String filename) {
        String lowerName = filename.toLowerCase();

        if (lowerName.contains("kraken") && lowerName.contains("report")) {
            return "Kraken Report";
        } else if (lowerName.contains("krona")) {
            return "Krona Report";
        } else if (lowerName.contains("count")) {
            return "Count Data";
        } else if (lowerName.contains("top")) {
            return "Top Results";
        } else if (lowerName.contains("combine")) {
            return "Combined Results";
        } else if (lowerName.contains("multiqc")) {
            return "MultiQC Report";
        } else if (lowerName.contains("taxtriage")) {
            return "TaxTriage Analysis";
        } else if (lowerName.contains("assembly")) {
            return "Assembly Summary";
        } else if (lowerName.contains("complete")) {
            return "Completion Report";
        } else {
            return "Analysis Result"; // Default generic type
        }
    }
}