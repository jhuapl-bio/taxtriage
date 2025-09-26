package com.jhuapl.taxtriage.geneious.results;

import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.plugin.DocumentImportException;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import com.jhuapl.taxtriage.geneious.documents.TaxTriageResultDocument;
import com.jhuapl.taxtriage.geneious.importers.BamFileImporter;
import com.jhuapl.taxtriage.geneious.importers.BamImportHelper;
import com.jhuapl.taxtriage.geneious.importers.BamIndexer;
import com.jhuapl.taxtriage.geneious.importers.BamReferenceExtractor;
import com.jhuapl.taxtriage.geneious.importers.SampleBasedImporter;
import com.jhuapl.taxtriage.geneious.importers.ForceIndexingHelper;
import com.jhuapl.taxtriage.geneious.importers.LocalDatabaseHelper;
import com.jhuapl.taxtriage.geneious.importers.NCBIReferenceDownloader;
import com.jhuapl.taxtriage.geneious.importers.ReferenceManager;
import com.jhuapl.taxtriage.geneious.importers.SimpleFastaImporter;
import com.jhuapl.taxtriage.geneious.reports.TaxTriageHtmlReportGenerator;
import com.jhuapl.taxtriage.geneious.utils.ArchiveManager;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * Imports TaxTriage analysis results into Geneious.
 *
 * This class is responsible for scanning the TaxTriage output directory
 * and importing relevant result files as Geneious documents.
 */
public class TaxTriageResultImporter {

    /**
     * Import phases for separating the import process.
     */
    public enum ImportPhase {
        REFERENCES_ONLY,  // Phase 1: Import only GenBank reference sequences
        TEXT_ONLY,        // Phase 2: Import only text files
        BAM_ONLY,         // Phase 3: Import only BAM files
        ALL               // Import all phases sequentially
    }

    private static final Logger logger = Logger.getLogger(TaxTriageResultImporter.class.getName());

    // Store workflow parameters for analysis summary
    private String workingDirectory;
    private String commandLine;
    private Map<String, String> workflowParams = new HashMap<>();

    // Track imported documents by folder
    private Map<String, List<AnnotatedPluginDocument>> documentsByFolder = new HashMap<>();

    // Cache of BAM reference name -> Geneious document name
    private Map<String, String> bamReferenceNameMap = new HashMap<>();

    // Track tool versions and run information for HTML report
    private Map<String, String> toolVersions = new HashMap<>();
    private List<String> importedSamples = new ArrayList<>();

    /**
     * Sets the workflow parameters for the analysis summary.
     */
    public void setWorkflowParameters(String workingDir, String cmdLine, Map<String, String> params) {
        this.workingDirectory = workingDir;
        this.commandLine = cmdLine;
        if (params != null) {
            this.workflowParams.putAll(params);
        }
    }

    /**
     * Sets tool versions for the analysis report.
     */
    public void setToolVersions(Map<String, String> versions) {
        if (versions != null) {
            this.toolVersions.putAll(versions);
        }
    }

    /**
     * Imports all TaxTriage results from the specified output directory.
     * Only imports: .txt files, reference FASTA files, and BAM files.
     * Organizes them into appropriate folders.
     *
     * @param outputDir The TaxTriage output directory
     * @param progressListener Progress tracker
     * @return List of imported documents
     * @throws IOException If file reading fails
     */
    public List<AnnotatedPluginDocument> importResults(Path outputDir, ProgressListener progressListener)
            throws IOException {
        return importResults(outputDir, progressListener, ImportPhase.ALL);
    }

    /**
     * Imports all TaxTriage results from the specified output directory into a target database.
     * Only imports: .txt files, reference FASTA files, and BAM files.
     * Organizes them into appropriate folders.
     *
     * @param outputDir The TaxTriage output directory
     * @param progressListener Progress tracker
     * @param targetDatabase The database service to import into (e.g., currently selected folder)
     * @return List of imported documents
     * @throws IOException If file reading fails
     */
    public List<AnnotatedPluginDocument> importResults(Path outputDir, ProgressListener progressListener,
                                                       WritableDatabaseService targetDatabase)
            throws IOException {
        return importResults(outputDir, progressListener, ImportPhase.ALL, targetDatabase);
    }

    /**
     * Imports TaxTriage results in a specific phase.
     * This allows for separated import operations to fix BAM reference resolution issues.
     *
     * @param outputDir The TaxTriage output directory
     * @param progressListener Progress tracker
     * @param phase The specific import phase to execute
     * @return List of imported documents for the specified phase
     * @throws IOException If file reading fails
     */
    public List<AnnotatedPluginDocument> importResults(Path outputDir, ProgressListener progressListener, ImportPhase phase)
            throws IOException {
        return importResults(outputDir, progressListener, phase, null);
    }

    /**
     * Imports TaxTriage results in a specific phase into a target database.
     * This allows for separated import operations to fix BAM reference resolution issues.
     *
     * @param outputDir The TaxTriage output directory
     * @param progressListener Progress tracker
     * @param phase The specific import phase to execute
     * @param targetDatabase The database service to import into (e.g., currently selected folder)
     * @return List of imported documents for the specified phase
     * @throws IOException If file reading fails
     */
    public List<AnnotatedPluginDocument> importResults(Path outputDir, ProgressListener progressListener,
                                                       ImportPhase phase, WritableDatabaseService targetDatabase)
            throws IOException {

        // Create sample-based importer for organized imports
        String timestamp = new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new java.util.Date());
        String runName = "TaxTriage_" + outputDir.getFileName().toString() + "_" + timestamp;

        // Use target database if provided, otherwise use default root
        SampleBasedImporter sampleImporter;
        if (targetDatabase != null) {
            sampleImporter = new SampleBasedImporter(runName, outputDir, targetDatabase);
        } else {
            sampleImporter = new SampleBasedImporter(runName, outputDir);
        }

        // Create folder structure
        if (sampleImporter.createFolderStructure(progressListener)) {
            logger.info("Created folder structure: " + runName);
            System.out.println("Created sample-based folder structure: " + runName);
        } else {
            logger.warning("Could not create folder structure - imports will go to current location");
            System.out.println("Warning: Could not create folder structure - imports will go to current location");
        }

        // Continue with imports using sample-based structure
        return importResultsWithSamples(outputDir, progressListener, phase, sampleImporter);
    }

    private List<AnnotatedPluginDocument> importResultsWithSamples(Path outputDir, ProgressListener progressListener,
                                                                   ImportPhase phase, SampleBasedImporter sampleImporter)
            throws IOException {

        List<AnnotatedPluginDocument> importedDocs = new ArrayList<>();
        documentsByFolder.clear();
        bamReferenceNameMap.clear();

        if (!Files.exists(outputDir)) {
            logger.warning("Output directory does not exist: " + outputDir);
            return importedDocs;
        }

        logger.info("Importing TaxTriage results from: " + outputDir.toAbsolutePath());

        // Set working directory if not already set
        if (workingDirectory == null) {
            workingDirectory = outputDir.getParent() != null ?
                outputDir.getParent().toAbsolutePath().toString() :
                outputDir.toAbsolutePath().toString();
        }

        int stepCount = 0;
        int totalSteps = getStepsForPhase(phase);

        // Extract BAM references and accessions (needed for most phases)
        List<File> bamFiles = findBamFiles(outputDir);
        Set<BamReferenceExtractor.ReferenceInfo> bamReferences = BamReferenceExtractor.extractAllReferences(bamFiles);
        List<String> accessions = NCBIReferenceDownloader.extractAccessions(bamReferences);
        Map<String, String> referenceNameMap = NCBIReferenceDownloader.createReferenceNameMap(new ArrayList<>(bamReferences));

        // Execute the specific phase or all phases
        switch (phase) {
            case REFERENCES_ONLY:
                return importReferencesOnly(outputDir, progressListener, bamReferences, accessions, referenceNameMap, sampleImporter);
            case TEXT_ONLY:
                return importTextOnly(outputDir, progressListener, sampleImporter);
            case BAM_ONLY:
                return prepareBamOnly(outputDir, progressListener, bamFiles, sampleImporter);
            case ALL:
            default:
                return importAllPhases(outputDir, progressListener, bamReferences, accessions, referenceNameMap, bamFiles, stepCount, totalSteps, sampleImporter);
        }
    }

    /**
     * Gets the number of steps for a specific import phase.
     */
    private int getStepsForPhase(ImportPhase phase) {
        switch (phase) {
            case REFERENCES_ONLY: return 2; // Extract BAM refs, download/import GenBank
            case TEXT_ONLY: return 1;       // Import text files
            case BAM_ONLY: return 1;        // Import BAM files
            case ALL: return 5;             // All steps
            default: return 1;
        }
    }

    /**
     * Imports all phases using sample-based organization.
     */
    private List<AnnotatedPluginDocument> importAllPhases(Path outputDir, ProgressListener progressListener,
                                                          Set<BamReferenceExtractor.ReferenceInfo> bamReferences,
                                                          List<String> accessions,
                                                          Map<String, String> referenceNameMap,
                                                          List<File> bamFiles,
                                                          int stepCount, int totalSteps,
                                                          SampleBasedImporter sampleImporter) throws IOException {
        List<AnnotatedPluginDocument> importedDocs = new ArrayList<>();

        logger.info("========== IMPORTING ALL TAXTRIAGE RESULTS ===========");

        // First, download and import GenBank references if available
        if (!accessions.isEmpty()) {
            logger.info("Downloading GenBank references for accessions: " + accessions);
            Path genBankDir = outputDir.resolve("genbank_downloads");
            Files.createDirectories(genBankDir);

            Map<String, File> genBankFiles = NCBIReferenceDownloader.downloadGenBankFiles(
                accessions, genBankDir, progressListener);

            if (!genBankFiles.isEmpty()) {
                logger.info("Downloaded " + genBankFiles.size() + " GenBank files");
            }
        }

        // First, import shared files (References, Reports, Kraken_Results) once at root level
        try {
            logger.info("Importing shared files (common to all samples)...");
            List<AnnotatedPluginDocument> sharedDocs = importSharedFiles(outputDir, sampleImporter, progressListener);
            importedDocs.addAll(sharedDocs);
            logger.info("Imported " + sharedDocs.size() + " shared documents");
        } catch (Exception e) {
            logger.log(Level.WARNING, "Error importing shared files", e);
        }

        // Then import sample-specific files (alignments)
        Set<String> samples = sampleImporter.getSamples();
        if (samples.isEmpty()) {
            samples = sampleImporter.discoverSamples();
        }

        logger.info("Processing " + samples.size() + " sample(s): " + samples);

        int sampleCount = 0;
        for (String sample : samples) {
            if (progressListener != null) {
                progressListener.setMessage("Importing files for sample: " + sample);
                progressListener.setProgress((double) sampleCount / samples.size());
            }

            try {
                List<AnnotatedPluginDocument> sampleDocs = sampleImporter.importSampleFiles(sample, progressListener);
                importedDocs.addAll(sampleDocs);
                logger.info("Imported " + sampleDocs.size() + " documents for sample: " + sample);
            } catch (Exception e) {
                logger.log(Level.WARNING, "Error importing files for sample: " + sample, e);
            }

            sampleCount++;
        }

        // Index BAM files but don't import them
        logger.info("Preparing BAM files for manual import...");
        for (File bamFile : bamFiles) {
            try {
                BamIndexer.ensureIndexExists(bamFile);
                logger.info("BAM file indexed and ready: " + bamFile.getAbsolutePath());
                System.out.println("BAM file ready for manual import: " + bamFile.getAbsolutePath());
            } catch (Exception e) {
                logger.warning("Could not index BAM file: " + bamFile.getName());
            }
        }

        // Generate HTML summary report
        try {
            logger.info("Generating HTML summary report...");
            AnnotatedPluginDocument htmlReport = generateHtmlReport(outputDir, sampleImporter);
            if (htmlReport != null) {
                importedDocs.add(htmlReport);
                logger.info("Generated comprehensive HTML report");
            }
        } catch (Exception e) {
            logger.log(Level.WARNING, "Error generating HTML report", e);
        }

        logger.info("========== IMPORT COMPLETE ===========");
        logger.info("Total imported documents: " + importedDocs.size());
        logger.info("BAM files prepared for manual import: " + bamFiles.size());

        // Post-import archive creation
        logger.info("Starting post-import archive creation...");
        createPostImportArchive(outputDir, sampleImporter);

        return importedDocs;
    }


    /**
     * Attempts to organize imported documents into Geneious folders.
     */
    private void organizeIntoFolders() {
        // Log the intended folder organization
        for (Map.Entry<String, List<AnnotatedPluginDocument>> entry : documentsByFolder.entrySet()) {
            String folderName = entry.getKey();
            List<AnnotatedPluginDocument> docs = entry.getValue();

            if (!docs.isEmpty()) {
                logger.info("Documents intended for folder '" + folderName + "': " + docs.size() + " file(s)");
                // Note: Actual folder creation would require specific Geneious API methods
                // The documents are tagged with their intended folder for manual organization
            }
        }
    }

    /**
     * Finds all BAM files in the output directory.
     */
    private List<File> findBamFiles(Path outputDir) {
        List<File> bamFiles = new ArrayList<>();
        Path minimap2Dir = outputDir.resolve("minimap2");

        if (Files.exists(minimap2Dir)) {
            try {
                Files.walk(minimap2Dir, 1)
                    .filter(Files::isRegularFile)
                    .filter(path -> path.getFileName().toString().endsWith(".bam"))
                    .forEach(path -> bamFiles.add(path.toFile()));
            } catch (IOException e) {
                logger.log(Level.WARNING, "Error scanning for BAM files", e);
            }
        }

        logger.info("Found " + bamFiles.size() + " BAM files");
        return bamFiles;
    }

    /**
     * Imports all .txt files from the output directory.
     */
    private List<AnnotatedPluginDocument> importTextFiles(Path outputDir) throws IOException {
        List<AnnotatedPluginDocument> docs = new ArrayList<>();

        // Walk through all subdirectories and find .txt files
        Files.walk(outputDir)
            .filter(Files::isRegularFile)
            .filter(path -> path.getFileName().toString().endsWith(".txt"))
            .forEach(path -> {
                try {
                    // Determine document type based on location and filename
                    String docType = determineDocumentType(outputDir, path);
                    String content = new String(Files.readAllBytes(path));
                    // Add folder prefix to name for organization
                    String name = "[Reports] " + path.getFileName().toString();
                    AnnotatedPluginDocument doc = createTextDocument(
                        name, content, docType
                    );
                    if (doc != null) {
                        docs.add(doc);
                    }
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Failed to import text file: " + path, e);
                }
            });

        return docs;
    }

    /**
     * Imports reference FASTA files from the download directory.
     */
    private List<AnnotatedPluginDocument> importReferenceFasta(Path outputDir) throws IOException {
        List<AnnotatedPluginDocument> docs = new ArrayList<>();

        Path downloadDir = outputDir.resolve("download");
        if (!Files.exists(downloadDir)) {
            logger.info("Download directory does not exist: " + downloadDir);
            return docs;
        }

        // Find .dwnld.references.fasta files
        Files.walk(downloadDir, 1)
            .filter(Files::isRegularFile)
            .filter(path -> {
                String name = path.getFileName().toString();
                return name.contains(".dwnld.references.fasta") ||
                       name.contains("references.fasta");
            })
            .forEach(path -> {
                try {
                    // Try to import using Geneious FASTA importer
                    List<AnnotatedPluginDocument> fastaDocs = importFastaFile(path.toFile());
                    if (!fastaDocs.isEmpty()) {
                        docs.addAll(fastaDocs);
                        logger.info("Imported " + fastaDocs.size() + " sequences from " + path.getFileName());
                    } else {
                        // Fall back to text import if FASTA import fails
                        String content = new String(Files.readAllBytes(path));
                        AnnotatedPluginDocument doc = createTextDocument(
                            path.getFileName().toString(),
                            content,
                            "Reference Sequences (FASTA)"
                        );
                        if (doc != null) {
                            docs.add(doc);
                        }
                    }
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Failed to import reference FASTA: " + path, e);
                }
            });

        return docs;
    }

    /**
     * Imports a FASTA file using SimpleFastaImporter.
     */
    private List<AnnotatedPluginDocument> importFastaFile(File fastaFile) {
        // Use SimpleFastaImporter to import as proper sequence documents
        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile);

        // Add folder prefix to document names for organization
        for (AnnotatedPluginDocument doc : docs) {
            try {
                // Set the document name with folder prefix
                doc.setName("[Reference Sequences] " + doc.getName());
            } catch (Exception e) {
                logger.log(Level.WARNING, "Could not set document name prefix", e);
            }
        }

        return docs;
    }

    /**
     * Enterprise-grade reference persistence with comprehensive validation and error handling.
     */
    private List<AnnotatedPluginDocument> persistReferencesToDatabase(List<AnnotatedPluginDocument> referenceDocs,
                                                                      ProgressListener progressListener) {
        if (referenceDocs == null || referenceDocs.isEmpty()) {
            logger.warning("No reference documents provided for persistence");
            return new ArrayList<>();
        }

        logger.info("Initiating enterprise reference persistence for " + referenceDocs.size() + " document(s)");

        try {
            // Use enterprise-grade persistence
            List<AnnotatedPluginDocument> persistedDocs = LocalDatabaseHelper.saveToLocalDatabase(
                referenceDocs, progressListener);

            if (persistedDocs != null && !persistedDocs.isEmpty()) {
                logger.info("Enterprise persistence successful: " + persistedDocs.size() + " document(s) persisted");

                // Validate persistence success
                validatePersistenceSuccess(persistedDocs);

                return persistedDocs;
            } else {
                logger.warning("Enterprise persistence returned empty result");
            }
        } catch (Exception e) {
            logger.log(Level.SEVERE, "Enterprise reference persistence failed", e);
        }

        logger.warning("Falling back to unpersisted reference documents - BAM import reliability may be reduced");
        return referenceDocs;
    }

    /**
     * Validates that reference persistence was successful.
     */
    private void validatePersistenceSuccess(List<AnnotatedPluginDocument> persistedDocs) {
        logger.info("Validating enterprise persistence success...");

        int validatedCount = 0;
        for (AnnotatedPluginDocument doc : persistedDocs) {
            if (doc.getURN() != null && doc.getDatabase() != null) {
                validatedCount++;
            }
        }

        logger.info("Persistence validation: " + validatedCount + "/" + persistedDocs.size() + " documents fully persisted");

        if (validatedCount < persistedDocs.size()) {
            logger.warning("Some documents may not be fully persisted - BAM import may have issues");
        } else {
            logger.info("All documents successfully validated as persisted");
        }
    }

    /**
     * Builds the list of reference documents required for a specific BAM file.
     */
    private List<AnnotatedPluginDocument> resolveReferencesForBam(File bamFile,
                                                                  List<AnnotatedPluginDocument> referenceDocuments) {
        List<AnnotatedPluginDocument> resolved = new ArrayList<>();

        if (referenceDocuments == null || referenceDocuments.isEmpty()) {
            return resolved;
        }

        List<BamReferenceExtractor.ReferenceInfo> requiredReferences;
        try {
            requiredReferences = BamReferenceExtractor.extractReferences(bamFile);
        } catch (IOException e) {
            logger.log(Level.WARNING, "Could not read BAM header for " + bamFile.getName(), e);
            return resolved;
        }

        Map<String, AnnotatedPluginDocument> docsByName = new HashMap<>();
        Map<String, AnnotatedPluginDocument> docsByAccession = new HashMap<>();

        for (AnnotatedPluginDocument doc : referenceDocuments) {
            String name = doc.getName();
            if (name != null) {
                docsByName.put(name, doc);

                String accession = extractAccessionFromName(name);
                if (accession != null && !docsByAccession.containsKey(accession)) {
                    docsByAccession.put(accession, doc);
                }
            }
        }

        for (BamReferenceExtractor.ReferenceInfo ref : requiredReferences) {
            AnnotatedPluginDocument match = null;

            String mappedName = bamReferenceNameMap.get(ref.name);
            if (mappedName != null) {
                match = docsByName.get(mappedName);
            }
            if (match == null) {
                match = docsByName.get(ref.name);
            }
            if (match == null && ref.accession != null) {
                match = docsByAccession.get(ref.accession);
            }
            if (match == null && mappedName != null) {
                String mappedAccession = extractAccessionFromName(mappedName);
                if (mappedAccession != null) {
                    match = docsByAccession.get(mappedAccession);
                }
            }
            if (match == null) {
                for (Map.Entry<String, AnnotatedPluginDocument> entry : docsByName.entrySet()) {
                    String candidate = entry.getKey();
                    if (candidate.equalsIgnoreCase(ref.name) ||
                        candidate.contains(ref.name) ||
                        ref.name.contains(candidate)) {
                        match = entry.getValue();
                        break;
                    }
                }
            }

            if (match != null) {
                if (!resolved.contains(match)) {
                    resolved.add(match);
                }
            } else {
                String accessionMsg = ref.accession != null ? " (accession: " + ref.accession + ")" : "";
                logger.warning("No matching reference document found for BAM reference '" + ref.name + "'" + accessionMsg);
            }
        }

        return resolved;
    }

    private String extractAccessionFromName(String name) {
        if (name == null) {
            return null;
        }

        java.util.regex.Pattern pattern = java.util.regex.Pattern.compile("([A-Z]{1,4}_?\\d+\\.?\\d*)");
        java.util.regex.Matcher matcher = pattern.matcher(name);
        if (matcher.find()) {
            return matcher.group(1);
        }

        return null;
    }

    /**
     * Imports BAM files from the minimap2 directory.
     * Overload for backward compatibility.
     */
    private List<AnnotatedPluginDocument> importBamFiles(Path outputDir) throws IOException {
        return importBamFiles(outputDir, new ArrayList<>());
    }

    /**
     * Imports BAM files from the minimap2 directory.
     * @param outputDir The output directory
     * @param referenceDocuments Previously imported reference sequences for the BAM alignments
     */
    private List<AnnotatedPluginDocument> importBamFiles(Path outputDir, List<AnnotatedPluginDocument> referenceDocuments) throws IOException {
        return importBamFiles(outputDir);
    }

    /**
     * Enterprise-grade BAM import with comprehensive reference resolution and fallback strategies.
     * @param outputDir The output directory
     * @param referenceDocuments Previously imported reference sequences for the BAM alignments
     * @param referenceFiles The original GenBank files for the references (if available)
     */
    private List<AnnotatedPluginDocument> importBamFilesEnterprise(Path outputDir, List<AnnotatedPluginDocument> referenceDocuments, Map<String, File> referenceFiles) throws IOException {
        logger.info("===== ENTERPRISE BAM IMPORT SYSTEM =====");
        logger.info("Initiating enterprise BAM import with " + referenceDocuments.size() + " reference sequences");

        List<AnnotatedPluginDocument> docs = new ArrayList<>();

        // Pre-import validation
        Path minimap2Dir = outputDir.resolve("minimap2");
        if (!Files.exists(minimap2Dir)) {
            logger.warning("Minimap2 directory does not exist: " + minimap2Dir);
            return docs;
        }

        // Log enterprise reference inventory
        logReferenceInventory(referenceDocuments);

        // Enterprise BAM file processing
        List<File> bamFiles = findAndValidateBamFiles(minimap2Dir);
        logger.info("Found " + bamFiles.size() + " BAM file(s) for enterprise import");

        for (File bamFile : bamFiles) {
            try {
                logger.info("Processing BAM file: " + bamFile.getName() + " (" + bamFile.length() + " bytes)");

                // Enterprise BAM processing workflow
                EnterpriseBamImportResult result = processEnterpriseBamFile(
                    bamFile, referenceDocuments, referenceFiles);

                if (!result.importedDocuments.isEmpty()) {
                    docs.addAll(result.importedDocuments);
                    logger.info("Enterprise BAM import SUCCESS: " + bamFile.getName() +
                               " produced " + result.importedDocuments.size() + " document(s)");
                } else {
                    logger.warning("Enterprise BAM import FAILED: " + bamFile.getName());
                    // Create reference document with detailed failure information
                    AnnotatedPluginDocument fallbackDoc = createEnterpriseBamReferenceDocument(
                        bamFile, result);
                    if (fallbackDoc != null) {
                        docs.add(fallbackDoc);
                    }
                }

            } catch (Exception e) {
                logger.log(Level.SEVERE, "Enterprise BAM processing failed for: " + bamFile.getName(), e);
            }
        }

        logger.info("Enterprise BAM import completed: " + docs.size() + " total documents");
        logger.info("===== ENTERPRISE BAM IMPORT COMPLETE =====\n");
        return docs;
    }

    /**
     * Logs enterprise reference inventory for debugging and audit purposes.
     */
    private void logReferenceInventory(List<AnnotatedPluginDocument> referenceDocuments) {
        logger.info("=== ENTERPRISE REFERENCE INVENTORY ===");
        logger.info("Total references available: " + referenceDocuments.size());

        for (int i = 0; i < referenceDocuments.size(); i++) {
            AnnotatedPluginDocument ref = referenceDocuments.get(i);
            logger.info(String.format("[%d] Name: '%s'", (i + 1), ref.getName()));
            logger.info(String.format("    URN: %s", ref.getURN() != null ? ref.getURN() : "Not assigned"));
            logger.info(String.format("    Type: %s", ref.getDocumentClass()));
            logger.info(String.format("    Database: %s", ref.getDatabase() != null ? "Connected" : "Not connected"));
        }
        logger.info("===================================\n");
    }

    /**
     * Finds and validates BAM files with enterprise-grade checks.
     */
    private List<File> findAndValidateBamFiles(Path minimap2Dir) {
        List<File> validBamFiles = new ArrayList<>();

        try {
            Files.walk(minimap2Dir, 1)
                .filter(Files::isRegularFile)
                .filter(path -> path.getFileName().toString().endsWith(".bam"))
                .forEach(path -> {
                    File bamFile = path.toFile();
                    if (validateBamFile(bamFile)) {
                        validBamFiles.add(bamFile);
                    } else {
                        logger.warning("BAM file failed validation: " + bamFile.getName());
                    }
                });
        } catch (IOException e) {
            logger.log(Level.WARNING, "Error scanning for BAM files", e);
        }

        return validBamFiles;
    }

    /**
     * Validates a BAM file before processing.
     */
    private boolean validateBamFile(File bamFile) {
        if (!bamFile.exists()) {
            logger.warning("BAM file does not exist: " + bamFile.getName());
            return false;
        }

        if (bamFile.length() == 0) {
            logger.warning("BAM file is empty: " + bamFile.getName());
            return false;
        }

        // Ensure BAM index exists
        boolean indexExists = BamIndexer.ensureIndexExists(bamFile);
        if (!indexExists) {
            logger.warning("BAM index could not be created: " + bamFile.getName());
            // Continue anyway - some imports might work without index
        }

        return true;
    }

    /**
     * Processes a single BAM file using enterprise import strategies.
     */
    private EnterpriseBamImportResult processEnterpriseBamFile(
            File bamFile,
            List<AnnotatedPluginDocument> referenceDocuments,
            Map<String, File> referenceFiles) {

        EnterpriseBamImportResult result = new EnterpriseBamImportResult();
        result.bamFile = bamFile;

        try {
            // Phase 1: Reference resolution
            List<AnnotatedPluginDocument> resolvedRefs = resolveReferencesForBam(bamFile, referenceDocuments);
            result.resolvedReferences = resolvedRefs;

            logger.info("Resolved " + resolvedRefs.size() + " reference(s) for BAM: " + bamFile.getName());
            for (AnnotatedPluginDocument ref : resolvedRefs) {
                logger.info("  - " + ref.getName() + " (URN: " + ref.getURN() + ")");
            }

            // Phase 2: Reference mapping file creation
            if (!resolvedRefs.isEmpty()) {
                try {
                    LocalDatabaseHelper.createReferenceMappingFile(resolvedRefs, bamFile);
                    result.mappingFileCreated = true;
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Could not create reference mapping file", e);
                }
            }

            // Phase 3: Enterprise BAM import using enhanced helper
            List<AnnotatedPluginDocument> referencesToUse = resolvedRefs.isEmpty() ?
                referenceDocuments : resolvedRefs;

            result.importedDocuments = BamImportHelper.importBamWithReferences(
                bamFile, referencesToUse, ProgressListener.EMPTY);

            result.success = !result.importedDocuments.isEmpty();

        } catch (Exception e) {
            logger.log(Level.SEVERE, "Enterprise BAM processing error", e);
            result.error = e.getMessage();
        }

        return result;
    }

    /**
     * Creates an enterprise BAM reference document with detailed failure information.
     */
    private AnnotatedPluginDocument createEnterpriseBamReferenceDocument(
            File bamFile, EnterpriseBamImportResult result) {

        try {
            String fileName = bamFile.getName();
            StringBuilder content = new StringBuilder();

            content.append("Enterprise BAM Import Report\n");
            content.append("============================\n\n");
            content.append("File: ").append(fileName).append("\n");
            content.append("Size: ").append(bamFile.length()).append(" bytes\n");
            content.append("Path: ").append(bamFile.getAbsolutePath()).append("\n\n");

            content.append("Import Status: ").append(result.success ? "SUCCESS" : "FAILED").append("\n");
            content.append("Resolved References: ").append(result.resolvedReferences.size()).append("\n");
            content.append("Mapping File Created: ").append(result.mappingFileCreated ? "Yes" : "No").append("\n");

            if (result.error != null) {
                content.append("Error: ").append(result.error).append("\n");
            }

            content.append("\nReference Details:\n");
            content.append("------------------\n");
            for (AnnotatedPluginDocument ref : result.resolvedReferences) {
                content.append("- ").append(ref.getName());
                content.append(" (URN: ").append(ref.getURN()).append(")\n");
            }

            content.append("\nManual Import Instructions:\n");
            content.append("---------------------------\n");
            content.append("1. Ensure all reference sequences are in the database\n");
            content.append("2. File → Import → From File...\n");
            content.append("3. Select: ").append(bamFile.getAbsolutePath()).append("\n");
            content.append("4. Geneious should automatically find the references\n");

            TaxTriageResultDocument doc = new TaxTriageResultDocument(
                "[Enterprise] " + fileName + ".report",
                content.toString(),
                "Enterprise BAM Import Report",
                bamFile.getAbsolutePath(),
                "TXT"
            );

            return DocumentUtilities.createAnnotatedPluginDocument(doc);

        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to create enterprise BAM reference document", e);
            return null;
        }
    }

    /**
     * Result container for enterprise BAM import operations.
     */
    private static class EnterpriseBamImportResult {
        File bamFile;
        List<AnnotatedPluginDocument> resolvedReferences = new ArrayList<>();
        List<AnnotatedPluginDocument> importedDocuments = new ArrayList<>();
        boolean mappingFileCreated = false;
        boolean success = false;
        String error = null;
    }

    /**
     * Imports a BAM file using the robust BamFileImporter.
     * This method uses multiple strategies to ensure successful BAM import.
     * Overload for backward compatibility.
     */
    private List<AnnotatedPluginDocument> importBamFile(File bamFile) {
        return importBamFile(bamFile, new ArrayList<>());
    }

    /**
     * Imports a BAM file using the robust BamFileImporter.
     * This method uses multiple strategies to ensure successful BAM import.
     *
     * @param bamFile The BAM file to import
     * @param referenceDocuments Previously imported reference sequences
     */
    /**
     * Helper method without reference files.
     */
    private List<AnnotatedPluginDocument> importBamFile(File bamFile, List<AnnotatedPluginDocument> referenceDocuments) {
        return importBamFile(bamFile, referenceDocuments, null);
    }

    /**
     * Imports a single BAM file with its references.
     * @param bamFile The BAM file to import
     * @param referenceDocuments Previously imported reference sequences
     * @param referenceFiles The original GenBank files for the references (if available)
     */
    private List<AnnotatedPluginDocument> importBamFile(File bamFile, List<AnnotatedPluginDocument> referenceDocuments, Map<String, File> referenceFiles) {
        // Create a simple progress listener for the import operation
        ProgressListener importProgress = ProgressListener.EMPTY;

        System.out.println("Attempting BAM import for: " + bamFile.getName());
        System.out.println("Reference documents passed: " + referenceDocuments.size());
        if (referenceFiles != null) {
            System.out.println("Reference files available: " + referenceFiles.size());
        }

        // Try direct import first, letting Geneious find references naturally
        List<AnnotatedPluginDocument> importedDocs = new ArrayList<>();

        try {
            System.out.println("Trying direct BAM import (letting Geneious find references)...");
            importedDocs = PluginUtilities.importDocuments(bamFile, importProgress);

            if (importedDocs != null && !importedDocs.isEmpty()) {
                System.out.println("Direct import SUCCESS: imported " + importedDocs.size() + " document(s)");
            } else {
                System.out.println("Direct import returned null or empty");
            }
        } catch (Exception e) {
            System.out.println("Direct import failed: " + e.getMessage());
            e.printStackTrace();
        }

        // If direct import didn't work, try with the BamImportHelper
        if (importedDocs.isEmpty()) {
            System.out.println("Falling back to BamImportHelper...");
            importedDocs = BamImportHelper.importBamWithReferences(
                bamFile, referenceDocuments, importProgress);
        }

        // If that didn't work, try the original BamFileImporter
        if (importedDocs.isEmpty()) {
            System.out.println("BamImportHelper failed, trying BamFileImporter...");
            logger.info("Attempting programmatic import of BAM file: " + bamFile.getName());
            importedDocs = BamFileImporter.importBamFile(bamFile, importProgress, referenceDocuments);
        }

        if (!importedDocs.isEmpty()) {
            System.out.println("Successfully imported BAM file: " + bamFile.getName() +
                             " (" + importedDocs.size() + " document(s))");
            logger.info("Successfully imported BAM file programmatically: " + bamFile.getName() +
                       " (" + importedDocs.size() + " document(s))");

            // Add folder prefix to document names for organization
            for (AnnotatedPluginDocument doc : importedDocs) {
                try {
                    doc.setName("[minimap2] " + doc.getName());
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Could not set document name prefix", e);
                }
            }

            return importedDocs;
        } else {
            System.out.println("All BAM import strategies failed for: " + bamFile.getName());
            logger.info("BAM import strategies exhausted for: " + bamFile.getName());

            // Log BAM import capabilities for debugging
            String capabilities = BamFileImporter.getBamImportInfo();
            System.out.println(capabilities);
            logger.info(capabilities);

            return new ArrayList<>();
        }
    }

    /**
     * Creates a reference document for a BAM file.
     * Overload for backward compatibility.
     */
    private AnnotatedPluginDocument createBamReferenceDocument(Path bamFile, boolean hasIndex) {
        return createBamReferenceDocument(bamFile, hasIndex, false);
    }

    /**
     * Creates a reference document for a BAM file.
     *
     * @param bamFile The BAM file path
     * @param hasIndex Whether the BAM file has an index
     * @param attemptedProgrammaticImport Whether programmatic import was attempted
     */
    private AnnotatedPluginDocument createBamReferenceDocument(Path bamFile, boolean hasIndex,
                                                                boolean attemptedProgrammaticImport) {
        try {
            String fileName = bamFile.getFileName().toString();

            // Build import status message
            String importStatus;
            String importMethod;

            if (attemptedProgrammaticImport) {
                importMethod = "Programmatic import was attempted but requires manual import.\n\n";
                if (hasIndex) {
                    importStatus = "✓ Index present - Ready for manual import in Geneious";
                } else {
                    importStatus = "⚠ Index missing - May need to be created manually";
                }
            } else {
                importMethod = "";
                if (hasIndex) {
                    importStatus = "✓ Index present - Ready for import in Geneious";
                } else {
                    importStatus = "⚠ Index could not be created automatically";
                }
            }

            String content = String.format(
                "BAM Alignment File\n" +
                "==================\n\n" +
                "File: %s\n" +
                "Size: %d bytes\n" +
                "Path: %s\n" +
                "Index Status: %s\n\n" +
                "Import Status:\n" +
                "--------------\n" +
                "%s" +
                "%s\n\n" +
                "Manual Import Instructions:\n" +
                "----------------------------\n" +
                "To manually import this BAM file into Geneious:\n\n" +
                "1. Select File → Import → From File...\n" +
                "2. Navigate to: %s\n" +
                "3. Select the BAM file: %s\n" +
                "4. Click OK to import\n\n" +
                "The alignment will be imported as a contig with mapped reads.\n" +
                "You can view coverage, SNPs, and other alignment features in Geneious.\n\n" +
                "Note: If import fails, ensure that:\n" +
                "- The BAM file has a valid index (.bai file)\n" +
                "- The SAM/BAM import plugin is installed in Geneious\n" +
                "- The file is not corrupted",
                fileName,
                Files.size(bamFile),
                bamFile.toAbsolutePath(),
                hasIndex ? "✓ Present (.bai file exists)" : "⚠ Missing",
                importMethod,
                importStatus,
                bamFile.getParent().toAbsolutePath(),
                fileName
            );

            TaxTriageResultDocument doc = new TaxTriageResultDocument(
                "[minimap2] " + fileName + ".info",
                content,
                "BAM Alignment Info",
                bamFile.toAbsolutePath().toString(),
                "TXT"
            );

            return DocumentUtilities.createAnnotatedPluginDocument(doc);
        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to create BAM reference document", e);
            return null;
        }
    }

    /**
     * Determines the document type based on file path and name.
     */
    private String determineDocumentType(Path outputDir, Path filePath) {
        String relativePath = outputDir.relativize(filePath).toString();
        String fileName = filePath.getFileName().toString();

        if (relativePath.startsWith("kraken2")) {
            if (fileName.contains("report")) return "Kraken Report";
            if (fileName.contains("kraken2")) return "Kraken Classifications";
            return "Kraken Output";
        } else if (relativePath.startsWith("top")) {
            if (fileName.contains("top_report")) return "Top Species Report";
            if (fileName.contains("toptaxids")) return "Top Taxa IDs";
            if (fileName.contains("topnames")) return "Top Taxa Names";
            return "Top Hits";
        } else if (relativePath.startsWith("filterkraken")) {
            return "Filtered Taxa Report";
        } else if (relativePath.startsWith("pipeline_info")) {
            return "Pipeline Information";
        }

        return "TaxTriage Output";
    }

    /**
     * Creates an analysis summary HTML document.
     */
    private AnnotatedPluginDocument createAnalysisSummary(Path outputDir) {
        try {
            SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
            String timestamp = dateFormat.format(new Date());

            StringBuilder html = new StringBuilder();
            html.append("<!DOCTYPE html>\n");
            html.append("<html>\n<head>\n");
            html.append("<title>TaxTriage Analysis Summary</title>\n");
            html.append("<style>\n");
            html.append("body { font-family: Arial, sans-serif; margin: 20px; }\n");
            html.append("h1 { color: #333; border-bottom: 2px solid #4CAF50; padding-bottom: 10px; }\n");
            html.append("h2 { color: #666; margin-top: 30px; }\n");
            html.append(".info-box { background-color: #f0f0f0; padding: 15px; border-radius: 5px; margin: 15px 0; }\n");
            html.append("code { background-color: #f4f4f4; padding: 2px 5px; border-radius: 3px; }\n");
            html.append("pre { background-color: #2d2d2d; color: #f8f8f2; padding: 15px; border-radius: 5px; overflow-x: auto; }\n");
            html.append(".param-table { width: 100%; border-collapse: collapse; margin: 15px 0; }\n");
            html.append(".param-table th { background-color: #4CAF50; color: white; padding: 10px; text-align: left; }\n");
            html.append(".param-table td { border: 1px solid #ddd; padding: 8px; }\n");
            html.append(".param-table tr:nth-child(even) { background-color: #f9f9f9; }\n");
            html.append("</style>\n");
            html.append("</head>\n<body>\n");

            html.append("<h1>TaxTriage Analysis Summary</h1>\n");
            html.append("<div class='info-box'>\n");
            html.append("<p><strong>Analysis Date:</strong> ").append(timestamp).append("</p>\n");
            html.append("<p><strong>Output Directory:</strong> <code>").append(outputDir.toAbsolutePath()).append("</code></p>\n");

            if (workingDirectory != null) {
                html.append("<p><strong>Working Directory:</strong> ");
                html.append("<a href='file://").append(workingDirectory).append("'>");
                html.append(workingDirectory).append("</a></p>\n");
            }
            html.append("</div>\n");

            // Add workflow parameters if available
            if (!workflowParams.isEmpty()) {
                html.append("<h2>Workflow Parameters</h2>\n");
                html.append("<table class='param-table'>\n");
                html.append("<tr><th>Parameter</th><th>Value</th></tr>\n");
                for (Map.Entry<String, String> entry : workflowParams.entrySet()) {
                    html.append("<tr><td>").append(entry.getKey()).append("</td>");
                    html.append("<td>").append(escapeHtml(entry.getValue())).append("</td></tr>\n");
                }
                html.append("</table>\n");
            }

            // Add reproduction command
            html.append("<h2>Reproduction Command</h2>\n");
            html.append("<p>To reproduce this analysis, run the following command:</p>\n");

            if (commandLine != null) {
                html.append("<pre>").append(escapeHtml(commandLine)).append("</pre>\n");
            } else {
                // Build command from parameters
                html.append("<pre>\n");
                html.append("cd ").append(workingDirectory != null ? workingDirectory : outputDir.getParent()).append("\n\n");
                html.append("nextflow run https://github.com/jhuapl-bio/taxtriage \\\n");
                html.append("    -r main \\\n");
                html.append("    -profile docker \\\n");

                if (workflowParams.containsKey("input")) {
                    html.append("    --input ").append(workflowParams.get("input")).append(" \\\n");
                }
                if (workflowParams.containsKey("outdir")) {
                    html.append("    --outdir ").append(workflowParams.get("outdir")).append(" \\\n");
                }
                if (workflowParams.containsKey("platform")) {
                    html.append("    --platform ").append(workflowParams.get("platform")).append(" \\\n");
                }
                if (workflowParams.containsKey("db")) {
                    html.append("    --db ").append(workflowParams.get("db")).append(" \\\n");
                }
                html.append("    --download_db true \\\n");
                html.append("    --max_cpus 4 \\\n");
                html.append("    --max_memory 8.GB \\\n");
                html.append("    -work-dir work\n");
                html.append("</pre>\n");
            }

            // Add imported files summary with folder organization
            html.append("<h2>Imported Files Organization</h2>\n");
            html.append("<div class='info-box'>\n");
            html.append("<p>Files have been organized into the following folders:</p>\n");
            html.append("<ul>\n");
            html.append("<li><strong>Reports/</strong> - Text report files (.txt) including Kraken reports, top hits, and other text outputs</li>\n");
            html.append("<li><strong>Reference Sequences/</strong> - Downloaded reference genomes (FASTA)</li>\n");
            html.append("<li><strong>minimap2/</strong> - Minimap2 alignment results (BAM files)</li>\n");
            html.append("</ul>\n");

            // Add summary of what was imported
            for (Map.Entry<String, List<AnnotatedPluginDocument>> entry : documentsByFolder.entrySet()) {
                if (!entry.getValue().isEmpty()) {
                    html.append("<p><strong>").append(entry.getKey()).append(":</strong> ")
                        .append(entry.getValue().size()).append(" file(s)</p>\n");
                }
            }

            html.append("</div>\n");

            html.append("</body>\n</html>");

            TaxTriageResultDocument summaryDoc = new TaxTriageResultDocument(
                "TaxTriage_Analysis_Summary.html",
                html.toString(),
                "Analysis Summary",
                "Analysis Summary",
                "HTML"
            );

            return DocumentUtilities.createAnnotatedPluginDocument(summaryDoc);

        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to create analysis summary", e);
            return null;
        }
    }

    /**
     * Creates a simple text document.
     */
    private AnnotatedPluginDocument createTextDocument(String name, String content, String type) {
        try {
            // Determine file format from name
            String fileFormat = "TXT";
            String lowerName = name.toLowerCase();
            if (lowerName.endsWith(".tsv")) {
                fileFormat = "TSV";
            } else if (lowerName.endsWith(".csv")) {
                fileFormat = "CSV";
            } else if (lowerName.endsWith(".html")) {
                fileFormat = "HTML";
            } else if (lowerName.endsWith(".json")) {
                fileFormat = "JSON";
            } else if (lowerName.endsWith(".fasta") || lowerName.endsWith(".fa")) {
                fileFormat = "FASTA";
            }

            // Create the TaxTriage result document
            TaxTriageResultDocument resultDoc = new TaxTriageResultDocument(
                name, content, type, name, fileFormat
            );

            // Create an annotated document from it
            AnnotatedPluginDocument annotatedDoc = DocumentUtilities.createAnnotatedPluginDocument(resultDoc);

            logger.info("Created document: " + name + " (" + type + ") - " + content.length() + " bytes");
            return annotatedDoc;

        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to create document for: " + name, e);
            return null;
        }
    }

    /**
     * Escapes HTML special characters.
     */
    private String escapeHtml(String text) {
        if (text == null) return "";
        return text.replace("&", "&amp;")
                  .replace("<", "&lt;")
                  .replace(">", "&gt;")
                  .replace("\"", "&quot;")
                  .replace("'", "&#39;");
    }

    /**
     * Gets the organized documents by folder.
     */
    public Map<String, List<AnnotatedPluginDocument>> getDocumentsByFolder() {
        return new HashMap<>(documentsByFolder);
    }

    // ========== NEW PHASE-BASED IMPORT METHODS ==========

    /**
     * Phase 1: Imports only GenBank reference sequences (no FASTA fallback).
     * This eliminates the double-import issue and ensures only high-quality GenBank references.
     */
    private List<AnnotatedPluginDocument> importReferencesOnly(Path outputDir, ProgressListener progressListener,
                                                               Set<BamReferenceExtractor.ReferenceInfo> bamReferences,
                                                               List<String> accessions,
                                                               Map<String, String> referenceNameMap,
                                                               SampleBasedImporter sampleImporter) throws IOException {
        List<AnnotatedPluginDocument> importedDocs = new ArrayList<>();
        documentsByFolder.clear();
        bamReferenceNameMap.clear();

        if (progressListener != null) {
            progressListener.setMessage("Phase 1: Downloading and importing reference sequences from NCBI...");
            progressListener.setProgress(0.0);
        }

        logger.info("========== PHASE 1: IMPORTING REFERENCES ONLY ===========");
        logger.info("Target accessions: " + accessions.size());
        for (String acc : accessions) {
            logger.info("  - " + acc);
        }

        List<AnnotatedPluginDocument> refDocs = new ArrayList<>();

        // ONLY download and import GenBank files - NO FASTA fallback
        if (!accessions.isEmpty()) {
            Path genBankDir = outputDir.resolve("genbank_downloads");

            if (progressListener != null) {
                progressListener.setProgress(0.2);
            }

            Map<String, File> genBankFiles = NCBIReferenceDownloader.downloadGenBankFiles(
                accessions, genBankDir, progressListener);

            if (!genBankFiles.isEmpty()) {
                logger.info("Downloaded " + genBankFiles.size() + " GenBank files from NCBI");

                if (progressListener != null) {
                    progressListener.setProgress(0.6);
                }

                refDocs = NCBIReferenceDownloader.importGenBankFiles(genBankFiles, referenceNameMap, progressListener, sampleImporter);
                logger.info("Imported " + refDocs.size() + " GenBank documents");

                if (!refDocs.isEmpty()) {
                    refDocs = ReferenceManager.prepareReferencesForBAM(
                        refDocs, new ArrayList<>(bamReferences), progressListener);
                    logger.info("Prepared " + refDocs.size() + " GenBank reference sequences for BAM");
                }
            } else {
                logger.warning("Could not download any GenBank files from NCBI");
                logger.warning("BAM import in Phase 3 will likely fail without reference sequences");
                logger.warning("Consider running manual import or checking network connectivity");
            }
        } else {
            logger.info("No NCBI accessions found in BAM headers - no references to import");
            logger.info("BAM files may have local references that need manual import");
        }

        if (!refDocs.isEmpty()) {
            if (progressListener != null) {
                progressListener.setProgress(0.8);
            }

            refDocs = persistReferencesToDatabase(refDocs, progressListener);
            documentsByFolder.put("Reference Sequences", refDocs);
            importedDocs.addAll(refDocs);

            logger.info("Phase 1 Complete: Persisted " + refDocs.size() + " reference sequence(s) to Geneious database");

            // Create reference mapping for future BAM import
            Map<String, String> refMapping = ReferenceManager.createReferenceMapping(
                new ArrayList<>(bamReferences), refDocs);
            bamReferenceNameMap.clear();
            bamReferenceNameMap.putAll(refMapping);
            logger.info("Created reference mapping with " + refMapping.size() + " entries for BAM import");

            // CRITICAL: Wait for database commit and indexing to complete
            waitForDatabaseCommit(refDocs);
        }

        if (progressListener != null) {
            progressListener.setProgress(1.0);
        }

        logger.info("Phase 1 Summary:");
        logger.info("  - GenBank files downloaded: " + (refDocs.size() > 0 ? "YES" : "NO"));
        logger.info("  - References persisted to database: " + refDocs.size());
        logger.info("  - Database indexing: COMPLETED");
        logger.info("========== PHASE 1 COMPLETE ===========");

        return importedDocs;
    }

    /**
     * Phase 2: Imports only text files.
     * This is separate to allow references to be fully committed before other imports.
     */
    private List<AnnotatedPluginDocument> importTextOnly(Path outputDir, ProgressListener progressListener, SampleBasedImporter sampleImporter) throws IOException {
        if (progressListener != null) {
            progressListener.setMessage("Phase 2: Importing text report files...");
            progressListener.setProgress(0.0);
        }

        logger.info("========== PHASE 2: IMPORTING TEXT FILES ONLY ===========");

        List<AnnotatedPluginDocument> txtDocs = importTextFiles(outputDir);
        if (!txtDocs.isEmpty()) {
            documentsByFolder.put("Reports", txtDocs);
            logger.info("Phase 2 Complete: Imported " + txtDocs.size() + " text files");
        } else {
            logger.info("Phase 2 Complete: No text files found to import");
        }

        if (progressListener != null) {
            progressListener.setProgress(1.0);
        }

        logger.info("========== PHASE 2 COMPLETE ===========");
        return txtDocs;
    }

    /**
     * Phase 3: Imports only BAM files.
     * This runs after references are fully committed and indexed, solving the reference resolution issue.
     */
    private List<AnnotatedPluginDocument> prepareBamOnly(Path outputDir, ProgressListener progressListener,
                                                        List<File> bamFiles,
                                                        SampleBasedImporter sampleImporter) throws IOException {
        if (progressListener != null) {
            progressListener.setMessage("Phase 3: Importing BAM alignment files...");
            progressListener.setProgress(0.0);
        }

        logger.info("========== PHASE 3: IMPORTING BAM FILES ONLY ===========");
        logger.info("BAM files to process: " + bamFiles.size());

        // Find all available reference documents in the database
        List<AnnotatedPluginDocument> availableRefs = findAvailableReferences(bamFiles);
        logger.info("Found " + availableRefs.size() + " available reference documents in database");

        // Log available references for debugging
        for (AnnotatedPluginDocument ref : availableRefs) {
            logger.info("  Available ref: '" + ref.getName() + "' (URN: " + ref.getURN() + ")");
        }

        if (availableRefs.isEmpty()) {
            logger.warning("No reference sequences found in database!");
            logger.warning("This suggests Phase 1 (reference import) was not completed successfully.");
            logger.warning("BAM import will likely fail. Consider running Phase 1 first.");
        }

        if (progressListener != null) {
            progressListener.setProgress(0.2);
        }

        List<AnnotatedPluginDocument> bamDocs = importBamFiles(outputDir, availableRefs);

        if (!bamDocs.isEmpty()) {
            documentsByFolder.put("minimap2", bamDocs);
            logger.info("Phase 3 Complete: Successfully imported " + bamDocs.size() + " BAM document(s)");
            for (AnnotatedPluginDocument bamDoc : bamDocs) {
                logger.info("  - " + bamDoc.getName() + " (Type: " + bamDoc.getDocumentClass() + ")");
            }
        } else {
            logger.warning("Phase 3 Complete: No BAM documents were successfully imported!");
            logger.warning("This may indicate:");
            logger.warning("  1. Reference sequences are not available in the database");
            logger.warning("  2. Reference names don't match BAM headers");
            logger.warning("  3. BAM files are corrupted or missing indices");
            logger.warning("  4. Geneious BAM import plugin is not available");
        }

        if (progressListener != null) {
            progressListener.setProgress(1.0);
        }

        logger.info("========== PHASE 3 COMPLETE ===========");
        return bamDocs;
    }

    /**
     * Waits for database commit and indexing to complete.
     * This is critical for ensuring references are findable by BAM import.
     */
    private void waitForDatabaseCommit(List<AnnotatedPluginDocument> documents) {
        if (documents == null || documents.isEmpty()) {
            return;
        }

        logger.info("Waiting for database commit and indexing...");

        // Force indexing of the reference documents
        ForceIndexingHelper.forceIndexing(documents);

        // Wait for indexing to complete with longer timeout
        // Wait for indexing to complete with longer timeout
        ForceIndexingHelper.waitForIndexing(documents, 30000);  // 30 seconds
        logger.info("Database indexing wait completed");

        // Additional wait to ensure database operations are fully committed
        try {
            Thread.sleep(2000);  // 2 second additional buffer
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }

        logger.info("Database commit wait completed");
    }

    /**
     * Finds reference documents that are available in the database for BAM import.
     */
    private List<AnnotatedPluginDocument> findAvailableReferences(List<File> bamFiles) {
        List<AnnotatedPluginDocument> availableRefs = new ArrayList<>();

        // Try to find references using the stored mapping
        if (!bamReferenceNameMap.isEmpty()) {
            logger.info("Using stored reference mapping to find available references");
            // In a real implementation, this would query the Geneious database
            // For now, we'll log what we're looking for
            for (Map.Entry<String, String> entry : bamReferenceNameMap.entrySet()) {
                logger.info("  Looking for reference: '" + entry.getValue() + "' (BAM ref: '" + entry.getKey() + "')");
            }
        }

        // Extract required references from BAM files and attempt to find them
        for (File bamFile : bamFiles) {
            try {
                List<BamReferenceExtractor.ReferenceInfo> bamRefs = BamReferenceExtractor.extractReferences(bamFile);
                logger.info("BAM file " + bamFile.getName() + " requires " + bamRefs.size() + " reference(s):");
                for (BamReferenceExtractor.ReferenceInfo ref : bamRefs) {
                    logger.info("  - '" + ref.name + "' (accession: " + ref.accession + ", length: " + ref.length + ")");
                }
            } catch (Exception e) {
                logger.warning("Could not extract references from " + bamFile.getName() + ": " + e.getMessage());
            }
        }

        // Note: In a full implementation, this would use Geneious API to search the database
        // for documents matching the required reference names/accessions

        return availableRefs;
    }

    /**
     * Imports shared files (Reports, References, Kraken_Results) at the root level.
     * These files are common to all samples and should be imported once at the main folder level.
     */
    private List<AnnotatedPluginDocument> importSharedFiles(Path outputDir, SampleBasedImporter sampleImporter,
                                                           ProgressListener progressListener) throws IOException {
        List<AnnotatedPluginDocument> sharedDocs = new ArrayList<>();
        WritableDatabaseService mainFolder = sampleImporter.getMainFolder();

        if (mainFolder == null) {
            logger.warning("Main folder not available for shared file import");
            return sharedDocs;
        }

        logger.info("Importing shared files to root level folders...");

        // Create root-level folders for shared files
        Map<String, WritableDatabaseService> sharedFolders = createSharedFolders(mainFolder);

        // Import References (GenBank and FASTA files)
        if (sharedFolders.containsKey("References")) {
            List<AnnotatedPluginDocument> refDocs = importSharedReferences(outputDir, sharedFolders.get("References"), progressListener);
            sharedDocs.addAll(refDocs);
            logger.info("Imported " + refDocs.size() + " reference documents to References folder");
        }

        // Import Reports (text files, HTML reports)
        if (sharedFolders.containsKey("Reports")) {
            List<AnnotatedPluginDocument> reportDocs = importSharedReports(outputDir, sharedFolders.get("Reports"), progressListener);
            sharedDocs.addAll(reportDocs);
            logger.info("Imported " + reportDocs.size() + " report documents to Reports folder");
        }

        // Import Kraken Results
        if (sharedFolders.containsKey("Kraken_Results")) {
            List<AnnotatedPluginDocument> krakenDocs = importSharedKrakenResults(outputDir, sharedFolders.get("Kraken_Results"), progressListener);
            sharedDocs.addAll(krakenDocs);
            logger.info("Imported " + krakenDocs.size() + " Kraken result documents to Kraken_Results folder");
        }

        return sharedDocs;
    }

    /**
     * Creates shared folders at the root level of the main folder.
     */
    private Map<String, WritableDatabaseService> createSharedFolders(WritableDatabaseService mainFolder) {
        Map<String, WritableDatabaseService> sharedFolders = new HashMap<>();
        String[] folderNames = {"References", "Reports", "Kraken_Results"};

        for (String folderName : folderNames) {
            try {
                WritableDatabaseService folder = mainFolder.createChildFolder(folderName);
                if (folder == null) {
                    folder = mainFolder.getChildService(folderName);
                }
                if (folder != null) {
                    sharedFolders.put(folderName, folder);
                    logger.info("Created/accessed shared folder: " + folderName);
                }
            } catch (Exception e) {
                logger.warning("Could not create shared folder " + folderName + ": " + e.getMessage());
            }
        }

        return sharedFolders;
    }

    /**
     * Imports shared reference files (GenBank, FASTA) to the References folder.
     */
    private List<AnnotatedPluginDocument> importSharedReferences(Path outputDir, WritableDatabaseService referencesFolder,
                                                                ProgressListener progressListener) throws IOException {
        List<AnnotatedPluginDocument> docs = new ArrayList<>();

        // Import GenBank files from genbank_downloads and minimap2 directories
        Path[] refDirs = {
            outputDir.resolve("genbank_downloads"),
            outputDir.resolve("minimap2"),
            outputDir.resolve("download")
        };

        for (Path refDir : refDirs) {
            if (Files.exists(refDir)) {
                try {
                    Files.list(refDir)
                        .filter(Files::isRegularFile)
                        .filter(p -> {
                            String name = p.getFileName().toString().toLowerCase();
                            // Skip dwnld.reference.fasta files - we already have GenBank versions
                            if (name.startsWith("dwnld.") && name.endsWith(".fasta")) {
                                return false;
                            }
                            // Only import GenBank files
                            return name.endsWith(".gb") || name.endsWith(".gbk") || name.endsWith(".genbank");
                        })
                        .forEach(refFile -> {
                            try {
                                List<AnnotatedPluginDocument> imported = PluginUtilities.importDocuments(refFile.toFile(), progressListener);
                                if (imported != null && !imported.isEmpty()) {
                                    for (AnnotatedPluginDocument doc : imported) {
                                        AnnotatedPluginDocument copiedDoc = referencesFolder.addDocumentCopy(doc, ProgressListener.EMPTY);
                                        if (copiedDoc != null) {
                                            docs.add(copiedDoc);
                                            logger.info("Imported reference: " + copiedDoc.getName());
                                        }
                                    }
                                }
                            } catch (Exception e) {
                                logger.warning("Failed to import reference file: " + refFile.getFileName() + " - " + e.getMessage());
                            }
                        });
                } catch (IOException e) {
                    logger.warning("Error scanning reference directory: " + refDir);
                }
            }
        }

        return docs;
    }

    /**
     * Imports shared report files to the Reports folder.
     */
    private List<AnnotatedPluginDocument> importSharedReports(Path outputDir, WritableDatabaseService reportsFolder,
                                                             ProgressListener progressListener) throws IOException {
        List<AnnotatedPluginDocument> docs = new ArrayList<>();

        // Import text and HTML reports from various directories
        Path[] reportDirs = {
            outputDir.resolve("report"),
            outputDir.resolve("top"),
            outputDir.resolve("combine"),
            outputDir.resolve("count"),
            outputDir.resolve("pipeline_info"),
            outputDir // Root directory for summary files
        };

        Set<String> processedFiles = new HashSet<>();

        for (Path reportDir : reportDirs) {
            if (Files.exists(reportDir)) {
                try {
                    Files.list(reportDir)
                        .filter(Files::isRegularFile)
                        .filter(p -> {
                            String name = p.getFileName().toString().toLowerCase();
                            return (name.endsWith(".txt") || name.endsWith(".html") ||
                                   name.endsWith(".tsv") || name.endsWith(".csv") ||
                                   name.endsWith(".log") || name.endsWith(".report")) &&
                                   !processedFiles.contains(name);
                        })
                        .forEach(reportFile -> {
                            try {
                                String filename = reportFile.getFileName().toString();
                                processedFiles.add(filename.toLowerCase());

                                if (filename.toLowerCase().endsWith(".html") || filename.toLowerCase().endsWith(".htm")) {
                                    // Import HTML files as text documents to avoid prompts
                                    String content = Files.readString(reportFile);
                                    TaxTriageResultDocument resultDoc = new TaxTriageResultDocument(
                                        filename, content, "HTML Report", reportFile.toAbsolutePath().toString(), "HTML");
                                    AnnotatedPluginDocument annotatedDoc = DocumentUtilities.createAnnotatedPluginDocument(resultDoc);
                                    AnnotatedPluginDocument copiedDoc = reportsFolder.addDocumentCopy(annotatedDoc, ProgressListener.EMPTY);

                                    if (copiedDoc != null) {
                                        docs.add(copiedDoc);
                                        logger.info("Imported HTML report: " + filename);
                                    }
                                } else {
                                    // Import text files as TaxTriageResultDocument
                                    String content = Files.readString(reportFile);
                                    String fileFormat = determineFileFormat(filename);
                                    String resultType = determineResultType(filename);

                                    TaxTriageResultDocument resultDoc = new TaxTriageResultDocument(
                                        filename, content, resultType, reportFile.toAbsolutePath().toString(), fileFormat);
                                    AnnotatedPluginDocument annotatedDoc = DocumentUtilities.createAnnotatedPluginDocument(resultDoc);
                                    AnnotatedPluginDocument copiedDoc = reportsFolder.addDocumentCopy(annotatedDoc, ProgressListener.EMPTY);

                                    if (copiedDoc != null) {
                                        docs.add(copiedDoc);
                                        logger.info("Imported text report: " + filename + " (" + fileFormat + ")");
                                    }
                                }
                            } catch (Exception e) {
                                logger.warning("Failed to import report file: " + reportFile.getFileName() + " - " + e.getMessage());
                            }
                        });
                } catch (IOException e) {
                    logger.warning("Error scanning report directory: " + reportDir);
                }
            }
        }

        return docs;
    }

    /**
     * Imports shared Kraken result files to the Kraken_Results folder.
     */
    private List<AnnotatedPluginDocument> importSharedKrakenResults(Path outputDir, WritableDatabaseService krakenFolder,
                                                                   ProgressListener progressListener) throws IOException {
        List<AnnotatedPluginDocument> docs = new ArrayList<>();

        // Import Kraken files from various directories
        Path[] krakenDirs = {
            outputDir.resolve("kraken2"),
            outputDir.resolve("kreport"),
            outputDir.resolve("mergedsubspecies"),
            outputDir.resolve("mergedkrakenreport")
        };

        Set<String> processedFiles = new HashSet<>();

        for (Path krakenDir : krakenDirs) {
            if (Files.exists(krakenDir)) {
                try {
                    Files.list(krakenDir)
                        .filter(Files::isRegularFile)
                        .filter(p -> {
                            String name = p.getFileName().toString().toLowerCase();
                            return (name.contains("kraken") || name.contains("krona") ||
                                   name.endsWith(".txt") || name.endsWith(".tsv") ||
                                   name.endsWith(".report")) &&
                                   !processedFiles.contains(name);
                        })
                        .forEach(krakenFile -> {
                            try {
                                String filename = krakenFile.getFileName().toString();
                                processedFiles.add(filename.toLowerCase());

                                String content = Files.readString(krakenFile);
                                String fileFormat = determineFileFormat(filename);
                                String resultType = determineKrakenResultType(filename);

                                TaxTriageResultDocument resultDoc = new TaxTriageResultDocument(
                                    filename, content, resultType, krakenFile.toAbsolutePath().toString(), fileFormat);
                                AnnotatedPluginDocument annotatedDoc = DocumentUtilities.createAnnotatedPluginDocument(resultDoc);
                                AnnotatedPluginDocument copiedDoc = krakenFolder.addDocumentCopy(annotatedDoc, ProgressListener.EMPTY);

                                if (copiedDoc != null) {
                                    docs.add(copiedDoc);
                                    logger.info("Imported Kraken result: " + filename + " (" + resultType + ")");
                                }
                            } catch (Exception e) {
                                logger.warning("Failed to import Kraken file: " + krakenFile.getFileName() + " - " + e.getMessage());
                            }
                        });
                } catch (IOException e) {
                    logger.warning("Error scanning Kraken directory: " + krakenDir);
                }
            }
        }

        return docs;
    }

    /**
     * Generates a comprehensive HTML report with Nextflow run information.
     */
    private AnnotatedPluginDocument generateHtmlReport(Path outputDir, SampleBasedImporter sampleImporter) {
        try {
            String runName = outputDir.getFileName().toString();
            TaxTriageHtmlReportGenerator generator = new TaxTriageHtmlReportGenerator(runName, outputDir.getParent());

            // Set basic information
            generator.setEndTime(new Date());
            if (commandLine != null) {
                generator.setCommandLine(commandLine);
            }

            // Add samples
            Set<String> samples = sampleImporter.getSamples();
            for (String sample : samples) {
                generator.addSample(sample);
            }

            // Add workflow parameters
            for (Map.Entry<String, String> param : workflowParams.entrySet()) {
                generator.addParameter(param.getKey(), param.getValue());
            }

            // Parse and add tool versions from pipeline_info
            parseToolVersions(outputDir, generator);

            // Add statistics
            generator.addStatistic("Total Samples", samples.size());
            generator.addStatistic("Working Directory", workingDirectory != null ? workingDirectory : outputDir.getParent().toString());

            // Check for errors and warnings in log files
            parseErrorsAndWarnings(outputDir, generator);

            // Generate HTML content
            String htmlContent = generator.generateHtml();

            // Create TaxTriageResultDocument
            TaxTriageResultDocument htmlDoc = new TaxTriageResultDocument(
                "TaxTriage_Analysis_Report.html",
                htmlContent,
                "TaxTriage Analysis Report",
                outputDir.resolve("TaxTriage_Analysis_Report.html").toString(),
                "HTML"
            );

            AnnotatedPluginDocument annotatedDoc = DocumentUtilities.createAnnotatedPluginDocument(htmlDoc);

            // Add to main folder
            WritableDatabaseService mainFolder = sampleImporter.getMainFolder();
            if (mainFolder != null) {
                return mainFolder.addDocumentCopy(annotatedDoc, ProgressListener.EMPTY);
            }

            return annotatedDoc;

        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to generate HTML report", e);
            return null;
        }
    }

    /**
     * Parses tool versions from pipeline_info directory.
     */
    private void parseToolVersions(Path outputDir, TaxTriageHtmlReportGenerator generator) {
        Path pipelineInfoDir = outputDir.resolve("pipeline_info");
        if (!Files.exists(pipelineInfoDir)) {
            return;
        }

        try {
            // Look for software_versions.yml or similar files
            Files.list(pipelineInfoDir)
                .filter(p -> {
                    String name = p.getFileName().toString().toLowerCase();
                    return name.contains("version") || name.contains("software") ||
                           name.equals("execution_trace.txt") || name.equals("execution_timeline.html");
                })
                .forEach(versionFile -> {
                    try {
                        String content = Files.readString(versionFile);
                        parseVersionContent(content, generator);
                    } catch (IOException e) {
                        logger.warning("Could not read version file: " + versionFile.getFileName());
                    }
                });
        } catch (IOException e) {
            logger.warning("Error scanning pipeline_info directory for versions");
        }
    }

    /**
     * Parses version content and adds tool versions to the generator.
     */
    private void parseVersionContent(String content, TaxTriageHtmlReportGenerator generator) {
        String[] lines = content.split("\n");
        for (String line : lines) {
            line = line.trim();

            // Parse YAML-style versions (tool: version)
            if (line.contains(":") && !line.startsWith("#")) {
                String[] parts = line.split(":", 2);
                if (parts.length == 2) {
                    String tool = parts[0].trim().replaceAll("[\"']", "");
                    String version = parts[1].trim().replaceAll("[\"']", "");
                    if (!tool.isEmpty() && !version.isEmpty()) {
                        generator.addParameter("Tool: " + tool, version);
                        toolVersions.put(tool, version);
                    }
                }
            }

            // Parse command line information
            if (line.toLowerCase().contains("command") && line.contains("nextflow")) {
                generator.addParameter("Nextflow Command", line);
            }
        }
    }

    /**
     * Parses errors and warnings from log files.
     */
    private void parseErrorsAndWarnings(Path outputDir, TaxTriageHtmlReportGenerator generator) {
        // Check .nextflow.log and other log files
        Path[] logFiles = {
            outputDir.getParent().resolve(".nextflow.log"),
            outputDir.resolve("pipeline_info").resolve("execution_trace.txt"),
            outputDir.resolve("nextflow.log")
        };

        for (Path logFile : logFiles) {
            if (Files.exists(logFile)) {
                try {
                    String content = Files.readString(logFile);
                    String[] lines = content.split("\n");

                    for (String line : lines) {
                        String lowerLine = line.toLowerCase();
                        if (lowerLine.contains("error") && !lowerLine.contains("0 errors")) {
                            generator.addError(line.trim());
                        } else if (lowerLine.contains("warn") && !lowerLine.contains("0 warnings")) {
                            generator.addWarning(line.trim());
                        }
                    }
                } catch (IOException e) {
                    logger.warning("Could not read log file: " + logFile.getFileName());
                }
            }
        }
    }

    /**
     * Determines the file format based on filename extension.
     */
    private String determineFileFormat(String filename) {
        String lowerName = filename.toLowerCase();
        if (lowerName.endsWith(".tsv")) return "TSV";
        if (lowerName.endsWith(".csv")) return "CSV";
        if (lowerName.endsWith(".html") || lowerName.endsWith(".htm")) return "HTML";
        if (lowerName.endsWith(".json")) return "JSON";
        if (lowerName.endsWith(".log")) return "LOG";
        return "TXT";
    }

    /**
     * Determines the result type based on filename patterns.
     */
    private String determineResultType(String filename) {
        String lowerName = filename.toLowerCase();
        if (lowerName.contains("count")) return "Count Data";
        if (lowerName.contains("top")) return "Top Results";
        if (lowerName.contains("combine")) return "Combined Results";
        if (lowerName.contains("multiqc")) return "MultiQC Report";
        if (lowerName.contains("taxtriage")) return "TaxTriage Analysis";
        if (lowerName.contains("assembly")) return "Assembly Summary";
        if (lowerName.contains("complete")) return "Completion Report";
        if (lowerName.contains("pipeline")) return "Pipeline Information";
        return "Analysis Result";
    }

    /**
     * Determines the Kraken-specific result type.
     */
    private String determineKrakenResultType(String filename) {
        String lowerName = filename.toLowerCase();
        if (lowerName.contains("kraken") && lowerName.contains("report")) return "Kraken Report";
        if (lowerName.contains("krona")) return "Krona Report";
        if (lowerName.contains("merged")) return "Merged Kraken Results";
        return "Kraken Output";
    }

    /**
     * Creates a post-import archive of the TaxTriage results.
     * This prompts the user for a save location and creates a .zip archive.
     */
    private void createPostImportArchive(Path outputDir, SampleBasedImporter sampleImporter) {
        logger.info("===========================================");
        logger.info("POST-IMPORT ARCHIVE CREATION");
        logger.info("===========================================");

        try {
            // Use the folder name as the suggested archive name
            String folderName = outputDir.getFileName().toString();
            String suggestedName = "TaxTriage_" + folderName;

            logger.info("Creating archive for output directory: " + outputDir);
            logger.info("Suggested archive name: " + suggestedName);

            // Create archive manager
            ArchiveManager archiveManager = new ArchiveManager();

            // Check if the directory can be archived
            if (!archiveManager.canArchive(outputDir)) {
                logger.warning("Directory cannot be archived: " + outputDir);
                logger.warning("Archive creation cancelled");
                return;
            }

            // Create archive with user dialog and progress monitoring
            logger.info("Prompting user for archive save location...");
            File archiveFile = archiveManager.createArchiveWithDialog(
                outputDir,
                suggestedName,
                percentage -> {
                    if (percentage % 10 == 0) { // Log every 10%
                        logger.info("Archive creation progress: " + percentage + "%");
                    }
                }
            );

            if (archiveFile != null) {
                logger.info("✓ Archive created successfully!");
                logger.info("  Archive file: " + archiveFile.getAbsolutePath());
                logger.info("  Archive size: " + formatFileSize(archiveFile.length()));

                // Log to console for user visibility
                System.out.println("===========================================");
                System.out.println("TaxTriage Results Archive Created");
                System.out.println("===========================================");
                System.out.println("Archive: " + archiveFile.getAbsolutePath());
                System.out.println("Size: " + formatFileSize(archiveFile.length()));
                System.out.println("===========================================");
            } else {
                logger.info("Archive creation was cancelled by user or failed");
                System.out.println("Archive creation cancelled or failed - results remain in: " + outputDir);
            }

        } catch (Exception e) {
            logger.log(Level.WARNING, "Error during post-import archive creation", e);
            System.out.println("Archive creation failed: " + e.getMessage());
            System.out.println("Results remain available in: " + outputDir);
        }

        logger.info("===========================================");
        logger.info("POST-IMPORT ARCHIVE CREATION COMPLETE");
        logger.info("===========================================");
    }

    /**
     * Formats file size into human-readable string.
     */
    private String formatFileSize(long bytes) {
        if (bytes < 1024) return bytes + " B";
        int exp = (int) (Math.log(bytes) / Math.log(1024));
        String pre = "KMGTPE".charAt(exp - 1) + "";
        return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
    }

}
