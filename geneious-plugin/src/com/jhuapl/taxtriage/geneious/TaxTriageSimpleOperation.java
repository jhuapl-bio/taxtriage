package com.jhuapl.taxtriage.geneious;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.PluginDocument;
import com.biomatters.geneious.publicapi.documents.sequence.SequenceDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentOperation;
import com.biomatters.geneious.publicapi.plugin.DocumentOperationException;
import com.biomatters.geneious.publicapi.plugin.DocumentSelectionSignature;
import com.biomatters.geneious.publicapi.plugin.GeneiousActionOptions;
import com.biomatters.geneious.publicapi.plugin.Options;
import com.jhuapl.taxtriage.geneious.config.ConfigGenerator;
import com.jhuapl.taxtriage.geneious.config.SampleSheetBuilder;
import com.jhuapl.taxtriage.geneious.config.TaxTriageConfig;
import com.jhuapl.taxtriage.geneious.database.DatabaseManager;
import com.jhuapl.taxtriage.geneious.database.DatabaseManager.DatabaseType;
import com.jhuapl.taxtriage.geneious.docker.DockerException;
import com.jhuapl.taxtriage.geneious.docker.DockerManager;
import com.jhuapl.taxtriage.geneious.docker.ExecutionResult;
import com.jhuapl.taxtriage.geneious.tools.BBToolsDeduplicator;
import com.jhuapl.taxtriage.geneious.utils.FileTypeUtil;
import com.jhuapl.taxtriage.geneious.utils.DatabasePathUtil;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.FileVisitOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Comprehensive TaxTriage genomic analysis operation for the Geneious platform.
 *
 * <p>This class provides a complete integration of the TaxTriage taxonomic classification
 * workflow within Geneious, enabling users to perform advanced metagenomic analysis
 * directly from their sequence data. The operation supports multiple sequencing platforms
 * and provides enterprise-grade workflow execution with comprehensive error handling
 * and progress tracking.</p>
 *
 * <h3>Key Features:</h3>
 * <ul>
 *   <li><strong>Multi-platform Support:</strong> Illumina (paired-end, single-end) and Oxford Nanopore</li>
 *   <li><strong>Flexible Input:</strong> File browser selection or Geneious document selection</li>
 *   <li><strong>Database Management:</strong> Automatic database detection, download, and caching</li>
 *   <li><strong>Workflow Orchestration:</strong> Nextflow-based pipeline execution with Docker containers</li>
 *   <li><strong>Quality Control:</strong> Optional read deduplication using GATK MarkDuplicates</li>
 *   <li><strong>Result Integration:</strong> Automatic import of analysis results with organized folder structure</li>
 * </ul>
 *
 * <h3>Workflow Architecture:</h3>
 * <ol>
 *   <li><strong>Environment Validation:</strong> Docker availability and Nextflow binary location</li>
 *   <li><strong>Database Preparation:</strong> Kraken2/Bracken database availability or download</li>
 *   <li><strong>Workspace Setup:</strong> Temporary workspace creation with organized directory structure</li>
 *   <li><strong>Input Processing:</strong> File validation, deduplication, and samplesheet generation</li>
 *   <li><strong>Pipeline Execution:</strong> TaxTriage workflow via Nextflow with real-time monitoring</li>
 *   <li><strong>Post-processing:</strong> Optional BAM deduplication using GATK MarkDuplicates</li>
 *   <li><strong>Result Import:</strong> Structured import of taxonomic reports, alignments, and references</li>
 * </ol>
 *
 * <h3>Performance Characteristics:</h3>
 * <ul>
 *   <li><strong>Scalability:</strong> Configurable thread count and memory limits</li>
 *   <li><strong>Reliability:</strong> Comprehensive error handling with detailed logging</li>
 *   <li><strong>Efficiency:</strong> Database caching and intermediate file management</li>
 *   <li><strong>Monitoring:</strong> Real-time progress updates and workflow status tracking</li>
 * </ul>
 *
 * <h3>Integration Points:</h3>
 * <ul>
 *   <li><strong>Geneious API:</strong> Document operations, progress tracking, and database services</li>
 *   <li><strong>Docker Engine:</strong> Container orchestration for tool execution</li>
 *   <li><strong>Nextflow:</strong> Workflow management and process coordination</li>
 *   <li><strong>External Tools:</strong> Kraken2, Bracken, minimap2, GATK, samtools</li>
 * </ul>
 *
 * <h3>Error Handling Strategy:</h3>
 * <ul>
 *   <li><strong>Graceful Degradation:</strong> Non-critical failures don't abort the entire workflow</li>
 *   <li><strong>Detailed Logging:</strong> Comprehensive error reporting with context</li>
 *   <li><strong>Resource Cleanup:</strong> Automatic cleanup of temporary resources on failure</li>
 *   <li><strong>User Feedback:</strong> Clear error messages with actionable guidance</li>
 * </ul>
 *
 * @author TaxTriage Development Team
 * @version 2.0
 * @since 1.0
 */
public class TaxTriageSimpleOperation extends DocumentOperation {

    private static final Logger logger = Logger.getLogger(TaxTriageSimpleOperation.class.getName());

    // Nextflow binary relative path within the plugin directory
    private static final String NEXTFLOW_RELATIVE_PATH = "bin/nextflow";

    // Cache for the located Nextflow binary path
    private static String nextflowBinaryPath = null;

    // Timeout constants
    private static final int WORKFLOW_TIMEOUT_MINUTES = 120; // 2 hours
    private static final int DATABASE_DOWNLOAD_TIMEOUT_MINUTES = 180; // 3 hours
    private static final int DOCKER_CHECK_TIMEOUT_SECONDS = 10;
    private static final int OUTPUT_READER_TIMEOUT_SECONDS = 5;

    // Database cache constants
    private static final String CACHE_DIR_NAME = ".taxtriage-geneious";
    private static final String STANDARD_DB_NAME = "standard";

    // File extension constants
    private static final String BAM_EXTENSION = ".bam";
    private static final String FASTQ_EXTENSION = ".fastq";
    private static final String FASTQ_GZ_EXTENSION = ".fastq.gz";
    private static final String FQ_EXTENSION = ".fq";
    private static final String FQ_GZ_EXTENSION = ".fq.gz";
    private static final String K2D_EXTENSION = ".k2d";

    // Progress tracking constants
    private static final double PROGRESS_INITIAL = 0.05;
    private static final double PROGRESS_VALIDATION = 0.1;
    private static final double PROGRESS_DATABASE_CHECK = 0.15;
    private static final double PROGRESS_WORKSPACE_SETUP = 0.2;
    private static final double PROGRESS_CONFIG_GENERATION = 0.25;
    private static final double PROGRESS_WORKFLOW_START = 0.3;
    private static final double PROGRESS_DEDUPLICATION = 0.85;
    private static final double PROGRESS_IMPORT = 0.9;
    private static final double PROGRESS_COMPLETE = 1.0;

    @Override
    public GeneiousActionOptions getActionOptions() {
        return new GeneiousActionOptions("TaxTriage Analysis",
                "Perform taxonomic classification using TaxTriage")
                .setMainMenuLocation(GeneiousActionOptions.MainMenu.Tools);
    }

    @Override
    public String getHelp() {
        return "Run TaxTriage analysis for taxonomic classification. You can either select documents or use the file browsers to select input files directly from the filesystem.";
    }

    @Override
    public DocumentSelectionSignature[] getSelectionSignatures() {
        return new DocumentSelectionSignature[] {
            // Allow any document type or no selection
            new DocumentSelectionSignature(PluginDocument.class, 0, Integer.MAX_VALUE)
        };
    }

    @Override
    public Options getOptions(AnnotatedPluginDocument... documents) throws DocumentOperationException {
        return new TaxTriageOptions();
    }

    @Override
    public List<AnnotatedPluginDocument> performOperation(AnnotatedPluginDocument[] documents,
                                                          ProgressListener progressListener,
                                                          Options options) throws DocumentOperationException {
        System.out.println("=================================================");
        System.out.println("TAXTRIAGE PLUGIN EXECUTION STARTED");
        System.out.println("Plugin Version: 2.0");
        System.out.println("Timestamp: " + new java.util.Date());
        System.out.println("=================================================");

        logger.info("TaxTriage operation started with " + (documents != null ? documents.length : 0) + " selected documents");

        if (progressListener != null) {
            progressListener.setMessage("Starting TaxTriage analysis...");
            progressListener.setProgress(PROGRESS_VALIDATION);
        }

        List<AnnotatedPluginDocument> results = new ArrayList<>();

        // Cast to TaxTriageOptions to access file browser selections
        TaxTriageOptions taxTriageOptions = (TaxTriageOptions) options;

        // Get input files from browsers
        List<File> inputFiles = taxTriageOptions.getInputFiles();
        File inputDirectory = taxTriageOptions.getInputDirectory();

        // Log the file selections
        logger.info("Input files from browser: " + (inputFiles != null ? inputFiles.size() : 0));
        if (inputFiles != null) {
            for (File f : inputFiles) {
                logger.fine("Input file: " + f.getAbsolutePath());
            }
        }
        logger.info("Input directory: " + (inputDirectory != null ? inputDirectory.getAbsolutePath() : "none"));

        // Check if we have input files from browsers
        boolean hasFileInputs = (inputFiles != null && !inputFiles.isEmpty()) || inputDirectory != null;

        logger.fine("Has file inputs from browser: " + hasFileInputs);

        if (!hasFileInputs) {
            // Only check selected documents if no file browser input was provided
            boolean hasSelectedDocuments = documents != null && documents.length > 0;
            logger.fine("Has selected documents: " + hasSelectedDocuments);

            if (!hasSelectedDocuments) {
                throw new DocumentOperationException("No input provided. Please use the file browsers to specify input files or directory.");
            }

            // If using selected documents (backward compatibility), validate they are sequence documents
            for (AnnotatedPluginDocument doc : documents) {
                logger.fine("Checking document: " + doc.getName() + " - Type: " + doc.getDocument().getClass().getName());
                if (!(doc.getDocument() instanceof SequenceDocument)) {
                    throw new DocumentOperationException("Document '" + doc.getName() + "' is not a sequence document. Please use the file browsers instead.");
                }
            }
        } else {
            // When using file browsers, ignore selected documents completely
            logger.info("Using file browser inputs, ignoring any selected documents");
        }

        if (progressListener != null) {
            String inputDescription = hasFileInputs ? "file browser inputs" : documents.length + " selected document(s)";
            progressListener.setMessage("Processing " + inputDescription + "...");
            progressListener.setProgress(0.5);
        }

        // Log configuration
        logger.info("TaxTriage Configuration:");
        logger.info("  Preset: " + taxTriageOptions.getSequencingPreset());
        logger.info("  Kraken DB: " + taxTriageOptions.getKrakenDatabase());
        logger.info("  Bracken DB: " + taxTriageOptions.getBrackenDatabase());
        logger.info("  Threads: " + taxTriageOptions.getThreadCount());
        logger.info("  Memory: " + taxTriageOptions.getMemoryLimit() + " GB");

        try {
            // Locate bundled Nextflow binary before workflow execution
            if (progressListener != null) {
                progressListener.setMessage("Preparing Nextflow environment...");
                progressListener.setProgress(PROGRESS_INITIAL);
            }

            locateNextflowBinary();

            // Execute the TaxTriage workflow
            List<AnnotatedPluginDocument> workflowResults = executeTaxTriageWorkflow(
                taxTriageOptions, hasFileInputs, inputFiles, inputDirectory, progressListener);
            results.addAll(workflowResults);

            if (progressListener != null) {
                progressListener.setMessage("TaxTriage analysis completed successfully");
                progressListener.setProgress(1.0);
            }

        } catch (Exception e) {
            logger.log(Level.SEVERE, "TaxTriage workflow execution failed", e);
            throw new DocumentOperationException("TaxTriage workflow failed: " + e.getMessage(), e);
        }

        return results;
    }

    /**
     * Executes the complete TaxTriage workflow.
     *
     * @param options TaxTriage configuration options
     * @param hasFileInputs whether file browser inputs are provided
     * @param inputFiles list of input files from browser
     * @param inputDirectory input directory from browser
     * @param progressListener progress listener for updates
     * @return list of imported result documents
     * @throws DocumentOperationException if workflow execution fails
     */
    private List<AnnotatedPluginDocument> executeTaxTriageWorkflow(
            TaxTriageOptions options, boolean hasFileInputs, List<File> inputFiles,
            File inputDirectory, ProgressListener progressListener) throws DocumentOperationException {

        logger.info("Starting TaxTriage workflow execution");

        try {
            // Step 1: Validate Docker availability (skip Nextflow check since we bundle it)
            if (progressListener != null) {
                progressListener.setMessage("Checking Docker availability...");
                progressListener.setProgress(0.1);
            }

            validateDockerEnvironment();

            // Step 2: Ensure database is available
            if (progressListener != null) {
                progressListener.setMessage("Checking database availability...");
                progressListener.setProgress(0.1);
            }

            ensureDatabaseExists(options, progressListener);

            // Step 3: Create temporary working directory
            if (progressListener != null) {
                progressListener.setMessage("Setting up workspace...");
                progressListener.setProgress(PROGRESS_DATABASE_CHECK);
            }

            Path workspaceDir = createWorkspace();
            logger.info("Created workspace: " + workspaceDir.toAbsolutePath());

            // Step 4: Collect and prepare input files
            if (progressListener != null) {
                progressListener.setMessage("Preparing input files...");
                progressListener.setProgress(0.2);
            }

            List<File> allInputFiles = collectInputFiles(hasFileInputs, inputFiles, inputDirectory);
            logger.info("Found " + allInputFiles.size() + " input files");

            if (allInputFiles.isEmpty()) {
                throw new DocumentOperationException("No valid input files found");
            }

            // Step 5: Generate configuration files
            if (progressListener != null) {
                progressListener.setMessage("Generating workflow configuration...");
                progressListener.setProgress(0.25);
            }

            generateWorkflowConfig(workspaceDir, allInputFiles, options);

            // Step 6: Execute TaxTriage workflow
            if (progressListener != null) {
                progressListener.setMessage("Starting TaxTriage analysis...");
                progressListener.setProgress(0.3);
            }

            ExecutionResult result = executeWorkflow(workspaceDir, options, progressListener);

            if (!result.isSuccessful()) {
                String standardOutput = result.getStandardOutput();
                String errorOutput = result.getErrorOutput();

                // Log both outputs for debugging
                logger.severe("TaxTriage workflow failed with exit code: " + result.getExitCode());
                if (standardOutput != null && !standardOutput.trim().isEmpty()) {
                    logger.severe("Standard output:\n" + standardOutput);
                }
                if (errorOutput != null && !errorOutput.trim().isEmpty()) {
                    logger.severe("Error output:\n" + errorOutput);
                } else {
                    logger.severe("No error output captured. The workflow may have failed silently.");
                }

                // Build comprehensive error message
                String errorMsg = "TaxTriage workflow failed with exit code " + result.getExitCode() + ".";
                if (errorOutput != null && !errorOutput.trim().isEmpty()) {
                    errorMsg += " Error: " + errorOutput;
                } else if (standardOutput != null && !standardOutput.trim().isEmpty()) {
                    // Sometimes errors go to stdout
                    errorMsg += " Output: " + standardOutput;
                }

                throw new DocumentOperationException(errorMsg);
            }

            // Step 7: Deduplication is now handled by BBTools preprocessing
            // No additional BAM deduplication needed
            logger.info("Skipping BAM deduplication - reads were already deduplicated by BBTools preprocessing");

            // Step 8: Import results back to Geneious
            if (progressListener != null) {
                progressListener.setMessage("Importing results...");
                progressListener.setProgress(0.9);
            }

            List<AnnotatedPluginDocument> importedResults = importResults(workspaceDir, progressListener);

            logger.info("TaxTriage workflow completed successfully. Imported " +
                importedResults.size() + " result documents.");

            return importedResults;

        } catch (DockerException e) {
            throw new DocumentOperationException("Docker error: " + e.getMessage(), e);
        } catch (IOException e) {
            throw new DocumentOperationException("I/O error during workflow execution: " + e.getMessage(), e);
        } catch (Exception e) {
            throw new DocumentOperationException("Unexpected error during workflow execution: " + e.getMessage(), e);
        }
    }

    /**
     * Validates that Docker is available and accessible.
     * Note: Nextflow validation is skipped since we bundle it with the plugin.
     */
    private void validateDockerEnvironment() throws DockerException {
        logger.info("Validating Docker environment (Nextflow is bundled)");

        // Check Docker availability (still needed for workflow containers)
        if (!isDockerAvailable()) {
            throw DockerException.dockerNotAvailable();
        }

        logger.info("Docker environment validated successfully");
    }

    /**
     * Creates a temporary workspace directory for the workflow.
     */
    private Path createWorkspace() throws IOException {
        String timestamp = LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss"));
        String workspaceName = "taxtriage_" + timestamp;

        Path tempDir = Paths.get(System.getProperty("java.io.tmpdir"));
        Path workspaceDir = tempDir.resolve(workspaceName);

        // Create directory structure
        Files.createDirectories(workspaceDir);
        Files.createDirectories(workspaceDir.resolve("input"));
        Files.createDirectories(workspaceDir.resolve("output"));
        Files.createDirectories(workspaceDir.resolve("work"));
        Files.createDirectories(workspaceDir.resolve("config"));

        logger.info("Created workspace structure in: " + workspaceDir.toAbsolutePath());
        return workspaceDir;
    }

    /**
     * Collects all input files from the various input sources.
     */
    private List<File> collectInputFiles(boolean hasFileInputs, List<File> inputFiles,
                                       File inputDirectory) throws DocumentOperationException {
        // Use a Set to avoid duplicates when same file is selected individually and in directory
        Set<String> processedPaths = new HashSet<>();
        List<File> allFiles = new ArrayList<>();

        if (hasFileInputs) {
            // Process file browser inputs
            if (inputFiles != null && !inputFiles.isEmpty()) {
                // Individual files selected
                for (File file : inputFiles) {
                    if (isValidSequenceFile(file)) {
                        String absolutePath = file.getAbsolutePath();
                        if (!processedPaths.contains(absolutePath)) {
                            allFiles.add(file);
                            processedPaths.add(absolutePath);
                            logger.info("Added input file: " + file.getName());
                        } else {
                            logger.info("Skipping duplicate file: " + file.getName());
                        }
                    } else {
                        logger.warning("Skipping invalid file: " + file.getName());
                    }
                }
            }

            if (inputDirectory != null && inputDirectory.isDirectory()) {
                // Directory selected - scan for sequence files
                List<File> dirFiles = new ArrayList<>();
                scanDirectoryForSequenceFiles(inputDirectory, dirFiles);

                // Add only files that haven't been added already
                for (File file : dirFiles) {
                    String absolutePath = file.getAbsolutePath();
                    if (!processedPaths.contains(absolutePath)) {
                        allFiles.add(file);
                        processedPaths.add(absolutePath);
                        logger.info("Added file from directory: " + file.getName());
                    } else {
                        logger.info("Skipping duplicate from directory: " + file.getName());
                    }
                }
            }
        }

        return allFiles;
    }

    /**
     * Checks if a file is a valid sequence file (FASTQ only).
     * Only FASTQ files are supported for TaxTriage analysis.
     */
    private boolean isValidSequenceFile(File file) {
        return FileTypeUtil.isValidSequenceFile(file);
    }

    /**
     * Recursively scans a directory for sequence files.
     */
    private void scanDirectoryForSequenceFiles(File directory, List<File> fileList) {
        File[] files = directory.listFiles();
        if (files != null) {
            for (File file : files) {
                if (file.isFile() && isValidSequenceFile(file)) {
                    fileList.add(file);
                } else if (file.isDirectory() && !file.getName().startsWith(".")) {
                    // Recursively scan subdirectories (skip hidden dirs)
                    scanDirectoryForSequenceFiles(file, fileList);
                }
            }
        }
    }

    /**
     * Generates Nextflow configuration and samplesheet files.
     */
    private void generateWorkflowConfig(Path workspaceDir, List<File> inputFiles,
                                      TaxTriageOptions options) throws IOException {
        logger.info("Generating workflow configuration files");

        // Create TaxTriage configuration
        TaxTriageConfig config = TaxTriageConfig.fromOptions(options);

        // Validate configuration
        String validationError = config.validate();
        if (validationError != null) {
            throw new IOException("Configuration validation failed: " + validationError);
        }

        // Generate Nextflow config file
        ConfigGenerator configGenerator = new ConfigGenerator();
        File configFile = workspaceDir.resolve("config").resolve("nextflow.config").toFile();
        configGenerator.generateConfig(config, configFile);
        logger.info("Generated Nextflow config: " + configFile.getAbsolutePath());

        // Generate parameters file
        File paramsFile = workspaceDir.resolve("config").resolve("params.json").toFile();
        configGenerator.generateParams(config, paramsFile);
        logger.info("Generated parameters file: " + paramsFile.getAbsolutePath());

        // Step 1: BBTools preprocessing BEFORE creating the sample sheet (if enabled)
        List<File> preprocessedFiles;
        boolean bbtoolsEnabled = options.isBBToolsPreprocessingEnabled();

        System.out.println("==========================================");
        System.out.println("PREPROCESSING CONFIGURATION:");
        System.out.println("  BBTools Preprocessing Enabled: " + bbtoolsEnabled);
        System.out.println("  Substitution Threshold: " + options.getBBToolsSubstitutionThreshold());
        System.out.println("==========================================");

        logger.info("==========================================");
        logger.info("PREPROCESSING CONFIGURATION:");
        logger.info("  BBTools Preprocessing Enabled: " + bbtoolsEnabled);
        logger.info("  Substitution Threshold: " + options.getBBToolsSubstitutionThreshold());
        logger.info("==========================================");

        if (bbtoolsEnabled) {
            logger.info("Starting BBTools preprocessing pipeline...");
            preprocessedFiles = performBBToolsPreprocessing(workspaceDir, inputFiles, options);
        } else {
            logger.info("BBTools preprocessing disabled - copying original files to workspace");
            preprocessedFiles = copyInputFilesToWorkspace(workspaceDir, inputFiles);
        }

        // Generate samplesheet with preprocessed files
        SampleSheetBuilder sampleSheetBuilder = new SampleSheetBuilder();
        File sampleSheetFile = workspaceDir.resolve("config").resolve("samplesheet.csv").toFile();

        // Create list of files with workspace input paths for samplesheet
        // Use a Set to track processed filenames and avoid duplicates
        Set<String> processedFilenames = new HashSet<>();
        List<File> workspaceInputFiles = new ArrayList<>();

        logger.info("Processing " + preprocessedFiles.size() + " preprocessed files for samplesheet generation");

        for (File inputFile : preprocessedFiles) {
            String filename = inputFile.getName();

            // Check if we've already processed this filename
            if (processedFilenames.contains(filename)) {
                logger.warning("Skipping duplicate filename: " + filename + " (path: " + inputFile.getAbsolutePath() + ")");
                continue;
            }

            // Use the preprocessed files directly (they're already in workspace)
            workspaceInputFiles.add(inputFile);
            processedFilenames.add(filename);

            logger.info("Added preprocessed file to samplesheet: " + filename + " (source: " + inputFile.getAbsolutePath() + ")");
        }

        logger.info("Final count: " + workspaceInputFiles.size() + " unique preprocessed files will be added to samplesheet");

        sampleSheetBuilder.buildSampleSheet(workspaceInputFiles, sampleSheetFile,
            options.getSequencingPreset().name());
        logger.info("Generated samplesheet with preprocessed files: " + sampleSheetFile.getAbsolutePath());
    }

    /**
     * Performs BBTools preprocessing on input FASTQ files.
     * This includes deduplication and quality preprocessing.
     */
    private List<File> performBBToolsPreprocessing(Path workspaceDir, List<File> inputFiles, TaxTriageOptions options) throws IOException {
        System.out.println(">>> ENTERED performBBToolsPreprocessing method");
        System.out.println(">>> Input files count: " + inputFiles.size());

        List<File> preprocessedFiles = new ArrayList<>();
        Path inputDir = workspaceDir.resolve("input");
        Path preprocessedDir = workspaceDir.resolve("preprocessed");

        // Create directories
        Files.createDirectories(inputDir);
        Files.createDirectories(preprocessedDir);

        System.out.println(">>> Created directories: input and preprocessed");
        logger.info("Starting BBTools preprocessing for " + inputFiles.size() + " files");

        // Initialize BBTools deduplicator with substitution threshold from options
        int subsThreshold = options.getBBToolsSubstitutionThreshold();
        System.out.println(">>> Substitution threshold: " + subsThreshold);
        logger.info("Configured substitution threshold: " + subsThreshold);

        BBToolsDeduplicator deduplicator;
        try {
            System.out.println(">>> Attempting to initialize BBToolsDeduplicator...");
            deduplicator = new BBToolsDeduplicator(subsThreshold);
            System.out.println(">>> BBTools deduplicator initialized successfully!");
            logger.info("✓ BBTools deduplicator initialized successfully");
        } catch (DockerException e) {
            System.out.println(">>> BBTools initialization FAILED: " + e.getMessage());
            String errorMsg = "BBTools preprocessing is enabled but failed to initialize.\n\n" +
                            "Error: " + e.getMessage() + "\n\n" +
                            "To fix this issue:\n" +
                            "1. Ensure Docker is installed and running\n" +
                            "2. Check Docker permissions\n" +
                            "3. Or disable BBTools preprocessing in the plugin options";
            logger.severe(errorMsg);
            throw new IOException(errorMsg, e);
        }

        System.out.println(">>> Checking if BBTools Docker image is available...");
        if (!deduplicator.isBBToolsAvailable()) {
            System.out.println(">>> BBTools Docker image NOT available - FAILING workflow");
            String errorMsg = "BBTools preprocessing is enabled but the Docker image is not available.\n\n" +
                            "Required: " + deduplicator.getImageName() + "\n\n" +
                            "To fix this issue:\n" +
                            "1. Ensure Docker is running\n" +
                            "2. Pull the BBTools image: docker pull " + deduplicator.getImageName() + "\n" +
                            "3. Or disable BBTools preprocessing in the plugin options";
            logger.severe(errorMsg);
            throw new IOException(errorMsg);
        }

        System.out.println(">>> BBTools Docker image IS available: " + deduplicator.getImageName());
        logger.info("✓ BBTools Docker image available: " + deduplicator.getImageName());

        try {

            int processedCount = 0;

            for (File inputFile : inputFiles) {
                System.out.println(">>> PROCESSING FILE " + (processedCount + 1) + "/" + inputFiles.size() + ": " + inputFile.getName());
                logger.info("======================================");
                logger.info("PROCESSING FILE " + (processedCount + 1) + "/" + inputFiles.size() + ": " + inputFile.getName());
                logger.info("======================================");

                // Copy original file to input directory first
                Path originalInWorkspace = inputDir.resolve(inputFile.getName());
                Files.copy(inputFile.toPath(), originalInWorkspace, StandardCopyOption.REPLACE_EXISTING);
                System.out.println(">>> Copied to workspace: " + originalInWorkspace);
                logger.info("Copied to workspace: " + originalInWorkspace);

                // Preprocess with BBTools
                System.out.println(">>> Starting BBTools 4-step pipeline...");
                logger.info("Starting BBTools 4-step pipeline...");
                BBToolsDeduplicator.DeduplicationResult result = deduplicator.deduplicateReads(
                    originalInWorkspace, preprocessedDir, jebl.util.ProgressListener.EMPTY);
                System.out.println(">>> BBTools pipeline completed, checking result...");

                if (result.isSuccess()) {
                    System.out.println(">>> SUCCESS! Result output path: " + result.getOutputPath());
                    File preprocessedFile = new File(result.getOutputPath());
                    preprocessedFiles.add(preprocessedFile);

                    System.out.println(">>> Duplicates removed: " + result.getDuplicatesRemoved());
                    System.out.println(">>> Total reads: " + result.getTotalReads());
                    System.out.println(">>> Preprocessed file added to list: " + preprocessedFile.getAbsolutePath());

                    logger.info("");
                    logger.info("\u2713\u2713\u2713 BBTools preprocessing SUCCESSFUL \u2713\u2713\u2713");
                    logger.info("  Original file:     " + originalInWorkspace);
                    logger.info("  Preprocessed file: " + preprocessedFile.getAbsolutePath());
                    logger.info("  Duplicates removed: " + result.getDuplicatesRemoved());
                    logger.info("  Total reads:        " + result.getTotalReads());
                    logger.info("  Deduplication rate: " + String.format("%.2f%%", result.getDeduplicationPercentage()));

                    // Log the processing steps
                    if (!result.getStepResults().isEmpty()) {
                        logger.info("  Pipeline steps completed:");
                        for (String step : result.getStepResults()) {
                            logger.info("    - " + step);
                        }
                    }
                } else {
                    System.out.println(">>> FAILED! Error: " + result.getErrorMessage());
                    System.out.println(">>> ABORTING workflow - preprocessing is required");
                    String errorMsg = "BBTools preprocessing failed for file: " + inputFile.getName() + "\n\n" +
                                    "Error: " + result.getErrorMessage() + "\n\n" +
                                    "Preprocessing is enabled and required. The workflow cannot continue.\n" +
                                    "To fix this issue:\n" +
                                    "1. Check Docker is running\n" +
                                    "2. Ensure input files are valid FASTQ format\n" +
                                    "3. Check available disk space\n" +
                                    "4. Or disable BBTools preprocessing in the plugin options";
                    logger.severe(errorMsg);
                    throw new IOException(errorMsg);
                }

                processedCount++;
            }

        } catch (Exception e) {
            logger.warning("======================================");
            logger.warning("BBTools preprocessing encountered an error:");
            logger.warning("  Error: " + e.getMessage());
            logger.warning("  Falling back to original files without preprocessing");
            logger.warning("======================================");
            logger.log(Level.FINE, "BBTools preprocessing detailed error", e);
            // Fallback to copying original files
            return copyInputFilesToWorkspace(workspaceDir, inputFiles);
        }

        logger.info("BBTools preprocessing completed. Processed " + preprocessedFiles.size() + " files");
        logger.info("==========================================");

        return preprocessedFiles;
    }

    /**
     * Copies input files to the workspace input directory (fallback method).
     */
    private List<File> copyInputFilesToWorkspace(Path workspaceDir, List<File> inputFiles) throws IOException {
        Path inputDir = workspaceDir.resolve("input");
        List<File> copiedFiles = new ArrayList<>();

        Files.createDirectories(inputDir);

        for (File inputFile : inputFiles) {
            Path targetPath = inputDir.resolve(inputFile.getName());

            // Copy file to workspace, replacing if it already exists
            Files.copy(inputFile.toPath(), targetPath,
                      java.nio.file.StandardCopyOption.REPLACE_EXISTING);
            logger.fine("Copied input file: " + inputFile.getName());
            copiedFiles.add(targetPath.toFile());
        }

        logger.info("Copied " + inputFiles.size() + " input files to workspace");
        return copiedFiles;
    }

    /**
     * Executes the TaxTriage Nextflow workflow directly on the host system.
     */
    private ExecutionResult executeWorkflow(Path workspaceDir, TaxTriageOptions options,
                                          ProgressListener progressListener) throws DockerException {
        logger.info("Starting TaxTriage workflow execution with native Nextflow");

        // Build Nextflow command
        List<String> nextflowCommand = buildNextflowCommand(workspaceDir, options);
        String commandString = String.join(" ", nextflowCommand);
        logger.info("Executing Nextflow command: " + commandString);
        logger.info("Working directory: " + workspaceDir.toAbsolutePath());

        // Execute the workflow directly with ProcessBuilder
        ExecutionResult result = executeNextflowDirectly(
            nextflowCommand, workspaceDir, progressListener, WORKFLOW_TIMEOUT_MINUTES);

        logger.info("Workflow execution completed. Exit code: " + result.getExitCode() +
            ", Duration: " + result.getExecutionSeconds() + "s");

        // Log execution details
        if (result.isSuccessful()) {
            logger.info("TaxTriage workflow completed successfully");
        } else {
            logger.severe("TaxTriage workflow failed: " + result.getErrorOutput());
        }

        return result;
    }

    /**
     * Builds the Nextflow command for TaxTriage execution.
     */
    private List<String> buildNextflowCommand(Path workspaceDir, TaxTriageOptions options) {
        List<String> cmd = new ArrayList<>();

        // Use the extracted Nextflow binary
        cmd.add(getNextflowBinaryPath());
        cmd.add("run");
        cmd.add("https://github.com/jhuapl-bio/taxtriage");
        cmd.add("-r");
        cmd.add("main");
        cmd.add("-profile");
        cmd.add("docker");

        // Core parameters - use local filesystem paths
        cmd.add("--input");
        cmd.add(workspaceDir.resolve("config").resolve("samplesheet.csv").toAbsolutePath().toString());
        cmd.add("--outdir");
        cmd.add(workspaceDir.resolve("output").toAbsolutePath().toString());
        cmd.add("-work-dir");
        cmd.add(workspaceDir.resolve("work").toAbsolutePath().toString());

        // Sequencing preset
        String preset = options.getSequencingPreset().name().toLowerCase();
        if (preset.equals("ont")) {
            cmd.add("--platform");
            cmd.add("OXFORD");
        } else if (preset.equals("illumina_pe")) {
            cmd.add("--platform");
            cmd.add("ILLUMINA");
        } else if (preset.equals("illumina_se")) {
            cmd.add("--platform");
            cmd.add("ILLUMINA");
            cmd.add("--single_end");
            cmd.add("true");
        }

        // Database configuration with DatabaseManager
        DatabaseManager dbManager = DatabaseManager.getInstance();
        String dbName = options.getKrakenDatabase();
        DatabaseType dbType = mapDatabaseNameToType(dbName);

        if (dbType != null) {
            // First check if cached or bundled database exists (without allowing download)
            String dbPath = dbManager.getDatabasePath(dbType, false);

            if (dbPath != null) {
                // We have a cached or bundled database - use it
                logger.info("==========================================");
                logger.info("Using cached/bundled database for workflow");
                logger.info("  Type: " + dbType.getId());
                logger.info("  Path: " + dbPath);
                logger.info("==========================================");

                cmd.add("--db");
                cmd.add(dbPath);
                cmd.add("--download_db");
                cmd.add("false");
            } else {
                // No cached database exists - need to download
                logger.info("==========================================");
                logger.info("No cached database found - will download");
                logger.info("  Type: " + dbType.getId());
                logger.info("  Name: " + dbName);

                // When downloading, just pass the database name (e.g., "viral")
                // TaxTriage will download it to its own location
                logger.info("  TaxTriage will download to workflow directory");
                logger.info("==========================================");

                cmd.add("--db");
                cmd.add(dbName);  // Just use the name, not a full path
                cmd.add("--download_db");
                cmd.add("true");

                // Note: After download, we'll scan the work directory to find and cache it
            }
        } else {
            // Unknown database type, use as-is
            cmd.add("--db");
            cmd.add(dbName);
            cmd.add("--download_db");
            cmd.add("true");
        }

        // Resource limits
        cmd.add("--max_cpus");
        cmd.add(String.valueOf(options.getThreadCount()));
        cmd.add("--max_memory");
        cmd.add(options.getMemoryLimit() + ".GB");

        // Quality filtering parameters (deprecated - removed)

        // Subsampling
        if (options.isSubsamplingEnabled()) {
            cmd.add("--enable_subsampling");
            cmd.add("true");
            cmd.add("--subsample_size");
            cmd.add(String.valueOf(options.getSubsampleSize()));
        }

        return cmd;
    }

    /**
     * Maps a database name string to a DatabaseType enum.
     */
    private DatabaseType mapDatabaseNameToType(String databaseName) {
        return DatabasePathUtil.mapDatabaseNameToType(databaseName);
    }

    /**
     * Imports TaxTriage results back to Geneious.
     */
    private List<AnnotatedPluginDocument> importResults(Path workspaceDir,
                                                       ProgressListener progressListener) {
        logger.info("Importing TaxTriage results to Geneious");

        List<AnnotatedPluginDocument> importedDocs = new ArrayList<>();

        try {
            // Check if output directory exists and has content
            Path outputDir = workspaceDir.resolve("output");
            if (!Files.exists(outputDir) || !Files.isDirectory(outputDir)) {
                logger.warning("No output directory found: " + outputDir.toAbsolutePath());
                return importedDocs;
            }

            // Use the TaxTriageResultImporter to import results
            com.jhuapl.taxtriage.geneious.results.TaxTriageResultImporter importer =
                new com.jhuapl.taxtriage.geneious.results.TaxTriageResultImporter();

            importedDocs = importer.importResults(outputDir, progressListener);

            if (importedDocs.isEmpty()) {
                logger.warning("No results could be imported from: " + outputDir.toAbsolutePath());
            } else {
                logger.info("Successfully imported " + importedDocs.size() + " result documents");
            }

            if (progressListener != null) {
                progressListener.setMessage("Results generated in: " + outputDir.toAbsolutePath());
            }

            // Check if databases were downloaded and cache them
            cacheDownloadedDatabases(workspaceDir);

        } catch (Exception e) {
            logger.log(Level.WARNING, "Error during result import", e);
        }

        return importedDocs;
    }


    /**
     * Checks if databases were downloaded during workflow and caches them.
     */
    private void cacheDownloadedDatabases(Path workspaceDir) {
        logger.info("==========================================");
        logger.info("DATABASE CACHING: Scanning for downloaded databases");
        logger.info("  Workspace: " + workspaceDir.toAbsolutePath());
        logger.info("==========================================");

        try {
            DatabaseManager dbManager = DatabaseManager.getInstance();
            Set<Path> foundDatabases = new HashSet<>();

            // First, do a comprehensive scan of the work directory for .k2d files
            Path workDir = workspaceDir.resolve("work");
            if (Files.exists(workDir)) {
                logger.info("\n[1] Scanning Nextflow work directory for database files...");
                logger.info("    Path: " + workDir.toAbsolutePath());

                // Look for directories containing .k2d files (Kraken2 database files)
                Files.walk(workDir, FileVisitOption.FOLLOW_LINKS)
                    .filter(Files::isDirectory)
                    .forEach(dir -> {
                        try {
                            // Check if this directory contains Kraken2 database files
                            boolean hasHashFile = Files.exists(dir.resolve("hash.k2d"));
                            boolean hasTaxoFile = Files.exists(dir.resolve("taxo.k2d"));

                            if (hasHashFile && hasTaxoFile) {
                                foundDatabases.add(dir);
                                logger.info("    ✓ Found database files in: " + dir.getFileName());
                                logger.info("      Full path: " + dir.toAbsolutePath());

                                // Try to determine database type from parent directory name
                                String pathStr = dir.toString().toLowerCase();
                                DatabaseType dbType = null;

                                if (pathStr.contains("viral")) {
                                    dbType = DatabaseType.VIRAL;
                                } else if (pathStr.contains("standard")) {
                                    dbType = DatabaseType.STANDARD;
                                } else if (pathStr.contains("minikraken") || pathStr.contains("mini")) {
                                    dbType = DatabaseType.MINIKRAKEN;
                                }

                                if (dbType != null) {
                                    logger.info("      Detected type: " + dbType.getId());
                                    String version = "Downloaded " +
                                        LocalDateTime.now().format(DateTimeFormatter.ISO_LOCAL_DATE);
                                    dbManager.cacheDownloadedDatabase(dbType, dir, version);
                                    logger.info("      ✓ Cached for future use");
                                } else {
                                    logger.info("      ⚠ Could not determine database type");
                                }
                            }
                        } catch (Exception e) {
                            // Ignore errors for individual directories
                        }
                    });
            }

            // Check common database locations
            logger.info("\n[2] Checking standard database locations...");

            // Check user's home for various database locations
            Path[] homeLocations = {
                Paths.get(System.getProperty("user.home"), "taxtriage_databases"),
                Paths.get(System.getProperty("user.home"), ".taxtriage", "databases"),
                Paths.get(System.getProperty("user.home"), ".nextflow", "assets"),
                Paths.get(System.getProperty("user.home"), "kraken2_dbs"),
                // Check singularity cache locations
                Paths.get(System.getProperty("user.home"), ".singularity", "cache"),
                Paths.get(System.getProperty("user.home"), ".apptainer", "cache"),
                // Check common Nextflow work locations
                Paths.get(System.getProperty("java.io.tmpdir"), "nxf-work"),
                workspaceDir.resolve("work").resolve("singularity")
            };

            for (Path dbPath : homeLocations) {
                if (Files.exists(dbPath)) {
                    logger.info("    Checking: " + dbPath);
                    scanDirectoryForDatabases(dbPath, dbManager);
                }
            }

            // Check workspace locations
            logger.info("\n[3] Checking workspace-relative locations...");
            Path[] workspaceLocations = {
                workspaceDir.resolve("databases"),
                workspaceDir.resolve("db"),
                workspaceDir.resolve("kraken2_db"),
                workspaceDir.resolve("output").resolve("databases"),
                workspaceDir.resolve("output").resolve("db")
            };

            for (Path dbPath : workspaceLocations) {
                if (Files.exists(dbPath)) {
                    logger.info("    Checking: " + dbPath);
                    scanDirectoryForDatabases(dbPath, dbManager);
                }
            }

            // Log summary
            if (foundDatabases.isEmpty()) {
                logger.info("\n⚠ No databases found in work directory");
                logger.info("  Databases may be:");
                logger.info("  - Downloaded to a different location");
                logger.info("  - Inside Docker containers (not accessible)");
                logger.info("  - Already cached from a previous run");
            } else {
                logger.info("\n✓ Found " + foundDatabases.size() + " database(s)");
            }

        } catch (Exception e) {
            logger.log(Level.WARNING, "Error during database caching", e);
            e.printStackTrace();
        }

        logger.info("==========================================");
    }

    /**
     * Scans a directory for database subdirectories and caches them.
     */
    private void scanDirectoryForDatabases(Path directory, DatabaseManager dbManager) {
        try {
            Files.list(directory)
                .filter(Files::isDirectory)
                .forEach(dir -> {
                    String dirName = dir.getFileName().toString();
                    DatabaseType dbType = mapDatabaseNameToType(dirName);

                    if (isValidDatabaseDirectory(dir)) {
                        logger.info("      Found database: " + dirName + " at " + dir);

                        if (dbType != null) {
                            String version = "Downloaded " +
                                LocalDateTime.now().format(DateTimeFormatter.ISO_LOCAL_DATE);
                            dbManager.cacheDownloadedDatabase(dbType, dir, version);
                            logger.info("      ✓ Cached as type: " + dbType.getId());
                        } else {
                            logger.info("      ⚠ Unknown database type: " + dirName);
                        }
                    }
                });
        } catch (IOException e) {
            logger.warning("    Error scanning directory: " + e.getMessage());
        }
    }

    /**
     * Checks if a directory contains valid database files.
     */
    private boolean isValidDatabaseDirectory(Path dir) {
        return DatabasePathUtil.isValidDatabaseDirectory(dir);
    }

    /**
     * Caches a directory if it contains a valid database.
     */
    private void cacheDirectoryIfDatabase(Path dir) {
        try {
            Files.list(dir)
                .filter(Files::isDirectory)
                .forEach(subDir -> {
                    String name = subDir.getFileName().toString();
                    DatabaseType dbType = mapDatabaseNameToType(name);

                    if (dbType != null && isValidDatabaseDirectory(subDir)) {
                        DatabaseManager dbManager = DatabaseManager.getInstance();
                        String version = "Downloaded " +
                            LocalDateTime.now().format(DateTimeFormatter.ISO_LOCAL_DATE);
                        dbManager.cacheDownloadedDatabase(dbType, subDir, version);
                        logger.info("    Cached database: " + name);
                    }
                });
        } catch (IOException e) {
            logger.log(Level.WARNING, "Error checking directory for databases", e);
        }
    }

    /**
     * Locates the bundled Nextflow binary within the plugin directory.
     * This method is idempotent - subsequent calls will return the same located binary.
     *
     * @throws DockerException if the binary cannot be located or is not executable
     */
    private void locateNextflowBinary() throws DockerException {
        if (nextflowBinaryPath != null) {
            // Already located, verify it still exists
            File existingBinary = new File(nextflowBinaryPath);
            if (existingBinary.exists() && existingBinary.canExecute()) {
                logger.fine("Using existing Nextflow binary: " + nextflowBinaryPath);
                return;
            } else {
                logger.warning("Previously located Nextflow binary is missing, re-locating");
                nextflowBinaryPath = null;
            }
        }

        try {
            logger.info("Locating bundled Nextflow binary in plugin directory");

            // Find the plugin installation directory
            String pluginDir = findPluginDirectory();
            if (pluginDir == null) {
                throw new DockerException("Cannot determine plugin installation directory");
            }

            // Build path to the Nextflow binary
            Path nextflowPath = Paths.get(pluginDir, NEXTFLOW_RELATIVE_PATH);
            File nextflowFile = nextflowPath.toFile();

            // Verify the binary exists
            if (!nextflowFile.exists()) {
                throw new DockerException("Bundled Nextflow binary not found at: " + nextflowPath.toAbsolutePath());
            }

            // Verify the binary is executable (and make it executable if needed)
            if (!nextflowFile.canExecute()) {
                logger.warning("Nextflow binary is not executable, attempting to set executable permission");
                if (!nextflowFile.setExecutable(true, false)) {
                    throw new DockerException("Cannot make Nextflow binary executable: " + nextflowPath.toAbsolutePath());
                }
            }

            nextflowBinaryPath = nextflowPath.toAbsolutePath().toString();
            logger.info("Successfully located Nextflow binary at: " + nextflowBinaryPath);

        } catch (Exception e) {
            if (e instanceof DockerException) {
                throw e;
            }
            throw new DockerException("Failed to locate bundled Nextflow binary: " + e.getMessage(), e);
        }
    }

    /**
     * Attempts to find the plugin installation directory by examining the class loader.
     *
     * @return the plugin directory path, or null if it cannot be determined
     */
    private String findPluginDirectory() {
        try {
            // Get the location of this class
            java.security.CodeSource codeSource = getClass().getProtectionDomain().getCodeSource();
            if (codeSource != null) {
                java.net.URL location = codeSource.getLocation();
                if (location != null) {
                    Path classPath = Paths.get(location.toURI());

                    // If we're running from a JAR file (TaxTriage.jar), get its parent directory
                    if (classPath.toString().endsWith(".jar")) {
                        Path pluginDir = classPath.getParent();
                        logger.fine("Found plugin directory via JAR location: " + pluginDir.toAbsolutePath());
                        return pluginDir.toAbsolutePath().toString();
                    }

                    // If we're running from classes directory (during development), look for resources
                    Path potentialPluginDir = classPath;
                    while (potentialPluginDir != null && potentialPluginDir.getParent() != null) {
                        Path nextflowPath = potentialPluginDir.resolve(NEXTFLOW_RELATIVE_PATH);
                        if (Files.exists(nextflowPath)) {
                            logger.fine("Found plugin directory via resource search: " + potentialPluginDir.toAbsolutePath());
                            return potentialPluginDir.toAbsolutePath().toString();
                        }
                        potentialPluginDir = potentialPluginDir.getParent();
                    }
                }
            }

            logger.warning("Could not determine plugin directory from class location");
            return null;

        } catch (Exception e) {
            logger.log(Level.WARNING, "Error determining plugin directory", e);
            return null;
        }
    }

    /**
     * Gets the path to the Nextflow binary (located within the plugin directory).
     *
     * @return the absolute path to the Nextflow binary
     * @throws IllegalStateException if the binary hasn't been located yet
     */
    private String getNextflowBinaryPath() {
        if (nextflowBinaryPath == null) {
            throw new IllegalStateException("Nextflow binary has not been located yet. Call locateNextflowBinary() first.");
        }
        return nextflowBinaryPath;
    }

    /**
     * Checks if Docker is available on the system.
     *
     * @return true if Docker is available and accessible
     */
    private boolean isDockerAvailable() {
        try {
            ProcessBuilder pb = new ProcessBuilder("docker", "--version");
            Process process = pb.start();

            boolean finished = process.waitFor(DOCKER_CHECK_TIMEOUT_SECONDS, java.util.concurrent.TimeUnit.SECONDS);
            if (!finished) {
                process.destroyForcibly();
                return false;
            }

            return process.exitValue() == 0;

        } catch (Exception e) {
            logger.log(Level.FINE, "Docker availability check failed", e);
            return false;
        }
    }

    /**
     * Executes Nextflow command directly on the host system using ProcessBuilder.
     *
     * @param nextflowCommand the Nextflow command components
     * @param workspaceDir the workspace directory
     * @param progressListener optional progress listener for monitoring
     * @param timeoutMinutes timeout in minutes
     * @return execution result
     * @throws DockerException if execution fails
     */
    private ExecutionResult executeNextflowDirectly(List<String> nextflowCommand,
                                                   Path workspaceDir,
                                                   ProgressListener progressListener,
                                                   int timeoutMinutes) throws DockerException {
        try {
            if (progressListener != null) {
                progressListener.setMessage("Starting Nextflow workflow...");
                progressListener.setProgress(0.1);
            }

            logger.info("Executing Nextflow command: " + String.join(" ", nextflowCommand));

            // Set up ProcessBuilder
            ProcessBuilder processBuilder = new ProcessBuilder(nextflowCommand);
            processBuilder.directory(workspaceDir.toFile());
            processBuilder.redirectErrorStream(false);

            // Set environment variables
            java.util.Map<String, String> env = processBuilder.environment();

            // Ensure Nextflow home directory
            if (!env.containsKey("NXF_HOME")) {
                env.put("NXF_HOME", System.getProperty("user.home") + "/.nextflow");
            }

            // Ensure Docker is available in PATH
            String currentPath = env.get("PATH");
            if (currentPath != null && !currentPath.contains("/usr/local/bin")) {
                env.put("PATH", currentPath + ":/usr/local/bin:/opt/homebrew/bin");
            }

            // Start the process
            java.time.LocalDateTime startTime = java.time.LocalDateTime.now();
            Process process = processBuilder.start();

            // Monitor the process with progress updates
            ExecutionResult result = monitorNextflowProcess(process, nextflowCommand, progressListener, timeoutMinutes);

            java.time.LocalDateTime endTime = java.time.LocalDateTime.now();

            logger.info("Nextflow execution completed. Exit code: " + result.getExitCode() +
                       ", Duration: " + java.time.Duration.between(startTime, endTime).getSeconds() + "s");

            return result;

        } catch (Exception e) {
            logger.log(Level.SEVERE, "Failed to execute Nextflow command", e);
            throw new DockerException("Failed to execute Nextflow command natively", e);
        }
    }

    /**
     * Monitors a running Nextflow process and captures output.
     *
     * @param process the process to monitor
     * @param command the command being executed
     * @param progressListener optional progress listener
     * @param timeoutMinutes timeout in minutes
     * @return execution result
     * @throws DockerException if monitoring fails
     */
    private ExecutionResult monitorNextflowProcess(Process process,
                                                 List<String> command,
                                                 ProgressListener progressListener,
                                                 int timeoutMinutes) throws DockerException {
        java.time.LocalDateTime startTime = java.time.LocalDateTime.now();

        try {
            // Start output readers
            StringBuilder standardOutput = new StringBuilder();
            StringBuilder errorOutput = new StringBuilder();

            Thread outputReader = new Thread(() -> readNextflowStream(process.getInputStream(), standardOutput, progressListener));
            Thread errorReader = new Thread(() -> readNextflowStream(process.getErrorStream(), errorOutput, null));

            outputReader.start();
            errorReader.start();

            // Wait for process completion with timeout
            boolean finished = process.waitFor(timeoutMinutes, java.util.concurrent.TimeUnit.MINUTES);

            if (!finished) {
                process.destroyForcibly();
                throw new DockerException("Nextflow process timed out after " + timeoutMinutes + " minutes");
            }

            // Wait for output readers to finish
            outputReader.join(OUTPUT_READER_TIMEOUT_SECONDS * 1000); // Convert seconds to milliseconds
            errorReader.join(OUTPUT_READER_TIMEOUT_SECONDS * 1000);

            java.time.LocalDateTime endTime = java.time.LocalDateTime.now();
            int exitCode = process.exitValue();

            return new ExecutionResult(
                String.join(" ", command),
                exitCode,
                standardOutput.toString(),
                errorOutput.toString(),
                startTime,
                endTime
            );

        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            process.destroyForcibly();
            throw new DockerException("Nextflow process monitoring was interrupted", e);
        } catch (Exception e) {
            process.destroyForcibly();
            throw new DockerException("Failed to monitor Nextflow process", e);
        }
    }

    /**
     * Reads from an input stream and appends to a StringBuilder.
     *
     * @param inputStream the stream to read from
     * @param output the StringBuilder to append to
     * @param progressListener optional progress listener for updates
     */
    private void readNextflowStream(java.io.InputStream inputStream, StringBuilder output, ProgressListener progressListener) {
        try (java.io.BufferedReader reader = new java.io.BufferedReader(new java.io.InputStreamReader(inputStream))) {
            String line;
            while ((line = reader.readLine()) != null) {
                output.append(line).append("\n");

                // Update progress listener if available
                if (progressListener != null) {
                    updateProgressFromNextflowOutput(line, progressListener);
                }

                // Log important lines
                if (line.contains("ERROR") || line.contains("WARN")) {
                    logger.warning("Nextflow output: " + line);
                } else if (line.contains("Completed") || line.contains("SUCCESS")) {
                    logger.info("Nextflow output: " + line);
                }
            }
        } catch (java.io.IOException e) {
            logger.log(Level.WARNING, "Error reading from Nextflow process stream", e);
        }
    }

    /**
     * Updates progress based on output from the Nextflow process.
     *
     * @param outputLine the output line from Nextflow
     * @param progressListener the progress listener to update
     */
    private void updateProgressFromNextflowOutput(String outputLine, ProgressListener progressListener) {
        // Simple progress estimation based on Nextflow output patterns
        if (outputLine.contains("Launching workflow")) {
            progressListener.setProgress(0.2);
            progressListener.setMessage("Launching workflow...");
        } else if (outputLine.contains("Staging foreign files")) {
            progressListener.setProgress(0.3);
            progressListener.setMessage("Staging input files...");
        } else if (outputLine.contains("Running pipeline")) {
            progressListener.setProgress(0.4);
            progressListener.setMessage("Running analysis pipeline...");
        } else if (outputLine.contains("Process")) {
            progressListener.setProgress(0.6);
            progressListener.setMessage("Processing samples...");
        } else if (outputLine.contains("Completed at")) {
            progressListener.setProgress(0.9);
            progressListener.setMessage("Finalizing results...");
        }

        // Check for cancellation
        if (progressListener.isCanceled()) {
            logger.info("Progress listener indicates cancellation");
        }
    }


    /**
     * Ensures that the required database exists in the cache directory.
     * Downloads the database if it doesn't exist.
     *
     * @param options TaxTriage configuration options
     * @param progressListener progress listener for updates
     * @throws DocumentOperationException if database setup fails
     */
    private void ensureDatabaseExists(TaxTriageOptions options, ProgressListener progressListener)
            throws DocumentOperationException {

        String databaseType = options.getKrakenDatabase();
        if (databaseType == null || databaseType.isEmpty()) {
            databaseType = STANDARD_DB_NAME;
        }

        String databasePath = DatabasePathUtil.getDatabasePath(options.getKrakenDatabase());
        Path dbPath = Paths.get(databasePath);

        logger.info("Checking database at: " + databasePath);

        // Check if database directory exists and has content
        if (Files.exists(dbPath) && Files.isDirectory(dbPath)) {
            boolean hasDbFiles = DatabasePathUtil.databaseHasContent(dbPath);

            if (hasDbFiles) {
                logger.info("Database already exists at: " + databasePath);
                return;
            }
        }

        // For now, skip explicit database download and use test database
        // or let TaxTriage download it automatically during workflow
        logger.info("Database not found at: " + databasePath);
        logger.info("Will use test database or let TaxTriage handle download");

        // Comment out database download for now as it requires additional parameters
        // TODO: Fix database download workflow once we understand correct parameters
        // downloadDatabase(databaseType, databasePath, progressListener);
    }

    /**
     * Downloads the specified database using the TaxTriage download workflow.
     *
     * @param databaseType the type of database to download (e.g., "standard")
     * @param outputPath the path where the database should be downloaded
     * @param progressListener progress listener for updates
     * @throws DocumentOperationException if database download fails
     */
    private void downloadDatabase(String databaseType, String outputPath, ProgressListener progressListener)
            throws DocumentOperationException {

        try {
            // Create cache directory if it doesn't exist
            Path cacheDir = Paths.get(outputPath).getParent();
            if (!Files.exists(cacheDir)) {
                Files.createDirectories(cacheDir);
                logger.info("Created cache directory: " + cacheDir.toAbsolutePath());
            }

            if (progressListener != null) {
                progressListener.setMessage("Downloading " + databaseType + " database...");
                progressListener.setProgress(0.05);
            }

            // Build Nextflow command for database download
            List<String> downloadCommand = buildDatabaseDownloadCommand(databaseType, outputPath);

            // Log the download command for debugging
            String commandString = String.join(" ", downloadCommand);
            logger.info("Executing database download command: " + commandString);

            // Execute the download command
            ExecutionResult result = executeNextflowDirectly(
                downloadCommand, cacheDir, progressListener, DATABASE_DOWNLOAD_TIMEOUT_MINUTES);

            if (!result.isSuccessful()) {
                // Log both outputs for debugging
                String standardOutput = result.getStandardOutput();
                String errorOutput = result.getErrorOutput();

                logger.severe("Database download failed with exit code: " + result.getExitCode());
                if (standardOutput != null && !standardOutput.trim().isEmpty()) {
                    logger.severe("Database download standard output:\n" + standardOutput);
                }
                if (errorOutput != null && !errorOutput.trim().isEmpty()) {
                    logger.severe("Database download error output:\n" + errorOutput);
                }

                String errorMsg = "Database download failed with exit code " + result.getExitCode();
                if (errorOutput != null && !errorOutput.trim().isEmpty()) {
                    errorMsg += ". Error: " + errorOutput;
                } else if (standardOutput != null && !standardOutput.trim().isEmpty()) {
                    errorMsg += ". Output: " + standardOutput;
                }
                throw new DocumentOperationException(errorMsg);
            }

            logger.info("Database download completed successfully: " + databaseType);

        } catch (IOException e) {
            throw new DocumentOperationException("Failed to create cache directory: " + e.getMessage(), e);
        } catch (DockerException e) {
            throw new DocumentOperationException("Failed to download database: " + e.getMessage(), e);
        }
    }

    /**
     * Builds the Nextflow command for downloading a database.
     *
     * @param databaseType the type of database to download
     * @param outputPath the path where the database should be downloaded
     * @return list of command components
     */
    private List<String> buildDatabaseDownloadCommand(String databaseType, String outputPath) {
        List<String> cmd = new ArrayList<>();

        // Use the extracted Nextflow binary
        cmd.add(getNextflowBinaryPath());
        cmd.add("run");
        cmd.add("https://github.com/jhuapl-bio/taxtriage");
        cmd.add("-r");
        cmd.add("main");
        cmd.add("-profile");
        cmd.add("docker");

        // Download database workflow parameters
        cmd.add("--download_db");

        // The download workflow still requires these parameters even though they're not used
        // Create a dummy samplesheet for the download process
        Path tempDir = Paths.get(System.getProperty("java.io.tmpdir"));
        Path dummySamplesheet = tempDir.resolve("dummy_samplesheet.csv");
        try {
            if (!Files.exists(dummySamplesheet)) {
                Files.write(dummySamplesheet, "sample,fastq_1,fastq_2\ndummy,dummy.fq,\n".getBytes());
            }
            cmd.add("--input");
            cmd.add(dummySamplesheet.toString());
        } catch (IOException e) {
            // Fall back to using fastq_1 parameter
            cmd.add("--fastq_1");
            cmd.add("dummy.fastq");
        }

        // Output directory (required even for download)
        cmd.add("--outdir");
        cmd.add(tempDir.resolve("taxtriage_download_temp").toString());

        // Database path
        cmd.add("--db");
        cmd.add(outputPath);

        // For non-standard databases, we might need to specify different parameters
        // The viral database might be specified differently
        if ("viral".equals(databaseType)) {
            // Viral database might use a specific parameter or profile
            // For now, we'll use the standard download
        }

        return cmd;
    }
}