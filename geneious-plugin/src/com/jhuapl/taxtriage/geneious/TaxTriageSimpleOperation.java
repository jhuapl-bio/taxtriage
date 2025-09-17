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
import com.jhuapl.taxtriage.geneious.docker.DockerException;
import com.jhuapl.taxtriage.geneious.docker.DockerManager;
import com.jhuapl.taxtriage.geneious.docker.ExecutionResult;
import com.jhuapl.taxtriage.geneious.importer.ImportResult;
import com.jhuapl.taxtriage.geneious.importer.ResultImporter;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Simplified TaxTriage operation that appears in the Tools menu.
 *
 * This class implements the complete TaxTriage workflow execution including:
 * - File selection and validation
 * - Docker environment setup
 * - Nextflow samplesheet and config generation
 * - TaxTriage workflow execution via Docker
 * - Result import back to Geneious
 */
public class TaxTriageSimpleOperation extends DocumentOperation {

    private static final Logger logger = Logger.getLogger(TaxTriageSimpleOperation.class.getName());

    // Nextflow binary relative path within the plugin directory
    private static final String NEXTFLOW_RELATIVE_PATH = "bin/nextflow";

    // Cache for the located Nextflow binary path
    private static String nextflowBinaryPath = null;

    // Timeout for workflow execution (2 hours)
    private static final int WORKFLOW_TIMEOUT_MINUTES = 120;

    // Database cache directory
    private static final String CACHE_DIR_NAME = ".taxtriage-geneious";
    private static final String STANDARD_DB_NAME = "standard";

    // Database download timeout (3 hours)
    private static final int DATABASE_DOWNLOAD_TIMEOUT_MINUTES = 180;

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
        System.out.println("[TaxTriage] performOperation called");
        System.out.println("[TaxTriage] Number of selected documents: " + (documents != null ? documents.length : 0));

        if (progressListener != null) {
            progressListener.setMessage("Starting TaxTriage analysis...");
            progressListener.setProgress(0.1);
        }

        List<AnnotatedPluginDocument> results = new ArrayList<>();

        // Cast to TaxTriageOptions to access file browser selections
        TaxTriageOptions taxTriageOptions = (TaxTriageOptions) options;

        // Get input files from browsers
        List<File> inputFiles = taxTriageOptions.getInputFiles();
        File inputDirectory = taxTriageOptions.getInputDirectory();
        File outputDirectory = taxTriageOptions.getOutputDirectory();

        // Log the file selections
        System.out.println("[TaxTriage] Input files from browser: " + (inputFiles != null ? inputFiles.size() : 0));
        if (inputFiles != null) {
            for (File f : inputFiles) {
                System.out.println("[TaxTriage]   - File: " + f.getAbsolutePath());
            }
        }
        System.out.println("[TaxTriage] Input directory: " + (inputDirectory != null ? inputDirectory.getAbsolutePath() : "null"));
        System.out.println("[TaxTriage] Output directory: " + (outputDirectory != null ? outputDirectory.getAbsolutePath() : "null"));

        // Check if we have input files from browsers
        boolean hasFileInputs = (inputFiles != null && !inputFiles.isEmpty()) || inputDirectory != null;

        System.out.println("[TaxTriage] Has file inputs from browser: " + hasFileInputs);

        if (!hasFileInputs) {
            // Only check selected documents if no file browser input was provided
            boolean hasSelectedDocuments = documents != null && documents.length > 0;
            System.out.println("[TaxTriage] Has selected documents: " + hasSelectedDocuments);

            if (!hasSelectedDocuments) {
                throw new DocumentOperationException("No input provided. Please use the file browsers to specify input files or directory.");
            }

            // If using selected documents (backward compatibility), validate they are sequence documents
            for (AnnotatedPluginDocument doc : documents) {
                System.out.println("[TaxTriage] Checking document: " + doc.getName() + " - Type: " + doc.getDocument().getClass().getName());
                if (!(doc.getDocument() instanceof SequenceDocument)) {
                    throw new DocumentOperationException("Document '" + doc.getName() + "' is not a sequence document. Please use the file browsers instead.");
                }
            }
        } else {
            // When using file browsers, ignore selected documents completely
            System.out.println("[TaxTriage] Using file browser inputs, ignoring any selected documents");
        }

        if (progressListener != null) {
            String inputDescription = hasFileInputs ? "file browser inputs" : documents.length + " selected document(s)";
            progressListener.setMessage("Processing " + inputDescription + "...");
            progressListener.setProgress(0.5);
        }

        // Log configuration
        System.out.println("[TaxTriage] Configuration:");
        System.out.println("[TaxTriage]   Preset: " + taxTriageOptions.getSequencingPreset());
        System.out.println("[TaxTriage]   Kraken DB: " + taxTriageOptions.getKrakenDatabase());
        System.out.println("[TaxTriage]   Bracken DB: " + taxTriageOptions.getBrackenDatabase());
        System.out.println("[TaxTriage]   Threads: " + taxTriageOptions.getThreadCount());
        System.out.println("[TaxTriage]   Memory: " + taxTriageOptions.getMemoryLimit() + " GB");

        try {
            // Locate bundled Nextflow binary before workflow execution
            if (progressListener != null) {
                progressListener.setMessage("Preparing Nextflow environment...");
                progressListener.setProgress(0.05);
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
                progressListener.setProgress(0.15);
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
                    System.err.println("[TaxTriage] Standard output:\n" + standardOutput);
                }
                if (errorOutput != null && !errorOutput.trim().isEmpty()) {
                    logger.severe("Error output:\n" + errorOutput);
                    System.err.println("[TaxTriage] Error output:\n" + errorOutput);
                } else {
                    logger.severe("No error output captured. The workflow may have failed silently.");
                    System.err.println("[TaxTriage] No error output captured.");
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

            // Step 7: Import results back to Geneious
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
        if (!file.exists() || !file.isFile()) {
            return false;
        }

        String name = file.getName().toLowerCase();

        // Only accept FASTQ files, not FASTA files
        // Also filter out specific files that shouldn't be included
        if (name.contains(".fastp.") || name.contains(".dwnld.") || name.contains("references")) {
            return false;
        }

        return name.endsWith(".fastq") || name.endsWith(".fq") ||
               name.endsWith(".fastq.gz") || name.endsWith(".fq.gz");
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

        // Generate samplesheet with local filesystem paths
        SampleSheetBuilder sampleSheetBuilder = new SampleSheetBuilder();
        File sampleSheetFile = workspaceDir.resolve("config").resolve("samplesheet.csv").toFile();

        // Create list of files with workspace input paths for samplesheet
        // Use a Set to track processed filenames and avoid duplicates
        Set<String> processedFilenames = new HashSet<>();
        List<File> workspaceInputFiles = new ArrayList<>();

        logger.info("Processing " + inputFiles.size() + " input files for samplesheet generation");

        for (File inputFile : inputFiles) {
            String filename = inputFile.getName();

            // Check if we've already processed this filename
            if (processedFilenames.contains(filename)) {
                logger.warning("Skipping duplicate filename: " + filename + " (path: " + inputFile.getAbsolutePath() + ")");
                continue;
            }

            // The files are copied to the workspace input directory
            File workspaceFile = workspaceDir.resolve("input").resolve(filename).toFile();
            workspaceInputFiles.add(workspaceFile);
            processedFilenames.add(filename);

            logger.info("Added file to samplesheet: " + filename + " (source: " + inputFile.getAbsolutePath() + ")");
        }

        logger.info("Final count: " + workspaceInputFiles.size() + " unique files will be added to samplesheet");

        sampleSheetBuilder.buildSampleSheet(workspaceInputFiles, sampleSheetFile,
            options.getSequencingPreset().name());
        logger.info("Generated samplesheet: " + sampleSheetFile.getAbsolutePath());

        // Copy input files to workspace (or create symlinks)
        copyInputFilesToWorkspace(workspaceDir, inputFiles);
    }

    /**
     * Copies input files to the workspace input directory.
     */
    private void copyInputFilesToWorkspace(Path workspaceDir, List<File> inputFiles) throws IOException {
        Path inputDir = workspaceDir.resolve("input");

        for (File inputFile : inputFiles) {
            Path targetPath = inputDir.resolve(inputFile.getName());

            // Copy file to workspace, replacing if it already exists
            Files.copy(inputFile.toPath(), targetPath,
                      java.nio.file.StandardCopyOption.REPLACE_EXISTING);
            logger.fine("Copied input file: " + inputFile.getName());
        }

        logger.info("Copied " + inputFiles.size() + " input files to workspace");
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
        logger.info("Executing command: " + commandString);
        System.out.println("[TaxTriage] Executing Nextflow command: " + commandString);
        System.out.println("[TaxTriage] Working directory: " + workspaceDir.toAbsolutePath());

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

        // Database configuration
        // Use database name (e.g., "viral", "standard-8", "flukraken2")
        cmd.add("--db");
        cmd.add(options.getKrakenDatabase());

        // Enable automatic database download
        cmd.add("--download_db");
        cmd.add("true");

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

        } catch (Exception e) {
            logger.log(Level.WARNING, "Error during result import", e);
        }

        return importedDocs;
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

            boolean finished = process.waitFor(10, java.util.concurrent.TimeUnit.SECONDS);
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
            outputReader.join(5000); // 5 second timeout for output reading
            errorReader.join(5000);

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
     * Gets the path to the cached database directory for the specified database type.
     *
     * @param options TaxTriage configuration options
     * @return absolute path to the database directory
     */
    private String getDatabasePath(TaxTriageOptions options) {
        String userHome = System.getProperty("user.home");
        Path cacheDir = Paths.get(userHome, CACHE_DIR_NAME);

        String databaseType = options.getKrakenDatabase();
        if (databaseType == null || databaseType.isEmpty()) {
            databaseType = STANDARD_DB_NAME;
        }

        Path databasePath = cacheDir.resolve(databaseType);
        return databasePath.toAbsolutePath().toString();
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

        String databasePath = getDatabasePath(options);
        Path dbPath = Paths.get(databasePath);

        logger.info("Checking database at: " + databasePath);

        // Check if database directory exists and has content
        if (Files.exists(dbPath) && Files.isDirectory(dbPath)) {
            try {
                // Check if directory has database files (any .k2d files indicate a valid Kraken database)
                boolean hasDbFiles = Files.walk(dbPath)
                    .anyMatch(path -> path.toString().endsWith(".k2d") ||
                                    path.toString().endsWith(".kmer_distrib") ||
                                    path.toString().contains("taxonomy"));

                if (hasDbFiles) {
                    logger.info("Database already exists at: " + databasePath);
                    return;
                }
            } catch (IOException e) {
                logger.warning("Error checking existing database: " + e.getMessage());
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
            System.out.println("[TaxTriage] Downloading database with command: " + commandString);

            // Execute the download command
            ExecutionResult result = executeNextflowDirectly(
                downloadCommand, cacheDir, progressListener, DATABASE_DOWNLOAD_TIMEOUT_MINUTES);

            if (!result.isSuccessful()) {
                // Log both outputs for debugging
                String standardOutput = result.getStandardOutput();
                String errorOutput = result.getErrorOutput();

                logger.severe("Database download failed with exit code: " + result.getExitCode());
                if (standardOutput != null && !standardOutput.trim().isEmpty()) {
                    logger.severe("Standard output:\n" + standardOutput);
                    System.err.println("[TaxTriage] Database download standard output:\n" + standardOutput);
                }
                if (errorOutput != null && !errorOutput.trim().isEmpty()) {
                    logger.severe("Error output:\n" + errorOutput);
                    System.err.println("[TaxTriage] Database download error output:\n" + errorOutput);
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