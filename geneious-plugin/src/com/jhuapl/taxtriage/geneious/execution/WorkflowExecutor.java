package com.jhuapl.taxtriage.geneious.execution;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.sequence.SequenceDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentOperationException;
import com.jhuapl.taxtriage.geneious.config.ConfigGenerator;
import com.jhuapl.taxtriage.geneious.config.SampleSheetBuilder;
import com.jhuapl.taxtriage.geneious.config.TaxTriageConfig;
import com.jhuapl.taxtriage.geneious.docker.DockerException;
import com.jhuapl.taxtriage.geneious.docker.DockerManager;
import com.jhuapl.taxtriage.geneious.docker.ExecutionResult;
import jebl.util.ProgressListener;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Orchestrates the complete TaxTriage workflow execution.
 *
 * This class manages the end-to-end execution of TaxTriage workflows, including:
 * - Exporting Geneious documents to FASTQ files
 * - Generating configuration and samplesheet files
 * - Executing Nextflow workflows via Docker
 * - Monitoring progress and handling errors
 * - Cleaning up temporary files
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class WorkflowExecutor {

    private static final Logger logger = Logger.getLogger(WorkflowExecutor.class.getName());

    /** Default timeout for workflow execution in minutes */
    private static final int DEFAULT_TIMEOUT_MINUTES = 120;

    /** Progress weights for different execution phases */
    private static final double EXPORT_WEIGHT = 0.15;
    private static final double CONFIG_WEIGHT = 0.05;
    private static final double SAMPLESHEET_WEIGHT = 0.05;
    private static final double EXECUTION_WEIGHT = 0.70;
    private static final double IMPORT_WEIGHT = 0.05;

    private final DockerManager dockerManager;
    private final ConfigGenerator configGenerator;
    private final SampleSheetBuilder sampleSheetBuilder;
    private final SequenceExporter sequenceExporter;

    /**
     * Creates a new WorkflowExecutor with default components.
     *
     * @throws DockerException if Docker is not available
     */
    public WorkflowExecutor() throws DockerException {
        this(new DockerManager(), new ConfigGenerator(), new SampleSheetBuilder(), new SequenceExporter());
    }

    /**
     * Creates a new WorkflowExecutor with custom components for dependency injection.
     *
     * @param dockerManager manages Docker operations
     * @param configGenerator generates Nextflow configuration
     * @param sampleSheetBuilder builds sample sheets
     * @param sequenceExporter exports Geneious sequences to files
     */
    public WorkflowExecutor(DockerManager dockerManager,
                           ConfigGenerator configGenerator,
                           SampleSheetBuilder sampleSheetBuilder,
                           SequenceExporter sequenceExporter) {
        this.dockerManager = dockerManager;
        this.configGenerator = configGenerator;
        this.sampleSheetBuilder = sampleSheetBuilder;
        this.sequenceExporter = sequenceExporter;
    }

    /**
     * Executes the complete TaxTriage workflow.
     *
     * @param context the workflow context containing all necessary information
     * @param progressListener optional progress listener for tracking execution
     * @return CompletableFuture that completes when the workflow finishes
     */
    public CompletableFuture<WorkflowContext> executeWorkflow(WorkflowContext context,
                                                             ProgressListener progressListener) {
        return CompletableFuture.supplyAsync(() -> {
            try {
                context.setStartedAt(LocalDateTime.now());
                updateProgress(context, progressListener, WorkflowContext.ExecutionState.INITIALIZED,
                              "Starting TaxTriage workflow...", 0.0);

                // Phase 1: Export sequences
                exportSequences(context, progressListener);

                // Phase 2: Generate configuration
                generateConfiguration(context, progressListener);

                // Phase 3: Create samplesheet
                createSampleSheet(context, progressListener);

                // Phase 4: Execute workflow
                executeNextflowWorkflow(context, progressListener);

                // Phase 5: Import results
                importResults(context, progressListener);

                // Mark as completed
                context.markCompleted();
                if (progressListener != null) {
                    progressListener.setMessage("Workflow completed successfully");
                    progressListener.setProgress(1.0);
                }

                logger.info("TaxTriage workflow completed successfully: " + context.getWorkflowId());
                return context;

            } catch (Exception e) {
                logger.log(Level.SEVERE, "Workflow execution failed", e);
                context.markFailed(e, "Workflow execution failed: " + e.getMessage());

                if (progressListener != null) {
                    progressListener.setMessage("Workflow failed: " + e.getMessage());
                }

                // Clean up on failure
                try {
                    cleanupTemporaryFiles(context);
                } catch (Exception cleanupError) {
                    logger.log(Level.WARNING, "Failed to clean up temporary files", cleanupError);
                }

                return context;
            }
        });
    }

    /**
     * Cancels a running workflow.
     *
     * @param context the workflow context to cancel
     */
    public void cancelWorkflow(WorkflowContext context) {
        context.markCancelled();
        logger.info("Workflow cancellation requested: " + context.getWorkflowId());

        // Note: Actual process cancellation would need to be implemented
        // by tracking running processes and calling their destroy methods
    }

    /**
     * Exports Geneious sequence documents to FASTQ files.
     *
     * @param context the workflow context
     * @param progressListener optional progress listener
     * @throws DocumentOperationException if export fails
     */
    private void exportSequences(WorkflowContext context, ProgressListener progressListener)
            throws DocumentOperationException {
        updateProgress(context, progressListener, WorkflowContext.ExecutionState.EXPORTING_SEQUENCES,
                      "Exporting sequences to files...", 0.0);

        try {
            // Create input directory
            Files.createDirectories(context.getInputDirectory());

            // Export each document
            List<File> exportedFiles = new ArrayList<>();
            List<AnnotatedPluginDocument> inputDocs = context.getInputDocuments();

            for (int i = 0; i < inputDocs.size(); i++) {
                AnnotatedPluginDocument doc = inputDocs.get(i);

                if (!(doc.getDocument() instanceof SequenceDocument)) {
                    throw new DocumentOperationException("Document " + doc.getName() + " is not a sequence document");
                }

                SequenceDocument seqDoc = (SequenceDocument) doc.getDocument();
                File outputFile = context.getInputDirectory().resolve(sanitizeFileName(doc.getName()) + ".fastq").toFile();

                sequenceExporter.exportToFastq(seqDoc, outputFile);
                exportedFiles.add(outputFile);

                // Update progress
                double phaseProgress = (double) (i + 1) / inputDocs.size();
                double totalProgress = EXPORT_WEIGHT * phaseProgress;
                updateProgress(context, progressListener, WorkflowContext.ExecutionState.EXPORTING_SEQUENCES,
                              "Exported " + (i + 1) + "/" + inputDocs.size() + " sequences", totalProgress);
            }

            context.setExportedFiles(exportedFiles);
            logger.info("Exported " + exportedFiles.size() + " sequence files");

        } catch (Exception e) {
            throw new DocumentOperationException("Failed to export sequences", e);
        }
    }

    /**
     * Generates Nextflow configuration files.
     *
     * @param context the workflow context
     * @param progressListener optional progress listener
     * @throws DocumentOperationException if configuration generation fails
     */
    private void generateConfiguration(WorkflowContext context, ProgressListener progressListener)
            throws DocumentOperationException {
        updateProgress(context, progressListener, WorkflowContext.ExecutionState.GENERATING_CONFIG,
                      "Generating workflow configuration...", EXPORT_WEIGHT);

        try {
            // Create config directory
            Files.createDirectories(context.getConfigDirectory());

            TaxTriageConfig config = context.getConfig();

            // Generate nextflow.config
            File configFile = context.getConfigDirectory().resolve("nextflow.config").toFile();
            configGenerator.generateConfig(config, configFile);
            context.setConfigFile(configFile);

            // Generate params.json
            File paramsFile = context.getConfigDirectory().resolve("params.json").toFile();
            configGenerator.generateParams(config, paramsFile);
            context.setParamsFile(paramsFile);

            logger.info("Generated configuration files: " + configFile.getName() + ", " + paramsFile.getName());

        } catch (Exception e) {
            throw new DocumentOperationException("Failed to generate configuration", e);
        }
    }

    /**
     * Creates the samplesheet file.
     *
     * @param context the workflow context
     * @param progressListener optional progress listener
     * @throws DocumentOperationException if samplesheet creation fails
     */
    private void createSampleSheet(WorkflowContext context, ProgressListener progressListener)
            throws DocumentOperationException {
        updateProgress(context, progressListener, WorkflowContext.ExecutionState.CREATING_SAMPLESHEET,
                      "Creating samplesheet...", EXPORT_WEIGHT + CONFIG_WEIGHT);

        try {
            File sampleSheetFile = context.getConfigDirectory().resolve("samplesheet.csv").toFile();
            String preset = context.getOptions().getSequencingPreset().name();

            sampleSheetBuilder.buildSampleSheet(context.getExportedFiles(), sampleSheetFile, preset);
            context.setSampleSheetFile(sampleSheetFile);

            logger.info("Created samplesheet: " + sampleSheetFile.getName());

        } catch (Exception e) {
            throw new DocumentOperationException("Failed to create samplesheet", e);
        }
    }

    /**
     * Executes the Nextflow workflow in Docker.
     *
     * @param context the workflow context
     * @param progressListener optional progress listener
     * @throws DocumentOperationException if workflow execution fails
     */
    private void executeNextflowWorkflow(WorkflowContext context, ProgressListener progressListener)
            throws DocumentOperationException {
        updateProgress(context, progressListener, WorkflowContext.ExecutionState.EXECUTING_WORKFLOW,
                      "Executing TaxTriage workflow...", EXPORT_WEIGHT + CONFIG_WEIGHT + SAMPLESHEET_WEIGHT);

        try {
            // Create output directory
            Files.createDirectories(context.getOutputDirectory());

            // Build Nextflow command
            String nextflowCommand = buildNextflowCommand(context);

            logger.info("Executing Nextflow command: " + nextflowCommand);

            // Create a delegating progress listener that updates context and preserves detailed messages
            ProgressListener delegatingListener = new ProgressListener() {
                @Override
                protected void _setProgress(double progress) {
                    if (progressListener != null) {
                        double adjustedProgress = EXPORT_WEIGHT + CONFIG_WEIGHT + SAMPLESHEET_WEIGHT +
                                                (EXECUTION_WEIGHT * progress);
                        progressListener.setProgress(adjustedProgress);
                    }
                    context.setProgress(EXPORT_WEIGHT + CONFIG_WEIGHT + SAMPLESHEET_WEIGHT +
                                      (EXECUTION_WEIGHT * progress));
                }

                @Override
                protected void _setMessage(String message) {
                    // Pass through the detailed message from ExecutionMonitor
                    if (progressListener != null) {
                        progressListener.setMessage(message);
                    }
                    context.setCurrentMessage(message);

                    // Log process changes for debugging
                    if (message != null && message.contains("TaxTriage:")) {
                        logger.fine("Workflow progress: " + message);
                    }
                }

                @Override
                protected void _setIndeterminateProgress() {
                    if (progressListener != null) {
                        progressListener.setIndeterminateProgress();
                    }
                }

                @Override
                public boolean isCanceled() {
                    return progressListener != null && progressListener.isCanceled();
                }
            };

            updateProgress(context, progressListener, WorkflowContext.ExecutionState.MONITORING_PROGRESS,
                          "Starting workflow execution...", EXPORT_WEIGHT + CONFIG_WEIGHT + SAMPLESHEET_WEIGHT);

            // Execute the workflow
            ExecutionResult result = dockerManager.executeNextflowCommand(
                    nextflowCommand,
                    context.getInputDirectory(),
                    context.getOutputDirectory(),
                    context.getWorkingDirectory(),
                    delegatingListener
            );

            if (result.getExitCode() != 0) {
                throw new DocumentOperationException("Nextflow execution failed with exit code " +
                        result.getExitCode() + ": " + result.getErrorOutput());
            }

            logger.info("Nextflow workflow completed successfully");

        } catch (Exception e) {
            throw new DocumentOperationException("Failed to execute workflow", e);
        }
    }

    /**
     * Imports and processes workflow results.
     *
     * @param context the workflow context
     * @param progressListener optional progress listener
     * @throws DocumentOperationException if result import fails
     */
    private void importResults(WorkflowContext context, ProgressListener progressListener)
            throws DocumentOperationException {
        updateProgress(context, progressListener, WorkflowContext.ExecutionState.IMPORTING_RESULTS,
                      "Processing workflow results...", EXPORT_WEIGHT + CONFIG_WEIGHT + SAMPLESHEET_WEIGHT + EXECUTION_WEIGHT);

        try {
            Path outputDir = context.getOutputDirectory();
            List<File> outputFiles = new ArrayList<>();

            // Collect output files
            if (Files.exists(outputDir)) {
                Files.walk(outputDir)
                     .filter(Files::isRegularFile)
                     .forEach(path -> outputFiles.add(path.toFile()));
            }

            context.setOutputFiles(outputFiles);

            // Look for specific result files
            findResultFiles(context, outputDir);

            logger.info("Processed " + outputFiles.size() + " output files");

        } catch (Exception e) {
            throw new DocumentOperationException("Failed to import results", e);
        }
    }

    /**
     * Finds and sets specific result files in the workflow context.
     *
     * @param context the workflow context
     * @param outputDir the output directory to search
     */
    private void findResultFiles(WorkflowContext context, Path outputDir) {
        try {
            // Look for common TaxTriage output files
            Files.walk(outputDir)
                 .filter(Files::isRegularFile)
                 .forEach(path -> {
                     String fileName = path.getFileName().toString().toLowerCase();
                     if (fileName.contains("report") && fileName.endsWith(".html")) {
                         context.setWorkflowReport(path.toFile());
                     } else if (fileName.contains("taxonomy") || fileName.contains("classification")) {
                         context.setTaxonomyResults(path.toFile());
                     } else if (fileName.contains("quality") || fileName.contains("qc")) {
                         context.setQualityReport(path.toFile());
                     }
                 });
        } catch (Exception e) {
            logger.log(Level.WARNING, "Error finding result files", e);
        }
    }

    /**
     * Builds the Nextflow command string.
     *
     * @param context the workflow context
     * @return the Nextflow command string
     */
    private String buildNextflowCommand(WorkflowContext context) {
        StringBuilder cmd = new StringBuilder();
        cmd.append("nextflow run");

        // Add logging options for better output capture
        cmd.append(" -log /data/work/.nextflow.log");
        cmd.append(" -ansi-log false");  // Disable ANSI formatting

        // Enable trace file for progress monitoring
        // This creates output/pipeline_info/execution_trace_*.txt
        cmd.append(" -with-trace");

        // Add TaxTriage workflow (assuming it's bundled in the Docker image)
        cmd.append(" /opt/taxtriage");

        // Add configuration
        cmd.append(" -c ").append(context.getConfigFile().getAbsolutePath());

        // Add parameters
        cmd.append(" -params-file ").append(context.getParamsFile().getAbsolutePath());

        // Add samplesheet
        cmd.append(" --input ").append(context.getSampleSheetFile().getAbsolutePath());

        // Add output directory
        cmd.append(" --outdir ").append(context.getOutputDirectory().toAbsolutePath());

        // Add preset
        cmd.append(" --preset ").append(context.getOptions().getSequencingPreset().name().toLowerCase());

        return cmd.toString();
    }

    /**
     * Sanitizes a filename by removing invalid characters.
     *
     * @param name the original filename
     * @return sanitized filename
     */
    private String sanitizeFileName(String name) {
        if (name == null) {
            return "sequence";
        }
        return name.replaceAll("[^a-zA-Z0-9._-]", "_");
    }

    /**
     * Updates progress and state in both context and progress listener.
     *
     * @param context the workflow context
     * @param progressListener optional progress listener
     * @param state the new execution state
     * @param message the progress message
     * @param progress the progress value (0.0 to 1.0)
     */
    private void updateProgress(WorkflowContext context, ProgressListener progressListener,
                               WorkflowContext.ExecutionState state, String message, double progress) {
        context.updateState(state, message, progress);

        if (progressListener != null) {
            progressListener.setMessage(message);
            progressListener.setProgress(progress);
        }
    }

    /**
     * Cleans up temporary files and directories.
     *
     * @param context the workflow context
     */
    private void cleanupTemporaryFiles(WorkflowContext context) {
        try {
            // Note: In a real implementation, you might want to make cleanup optional
            // based on user preferences or debugging needs
            Path workingDir = context.getWorkingDirectory();
            if (Files.exists(workingDir)) {
                Files.walk(workingDir)
                     .sorted((a, b) -> b.compareTo(a)) // Delete files before directories
                     .forEach(path -> {
                         try {
                             Files.deleteIfExists(path);
                         } catch (Exception e) {
                             logger.log(Level.WARNING, "Failed to delete: " + path, e);
                         }
                     });
            }
            logger.info("Cleaned up temporary files for workflow: " + context.getWorkflowId());
        } catch (Exception e) {
            logger.log(Level.WARNING, "Error during cleanup", e);
        }
    }

    /**
     * Helper class for exporting Geneious sequences to FASTQ files.
     */
    public static class SequenceExporter {

        /**
         * Exports a sequence document to a FASTQ file.
         *
         * @param sequenceDoc the sequence document to export
         * @param outputFile the output FASTQ file
         * @throws IOException if export fails
         */
        public void exportToFastq(SequenceDocument sequenceDoc, File outputFile) throws IOException {
            // Ensure parent directory exists
            File parentDir = outputFile.getParentFile();
            if (parentDir != null && !parentDir.exists()) {
                if (!parentDir.mkdirs()) {
                    throw new IOException("Failed to create directory: " + parentDir.getAbsolutePath());
                }
            }

            try (BufferedWriter writer = Files.newBufferedWriter(outputFile.toPath(), StandardCharsets.UTF_8)) {
                // Write FASTQ format
                writer.write("@" + sequenceDoc.getName());
                writer.newLine();
                writer.write(sequenceDoc.getSequenceString());
                writer.newLine();
                writer.write("+");
                writer.newLine();

                // Generate quality scores (placeholder - real implementation would use actual quality)
                String quality = generateQualityString(sequenceDoc.getSequenceLength());
                writer.write(quality);
                writer.newLine();
            }
        }

        /**
         * Generates a quality string for FASTQ format.
         *
         * @param length the sequence length
         * @return quality string with default high quality scores
         */
        private String generateQualityString(int length) {
            StringBuilder quality = new StringBuilder();
            for (int i = 0; i < length; i++) {
                quality.append("I"); // ASCII 73 = Phred 40 (high quality)
            }
            return quality.toString();
        }
    }
}