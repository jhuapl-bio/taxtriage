package com.jhuapl.taxtriage.geneious.execution;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentOperationException;
import com.jhuapl.taxtriage.geneious.TaxTriageOptions;
import com.jhuapl.taxtriage.geneious.config.TaxTriageConfig;
import com.jhuapl.taxtriage.geneious.docker.DockerException;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * High-level API for running TaxTriage workflows.
 *
 * This class provides a simplified interface for executing TaxTriage workflows,
 * handling the coordination between Geneious documents, workflow execution,
 * and result processing. It integrates with Geneious ProgressListener for
 * user feedback and manages the complete workflow lifecycle.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class TaxTriageWorkflowRunner {

    private static final Logger logger = Logger.getLogger(TaxTriageWorkflowRunner.class.getName());

    /** Default workflow timeout in minutes */
    private static final int DEFAULT_TIMEOUT_MINUTES = 120;

    /** Temporary directory prefix for workflow execution */
    private static final String TEMP_DIR_PREFIX = "taxtriage_workflow_";

    private final WorkflowExecutor workflowExecutor;
    private final Path tempBaseDirectory;

    /**
     * Creates a new TaxTriageWorkflowRunner with default settings.
     *
     * @throws DockerException if Docker is not available
     * @throws IOException if temporary directory creation fails
     */
    public TaxTriageWorkflowRunner() throws DockerException, IOException {
        this(new WorkflowExecutor(), createDefaultTempDirectory());
    }

    /**
     * Creates a new TaxTriageWorkflowRunner with custom components.
     *
     * @param workflowExecutor the workflow executor to use
     * @param tempBaseDirectory the base directory for temporary files
     */
    public TaxTriageWorkflowRunner(WorkflowExecutor workflowExecutor, Path tempBaseDirectory) {
        this.workflowExecutor = workflowExecutor;
        this.tempBaseDirectory = tempBaseDirectory;
    }

    /**
     * Runs a TaxTriage workflow with the given documents and options.
     *
     * @param documents the input Geneious documents
     * @param options the TaxTriage options
     * @param progressListener optional progress listener for tracking
     * @return the workflow context containing execution results
     * @throws DocumentOperationException if the workflow fails
     */
    public WorkflowContext runWorkflow(AnnotatedPluginDocument[] documents,
                                      TaxTriageOptions options,
                                      ProgressListener progressListener) throws DocumentOperationException {
        return runWorkflow(documents, options, progressListener, DEFAULT_TIMEOUT_MINUTES);
    }

    /**
     * Runs a TaxTriage workflow with the given documents, options, and timeout.
     *
     * @param documents the input Geneious documents
     * @param options the TaxTriage options
     * @param progressListener optional progress listener for tracking
     * @param timeoutMinutes the maximum execution time in minutes
     * @return the workflow context containing execution results
     * @throws DocumentOperationException if the workflow fails
     */
    public WorkflowContext runWorkflow(AnnotatedPluginDocument[] documents,
                                      TaxTriageOptions options,
                                      ProgressListener progressListener,
                                      int timeoutMinutes) throws DocumentOperationException {
        // Validate inputs
        validateInputs(documents, options);

        WorkflowContext context = null;
        try {
            // Create workflow context
            context = createWorkflowContext(documents, options);

            if (progressListener != null) {
                progressListener.setMessage("Initializing TaxTriage workflow...");
                progressListener.setProgress(0.0);
            }

            logger.info("Starting TaxTriage workflow: " + context.getWorkflowId());

            // Execute workflow asynchronously
            CompletableFuture<WorkflowContext> workflowFuture = workflowExecutor.executeWorkflow(context, progressListener);

            // Wait for completion with timeout
            WorkflowContext result;
            try {
                result = workflowFuture.get(timeoutMinutes, TimeUnit.MINUTES);
            } catch (TimeoutException e) {
                workflowExecutor.cancelWorkflow(context);
                throw new DocumentOperationException("Workflow execution timed out after " + timeoutMinutes + " minutes");
            } catch (ExecutionException e) {
                Throwable cause = e.getCause();
                if (cause instanceof DocumentOperationException) {
                    throw (DocumentOperationException) cause;
                } else {
                    throw new DocumentOperationException("Workflow execution failed", cause);
                }
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                workflowExecutor.cancelWorkflow(context);
                throw new DocumentOperationException("Workflow execution was interrupted");
            }

            // Check for cancellation by user
            if (progressListener != null && progressListener.isCanceled()) {
                workflowExecutor.cancelWorkflow(result);
                throw new DocumentOperationException("Workflow was cancelled by user");
            }

            // Validate results
            if (!result.isSuccessful()) {
                String errorMessage = "Workflow failed";
                if (result.getLastError() != null) {
                    errorMessage = result.getLastError().getMessage();
                }
                throw new DocumentOperationException(errorMessage, result.getLastError());
            }

            logger.info("TaxTriage workflow completed successfully: " + result.getWorkflowId());
            return result;

        } catch (Exception e) {
            if (context != null) {
                context.markFailed(e, "Workflow execution failed");
            }

            if (e instanceof DocumentOperationException) {
                throw e;
            } else {
                throw new DocumentOperationException("Unexpected error during workflow execution", e);
            }
        }
    }

    /**
     * Runs a workflow asynchronously and returns a future.
     *
     * @param documents the input Geneious documents
     * @param options the TaxTriage options
     * @param progressListener optional progress listener for tracking
     * @return CompletableFuture that will contain the workflow results
     */
    public CompletableFuture<WorkflowContext> runWorkflowAsync(AnnotatedPluginDocument[] documents,
                                                              TaxTriageOptions options,
                                                              ProgressListener progressListener) {
        return CompletableFuture.supplyAsync(() -> {
            try {
                return runWorkflow(documents, options, progressListener);
            } catch (DocumentOperationException e) {
                throw new RuntimeException(e);
            }
        });
    }

    /**
     * Creates a workflow context from input documents and options.
     *
     * @param documents the input documents
     * @param options the TaxTriage options
     * @return the created workflow context
     * @throws DocumentOperationException if context creation fails
     */
    private WorkflowContext createWorkflowContext(AnnotatedPluginDocument[] documents,
                                                 TaxTriageOptions options) throws DocumentOperationException {
        try {
            // Create temporary working directory
            Path workingDirectory = createWorkingDirectory();

            // Convert options to configuration
            TaxTriageConfig config = TaxTriageConfig.fromOptions(options);

            // Validate configuration
            String validationError = config.validate();
            if (validationError != null) {
                throw new DocumentOperationException("Invalid configuration: " + validationError);
            }

            // Create context
            List<AnnotatedPluginDocument> docList = Arrays.asList(documents);
            return new WorkflowContext(docList, options, config, workingDirectory);

        } catch (Exception e) {
            throw new DocumentOperationException("Failed to create workflow context", e);
        }
    }

    /**
     * Creates a unique working directory for the workflow.
     *
     * @return the path to the working directory
     * @throws IOException if directory creation fails
     */
    private Path createWorkingDirectory() throws IOException {
        // Create timestamp-based directory name
        String timestamp = LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss_SSS"));
        String dirName = TEMP_DIR_PREFIX + timestamp;

        Path workingDir = tempBaseDirectory.resolve(dirName);
        Files.createDirectories(workingDir);

        logger.info("Created working directory: " + workingDir.toAbsolutePath());
        return workingDir;
    }

    /**
     * Validates input documents and options.
     *
     * @param documents the input documents to validate
     * @param options the options to validate
     * @throws DocumentOperationException if validation fails
     */
    private void validateInputs(AnnotatedPluginDocument[] documents, TaxTriageOptions options)
            throws DocumentOperationException {
        if (documents == null || documents.length == 0) {
            throw new DocumentOperationException("No input documents provided");
        }

        if (options == null) {
            throw new DocumentOperationException("No options provided");
        }

        // Validate each document
        for (AnnotatedPluginDocument doc : documents) {
            if (doc == null) {
                throw new DocumentOperationException("Null document in input array");
            }
            if (doc.getDocument() == null) {
                throw new DocumentOperationException("Document " + doc.getName() + " has null content");
            }
        }

        logger.info("Validated " + documents.length + " input documents");
    }

    /**
     * Gets the temporary base directory used for workflow execution.
     *
     * @return the temporary base directory path
     */
    public Path getTempBaseDirectory() {
        return tempBaseDirectory;
    }

    /**
     * Estimates the execution time for a workflow based on input size and options.
     *
     * @param documents the input documents
     * @param options the TaxTriage options
     * @return estimated execution time in minutes
     */
    public int estimateExecutionTime(AnnotatedPluginDocument[] documents, TaxTriageOptions options) {
        if (documents == null || documents.length == 0) {
            return 5; // Minimum time
        }

        // Base time for workflow overhead
        int baseTimeMinutes = 10;

        // Additional time per document (rough estimate)
        int timePerDocument = 2;

        // Additional time based on preset complexity
        int presetMultiplier = 1;
        if (options != null) {
            switch (options.getSequencingPreset()) {
                case ONT:
                    presetMultiplier = 2; // Long reads take longer
                    break;
                case ILLUMINA_PE:
                    presetMultiplier = 2; // Paired-end processing
                    break;
                case ILLUMINA_SE:
                    presetMultiplier = 1; // Single-end is fastest
                    break;
            }
        }

        int estimatedTime = baseTimeMinutes + (documents.length * timePerDocument * presetMultiplier);

        // Cap at reasonable maximum
        return Math.min(estimatedTime, 240); // 4 hours max
    }

    /**
     * Checks if the system requirements are met for running workflows.
     *
     * @return system status information
     */
    public SystemStatus checkSystemRequirements() {
        SystemStatus status = new SystemStatus();

        try {
            // Check Docker availability
            WorkflowExecutor tempExecutor = new WorkflowExecutor();
            status.dockerAvailable = true;
            status.dockerMessage = "Docker is available and running";
        } catch (DockerException e) {
            status.dockerAvailable = false;
            status.dockerMessage = "Docker not available: " + e.getMessage();
        }

        // Check disk space
        try {
            long freeSpace = tempBaseDirectory.toFile().getFreeSpace();
            status.availableDiskSpaceGB = freeSpace / (1024L * 1024L * 1024L);
            status.sufficientDiskSpace = status.availableDiskSpaceGB >= 5; // Require at least 5GB
        } catch (Exception e) {
            status.sufficientDiskSpace = false;
            status.diskSpaceMessage = "Unable to check disk space: " + e.getMessage();
        }

        // Check memory
        Runtime runtime = Runtime.getRuntime();
        long maxMemoryMB = runtime.maxMemory() / (1024L * 1024L);
        status.availableMemoryMB = maxMemoryMB;
        status.sufficientMemory = maxMemoryMB >= 2048; // Require at least 2GB

        status.overall = status.dockerAvailable && status.sufficientDiskSpace && status.sufficientMemory;

        return status;
    }

    /**
     * Creates the default temporary directory for workflows.
     *
     * @return the default temporary directory path
     * @throws IOException if directory creation fails
     */
    private static Path createDefaultTempDirectory() throws IOException {
        String tempDirProperty = System.getProperty("java.io.tmpdir");
        Path systemTempDir = Paths.get(tempDirProperty);
        Path taxTriageTempDir = systemTempDir.resolve("taxtriage");

        Files.createDirectories(taxTriageTempDir);
        return taxTriageTempDir;
    }

    /**
     * Represents system status information for workflow execution.
     */
    public static class SystemStatus {
        public boolean overall;
        public boolean dockerAvailable;
        public String dockerMessage;
        public boolean sufficientDiskSpace;
        public long availableDiskSpaceGB;
        public String diskSpaceMessage;
        public boolean sufficientMemory;
        public long availableMemoryMB;

        /**
         * Gets a human-readable summary of the system status.
         *
         * @return status summary
         */
        public String getSummary() {
            StringBuilder summary = new StringBuilder();
            summary.append("System Status: ").append(overall ? "Ready" : "Not Ready").append("\n");
            summary.append("Docker: ").append(dockerAvailable ? "Available" : "Not Available");
            if (dockerMessage != null) {
                summary.append(" - ").append(dockerMessage);
            }
            summary.append("\n");
            summary.append("Disk Space: ").append(availableDiskSpaceGB).append(" GB available");
            if (!sufficientDiskSpace) {
                summary.append(" (Insufficient - need at least 5 GB)");
            }
            summary.append("\n");
            summary.append("Memory: ").append(availableMemoryMB).append(" MB available");
            if (!sufficientMemory) {
                summary.append(" (Insufficient - need at least 2 GB)");
            }
            return summary.toString();
        }
    }
}