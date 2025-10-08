package com.jhuapl.taxtriage.geneious.execution;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.jhuapl.taxtriage.geneious.TaxTriageOptions;
import com.jhuapl.taxtriage.geneious.config.TaxTriageConfig;

import java.io.File;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.util.List;
import java.util.UUID;

/**
 * Stores workflow state and configuration for TaxTriage execution.
 *
 * This class maintains all the context information needed during workflow execution,
 * including input/output directories, configuration, execution status, and results.
 * It provides a centralized way to track workflow state across different execution phases.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class WorkflowContext {

    /** Possible workflow execution states */
    public enum ExecutionState {
        INITIALIZED,
        EXPORTING_SEQUENCES,
        GENERATING_CONFIG,
        CREATING_SAMPLESHEET,
        EXECUTING_WORKFLOW,
        MONITORING_PROGRESS,
        IMPORTING_RESULTS,
        COMPLETED,
        FAILED,
        CANCELLED
    }

    // Workflow identification
    private final String workflowId;
    private final LocalDateTime createdAt;
    private LocalDateTime startedAt;
    private LocalDateTime completedAt;

    // Input data
    private final List<AnnotatedPluginDocument> inputDocuments;
    private final TaxTriageOptions options;
    private final TaxTriageConfig config;

    // Execution directories
    private final Path workingDirectory;
    private final Path inputDirectory;
    private final Path outputDirectory;
    private final Path configDirectory;

    // Generated files
    private File sampleSheetFile;
    private File configFile;
    private File paramsFile;
    private List<File> exportedFiles;

    // Execution state
    private ExecutionState state;
    private String currentMessage;
    private double progress;
    private Exception lastError;

    // Results
    private File workflowReport;
    private File taxonomyResults;
    private File qualityReport;
    private List<File> outputFiles;

    /**
     * Creates a new workflow context.
     *
     * @param inputDocuments the input Geneious documents
     * @param options the TaxTriage options
     * @param config the TaxTriage configuration
     * @param workingDirectory the working directory for temporary files
     */
    public WorkflowContext(List<AnnotatedPluginDocument> inputDocuments,
                          TaxTriageOptions options,
                          TaxTriageConfig config,
                          Path workingDirectory) {
        this.workflowId = UUID.randomUUID().toString();
        this.createdAt = LocalDateTime.now();
        this.inputDocuments = inputDocuments;
        this.options = options;
        this.config = config;
        this.workingDirectory = workingDirectory;
        this.inputDirectory = workingDirectory.resolve("input");
        this.outputDirectory = workingDirectory.resolve("output");
        this.configDirectory = workingDirectory.resolve("config");
        this.state = ExecutionState.INITIALIZED;
        this.progress = 0.0;
        this.currentMessage = "Workflow initialized";
    }

    /**
     * Gets the unique workflow identifier.
     *
     * @return the workflow ID
     */
    public String getWorkflowId() {
        return workflowId;
    }

    /**
     * Gets the workflow creation timestamp.
     *
     * @return creation timestamp
     */
    public LocalDateTime getCreatedAt() {
        return createdAt;
    }

    /**
     * Gets the workflow start timestamp.
     *
     * @return start timestamp, or null if not started
     */
    public LocalDateTime getStartedAt() {
        return startedAt;
    }

    /**
     * Sets the workflow start timestamp.
     *
     * @param startedAt start timestamp
     */
    public void setStartedAt(LocalDateTime startedAt) {
        this.startedAt = startedAt;
    }

    /**
     * Gets the workflow completion timestamp.
     *
     * @return completion timestamp, or null if not completed
     */
    public LocalDateTime getCompletedAt() {
        return completedAt;
    }

    /**
     * Sets the workflow completion timestamp.
     *
     * @param completedAt completion timestamp
     */
    public void setCompletedAt(LocalDateTime completedAt) {
        this.completedAt = completedAt;
    }

    /**
     * Gets the input Geneious documents.
     *
     * @return input documents
     */
    public List<AnnotatedPluginDocument> getInputDocuments() {
        return inputDocuments;
    }

    /**
     * Gets the TaxTriage options.
     *
     * @return TaxTriage options
     */
    public TaxTriageOptions getOptions() {
        return options;
    }

    /**
     * Gets the TaxTriage configuration.
     *
     * @return TaxTriage configuration
     */
    public TaxTriageConfig getConfig() {
        return config;
    }

    /**
     * Gets the working directory.
     *
     * @return working directory path
     */
    public Path getWorkingDirectory() {
        return workingDirectory;
    }

    /**
     * Gets the input directory for exported sequences.
     *
     * @return input directory path
     */
    public Path getInputDirectory() {
        return inputDirectory;
    }

    /**
     * Gets the output directory for workflow results.
     *
     * @return output directory path
     */
    public Path getOutputDirectory() {
        return outputDirectory;
    }

    /**
     * Gets the configuration directory.
     *
     * @return configuration directory path
     */
    public Path getConfigDirectory() {
        return configDirectory;
    }

    /**
     * Gets the samplesheet file.
     *
     * @return samplesheet file, or null if not generated
     */
    public File getSampleSheetFile() {
        return sampleSheetFile;
    }

    /**
     * Sets the samplesheet file.
     *
     * @param sampleSheetFile the samplesheet file
     */
    public void setSampleSheetFile(File sampleSheetFile) {
        this.sampleSheetFile = sampleSheetFile;
    }

    /**
     * Gets the configuration file.
     *
     * @return configuration file, or null if not generated
     */
    public File getConfigFile() {
        return configFile;
    }

    /**
     * Sets the configuration file.
     *
     * @param configFile the configuration file
     */
    public void setConfigFile(File configFile) {
        this.configFile = configFile;
    }

    /**
     * Gets the parameters file.
     *
     * @return parameters file, or null if not generated
     */
    public File getParamsFile() {
        return paramsFile;
    }

    /**
     * Sets the parameters file.
     *
     * @param paramsFile the parameters file
     */
    public void setParamsFile(File paramsFile) {
        this.paramsFile = paramsFile;
    }

    /**
     * Gets the list of exported sequence files.
     *
     * @return exported files, or null if not exported
     */
    public List<File> getExportedFiles() {
        return exportedFiles;
    }

    /**
     * Sets the list of exported sequence files.
     *
     * @param exportedFiles the exported files
     */
    public void setExportedFiles(List<File> exportedFiles) {
        this.exportedFiles = exportedFiles;
    }

    /**
     * Gets the current execution state.
     *
     * @return execution state
     */
    public ExecutionState getState() {
        return state;
    }

    /**
     * Sets the execution state.
     *
     * @param state the new execution state
     */
    public void setState(ExecutionState state) {
        this.state = state;
    }

    /**
     * Gets the current progress message.
     *
     * @return current message
     */
    public String getCurrentMessage() {
        return currentMessage;
    }

    /**
     * Sets the current progress message.
     *
     * @param currentMessage the progress message
     */
    public void setCurrentMessage(String currentMessage) {
        this.currentMessage = currentMessage;
    }

    /**
     * Gets the current progress (0.0 to 1.0).
     *
     * @return progress value
     */
    public double getProgress() {
        return progress;
    }

    /**
     * Sets the current progress.
     *
     * @param progress the progress value (0.0 to 1.0)
     */
    public void setProgress(double progress) {
        this.progress = Math.max(0.0, Math.min(1.0, progress));
    }

    /**
     * Gets the last error that occurred.
     *
     * @return last error, or null if none
     */
    public Exception getLastError() {
        return lastError;
    }

    /**
     * Sets the last error.
     *
     * @param lastError the error that occurred
     */
    public void setLastError(Exception lastError) {
        this.lastError = lastError;
    }

    /**
     * Gets the workflow report file.
     *
     * @return workflow report file, or null if not available
     */
    public File getWorkflowReport() {
        return workflowReport;
    }

    /**
     * Sets the workflow report file.
     *
     * @param workflowReport the workflow report file
     */
    public void setWorkflowReport(File workflowReport) {
        this.workflowReport = workflowReport;
    }

    /**
     * Gets the taxonomy results file.
     *
     * @return taxonomy results file, or null if not available
     */
    public File getTaxonomyResults() {
        return taxonomyResults;
    }

    /**
     * Sets the taxonomy results file.
     *
     * @param taxonomyResults the taxonomy results file
     */
    public void setTaxonomyResults(File taxonomyResults) {
        this.taxonomyResults = taxonomyResults;
    }

    /**
     * Gets the quality report file.
     *
     * @return quality report file, or null if not available
     */
    public File getQualityReport() {
        return qualityReport;
    }

    /**
     * Sets the quality report file.
     *
     * @param qualityReport the quality report file
     */
    public void setQualityReport(File qualityReport) {
        this.qualityReport = qualityReport;
    }

    /**
     * Gets all output files generated by the workflow.
     *
     * @return list of output files, or null if not available
     */
    public List<File> getOutputFiles() {
        return outputFiles;
    }

    /**
     * Sets the output files generated by the workflow.
     *
     * @param outputFiles the output files
     */
    public void setOutputFiles(List<File> outputFiles) {
        this.outputFiles = outputFiles;
    }

    /**
     * Checks if the workflow is in a terminal state.
     *
     * @return true if workflow is completed, failed, or cancelled
     */
    public boolean isTerminal() {
        return state == ExecutionState.COMPLETED ||
               state == ExecutionState.FAILED ||
               state == ExecutionState.CANCELLED;
    }

    /**
     * Checks if the workflow completed successfully.
     *
     * @return true if workflow completed successfully
     */
    public boolean isSuccessful() {
        return state == ExecutionState.COMPLETED && lastError == null;
    }

    /**
     * Updates the workflow state and progress.
     *
     * @param newState the new execution state
     * @param message the progress message
     * @param progress the progress value (0.0 to 1.0)
     */
    public void updateState(ExecutionState newState, String message, double progress) {
        this.state = newState;
        this.currentMessage = message;
        setProgress(progress);
    }

    /**
     * Marks the workflow as failed with an error.
     *
     * @param error the error that caused the failure
     * @param message the failure message
     */
    public void markFailed(Exception error, String message) {
        this.state = ExecutionState.FAILED;
        this.lastError = error;
        this.currentMessage = message != null ? message : "Workflow failed";
        this.completedAt = LocalDateTime.now();
    }

    /**
     * Marks the workflow as completed successfully.
     */
    public void markCompleted() {
        this.state = ExecutionState.COMPLETED;
        this.currentMessage = "Workflow completed successfully";
        this.progress = 1.0;
        this.completedAt = LocalDateTime.now();
    }

    /**
     * Marks the workflow as cancelled.
     */
    public void markCancelled() {
        this.state = ExecutionState.CANCELLED;
        this.currentMessage = "Workflow cancelled";
        this.completedAt = LocalDateTime.now();
    }

    /**
     * Gets a summary string of the workflow state.
     *
     * @return workflow summary
     */
    public String getSummary() {
        StringBuilder summary = new StringBuilder();
        summary.append("Workflow ID: ").append(workflowId).append("\n");
        summary.append("State: ").append(state).append("\n");
        summary.append("Progress: ").append(String.format("%.1f%%", progress * 100)).append("\n");
        summary.append("Message: ").append(currentMessage).append("\n");
        summary.append("Created: ").append(createdAt).append("\n");

        if (startedAt != null) {
            summary.append("Started: ").append(startedAt).append("\n");
        }

        if (completedAt != null) {
            summary.append("Completed: ").append(completedAt).append("\n");
        }

        if (lastError != null) {
            summary.append("Last Error: ").append(lastError.getMessage()).append("\n");
        }

        return summary.toString();
    }
}