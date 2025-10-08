package com.jhuapl.taxtriage.geneious.docker;

import jebl.util.ProgressListener;

import java.nio.file.Path;
import java.time.LocalDateTime;

/**
 * Mock implementation of Docker operations for testing without requiring actual Docker.
 *
 * This class simulates Docker behavior for unit tests, allowing complete testing
 * of the TaxTriage plugin functionality in environments where Docker is not available.
 *
 * Note: This is a standalone mock that does not extend DockerManager to avoid
 * Docker validation requirements during test initialization.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class MockDockerManager {

    private boolean dockerAvailable = true;
    private boolean shouldFailExecution = false;
    private String mockOutput = "Mock Nextflow execution completed successfully";
    private String mockError = "";
    private int mockExitCode = 0;
    private long executionDelayMs = 100;
    private final String dockerImage;

    /**
     * Creates a new MockDockerManager with default successful behavior.
     */
    public MockDockerManager() {
        this.dockerImage = "mock-docker";
    }

    /**
     * Creates a MockDockerManager with custom docker image name.
     *
     * @param dockerImage the mock docker image name
     */
    public MockDockerManager(String dockerImage) {
        this.dockerImage = dockerImage;
    }

    public boolean isDockerAvailable() {
        return dockerAvailable;
    }

    public void pullImageIfNeeded() throws DockerException {
        if (!dockerAvailable) {
            throw new DockerException("Mock Docker is not available");
        }
        // Mock image pull - do nothing
    }

    public ExecutionResult executeNextflowCommand(
            String nextflowCommand,
            Path inputDirectory,
            Path outputDirectory,
            Path workingDirectory,
            ProgressListener progressListener) throws DockerException {

        if (!dockerAvailable) {
            throw new DockerException("Mock Docker is not available");
        }

        if (shouldFailExecution) {
            throw new DockerException("Mock execution failure: " + mockError);
        }

        if (progressListener != null) {
            progressListener.setMessage("Mock: Starting Nextflow execution");
            progressListener.setProgress(0.0);

            // Simulate progress updates
            simulateProgress(progressListener);

            progressListener.setMessage("Mock: Execution completed");
            progressListener.setProgress(1.0);
        }

        LocalDateTime now = LocalDateTime.now();
        return new ExecutionResult(nextflowCommand, mockExitCode, mockOutput, mockError, now, now);
    }

    public ExecutionResult executeNextflowCommand(
            String nextflowCommand,
            Path inputDirectory,
            Path outputDirectory,
            Path workingDirectory,
            ProgressListener progressListener,
            int timeoutMinutes) throws DockerException {
        // Timeout is ignored in mock, just delegate to main method
        return executeNextflowCommand(nextflowCommand, inputDirectory, outputDirectory, workingDirectory, progressListener);
    }

    public String getDockerImage() {
        return dockerImage;
    }

    public boolean isImageAvailable() {
        return dockerAvailable;
    }

    /**
     * Sets whether Docker should be reported as available.
     *
     * @param available true to simulate Docker availability, false otherwise
     */
    public void setDockerAvailable(boolean available) {
        this.dockerAvailable = available;
    }

    /**
     * Sets whether execution should fail.
     *
     * @param shouldFail true to simulate execution failure, false for success
     */
    public void setShouldFailExecution(boolean shouldFail) {
        this.shouldFailExecution = shouldFail;
    }

    /**
     * Sets the mock output for successful executions.
     *
     * @param output the mock output string
     */
    public void setMockOutput(String output) {
        this.mockOutput = output;
    }

    /**
     * Sets the mock error output.
     *
     * @param error the mock error string
     */
    public void setMockError(String error) {
        this.mockError = error;
    }

    /**
     * Sets the mock exit code.
     *
     * @param exitCode the mock exit code
     */
    public void setMockExitCode(int exitCode) {
        this.mockExitCode = exitCode;
    }

    /**
     * Sets the execution delay for simulating processing time.
     *
     * @param delayMs delay in milliseconds
     */
    public void setExecutionDelay(long delayMs) {
        this.executionDelayMs = delayMs;
    }

    /**
     * Simulates progress updates during execution.
     *
     * @param progressListener the progress listener to update
     */
    private void simulateProgress(ProgressListener progressListener) {
        try {
            String[] stages = {
                "Mock: Initializing workflow",
                "Mock: Processing input files",
                "Mock: Running taxonomic classification",
                "Mock: Generating reports"
            };

            for (int i = 0; i < stages.length; i++) {
                if (progressListener != null) {
                    progressListener.setMessage(stages[i]);
                    progressListener.setProgress((double) (i + 1) / (stages.length + 1));
                }
                Thread.sleep(executionDelayMs / stages.length);
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
    }
}