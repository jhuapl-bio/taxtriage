package com.jhuapl.taxtriage.geneious.docker;

import jebl.util.ProgressListener;

import java.nio.file.Path;
import java.util.concurrent.CompletableFuture;

/**
 * Mock implementation of Docker operations for testing without requiring actual Docker.
 *
 * This class simulates Docker behavior for unit tests, allowing complete testing
 * of the TaxTriage plugin functionality in environments where Docker is not available.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class MockDockerManager extends DockerManager {

    private boolean dockerAvailable = true;
    private boolean shouldFailExecution = false;
    private String mockOutput = "Mock Nextflow execution completed successfully";
    private String mockError = "";
    private int mockExitCode = 0;
    private long executionDelayMs = 100;

    /**
     * Creates a new MockDockerManager with default successful behavior.
     *
     * @throws DockerException never thrown in mock implementation
     */
    public MockDockerManager() throws DockerException {
        super(new MockVolumeMapper(), new MockExecutionMonitor(), "mock-docker");
    }

    /**
     * Creates a MockDockerManager with custom mock components.
     *
     * @param volumeMapper mock volume mapper
     * @param executionMonitor mock execution monitor
     * @throws DockerException never thrown in mock implementation
     */
    public MockDockerManager(VolumeMapper volumeMapper, ExecutionMonitor executionMonitor) throws DockerException {
        super(volumeMapper, executionMonitor, "mock-docker");
    }

    @Override
    public boolean isDockerAvailable() {
        return dockerAvailable;
    }

    @Override
    public void pullImageIfNeeded(String imageName, ProgressListener progressListener) throws DockerException {
        if (!dockerAvailable) {
            throw new DockerException("Mock Docker is not available");
        }

        if (progressListener != null) {
            progressListener.setMessage("Mock: Checking image " + imageName);
            try {
                Thread.sleep(50);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }
            progressListener.setMessage("Mock: Image " + imageName + " is ready");
        }
    }

    @Override
    public ExecutionResult executeNextflowCommand(
            String imageName,
            String nextflowCommand,
            Path inputDirectory,
            Path outputDirectory,
            Path workDirectory,
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

        return new ExecutionResult(mockExitCode, mockOutput, mockError);
    }

    @Override
    public void stopContainer(String containerId) throws DockerException {
        if (!dockerAvailable) {
            throw new DockerException("Mock Docker is not available");
        }
        // Mock container stop - do nothing
    }

    @Override
    public void removeContainer(String containerId) throws DockerException {
        if (!dockerAvailable) {
            throw new DockerException("Mock Docker is not available");
        }
        // Mock container removal - do nothing
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

    /**
     * Mock VolumeMapper implementation for testing.
     */
    public static class MockVolumeMapper extends VolumeMapper {

        public MockVolumeMapper() {
            super("mock-platform");
        }

        @Override
        public boolean validatePermissions(Path path) {
            return true; // Always return true for mock
        }
    }

    /**
     * Mock ExecutionMonitor implementation for testing.
     */
    public static class MockExecutionMonitor extends ExecutionMonitor {

        @Override
        public CompletableFuture<ExecutionResult> monitorExecution(Process process, ProgressListener progressListener) {
            return CompletableFuture.supplyAsync(() -> {
                if (progressListener != null) {
                    progressListener.setMessage("Mock: Monitoring execution");
                    progressListener.setProgress(1.0);
                }
                return new ExecutionResult(0, "Mock execution output", "");
            });
        }

        @Override
        public double parseProgress(String line) {
            if (line != null && line.contains("Mock")) {
                return 0.5; // Return mock progress
            }
            return super.parseProgress(line);
        }
    }
}