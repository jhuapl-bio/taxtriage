package com.jhuapl.taxtriage.geneious.docker;

import jebl.util.ProgressListener;

import java.io.*;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Manages Docker operations for TaxTriage workflow execution.
 *
 * This class provides a high-level interface for executing Docker containers,
 * specifically for running TaxTriage Nextflow workflows. It handles container
 * lifecycle, volume mounting, and output capturing.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class DockerManager {

    private static final Logger logger = Logger.getLogger(DockerManager.class.getName());

    /** Default TaxTriage Docker image name */
    private static final String DEFAULT_IMAGE = "taxtriage:latest";

    /** Default timeout for Docker operations in minutes */
    private static final int DEFAULT_TIMEOUT_MINUTES = 120;

    private final String dockerImage;

    /**
     * Creates a new DockerManager with the default TaxTriage image.
     *
     * @throws DockerException if Docker is not available
     */
    public DockerManager() throws DockerException {
        this(DEFAULT_IMAGE);
    }

    /**
     * Creates a new DockerManager with a custom Docker image.
     *
     * @param dockerImage the Docker image to use
     * @throws DockerException if Docker is not available
     */
    public DockerManager(String dockerImage) throws DockerException {
        this.dockerImage = dockerImage;
        validateDockerAvailability();
    }

    /**
     * Executes a Nextflow command in a Docker container.
     *
     * @param nextflowCommand the Nextflow command to execute
     * @param inputDirectory the input directory to mount
     * @param outputDirectory the output directory to mount
     * @param workingDirectory the working directory to mount
     * @param progressListener optional progress listener for monitoring
     * @return execution result
     * @throws DockerException if execution fails
     */
    public ExecutionResult executeNextflowCommand(String nextflowCommand,
                                                 Path inputDirectory,
                                                 Path outputDirectory,
                                                 Path workingDirectory,
                                                 ProgressListener progressListener) throws DockerException {
        return executeNextflowCommand(nextflowCommand, inputDirectory, outputDirectory, workingDirectory,
                                    progressListener, DEFAULT_TIMEOUT_MINUTES);
    }

    /**
     * Executes a Nextflow command in a Docker container with custom timeout.
     *
     * @param nextflowCommand the Nextflow command to execute
     * @param inputDirectory the input directory to mount
     * @param outputDirectory the output directory to mount
     * @param workingDirectory the working directory to mount
     * @param progressListener optional progress listener for monitoring
     * @param timeoutMinutes the timeout in minutes
     * @return execution result
     * @throws DockerException if execution fails
     */
    public ExecutionResult executeNextflowCommand(String nextflowCommand,
                                                 Path inputDirectory,
                                                 Path outputDirectory,
                                                 Path workingDirectory,
                                                 ProgressListener progressListener,
                                                 int timeoutMinutes) throws DockerException {
        try {
            // Build the Docker command
            List<String> dockerCommand = buildDockerCommand(nextflowCommand, inputDirectory,
                                                           outputDirectory, workingDirectory);

            if (progressListener != null) {
                progressListener.setMessage("Starting Docker container...");
                progressListener.setProgress(0.1);
            }

            logger.info("Executing Docker command: " + String.join(" ", dockerCommand));

            // Execute the command
            LocalDateTime startTime = LocalDateTime.now();
            ProcessBuilder processBuilder = new ProcessBuilder(dockerCommand);
            processBuilder.redirectErrorStream(false);

            Process process = processBuilder.start();

            // Monitor the process with progress updates
            ExecutionResult result = monitorProcess(process, dockerCommand, progressListener, timeoutMinutes);

            LocalDateTime endTime = LocalDateTime.now();

            logger.info("Docker execution completed. Exit code: " + result.getExitCode() +
                       ", Duration: " + java.time.Duration.between(startTime, endTime).getSeconds() + "s");

            return result;

        } catch (Exception e) {
            logger.log(Level.SEVERE, "Failed to execute Docker command", e);
            throw new DockerException("Failed to execute Nextflow command in Docker", e);
        }
    }

    /**
     * Builds the Docker command with proper volume mounts and options.
     *
     * @param nextflowCommand the Nextflow command to execute
     * @param inputDirectory the input directory
     * @param outputDirectory the output directory
     * @param workingDirectory the working directory
     * @return list of command components
     */
    private List<String> buildDockerCommand(String nextflowCommand,
                                          Path inputDirectory,
                                          Path outputDirectory,
                                          Path workingDirectory) {
        List<String> command = new ArrayList<>();

        // Base Docker command
        command.add("docker");
        command.add("run");
        command.add("--rm"); // Remove container after execution
        command.add("-i"); // Interactive mode

        // Volume mounts
        command.add("-v");
        command.add(inputDirectory.toAbsolutePath() + ":/data/input:ro"); // Read-only input

        command.add("-v");
        command.add(outputDirectory.toAbsolutePath() + ":/data/output:rw"); // Read-write output

        command.add("-v");
        command.add(workingDirectory.toAbsolutePath() + ":/data/work:rw"); // Read-write work

        // Set working directory
        command.add("-w");
        command.add("/data/work");

        // Don't add user mapping - it causes issues when username doesn't exist in container
        // Files will be created as root but this is acceptable for temporary workflow files

        // Environment variables
        command.add("-e");
        command.add("NXF_HOME=/tmp/.nextflow");

        // Docker image
        command.add(dockerImage);

        // The command to execute (using shell to handle complex commands)
        command.add("sh");
        command.add("-c");
        command.add(nextflowCommand);

        return command;
    }

    /**
     * Monitors a running process and captures output.
     *
     * @param process the process to monitor
     * @param command the command being executed
     * @param progressListener optional progress listener
     * @param timeoutMinutes timeout in minutes
     * @return execution result
     * @throws DockerException if monitoring fails
     */
    private ExecutionResult monitorProcess(Process process,
                                         List<String> command,
                                         ProgressListener progressListener,
                                         int timeoutMinutes) throws DockerException {
        LocalDateTime startTime = LocalDateTime.now();

        try {
            // Start output readers
            StringBuilder standardOutput = new StringBuilder();
            StringBuilder errorOutput = new StringBuilder();

            Thread outputReader = new Thread(() -> readStream(process.getInputStream(), standardOutput, progressListener));
            Thread errorReader = new Thread(() -> readStream(process.getErrorStream(), errorOutput, null));

            outputReader.start();
            errorReader.start();

            // Wait for process completion with timeout
            boolean finished = process.waitFor(timeoutMinutes, TimeUnit.MINUTES);

            if (!finished) {
                process.destroyForcibly();
                throw new DockerException("Docker process timed out after " + timeoutMinutes + " minutes");
            }

            // Wait for output readers to finish
            outputReader.join(5000); // 5 second timeout for output reading
            errorReader.join(5000);

            LocalDateTime endTime = LocalDateTime.now();
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
            throw new DockerException("Docker process monitoring was interrupted", e);
        } catch (Exception e) {
            process.destroyForcibly();
            throw new DockerException("Failed to monitor Docker process", e);
        }
    }

    /**
     * Reads from an input stream and appends to a StringBuilder.
     *
     * @param inputStream the stream to read from
     * @param output the StringBuilder to append to
     * @param progressListener optional progress listener for updates
     */
    private void readStream(InputStream inputStream, StringBuilder output, ProgressListener progressListener) {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream))) {
            String line;
            while ((line = reader.readLine()) != null) {
                output.append(line).append("\n");

                // Update progress listener if available
                if (progressListener != null) {
                    updateProgressFromOutput(line, progressListener);
                }

                // Log important lines
                if (line.contains("ERROR") || line.contains("WARN")) {
                    logger.warning("Docker output: " + line);
                } else if (line.contains("Completed") || line.contains("SUCCESS")) {
                    logger.info("Docker output: " + line);
                }
            }
        } catch (IOException e) {
            logger.log(Level.WARNING, "Error reading from Docker process stream", e);
        }
    }

    /**
     * Updates progress based on output from the Docker process.
     *
     * @param outputLine the output line from Docker
     * @param progressListener the progress listener to update
     */
    private void updateProgressFromOutput(String outputLine, ProgressListener progressListener) {
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
     * Validates that Docker is available and accessible.
     *
     * @throws DockerException if Docker is not available
     */
    private void validateDockerAvailability() throws DockerException {
        try {
            // Test Docker availability with simple version command
            ProcessBuilder pb = new ProcessBuilder("docker", "--version");
            Process process = pb.start();

            boolean finished = process.waitFor(10, TimeUnit.SECONDS);
            if (!finished) {
                process.destroyForcibly();
                throw DockerException.dockerNotAvailable();
            }

            int exitCode = process.exitValue();
            if (exitCode != 0) {
                throw DockerException.dockerNotAvailable();
            }

            logger.info("Docker availability confirmed");

        } catch (IOException e) {
            throw new DockerException("Docker command not found", e);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new DockerException("Docker availability check was interrupted", e);
        }
    }

    /**
     * Checks if the TaxTriage Docker image is available.
     *
     * @return true if the image is available
     */
    public boolean isImageAvailable() {
        try {
            ProcessBuilder pb = new ProcessBuilder("docker", "images", "-q", dockerImage);
            Process process = pb.start();

            boolean finished = process.waitFor(30, TimeUnit.SECONDS);
            if (!finished) {
                process.destroyForcibly();
                return false;
            }

            // If exit code is 0 and there's output, image exists
            if (process.exitValue() == 0) {
                try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                    String line = reader.readLine();
                    return line != null && !line.trim().isEmpty();
                }
            }

            return false;

        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to check Docker image availability", e);
            return false;
        }
    }

    /**
     * Checks if Docker is available on the system.
     *
     * @return true if Docker is available and accessible
     */
    public boolean isDockerAvailable() {
        try {
            ProcessBuilder pb = new ProcessBuilder("docker", "--version");
            Process process = pb.start();

            boolean finished = process.waitFor(10, TimeUnit.SECONDS);
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
     * Gets the Docker image name being used.
     *
     * @return the Docker image name
     */
    public String getDockerImage() {
        return dockerImage;
    }

    /**
     * Pulls the TaxTriage Docker image if not available.
     *
     * @throws DockerException if pull fails
     */
    public void pullImageIfNeeded() throws DockerException {
        if (!isImageAvailable()) {
            logger.info("TaxTriage Docker image not found, attempting to pull...");
            pullImage();
        } else {
            logger.info("TaxTriage Docker image is available");
        }
    }

    /**
     * Pulls the TaxTriage Docker image.
     *
     * @throws DockerException if pull fails
     */
    private void pullImage() throws DockerException {
        try {
            ProcessBuilder pb = new ProcessBuilder("docker", "pull", dockerImage);
            Process process = pb.start();

            boolean finished = process.waitFor(10, TimeUnit.MINUTES); // 10 minute timeout for pull
            if (!finished) {
                process.destroyForcibly();
                throw new DockerException("Docker image pull timed out");
            }

            int exitCode = process.exitValue();
            if (exitCode != 0) {
                try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getErrorStream()))) {
                    String errorOutput = reader.lines().collect(
                        java.util.stream.Collectors.joining("\n")
                    );
                    throw DockerException.imageNotFound(dockerImage + ". Error: " + errorOutput);
                }
            }

            logger.info("Successfully pulled Docker image: " + dockerImage);

        } catch (IOException e) {
            throw new DockerException("Failed to execute docker pull command", e);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new DockerException("Docker image pull was interrupted", e);
        }
    }
}