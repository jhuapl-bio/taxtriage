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
    private final String dockerCommand;

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
        this.dockerCommand = findDockerCommand();
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

            String fullCommand = String.join(" ", dockerCommand);
            logger.info("Executing Docker command: " + fullCommand);
            System.out.println("========================================");
            System.out.println("[TaxTriage] DOCKER COMMAND:");
            System.out.println(fullCommand);
            System.out.println("========================================");
            System.out.flush();

            // Execute the command
            LocalDateTime startTime = LocalDateTime.now();
            ProcessBuilder processBuilder = new ProcessBuilder(dockerCommand);
            processBuilder.redirectErrorStream(false);

            Process process = processBuilder.start();

            System.out.println("[TaxTriage] Docker process started, PID (if available): " + process.pid());
            System.out.flush();

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
        command.add(dockerCommand);
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

        // Environment variables for Nextflow
        command.add("-e");
        command.add("NXF_HOME=/tmp/.nextflow");

        // Force plain text output (disable ANSI colors and terminal features)
        command.add("-e");
        command.add("NXF_ANSI_LOG=false");

        command.add("-e");
        command.add("NXF_ANSI_SUMMARY=false");

        // Set log level to info for more output
        command.add("-e");
        command.add("NXF_LOG_LEVEL=info");

        // Ensure output is unbuffered
        command.add("-e");
        command.add("PYTHONUNBUFFERED=1");

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
     * Uses ExecutionMonitor for detailed progress tracking.
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
            // Use ExecutionMonitor for enhanced progress tracking
            ExecutionMonitor monitor = new ExecutionMonitor();

            if (progressListener != null) {
                progressListener.setMessage("TaxTriage: Initializing workflow...");
            }

            // Monitor the execution with the ExecutionMonitor
            java.util.concurrent.CompletableFuture<ExecutionResult> monitorFuture =
                monitor.monitorExecution(process, progressListener);

            // Wait for completion with timeout
            try {
                ExecutionResult result = monitorFuture.get(timeoutMinutes, TimeUnit.MINUTES);
                return result;
            } catch (java.util.concurrent.TimeoutException e) {
                monitor.requestCancellation();
                process.destroyForcibly();
                throw new DockerException("Docker process timed out after " + timeoutMinutes + " minutes");
            }

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
     * Validates that Docker is available and accessible.
     *
     * @throws DockerException if Docker is not available
     */
    private void validateDockerAvailability() throws DockerException {
        logger.info("=== Docker Validation Starting ===");
        logger.info("Docker command path: " + dockerCommand);

        try {
            // Test Docker daemon availability with 'docker ps' command
            // This verifies both that docker exists AND that the daemon is running
            logger.info("Executing: " + dockerCommand + " ps");
            ProcessBuilder pb = new ProcessBuilder(dockerCommand, "ps");
            Process process = pb.start();

            boolean finished = process.waitFor(10, TimeUnit.SECONDS);
            if (!finished) {
                process.destroyForcibly();
                String error = "Docker daemon is not responding. Please ensure Docker Desktop is running.";
                logger.severe("Docker validation failed: " + error);
                throw new DockerException(error);
            }

            int exitCode = process.exitValue();
            logger.info("Docker ps exit code: " + exitCode);

            if (exitCode != 0) {
                // Capture both stdout and stderr for better diagnostics
                StringBuilder errorOutput = new StringBuilder();
                StringBuilder stdOutput = new StringBuilder();

                try (BufferedReader errReader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
                     BufferedReader outReader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                    String line;
                    while ((line = errReader.readLine()) != null) {
                        errorOutput.append(line).append("\n");
                        logger.warning("Docker stderr: " + line);
                    }
                    while ((line = outReader.readLine()) != null) {
                        stdOutput.append(line).append("\n");
                        logger.info("Docker stdout: " + line);
                    }
                }

                String errorMsg = "Docker daemon is not running. Please start Docker Desktop.\n" +
                                 "Command: " + dockerCommand + " ps\n" +
                                 "Exit code: " + exitCode + "\n" +
                                 "Error output: " + errorOutput.toString() +
                                 "Standard output: " + stdOutput.toString();
                logger.severe("Docker validation failed: " + errorMsg);
                throw new DockerException(errorMsg);
            }

            logger.info("Docker daemon availability confirmed at: " + dockerCommand);
            logger.info("=== Docker Validation Complete - SUCCESS ===");

        } catch (IOException e) {
            String errorMsg = "Docker command not found at: " + dockerCommand +
                             ". Please ensure Docker is installed and in your PATH.";
            logger.log(Level.SEVERE, "Docker validation failed with IOException: " + errorMsg, e);
            throw new DockerException(errorMsg, e);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            String errorMsg = "Docker availability check was interrupted";
            logger.log(Level.SEVERE, "Docker validation interrupted", e);
            throw new DockerException(errorMsg, e);
        }
    }

    /**
     * Finds the Docker command executable, checking common installation locations.
     *
     * @return the path to the docker executable
     */
    private String findDockerCommand() {
        logger.info("=== Searching for Docker executable ===");

        // Log current PATH for debugging
        String path = System.getenv("PATH");
        logger.info("System PATH: " + (path != null ? path : "null"));

        // Try standard 'docker' command first (will use PATH)
        logger.info("Checking if 'docker' command is available in PATH...");
        if (isCommandAvailable("docker")) {
            logger.info("Found 'docker' in PATH");
            return "docker";
        }
        logger.info("'docker' not found in PATH, checking common installation locations...");

        // Try common installation paths on macOS and Linux
        String[] commonPaths = {
            "/usr/local/bin/docker",           // macOS Homebrew, Linux standard
            "/usr/bin/docker",                 // Linux standard
            "/opt/homebrew/bin/docker",        // macOS Apple Silicon Homebrew
            "/Applications/Docker.app/Contents/Resources/bin/docker"  // macOS Docker Desktop
        };

        for (String checkPath : commonPaths) {
            logger.info("Checking: " + checkPath);
            File dockerFile = new File(checkPath);
            if (dockerFile.exists()) {
                logger.info("  File exists: " + dockerFile.exists());
                logger.info("  Is executable: " + dockerFile.canExecute());
                if (dockerFile.canExecute()) {
                    logger.info("SUCCESS: Found Docker at: " + checkPath);
                    return checkPath;
                }
            } else {
                logger.info("  File does not exist");
            }
        }

        // Default to 'docker' and let it fail with a clear error
        logger.warning("Docker not found in any common locations. Will try 'docker' command and may fail.");
        logger.warning("If Docker is installed in a non-standard location, consider adding it to PATH.");
        return "docker";
    }

    /**
     * Checks if a command is available in the system PATH.
     *
     * @param command the command to check
     * @return true if the command is available
     */
    private boolean isCommandAvailable(String command) {
        try {
            logger.fine("Testing command availability: " + command + " --version");
            ProcessBuilder pb = new ProcessBuilder(command, "--version");
            Process process = pb.start();
            boolean finished = process.waitFor(5, TimeUnit.SECONDS);
            if (!finished) {
                logger.fine("Command timed out: " + command);
                process.destroyForcibly();
                return false;
            }
            int exitCode = process.exitValue();
            logger.fine("Command " + command + " exit code: " + exitCode);
            return exitCode == 0;
        } catch (Exception e) {
            logger.fine("Command not available: " + command + " - " + e.getMessage());
            return false;
        }
    }

    /**
     * Checks if the TaxTriage Docker image is available.
     *
     * @return true if the image is available
     */
    public boolean isImageAvailable() {
        try {
            ProcessBuilder pb = new ProcessBuilder(dockerCommand, "images", "-q", dockerImage);
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
            ProcessBuilder pb = new ProcessBuilder(dockerCommand, "ps");
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
     * Gets the Docker command path being used.
     *
     * @return the Docker command path
     */
    public String getDockerCommand() {
        return dockerCommand;
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
            ProcessBuilder pb = new ProcessBuilder(dockerCommand, "pull", dockerImage);
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