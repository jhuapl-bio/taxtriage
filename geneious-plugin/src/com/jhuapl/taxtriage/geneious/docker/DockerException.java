package com.jhuapl.taxtriage.geneious.docker;

/**
 * Exception thrown when Docker operations fail.
 *
 * This exception is used to wrap Docker-related errors such as:
 * - Docker not being installed or available
 * - Container execution failures
 * - Image pulling problems
 * - Volume mounting issues
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class DockerException extends Exception {

    private final String dockerCommand;
    private final int exitCode;

    /**
     * Creates a new DockerException with a message.
     *
     * @param message the error message
     */
    public DockerException(String message) {
        super(message);
        this.dockerCommand = null;
        this.exitCode = -1;
    }

    /**
     * Creates a new DockerException with a message and cause.
     *
     * @param message the error message
     * @param cause the underlying cause
     */
    public DockerException(String message, Throwable cause) {
        super(message, cause);
        this.dockerCommand = null;
        this.exitCode = -1;
    }

    /**
     * Creates a new DockerException with a cause.
     *
     * @param cause the underlying cause
     */
    public DockerException(Throwable cause) {
        super(cause);
        this.dockerCommand = null;
        this.exitCode = -1;
    }

    /**
     * Creates a new DockerException for a failed Docker command.
     *
     * @param message the error message
     * @param dockerCommand the Docker command that failed
     * @param exitCode the exit code from the command
     */
    public DockerException(String message, String dockerCommand, int exitCode) {
        super(message);
        this.dockerCommand = dockerCommand;
        this.exitCode = exitCode;
    }

    /**
     * Creates a new DockerException for a failed Docker command with cause.
     *
     * @param message the error message
     * @param dockerCommand the Docker command that failed
     * @param exitCode the exit code from the command
     * @param cause the underlying cause
     */
    public DockerException(String message, String dockerCommand, int exitCode, Throwable cause) {
        super(message, cause);
        this.dockerCommand = dockerCommand;
        this.exitCode = exitCode;
    }

    /**
     * Gets the Docker command that failed.
     *
     * @return the Docker command, or null if not applicable
     */
    public String getDockerCommand() {
        return dockerCommand;
    }

    /**
     * Gets the exit code from the failed command.
     *
     * @return the exit code, or -1 if not applicable
     */
    public int getExitCode() {
        return exitCode;
    }

    /**
     * Returns a detailed error message including command and exit code if available.
     *
     * @return detailed error message
     */
    @Override
    public String getMessage() {
        String baseMessage = super.getMessage();
        if (baseMessage == null) {
            baseMessage = "";
        }
        StringBuilder msg = new StringBuilder(baseMessage);

        if (dockerCommand != null) {
            if (msg.length() > 0) {
                msg.append(" ");
            }
            msg.append("(Command: ").append(dockerCommand).append(")");
        }

        if (exitCode != -1) {
            if (msg.length() > 0) {
                msg.append(" ");
            }
            msg.append("(Exit code: ").append(exitCode).append(")");
        }

        return msg.toString();
    }

    /**
     * Creates a DockerException for when Docker is not available.
     *
     * @return configured DockerException
     */
    public static DockerException dockerNotAvailable() {
        return new DockerException("Docker is not available. Please ensure Docker is installed and running.");
    }

    /**
     * Creates a DockerException for when a Docker image is not found.
     *
     * @param imageName the name of the missing image
     * @return configured DockerException
     */
    public static DockerException imageNotFound(String imageName) {
        return new DockerException("Docker image not found: " + imageName +
                                  ". Please ensure the TaxTriage Docker image is available.");
    }

    /**
     * Creates a DockerException for container execution failures.
     *
     * @param containerName the name of the container
     * @param exitCode the exit code
     * @param errorOutput the error output from the container
     * @return configured DockerException
     */
    public static DockerException executionFailed(String containerName, int exitCode, String errorOutput) {
        String message = "Container execution failed: " + containerName;
        if (errorOutput != null && !errorOutput.trim().isEmpty()) {
            message += "\nError output: " + errorOutput;
        }
        return new DockerException(message, "docker run " + containerName, exitCode);
    }

    /**
     * Creates a DockerException for volume mounting failures.
     *
     * @param volumePath the path that failed to mount
     * @return configured DockerException
     */
    public static DockerException volumeMountFailed(String volumePath) {
        return new DockerException("Failed to mount volume: " + volumePath +
                                  ". Please check that the path exists and is accessible.");
    }
}