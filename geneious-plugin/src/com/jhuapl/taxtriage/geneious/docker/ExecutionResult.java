package com.jhuapl.taxtriage.geneious.docker;

import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Represents the result of a Docker container execution.
 *
 * This class encapsulates all information about a completed Docker execution,
 * including exit code, output, timing information, and any errors that occurred.
 * It provides methods to analyze the execution results and determine success/failure.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class ExecutionResult {

    private final String command;
    private final int exitCode;
    private final String standardOutput;
    private final String errorOutput;
    private final LocalDateTime startTime;
    private final LocalDateTime endTime;
    private final Duration executionTime;
    private final List<String> outputLines;
    private final List<String> errorLines;

    /**
     * Creates a new ExecutionResult.
     *
     * @param command the Docker command that was executed
     * @param exitCode the exit code from the execution
     * @param standardOutput the standard output
     * @param errorOutput the error output
     * @param startTime when execution started
     * @param endTime when execution completed
     */
    public ExecutionResult(String command, int exitCode, String standardOutput, String errorOutput,
                          LocalDateTime startTime, LocalDateTime endTime) {
        this.command = command;
        this.exitCode = exitCode;
        this.standardOutput = standardOutput != null ? standardOutput : "";
        this.errorOutput = errorOutput != null ? errorOutput : "";
        this.startTime = startTime;
        this.endTime = endTime;
        this.executionTime = Duration.between(startTime, endTime);

        // Parse output into lines for easier processing
        this.outputLines = parseLines(this.standardOutput);
        this.errorLines = parseLines(this.errorOutput);
    }

    /**
     * Creates a successful ExecutionResult with minimal information.
     *
     * @param command the command that was executed
     * @param output the output from the command
     * @return successful ExecutionResult
     */
    public static ExecutionResult success(String command, String output) {
        LocalDateTime now = LocalDateTime.now();
        return new ExecutionResult(command, 0, output, "", now, now);
    }

    /**
     * Creates a failed ExecutionResult.
     *
     * @param command the command that was executed
     * @param exitCode the non-zero exit code
     * @param errorOutput the error output
     * @return failed ExecutionResult
     */
    public static ExecutionResult failure(String command, int exitCode, String errorOutput) {
        LocalDateTime now = LocalDateTime.now();
        return new ExecutionResult(command, exitCode, "", errorOutput, now, now);
    }

    /**
     * Creates a failed ExecutionResult from an exception.
     *
     * @param command the command that was attempted
     * @param exception the exception that occurred
     * @return failed ExecutionResult
     */
    public static ExecutionResult fromException(String command, Exception exception) {
        LocalDateTime now = LocalDateTime.now();
        String errorMessage = exception.getMessage();
        if (errorMessage == null) {
            errorMessage = exception.getClass().getSimpleName();
        }
        return new ExecutionResult(command, -1, "", errorMessage, now, now);
    }

    /**
     * Parses text into individual lines.
     *
     * @param text the text to parse
     * @return list of lines
     */
    private List<String> parseLines(String text) {
        if (text == null || text.isEmpty()) {
            return Collections.emptyList();
        }

        List<String> lines = new ArrayList<>();
        String[] splitLines = text.split("\r?\n");
        for (String line : splitLines) {
            lines.add(line);
        }
        return Collections.unmodifiableList(lines);
    }

    /**
     * Gets the Docker command that was executed.
     *
     * @return the command
     */
    public String getCommand() {
        return command;
    }

    /**
     * Gets the exit code from the execution.
     *
     * @return the exit code (0 indicates success)
     */
    public int getExitCode() {
        return exitCode;
    }

    /**
     * Gets the standard output from the execution.
     *
     * @return the standard output
     */
    public String getStandardOutput() {
        return standardOutput;
    }

    /**
     * Gets the error output from the execution.
     *
     * @return the error output
     */
    public String getErrorOutput() {
        return errorOutput;
    }

    /**
     * Gets the execution start time.
     *
     * @return the start time
     */
    public LocalDateTime getStartTime() {
        return startTime;
    }

    /**
     * Gets the execution end time.
     *
     * @return the end time
     */
    public LocalDateTime getEndTime() {
        return endTime;
    }

    /**
     * Gets the total execution time.
     *
     * @return the execution duration
     */
    public Duration getExecutionTime() {
        return executionTime;
    }

    /**
     * Gets the standard output as individual lines.
     *
     * @return list of output lines
     */
    public List<String> getOutputLines() {
        return outputLines;
    }

    /**
     * Gets the error output as individual lines.
     *
     * @return list of error lines
     */
    public List<String> getErrorLines() {
        return errorLines;
    }

    /**
     * Checks if the execution was successful (exit code 0).
     *
     * @return true if successful
     */
    public boolean isSuccessful() {
        return exitCode == 0;
    }

    /**
     * Checks if the execution was successful (exit code 0).
     * Alias for isSuccessful() for backward compatibility.
     *
     * @return true if successful
     */
    public boolean isSuccess() {
        return isSuccessful();
    }

    /**
     * Checks if the execution failed (non-zero exit code).
     *
     * @return true if failed
     */
    public boolean isFailed() {
        return exitCode != 0;
    }

    /**
     * Checks if there is any error output.
     *
     * @return true if there is error output
     */
    public boolean hasErrorOutput() {
        return errorOutput != null && !errorOutput.trim().isEmpty();
    }

    /**
     * Gets the execution time in seconds.
     *
     * @return execution time in seconds
     */
    public long getExecutionSeconds() {
        return executionTime.getSeconds();
    }

    /**
     * Gets the execution time in milliseconds.
     *
     * @return execution time in milliseconds
     */
    public long getExecutionMillis() {
        return executionTime.toMillis();
    }

    /**
     * Searches for a specific pattern in the output lines.
     *
     * @param pattern the pattern to search for
     * @return true if the pattern is found in any output line
     */
    public boolean containsInOutput(String pattern) {
        return outputLines.stream().anyMatch(line -> line.contains(pattern));
    }

    /**
     * Searches for a specific pattern in the error lines.
     *
     * @param pattern the pattern to search for
     * @return true if the pattern is found in any error line
     */
    public boolean containsInError(String pattern) {
        return errorLines.stream().anyMatch(line -> line.contains(pattern));
    }

    /**
     * Gets the first line of output that matches a pattern.
     *
     * @param pattern the pattern to search for
     * @return the first matching line, or null if not found
     */
    public String findInOutput(String pattern) {
        return outputLines.stream()
                .filter(line -> line.contains(pattern))
                .findFirst()
                .orElse(null);
    }

    /**
     * Gets the first line of error output that matches a pattern.
     *
     * @param pattern the pattern to search for
     * @return the first matching error line, or null if not found
     */
    public String findInError(String pattern) {
        return errorLines.stream()
                .filter(line -> line.contains(pattern))
                .findFirst()
                .orElse(null);
    }

    /**
     * Creates a summary string of the execution result.
     *
     * @return execution summary
     */
    public String getSummary() {
        StringBuilder summary = new StringBuilder();
        summary.append("Command: ").append(command).append("\n");
        summary.append("Exit Code: ").append(exitCode).append(" (").append(isSuccessful() ? "SUCCESS" : "FAILED").append(")\n");
        summary.append("Execution Time: ").append(getExecutionSeconds()).append(" seconds\n");
        summary.append("Output Lines: ").append(outputLines.size()).append("\n");
        summary.append("Error Lines: ").append(errorLines.size()).append("\n");

        if (hasErrorOutput()) {
            summary.append("Error Output: ").append(errorOutput.substring(0, Math.min(200, errorOutput.length())));
            if (errorOutput.length() > 200) {
                summary.append("...");
            }
        }

        return summary.toString();
    }

    @Override
    public String toString() {
        return "ExecutionResult{" +
               "command='" + command + '\'' +
               ", exitCode=" + exitCode +
               ", successful=" + isSuccessful() +
               ", executionTime=" + getExecutionSeconds() + "s" +
               ", outputLines=" + outputLines.size() +
               ", errorLines=" + errorLines.size() +
               '}';
    }
}