package com.jhuapl.taxtriage.geneious.docker;

import jebl.util.ProgressListener;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Monitors the execution of Nextflow processes and provides progress tracking,
 * error detection, and cancellation support for TaxTriage workflows.
 *
 * This class parses Nextflow output in real-time to extract progress information,
 * detect completion status, and handle error conditions during workflow execution.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class ExecutionMonitor {

    private static final Logger logger = Logger.getLogger(ExecutionMonitor.class.getName());

    /** Pattern to match Nextflow progress lines */
    private static final Pattern PROGRESS_PATTERN = Pattern.compile(
            "\\[([^\\]]+)\\]\\s+process\\s+>\\s+([^\\s]+)\\s+\\[(\\d+)\\s+of\\s+(\\d+)\\]");

    /** Pattern to match Nextflow completion status */
    private static final Pattern COMPLETION_PATTERN = Pattern.compile(
            "Pipeline completed at:\\s+(.+)");

    /** Pattern to match Nextflow error messages */
    private static final Pattern ERROR_PATTERN = Pattern.compile(
            "(ERROR|WARN|Exception|Failed|Error).*", Pattern.CASE_INSENSITIVE);

    /** Pattern to match workflow execution summary */
    private static final Pattern SUMMARY_PATTERN = Pattern.compile(
            "Completed at:\\s+(.+)|Duration:\\s+(.+)|CPU hours:\\s+(.+)|Succeeded:\\s+(\\d+)");

    private final AtomicBoolean cancelled = new AtomicBoolean(false);

    /**
     * Monitors the execution of a process and provides progress updates.
     *
     * @param process the process to monitor
     * @param progressListener optional progress listener for updates
     * @return CompletableFuture containing the execution result
     */
    public CompletableFuture<ExecutionResult> monitorExecution(Process process, ProgressListener progressListener) {
        return CompletableFuture.supplyAsync(() -> {
            StringBuilder outputBuilder = new StringBuilder();
            StringBuilder errorBuilder = new StringBuilder();

            try {
                // Monitor stdout and stderr concurrently
                CompletableFuture<Void> stdoutMonitor = monitorStream(
                        process.getInputStream(), outputBuilder, progressListener, false);

                CompletableFuture<Void> stderrMonitor = monitorStream(
                        process.getErrorStream(), errorBuilder, progressListener, true);

                // Wait for process completion or cancellation
                while (process.isAlive() && !cancelled.get()) {
                    try {
                        Thread.sleep(100);
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                        break;
                    }
                }

                if (cancelled.get()) {
                    process.destroyForcibly();
                    throw new RuntimeException("Execution was cancelled");
                }

                // Wait for streams to finish reading
                stdoutMonitor.join();
                stderrMonitor.join();

                int exitCode = process.waitFor();
                String output = outputBuilder.toString();
                String error = errorBuilder.toString();

                logger.info("Process completed with exit code: " + exitCode);

                if (progressListener != null) {
                    if (exitCode == 0) {
                        progressListener.setMessage("Workflow completed successfully");
                        progressListener.setProgress(1.0);
                    } else {
                        progressListener.setMessage("Workflow failed with exit code: " + exitCode);
                    }
                }

                return new ExecutionResult("nextflow-command", exitCode, output, error,
                        java.time.LocalDateTime.now().minusSeconds(1), java.time.LocalDateTime.now());

            } catch (Exception e) {
                logger.log(Level.SEVERE, "Error monitoring process execution", e);
                return new ExecutionResult("nextflow-command", -1, outputBuilder.toString(),
                        errorBuilder.toString() + "\nMonitoring error: " + e.getMessage(),
                        java.time.LocalDateTime.now().minusSeconds(1), java.time.LocalDateTime.now());
            }
        });
    }

    /**
     * Requests cancellation of the currently monitored execution.
     */
    public void requestCancellation() {
        cancelled.set(true);
        logger.info("Execution cancellation requested");
    }

    /**
     * Checks if a cancellation has been requested.
     *
     * @return true if cancellation was requested, false otherwise
     */
    public boolean isCancellationRequested() {
        return cancelled.get();
    }

    /**
     * Resets the cancellation state for reuse of this monitor.
     */
    public void reset() {
        cancelled.set(false);
    }

    /**
     * Parses Nextflow output to extract progress information.
     *
     * @param line the output line to parse
     * @return progress value between 0.0 and 1.0, or -1 if no progress found
     */
    public double parseProgress(String line) {
        if (line == null || line.isEmpty()) {
            return -1;
        }

        Matcher matcher = PROGRESS_PATTERN.matcher(line);
        if (matcher.find()) {
            try {
                int completed = Integer.parseInt(matcher.group(3));
                int total = Integer.parseInt(matcher.group(4));

                if (total > 0) {
                    double progress = (double) completed / total;
                    logger.fine("Parsed progress: " + completed + "/" + total + " = " + progress);
                    return progress;
                }
            } catch (NumberFormatException e) {
                logger.warning("Failed to parse progress numbers from: " + line);
            }
        }

        // Check for completion
        if (COMPLETION_PATTERN.matcher(line).find()) {
            logger.info("Detected workflow completion");
            return 1.0;
        }

        return -1;
    }

    /**
     * Checks if a line contains an error message.
     *
     * @param line the line to check
     * @return true if the line contains an error, false otherwise
     */
    public boolean isErrorLine(String line) {
        if (line == null || line.isEmpty()) {
            return false;
        }

        return ERROR_PATTERN.matcher(line).find();
    }

    /**
     * Extracts the current process name from a Nextflow progress line.
     *
     * @param line the progress line to parse
     * @return the process name, or null if not found
     */
    public String extractProcessName(String line) {
        if (line == null || line.isEmpty()) {
            return null;
        }

        Matcher matcher = PROGRESS_PATTERN.matcher(line);
        if (matcher.find()) {
            return matcher.group(2);
        }

        return null;
    }

    /**
     * Checks if a line indicates workflow completion.
     *
     * @param line the line to check
     * @return true if the line indicates completion, false otherwise
     */
    public boolean isCompletionLine(String line) {
        if (line == null || line.isEmpty()) {
            return false;
        }

        return COMPLETION_PATTERN.matcher(line).find() ||
               line.contains("Pipeline completed") ||
               SUMMARY_PATTERN.matcher(line).find();
    }

    /**
     * Monitors an input stream and processes output lines.
     *
     * @param inputStream the input stream to monitor
     * @param outputBuilder builder to accumulate output
     * @param progressListener optional progress listener
     * @param isErrorStream true if this is an error stream
     * @return CompletableFuture for the monitoring task
     */
    private CompletableFuture<Void> monitorStream(
            java.io.InputStream inputStream,
            StringBuilder outputBuilder,
            ProgressListener progressListener,
            boolean isErrorStream) {

        return CompletableFuture.runAsync(() -> {
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream))) {
                String line;
                while ((line = reader.readLine()) != null && !cancelled.get()) {
                    outputBuilder.append(line).append("\n");

                    if (isErrorStream) {
                        if (isErrorLine(line)) {
                            logger.warning("Error detected: " + line);
                        }
                    } else {
                        processOutputLine(line, progressListener);
                    }
                }
            } catch (IOException e) {
                logger.log(Level.WARNING, "Error reading from " +
                          (isErrorStream ? "error" : "output") + " stream", e);
            }
        });
    }

    /**
     * Processes a single output line for progress and status updates.
     *
     * @param line the output line to process
     * @param progressListener optional progress listener
     */
    private void processOutputLine(String line, ProgressListener progressListener) {
        if (progressListener == null) {
            return;
        }

        // Update progress if found
        double progress = parseProgress(line);
        if (progress >= 0) {
            progressListener.setProgress(progress);
        }

        // Update message based on content
        String processName = extractProcessName(line);
        if (processName != null) {
            progressListener.setMessage("Running process: " + processName);
        } else if (isCompletionLine(line)) {
            progressListener.setMessage("Workflow completed");
        } else if (line.contains("Launching")) {
            progressListener.setMessage("Launching workflow...");
        } else if (line.contains("Pulling")) {
            progressListener.setMessage("Pulling container images...");
        } else if (line.contains("executor")) {
            progressListener.setMessage("Initializing executor...");
        }

        logger.fine("Processed output line: " + line);
    }
}