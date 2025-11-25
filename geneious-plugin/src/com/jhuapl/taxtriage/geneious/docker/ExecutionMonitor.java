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

    /** Pattern to match Nextflow progress lines (old format with "process >") */
    private static final Pattern PROGRESS_PATTERN = Pattern.compile(
            "\\[([^\\]]+)\\]\\s+process\\s+>\\s+([^\\s]+)\\s+\\[(\\d+)\\s+of\\s+(\\d+)\\]");

    /** Pattern to match Nextflow process submission (executor log) */
    private static final Pattern EXECUTOR_PATTERN = Pattern.compile(
            "\\[([^\\]]+)\\]\\s+Submitted\\s+process\\s+>\\s+([^\\s(]+)");

    /** Pattern to match Nextflow process start */
    private static final Pattern PROCESS_START_PATTERN = Pattern.compile(
            "executor\\s+>\\s+([^\\s(]+).*started");

    /** Pattern to match modern Nextflow process lines with [hash] PREFIX:PROCESS_NAME (sample) [%] X of Y format */
    private static final Pattern MODERN_PROCESS_PATTERN = Pattern.compile(
            "\\[([^\\]]+)\\]\\s+([A-Z0-9_:…]+)\\s+\\(([^)]+)\\)\\s+\\[([^\\]]+)\\]\\s+(\\d+)\\s+of\\s+(\\d+)");

    /** Pattern to match executor summary lines like "executor >  local (26)" */
    private static final Pattern EXECUTOR_SUMMARY_PATTERN = Pattern.compile(
            "executor\\s+>\\s+\\w+\\s+\\((\\d+)\\)");

    /** Pattern to match Nextflow completion status */
    private static final Pattern COMPLETION_PATTERN = Pattern.compile(
            "(Pipeline completed|Completed) at:\\s*(.+)");

    /** Pattern to match Nextflow error messages */
    private static final Pattern ERROR_PATTERN = Pattern.compile(
            "(ERROR|WARN|Exception|Failed|Error).*", Pattern.CASE_INSENSITIVE);

    /** Pattern to match workflow execution summary */
    private static final Pattern SUMMARY_PATTERN = Pattern.compile(
            "Completed at:\\s+(.+)|Duration:\\s+(.+)|CPU hours:\\s+(.+)|Succeeded:\\s+(\\d+)");

    private String currentProcessName = null;
    private String lastReportedProcess = null;

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
            System.out.println("========================================");
            System.out.println("[TaxTriage] ExecutionMonitor starting...");
            System.out.println("[TaxTriage] ALL Nextflow output will be displayed below");
            System.out.println("========================================");
            System.out.flush();

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

                System.out.println("========================================");
                System.out.println("[TaxTriage] Process completed with exit code: " + exitCode);
                System.out.println("========================================");
                System.out.flush();

                logger.info("Process completed with exit code: " + exitCode);

                if (progressListener != null) {
                    if (exitCode == 0) {
                        progressListener.setMessage("Workflow completed successfully");
                        progressListener.setProgress(1.0);
                        System.out.println("[TaxTriage Progress] >>> Workflow completed successfully");
                    } else {
                        String msg = "Workflow failed with exit code: " + exitCode;
                        progressListener.setMessage(msg);
                        System.out.println("[TaxTriage Progress] >>> " + msg);
                    }
                    System.out.flush();
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
        resetProcessTracking();
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

        // Try modern format first: [hash] PROCESS (sample) [100%] 1 of 1
        Matcher modernMatcher = MODERN_PROCESS_PATTERN.matcher(line);
        if (modernMatcher.find()) {
            try {
                int completed = Integer.parseInt(modernMatcher.group(5));
                int total = Integer.parseInt(modernMatcher.group(6));

                if (total > 0) {
                    double progress = (double) completed / total;
                    logger.fine("Parsed modern progress: " + completed + "/" + total + " = " + progress);
                    return progress;
                }
            } catch (NumberFormatException e) {
                logger.warning("Failed to parse modern progress numbers from: " + line);
            }
        }

        // Try old format: [hash] process > NAME [1 of 4]
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

        // Try modern format first: [hash] NFCORE_TAXTRIAGE:PROCESS_NAME (sample) [%] X of Y
        Matcher modernMatcher = MODERN_PROCESS_PATTERN.matcher(line);
        if (modernMatcher.find()) {
            String processName = modernMatcher.group(2); // Full process path
            currentProcessName = processName;
            return processName;
        }

        // Try to match standard progress pattern
        Matcher matcher = PROGRESS_PATTERN.matcher(line);
        if (matcher.find()) {
            String processName = matcher.group(2);
            currentProcessName = processName;
            return processName;
        }

        // Try to match executor submission pattern
        matcher = EXECUTOR_PATTERN.matcher(line);
        if (matcher.find()) {
            String processName = matcher.group(2);
            currentProcessName = processName;
            return processName;
        }

        // Try to match process start pattern
        matcher = PROCESS_START_PATTERN.matcher(line);
        if (matcher.find()) {
            String processName = matcher.group(1);
            currentProcessName = processName;
            return processName;
        }

        return null;
    }

    /**
     * Gets the current process name being executed.
     *
     * @return the current process name, or null if none
     */
    public String getCurrentProcessName() {
        return currentProcessName;
    }

    /**
     * Resets the current process tracking state.
     */
    public void resetProcessTracking() {
        currentProcessName = null;
        lastReportedProcess = null;
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
            String streamType = isErrorStream ? "STDERR" : "STDOUT";
            System.out.println("[TaxTriage] Starting to monitor " + streamType + "...");
            System.out.flush();

            try (BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream))) {
                String line;
                int lineCount = 0;
                while ((line = reader.readLine()) != null && !cancelled.get()) {
                    lineCount++;
                    outputBuilder.append(line).append("\n");

                    if (isErrorStream) {
                        // Always print error stream to console
                        System.err.println("[TaxTriage Error] " + line);
                        System.err.flush();
                        if (isErrorLine(line)) {
                            logger.warning("Error detected: " + line);
                        }
                    } else {
                        processOutputLine(line, progressListener);
                    }
                }
                System.out.println("[TaxTriage] Finished monitoring " + streamType + " (" + lineCount + " lines)");
                System.out.flush();
            } catch (IOException e) {
                System.err.println("[TaxTriage] ERROR reading from " + streamType + ": " + e.getMessage());
                System.err.flush();
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
        // ALWAYS print ALL output to console for complete visibility
        // This ensures we see everything Nextflow produces
        System.out.println("[TaxTriage Nextflow] " + line);
        System.out.flush(); // Ensure immediate output

        if (progressListener == null) {
            return;
        }

        // Update progress if found
        double progress = parseProgress(line);
        if (progress >= 0) {
            progressListener.setProgress(progress);
            System.out.println("[TaxTriage Progress] Progress: " + (int)(progress * 100) + "%");
            System.out.flush();
        }

        // Update message based on content
        String processName = extractProcessName(line);
        if (processName != null && !processName.equals(lastReportedProcess)) {
            // Format the process name for better readability
            String formattedName = formatProcessName(processName);
            String message = "TaxTriage: " + formattedName;
            progressListener.setMessage(message);
            lastReportedProcess = processName;
            logger.info("Started Nextflow process: " + processName);
            System.out.println("[TaxTriage Progress] >>> " + message);
            System.out.flush();
        } else if (isCompletionLine(line)) {
            String message = "TaxTriage: Workflow completed";
            progressListener.setMessage(message);
            System.out.println("[TaxTriage Progress] >>> " + message);
            System.out.flush();
        } else if (line.contains("Launching")) {
            String message = "TaxTriage: Launching workflow...";
            progressListener.setMessage(message);
            System.out.println("[TaxTriage Progress] >>> " + message);
            System.out.flush();
        } else if (line.contains("Pulling")) {
            String message = "TaxTriage: Pulling container images...";
            progressListener.setMessage(message);
            System.out.println("[TaxTriage Progress] >>> " + message);
            System.out.flush();
        } else if (line.contains("Staging foreign files")) {
            String message = "TaxTriage: Staging input files...";
            progressListener.setMessage(message);
            System.out.println("[TaxTriage Progress] >>> " + message);
            System.out.flush();
        } else if (line.contains("executor") && !line.contains("process")) {
            String message = "TaxTriage: Initializing executor...";
            progressListener.setMessage(message);
            System.out.println("[TaxTriage Progress] >>> " + message);
            System.out.flush();
        }

        logger.fine("Processed output line: " + line);
    }

    /**
     * Formats a Nextflow process name for better human readability.
     * Converts technical process names like "FASTQC" or "KRAKEN2_KRAKEN2"
     * into more readable forms like "FastQC" or "Kraken2".
     *
     * @param processName the raw process name from Nextflow
     * @return formatted process name
     */
    private String formatProcessName(String processName) {
        if (processName == null || processName.isEmpty()) {
            return processName;
        }

        // Remove common prefixes if they appear (e.g., "NFCORE_TAXTRIAGE:")
        String name = processName;
        if (name.contains(":")) {
            String[] parts = name.split(":");
            name = parts[parts.length - 1]; // Take the last part
        }

        // Handle duplicated names (e.g., "KRAKEN2_KRAKEN2" -> "KRAKEN2")
        String[] parts = name.split("_");
        if (parts.length == 2 && parts[0].equals(parts[1])) {
            name = parts[0];
        }

        // Convert from SCREAMING_SNAKE_CASE to Title Case
        String[] words = name.split("_");
        StringBuilder formatted = new StringBuilder();
        for (int i = 0; i < words.length; i++) {
            if (i > 0) {
                formatted.append(" ");
            }
            String word = words[i].toLowerCase();
            if (word.length() > 0) {
                // Keep known acronyms in uppercase
                if (isKnownAcronym(word)) {
                    formatted.append(word.toUpperCase());
                } else {
                    formatted.append(Character.toUpperCase(word.charAt(0)));
                    formatted.append(word.substring(1));
                }
            }
        }

        return formatted.toString();
    }

    /**
     * Checks if a word is a known bioinformatics acronym that should stay uppercase.
     *
     * @param word the word to check
     * @return true if it's a known acronym
     */
    private boolean isKnownAcronym(String word) {
        String upper = word.toUpperCase();
        return upper.equals("QC") || upper.equals("DNA") || upper.equals("RNA") ||
               upper.equals("BWA") || upper.equals("BAM") || upper.equals("SAM") ||
               upper.equals("VCF") || upper.equals("GFF") || upper.equals("GTF") ||
               upper.equals("NCBI") || upper.equals("BLAST");
    }
}