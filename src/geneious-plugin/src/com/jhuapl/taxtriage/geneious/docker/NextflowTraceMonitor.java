package com.jhuapl.taxtriage.geneious.docker;

import jebl.util.ProgressListener;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.logging.Logger;

/**
 * Monitors Nextflow execution by reading the execution_trace*.txt file in real-time.
 * This provides accurate progress tracking based on actual process completion.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class NextflowTraceMonitor implements Runnable {

    private static final Logger logger = Logger.getLogger(NextflowTraceMonitor.class.getName());

    private final Path outputDirectory;
    private final ProgressListener progressListener;
    private final AtomicBoolean running = new AtomicBoolean(true);
    private final Set<String> completedProcesses = new HashSet<>();
    private final Map<String, ProcessInfo> processMap = new LinkedHashMap<>();

    private long lastFilePosition = 0;
    private int totalProcesses = 0;
    private int completedCount = 0;
    private String currentProcess = null;

    /**
     * Information about a Nextflow process from the trace file.
     */
    private static class ProcessInfo {
        String taskId;
        String processName;
        String status;
        String tag;

        ProcessInfo(String taskId, String processName, String status, String tag) {
            this.taskId = taskId;
            this.processName = processName;
            this.status = status;
            this.tag = tag;
        }
    }

    public NextflowTraceMonitor(Path outputDirectory, ProgressListener progressListener) {
        this.outputDirectory = outputDirectory;
        this.progressListener = progressListener;
    }

    @Override
    public void run() {
        System.out.println("========================================");
        System.out.println("[TaxTriage] NextflowTraceMonitor started");
        System.out.println("[TaxTriage] Monitoring: " + outputDirectory.resolve("pipeline_info"));
        System.out.println("========================================");
        System.out.flush();

        try {
            while (running.get()) {
                try {
                    Path traceFile = findTraceFile();
                    if (traceFile != null && Files.exists(traceFile)) {
                        readTraceUpdates(traceFile);
                    }

                    // Poll every 500ms
                    Thread.sleep(500);
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    break;
                } catch (Exception e) {
                    logger.warning("Error reading trace file: " + e.getMessage());
                }
            }
        } finally {
            System.out.println("========================================");
            System.out.println("[TaxTriage] NextflowTraceMonitor stopped");
            System.out.println("[TaxTriage] Completed processes: " + completedCount + "/" + totalProcesses);
            System.out.println("========================================");
            System.out.flush();
        }
    }

    /**
     * Finds the execution trace file in the pipeline_info directory.
     */
    private Path findTraceFile() {
        try {
            Path pipelineInfo = outputDirectory.resolve("pipeline_info");
            if (!Files.exists(pipelineInfo)) {
                return null;
            }

            // Look for execution_trace*.txt files
            File[] traceFiles = pipelineInfo.toFile().listFiles((dir, name) ->
                name.startsWith("execution_trace") && name.endsWith(".txt"));

            if (traceFiles != null && traceFiles.length > 0) {
                // Use the most recent one
                Arrays.sort(traceFiles, Comparator.comparingLong(File::lastModified).reversed());
                return traceFiles[0].toPath();
            }
        } catch (Exception e) {
            logger.warning("Error finding trace file: " + e.getMessage());
        }
        return null;
    }

    /**
     * Reads new lines from the trace file since last read.
     */
    private void readTraceUpdates(Path traceFile) {
        try {
            long fileSize = Files.size(traceFile);
            if (fileSize < lastFilePosition) {
                // File was truncated/replaced, start over
                lastFilePosition = 0;
                completedProcesses.clear();
                processMap.clear();
            }

            if (fileSize == lastFilePosition) {
                // No new data
                return;
            }

            try (RandomAccessFile raf = new RandomAccessFile(traceFile.toFile(), "r")) {
                raf.seek(lastFilePosition);

                String line;
                boolean isFirstLine = (lastFilePosition == 0);

                while ((line = raf.readLine()) != null) {
                    if (isFirstLine) {
                        // Skip header line
                        isFirstLine = false;
                        continue;
                    }

                    processTraceLine(line);
                }

                lastFilePosition = raf.getFilePointer();
            }

        } catch (IOException e) {
            logger.warning("Error reading trace file: " + e.getMessage());
        }
    }

    /**
     * Processes a single line from the trace file.
     * Format: task_id hash native_id name status exit submit duration realtime %cpu %mem rss vmem ...
     */
    private void processTraceLine(String line) {
        if (line == null || line.trim().isEmpty()) {
            return;
        }

        try {
            String[] fields = line.split("\t");
            if (fields.length < 4) {
                return;
            }

            String taskId = fields[0].trim();
            String processName = fields[3].trim();
            String status = fields.length > 4 ? fields[4].trim() : "";
            String tag = fields.length > 13 ? fields[13].trim() : "";

            ProcessInfo info = new ProcessInfo(taskId, processName, status, tag);
            processMap.put(taskId, info);

            // Count total unique processes
            totalProcesses = processMap.size();

            // Check if completed
            if ("COMPLETED".equals(status) && !completedProcesses.contains(taskId)) {
                completedProcesses.add(taskId);
                completedCount++;

                // Update progress
                updateProgress(processName, tag);

                System.out.println("[TaxTriage Trace] Completed: " + processName +
                    (tag.isEmpty() ? "" : " (" + tag + ")") +
                    " - " + completedCount + "/" + totalProcesses);
                System.out.flush();
            } else if ("RUNNING".equals(status)) {
                if (currentProcess == null || !currentProcess.equals(processName)) {
                    currentProcess = processName;
                    updateProgress(processName, tag);

                    System.out.println("[TaxTriage Trace] Running: " + processName +
                        (tag.isEmpty() ? "" : " (" + tag + ")"));
                    System.out.flush();
                }
            }

        } catch (Exception e) {
            logger.warning("Error parsing trace line: " + e.getMessage());
        }
    }

    /**
     * Updates the progress listener with current process and percentage.
     */
    private void updateProgress(String processName, String tag) {
        if (progressListener == null) {
            return;
        }

        // Format process name for display
        String displayName = formatProcessName(processName);

        // Create progress message
        String message = "TaxTriage: " + displayName;
        if (tag != null && !tag.isEmpty() && !"-".equals(tag)) {
            message += " (" + tag + ")";
        }

        progressListener.setMessage(message);

        // Calculate progress percentage
        if (totalProcesses > 0) {
            double progress = (double) completedCount / totalProcesses;
            progressListener.setProgress(progress);

            System.out.println("[TaxTriage Progress] >>> " + message +
                " - " + (int)(progress * 100) + "%");
            System.out.flush();
        }
    }

    /**
     * Formats a Nextflow process name for better readability.
     * Converts "NFCORE_TAXTRIAGE:CLASSIFIER:KRAKEN2_KRAKEN2" to "Kraken2"
     */
    private String formatProcessName(String processName) {
        if (processName == null || processName.isEmpty()) {
            return "Processing";
        }

        // Remove prefix and get last component
        String name = processName;
        if (name.contains(":")) {
            String[] parts = name.split(":");
            name = parts[parts.length - 1];
        }

        // Handle duplicated names (e.g., "KRAKEN2_KRAKEN2" -> "KRAKEN2")
        String[] parts = name.split("_");
        if (parts.length == 2 && parts[0].equals(parts[1])) {
            name = parts[0];
        }

        // Convert to Title Case
        String[] words = name.split("_");
        StringBuilder formatted = new StringBuilder();
        for (int i = 0; i < words.length; i++) {
            if (i > 0) formatted.append(" ");
            String word = words[i].toLowerCase();
            if (word.length() > 0) {
                // Keep known acronyms uppercase
                if (isKnownAcronym(word)) {
                    formatted.append(word.toUpperCase());
                } else {
                    formatted.append(Character.toUpperCase(word.charAt(0)));
                    if (word.length() > 1) {
                        formatted.append(word.substring(1));
                    }
                }
            }
        }

        return formatted.toString();
    }

    /**
     * Checks if a word is a known bioinformatics acronym.
     */
    private boolean isKnownAcronym(String word) {
        String upper = word.toUpperCase();
        return upper.equals("QC") || upper.equals("DNA") || upper.equals("RNA") ||
               upper.equals("BWA") || upper.equals("BAM") || upper.equals("SAM") ||
               upper.equals("VCF") || upper.equals("GFF") || upper.equals("GTF") ||
               upper.equals("NCBI") || upper.equals("BLAST") || upper.equals("FASTQC");
    }

    /**
     * Stops the monitor.
     */
    public void stop() {
        running.set(false);
    }

    /**
     * Gets the current completion count.
     */
    public int getCompletedCount() {
        return completedCount;
    }

    /**
     * Gets the total process count.
     */
    public int getTotalProcesses() {
        return totalProcesses;
    }
}
