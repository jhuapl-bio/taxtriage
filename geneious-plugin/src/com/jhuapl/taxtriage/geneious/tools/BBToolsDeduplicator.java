package com.jhuapl.taxtriage.geneious.tools;

import com.jhuapl.taxtriage.geneious.docker.DockerException;
import com.jhuapl.taxtriage.geneious.docker.DockerManager;
import com.jhuapl.taxtriage.geneious.docker.ExecutionResult;
import jebl.util.ProgressListener;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Handles read deduplication and preprocessing using BBTools suite via Docker.
 *
 * <p>This class implements the complete BBTools preprocessing pipeline for FASTQ files:
 * <ol>
 *   <li><strong>Split interleaved FASTQ:</strong> Separate R1 and R2 reads from interleaved format</li>
 *   <li><strong>Deduplicate reads:</strong> Remove duplicate reads with configurable similarity threshold</li>
 *   <li><strong>Clean read names:</strong> Standardize read naming for downstream compatibility</li>
 *   <li><strong>Re-interleave:</strong> Combine R1 and R2 back into interleaved format</li>
 * </ol>
 *
 * <h3>BBTools Commands Used:</h3>
 * <ul>
 *   <li><strong>reformat.sh:</strong> Split and re-interleave FASTQ files</li>
 *   <li><strong>clumpify.sh:</strong> Remove duplicate reads with sequence similarity</li>
 *   <li><strong>rename.sh:</strong> Clean and standardize read names</li>
 * </ul>
 *
 * <h3>Docker Integration:</h3>
 * <ul>
 *   <li><strong>Container:</strong> staphb/bbtools:39.01</li>
 *   <li><strong>Volume mounting:</strong> Input/output directories mounted for processing</li>
 *   <li><strong>Progress tracking:</strong> Real-time monitoring of processing steps</li>
 *   <li><strong>Error handling:</strong> Graceful cleanup and detailed error reporting</li>
 * </ul>
 *
 * @author TaxTriage Development Team
 * @version 2.0
 * @since 2.0
 */
public class BBToolsDeduplicator {

    private static final Logger logger = Logger.getLogger(BBToolsDeduplicator.class.getName());

    /** BBTools Docker image */
    private static final String BBTOOLS_IMAGE = "staphb/bbtools:39.01";

    /** Default timeout for each BBTools operation in minutes */
    private static final int DEFAULT_TIMEOUT_MINUTES = 30;

    /** Default substitution threshold for clumpify deduplication */
    private static final int DEFAULT_SUBS_THRESHOLD = 5;

    /**
     * Result of a BBTools deduplication operation.
     */
    public static class DeduplicationResult {
        private final boolean success;
        private final String outputPath;
        private final String errorMessage;
        private final long duplicatesRemoved;
        private final long totalReads;
        private final List<String> stepResults;

        public DeduplicationResult(boolean success, String outputPath, String errorMessage,
                                 long duplicatesRemoved, long totalReads, List<String> stepResults) {
            this.success = success;
            this.outputPath = outputPath;
            this.errorMessage = errorMessage;
            this.duplicatesRemoved = duplicatesRemoved;
            this.totalReads = totalReads;
            this.stepResults = new ArrayList<>(stepResults != null ? stepResults : new ArrayList<>());
        }

        public boolean isSuccess() { return success; }
        public String getOutputPath() { return outputPath; }
        public String getErrorMessage() { return errorMessage; }
        public long getDuplicatesRemoved() { return duplicatesRemoved; }
        public long getTotalReads() { return totalReads; }
        public List<String> getStepResults() { return stepResults; }

        /**
         * Gets the deduplication percentage.
         * @return percentage of reads removed as duplicates
         */
        public double getDeduplicationPercentage() {
            if (totalReads == 0) return 0.0;
            return (double) duplicatesRemoved / totalReads * 100.0;
        }

        @Override
        public String toString() {
            return String.format("BBToolsDeduplicationResult{success=%s, output='%s', " +
                               "duplicatesRemoved=%d, totalReads=%d (%.2f%% removed)}",
                               success, outputPath, duplicatesRemoved, totalReads, getDeduplicationPercentage());
        }
    }

    private final DockerManager dockerManager;
    private final int substitutionThreshold;

    /**
     * Creates a new BBToolsDeduplicator instance with default settings.
     * @throws DockerException if Docker is not available
     */
    public BBToolsDeduplicator() throws DockerException {
        this(DEFAULT_SUBS_THRESHOLD);
    }

    /**
     * Creates a new BBToolsDeduplicator instance with custom substitution threshold.
     * @param subsThreshold substitution threshold for clumpify deduplication (0-10)
     * @throws DockerException if Docker is not available
     */
    public BBToolsDeduplicator(int subsThreshold) throws DockerException {
        this.dockerManager = new DockerManager(BBTOOLS_IMAGE);
        this.substitutionThreshold = Math.max(0, Math.min(10, subsThreshold)); // Clamp to 0-10
        validateBBToolsAvailability();
    }

    /**
     * Validates that BBTools Docker image is available.
     * @throws DockerException if BBTools image is not accessible
     */
    private void validateBBToolsAvailability() throws DockerException {
        if (!dockerManager.isDockerAvailable()) {
            throw new DockerException("Docker is not available for BBTools operations");
        }

        logger.info("BBTools Docker integration validated successfully");
    }

    /**
     * Deduplicates reads in an interleaved FASTQ file using the complete BBTools pipeline.
     *
     * @param inputFile Path to the input interleaved FASTQ file
     * @param outputDir Directory where deduplicated FASTQ will be saved
     * @param progressListener Optional progress listener for monitoring
     * @return DeduplicationResult containing the path to deduplicated FASTQ or error message
     */
    public DeduplicationResult deduplicateReads(Path inputFile, Path outputDir, ProgressListener progressListener) {
        return deduplicateReads(inputFile, outputDir, this.substitutionThreshold, progressListener);
    }

    /**
     * Deduplicates reads in an interleaved FASTQ file using the complete BBTools pipeline with custom threshold.
     *
     * @param inputFile Path to the input interleaved FASTQ file
     * @param outputDir Directory where deduplicated FASTQ will be saved
     * @param subsThreshold Substitution threshold for clumpify (default: 5)
     * @param progressListener Optional progress listener for monitoring
     * @return DeduplicationResult containing the path to deduplicated FASTQ or error message
     */
    public DeduplicationResult deduplicateReads(Path inputFile, Path outputDir, int subsThreshold,
                                               ProgressListener progressListener) {
        logger.info("==========================================");
        logger.info("BBTOOLS DEDUPLICATION: Starting read deduplication pipeline");
        logger.info("  Input: " + inputFile);
        logger.info("  Output dir: " + outputDir);
        logger.info("  Substitution threshold: " + subsThreshold);
        logger.info("  Docker image: " + BBTOOLS_IMAGE);
        logger.info("==========================================");

        List<String> stepResults = new ArrayList<>();

        // Validate inputs
        if (!Files.exists(inputFile)) {
            String error = "Input FASTQ file does not exist: " + inputFile;
            logger.warning(error);
            return new DeduplicationResult(false, null, error, 0, 0, stepResults);
        }

        // Create output directory
        try {
            Files.createDirectories(outputDir);
        } catch (IOException e) {
            String error = "Failed to create output directory: " + e.getMessage();
            logger.log(Level.WARNING, error, e);
            return new DeduplicationResult(false, null, error, 0, 0, stepResults);
        }

        // Generate file names
        String baseName = getBaseName(inputFile.getFileName().toString());
        Path workDir = outputDir.resolve("bbtools_work");

        try {
            Files.createDirectories(workDir);

            // Copy input file to workDir so all files are in the same directory for Docker volume mount
            Path inputInWorkDir = workDir.resolve(inputFile.getFileName());
            Files.copy(inputFile, inputInWorkDir, java.nio.file.StandardCopyOption.REPLACE_EXISTING);
            logger.info("Copied input file to work directory: " + inputInWorkDir);

            // Define intermediate and final file paths (all in workDir)
            Path splitR1 = workDir.resolve(baseName + ".R1.fq.gz");
            Path splitR2 = workDir.resolve(baseName + ".R2.fq.gz");
            Path dedupR1 = workDir.resolve(baseName + "_dedupe.R1.fq.gz");
            Path dedupR2 = workDir.resolve(baseName + "_dedupe.R2.fq.gz");
            Path finalOutput = workDir.resolve(baseName + "_dedupe.fq.gz");

            // Progress tracking
            if (progressListener != null) {
                progressListener.setMessage("Starting BBTools deduplication pipeline...");
                progressListener.setProgress(0.1);
            }

            // Step 1: Split interleaved FASTQ
            logger.info("\n[Step 1/4] Splitting interleaved FASTQ...");
            ExecutionResult step1Result = executeReformat(
                "Split interleaved reads",
                inputInWorkDir, splitR1, splitR2, workDir, progressListener
            );

            if (!step1Result.isSuccessful()) {
                String error = "Failed to split interleaved FASTQ: " + step1Result.getErrorOutput();
                logger.warning(error);
                return new DeduplicationResult(false, null, error, 0, 0, stepResults);
            }
            stepResults.add("Split: " + step1Result.getSummary());
            logger.info("  ✓ Splitting complete: " + splitR1.getFileName() + ", " + splitR2.getFileName());

            if (progressListener != null) {
                progressListener.setProgress(0.3);
                progressListener.setMessage("Deduplicating reads...");
            }

            // Step 2: Deduplicate with clumpify
            logger.info("\n[Step 2/4] Deduplicating reads with clumpify...");
            ExecutionResult step2Result = executeClumpify(
                splitR1, splitR2, dedupR1, dedupR2, subsThreshold, workDir, progressListener
            );

            if (!step2Result.isSuccessful()) {
                String error = "Failed to deduplicate reads: " + step2Result.getErrorOutput();
                logger.warning(error);
                return new DeduplicationResult(false, null, error, 0, 0, stepResults);
            }
            stepResults.add("Deduplicate: " + step2Result.getSummary());

            // Parse deduplication statistics
            long duplicatesRemoved = parseDuplicateCount(step2Result.getStandardOutput());
            long totalReads = parseReadCount(step2Result.getStandardOutput());

            logger.info("  ✓ Deduplication complete: " + duplicatesRemoved + " duplicates removed");

            if (progressListener != null) {
                progressListener.setProgress(0.6);
                progressListener.setMessage("Re-interleaving reads...");
            }

            // Step 3: Re-interleave with underscore replacement for copy counts
            logger.info("\n[Step 3/3] Re-interleaving reads and formatting names...");
            ExecutionResult step3Result = executeReformatInterleaveWithUnderscore(
                "Re-interleave deduplicated reads",
                dedupR1, dedupR2, finalOutput, workDir, progressListener
            );

            if (!step3Result.isSuccessful()) {
                String error = "Failed to re-interleave reads: " + step3Result.getErrorOutput();
                logger.warning(error);
                return new DeduplicationResult(false, null, error, 0, 0, stepResults);
            }
            stepResults.add("Re-interleave: " + step3Result.getSummary());
            logger.info("  ✓ Re-interleaving complete: " + finalOutput.getFileName());

            if (progressListener != null) {
                progressListener.setProgress(1.0);
                progressListener.setMessage("BBTools deduplication complete");
            }

            // Log summary
            logger.info("\n==========================================");
            logger.info("BBTOOLS DEDUPLICATION COMPLETE");
            logger.info("  Input: " + inputFile.getFileName());
            logger.info("  Output: " + finalOutput);
            logger.info("  Total reads processed: " + totalReads);
            logger.info("  Duplicates removed: " + duplicatesRemoved);
            logger.info("  Deduplication rate: " + String.format("%.2f%%",
                       totalReads > 0 ? (double) duplicatesRemoved / totalReads * 100.0 : 0.0));
            logger.info("  Working directory: " + workDir);
            logger.info("==========================================\n");

            return new DeduplicationResult(true, finalOutput.toString(), null,
                                         duplicatesRemoved, totalReads, stepResults);

        } catch (Exception e) {
            logger.log(Level.SEVERE, "BBTools deduplication failed with exception", e);
            String error = "Deduplication failed: " + e.getMessage();
            return new DeduplicationResult(false, null, error, 0, 0, stepResults);
        }
    }

    /**
     * Executes reformat.sh for splitting or re-interleaving FASTQ files.
     */
    private ExecutionResult executeReformat(String operation, Path input, Path output1, Path output2,
                                          Path workDir, ProgressListener progressListener) throws DockerException {
        List<String> command = new ArrayList<>();
        command.add("reformat.sh");
        command.add("in=/data/" + input.getFileName());
        command.add("out=/data/" + output1.getFileName());
        if (output2 != null) {
            command.add("out2=/data/" + output2.getFileName());
        }
        command.add("ow=t"); // Overwrite existing files

        return executeDockerCommandDirect(operation, command, workDir, progressListener);
    }

    /**
     * Executes reformat.sh for re-interleaving two FASTQ files into one.
     */
    private ExecutionResult executeReformatInterleave(String operation, Path input1, Path input2, Path output,
                                                     Path workDir, ProgressListener progressListener) throws DockerException {
        List<String> command = new ArrayList<>();
        command.add("reformat.sh");
        command.add("in=/data/" + input1.getFileName());
        command.add("in2=/data/" + input2.getFileName());
        command.add("out=/data/" + output.getFileName());
        command.add("ow=t"); // Overwrite existing files

        return executeDockerCommandDirect(operation, command, workDir, progressListener);
    }

    /**
     * Executes reformat.sh for re-interleaving two FASTQ files into one with underscore replacement.
     * This replaces spaces in read names with underscores, converting " copies=X" to "_copies=X".
     */
    private ExecutionResult executeReformatInterleaveWithUnderscore(String operation, Path input1, Path input2, Path output,
                                                                   Path workDir, ProgressListener progressListener) throws DockerException {
        List<String> command = new ArrayList<>();
        command.add("reformat.sh");
        command.add("in=/data/" + input1.getFileName());
        command.add("in2=/data/" + input2.getFileName());
        command.add("out=/data/" + output.getFileName());
        command.add("underscore=t"); // Replace spaces with underscores in read names
        command.add("ow=t"); // Overwrite existing files

        return executeDockerCommandDirect(operation, command, workDir, progressListener);
    }

    /**
     * Executes clumpify.sh for read deduplication.
     */
    private ExecutionResult executeClumpify(Path input1, Path input2, Path output1, Path output2,
                                          int subsThreshold, Path workDir, ProgressListener progressListener) throws DockerException {
        List<String> command = new ArrayList<>();
        command.add("clumpify.sh");
        command.add("in=/data/" + input1.getFileName());
        command.add("in2=/data/" + input2.getFileName());
        command.add("out=/data/" + output1.getFileName());
        command.add("out2=/data/" + output2.getFileName());
        command.add("dedupe");
        command.add("subs=" + subsThreshold);
        command.add("addcount=t"); // Add count information to read names
        command.add("ow=t"); // Overwrite existing files

        return executeDockerCommandDirect("Deduplicate reads", command, workDir, progressListener);
    }

    /**
     * Executes rename.sh for cleaning read names.
     */
    private ExecutionResult executeRename(Path input, Path output, Path workDir, ProgressListener progressListener) throws DockerException {
        List<String> command = new ArrayList<>();
        command.add("rename.sh");
        command.add("in=/data/" + input.getFileName());
        command.add("out=/data/" + output.getFileName());
        command.add("ow=t"); // Overwrite existing files

        return executeDockerCommandDirect("Clean read names", command, workDir, progressListener);
    }

    /**
     * Executes a BBTools command directly in Docker container with proper volume mounting.
     */
    private ExecutionResult executeDockerCommandDirect(String operation, List<String> command, Path workDir,
                                                      ProgressListener progressListener) throws DockerException {
        logger.info("  Executing: " + operation);
        logger.info("  Command: " + String.join(" ", command));

        try {
            // Build Docker command with proper volume mounting
            List<String> dockerCmd = new ArrayList<>();
            dockerCmd.add("docker");
            dockerCmd.add("run");
            dockerCmd.add("--rm");
            dockerCmd.add("-v");
            dockerCmd.add(workDir.toAbsolutePath() + ":/data");
            dockerCmd.add(BBTOOLS_IMAGE);
            dockerCmd.addAll(command);

            logger.info("  Docker command: " + String.join(" ", dockerCmd));

            // Execute the command
            ProcessBuilder pb = new ProcessBuilder(dockerCmd);
            pb.redirectErrorStream(false);
            Process process = pb.start();

            // Capture output
            StringBuilder stdout = new StringBuilder();
            StringBuilder stderr = new StringBuilder();

            try (BufferedReader stdoutReader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                 BufferedReader stderrReader = new BufferedReader(new InputStreamReader(process.getErrorStream()))) {

                String line;
                while ((line = stdoutReader.readLine()) != null) {
                    stdout.append(line).append("\n");
                    logger.fine(line);
                }

                while ((line = stderrReader.readLine()) != null) {
                    stderr.append(line).append("\n");
                    logger.fine(line);
                }
            }

            // Wait for completion
            boolean finished = process.waitFor(DEFAULT_TIMEOUT_MINUTES, TimeUnit.MINUTES);
            if (!finished) {
                process.destroy();
                throw new DockerException("BBTools command timed out after " + DEFAULT_TIMEOUT_MINUTES + " minutes");
            }

            int exitCode = process.exitValue();
            return new com.jhuapl.taxtriage.geneious.docker.ExecutionResult(
                String.join(" ", dockerCmd), exitCode, stdout.toString(), stderr.toString(),
                java.time.LocalDateTime.now(), java.time.LocalDateTime.now()
            );

        } catch (Exception e) {
            logger.log(Level.SEVERE, "Failed to execute BBTools command", e);
            throw new DockerException("Failed to execute BBTools command: " + e.getMessage(), e);
        }
    }

    /**
     * Parses the duplicate count from clumpify output.
     */
    private long parseDuplicateCount(String output) {
        if (output == null) return 0;

        // Look for patterns like "Duplicates found: 12345"
        Pattern pattern = Pattern.compile("(?i)duplicates?\\s+(?:found|removed):\\s*(\\d+)");
        Matcher matcher = pattern.matcher(output);
        if (matcher.find()) {
            try {
                return Long.parseLong(matcher.group(1));
            } catch (NumberFormatException e) {
                logger.fine("Could not parse duplicate count: " + matcher.group(1));
            }
        }

        // Alternative pattern: "Removed X duplicates"
        pattern = Pattern.compile("(?i)removed\\s+(\\d+)\\s+duplicates?");
        matcher = pattern.matcher(output);
        if (matcher.find()) {
            try {
                return Long.parseLong(matcher.group(1));
            } catch (NumberFormatException e) {
                logger.fine("Could not parse removed duplicates count: " + matcher.group(1));
            }
        }

        logger.fine("No duplicate count pattern found in output");
        return 0;
    }

    /**
     * Parses the total read count from clumpify output.
     */
    private long parseReadCount(String output) {
        if (output == null) return 0;

        // Look for patterns like "Reads processed: 67890"
        Pattern pattern = Pattern.compile("(?i)reads?\\s+(?:processed|input):\\s*(\\d+)");
        Matcher matcher = pattern.matcher(output);
        if (matcher.find()) {
            try {
                return Long.parseLong(matcher.group(1));
            } catch (NumberFormatException e) {
                logger.fine("Could not parse read count: " + matcher.group(1));
            }
        }

        // Alternative pattern: "Input reads: X"
        pattern = Pattern.compile("(?i)input\\s+reads?:\\s*(\\d+)");
        matcher = pattern.matcher(output);
        if (matcher.find()) {
            try {
                return Long.parseLong(matcher.group(1));
            } catch (NumberFormatException e) {
                logger.fine("Could not parse input reads count: " + matcher.group(1));
            }
        }

        logger.fine("No read count pattern found in output");
        return 0;
    }

    /**
     * Gets the base name of a file without extensions.
     */
    private String getBaseName(String filename) {
        if (filename == null || filename.isEmpty()) {
            return "processed";
        }

        // Remove common extensions
        String baseName = filename;
        if (baseName.toLowerCase().endsWith(".gz")) {
            baseName = baseName.substring(0, baseName.length() - 3);
        }
        if (baseName.toLowerCase().endsWith(".fq") || baseName.toLowerCase().endsWith(".fastq")) {
            int lastDot = baseName.lastIndexOf('.');
            if (lastDot > 0) {
                baseName = baseName.substring(0, lastDot);
            }
        }

        return baseName;
    }

    /**
     * Checks if BBTools is available via Docker.
     * @return true if BBTools Docker image is accessible
     */
    public boolean isBBToolsAvailable() {
        try {
            return dockerManager.isDockerAvailable() && dockerManager.isImageAvailable();
        } catch (Exception e) {
            logger.log(Level.FINE, "BBTools availability check failed", e);
            return false;
        }
    }

    /**
     * Gets the BBTools Docker image name being used.
     * @return the Docker image name
     */
    public String getImageName() {
        return BBTOOLS_IMAGE;
    }

    /**
     * Gets the substitution threshold being used for deduplication.
     * @return the substitution threshold (0-10)
     */
    public int getSubstitutionThreshold() {
        return substitutionThreshold;
    }
}