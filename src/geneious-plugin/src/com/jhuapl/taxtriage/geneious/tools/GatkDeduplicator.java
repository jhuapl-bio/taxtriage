package com.jhuapl.taxtriage.geneious.tools;

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

/**
 * Handles deduplication of mapped reads in BAM files using GATK MarkDuplicates via Docker.
 *
 * This class uses GATK MarkDuplicates for more sophisticated duplicate detection:
 * 1. Sort BAM by coordinate (samtools sort)
 * 2. Remove duplicates using GATK MarkDuplicates via Docker
 * 3. Index the deduplicated BAM (samtools index)
 *
 * GATK provides better handling of paired-end reads and generates detailed metrics.
 */
public class GatkDeduplicator {

    private static final Logger logger = Logger.getLogger(GatkDeduplicator.class.getName());
    private static final int PROCESS_TIMEOUT_MINUTES = 30;
    private static final String GATK_DOCKER_IMAGE = "broadinstitute/gatk:4.3.0.0";
    private static final String dockerCommand = findDockerCommand();

    /**
     * Result of a deduplication operation.
     */
    public static class DeduplicationResult {
        public final boolean success;
        public final String outputPath;
        public final String errorMessage;
        public final long duplicatesRemoved;

        public DeduplicationResult(boolean success, String outputPath, String errorMessage, long duplicatesRemoved) {
            this.success = success;
            this.outputPath = outputPath;
            this.errorMessage = errorMessage;
            this.duplicatesRemoved = duplicatesRemoved;
        }
    }

    /**
     * Check if Docker is available and GATK image is accessible.
     */
    public static boolean isGatkAvailable() {
        logger.info("=== GATK Availability Check ===");
        logger.info("Docker command: " + dockerCommand);

        try {
            // Check if Docker is running - use 'docker ps' to verify daemon is accessible
            logger.info("Checking Docker daemon with: " + dockerCommand + " ps");
            Process dockerCheck = new ProcessBuilder(dockerCommand, "ps")
                .redirectErrorStream(true)
                .start();

            boolean finished = dockerCheck.waitFor(5, TimeUnit.SECONDS);
            int exitCode = finished ? dockerCheck.exitValue() : -1;
            logger.info("Docker ps - finished: " + finished + ", exit code: " + exitCode);

            boolean dockerAvailable = finished && exitCode == 0;
            if (!dockerAvailable) {
                // Capture error output
                try (BufferedReader reader = new BufferedReader(new InputStreamReader(dockerCheck.getInputStream()))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        logger.warning("Docker check output: " + line);
                    }
                }
                logger.warning("Docker daemon is not available or not running");
                return false;
            }
            logger.info("Docker daemon is available");

            // Check if GATK image exists or can be pulled
            logger.info("Checking for GATK Docker image: " + GATK_DOCKER_IMAGE);
            Process imageCheck = new ProcessBuilder(dockerCommand, "images", "-q", GATK_DOCKER_IMAGE)
                .redirectErrorStream(true)
                .start();

            BufferedReader reader = new BufferedReader(new InputStreamReader(imageCheck.getInputStream()));
            String imageId = reader.readLine();
            reader.close();

            if (imageId == null || imageId.trim().isEmpty()) {
                logger.info("GATK image not found locally, will be pulled on first use");
                // Image will be pulled automatically when we try to run it
            } else {
                logger.info("GATK image found: " + imageId);
            }

            logger.info("=== GATK Availability Check Complete - SUCCESS ===");
            return true;
        } catch (Exception e) {
            logger.log(Level.SEVERE, "Error checking GATK availability", e);
            return false;
        }
    }

    /**
     * Deduplicate reads in a BAM file using GATK MarkDuplicates.
     *
     * @param inputBam Path to the input BAM file
     * @param outputDir Directory where deduplicated BAM will be saved
     * @param threads Number of threads to use (not used by GATK MarkDuplicates)
     * @return DeduplicationResult containing the path to deduplicated BAM or error message
     */
    public DeduplicationResult deduplicateBam(Path inputBam, Path outputDir, int threads) {
        logger.info("==========================================");
        logger.info("GATK DEDUPLICATION: Starting BAM deduplication");
        logger.info("  Input: " + inputBam);
        logger.info("  Output dir: " + outputDir);
        logger.info("  Docker image: " + GATK_DOCKER_IMAGE);
        logger.info("==========================================");

        // Check for null inputs
        if (inputBam == null) {
            String error = "Input BAM path is null";
            logger.warning(error);
            return new DeduplicationResult(false, null, error, 0);
        }

        if (outputDir == null) {
            String error = "Output directory path is null";
            logger.warning(error);
            return new DeduplicationResult(false, null, error, 0);
        }

        // Check if Docker and GATK are available
        if (!isGatkAvailable()) {
            String error = "Docker or GATK image is not available. Cannot perform deduplication.";
            logger.warning(error);
            return new DeduplicationResult(false, null, error, 0);
        }

        // Check if input file exists
        if (!Files.exists(inputBam)) {
            String error = "Input BAM file does not exist: " + inputBam;
            logger.warning(error);
            return new DeduplicationResult(false, null, error, 0);
        }

        // Create output directory if it doesn't exist
        try {
            Files.createDirectories(outputDir);
        } catch (IOException e) {
            String error = "Failed to create output directory: " + e.getMessage();
            logger.log(Level.WARNING, error, e);
            return new DeduplicationResult(false, null, error, 0);
        }

        // Generate output file names
        String baseName = inputBam.getFileName().toString();
        if (baseName.endsWith(".bam")) {
            baseName = baseName.substring(0, baseName.length() - 4);
        }

        Path tempDir = outputDir.resolve("gatk_dedup_temp");
        Path sortedBam = tempDir.resolve(baseName + ".sorted.bam");
        Path deduplicated = outputDir.resolve(baseName + ".dedup.bam");
        Path metricsFile = outputDir.resolve(baseName + ".dedup.metrics.txt");

        try {
            Files.createDirectories(tempDir);

            // Step 1: Sort BAM by coordinate using samtools
            logger.info("\n[Step 1/3] Sorting BAM by coordinate...");
            if (!runCommand(new String[]{
                "samtools", "sort",
                "-@", String.valueOf(threads),
                "-o", sortedBam.toString(),
                inputBam.toString()
            })) {
                return new DeduplicationResult(false, null, "Failed to sort BAM file", 0);
            }
            logger.info("  ✓ Sorting complete: " + sortedBam.getFileName());

            // Step 2: Remove duplicates using GATK MarkDuplicates via Docker
            logger.info("\n[Step 2/3] Removing duplicates with GATK MarkDuplicates...");
            logger.info("  This may take a while if the Docker image needs to be downloaded...");

            // Get absolute paths for Docker volume mounting
            Path absoluteOutputDir = outputDir.toAbsolutePath();
            Path absoluteTempDir = tempDir.toAbsolutePath();

            // Build Docker command
            List<String> dockerCmd = new ArrayList<>();
            dockerCmd.add(dockerCommand);
            dockerCmd.add("run");
            dockerCmd.add("--rm");  // Remove container after execution
            dockerCmd.add("-v");
            dockerCmd.add(absoluteOutputDir + ":/output");  // Mount output directory
            dockerCmd.add("-v");
            dockerCmd.add(absoluteTempDir + ":/temp");  // Mount temp directory
            dockerCmd.add(GATK_DOCKER_IMAGE);
            dockerCmd.add("gatk");
            dockerCmd.add("MarkDuplicates");
            dockerCmd.add("-I");
            dockerCmd.add("/temp/" + sortedBam.getFileName());
            dockerCmd.add("-O");
            dockerCmd.add("/output/" + deduplicated.getFileName());
            dockerCmd.add("-M");
            dockerCmd.add("/output/" + metricsFile.getFileName());
            dockerCmd.add("--REMOVE_DUPLICATES");
            dockerCmd.add("true");
            dockerCmd.add("--ASSUME_SORT_ORDER");
            dockerCmd.add("coordinate");
            dockerCmd.add("--CREATE_INDEX");
            dockerCmd.add("true");  // This creates the .bai index automatically

            logger.info("  Running Docker command: " + String.join(" ", dockerCmd));

            if (!runCommand(dockerCmd.toArray(new String[0]))) {
                return new DeduplicationResult(false, null, "Failed to run GATK MarkDuplicates", 0);
            }
            logger.info("  ✓ Duplicate removal complete");

            // Step 3: Parse metrics to get duplicate count
            long duplicatesRemoved = parseGatkMetrics(metricsFile);

            // Keep intermediate files for debugging
            logger.info("\nIntermediate files kept in: " + tempDir);
            logger.info("  Sorted BAM: " + sortedBam.getFileName());

            // Log summary
            logger.info("\n==========================================");
            logger.info("GATK DEDUPLICATION COMPLETE");
            logger.info("  Output: " + deduplicated);
            logger.info("  Duplicates removed: " + duplicatesRemoved);
            logger.info("  Metrics saved to: " + metricsFile);
            logger.info("  Index created: " + deduplicated + ".bai");
            logger.info("==========================================\n");

            return new DeduplicationResult(true, deduplicated.toString(), null, duplicatesRemoved);

        } catch (Exception e) {
            logger.log(Level.SEVERE, "GATK deduplication failed with exception", e);

            // Keep all files for debugging
            logger.info("Keeping intermediate files for debugging in: " + tempDir);

            return new DeduplicationResult(false, null, "Deduplication failed: " + e.getMessage(), 0);
        }
    }

    /**
     * Run a command and wait for completion.
     */
    private boolean runCommand(String[] command) {
        try {
            logger.fine("Running command: " + String.join(" ", command));

            ProcessBuilder pb = new ProcessBuilder(command);
            Process process = pb.start();

            // Read stdout
            BufferedReader stdoutReader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            BufferedReader stderrReader = new BufferedReader(new InputStreamReader(process.getErrorStream()));

            // Read output in separate threads to prevent blocking
            Thread stdoutThread = new Thread(() -> {
                try {
                    String line;
                    while ((line = stdoutReader.readLine()) != null) {
                        logger.fine("  > " + line);
                    }
                } catch (IOException e) {
                    logger.fine("Error reading stdout: " + e.getMessage());
                }
            });

            Thread stderrThread = new Thread(() -> {
                try {
                    String line;
                    while ((line = stderrReader.readLine()) != null) {
                        logger.fine("  > " + line);
                    }
                } catch (IOException e) {
                    logger.fine("Error reading stderr: " + e.getMessage());
                }
            });

            stdoutThread.start();
            stderrThread.start();

            // Wait for process to complete
            boolean finished = process.waitFor(PROCESS_TIMEOUT_MINUTES, TimeUnit.MINUTES);

            if (!finished) {
                process.destroyForcibly();
                logger.warning("Command timed out after " + PROCESS_TIMEOUT_MINUTES + " minutes");
                return false;
            }

            // Wait for output threads to finish
            stdoutThread.join(1000);
            stderrThread.join(1000);

            int exitCode = process.exitValue();
            if (exitCode != 0) {
                logger.warning("Command failed with exit code: " + exitCode);
                return false;
            }

            return true;

        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to run command: " + String.join(" ", command), e);
            return false;
        }
    }

    /**
     * Parse the GATK metrics file to extract duplicate count.
     */
    private long parseGatkMetrics(Path metricsFile) {
        if (!Files.exists(metricsFile)) {
            logger.warning("Metrics file not found: " + metricsFile);
            return 0;
        }

        try {
            List<String> lines = Files.readAllLines(metricsFile);
            boolean inMetricsSection = false;
            int duplicateIndex = -1;

            for (String line : lines) {
                line = line.trim();

                // Skip comments and empty lines
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }

                // Look for METRICS CLASS header
                if (line.startsWith("METRICS CLASS")) {
                    inMetricsSection = true;
                    continue;
                }

                if (inMetricsSection) {
                    // Find the header line
                    if (line.contains("UNPAIRED_READS_EXAMINED") || line.contains("READ_PAIRS_EXAMINED")) {
                        String[] headers = line.split("\t");
                        for (int i = 0; i < headers.length; i++) {
                            if (headers[i].equals("READ_PAIR_DUPLICATES") ||
                                headers[i].equals("UNPAIRED_READ_DUPLICATES")) {
                                duplicateIndex = i;
                            }
                        }
                    } else if (duplicateIndex >= 0) {
                        // This should be the data line
                        String[] values = line.split("\t");
                        if (values.length > duplicateIndex) {
                            try {
                                return Long.parseLong(values[duplicateIndex]);
                            } catch (NumberFormatException e) {
                                logger.fine("Could not parse duplicate count: " + values[duplicateIndex]);
                            }
                        }
                        break;
                    }
                }
            }
        } catch (Exception e) {
            logger.fine("Could not parse metrics file: " + e.getMessage());
        }

        return 0;
    }

    /**
     * Process all BAM files in a directory.
     */
    public List<DeduplicationResult> deduplicateAllBams(Path directory, int threads) {
        List<DeduplicationResult> results = new ArrayList<>();

        if (!Files.exists(directory) || !Files.isDirectory(directory)) {
            logger.warning("Directory does not exist or is not a directory: " + directory);
            return results;
        }

        try {
            // Find all BAM files in the directory
            Files.walk(directory)
                .filter(Files::isRegularFile)
                .filter(p -> p.getFileName().toString().endsWith(".bam"))
                .filter(p -> !p.getFileName().toString().contains(".dedup.")) // Skip already deduplicated
                .forEach(bamFile -> {
                    logger.info("Processing BAM file: " + bamFile.getFileName());

                    // Create output directory for this sample
                    Path outputDir = bamFile.getParent();

                    // Deduplicate the BAM
                    DeduplicationResult result = deduplicateBam(bamFile, outputDir, threads);
                    results.add(result);

                    if (result.success) {
                        logger.info("Successfully deduplicated: " + bamFile.getFileName());
                    } else {
                        logger.warning("Failed to deduplicate: " + bamFile.getFileName() + " - " + result.errorMessage);
                    }
                });

        } catch (IOException e) {
            logger.log(Level.SEVERE, "Error scanning directory for BAM files", e);
        }

        return results;
    }

    /**
     * Finds the Docker command executable, checking common installation locations.
     *
     * @return the path to the docker executable
     */
    private static String findDockerCommand() {
        logger.info("=== GATK: Searching for Docker executable ===");

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
    private static boolean isCommandAvailable(String command) {
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
}