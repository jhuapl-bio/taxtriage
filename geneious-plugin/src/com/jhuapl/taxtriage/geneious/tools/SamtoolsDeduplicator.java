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
 * Handles deduplication of mapped reads in BAM files using samtools markdup.
 *
 * This class implements the workflow with samtools markdup including preprocessing:
 * 1. Sort by name (samtools sort -n)
 * 2. Fix mate information (samtools fixmate -m)
 * 3. Sort by coordinate (samtools sort)
 * 4. Mark and remove duplicates (samtools markdup -r -S)
 * 5. Index the deduplicated BAM (samtools index)
 *
 * The -S flag ensures that supplementary alignments of duplicate reads are also marked/removed.
 * All intermediate files are kept for debugging purposes.
 */
public class SamtoolsDeduplicator {

    private static final Logger logger = Logger.getLogger(SamtoolsDeduplicator.class.getName());
    private static final int PROCESS_TIMEOUT_MINUTES = 30;

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
     * Check if samtools is available in the system PATH.
     */
    public static boolean isSamtoolsAvailable() {
        try {
            Process process = new ProcessBuilder("samtools", "--version")
                .redirectErrorStream(true)
                .start();

            boolean finished = process.waitFor(5, TimeUnit.SECONDS);
            if (finished && process.exitValue() == 0) {
                try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                    String line = reader.readLine();
                    if (line != null && line.contains("samtools")) {
                        logger.info("Samtools available: " + line);
                        return true;
                    }
                }
            }
        } catch (Exception e) {
            logger.log(Level.FINE, "Samtools not available", e);
        }
        return false;
    }

    /**
     * Deduplicate reads in a BAM file.
     *
     * @param inputBam Path to the input BAM file
     * @param outputDir Directory where deduplicated BAM will be saved
     * @param threads Number of threads to use
     * @return DeduplicationResult containing the path to deduplicated BAM or error message
     */
    public DeduplicationResult deduplicateBam(Path inputBam, Path outputDir, int threads) {
        logger.info("==========================================");
        logger.info("DEDUPLICATION: Starting BAM deduplication");
        logger.info("  Input: " + inputBam);
        logger.info("  Output dir: " + outputDir);
        logger.info("  Threads: " + threads);
        logger.info("==========================================");

        // Check if samtools is available
        if (!isSamtoolsAvailable()) {
            String error = "Samtools is not installed or not in PATH. Cannot perform deduplication.";
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

        // Generate output file names - keep all intermediate files for debugging
        String baseName = inputBam.getFileName().toString();
        if (baseName.endsWith(".bam")) {
            baseName = baseName.substring(0, baseName.length() - 4);
        }

        Path tempDir = outputDir.resolve("dedup_intermediate");
        Path nameSorted = tempDir.resolve(baseName + ".namesorted.bam");
        Path fixmated = tempDir.resolve(baseName + ".fixmate.bam");
        Path coordSorted = tempDir.resolve(baseName + ".sorted.bam");
        Path deduplicated = outputDir.resolve(baseName + ".dedup.bam");
        Path statsFile = outputDir.resolve(baseName + ".dedup.stats");

        try {
            Files.createDirectories(tempDir);

            // Step 1: Sort by name
            logger.info("\n[Step 1/5] Sorting BAM by read name...");
            if (!runSamtoolsCommand(new String[]{
                "samtools", "sort",
                "-n",  // Sort by name
                "-@", String.valueOf(threads),
                "-o", nameSorted.toString(),
                inputBam.toString()
            })) {
                return new DeduplicationResult(false, null, "Failed to sort BAM by name", 0);
            }
            logger.info("  ✓ Name sorting complete: " + nameSorted.getFileName());

            // Step 2: Fix mate information
            logger.info("\n[Step 2/5] Fixing mate information...");
            if (!runSamtoolsCommand(new String[]{
                "samtools", "fixmate",
                "-m",  // Add mate score tag for markdup
                "-@", String.valueOf(threads),
                nameSorted.toString(),
                fixmated.toString()
            })) {
                return new DeduplicationResult(false, null, "Failed to fix mate information", 0);
            }
            logger.info("  ✓ Mate fixing complete: " + fixmated.getFileName());

            // Step 3: Sort by coordinate
            logger.info("\n[Step 3/5] Sorting BAM by coordinate...");
            if (!runSamtoolsCommand(new String[]{
                "samtools", "sort",
                "-@", String.valueOf(threads),
                "-o", coordSorted.toString(),
                fixmated.toString()
            })) {
                return new DeduplicationResult(false, null, "Failed to sort BAM by coordinate", 0);
            }
            logger.info("  ✓ Coordinate sorting complete: " + coordSorted.getFileName());

            // Step 4: Mark and remove duplicates with -S flag for supplementary alignments
            logger.info("\n[Step 4/5] Marking and removing duplicates...");
            logger.info("  Using flags: -r (remove duplicates) -S (mark supplementary alignments)");
            if (!runSamtoolsCommand(new String[]{
                "samtools", "markdup",
                "-r",  // Remove duplicates (not just mark)
                "-S",  // Mark supplementary alignments of duplicates
                "-s",  // Print stats to stderr
                "-@", String.valueOf(threads),
                coordSorted.toString(),
                deduplicated.toString()
            }, statsFile)) {
                return new DeduplicationResult(false, null, "Failed to mark/remove duplicates", 0);
            }
            logger.info("  ✓ Duplicate removal complete");

            // Step 5: Index the deduplicated BAM
            logger.info("\n[Step 5/5] Indexing deduplicated BAM...");
            if (!runSamtoolsCommand(new String[]{
                "samtools", "index",
                "-@", String.valueOf(threads),
                deduplicated.toString()
            })) {
                logger.warning("Failed to index deduplicated BAM (non-critical)");
            } else {
                logger.info("  ✓ Indexing complete");
            }

            // Keep intermediate files for debugging - do not delete
            logger.info("Intermediate files kept in: " + tempDir);

            // Parse stats to get duplicate count
            long duplicatesRemoved = parseStatsFile(statsFile);

            // Log summary
            logger.info("\n==========================================");
            logger.info("DEDUPLICATION COMPLETE");
            logger.info("  Output: " + deduplicated);
            logger.info("  Duplicates removed: " + duplicatesRemoved);
            logger.info("  Stats saved to: " + statsFile);
            logger.info("==========================================\n");

            return new DeduplicationResult(true, deduplicated.toString(), null, duplicatesRemoved);

        } catch (Exception e) {
            logger.log(Level.SEVERE, "Deduplication failed with exception", e);

            // Keep all files for debugging - don't clean up
            logger.info("Keeping intermediate files for debugging in: " + tempDir);

            return new DeduplicationResult(false, null, "Deduplication failed: " + e.getMessage(), 0);
        }
    }

    /**
     * Run a samtools command.
     */
    private boolean runSamtoolsCommand(String[] command) {
        return runSamtoolsCommand(command, null);
    }

    /**
     * Run a samtools command with optional stats file.
     */
    private boolean runSamtoolsCommand(String[] command, Path statsFile) {
        try {
            logger.fine("Running command: " + String.join(" ", command));

            ProcessBuilder pb = new ProcessBuilder(command);

            // If stats file is provided, redirect stderr to it
            if (statsFile != null) {
                pb.redirectError(statsFile.toFile());
            }

            Process process = pb.start();

            // Read stdout
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    logger.fine("  > " + line);
                }
            }

            // Wait for process to complete
            boolean finished = process.waitFor(PROCESS_TIMEOUT_MINUTES, TimeUnit.MINUTES);

            if (!finished) {
                process.destroyForcibly();
                logger.warning("Command timed out after " + PROCESS_TIMEOUT_MINUTES + " minutes");
                return false;
            }

            int exitCode = process.exitValue();
            if (exitCode != 0) {
                logger.warning("Command failed with exit code: " + exitCode);

                // Read error output if not redirected
                if (statsFile == null) {
                    try (BufferedReader errorReader = new BufferedReader(new InputStreamReader(process.getErrorStream()))) {
                        String line;
                        while ((line = errorReader.readLine()) != null) {
                            logger.warning("  ERROR: " + line);
                        }
                    }
                }
                return false;
            }

            return true;

        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to run command: " + String.join(" ", command), e);
            return false;
        }
    }

    /**
     * Parse the stats file to extract duplicate count.
     */
    private long parseStatsFile(Path statsFile) {
        if (!Files.exists(statsFile)) {
            return 0;
        }

        try {
            List<String> lines = Files.readAllLines(statsFile);
            for (String line : lines) {
                // Look for line like: "DUPLICATE TOTAL: 12345"
                if (line.contains("DUPLICATE") && line.contains("TOTAL")) {
                    String[] parts = line.split(":");
                    if (parts.length > 1) {
                        String number = parts[1].trim().split("\\s+")[0];
                        return Long.parseLong(number);
                    }
                }
            }
        } catch (Exception e) {
            logger.fine("Could not parse stats file: " + e.getMessage());
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

                        // Optionally rename original file to indicate it's the original
                        try {
                            Path originalBackup = Paths.get(bamFile.toString() + ".original");
                            Files.move(bamFile, originalBackup);
                            logger.info("Original file backed up to: " + originalBackup.getFileName());
                        } catch (IOException e) {
                            logger.warning("Could not backup original file: " + e.getMessage());
                        }
                    } else {
                        logger.warning("Failed to deduplicate: " + bamFile.getFileName() + " - " + result.errorMessage);
                    }
                });

        } catch (IOException e) {
            logger.log(Level.SEVERE, "Error scanning directory for BAM files", e);
        }

        return results;
    }
}