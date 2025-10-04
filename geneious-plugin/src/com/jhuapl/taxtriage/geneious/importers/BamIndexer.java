package com.jhuapl.taxtriage.geneious.importers;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Utility class for creating BAM indices using samtools.
 */
public class BamIndexer {

    private static final Logger logger = Logger.getLogger(BamIndexer.class.getName());
    // Thread-safe publication of static state using volatile
    private static volatile String samtoolsPath = null;
    private static volatile Boolean dockerAvailable = null;

    static {
        // Try to find samtools in PATH only (no hardcoded paths)
        if (isSamtoolsAvailable("samtools")) {
            samtoolsPath = "samtools";
            logger.info("Found samtools in PATH");
            System.out.println("BamIndexer: Found native samtools in PATH");
        } else {
            logger.info("Native samtools not found in PATH");
            System.out.println("BamIndexer: Native samtools not found in PATH");
        }

        // Check if Docker is available as fallback
        dockerAvailable = isDockerAvailable();
        if (dockerAvailable) {
            System.out.println("BamIndexer: Docker is available for samtools operations");
        } else {
            System.out.println("BamIndexer: Docker is not available");
        }
    }

    /**
     * Check if Docker is available.
     */
    private static boolean isDockerAvailable() {
        try {
            Process p = Runtime.getRuntime().exec(new String[]{"docker", "version"});
            p.waitFor();
            return p.exitValue() == 0;
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Check if samtools is available at the given path.
     */
    private static boolean isSamtoolsAvailable(String path) {
        try {
            Process p = Runtime.getRuntime().exec(new String[]{path, "--version"});
            p.waitFor();
            return p.exitValue() == 0;
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Creates an index for a BAM file if it doesn't already exist.
     *
     * @param bamFile The BAM file to index
     * @return true if index was created or already exists, false on error
     */
    public static boolean ensureIndexExists(File bamFile) {
        if (bamFile == null) {
            logger.warning("BAM file is null");
            return false;
        }

        if (!bamFile.exists()) {
            logger.warning("BAM file does not exist: " + bamFile);
            return false;
        }

        if (!bamFile.isFile()) {
            logger.warning("Path is not a file: " + bamFile);
            return false;
        }

        // Check if index already exists (.bai or .bam.bai)
        File baiFile = new File(bamFile.getAbsolutePath() + ".bai");
        File bamBaiFile = new File(bamFile.getAbsolutePath().replace(".bam", ".bam.bai"));

        if (baiFile.exists() || bamBaiFile.exists()) {
            logger.info("BAM index already exists for: " + bamFile.getName());
            return true;
        }

        // Try native samtools first
        if (samtoolsPath != null) {
            if (createIndexWithNativeSamtools(bamFile)) {
                return true;
            }
        }

        // Fall back to Docker samtools
        if (dockerAvailable != null && dockerAvailable) {
            logger.info("Trying Docker samtools for indexing");
            return createIndexWithDockerSamtools(bamFile);
        }

        logger.warning("Cannot create BAM index - neither samtools nor Docker available");
        return false;
    }

    /**
     * Creates a BAM index using native samtools.
     */
    private static boolean createIndexWithNativeSamtools(File bamFile) {
        try {
            logger.info("Creating BAM index using native samtools for: " + bamFile.getName());

            String[] command = {samtoolsPath, "index", bamFile.getAbsolutePath()};

            ProcessBuilder pb = new ProcessBuilder(command);
            pb.redirectErrorStream(true);
            Process process = pb.start();

            // Read output
            StringBuilder output = new StringBuilder();
            try (BufferedReader reader = new BufferedReader(
                new InputStreamReader(process.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    output.append(line).append("\n");
                }
            }

            int exitCode = process.waitFor();

            if (exitCode == 0) {
                logger.info("Successfully created BAM index for: " + bamFile.getName());
                return true;
            } else {
                logger.warning("Failed to create BAM index. Exit code: " + exitCode +
                              "\nOutput: " + output.toString());
                return false;
            }

        } catch (Exception e) {
            logger.log(Level.SEVERE, "Error creating BAM index with native samtools", e);
            return false;
        }
    }

    /**
     * Creates a BAM index using Docker samtools.
     */
    private static boolean createIndexWithDockerSamtools(File bamFile) {
        try {
            logger.info("Creating BAM index using Docker samtools for: " + bamFile.getName());
            System.out.println("Creating BAM index using Docker samtools...");

            // Get parent directory for volume mounting
            File bamDir = bamFile.getParentFile();
            String containerBamPath = "/data/" + bamFile.getName();

            // Pull the image if needed
            Process pullProcess = Runtime.getRuntime().exec(new String[]{
                "docker", "pull", "-q", "staphb/samtools:1.19"
            });
            pullProcess.waitFor(30, java.util.concurrent.TimeUnit.SECONDS);

            // Run samtools index in Docker
            String[] dockerCommand = {
                "docker", "run", "--rm",
                "-v", bamDir.getAbsolutePath() + ":/data",
                "staphb/samtools:1.19",
                "samtools", "index", containerBamPath
            };

            ProcessBuilder pb = new ProcessBuilder(dockerCommand);
            pb.redirectErrorStream(true);
            Process process = pb.start();

            // Read output
            StringBuilder output = new StringBuilder();
            try (BufferedReader reader = new BufferedReader(
                new InputStreamReader(process.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    output.append(line).append("\n");
                    System.out.println("  Docker: " + line);
                }
            }

            boolean finished = process.waitFor(120, java.util.concurrent.TimeUnit.SECONDS);

            if (!finished) {
                logger.warning("Docker samtools index timed out");
                process.destroyForcibly();
                return false;
            }

            int exitCode = process.exitValue();

            if (exitCode == 0) {
                logger.info("Successfully created BAM index using Docker for: " + bamFile.getName());
                System.out.println("✓ Successfully created BAM index using Docker");
                return true;
            } else {
                logger.warning("Failed to create BAM index with Docker. Exit code: " + exitCode +
                              "\nOutput: " + output.toString());
                System.out.println("✗ Docker samtools index failed with exit code: " + exitCode);
                return false;
            }

        } catch (Exception e) {
            logger.log(Level.SEVERE, "Error creating BAM index with Docker samtools", e);
            System.out.println("✗ Error creating BAM index with Docker: " + e.getMessage());
            return false;
        }
    }

    /**
     * Checks if a BAM file has an index.
     */
    public static boolean hasIndex(File bamFile) {
        File baiFile = new File(bamFile.getAbsolutePath() + ".bai");
        File bamBaiFile = new File(bamFile.getAbsolutePath().replace(".bam", ".bam.bai"));
        return baiFile.exists() || bamBaiFile.exists();
    }

    /**
     * Returns whether samtools is available.
     */
    public static boolean isSamtoolsAvailable() {
        return samtoolsPath != null;
    }
}