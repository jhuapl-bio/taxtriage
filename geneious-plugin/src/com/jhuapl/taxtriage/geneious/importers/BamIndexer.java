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
    private static String samtoolsPath = null;

    static {
        // Try to find samtools in common locations
        String[] possiblePaths = {
            "/usr/local/bin/samtools",
            "/usr/bin/samtools",
            "/opt/homebrew/bin/samtools",
            "/Users/dho/miniforge3/bin/samtools",
            "samtools" // Try PATH
        };

        for (String path : possiblePaths) {
            if (isSamtoolsAvailable(path)) {
                samtoolsPath = path;
                logger.info("Found samtools at: " + path);
                break;
            }
        }

        if (samtoolsPath == null) {
            logger.warning("Samtools not found in common locations");
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
        if (!bamFile.exists()) {
            logger.warning("BAM file does not exist: " + bamFile);
            return false;
        }

        // Check if index already exists (.bai or .bam.bai)
        File baiFile = new File(bamFile.getAbsolutePath() + ".bai");
        File bamBaiFile = new File(bamFile.getAbsolutePath().replace(".bam", ".bam.bai"));

        if (baiFile.exists() || bamBaiFile.exists()) {
            logger.info("BAM index already exists for: " + bamFile.getName());
            return true;
        }

        // Try to create index using samtools
        if (samtoolsPath != null) {
            return createIndexWithSamtools(bamFile);
        }

        logger.warning("Cannot create BAM index - samtools not available");
        return false;
    }

    /**
     * Creates a BAM index using samtools.
     */
    private static boolean createIndexWithSamtools(File bamFile) {
        try {
            logger.info("Creating BAM index for: " + bamFile.getName());

            String[] command = {samtoolsPath, "index", bamFile.getAbsolutePath()};

            ProcessBuilder pb = new ProcessBuilder(command);
            pb.redirectErrorStream(true);
            Process process = pb.start();

            // Read output
            BufferedReader reader = new BufferedReader(
                new InputStreamReader(process.getInputStream())
            );
            String line;
            StringBuilder output = new StringBuilder();
            while ((line = reader.readLine()) != null) {
                output.append(line).append("\n");
            }
            reader.close();

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
            logger.log(Level.SEVERE, "Error creating BAM index", e);
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