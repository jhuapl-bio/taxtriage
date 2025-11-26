package com.jhuapl.taxtriage.geneious.importers;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Fixes GenBank files to ensure LOCUS and ACCESSION fields match the VERSION field.
 * This is necessary for proper reference matching in Geneious when importing BAM files.
 */
public class GenBankFileFixer {

    private static final Logger logger = Logger.getLogger(GenBankFileFixer.class.getName());

    /**
     * Fixes a GenBank file to ensure LOCUS and ACCESSION match VERSION.
     *
     * @param gbFile The GenBank file to fix
     * @return true if the file was successfully fixed
     */
    public static boolean fixGenBankFile(File gbFile) {
        if (!gbFile.exists()) {
            logger.warning("GenBank file does not exist: " + gbFile.getAbsolutePath());
            return false;
        }

        try {
            logger.info("Fixing GenBank file: " + gbFile.getName());

            // Read the file
            List<String> lines = Files.readAllLines(gbFile.toPath(), StandardCharsets.UTF_8);
            if (lines.isEmpty()) {
                logger.warning("GenBank file is empty: " + gbFile.getName());
                return false;
            }

            // Find the VERSION field value
            String versionValue = null;
            int versionLineIndex = -1;
            Pattern versionPattern = Pattern.compile("^VERSION\\s+(\\S+)");

            for (int i = 0; i < lines.size(); i++) {
                Matcher matcher = versionPattern.matcher(lines.get(i));
                if (matcher.find()) {
                    versionValue = matcher.group(1);
                    versionLineIndex = i;
                    break;
                }
            }

            if (versionValue == null) {
                logger.warning("Could not find VERSION field in GenBank file: " + gbFile.getName());
                return false;
            }

            logger.info("Found VERSION: " + versionValue);

            // Fix LOCUS line (first line)
            if (!lines.isEmpty()) {
                String locusLine = lines.get(0);
                if (locusLine.startsWith("LOCUS")) {
                    // Replace the identifier in the LOCUS line with the version value
                    // LOCUS line format: LOCUS       identifier     length bp ...
                    String[] parts = locusLine.split("\\s+");
                    if (parts.length >= 3) {
                        // Reconstruct the LOCUS line with the VERSION value
                        StringBuilder newLocus = new StringBuilder();
                        newLocus.append("LOCUS       ").append(versionValue);

                        // Calculate spacing to maintain alignment
                        int spacesNeeded = 16 - versionValue.length(); // Standard GenBank format
                        if (spacesNeeded < 1) spacesNeeded = 1;
                        for (int i = 0; i < spacesNeeded; i++) {
                            newLocus.append(" ");
                        }

                        // Add the rest of the line (length, molecule type, etc.)
                        boolean foundLength = false;
                        for (int i = 2; i < parts.length; i++) {
                            if (!foundLength && parts[i].matches("\\d+")) {
                                foundLength = true;
                            }
                            if (foundLength) {
                                if (i > 2) newLocus.append(" ");
                                newLocus.append(parts[i]);
                            }
                        }

                        lines.set(0, newLocus.toString());
                        logger.info("Fixed LOCUS line: " + lines.get(0).substring(0, Math.min(50, lines.get(0).length())) + "...");
                    }
                }
            }

            // Fix ACCESSION line
            boolean accessionFixed = false;
            Pattern accessionPattern = Pattern.compile("^ACCESSION\\s+");
            for (int i = 0; i < lines.size(); i++) {
                if (accessionPattern.matcher(lines.get(i)).find()) {
                    lines.set(i, "ACCESSION   " + versionValue);
                    accessionFixed = true;
                    logger.info("Fixed ACCESSION line: " + lines.get(i));
                    break;
                }
            }

            if (!accessionFixed) {
                // If ACCESSION line doesn't exist, add it after DEFINITION
                for (int i = 0; i < lines.size(); i++) {
                    if (lines.get(i).startsWith("DEFINITION")) {
                        // Find the end of the DEFINITION (it can span multiple lines)
                        int j = i + 1;
                        while (j < lines.size() && lines.get(j).startsWith("            ")) {
                            j++;
                        }
                        // Insert ACCESSION line
                        lines.add(j, "ACCESSION   " + versionValue);
                        accessionFixed = true;
                        logger.info("Added ACCESSION line: ACCESSION   " + versionValue);
                        break;
                    }
                }
            }

            // Write the fixed file back
            Files.write(gbFile.toPath(), lines, StandardCharsets.UTF_8);
            logger.info("Successfully fixed GenBank file: " + gbFile.getName());
            return true;

        } catch (IOException e) {
            logger.log(Level.SEVERE, "Error fixing GenBank file: " + gbFile.getName(), e);
            return false;
        }
    }

    /**
     * Fixes all GenBank files in a directory.
     *
     * @param directory The directory containing GenBank files
     * @return The number of files successfully fixed
     */
    public static int fixGenBankFilesInDirectory(Path directory) {
        if (!Files.exists(directory)) {
            logger.warning("Directory does not exist: " + directory);
            return 0;
        }

        int fixedCount = 0;

        try {
            List<File> gbFiles = new ArrayList<>();
            Files.list(directory)
                .filter(Files::isRegularFile)
                .filter(p -> {
                    String name = p.getFileName().toString().toLowerCase();
                    return name.endsWith(".gb") || name.endsWith(".gbk") || name.endsWith(".genbank");
                })
                .forEach(p -> gbFiles.add(p.toFile()));

            logger.info("Found " + gbFiles.size() + " GenBank files to fix in " + directory);

            for (File gbFile : gbFiles) {
                if (fixGenBankFile(gbFile)) {
                    fixedCount++;
                }
            }

            logger.info("Fixed " + fixedCount + " of " + gbFiles.size() + " GenBank files");

        } catch (IOException e) {
            logger.log(Level.SEVERE, "Error processing directory: " + directory, e);
        }

        return fixedCount;
    }

    /**
     * Validates that a GenBank file has matching LOCUS, ACCESSION, and VERSION fields.
     *
     * @param gbFile The GenBank file to validate
     * @return true if all fields match the VERSION value
     */
    public static boolean validateGenBankFile(File gbFile) {
        try {
            List<String> lines = Files.readAllLines(gbFile.toPath(), StandardCharsets.UTF_8);

            String locusId = null;
            String accessionId = null;
            String versionId = null;

            for (String line : lines) {
                if (line.startsWith("LOCUS")) {
                    String[] parts = line.split("\\s+");
                    if (parts.length >= 2) {
                        locusId = parts[1];
                    }
                } else if (line.startsWith("ACCESSION")) {
                    String[] parts = line.trim().split("\\s+");
                    if (parts.length >= 2) {
                        accessionId = parts[1];
                    }
                } else if (line.startsWith("VERSION")) {
                    String[] parts = line.trim().split("\\s+");
                    if (parts.length >= 2) {
                        versionId = parts[1];
                    }
                }

                // Stop after finding all three
                if (locusId != null && accessionId != null && versionId != null) {
                    break;
                }
            }

            boolean valid = versionId != null &&
                           versionId.equals(locusId) &&
                           versionId.equals(accessionId);

            if (!valid) {
                logger.info("GenBank file validation for " + gbFile.getName() + ":");
                logger.info("  LOCUS: " + locusId);
                logger.info("  ACCESSION: " + accessionId);
                logger.info("  VERSION: " + versionId);
                logger.info("  Valid: " + valid);
            }

            return valid;

        } catch (IOException e) {
            logger.log(Level.WARNING, "Error validating GenBank file: " + gbFile.getName(), e);
            return false;
        }
    }
}