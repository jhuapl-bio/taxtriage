package com.jhuapl.taxtriage.geneious.importers;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.zip.GZIPInputStream;

/**
 * Utility class for extracting reference sequence information from BAM files.
 * This class reads BAM headers to identify reference sequences without requiring
 * external dependencies like htsjdk.
 */
public class BamReferenceExtractor {

    private static final Logger logger = Logger.getLogger(BamReferenceExtractor.class.getName());

    // BAM magic number
    private static final byte[] BAM_MAGIC = {'B', 'A', 'M', 1};

    /**
     * Information about a reference sequence in a BAM file.
     */
    public static class ReferenceInfo {
        public final String name;
        public final int length;
        public final String accession;  // Extracted from name if present

        public ReferenceInfo(String name, int length) {
            this.name = name;
            this.length = length;
            // Try to extract accession from name (e.g., "NC_045512.2" or "NC_045512.2 Severe acute...")
            this.accession = extractAccession(name);
        }

        private String extractAccession(String name) {
            if (name == null) {
                System.out.println("        extractAccession: input is null");
                return null;
            }

            System.out.println("        extractAccession: trying to extract from '" + name + "'");

            // Common patterns for NCBI accessions
            String[] patterns = {
                "^([A-Z]{1,2}_\\d+\\.\\d+)",  // NC_045512.2
                "^([A-Z]{1,2}_\\d+)",          // NC_045512
                "^([A-Z]{2,4}\\d+\\.\\d+)",    // MN908947.3
                "^([A-Z]{2,4}\\d+)",           // MN908947
            };

            for (int i = 0; i < patterns.length; i++) {
                String pattern = patterns[i];
                java.util.regex.Pattern p = java.util.regex.Pattern.compile(pattern);
                java.util.regex.Matcher m = p.matcher(name);
                if (m.find()) {
                    String accession = m.group(1);
                    System.out.println("        extractAccession: ✓ MATCHED pattern " + (i+1) + " -> '" + accession + "'");
                    return accession;
                }
            }

            System.out.println("        extractAccession: no direct pattern match");

            // If name contains a space, check the first part
            int spaceIdx = name.indexOf(' ');
            if (spaceIdx > 0) {
                String firstPart = name.substring(0, spaceIdx);
                System.out.println("        extractAccession: has space, retrying with '" + firstPart + "'");
                return extractAccession(firstPart);
            }

            System.out.println("        extractAccession: ✗ FAILED - no accession found");
            return null;
        }

        @Override
        public String toString() {
            return String.format("ReferenceInfo[name=%s, length=%d, accession=%s]",
                                name, length, accession);
        }
    }

    /**
     * Extracts reference sequence information from a BAM file.
     *
     * @param bamFile The BAM file to read
     * @return List of reference sequences found in the BAM header
     * @throws IOException If the file cannot be read
     */
    public static List<ReferenceInfo> extractReferences(File bamFile) throws IOException {
        System.out.println("  DEBUG: BamReferenceExtractor.extractReferences() called for: " + bamFile.getName());
        List<ReferenceInfo> references = new ArrayList<>();

        // First try samtools if available
        System.out.println("  DEBUG: Attempting to use samtools for reference extraction...");
        List<ReferenceInfo> samtoolsRefs = extractReferencesUsingSamtools(bamFile);
        if (!samtoolsRefs.isEmpty()) {
            System.out.println("  DEBUG: Successfully extracted " + samtoolsRefs.size() + " reference(s) using samtools");
            return samtoolsRefs;
        }
        System.out.println("  DEBUG: Samtools extraction returned no results, falling back to direct BAM parsing");

        // Fall back to direct BAM parsing
        try (FileInputStream fis = new FileInputStream(bamFile);
             BufferedInputStream bis = new BufferedInputStream(fis)) {

            // Check BAM magic number
            byte[] magic = new byte[4];
            if (bis.read(magic) != 4) {
                throw new IOException("Could not read BAM magic number");
            }

            if (!Arrays.equals(magic, BAM_MAGIC)) {
                throw new IOException("Not a valid BAM file (invalid magic number)");
            }

            // Read header text length
            byte[] lengthBytes = new byte[4];
            if (bis.read(lengthBytes) != 4) {
                throw new IOException("Could not read header length");
            }
            int headerLength = ByteBuffer.wrap(lengthBytes).order(ByteOrder.LITTLE_ENDIAN).getInt();

            // Read and parse SAM header text
            if (headerLength > 0) {
                byte[] headerBytes = new byte[headerLength];
                if (bis.read(headerBytes) != headerLength) {
                    throw new IOException("Could not read complete header");
                }
                String headerText = new String(headerBytes);
                parseHeaderText(headerText, references);
            }

            // Read number of reference sequences
            if (bis.read(lengthBytes) != 4) {
                throw new IOException("Could not read reference count");
            }
            int refCount = ByteBuffer.wrap(lengthBytes).order(ByteOrder.LITTLE_ENDIAN).getInt();

            System.out.println("  DEBUG: BAM file contains " + refCount + " reference sequence(s)");
            logger.info("BAM file contains " + refCount + " reference sequences");

            // Read reference sequences
            for (int i = 0; i < refCount; i++) {
                // Read name length
                if (bis.read(lengthBytes) != 4) {
                    throw new IOException("Could not read reference name length");
                }
                int nameLength = ByteBuffer.wrap(lengthBytes).order(ByteOrder.LITTLE_ENDIAN).getInt();

                // Read name (null-terminated)
                byte[] nameBytes = new byte[nameLength];
                if (bis.read(nameBytes) != nameLength) {
                    throw new IOException("Could not read reference name");
                }
                String name = new String(nameBytes, 0, nameLength - 1); // Remove null terminator

                // Read sequence length
                if (bis.read(lengthBytes) != 4) {
                    throw new IOException("Could not read reference length");
                }
                int seqLength = ByteBuffer.wrap(lengthBytes).order(ByteOrder.LITTLE_ENDIAN).getInt();

                ReferenceInfo ref = new ReferenceInfo(name, seqLength);
                references.add(ref);
                System.out.println("    DEBUG: Found reference #" + (i+1) + ": name='" + ref.name + "', length=" + ref.length + ", accession=" + (ref.accession != null ? "'" + ref.accession + "'" : "null"));
                logger.fine("Found reference: " + ref);
            }

        } catch (Exception e) {
            System.out.println("  DEBUG: Exception during BAM parsing: " + e.getClass().getName() + ": " + e.getMessage());
            e.printStackTrace(System.out);
            logger.log(Level.WARNING, "Error extracting references from BAM file: " + bamFile.getName(), e);
            throw new IOException("Failed to extract references from BAM file", e);
        }

        System.out.println("  DEBUG: Successfully extracted " + references.size() + " references from BAM using direct parsing");
        return references;
    }

    /**
     * Parse SAM header text for additional reference information.
     */
    private static void parseHeaderText(String headerText, List<ReferenceInfo> references) {
        // Parse @SQ lines for sequence information
        String[] lines = headerText.split("\n");
        for (String line : lines) {
            if (line.startsWith("@SQ")) {
                // @SQ lines contain sequence dictionary information
                // Format: @SQ SN:ref_name LN:ref_length
                logger.fine("SAM header line: " + line);
            }
        }
    }

    /**
     * Attempts to use samtools to extract reference names if available.
     * This is more reliable than parsing the binary format directly.
     * First tries native samtools, then falls back to Docker if needed.
     */
    public static List<ReferenceInfo> extractReferencesUsingSamtools(File bamFile) {
        List<ReferenceInfo> references = new ArrayList<>();

        // First try native samtools if available
        if (BamIndexer.isSamtoolsAvailable()) {
            System.out.println("  DEBUG: Found native samtools, attempting to extract references...");
            references = extractReferencesWithNativeSamtools(bamFile);
            if (!references.isEmpty()) {
                return references;
            }
        } else {
            System.out.println("  DEBUG: Native samtools not available in PATH");
        }

        // Try Docker samtools as fallback
        System.out.println("  DEBUG: Attempting to use Docker samtools...");
        references = extractReferencesWithDockerSamtools(bamFile);
        if (!references.isEmpty()) {
            System.out.println("  DEBUG: Successfully extracted references using Docker samtools");
            return references;
        }

        System.out.println("  DEBUG: Docker samtools unavailable or failed");
        return references;
    }

    /**
     * Extracts references using native samtools installation.
     */
    private static List<ReferenceInfo> extractReferencesWithNativeSamtools(File bamFile) {
        List<ReferenceInfo> references = new ArrayList<>();
        try {
            ProcessBuilder pb = new ProcessBuilder("samtools", "view", "-H", bamFile.getAbsolutePath());
            Process process = pb.start();

            try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                String line;

                while ((line = reader.readLine()) != null) {
                    if (line.startsWith("@SQ")) {
                        ReferenceInfo ref = parseSQLine(line);
                        if (ref != null) {
                            references.add(ref);
                            System.out.println("    Found reference via native samtools: " + ref.name + " (length: " + ref.length + ", accession: " + ref.accession + ")");
                        }
                    }
                }
            }

            process.waitFor();

        } catch (Exception e) {
            logger.log(Level.WARNING, "Could not use native samtools to extract references", e);
            System.out.println("  DEBUG: Native samtools failed: " + e.getMessage());
        }

        return references;
    }

    /**
     * Extracts references using samtools in a Docker container.
     * Uses the staphb/samtools image which is reliable and well-maintained.
     */
    private static List<ReferenceInfo> extractReferencesWithDockerSamtools(File bamFile) {
        List<ReferenceInfo> references = new ArrayList<>();

        try {
            // Check if Docker is available
            Process dockerCheck = Runtime.getRuntime().exec(new String[]{"docker", "version"});
            dockerCheck.waitFor();
            if (dockerCheck.exitValue() != 0) {
                System.out.println("  DEBUG: Docker is not available");
                return references;
            }

            // Get parent directory for volume mounting
            File bamDir = bamFile.getParentFile();
            String containerBamPath = "/data/" + bamFile.getName();

            // Pull the image if needed (silent operation)
            System.out.println("  DEBUG: Ensuring samtools Docker image is available...");
            Process pullProcess = Runtime.getRuntime().exec(new String[]{
                "docker", "pull", "-q", "staphb/samtools:1.19"
            });
            pullProcess.waitFor(30, java.util.concurrent.TimeUnit.SECONDS);

            // Run samtools view -H in Docker
            String[] dockerCommand = {
                "docker", "run", "--rm",
                "-v", bamDir.getAbsolutePath() + ":/data:ro",
                "staphb/samtools:1.19",
                "samtools", "view", "-H", containerBamPath
            };

            System.out.println("  DEBUG: Running Docker samtools command...");
            ProcessBuilder pb = new ProcessBuilder(dockerCommand);
            pb.redirectErrorStream(true);
            Process process = pb.start();

            try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                String line;

                while ((line = reader.readLine()) != null) {
                    if (line.startsWith("@SQ")) {
                        ReferenceInfo ref = parseSQLine(line);
                        if (ref != null) {
                            references.add(ref);
                            System.out.println("    Found reference via Docker samtools: " + ref.name + " (length: " + ref.length + ", accession: " + ref.accession + ")");
                        }
                    }
                }
            }

            boolean finished = process.waitFor(60, java.util.concurrent.TimeUnit.SECONDS);

            if (!finished) {
                System.out.println("  DEBUG: Docker samtools timed out");
                process.destroyForcibly();
            }

        } catch (Exception e) {
            logger.log(Level.WARNING, "Could not use Docker samtools to extract references", e);
            System.out.println("  DEBUG: Docker samtools failed: " + e.getMessage());
        }

        return references;
    }

    /**
     * Parses an @SQ line from SAM/BAM header.
     */
    private static ReferenceInfo parseSQLine(String line) {
        String name = null;
        int length = 0;

        String[] parts = line.split("\t");
        for (String part : parts) {
            if (part.startsWith("SN:")) {
                name = part.substring(3);
            } else if (part.startsWith("LN:")) {
                try {
                    length = Integer.parseInt(part.substring(3));
                } catch (NumberFormatException e) {
                    // Ignore
                }
            }
        }

        if (name != null) {
            return new ReferenceInfo(name, length);
        }
        return null;
    }

    /**
     * Gets all unique reference sequences from multiple BAM files.
     */
    public static Set<ReferenceInfo> extractAllReferences(List<File> bamFiles) {
        Set<ReferenceInfo> allRefs = new LinkedHashSet<>();

        for (File bamFile : bamFiles) {
            try {
                // Try samtools first, fall back to direct parsing
                List<ReferenceInfo> refs = extractReferencesUsingSamtools(bamFile);
                if (refs.isEmpty()) {
                    refs = extractReferences(bamFile);
                }
                allRefs.addAll(refs);
            } catch (Exception e) {
                logger.log(Level.WARNING, "Could not extract references from: " + bamFile.getName(), e);
            }
        }

        return allRefs;
    }
}