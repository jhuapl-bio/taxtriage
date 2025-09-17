package com.jhuapl.taxtriage.geneious.importer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Detects file formats for TaxTriage output files and maps them to appropriate
 * Geneious document types.
 *
 * This class analyzes file extensions and content to determine the correct
 * format for importing into Geneious. It handles compressed files and
 * performs content-based detection when extension-based detection is insufficient.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class FileFormatDetector {

    private static final Logger logger = Logger.getLogger(FileFormatDetector.class.getName());

    /** Supported sequence file extensions */
    private static final Set<String> SEQUENCE_EXTENSIONS = new HashSet<>(Arrays.asList(
            ".fa", ".fasta", ".fas", ".fna", ".ffn", ".faa", ".frn",
            ".fq", ".fastq"
    ));

    /** Supported alignment file extensions */
    private static final Set<String> ALIGNMENT_EXTENSIONS = new HashSet<>(Arrays.asList(
            ".sam", ".bam", ".cram"
    ));

    /** Supported variant file extensions */
    private static final Set<String> VARIANT_EXTENSIONS = new HashSet<>(Arrays.asList(
            ".vcf", ".bcf"
    ));

    /** Supported report file extensions */
    private static final Set<String> REPORT_EXTENSIONS = new HashSet<>(Arrays.asList(
            ".html", ".htm", ".xml"
    ));

    /** Supported text file extensions */
    private static final Set<String> TEXT_EXTENSIONS = new HashSet<>(Arrays.asList(
            ".txt", ".tsv", ".csv", ".log", ".out", ".err", ".report", ".kreport",
            ".classification", ".taxonomy", ".summary", ".stats"
    ));

    /** Supported compressed file extensions */
    private static final Set<String> COMPRESSED_EXTENSIONS = new HashSet<>(Arrays.asList(
            ".gz", ".bz2", ".zip", ".tar"
    ));

    /** Known TaxTriage output file patterns */
    private static final Set<String> TAXTRIAGE_PATTERNS = new HashSet<>(Arrays.asList(
            "kraken", "bracken", "fastp", "multiqc", "pipeline_info",
            "classification", "taxonomy", "count", "report"
    ));

    /**
     * Detected file format information.
     */
    public static class FileFormat {
        private final FileType type;
        private final String subtype;
        private final boolean compressed;
        private final String compressionFormat;
        private final boolean importable;

        public FileFormat(FileType type, String subtype, boolean compressed,
                         String compressionFormat, boolean importable) {
            this.type = type;
            this.subtype = subtype;
            this.compressed = compressed;
            this.compressionFormat = compressionFormat;
            this.importable = importable;
        }

        public FileType getType() { return type; }
        public String getSubtype() { return subtype; }
        public boolean isCompressed() { return compressed; }
        public String getCompressionFormat() { return compressionFormat; }
        public boolean isImportable() { return importable; }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder(type.toString());
            if (subtype != null) {
                sb.append("/").append(subtype);
            }
            if (compressed) {
                sb.append(" (").append(compressionFormat).append(")");
            }
            return sb.toString();
        }
    }

    /**
     * Supported file types for import.
     */
    public enum FileType {
        SEQUENCE("Sequence"),
        ALIGNMENT("Alignment"),
        VARIANT("Variant"),
        REPORT("Report"),
        TEXT("Text"),
        KRAKEN_REPORT("Kraken Report"),
        UNKNOWN("Unknown");

        private final String displayName;

        FileType(String displayName) {
            this.displayName = displayName;
        }

        public String getDisplayName() {
            return displayName;
        }

        @Override
        public String toString() {
            return displayName;
        }
    }

    /**
     * Detects the format of a file based on its extension and content.
     *
     * @param file the file to analyze
     * @return detected file format
     */
    public FileFormat detectFormat(File file) {
        if (file == null || !file.exists() || !file.isFile()) {
            return new FileFormat(FileType.UNKNOWN, null, false, null, false);
        }

        String fileName = file.getName().toLowerCase();

        // Check for compression
        boolean compressed = false;
        String compressionFormat = null;
        String baseFileName = fileName;

        for (String ext : COMPRESSED_EXTENSIONS) {
            if (fileName.endsWith(ext)) {
                compressed = true;
                compressionFormat = ext.substring(1); // Remove dot
                baseFileName = fileName.substring(0, fileName.length() - ext.length());
                break;
            }
        }

        // Detect format based on extension
        FileFormat extensionFormat = detectByExtension(baseFileName, compressed, compressionFormat);

        // If extension detection is inconclusive, try content detection
        if (extensionFormat.getType() == FileType.UNKNOWN || extensionFormat.getType() == FileType.TEXT) {
            FileFormat contentFormat = detectByContent(file, baseFileName);
            if (contentFormat.getType() != FileType.UNKNOWN) {
                return new FileFormat(contentFormat.getType(), contentFormat.getSubtype(),
                                    compressed, compressionFormat, contentFormat.isImportable());
            }
        }

        return extensionFormat;
    }

    /**
     * Detects file format based on extension.
     *
     * @param fileName the file name (without compression extension)
     * @param compressed whether the file is compressed
     * @param compressionFormat the compression format
     * @return detected format
     */
    private FileFormat detectByExtension(String fileName, boolean compressed, String compressionFormat) {
        // Check for sequence files
        for (String ext : SEQUENCE_EXTENSIONS) {
            if (fileName.endsWith(ext)) {
                String subtype = ext.substring(1); // Remove dot
                return new FileFormat(FileType.SEQUENCE, subtype, compressed, compressionFormat, true);
            }
        }

        // Check for alignment files
        for (String ext : ALIGNMENT_EXTENSIONS) {
            if (fileName.endsWith(ext)) {
                String subtype = ext.substring(1);
                return new FileFormat(FileType.ALIGNMENT, subtype, compressed, compressionFormat, true);
            }
        }

        // Check for variant files
        for (String ext : VARIANT_EXTENSIONS) {
            if (fileName.endsWith(ext)) {
                String subtype = ext.substring(1);
                return new FileFormat(FileType.VARIANT, subtype, compressed, compressionFormat, true);
            }
        }

        // Check for report files
        for (String ext : REPORT_EXTENSIONS) {
            if (fileName.endsWith(ext)) {
                String subtype = ext.substring(1);
                return new FileFormat(FileType.REPORT, subtype, compressed, compressionFormat, true);
            }
        }

        // Check for Kraken reports
        if (fileName.contains("kraken") && (fileName.endsWith(".report") || fileName.endsWith(".kreport"))) {
            return new FileFormat(FileType.KRAKEN_REPORT, "kreport", compressed, compressionFormat, true);
        }

        // Check for text files
        for (String ext : TEXT_EXTENSIONS) {
            if (fileName.endsWith(ext)) {
                String subtype = ext.substring(1);
                return new FileFormat(FileType.TEXT, subtype, compressed, compressionFormat, true);
            }
        }

        // Check for TaxTriage-specific patterns
        for (String pattern : TAXTRIAGE_PATTERNS) {
            if (fileName.contains(pattern)) {
                return new FileFormat(FileType.TEXT, "taxtriage", compressed, compressionFormat, true);
            }
        }

        return new FileFormat(FileType.UNKNOWN, null, compressed, compressionFormat, false);
    }

    /**
     * Detects file format based on content analysis.
     *
     * @param file the file to analyze
     * @param fileName the file name for context
     * @return detected format
     */
    private FileFormat detectByContent(File file, String fileName) {
        try {
            // For small files, read a few lines to determine content
            if (file.length() > 10 * 1024 * 1024) { // Skip content detection for files > 10MB
                return new FileFormat(FileType.UNKNOWN, null, false, null, false);
            }

            String firstLines = readFirstLines(file, 10);
            if (firstLines == null || firstLines.trim().isEmpty()) {
                return new FileFormat(FileType.UNKNOWN, null, false, null, false);
            }

            // Check for FASTA format
            if (firstLines.startsWith(">")) {
                return new FileFormat(FileType.SEQUENCE, "fasta", false, null, true);
            }

            // Check for FASTQ format
            if (firstLines.startsWith("@") && firstLines.contains("\n+\n")) {
                return new FileFormat(FileType.SEQUENCE, "fastq", false, null, true);
            }

            // Check for SAM format
            if (firstLines.startsWith("@HD") || firstLines.startsWith("@SQ") || firstLines.startsWith("@RG")) {
                return new FileFormat(FileType.ALIGNMENT, "sam", false, null, true);
            }

            // Check for VCF format
            if (firstLines.startsWith("##fileformat=VCF")) {
                return new FileFormat(FileType.VARIANT, "vcf", false, null, true);
            }

            // Check for HTML content
            if (firstLines.toLowerCase().contains("<html") || firstLines.toLowerCase().contains("<!doctype html")) {
                return new FileFormat(FileType.REPORT, "html", false, null, true);
            }

            // Check for Kraken report format
            if (isKrakenReport(firstLines, fileName)) {
                return new FileFormat(FileType.KRAKEN_REPORT, "kreport", false, null, true);
            }

            // Check for TSV/CSV format
            if (firstLines.contains("\t") || firstLines.contains(",")) {
                String subtype = firstLines.contains("\t") ? "tsv" : "csv";
                return new FileFormat(FileType.TEXT, subtype, false, null, true);
            }

            // Default to text if it appears to be readable
            if (isPrintableText(firstLines)) {
                return new FileFormat(FileType.TEXT, "txt", false, null, true);
            }

        } catch (IOException e) {
            logger.log(Level.WARNING, "Failed to analyze file content: " + file.getName(), e);
        }

        return new FileFormat(FileType.UNKNOWN, null, false, null, false);
    }

    /**
     * Reads the first few lines of a file.
     *
     * @param file the file to read
     * @param maxLines maximum number of lines to read
     * @return the first lines as a string
     * @throws IOException if reading fails
     */
    private String readFirstLines(File file, int maxLines) throws IOException {
        StringBuilder content = new StringBuilder();
        int lineCount = 0;

        try (BufferedReader reader = Files.newBufferedReader(file.toPath())) {
            String line;
            while ((line = reader.readLine()) != null && lineCount < maxLines) {
                content.append(line).append("\n");
                lineCount++;
            }
        }

        return content.toString();
    }

    /**
     * Checks if content appears to be a Kraken report.
     *
     * @param content the file content
     * @param fileName the file name for context
     * @return true if it appears to be a Kraken report
     */
    private boolean isKrakenReport(String content, String fileName) {
        // Kraken reports have specific column structure
        String[] lines = content.split("\n");
        if (lines.length < 2) {
            return false;
        }

        // Check filename contains kraken
        if (fileName.contains("kraken")) {
            return true;
        }

        // Check for Kraken report structure (percentage, reads, direct reads, rank, taxid, name)
        for (String line : lines) {
            if (line.trim().isEmpty()) continue;

            String[] parts = line.split("\t");
            if (parts.length >= 6) {
                try {
                    // First column should be percentage
                    Float.parseFloat(parts[0]);
                    // Second and third should be integers (read counts)
                    Integer.parseInt(parts[1]);
                    Integer.parseInt(parts[2]);
                    // Fifth should be taxid (integer)
                    Integer.parseInt(parts[4]);
                    return true;
                } catch (NumberFormatException e) {
                    // Continue checking other lines
                }
            }
        }

        return false;
    }

    /**
     * Checks if content appears to be printable text.
     *
     * @param content the content to check
     * @return true if it appears to be text
     */
    private boolean isPrintableText(String content) {
        if (content == null || content.isEmpty()) {
            return false;
        }

        // Check if most characters are printable
        int printableCount = 0;
        int totalCount = 0;

        for (char c : content.toCharArray()) {
            totalCount++;
            if (Character.isWhitespace(c) || (c >= 32 && c <= 126)) {
                printableCount++;
            }
        }

        // Consider it text if at least 95% of characters are printable
        return totalCount > 0 && (double) printableCount / totalCount >= 0.95;
    }

    /**
     * Checks if a file should be imported based on its format.
     *
     * @param file the file to check
     * @return true if the file should be imported
     */
    public boolean shouldImport(File file) {
        FileFormat format = detectFormat(file);
        return format.isImportable();
    }

    /**
     * Gets a human-readable description of the file format.
     *
     * @param file the file to describe
     * @return format description
     */
    public String getFormatDescription(File file) {
        FileFormat format = detectFormat(file);
        return format.toString();
    }
}