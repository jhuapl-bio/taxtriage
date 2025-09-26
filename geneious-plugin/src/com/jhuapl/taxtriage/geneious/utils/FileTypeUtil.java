package com.jhuapl.taxtriage.geneious.utils;

import java.io.File;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Utility class for file type detection and validation in TaxTriage workflows.
 *
 * <p>This utility provides centralized file type checking logic to eliminate
 * code duplication across the TaxTriage plugin components. It supports
 * common genomic file formats and provides both case-insensitive and
 * strict validation methods.</p>
 *
 * <h3>Supported File Types:</h3>
 * <ul>
 *   <li><strong>FASTQ:</strong> .fastq, .fq, .fastq.gz, .fq.gz</li>
 *   <li><strong>FASTA:</strong> .fasta, .fa, .fas, .fna, .fasta.gz, .fa.gz</li>
 *   <li><strong>BAM:</strong> .bam, .bam.bai</li>
 *   <li><strong>SAM:</strong> .sam</li>
 *   <li><strong>GenBank:</strong> .gb, .gbk, .genbank</li>
 *   <li><strong>Text:</strong> .txt, .tsv, .csv, .log, .report</li>
 *   <li><strong>Kraken:</strong> .k2d, .kmer_distrib</li>
 * </ul>
 *
 * @author TaxTriage Development Team
 * @version 2.0
 * @since 2.0
 */
public final class FileTypeUtil {

    // FASTQ file extensions
    private static final Set<String> FASTQ_EXTENSIONS = new HashSet<>(Arrays.asList(
        ".fastq", ".fq", ".fastq.gz", ".fq.gz"
    ));

    // FASTA file extensions
    private static final Set<String> FASTA_EXTENSIONS = new HashSet<>(Arrays.asList(
        ".fasta", ".fa", ".fas", ".fna", ".fasta.gz", ".fa.gz"
    ));

    // BAM/SAM file extensions
    private static final Set<String> ALIGNMENT_EXTENSIONS = new HashSet<>(Arrays.asList(
        ".bam", ".sam", ".bam.bai"
    ));

    // GenBank file extensions
    private static final Set<String> GENBANK_EXTENSIONS = new HashSet<>(Arrays.asList(
        ".gb", ".gbk", ".genbank"
    ));

    // Text file extensions
    private static final Set<String> TEXT_EXTENSIONS = new HashSet<>(Arrays.asList(
        ".txt", ".tsv", ".csv", ".log", ".report"
    ));

    // Kraken database file extensions
    private static final Set<String> KRAKEN_DB_EXTENSIONS = new HashSet<>(Arrays.asList(
        ".k2d", ".kmer_distrib"
    ));

    // Excluded file patterns
    private static final Set<String> EXCLUDED_PATTERNS = new HashSet<>(Arrays.asList(
        ".fastp.", ".dwnld.", "references"
    ));

    /**
     * Private constructor to prevent instantiation.
     */
    private FileTypeUtil() {
        throw new UnsupportedOperationException("Utility class cannot be instantiated");
    }

    /**
     * Checks if a file is a valid FASTQ sequence file.
     *
     * @param file the file to check
     * @return true if the file is a valid FASTQ file
     */
    public static boolean isFastqFile(File file) {
        return isValidFile(file) && hasExtension(file.getName(), FASTQ_EXTENSIONS);
    }

    /**
     * Checks if a file is a valid FASTQ sequence file.
     *
     * @param path the path to check
     * @return true if the file is a valid FASTQ file
     */
    public static boolean isFastqFile(Path path) {
        return isFastqFile(path.toFile());
    }

    /**
     * Checks if a file is a valid FASTA sequence file.
     *
     * @param file the file to check
     * @return true if the file is a valid FASTA file
     */
    public static boolean isFastaFile(File file) {
        return isValidFile(file) && hasExtension(file.getName(), FASTA_EXTENSIONS);
    }

    /**
     * Checks if a file is a valid alignment file (BAM/SAM).
     *
     * @param file the file to check
     * @return true if the file is a valid alignment file
     */
    public static boolean isAlignmentFile(File file) {
        return isValidFile(file) && hasExtension(file.getName(), ALIGNMENT_EXTENSIONS);
    }

    /**
     * Checks if a file is a BAM file specifically.
     *
     * @param file the file to check
     * @return true if the file is a BAM file
     */
    public static boolean isBamFile(File file) {
        return isValidFile(file) && file.getName().toLowerCase().endsWith(".bam");
    }

    /**
     * Checks if a file is a GenBank file.
     *
     * @param file the file to check
     * @return true if the file is a GenBank file
     */
    public static boolean isGenBankFile(File file) {
        return isValidFile(file) && hasExtension(file.getName(), GENBANK_EXTENSIONS);
    }

    /**
     * Checks if a file is a text file.
     *
     * @param file the file to check
     * @return true if the file is a text file
     */
    public static boolean isTextFile(File file) {
        return isValidFile(file) && hasExtension(file.getName(), TEXT_EXTENSIONS);
    }

    /**
     * Checks if a file is a Kraken database file.
     *
     * @param file the file to check
     * @return true if the file is a Kraken database file
     */
    public static boolean isKrakenDatabaseFile(File file) {
        return isValidFile(file) && hasExtension(file.getName(), KRAKEN_DB_EXTENSIONS);
    }

    /**
     * Checks if a file is a valid sequence file (FASTQ only for TaxTriage).
     * Excludes specific patterns that shouldn't be included in analysis.
     *
     * @param file the file to check
     * @return true if the file is a valid sequence file for TaxTriage
     */
    public static boolean isValidSequenceFile(File file) {
        if (!isValidFile(file)) {
            return false;
        }

        String name = file.getName().toLowerCase();

        // Check for excluded patterns
        for (String pattern : EXCLUDED_PATTERNS) {
            if (name.contains(pattern)) {
                return false;
            }
        }

        // Only accept FASTQ files for TaxTriage analysis
        return isFastqFile(file);
    }

    /**
     * Checks if a file is a deduplicated BAM file.
     *
     * @param file the file to check
     * @return true if the file is a deduplicated BAM file
     */
    public static boolean isDeduplicatedBamFile(File file) {
        return isBamFile(file) && file.getName().contains(".dedup.");
    }

    /**
     * Checks if a file is an original (non-deduplicated) BAM file.
     *
     * @param file the file to check
     * @return true if the file is an original BAM file
     */
    public static boolean isOriginalBamFile(File file) {
        return isBamFile(file) && !isDeduplicatedBamFile(file);
    }

    /**
     * Gets the base name of a file without extension.
     *
     * @param filename the filename
     * @return the base name without extension
     */
    public static String getBaseName(String filename) {
        if (filename == null || filename.isEmpty()) {
            return "";
        }

        // Handle double extensions like .fastq.gz
        String lowerName = filename.toLowerCase();
        if (lowerName.endsWith(".gz")) {
            // Remove .gz first
            String withoutGz = filename.substring(0, filename.length() - 3);
            // Then remove the next extension
            int lastDot = withoutGz.lastIndexOf('.');
            if (lastDot > 0) {
                return withoutGz.substring(0, lastDot);
            }
            return withoutGz;
        }

        // Standard single extension
        int lastDot = filename.lastIndexOf('.');
        if (lastDot > 0) {
            return filename.substring(0, lastDot);
        }
        return filename;
    }

    /**
     * Checks if a file exists and is a regular file.
     *
     * @param file the file to check
     * @return true if the file exists and is a regular file
     */
    private static boolean isValidFile(File file) {
        return file != null && file.exists() && file.isFile();
    }

    /**
     * Checks if a filename has any of the specified extensions (case-insensitive).
     *
     * @param filename the filename to check
     * @param extensions the set of extensions to check against
     * @return true if the filename has one of the extensions
     */
    private static boolean hasExtension(String filename, Set<String> extensions) {
        if (filename == null) {
            return false;
        }

        String lowerName = filename.toLowerCase();
        return extensions.stream().anyMatch(lowerName::endsWith);
    }
}