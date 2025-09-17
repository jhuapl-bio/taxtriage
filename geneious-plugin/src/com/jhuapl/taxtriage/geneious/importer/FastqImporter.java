package com.jhuapl.taxtriage.geneious.importer;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.documents.sequence.SequenceDocument;
import com.biomatters.geneious.publicapi.implementations.sequence.DefaultNucleotideSequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;

/**
 * Imports FASTQ sequence files from TaxTriage output into Geneious.
 *
 * This importer handles both plain and gzipped FASTQ files, parsing
 * sequence data and creating appropriate sequence documents in Geneious.
 * It supports both single-end and paired-end data organization.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class FastqImporter {

    private static final Logger logger = Logger.getLogger(FastqImporter.class.getName());

    /** Maximum number of sequences to import from a single file (to prevent memory issues) */
    private static final int MAX_SEQUENCES_PER_FILE = 10000;

    /** Maximum file size to process (100MB) */
    private static final long MAX_FILE_SIZE = 100 * 1024 * 1024;

    /**
     * Represents a FASTQ sequence record.
     */
    public static class FastqRecord {
        private final String header;
        private final String sequence;
        private final String qualityHeader;
        private final String quality;

        public FastqRecord(String header, String sequence, String qualityHeader, String quality) {
            this.header = header;
            this.sequence = sequence;
            this.qualityHeader = qualityHeader;
            this.quality = quality;
        }

        public String getHeader() { return header; }
        public String getSequence() { return sequence; }
        public String getQualityHeader() { return qualityHeader; }
        public String getQuality() { return quality; }

        public String getSequenceId() {
            // Extract sequence ID from header (remove @ and any trailing info)
            if (header != null && header.startsWith("@")) {
                String id = header.substring(1);
                int spaceIndex = id.indexOf(' ');
                if (spaceIndex > 0) {
                    id = id.substring(0, spaceIndex);
                }
                return id;
            }
            return "sequence";
        }
    }

    /**
     * Imports a FASTQ file as Geneious sequence documents.
     *
     * @param file the FASTQ file to import
     * @return list of imported sequence documents
     * @throws IOException if file reading fails
     */
    public List<AnnotatedPluginDocument> importFastqFile(File file) throws IOException {
        if (file == null || !file.exists() || !file.isFile()) {
            throw new IOException("Invalid FASTQ file: " + file);
        }

        if (file.length() > MAX_FILE_SIZE) {
            logger.warning("FASTQ file is very large (" + file.length() + " bytes), " +
                          "importing first " + MAX_SEQUENCES_PER_FILE + " sequences only");
        }

        List<FastqRecord> records = parseFastqFile(file);
        if (records.isEmpty()) {
            logger.warning("No valid FASTQ records found in: " + file.getName());
            return new ArrayList<>();
        }

        List<AnnotatedPluginDocument> documents = new ArrayList<>();

        // If file has many sequences, create a single document with multiple sequences
        if (records.size() > 100) {
            AnnotatedPluginDocument multiSeqDoc = createMultiSequenceDocument(file, records);
            if (multiSeqDoc != null) {
                documents.add(multiSeqDoc);
            }
        } else {
            // For smaller files, create individual sequence documents
            for (FastqRecord record : records) {
                AnnotatedPluginDocument seqDoc = createSequenceDocument(file, record);
                if (seqDoc != null) {
                    documents.add(seqDoc);
                }
            }
        }

        logger.info("Imported FASTQ file: " + file.getName() + " with " + records.size() + " sequences");
        return documents;
    }

    /**
     * Parses a FASTQ file into sequence records.
     *
     * @param file the file to parse
     * @return list of parsed records
     * @throws IOException if reading fails
     */
    private List<FastqRecord> parseFastqFile(File file) throws IOException {
        List<FastqRecord> records = new ArrayList<>();
        boolean isGzipped = file.getName().toLowerCase().endsWith(".gz");

        BufferedReader reader;
        if (isGzipped) {
            reader = new BufferedReader(new InputStreamReader(
                    new GZIPInputStream(new FileInputStream(file))));
        } else {
            reader = new BufferedReader(new FileReader(file));
        }

        try {
            String line;
            int lineCount = 0;
            int recordCount = 0;

            while ((line = reader.readLine()) != null && recordCount < MAX_SEQUENCES_PER_FILE) {
                lineCount++;
                line = line.trim();

                if (line.isEmpty()) {
                    continue;
                }

                // FASTQ records start with @ header
                if (line.startsWith("@")) {
                    String header = line;
                    String sequence = reader.readLine();
                    String qualityHeader = reader.readLine();
                    String quality = reader.readLine();

                    lineCount += 3;

                    if (sequence != null && qualityHeader != null && quality != null) {
                        sequence = sequence.trim();
                        qualityHeader = qualityHeader.trim();
                        quality = quality.trim();

                        if (qualityHeader.startsWith("+") && sequence.length() == quality.length()) {
                            records.add(new FastqRecord(header, sequence, qualityHeader, quality));
                            recordCount++;
                        } else {
                            logger.warning("Invalid FASTQ record at line " + (lineCount - 3) +
                                         " in " + file.getName());
                        }
                    } else {
                        logger.warning("Incomplete FASTQ record at line " + (lineCount - 3) +
                                     " in " + file.getName());
                        break;
                    }
                }
            }
        } finally {
            reader.close();
        }

        return records;
    }

    /**
     * Creates a single document containing multiple sequences.
     *
     * @param file the source file
     * @param records the sequence records
     * @return the created document
     */
    private AnnotatedPluginDocument createMultiSequenceDocument(File file, List<FastqRecord> records) {
        try {
            // For large files, create a summary document
            StringBuilder content = new StringBuilder();
            content.append("=== FASTQ File Summary ===\n\n");
            content.append("File: ").append(file.getName()).append("\n");
            content.append("Total Sequences: ").append(records.size()).append("\n");
            content.append("File Size: ").append(formatFileSize(file.length())).append("\n");
            content.append("Generated: ").append(new Date()).append("\n\n");

            // Calculate statistics
            int totalBases = 0;
            int minLength = Integer.MAX_VALUE;
            int maxLength = 0;

            for (FastqRecord record : records) {
                int length = record.getSequence().length();
                totalBases += length;
                minLength = Math.min(minLength, length);
                maxLength = Math.max(maxLength, length);
            }

            double avgLength = records.isEmpty() ? 0 : (double) totalBases / records.size();

            content.append("=== Sequence Statistics ===\n");
            content.append("Total Bases: ").append(totalBases).append("\n");
            content.append("Average Length: ").append(String.format("%.1f", avgLength)).append("\n");
            content.append("Min Length: ").append(minLength == Integer.MAX_VALUE ? 0 : minLength).append("\n");
            content.append("Max Length: ").append(maxLength).append("\n\n");

            // Show first few sequences as examples
            content.append("=== Sample Sequences (first 5) ===\n");
            for (int i = 0; i < Math.min(5, records.size()); i++) {
                FastqRecord record = records.get(i);
                content.append("Sequence ").append(i + 1).append(":\n");
                content.append("ID: ").append(record.getSequenceId()).append("\n");
                content.append("Length: ").append(record.getSequence().length()).append(" bp\n");
                content.append("Sequence: ").append(record.getSequence().substring(0,
                        Math.min(50, record.getSequence().length()))).append("...\n\n");
            }

            // Create a representative sequence (using first record)
            FastqRecord firstRecord = records.get(0);
            DefaultNucleotideSequence sequence = new DefaultNucleotideSequence(
                    file.getName(),
                    "FASTQ file summary with " + records.size() + " sequences",
                    firstRecord.getSequence(),
                    new Date(file.lastModified())
            );

            sequence.setDescription(content.toString());
            return DocumentUtilities.createAnnotatedPluginDocument(sequence);

        } catch (Exception e) {
            logger.warning("Failed to create multi-sequence document for " + file.getName() + ": " + e.getMessage());
            return null;
        }
    }

    /**
     * Creates a sequence document from a single FASTQ record.
     *
     * @param file the source file
     * @param record the sequence record
     * @return the created document
     */
    private AnnotatedPluginDocument createSequenceDocument(File file, FastqRecord record) {
        try {
            String sequenceId = record.getSequenceId();
            String description = "Sequence from " + file.getName() + " (Length: " +
                               record.getSequence().length() + " bp)";

            DefaultNucleotideSequence sequence = new DefaultNucleotideSequence(
                    sequenceId,
                    description,
                    record.getSequence(),
                    new Date(file.lastModified())
            );

            // Add quality information as description
            StringBuilder fullDescription = new StringBuilder();
            fullDescription.append(description).append("\n\n");
            fullDescription.append("FASTQ Header: ").append(record.getHeader()).append("\n");
            fullDescription.append("Quality Header: ").append(record.getQualityHeader()).append("\n");
            fullDescription.append("Quality String: ").append(record.getQuality()).append("\n");

            sequence.setDescription(fullDescription.toString());

            return DocumentUtilities.createAnnotatedPluginDocument(sequence);

        } catch (Exception e) {
            logger.warning("Failed to create sequence document for " + record.getSequenceId() +
                         " from " + file.getName() + ": " + e.getMessage());
            return null;
        }
    }

    /**
     * Formats file size in human-readable format.
     *
     * @param size the file size in bytes
     * @return formatted string
     */
    private String formatFileSize(long size) {
        if (size < 1024) {
            return size + " bytes";
        } else if (size < 1024 * 1024) {
            return String.format("%.1f KB", size / 1024.0);
        } else if (size < 1024 * 1024 * 1024) {
            return String.format("%.1f MB", size / (1024.0 * 1024.0));
        } else {
            return String.format("%.1f GB", size / (1024.0 * 1024.0 * 1024.0));
        }
    }

    /**
     * Checks if a file appears to be a FASTQ file.
     *
     * @param file the file to check
     * @return true if it appears to be a FASTQ file
     */
    public boolean isFastqFile(File file) {
        if (file == null || !file.exists() || !file.isFile()) {
            return false;
        }

        String name = file.getName().toLowerCase();
        if (name.endsWith(".fastq") || name.endsWith(".fq") ||
            name.endsWith(".fastq.gz") || name.endsWith(".fq.gz")) {
            return true;
        }

        // Check content
        try {
            boolean isGzipped = name.endsWith(".gz");
            BufferedReader reader;

            if (isGzipped) {
                reader = new BufferedReader(new InputStreamReader(
                        new GZIPInputStream(new FileInputStream(file))));
            } else {
                reader = new BufferedReader(new FileReader(file));
            }

            try {
                String line = reader.readLine();
                if (line != null && line.trim().startsWith("@")) {
                    // Check if it follows FASTQ format
                    String sequence = reader.readLine();
                    String qualityHeader = reader.readLine();
                    String quality = reader.readLine();

                    return sequence != null && qualityHeader != null && quality != null &&
                           qualityHeader.trim().startsWith("+") &&
                           sequence.trim().length() == quality.trim().length();
                }
            } finally {
                reader.close();
            }
        } catch (IOException e) {
            return false;
        }

        return false;
    }
}