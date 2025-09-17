package com.jhuapl.taxtriage.geneious.importer;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.implementations.TextDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentOperationException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.logging.Logger;

/**
 * Imports Kraken taxonomic classification reports into Geneious.
 *
 * This importer processes Kraken and Bracken output files, parsing the
 * taxonomic classification data and creating structured documents that
 * can be viewed and analyzed within Geneious.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class KrakenReportImporter {

    private static final Logger logger = Logger.getLogger(KrakenReportImporter.class.getName());

    /**
     * Represents a single line in a Kraken report.
     */
    public static class KrakenEntry {
        private final double percentage;
        private final int totalReads;
        private final int directReads;
        private final String rank;
        private final int taxId;
        private final String name;

        public KrakenEntry(double percentage, int totalReads, int directReads,
                          String rank, int taxId, String name) {
            this.percentage = percentage;
            this.totalReads = totalReads;
            this.directReads = directReads;
            this.rank = rank;
            this.taxId = taxId;
            this.name = name;
        }

        public double getPercentage() { return percentage; }
        public int getTotalReads() { return totalReads; }
        public int getDirectReads() { return directReads; }
        public String getRank() { return rank; }
        public int getTaxId() { return taxId; }
        public String getName() { return name; }

        @Override
        public String toString() {
            return String.format("%.2f%%\t%d\t%d\t%s\t%d\t%s",
                    percentage, totalReads, directReads, rank, taxId, name);
        }
    }

    /**
     * Imports a Kraken report file as a Geneious document.
     *
     * @param file the Kraken report file
     * @return the imported document, or null if import fails
     * @throws IOException if file reading fails
     */
    public AnnotatedPluginDocument importKrakenReport(File file) throws IOException {
        if (file == null || !file.exists() || !file.isFile()) {
            throw new IOException("Invalid Kraken report file: " + file);
        }

        List<KrakenEntry> entries = parseKrakenReport(file);
        if (entries.isEmpty()) {
            logger.warning("No valid entries found in Kraken report: " + file.getName());
            return null;
        }

        // Create document content
        String content = generateReportContent(file, entries);
        String summary = generateSummary(entries);

        // Create the document
        TextDocument document;
        try {
            document = new TextDocument(file.getName(), content, TextDocument.Format.Plain);
        } catch (DocumentOperationException e) {
            throw new IOException("Failed to create TextDocument", e);
        }

        logger.info("Imported Kraken report: " + file.getName() + " with " + entries.size() + " entries");
        return DocumentUtilities.createAnnotatedPluginDocument(document);
    }

    /**
     * Parses a Kraken report file into structured entries.
     *
     * @param file the file to parse
     * @return list of parsed entries
     * @throws IOException if reading fails
     */
    private List<KrakenEntry> parseKrakenReport(File file) throws IOException {
        List<KrakenEntry> entries = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            String line;
            int lineNumber = 0;

            while ((line = reader.readLine()) != null) {
                lineNumber++;
                line = line.trim();

                if (line.isEmpty() || line.startsWith("#")) {
                    continue; // Skip empty lines and comments
                }

                try {
                    KrakenEntry entry = parseKrakenLine(line);
                    if (entry != null) {
                        entries.add(entry);
                    }
                } catch (Exception e) {
                    logger.warning("Failed to parse line " + lineNumber + " in " + file.getName() + ": " + line);
                }
            }
        }

        return entries;
    }

    /**
     * Parses a single line from a Kraken report.
     *
     * @param line the line to parse
     * @return parsed entry, or null if invalid
     */
    private KrakenEntry parseKrakenLine(String line) {
        String[] parts = line.split("\\t");
        if (parts.length < 6) {
            return null;
        }

        try {
            double percentage = Double.parseDouble(parts[0]);
            int totalReads = Integer.parseInt(parts[1]);
            int directReads = Integer.parseInt(parts[2]);
            String rank = parts[3].trim();
            int taxId = Integer.parseInt(parts[4]);
            String name = parts[5].trim();

            return new KrakenEntry(percentage, totalReads, directReads, rank, taxId, name);
        } catch (NumberFormatException e) {
            return null;
        }
    }

    /**
     * Generates formatted content for the report document.
     *
     * @param file the source file
     * @param entries the parsed entries
     * @return formatted content string
     */
    private String generateReportContent(File file, List<KrakenEntry> entries) {
        StringBuilder content = new StringBuilder();

        // Header
        content.append("=== Kraken Taxonomic Classification Report ===\n\n");
        content.append("File: ").append(file.getName()).append("\n");
        content.append("Generated: ").append(new Date()).append("\n");
        content.append("Total Entries: ").append(entries.size()).append("\n\n");

        // Summary statistics
        content.append("=== Summary Statistics ===\n");
        int totalClassified = entries.stream().mapToInt(KrakenEntry::getTotalReads).sum();
        content.append("Total Classified Reads: ").append(totalClassified).append("\n");

        // Top hits by percentage
        content.append("\n=== Top Classifications (by percentage) ===\n");
        content.append("Percentage\tTotal Reads\tDirect Reads\tRank\tTax ID\tName\n");
        content.append("----------\t-----------\t------------\t----\t------\t----\n");

        entries.stream()
                .filter(entry -> entry.getPercentage() > 0.1) // Only show significant hits
                .sorted((a, b) -> Double.compare(b.getPercentage(), a.getPercentage()))
                .limit(20) // Top 20
                .forEach(entry -> content.append(entry.toString()).append("\n"));

        // Species-level classifications
        content.append("\n=== Species-Level Classifications ===\n");
        content.append("Percentage\tTotal Reads\tDirect Reads\tTax ID\tName\n");
        content.append("----------\t-----------\t------------\t------\t----\n");

        entries.stream()
                .filter(entry -> "S".equals(entry.getRank()) && entry.getPercentage() > 0.01)
                .sorted((a, b) -> Double.compare(b.getPercentage(), a.getPercentage()))
                .limit(50)
                .forEach(entry -> content.append(String.format("%.2f%%\t%d\t%d\t%d\t%s\n",
                        entry.getPercentage(), entry.getTotalReads(), entry.getDirectReads(),
                        entry.getTaxId(), entry.getName())));

        // Full report
        content.append("\n=== Complete Classification Report ===\n");
        content.append("Percentage\tTotal Reads\tDirect Reads\tRank\tTax ID\tName\n");
        content.append("----------\t-----------\t------------\t----\t------\t----\n");

        for (KrakenEntry entry : entries) {
            content.append(entry.toString()).append("\n");
        }

        return content.toString();
    }

    /**
     * Generates a summary description for the document.
     *
     * @param entries the parsed entries
     * @return summary string
     */
    private String generateSummary(List<KrakenEntry> entries) {
        if (entries.isEmpty()) {
            return "Empty Kraken classification report";
        }

        int totalClassified = entries.stream().mapToInt(KrakenEntry::getTotalReads).sum();
        long speciesCount = entries.stream().filter(entry -> "S".equals(entry.getRank())).count();

        // Find top classification
        KrakenEntry topHit = entries.stream()
                .filter(entry -> entry.getPercentage() > 0.1)
                .max((a, b) -> Double.compare(a.getPercentage(), b.getPercentage()))
                .orElse(null);

        StringBuilder summary = new StringBuilder();
        summary.append("Kraken taxonomic classification report with ")
                .append(entries.size()).append(" total entries, ")
                .append(totalClassified).append(" classified reads, ")
                .append(speciesCount).append(" species detected.");

        if (topHit != null) {
            summary.append(" Top hit: ").append(topHit.getName())
                    .append(" (").append(String.format("%.1f%%", topHit.getPercentage())).append(")");
        }

        return summary.toString();
    }

    /**
     * Checks if a file appears to be a Kraken report.
     *
     * @param file the file to check
     * @return true if it appears to be a Kraken report
     */
    public boolean isKrakenReport(File file) {
        if (file == null || !file.exists() || !file.isFile()) {
            return false;
        }

        String name = file.getName().toLowerCase();
        if (name.contains("kraken") || name.contains("bracken") || name.endsWith(".kreport")) {
            return true;
        }

        // Check content structure
        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            String line;
            int checkedLines = 0;
            while ((line = reader.readLine()) != null && checkedLines < 10) {
                if (line.trim().isEmpty() || line.startsWith("#")) {
                    continue;
                }

                if (parseKrakenLine(line) != null) {
                    return true;
                }
                checkedLines++;
            }
        } catch (IOException e) {
            // If we can't read the file, assume it's not a Kraken report
            return false;
        }

        return false;
    }
}