package com.jhuapl.taxtriage.geneious.importer;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.implementations.TextDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentOperationException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.logging.Logger;

/**
 * Imports text files (TSV, CSV, TXT) from TaxTriage output into Geneious.
 *
 * This importer handles various text-based output files from TaxTriage,
 * including tabular data, log files, and summary reports. It provides
 * structured viewing and basic analysis of the content.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class TextFileImporter {

    private static final Logger logger = Logger.getLogger(TextFileImporter.class.getName());

    /** Maximum text file size to process (5MB) */
    private static final long MAX_TEXT_SIZE = 5 * 1024 * 1024;

    /** Maximum number of lines to preview in description */
    private static final int MAX_PREVIEW_LINES = 10;

    /**
     * Represents the structure of a text file.
     */
    public static class TextFileStructure {
        private final String fileType;
        private final boolean hasHeader;
        private final String delimiter;
        private final int columnCount;
        private final int lineCount;
        private final List<String> headers;

        public TextFileStructure(String fileType, boolean hasHeader, String delimiter,
                               int columnCount, int lineCount, List<String> headers) {
            this.fileType = fileType;
            this.hasHeader = hasHeader;
            this.delimiter = delimiter;
            this.columnCount = columnCount;
            this.lineCount = lineCount;
            this.headers = headers;
        }

        public String getFileType() { return fileType; }
        public boolean hasHeader() { return hasHeader; }
        public String getDelimiter() { return delimiter; }
        public int getColumnCount() { return columnCount; }
        public int getLineCount() { return lineCount; }
        public List<String> getHeaders() { return headers; }
    }

    /**
     * Imports a text file as a Geneious document.
     *
     * @param file the text file to import
     * @return the imported document, or null if import fails
     * @throws IOException if file reading fails
     */
    public AnnotatedPluginDocument importTextFile(File file) throws IOException {
        if (file == null || !file.exists() || !file.isFile()) {
            throw new IOException("Invalid text file: " + file);
        }

        if (file.length() > MAX_TEXT_SIZE) {
            logger.warning("Text file is very large (" + file.length() + " bytes): " + file.getName());
        }

        String content = Files.readString(file.toPath());
        TextFileStructure structure = analyzeStructure(file, content);

        // Create document
        String documentName = generateDocumentName(file, structure);
        String description = generateDescription(file, structure, content);

        TextDocument document;
        try {
            document = new TextDocument(documentName, content, TextDocument.Format.Plain);
        } catch (DocumentOperationException e) {
            throw new IOException("Failed to create TextDocument", e);
        }

        // TextDocument doesn't have setCreationDate method

        logger.info("Imported text file: " + file.getName() + " (Type: " + structure.getFileType() + ")");
        return DocumentUtilities.createAnnotatedPluginDocument(document);
    }

    /**
     * Analyzes the structure of a text file.
     *
     * @param file the file being analyzed
     * @param content the file content
     * @return structure information
     */
    private TextFileStructure analyzeStructure(File file, String content) {
        String[] lines = content.split("\n");
        if (lines.length == 0) {
            return new TextFileStructure("Empty File", false, null, 0, 0, new ArrayList<>());
        }

        // Detect file type from extension and content
        String fileName = file.getName().toLowerCase();
        String fileType = detectFileType(fileName, content);

        // Detect delimiter
        String delimiter = detectDelimiter(lines);

        // Count columns and analyze structure
        int columnCount = 0;
        boolean hasHeader = false;
        List<String> headers = new ArrayList<>();

        if (delimiter != null && lines.length > 0) {
            String[] firstLineParts = lines[0].split(delimiter, -1);
            columnCount = firstLineParts.length;

            // Check if first line looks like a header
            hasHeader = looksLikeHeader(firstLineParts, lines);

            if (hasHeader) {
                for (String part : firstLineParts) {
                    headers.add(part.trim());
                }
            }
        }

        return new TextFileStructure(fileType, hasHeader, delimiter, columnCount, lines.length, headers);
    }

    /**
     * Detects the file type based on name and content.
     *
     * @param fileName the file name
     * @param content the file content
     * @return detected file type
     */
    private String detectFileType(String fileName, String content) {
        if (fileName.endsWith(".csv")) {
            return "CSV Table";
        } else if (fileName.endsWith(".tsv")) {
            return "TSV Table";
        } else if (fileName.endsWith(".log")) {
            return "Log File";
        } else if (fileName.endsWith(".out")) {
            return "Output File";
        } else if (fileName.endsWith(".err")) {
            return "Error Log";
        } else if (fileName.contains("report") || fileName.contains("summary")) {
            return "Analysis Report";
        } else if (fileName.contains("stats") || fileName.contains("statistics")) {
            return "Statistics File";
        } else if (fileName.contains("classification")) {
            return "Classification Results";
        } else if (fileName.contains("taxonomy")) {
            return "Taxonomy File";
        } else if (fileName.contains("count")) {
            return "Count Data";
        } else if (content.contains("\t") && content.contains("\n")) {
            return "Tabular Data";
        } else if (content.contains(",") && content.contains("\n")) {
            return "Comma-Separated Data";
        } else {
            return "Text File";
        }
    }

    /**
     * Detects the delimiter used in the file.
     *
     * @param lines the file lines
     * @return detected delimiter, or null if none
     */
    private String detectDelimiter(String[] lines) {
        if (lines.length == 0) {
            return null;
        }

        // Check first few lines for consistent delimiters
        int tabCount = 0, commaCount = 0, semicolonCount = 0, pipeCount = 0;
        int linesToCheck = Math.min(5, lines.length);

        for (int i = 0; i < linesToCheck; i++) {
            tabCount += countOccurrences(lines[i], '\t');
            commaCount += countOccurrences(lines[i], ',');
            semicolonCount += countOccurrences(lines[i], ';');
            pipeCount += countOccurrences(lines[i], '|');
        }

        // Return the most common delimiter
        if (tabCount > commaCount && tabCount > semicolonCount && tabCount > pipeCount) {
            return "\t";
        } else if (commaCount > semicolonCount && commaCount > pipeCount) {
            return ",";
        } else if (semicolonCount > pipeCount) {
            return ";";
        } else if (pipeCount > 0) {
            return "\\|";
        }

        return null;
    }

    /**
     * Counts occurrences of a character in a string.
     *
     * @param str the string to search
     * @param ch the character to count
     * @return occurrence count
     */
    private int countOccurrences(String str, char ch) {
        int count = 0;
        for (int i = 0; i < str.length(); i++) {
            if (str.charAt(i) == ch) {
                count++;
            }
        }
        return count;
    }

    /**
     * Determines if the first line looks like a header.
     *
     * @param firstLineParts the parts of the first line
     * @param lines all lines in the file
     * @return true if it appears to be a header
     */
    private boolean looksLikeHeader(String[] firstLineParts, String[] lines) {
        if (lines.length < 2) {
            return false;
        }

        // Check if first line contains mostly text (not numbers)
        int textCount = 0;
        for (String part : firstLineParts) {
            if (!isNumeric(part.trim())) {
                textCount++;
            }
        }

        // If more than half the fields are non-numeric, it's likely a header
        return textCount > firstLineParts.length / 2;
    }

    /**
     * Checks if a string is numeric.
     *
     * @param str the string to check
     * @return true if numeric
     */
    private boolean isNumeric(String str) {
        if (str == null || str.isEmpty()) {
            return false;
        }
        try {
            Double.parseDouble(str);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    /**
     * Generates an appropriate document name.
     *
     * @param file the source file
     * @param structure the file structure
     * @return document name
     */
    private String generateDocumentName(File file, TextFileStructure structure) {
        String fileName = file.getName();

        // Remove common extensions
        String[] extensions = {".txt", ".tsv", ".csv", ".log", ".out", ".err", ".report"};
        for (String ext : extensions) {
            if (fileName.toLowerCase().endsWith(ext)) {
                fileName = fileName.substring(0, fileName.length() - ext.length());
                break;
            }
        }

        return fileName + " (" + structure.getFileType() + ")";
    }

    /**
     * Generates a description for the document.
     *
     * @param file the source file
     * @param structure the file structure
     * @param content the file content
     * @return description string
     */
    private String generateDescription(File file, TextFileStructure structure, String content) {
        StringBuilder description = new StringBuilder();
        description.append(structure.getFileType()).append(" from TaxTriage workflow. ");

        if (structure.getLineCount() > 0) {
            description.append(structure.getLineCount()).append(" lines");
            if (structure.getColumnCount() > 1) {
                description.append(", ").append(structure.getColumnCount()).append(" columns");
            }
            description.append(". ");
        }

        if (structure.hasHeader()) {
            description.append("Contains header row. ");
        }

        description.append("File size: ").append(formatFileSize(file.length())).append(".");

        return description.toString();
    }

    /**
     * Generates display content for the document.
     *
     * @param file the source file
     * @param structure the file structure
     * @param content the original content
     * @return formatted display content
     */
    private String generateDisplayContent(File file, TextFileStructure structure, String content) {
        StringBuilder display = new StringBuilder();

        // Header
        display.append("=== ").append(structure.getFileType()).append(" ===\n\n");
        display.append("File: ").append(file.getName()).append("\n");
        display.append("Size: ").append(formatFileSize(file.length())).append("\n");
        display.append("Lines: ").append(structure.getLineCount()).append("\n");
        display.append("Last Modified: ").append(new Date(file.lastModified())).append("\n");

        if (structure.getColumnCount() > 1) {
            display.append("Columns: ").append(structure.getColumnCount()).append("\n");
            if (structure.getDelimiter() != null) {
                String delimiterName = structure.getDelimiter().equals("\t") ? "Tab" :
                                     structure.getDelimiter().equals(",") ? "Comma" :
                                     structure.getDelimiter();
                display.append("Delimiter: ").append(delimiterName).append("\n");
            }
        }

        if (structure.hasHeader() && !structure.getHeaders().isEmpty()) {
            display.append("Headers: ").append(String.join(", ", structure.getHeaders())).append("\n");
        }

        display.append("\n=== Content ===\n");

        // Show preview of content
        String[] lines = content.split("\n");
        int linesToShow = Math.min(MAX_PREVIEW_LINES, lines.length);

        if (structure.getColumnCount() > 1 && structure.getDelimiter() != null) {
            // Format as table
            for (int i = 0; i < linesToShow; i++) {
                String[] parts = lines[i].split(structure.getDelimiter(), -1);
                for (int j = 0; j < parts.length; j++) {
                    if (j > 0) display.append("\t");
                    display.append(parts[j].trim());
                }
                display.append("\n");
            }
        } else {
            // Show as plain text
            for (int i = 0; i < linesToShow; i++) {
                display.append(lines[i]).append("\n");
            }
        }

        if (lines.length > linesToShow) {
            display.append("\n... (").append(lines.length - linesToShow).append(" more lines)\n");
        }

        display.append("\n=== Full Content ===\n");
        display.append(content);

        return display.toString();
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
        } else {
            return String.format("%.1f MB", size / (1024.0 * 1024.0));
        }
    }

    /**
     * Checks if a file appears to be a supported text file.
     *
     * @param file the file to check
     * @return true if it appears to be a text file
     */
    public boolean isTextFile(File file) {
        if (file == null || !file.exists() || !file.isFile()) {
            return false;
        }

        String name = file.getName().toLowerCase();
        String[] textExtensions = {".txt", ".tsv", ".csv", ".log", ".out", ".err",
                                  ".report", ".summary", ".stats", ".classification", ".taxonomy"};

        for (String ext : textExtensions) {
            if (name.endsWith(ext)) {
                return true;
            }
        }

        // Check if content appears to be text
        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            String firstLine = reader.readLine();
            if (firstLine != null) {
                // Check if it's printable text
                return isPrintableText(firstLine);
            }
        } catch (IOException e) {
            return false;
        }

        return false;
    }

    /**
     * Checks if a string contains mostly printable text.
     *
     * @param str the string to check
     * @return true if it appears to be text
     */
    private boolean isPrintableText(String str) {
        if (str == null || str.isEmpty()) {
            return false;
        }

        int printableCount = 0;
        for (char c : str.toCharArray()) {
            if (Character.isWhitespace(c) || (c >= 32 && c <= 126)) {
                printableCount++;
            }
        }

        // Consider it text if at least 95% of characters are printable
        return (double) printableCount / str.length() >= 0.95;
    }
}