package com.jhuapl.taxtriage.geneious.importer;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;

import java.io.File;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Tracks the results of importing TaxTriage output files into Geneious.
 *
 * This class maintains statistics about the import process, including
 * successfully imported files, failed imports, and folder structure creation.
 * It provides comprehensive reporting on the import operation results.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class ImportResult {

    private final LocalDateTime startTime;
    private LocalDateTime endTime;
    private final List<AnnotatedPluginDocument> importedDocuments;
    private final List<ImportError> errors;
    private final Map<String, Integer> fileTypeStats;
    private final List<String> createdFolders;
    private int totalFilesProcessed;
    private int successfulImports;
    private long totalBytesProcessed;

    /**
     * Creates a new import result tracker.
     */
    public ImportResult() {
        this.startTime = LocalDateTime.now();
        this.importedDocuments = new ArrayList<>();
        this.errors = new ArrayList<>();
        this.fileTypeStats = new HashMap<>();
        this.createdFolders = new ArrayList<>();
        this.totalFilesProcessed = 0;
        this.successfulImports = 0;
        this.totalBytesProcessed = 0;
    }

    /**
     * Records a successful import of a document.
     *
     * @param document the imported document
     * @param file the source file
     * @param fileType the detected file type
     */
    public void addImportedDocument(AnnotatedPluginDocument document, File file, String fileType) {
        if (document != null) {
            importedDocuments.add(document);
            successfulImports++;
            fileTypeStats.merge(fileType, 1, Integer::sum);
            if (file != null && file.exists()) {
                totalBytesProcessed += file.length();
            }
        }
        totalFilesProcessed++;
    }

    /**
     * Records an import error.
     *
     * @param file the file that failed to import
     * @param error the error that occurred
     * @param message additional error context
     */
    public void addError(File file, Exception error, String message) {
        errors.add(new ImportError(file, error, message));
        totalFilesProcessed++;
    }

    /**
     * Records the creation of a folder in Geneious.
     *
     * @param folderPath the path of the created folder
     */
    public void addCreatedFolder(String folderPath) {
        if (folderPath != null && !createdFolders.contains(folderPath)) {
            createdFolders.add(folderPath);
        }
    }

    /**
     * Marks the import operation as completed.
     */
    public void markCompleted() {
        this.endTime = LocalDateTime.now();
    }

    /**
     * Gets all successfully imported documents.
     *
     * @return list of imported documents
     */
    public List<AnnotatedPluginDocument> getImportedDocuments() {
        return new ArrayList<>(importedDocuments);
    }

    /**
     * Gets all import errors.
     *
     * @return list of import errors
     */
    public List<ImportError> getErrors() {
        return new ArrayList<>(errors);
    }

    /**
     * Gets statistics for imported file types.
     *
     * @return map of file type to count
     */
    public Map<String, Integer> getFileTypeStats() {
        return new HashMap<>(fileTypeStats);
    }

    /**
     * Gets the list of created folders.
     *
     * @return list of folder paths
     */
    public List<String> getCreatedFolders() {
        return new ArrayList<>(createdFolders);
    }

    /**
     * Gets the import start time.
     *
     * @return start timestamp
     */
    public LocalDateTime getStartTime() {
        return startTime;
    }

    /**
     * Gets the import end time.
     *
     * @return end timestamp, or null if not completed
     */
    public LocalDateTime getEndTime() {
        return endTime;
    }

    /**
     * Gets the total number of files processed.
     *
     * @return total files processed
     */
    public int getTotalFilesProcessed() {
        return totalFilesProcessed;
    }

    /**
     * Gets the number of successful imports.
     *
     * @return successful import count
     */
    public int getSuccessfulImports() {
        return successfulImports;
    }

    /**
     * Gets the number of failed imports.
     *
     * @return failed import count
     */
    public int getFailedImports() {
        return errors.size();
    }

    /**
     * Gets the total bytes processed during import.
     *
     * @return total bytes processed
     */
    public long getTotalBytesProcessed() {
        return totalBytesProcessed;
    }

    /**
     * Gets the import success rate as a percentage.
     *
     * @return success rate (0.0 to 100.0)
     */
    public double getSuccessRate() {
        if (totalFilesProcessed == 0) {
            return 0.0;
        }
        return (double) successfulImports / totalFilesProcessed * 100.0;
    }

    /**
     * Checks if the import operation was successful overall.
     *
     * @return true if there were no errors and at least one successful import
     */
    public boolean isSuccessful() {
        return errors.isEmpty() && successfulImports > 0;
    }

    /**
     * Checks if the import operation has completed.
     *
     * @return true if the import has finished
     */
    public boolean isCompleted() {
        return endTime != null;
    }

    /**
     * Generates a summary report of the import results.
     *
     * @return formatted summary string
     */
    public String getSummary() {
        StringBuilder summary = new StringBuilder();
        summary.append("=== TaxTriage Import Results ===\n\n");

        // Overall statistics
        summary.append("Import Summary:\n");
        summary.append("- Total Files Processed: ").append(totalFilesProcessed).append("\n");
        summary.append("- Successful Imports: ").append(successfulImports).append("\n");
        summary.append("- Failed Imports: ").append(getFailedImports()).append("\n");
        summary.append("- Success Rate: ").append(String.format("%.1f%%", getSuccessRate())).append("\n");
        summary.append("- Total Bytes Processed: ").append(formatBytes(totalBytesProcessed)).append("\n");
        summary.append("- Folders Created: ").append(createdFolders.size()).append("\n");

        // Timing
        if (startTime != null) {
            summary.append("- Start Time: ").append(startTime).append("\n");
        }
        if (endTime != null) {
            summary.append("- End Time: ").append(endTime).append("\n");
            summary.append("- Duration: ").append(formatDuration()).append("\n");
        }

        // File type breakdown
        if (!fileTypeStats.isEmpty()) {
            summary.append("\nFile Type Breakdown:\n");
            fileTypeStats.entrySet().stream()
                    .sorted(Map.Entry.<String, Integer>comparingByValue().reversed())
                    .forEach(entry -> summary.append("- ").append(entry.getKey())
                            .append(": ").append(entry.getValue()).append("\n"));
        }

        // Folder structure
        if (!createdFolders.isEmpty()) {
            summary.append("\nCreated Folders:\n");
            createdFolders.forEach(folder -> summary.append("- ").append(folder).append("\n"));
        }

        // Errors
        if (!errors.isEmpty()) {
            summary.append("\nErrors:\n");
            for (ImportError error : errors) {
                summary.append("- ").append(error.getFile() != null ? error.getFile().getName() : "Unknown")
                        .append(": ").append(error.getMessage()).append("\n");
            }
        }

        return summary.toString();
    }

    /**
     * Formats byte count in human-readable format.
     *
     * @param bytes the byte count
     * @return formatted string
     */
    private String formatBytes(long bytes) {
        if (bytes < 1024) {
            return bytes + " bytes";
        } else if (bytes < 1024 * 1024) {
            return String.format("%.1f KB", bytes / 1024.0);
        } else if (bytes < 1024 * 1024 * 1024) {
            return String.format("%.1f MB", bytes / (1024.0 * 1024.0));
        } else {
            return String.format("%.1f GB", bytes / (1024.0 * 1024.0 * 1024.0));
        }
    }

    /**
     * Formats the duration between start and end time.
     *
     * @return formatted duration string
     */
    private String formatDuration() {
        if (startTime == null || endTime == null) {
            return "Unknown";
        }

        long seconds = java.time.Duration.between(startTime, endTime).getSeconds();
        if (seconds < 60) {
            return seconds + " seconds";
        } else if (seconds < 3600) {
            return String.format("%d minutes, %d seconds", seconds / 60, seconds % 60);
        } else {
            return String.format("%d hours, %d minutes", seconds / 3600, (seconds % 3600) / 60);
        }
    }

    /**
     * Represents an import error.
     */
    public static class ImportError {
        private final File file;
        private final Exception exception;
        private final String message;
        private final LocalDateTime timestamp;

        /**
         * Creates a new import error.
         *
         * @param file the file that failed to import
         * @param exception the exception that occurred
         * @param message additional error context
         */
        public ImportError(File file, Exception exception, String message) {
            this.file = file;
            this.exception = exception;
            this.message = message;
            this.timestamp = LocalDateTime.now();
        }

        /**
         * Gets the file that failed to import.
         *
         * @return the file, or null if unknown
         */
        public File getFile() {
            return file;
        }

        /**
         * Gets the exception that occurred.
         *
         * @return the exception, or null if none
         */
        public Exception getException() {
            return exception;
        }

        /**
         * Gets the error message.
         *
         * @return the error message
         */
        public String getMessage() {
            return message != null ? message : (exception != null ? exception.getMessage() : "Unknown error");
        }

        /**
         * Gets the error timestamp.
         *
         * @return when the error occurred
         */
        public LocalDateTime getTimestamp() {
            return timestamp;
        }
    }
}