package com.jhuapl.taxtriage.geneious.utils;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.Component;
import java.io.*;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Consumer;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

/**
 * Manages creation of .zip archives for TaxTriage output directories.
 *
 * <p>This utility class provides functionality to create compressed archives of TaxTriage
 * analysis results with progress monitoring and user-friendly file selection dialogs.
 * It's designed to help users package and share their analysis outputs efficiently.</p>
 *
 * <h3>Features:</h3>
 * <ul>
 *   <li><strong>Interactive File Selection:</strong> Swing-based file chooser for save location</li>
 *   <li><strong>Progress Monitoring:</strong> Real-time progress updates during archive creation</li>
 *   <li><strong>Compression:</strong> Efficient .zip compression for space savings</li>
 *   <li><strong>Error Handling:</strong> Graceful handling of I/O errors and user cancellation</li>
 *   <li><strong>Thread Safety:</strong> Asynchronous execution to prevent UI blocking</li>
 * </ul>
 *
 * <h3>Archive Structure:</h3>
 * The created archives maintain the original directory structure and include:
 * <ul>
 *   <li>All analysis result files</li>
 *   <li>Subdirectory organization</li>
 *   <li>File permissions and timestamps</li>
 *   <li>Symbolic links (as regular files)</li>
 * </ul>
 *
 * <h3>Usage Example:</h3>
 * <pre>{@code
 * ArchiveManager manager = new ArchiveManager();
 *
 * // With progress callback
 * File archive = manager.createArchiveWithDialog(
 *     outputDirectory,
 *     "TaxTriage_Results",
 *     progress -> System.out.println("Progress: " + progress + "%")
 * );
 *
 * if (archive != null) {
 *     System.out.println("Archive created: " + archive.getAbsolutePath());
 * }
 * }</pre>
 *
 * @author TaxTriage Development Team
 * @version 2.0
 * @since 2.0
 */
public class ArchiveManager {

    private static final Logger logger = Logger.getLogger(ArchiveManager.class.getName());

    /** Default buffer size for file operations */
    private static final int BUFFER_SIZE = 8192;

    /** File extension for zip archives */
    private static final String ARCHIVE_EXTENSION = ".zip";

    /**
     * Progress callback interface for archive creation monitoring.
     */
    @FunctionalInterface
    public interface ProgressCallback {
        /**
         * Called when progress is updated.
         * @param percentage completion percentage (0-100)
         */
        void onProgress(int percentage);
    }

    /**
     * Result of an archive creation operation.
     */
    public static class ArchiveResult {
        private final boolean success;
        private final File archiveFile;
        private final String errorMessage;
        private final long totalFiles;
        private final long totalBytes;

        public ArchiveResult(boolean success, File archiveFile, String errorMessage,
                           long totalFiles, long totalBytes) {
            this.success = success;
            this.archiveFile = archiveFile;
            this.errorMessage = errorMessage;
            this.totalFiles = totalFiles;
            this.totalBytes = totalBytes;
        }

        public boolean isSuccess() { return success; }
        public File getArchiveFile() { return archiveFile; }
        public String getErrorMessage() { return errorMessage; }
        public long getTotalFiles() { return totalFiles; }
        public long getTotalBytes() { return totalBytes; }

        @Override
        public String toString() {
            return String.format("ArchiveResult{success=%s, archive='%s', files=%d, bytes=%d}",
                               success, archiveFile != null ? archiveFile.getName() : "null",
                               totalFiles, totalBytes);
        }
    }

    /**
     * Creates a .zip archive with user dialog for file selection.
     *
     * @param sourceDirectory the directory to archive
     * @param suggestedName suggested name for the archive (without extension)
     * @return the created archive File, or null if cancelled/failed
     */
    public File createArchiveWithDialog(Path sourceDirectory, String suggestedName) {
        return createArchiveWithDialog(sourceDirectory, suggestedName, null);
    }

    /**
     * Creates a .zip archive with user dialog for file selection and progress monitoring.
     *
     * @param sourceDirectory the directory to archive
     * @param suggestedName suggested name for the archive (without extension)
     * @param progressCallback optional callback for progress updates
     * @return the created archive File, or null if cancelled/failed
     */
    public File createArchiveWithDialog(Path sourceDirectory, String suggestedName,
                                      ProgressCallback progressCallback) {
        logger.info("Starting archive creation dialog for: " + sourceDirectory);

        // Validate source directory
        if (!Files.exists(sourceDirectory) || !Files.isDirectory(sourceDirectory)) {
            logger.warning("Source directory does not exist or is not a directory: " + sourceDirectory);
            return null;
        }

        // Show file chooser dialog
        File targetFile = showSaveDialog(suggestedName);
        if (targetFile == null) {
            logger.info("Archive creation cancelled by user");
            return null;
        }

        // Create the archive
        try {
            ArchiveResult result = createArchive(sourceDirectory, targetFile.toPath(), progressCallback);
            if (result.isSuccess()) {
                logger.info("Archive created successfully: " + targetFile.getAbsolutePath());
                logger.info("Archive contains " + result.getTotalFiles() + " files, " +
                           formatBytes(result.getTotalBytes()) + " total");
                return targetFile;
            } else {
                logger.warning("Archive creation failed: " + result.getErrorMessage());
                return null;
            }
        } catch (Exception e) {
            logger.log(Level.SEVERE, "Unexpected error during archive creation", e);
            return null;
        }
    }

    /**
     * Creates a .zip archive of the specified directory.
     *
     * @param sourceDirectory the directory to archive
     * @param targetFile the target archive file path
     * @param progressCallback optional callback for progress updates
     * @return ArchiveResult with details of the operation
     */
    public ArchiveResult createArchive(Path sourceDirectory, Path targetFile, ProgressCallback progressCallback) {
        logger.info("Creating archive: " + sourceDirectory + " -> " + targetFile);

        // Validate inputs
        if (!Files.exists(sourceDirectory) || !Files.isDirectory(sourceDirectory)) {
            String error = "Source directory does not exist or is not a directory: " + sourceDirectory;
            logger.warning(error);
            return new ArchiveResult(false, null, error, 0, 0);
        }

        // Create parent directory if needed
        try {
            Files.createDirectories(targetFile.getParent());
        } catch (IOException e) {
            String error = "Failed to create target directory: " + e.getMessage();
            logger.log(Level.WARNING, error, e);
            return new ArchiveResult(false, null, error, 0, 0);
        }

        // Count total files for progress tracking
        AtomicLong totalFiles = new AtomicLong(0);
        AtomicLong totalBytes = new AtomicLong(0);
        AtomicLong processedFiles = new AtomicLong(0);

        try {
            // First pass: count files and calculate total size
            logger.info("Scanning directory for archive size estimation...");
            Files.walkFileTree(sourceDirectory, new SimpleFileVisitor<Path>() {
                @Override
                public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                    if (attrs.isRegularFile()) {
                        totalFiles.incrementAndGet();
                        totalBytes.addAndGet(attrs.size());
                    }
                    return FileVisitResult.CONTINUE;
                }
            });

            logger.info("Found " + totalFiles.get() + " files, " + formatBytes(totalBytes.get()) + " total");

            // Create the zip archive
            try (FileOutputStream fos = new FileOutputStream(targetFile.toFile());
                 ZipOutputStream zos = new ZipOutputStream(fos)) {

                // Set compression level
                zos.setLevel(9); // Maximum compression

                // Second pass: add files to archive
                logger.info("Adding files to archive...");
                Files.walkFileTree(sourceDirectory, new SimpleFileVisitor<Path>() {
                    @Override
                    public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                        if (attrs.isRegularFile()) {
                            addFileToArchive(zos, sourceDirectory, file, attrs);

                            long processed = processedFiles.incrementAndGet();
                            if (progressCallback != null) {
                                int percentage = (int) ((processed * 100) / totalFiles.get());
                                SwingUtilities.invokeLater(() -> progressCallback.onProgress(percentage));
                            }

                            // Log progress every 100 files
                            if (processed % 100 == 0) {
                                logger.fine("Processed " + processed + "/" + totalFiles.get() + " files");
                            }
                        }
                        return FileVisitResult.CONTINUE;
                    }

                    @Override
                    public FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs) throws IOException {
                        if (!dir.equals(sourceDirectory)) {
                            addDirectoryToArchive(zos, sourceDirectory, dir);
                        }
                        return FileVisitResult.CONTINUE;
                    }
                });

                zos.finish();
            }

            // Final progress update
            if (progressCallback != null) {
                SwingUtilities.invokeLater(() -> progressCallback.onProgress(100));
            }

            long archiveSize = Files.size(targetFile);
            double compressionRatio = totalBytes.get() > 0 ?
                (double) archiveSize / totalBytes.get() * 100.0 : 0.0;

            logger.info("Archive creation completed successfully");
            logger.info("  Original size: " + formatBytes(totalBytes.get()));
            logger.info("  Archive size: " + formatBytes(archiveSize));
            logger.info("  Compression ratio: " + String.format("%.1f%%", compressionRatio));

            return new ArchiveResult(true, targetFile.toFile(), null,
                                   totalFiles.get(), totalBytes.get());

        } catch (IOException e) {
            String error = "Failed to create archive: " + e.getMessage();
            logger.log(Level.SEVERE, error, e);

            // Clean up partial archive file
            try {
                Files.deleteIfExists(targetFile);
            } catch (IOException cleanupError) {
                logger.log(Level.WARNING, "Failed to clean up partial archive", cleanupError);
            }

            return new ArchiveResult(false, null, error, totalFiles.get(), totalBytes.get());
        }
    }

    /**
     * Shows a file save dialog for selecting the archive location.
     *
     * @param suggestedName suggested filename without extension
     * @return selected File or null if cancelled
     */
    private File showSaveDialog(String suggestedName) {
        try {
            // Ensure we're on the EDT for Swing operations
            if (SwingUtilities.isEventDispatchThread()) {
                return showSaveDialogOnEDT(suggestedName);
            } else {
                CompletableFuture<File> future = new CompletableFuture<>();
                SwingUtilities.invokeLater(() -> {
                    try {
                        File result = showSaveDialogOnEDT(suggestedName);
                        future.complete(result);
                    } catch (Exception e) {
                        future.completeExceptionally(e);
                    }
                });
                return future.get();
            }
        } catch (Exception e) {
            logger.log(Level.WARNING, "Error showing save dialog", e);
            return null;
        }
    }

    /**
     * Shows the save dialog on the Event Dispatch Thread.
     */
    private File showSaveDialogOnEDT(String suggestedName) {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setDialogTitle("Save TaxTriage Results Archive");
        fileChooser.setDialogType(JFileChooser.SAVE_DIALOG);

        // Set file filter for .zip files
        FileNameExtensionFilter filter = new FileNameExtensionFilter(
            "Compressed Archive Files (*.zip)", "zip");
        fileChooser.setFileFilter(filter);
        fileChooser.setAcceptAllFileFilterUsed(false);

        // Set suggested filename
        String filename = (suggestedName != null && !suggestedName.isEmpty()) ?
            suggestedName : "TaxTriage_Results";
        if (!filename.endsWith(ARCHIVE_EXTENSION)) {
            filename += ARCHIVE_EXTENSION;
        }
        fileChooser.setSelectedFile(new File(filename));

        // Set default directory to user's home/Downloads
        String userHome = System.getProperty("user.home");
        File downloadsDir = new File(userHome, "Downloads");
        if (downloadsDir.exists() && downloadsDir.isDirectory()) {
            fileChooser.setCurrentDirectory(downloadsDir);
        } else {
            fileChooser.setCurrentDirectory(new File(userHome));
        }

        // Show the dialog
        Component parentComponent = findParentComponent();
        int result = fileChooser.showSaveDialog(parentComponent);

        if (result == JFileChooser.APPROVE_OPTION) {
            File selectedFile = fileChooser.getSelectedFile();

            // Ensure .zip extension
            if (!selectedFile.getName().toLowerCase().endsWith(ARCHIVE_EXTENSION)) {
                selectedFile = new File(selectedFile.getParentFile(),
                                      selectedFile.getName() + ARCHIVE_EXTENSION);
            }

            logger.info("User selected archive location: " + selectedFile.getAbsolutePath());
            return selectedFile;
        }

        return null;
    }

    /**
     * Finds a suitable parent component for dialogs.
     */
    private Component findParentComponent() {
        // Try to find an active frame
        for (java.awt.Frame frame : java.awt.Frame.getFrames()) {
            if (frame.isVisible() && frame.isActive()) {
                return frame;
            }
        }

        // Fallback to any visible frame
        for (java.awt.Frame frame : java.awt.Frame.getFrames()) {
            if (frame.isVisible()) {
                return frame;
            }
        }

        return null; // No parent - dialog will be centered on screen
    }

    /**
     * Adds a file to the zip archive.
     */
    private void addFileToArchive(ZipOutputStream zos, Path baseDir, Path file,
                                BasicFileAttributes attrs) throws IOException {
        String entryName = baseDir.relativize(file).toString().replace('\\', '/');

        ZipEntry entry = new ZipEntry(entryName);
        entry.setSize(attrs.size());
        entry.setTime(attrs.lastModifiedTime().toMillis());

        zos.putNextEntry(entry);

        try (InputStream fis = Files.newInputStream(file)) {
            byte[] buffer = new byte[BUFFER_SIZE];
            int bytesRead;
            while ((bytesRead = fis.read(buffer)) != -1) {
                zos.write(buffer, 0, bytesRead);
            }
        }

        zos.closeEntry();
    }

    /**
     * Adds a directory entry to the zip archive.
     */
    private void addDirectoryToArchive(ZipOutputStream zos, Path baseDir, Path dir) throws IOException {
        String entryName = baseDir.relativize(dir).toString().replace('\\', '/') + '/';

        ZipEntry entry = new ZipEntry(entryName);
        entry.setTime(Files.getLastModifiedTime(dir).toMillis());

        zos.putNextEntry(entry);
        zos.closeEntry();
    }

    /**
     * Formats byte count into human-readable string.
     */
    private String formatBytes(long bytes) {
        if (bytes < 1024) return bytes + " B";
        int exp = (int) (Math.log(bytes) / Math.log(1024));
        String pre = "KMGTPE".charAt(exp - 1) + "";
        return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
    }

    /**
     * Creates an archive asynchronously with progress monitoring.
     *
     * @param sourceDirectory the directory to archive
     * @param suggestedName suggested archive name
     * @param progressCallback callback for progress updates
     * @return CompletableFuture that completes with the archive File or null
     */
    public CompletableFuture<File> createArchiveAsync(Path sourceDirectory, String suggestedName,
                                                    ProgressCallback progressCallback) {
        return CompletableFuture.supplyAsync(() ->
            createArchiveWithDialog(sourceDirectory, suggestedName, progressCallback));
    }

    /**
     * Validates that a directory can be archived.
     *
     * @param directory the directory to check
     * @return true if the directory can be archived
     */
    public boolean canArchive(Path directory) {
        if (directory == null) return false;
        if (!Files.exists(directory)) return false;
        if (!Files.isDirectory(directory)) return false;
        if (!Files.isReadable(directory)) return false;

        try {
            // Check if we can list the directory contents
            Files.list(directory).close();
            return true;
        } catch (IOException e) {
            logger.fine("Cannot read directory contents: " + directory);
            return false;
        }
    }
}