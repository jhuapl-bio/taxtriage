package com.jhuapl.taxtriage.geneious.utils;

import com.jhuapl.taxtriage.geneious.utils.ArchiveManager.ArchiveResult;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.atomic.AtomicInteger;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for ArchiveManager (Group 4).
 *
 * Tests the archive creation functionality including:
 * - Archive creation and compression
 * - Progress callback handling
 * - Error handling for invalid inputs
 * - Directory validation
 * - Archive result structure
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class ArchiveManagerTest {

    @TempDir
    Path tempDir;

    private ArchiveManager archiveManager;

    @BeforeEach
    void setUp() {
        archiveManager = new ArchiveManager();
    }

    /**
     * GROUP 4: Test archive creation with empty directory
     */
    @Test
    void testCreateArchiveEmptyDirectory() throws Exception {
        Path emptyDir = tempDir.resolve("empty");
        Files.createDirectories(emptyDir);

        Path targetFile = tempDir.resolve("archive.zip");

        ArchiveResult result = archiveManager.createArchive(
            emptyDir, targetFile, null);

        assertNotNull(result, "Result should not be null");
        assertTrue(result.isSuccess(), "Should succeed with empty directory");
        assertEquals(0, result.getTotalFiles(), "Should have zero files");
        assertTrue(Files.exists(targetFile), "Archive file should be created");
    }

    /**
     * GROUP 4: Test archive creation with non-existent source
     */
    @Test
    void testCreateArchiveNonExistentSource() {
        Path nonExistent = tempDir.resolve("nonexistent");
        Path targetFile = tempDir.resolve("archive.zip");

        ArchiveResult result = archiveManager.createArchive(
            nonExistent, targetFile, null);

        assertNotNull(result, "Result should not be null");
        assertFalse(result.isSuccess(), "Should fail for non-existent source");
        assertNotNull(result.getErrorMessage(), "Should have error message");
    }

    /**
     * GROUP 4: Test archive creation with file as source
     */
    @Test
    void testCreateArchiveWithFileAsSource() throws Exception {
        Path sourceFile = tempDir.resolve("test.txt");
        Files.writeString(sourceFile, "test content");

        Path targetFile = tempDir.resolve("archive.zip");

        ArchiveResult result = archiveManager.createArchive(
            sourceFile, targetFile, null);

        assertNotNull(result, "Result should not be null");
        assertFalse(result.isSuccess(), "Should fail when source is a file");
    }

    /**
     * GROUP 4: Test archive creation with single file
     */
    @Test
    void testCreateArchiveWithSingleFile() throws Exception {
        Path sourceDir = tempDir.resolve("source");
        Files.createDirectories(sourceDir);
        Files.writeString(sourceDir.resolve("test.txt"), "test content");

        Path targetFile = tempDir.resolve("archive.zip");

        ArchiveResult result = archiveManager.createArchive(
            sourceDir, targetFile, null);

        assertNotNull(result, "Result should not be null");
        assertTrue(result.isSuccess(), "Should succeed with single file");
        assertEquals(1, result.getTotalFiles(), "Should have one file");
        assertTrue(Files.exists(targetFile), "Archive should be created");
    }

    /**
     * GROUP 4: Test archive creation with multiple files
     */
    @Test
    void testCreateArchiveWithMultipleFiles() throws Exception {
        Path sourceDir = tempDir.resolve("source");
        Files.createDirectories(sourceDir);

        // Create multiple files
        for (int i = 0; i < 10; i++) {
            Files.writeString(sourceDir.resolve("file" + i + ".txt"),
                "content " + i);
        }

        Path targetFile = tempDir.resolve("archive.zip");

        ArchiveResult result = archiveManager.createArchive(
            sourceDir, targetFile, null);

        assertNotNull(result, "Result should not be null");
        assertTrue(result.isSuccess(), "Should succeed with multiple files");
        assertEquals(10, result.getTotalFiles(), "Should have 10 files");
        assertTrue(result.getTotalBytes() > 0, "Should have non-zero bytes");
    }

    /**
     * GROUP 4: Test archive creation with subdirectories
     */
    @Test
    void testCreateArchiveWithSubdirectories() throws Exception {
        Path sourceDir = tempDir.resolve("source");
        Files.createDirectories(sourceDir);

        // Create subdirectory structure
        Path subDir1 = sourceDir.resolve("subdir1");
        Path subDir2 = sourceDir.resolve("subdir2");
        Files.createDirectories(subDir1);
        Files.createDirectories(subDir2);

        Files.writeString(subDir1.resolve("file1.txt"), "content 1");
        Files.writeString(subDir2.resolve("file2.txt"), "content 2");

        Path targetFile = tempDir.resolve("archive.zip");

        ArchiveResult result = archiveManager.createArchive(
            sourceDir, targetFile, null);

        assertNotNull(result, "Result should not be null");
        assertTrue(result.isSuccess(), "Should succeed with subdirectories");
        assertEquals(2, result.getTotalFiles(), "Should have 2 files");
    }

    /**
     * GROUP 4: Test progress callback
     */
    @Test
    void testProgressCallback() throws Exception {
        Path sourceDir = tempDir.resolve("source");
        Files.createDirectories(sourceDir);

        // Create multiple files to trigger progress updates
        for (int i = 0; i < 50; i++) {
            Files.writeString(sourceDir.resolve("file" + i + ".txt"),
                "content " + i);
        }

        Path targetFile = tempDir.resolve("archive.zip");

        AtomicInteger lastProgress = new AtomicInteger(0);
        AtomicInteger callCount = new AtomicInteger(0);

        ArchiveResult result = archiveManager.createArchive(
            sourceDir, targetFile, percentage -> {
                lastProgress.set(percentage);
                callCount.incrementAndGet();
            });

        assertTrue(result.isSuccess(), "Archive creation should succeed");
        assertTrue(callCount.get() > 0, "Progress callback should be called");
        // Progress may not reach exactly 100% due to async Swing updates
        assertTrue(lastProgress.get() >= 50,
            "Final progress should be at least 50% (was " + lastProgress.get() + "%)");
    }

    /**
     * GROUP 4: Test archive result structure
     */
    @Test
    void testArchiveResultStructure() {
        File archiveFile = new File("/path/to/archive.zip");
        ArchiveResult result = new ArchiveResult(
            true, archiveFile, null, 10, 1024);

        assertTrue(result.isSuccess());
        assertEquals(archiveFile, result.getArchiveFile());
        assertNull(result.getErrorMessage());
        assertEquals(10, result.getTotalFiles());
        assertEquals(1024, result.getTotalBytes());

        String str = result.toString();
        assertNotNull(str);
        assertTrue(str.contains("success=true"));
    }

    /**
     * GROUP 4: Test archive result with error
     */
    @Test
    void testArchiveResultWithError() {
        ArchiveResult result = new ArchiveResult(
            false, null, "Test error", 0, 0);

        assertFalse(result.isSuccess());
        assertNull(result.getArchiveFile());
        assertEquals("Test error", result.getErrorMessage());
    }

    /**
     * GROUP 4: Test canArchive validation
     */
    @Test
    void testCanArchiveValidation() throws Exception {
        // Valid directory
        Path validDir = tempDir.resolve("valid");
        Files.createDirectories(validDir);
        assertTrue(archiveManager.canArchive(validDir),
            "Should be able to archive valid directory");

        // Null path
        assertFalse(archiveManager.canArchive(null),
            "Should not archive null path");

        // Non-existent
        Path nonExistent = tempDir.resolve("nonexistent");
        assertFalse(archiveManager.canArchive(nonExistent),
            "Should not archive non-existent path");

        // File instead of directory
        Path file = tempDir.resolve("file.txt");
        Files.createFile(file);
        assertFalse(archiveManager.canArchive(file),
            "Should not archive file");
    }

    /**
     * GROUP 4: Test async archive creation
     */
    @Test
    void testCreateArchiveAsync() throws Exception {
        Path sourceDir = tempDir.resolve("source");
        Files.createDirectories(sourceDir);
        Files.writeString(sourceDir.resolve("test.txt"), "content");

        // Skip this test in headless mode (CI/automated testing) as it requires GUI
        boolean isHeadless = java.awt.GraphicsEnvironment.isHeadless();
        if (isHeadless) {
            System.out.println("Skipping testCreateArchiveAsync in headless mode");
            return;
        }

        // Note: This will try to show a dialog, which won't work in headless tests
        // We're testing that it doesn't throw an exception
        assertDoesNotThrow(() -> {
            CompletableFuture<File> future = archiveManager.createArchiveAsync(
                sourceDir, "test", null);
            assertNotNull(future, "Future should not be null");
        });
    }

    /**
     * GROUP 4: Test archive creation with null progress callback
     */
    @Test
    void testCreateArchiveNullCallback() throws Exception {
        Path sourceDir = tempDir.resolve("source");
        Files.createDirectories(sourceDir);
        Files.writeString(sourceDir.resolve("test.txt"), "content");

        Path targetFile = tempDir.resolve("archive.zip");

        ArchiveResult result = archiveManager.createArchive(
            sourceDir, targetFile, null);

        assertTrue(result.isSuccess(),
            "Should succeed with null progress callback");
    }

    /**
     * GROUP 4: Test archive compression ratio
     */
    @Test
    void testArchiveCompression() throws Exception {
        Path sourceDir = tempDir.resolve("source");
        Files.createDirectories(sourceDir);

        // Create file with compressible content
        String compressibleContent = "AAAA".repeat(1000);
        Files.writeString(sourceDir.resolve("compressible.txt"),
            compressibleContent);

        Path targetFile = tempDir.resolve("archive.zip");

        ArchiveResult result = archiveManager.createArchive(
            sourceDir, targetFile, null);

        assertTrue(result.isSuccess(), "Archive creation should succeed");
        long archiveSize = Files.size(targetFile);
        assertTrue(archiveSize < result.getTotalBytes(),
            "Archive should be compressed");
    }

    /**
     * GROUP 4: Test archive with large number of files
     */
    @Test
    void testArchiveWithManyFiles() throws Exception {
        Path sourceDir = tempDir.resolve("source");
        Files.createDirectories(sourceDir);

        // Create 200 small files
        for (int i = 0; i < 200; i++) {
            Files.writeString(sourceDir.resolve("file" + i + ".txt"),
                "content " + i);
        }

        Path targetFile = tempDir.resolve("archive.zip");

        ArchiveResult result = archiveManager.createArchive(
            sourceDir, targetFile, null);

        assertTrue(result.isSuccess(), "Should handle many files");
        assertEquals(200, result.getTotalFiles(), "Should archive all files");
    }

    /**
     * GROUP 4: Test archive with special characters in filenames
     */
    @Test
    void testArchiveWithSpecialCharacters() throws Exception {
        Path sourceDir = tempDir.resolve("source");
        Files.createDirectories(sourceDir);

        // Create files with special characters (OS-safe ones)
        Files.writeString(sourceDir.resolve("file-with-dash.txt"), "content");
        Files.writeString(sourceDir.resolve("file_with_underscore.txt"), "content");
        Files.writeString(sourceDir.resolve("file.with.dots.txt"), "content");

        Path targetFile = tempDir.resolve("archive.zip");

        ArchiveResult result = archiveManager.createArchive(
            sourceDir, targetFile, null);

        assertTrue(result.isSuccess(),
            "Should handle special characters in filenames");
        assertEquals(3, result.getTotalFiles(), "Should archive all files");
    }

    /**
     * GROUP 4: Test archive overwrites existing file
     */
    @Test
    void testArchiveOverwritesExisting() throws Exception {
        Path sourceDir = tempDir.resolve("source");
        Files.createDirectories(sourceDir);
        Files.writeString(sourceDir.resolve("test.txt"), "content");

        Path targetFile = tempDir.resolve("archive.zip");

        // Create first archive
        ArchiveResult result1 = archiveManager.createArchive(
            sourceDir, targetFile, null);
        assertTrue(result1.isSuccess(), "First archive should succeed");

        long firstSize = Files.size(targetFile);

        // Create second archive (should overwrite)
        Files.writeString(sourceDir.resolve("test2.txt"), "more content");
        ArchiveResult result2 = archiveManager.createArchive(
            sourceDir, targetFile, null);
        assertTrue(result2.isSuccess(), "Second archive should succeed");

        long secondSize = Files.size(targetFile);
        assertNotEquals(firstSize, secondSize,
            "Archive should be overwritten with different size");
    }
}
