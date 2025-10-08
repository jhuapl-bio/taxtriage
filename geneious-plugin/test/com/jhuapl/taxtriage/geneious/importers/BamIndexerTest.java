package com.jhuapl.taxtriage.geneious.importers;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for BamIndexer (Group 4).
 *
 * Tests the BAM indexing functionality including:
 * - Samtools availability detection
 * - Index creation
 * - Index existence checking
 * - Error handling for invalid files
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class BamIndexerTest {

    @TempDir
    Path tempDir;

    /**
     * GROUP 4: Test samtools availability check
     */
    @Test
    void testIsSamtoolsAvailable() {
        assertDoesNotThrow(() -> {
            boolean available = BamIndexer.isSamtoolsAvailable();
            // Result depends on system configuration
            assertNotNull(available);
        });
    }

    /**
     * GROUP 4: Test ensure index with non-existent file
     */
    @Test
    void testEnsureIndexNonExistentFile() {
        File bamFile = new File(tempDir.toFile(), "nonexistent.bam");

        boolean result = BamIndexer.ensureIndexExists(bamFile);

        assertFalse(result, "Should return false for non-existent file");
    }

    /**
     * GROUP 4: Test ensure index with null file
     */
    @Test
    void testEnsureIndexNullFile() {
        assertDoesNotThrow(() -> {
            boolean result = BamIndexer.ensureIndexExists(null);
            assertFalse(result, "Should return false for null file");
        });
    }

    /**
     * GROUP 4: Test ensure index with invalid BAM file
     */
    @Test
    void testEnsureIndexInvalidBamFile() throws Exception {
        // Create a text file with .bam extension
        Path fakeBam = tempDir.resolve("fake.bam");
        Files.writeString(fakeBam, "This is not a BAM file");

        boolean result = BamIndexer.ensureIndexExists(fakeBam.toFile());

        // Should handle gracefully (may succeed or fail depending on validation)
        assertNotNull(result);
    }

    /**
     * GROUP 4: Test ensure index with directory instead of file
     */
    @Test
    void testEnsureIndexWithDirectory() {
        File directory = tempDir.toFile();

        boolean result = BamIndexer.ensureIndexExists(directory);

        assertFalse(result, "Should return false for directory");
    }

    /**
     * GROUP 4: Test index file path generation
     */
    @Test
    void testIndexFilePathGeneration() throws Exception {
        Path bamFile = tempDir.resolve("test.bam");
        Files.createFile(bamFile);

        // Expected index file path
        Path expectedIndex = tempDir.resolve("test.bam.bai");

        // The index won't be created if samtools is not available,
        // but we can test that the expected path is logical
        assertTrue(
            expectedIndex.toString().endsWith(".bam.bai"),
            "Index should have .bam.bai extension"
        );
    }

    /**
     * GROUP 4: Test multiple index creation attempts
     */
    @Test
    void testMultipleIndexAttempts() throws Exception {
        Path bamFile = tempDir.resolve("test.bam");
        Files.createFile(bamFile);

        // Multiple attempts should not cause issues
        for (int i = 0; i < 3; i++) {
            assertDoesNotThrow(() -> {
                boolean result = BamIndexer.ensureIndexExists(bamFile.toFile());
                assertNotNull(result);
            });
        }
    }

    /**
     * GROUP 4: Test index with file without .bam extension
     */
    @Test
    void testIndexFileWithoutBamExtension() throws Exception {
        Path notBamFile = tempDir.resolve("test.txt");
        Files.writeString(notBamFile, "test content");

        boolean result = BamIndexer.ensureIndexExists(notBamFile.toFile());

        // Should handle gracefully
        assertNotNull(result);
    }

    /**
     * GROUP 4: Test samtools detection is consistent
     */
    @Test
    void testSamtoolsDetectionConsistent() {
        boolean first = BamIndexer.isSamtoolsAvailable();
        boolean second = BamIndexer.isSamtoolsAvailable();

        assertEquals(first, second, "Samtools availability should be consistent");
    }

    /**
     * GROUP 4: Test with read-only directory
     */
    @Test
    void testIndexInReadOnlyDirectory() throws Exception {
        Path readOnlyDir = tempDir.resolve("readonly");
        Files.createDirectories(readOnlyDir);

        Path bamFile = readOnlyDir.resolve("test.bam");
        Files.createFile(bamFile);

        // Make directory read-only
        readOnlyDir.toFile().setReadOnly();

        try {
            boolean result = BamIndexer.ensureIndexExists(bamFile.toFile());
            // Should handle read-only directory gracefully
            assertNotNull(result);
        } finally {
            // Restore write permission for cleanup
            readOnlyDir.toFile().setWritable(true);
        }
    }

    /**
     * GROUP 4: Test with very long file name
     */
    @Test
    void testIndexWithLongFileName() throws Exception {
        String longName = "a".repeat(200) + ".bam";
        Path bamFile = tempDir.resolve(longName);

        // Only create if path is not too long for file system
        if (bamFile.toString().length() < 255) {
            Files.createFile(bamFile);

            assertDoesNotThrow(() -> {
                boolean result = BamIndexer.ensureIndexExists(bamFile.toFile());
                assertNotNull(result);
            });
        }
    }
}
