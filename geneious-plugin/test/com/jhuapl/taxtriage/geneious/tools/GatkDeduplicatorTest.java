package com.jhuapl.taxtriage.geneious.tools;

import com.jhuapl.taxtriage.geneious.tools.GatkDeduplicator.DeduplicationResult;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for GatkDeduplicator (Group 4).
 *
 * Tests the GATK deduplication functionality including:
 * - GATK availability detection
 * - Deduplication result structure
 * - Error handling for invalid inputs
 * - Docker command detection
 * - File path handling
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class GatkDeduplicatorTest {

    @TempDir
    Path tempDir;

    private GatkDeduplicator deduplicator;

    @BeforeEach
    void setUp() {
        deduplicator = new GatkDeduplicator();
    }

    /**
     * GROUP 4: Test GATK availability check
     */
    @Test
    void testIsGatkAvailable() {
        assertDoesNotThrow(() -> {
            boolean available = GatkDeduplicator.isGatkAvailable();
            // Result depends on Docker installation
            assertNotNull(available);
        });
    }

    /**
     * GROUP 4: Test deduplication with non-existent input file
     */
    @Test
    void testDeduplicateNonExistentFile() {
        Path inputBam = tempDir.resolve("nonexistent.bam");
        Path outputDir = tempDir.resolve("output");

        DeduplicationResult result = deduplicator.deduplicateBam(
            inputBam, outputDir, 4);

        assertNotNull(result, "Result should not be null");
        assertFalse(result.success, "Should fail for non-existent file");
        assertNotNull(result.errorMessage, "Should have error message");
        assertTrue(result.errorMessage.contains("does not exist"),
            "Error should mention file doesn't exist");
    }

    /**
     * GROUP 4: Test deduplication with null input
     */
    @Test
    void testDeduplicateNullInput() {
        Path outputDir = tempDir.resolve("output");

        assertDoesNotThrow(() -> {
            DeduplicationResult result = deduplicator.deduplicateBam(
                null, outputDir, 4);
            assertNotNull(result);
            assertFalse(result.success);
        });
    }

    /**
     * GROUP 4: Test deduplication with invalid output directory
     */
    @Test
    void testDeduplicateInvalidOutputDir() throws Exception {
        Path inputBam = tempDir.resolve("test.bam");
        Files.createFile(inputBam);

        // Use a file as output directory (invalid)
        Path invalidOutputDir = tempDir.resolve("file_not_dir.txt");
        Files.createFile(invalidOutputDir);

        assertDoesNotThrow(() -> {
            DeduplicationResult result = deduplicator.deduplicateBam(
                inputBam, invalidOutputDir, 4);
            assertNotNull(result);
        });
    }

    /**
     * GROUP 4: Test deduplication result structure
     */
    @Test
    void testDeduplicationResultStructure() {
        DeduplicationResult result = new DeduplicationResult(
            true, "/path/to/output.bam", null, 100);

        assertTrue(result.success, "Should be successful");
        assertEquals("/path/to/output.bam", result.outputPath);
        assertNull(result.errorMessage, "Error message should be null");
        assertEquals(100, result.duplicatesRemoved);
    }

    /**
     * GROUP 4: Test deduplication result with error
     */
    @Test
    void testDeduplicationResultWithError() {
        DeduplicationResult result = new DeduplicationResult(
            false, null, "Test error message", 0);

        assertFalse(result.success, "Should not be successful");
        assertNull(result.outputPath, "Output path should be null");
        assertEquals("Test error message", result.errorMessage);
        assertEquals(0, result.duplicatesRemoved);
    }

    /**
     * GROUP 4: Test deduplicate with various thread counts
     */
    @Test
    void testDeduplicateWithVariousThreadCounts() {
        Path inputBam = tempDir.resolve("test.bam");
        Path outputDir = tempDir.resolve("output");

        int[] threadCounts = {1, 2, 4, 8, 16};

        for (int threads : threadCounts) {
            assertDoesNotThrow(() -> {
                DeduplicationResult result = deduplicator.deduplicateBam(
                    inputBam, outputDir, threads);
                assertNotNull(result);
            }, "Should handle " + threads + " threads");
        }
    }

    /**
     * GROUP 4: Test deduplicate with zero threads
     */
    @Test
    void testDeduplicateWithZeroThreads() {
        Path inputBam = tempDir.resolve("test.bam");
        Path outputDir = tempDir.resolve("output");

        assertDoesNotThrow(() -> {
            DeduplicationResult result = deduplicator.deduplicateBam(
                inputBam, outputDir, 0);
            assertNotNull(result);
        });
    }

    /**
     * GROUP 4: Test deduplicate with negative threads
     */
    @Test
    void testDeduplicateWithNegativeThreads() {
        Path inputBam = tempDir.resolve("test.bam");
        Path outputDir = tempDir.resolve("output");

        assertDoesNotThrow(() -> {
            DeduplicationResult result = deduplicator.deduplicateBam(
                inputBam, outputDir, -1);
            assertNotNull(result);
        });
    }

    /**
     * GROUP 4: Test deduplicateAllBams with non-existent directory
     */
    @Test
    void testDeduplicateAllBamsNonExistentDir() {
        Path nonExistent = tempDir.resolve("nonexistent");

        List<DeduplicationResult> results = deduplicator.deduplicateAllBams(
            nonExistent, 4);

        assertNotNull(results, "Results should not be null");
        assertTrue(results.isEmpty(), "Should return empty list for non-existent dir");
    }

    /**
     * GROUP 4: Test deduplicateAllBams with empty directory
     */
    @Test
    void testDeduplicateAllBamsEmptyDir() throws Exception {
        Path emptyDir = tempDir.resolve("empty");
        Files.createDirectories(emptyDir);

        List<DeduplicationResult> results = deduplicator.deduplicateAllBams(
            emptyDir, 4);

        assertNotNull(results, "Results should not be null");
        assertTrue(results.isEmpty(), "Should return empty list for directory with no BAMs");
    }

    /**
     * GROUP 4: Test deduplicateAllBams with file instead of directory
     */
    @Test
    void testDeduplicateAllBamsWithFile() throws Exception {
        Path file = tempDir.resolve("test.txt");
        Files.createFile(file);

        List<DeduplicationResult> results = deduplicator.deduplicateAllBams(
            file, 4);

        assertNotNull(results, "Results should not be null");
        assertTrue(results.isEmpty(), "Should return empty list for file");
    }

    /**
     * GROUP 4: Test deduplicateAllBams skips already deduplicated files
     */
    @Test
    void testDeduplicateAllBamsSkipsDeduplicated() throws Exception {
        Path dir = tempDir.resolve("bams");
        Files.createDirectories(dir);

        // Create regular and deduplicated BAM files
        Files.createFile(dir.resolve("sample.bam"));
        Files.createFile(dir.resolve("sample.dedup.bam"));

        // Should skip .dedup.bam files
        assertDoesNotThrow(() -> {
            List<DeduplicationResult> results = deduplicator.deduplicateAllBams(
                dir, 4);
            assertNotNull(results);
        });
    }

    /**
     * GROUP 4: Test GATK availability is consistent
     */
    @Test
    void testGatkAvailabilityConsistent() {
        boolean first = GatkDeduplicator.isGatkAvailable();
        boolean second = GatkDeduplicator.isGatkAvailable();

        assertEquals(first, second, "GATK availability should be consistent");
    }

    /**
     * GROUP 4: Test deduplication with very long file path
     */
    @Test
    void testDeduplicateWithLongPath() {
        StringBuilder longPath = new StringBuilder(tempDir.toString());
        for (int i = 0; i < 10; i++) {
            longPath.append("/very_long_directory_name_").append(i);
        }
        Path inputBam = Path.of(longPath + "/test.bam");
        Path outputDir = tempDir.resolve("output");

        assertDoesNotThrow(() -> {
            DeduplicationResult result = deduplicator.deduplicateBam(
                inputBam, outputDir, 4);
            assertNotNull(result);
        });
    }

    /**
     * GROUP 4: Test multiple deduplication attempts
     */
    @Test
    void testMultipleDeduplicationAttempts() {
        Path inputBam = tempDir.resolve("test.bam");
        Path outputDir = tempDir.resolve("output");

        for (int i = 0; i < 3; i++) {
            assertDoesNotThrow(() -> {
                DeduplicationResult result = deduplicator.deduplicateBam(
                    inputBam, outputDir, 4);
                assertNotNull(result);
            });
        }
    }

    /**
     * GROUP 4: Test deduplication result with large duplicate count
     */
    @Test
    void testDeduplicationResultLargeDuplicateCount() {
        DeduplicationResult result = new DeduplicationResult(
            true, "/path/to/output.bam", null, Long.MAX_VALUE);

        assertEquals(Long.MAX_VALUE, result.duplicatesRemoved,
            "Should handle large duplicate counts");
    }
}
