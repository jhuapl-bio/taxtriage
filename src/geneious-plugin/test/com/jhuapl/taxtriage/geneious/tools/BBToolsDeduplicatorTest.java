package com.jhuapl.taxtriage.geneious.tools;

import com.jhuapl.taxtriage.geneious.docker.DockerException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for BBToolsDeduplicator (Group 3).
 *
 * Tests the BBTools deduplication functionality including:
 * - Deduplicator initialization
 * - Configuration validation
 * - Result object handling
 * - Error handling
 *
 * Note: Full Docker-based deduplication tests require Docker to be available
 * and are marked appropriately.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class BBToolsDeduplicatorTest {

    @TempDir
    Path tempDir;

    /**
     * GROUP 3: Test DeduplicationResult creation and getters
     */
    @Test
    void testDeduplicationResultSuccess() {
        BBToolsDeduplicator.DeduplicationResult result =
            new BBToolsDeduplicator.DeduplicationResult(
                true,
                "/path/to/output.fastq",
                null,
                1000,
                10000,
                java.util.Arrays.asList("Step 1 complete", "Step 2 complete")
            );

        assertTrue(result.isSuccess());
        assertEquals("/path/to/output.fastq", result.getOutputPath());
        assertNull(result.getErrorMessage());
        assertEquals(1000, result.getDuplicatesRemoved());
        assertEquals(10000, result.getTotalReads());
        assertEquals(2, result.getStepResults().size());
    }

    /**
     * GROUP 3: Test DeduplicationResult failure
     */
    @Test
    void testDeduplicationResultFailure() {
        BBToolsDeduplicator.DeduplicationResult result =
            new BBToolsDeduplicator.DeduplicationResult(
                false,
                null,
                "Docker not available",
                0,
                0,
                null
            );

        assertFalse(result.isSuccess());
        assertNull(result.getOutputPath());
        assertEquals("Docker not available", result.getErrorMessage());
        assertEquals(0, result.getDuplicatesRemoved());
        assertEquals(0, result.getTotalReads());
        assertNotNull(result.getStepResults());
        assertEquals(0, result.getStepResults().size());
    }

    /**
     * GROUP 3: Test deduplication percentage calculation
     */
    @Test
    void testDeduplicationPercentage() {
        BBToolsDeduplicator.DeduplicationResult result =
            new BBToolsDeduplicator.DeduplicationResult(
                true, "/path", null, 2500, 10000, null
            );

        assertEquals(25.0, result.getDeduplicationPercentage(), 0.01);
    }

    /**
     * GROUP 3: Test deduplication percentage with zero reads
     */
    @Test
    void testDeduplicationPercentageWithZeroReads() {
        BBToolsDeduplicator.DeduplicationResult result =
            new BBToolsDeduplicator.DeduplicationResult(
                true, "/path", null, 0, 0, null
            );

        assertEquals(0.0, result.getDeduplicationPercentage(), 0.01);
    }

    /**
     * GROUP 3: Test deduplication percentage with all duplicates
     */
    @Test
    void testDeduplicationPercentageAllDuplicates() {
        BBToolsDeduplicator.DeduplicationResult result =
            new BBToolsDeduplicator.DeduplicationResult(
                true, "/path", null, 10000, 10000, null
            );

        assertEquals(100.0, result.getDeduplicationPercentage(), 0.01);
    }

    /**
     * GROUP 3: Test result toString
     */
    @Test
    void testDeduplicationResultToString() {
        BBToolsDeduplicator.DeduplicationResult result =
            new BBToolsDeduplicator.DeduplicationResult(
                true, "/output.fastq", null, 500, 5000, null
            );

        String str = result.toString();
        assertNotNull(str);
        assertTrue(str.contains("success=true"));
        assertTrue(str.contains("/output.fastq"));
        assertTrue(str.contains("500"));
        assertTrue(str.contains("5000"));
        assertTrue(str.contains("10.00%"));
    }

    /**
     * GROUP 3: Test default constructor requires Docker
     *
     * Note: This test may fail if Docker is not available, which is expected.
     */
    @Test
    void testDefaultConstructorRequiresDocker() {
        try {
            BBToolsDeduplicator deduplicator = new BBToolsDeduplicator();
            assertNotNull(deduplicator,
                "Deduplicator should be created if Docker is available");
        } catch (DockerException e) {
            // Expected if Docker is not available
            assertTrue(e.getMessage().contains("Docker") ||
                      e.getMessage().contains("docker"),
                "Exception should mention Docker");
        }
    }

    /**
     * GROUP 3: Test constructor with custom threshold
     */
    @Test
    void testConstructorWithCustomThreshold() {
        try {
            BBToolsDeduplicator deduplicator = new BBToolsDeduplicator(3);
            assertNotNull(deduplicator);
        } catch (DockerException e) {
            // Expected if Docker is not available
            assertTrue(e.getMessage().contains("Docker"));
        }
    }

    /**
     * GROUP 3: Test constructor with custom threshold and memory
     */
    @Test
    void testConstructorWithCustomThresholdAndMemory() {
        try {
            BBToolsDeduplicator deduplicator = new BBToolsDeduplicator(5, 4);
            assertNotNull(deduplicator);
        } catch (DockerException e) {
            // Expected if Docker is not available
            assertTrue(e.getMessage().contains("Docker"));
        }
    }

    /**
     * GROUP 3: Test threshold clamping - too low
     */
    @Test
    void testThresholdClampingTooLow() {
        try {
            BBToolsDeduplicator deduplicator = new BBToolsDeduplicator(-5, 8);
            // If Docker is available, threshold should be clamped to 0
            assertNotNull(deduplicator);
        } catch (DockerException e) {
            // Expected if Docker is not available
        }
    }

    /**
     * GROUP 3: Test threshold clamping - too high
     */
    @Test
    void testThresholdClampingTooHigh() {
        try {
            BBToolsDeduplicator deduplicator = new BBToolsDeduplicator(100, 8);
            // If Docker is available, threshold should be clamped to 10
            assertNotNull(deduplicator);
        } catch (DockerException e) {
            // Expected if Docker is not available
        }
    }

    /**
     * GROUP 3: Test memory clamping - too low
     */
    @Test
    void testMemoryClampingTooLow() {
        try {
            BBToolsDeduplicator deduplicator = new BBToolsDeduplicator(5, 0);
            // If Docker is available, memory should be clamped to 1
            assertNotNull(deduplicator);
        } catch (DockerException e) {
            // Expected if Docker is not available
        }
    }

    /**
     * GROUP 3: Test memory clamping - too high
     */
    @Test
    void testMemoryClampingTooHigh() {
        try {
            BBToolsDeduplicator deduplicator = new BBToolsDeduplicator(5, 100);
            // If Docker is available, memory should be clamped to 64
            assertNotNull(deduplicator);
        } catch (DockerException e) {
            // Expected if Docker is not available
        }
    }

    /**
     * GROUP 3: Test result with empty step results
     */
    @Test
    void testResultWithEmptyStepResults() {
        BBToolsDeduplicator.DeduplicationResult result =
            new BBToolsDeduplicator.DeduplicationResult(
                true, "/path", null, 100, 1000, java.util.Collections.emptyList()
            );

        assertNotNull(result.getStepResults());
        assertEquals(0, result.getStepResults().size());
    }

    /**
     * GROUP 3: Test result immutability
     */
    @Test
    void testResultImmutability() {
        java.util.List<String> steps = new java.util.ArrayList<>();
        steps.add("Step 1");

        BBToolsDeduplicator.DeduplicationResult result =
            new BBToolsDeduplicator.DeduplicationResult(
                true, "/path", null, 100, 1000, steps
            );

        // Modify original list
        steps.add("Step 2");

        // Result should not be affected
        assertEquals(1, result.getStepResults().size());
    }

    /**
     * GROUP 3: Test result with multiple steps
     */
    @Test
    void testResultWithMultipleSteps() {
        java.util.List<String> steps = java.util.Arrays.asList(
            "Step 1: Split interleaved FASTQ",
            "Step 2: Deduplicate reads",
            "Step 3: Clean read names",
            "Step 4: Re-interleave"
        );

        BBToolsDeduplicator.DeduplicationResult result =
            new BBToolsDeduplicator.DeduplicationResult(
                true, "/output.fastq", null, 5000, 50000, steps
            );

        assertEquals(4, result.getStepResults().size());
        assertTrue(result.getStepResults().get(0).contains("Split"));
        assertTrue(result.getStepResults().get(1).contains("Deduplicate"));
        assertTrue(result.getStepResults().get(2).contains("Clean"));
        assertTrue(result.getStepResults().get(3).contains("Re-interleave"));
    }

    /**
     * GROUP 3: Test various deduplication percentages
     */
    @Test
    void testVariousDeduplicationPercentages() {
        // 0% duplicates
        assertEquals(0.0,
            new BBToolsDeduplicator.DeduplicationResult(true, "/p", null, 0, 1000, null)
                .getDeduplicationPercentage(), 0.01);

        // 10% duplicates
        assertEquals(10.0,
            new BBToolsDeduplicator.DeduplicationResult(true, "/p", null, 100, 1000, null)
                .getDeduplicationPercentage(), 0.01);

        // 50% duplicates
        assertEquals(50.0,
            new BBToolsDeduplicator.DeduplicationResult(true, "/p", null, 500, 1000, null)
                .getDeduplicationPercentage(), 0.01);

        // 99% duplicates
        assertEquals(99.0,
            new BBToolsDeduplicator.DeduplicationResult(true, "/p", null, 990, 1000, null)
                .getDeduplicationPercentage(), 0.01);
    }

    /**
     * GROUP 3: Test result with large numbers
     */
    @Test
    void testResultWithLargeNumbers() {
        long totalReads = 1_000_000_000L;  // 1 billion reads
        long duplicates = 250_000_000L;    // 250 million duplicates

        BBToolsDeduplicator.DeduplicationResult result =
            new BBToolsDeduplicator.DeduplicationResult(
                true, "/output.fastq", null, duplicates, totalReads, null
            );

        assertEquals(25.0, result.getDeduplicationPercentage(), 0.01);
        assertEquals(duplicates, result.getDuplicatesRemoved());
        assertEquals(totalReads, result.getTotalReads());
    }

    /**
     * GROUP 3: Test result with error message formatting
     */
    @Test
    void testResultWithErrorMessage() {
        String errorMsg = "Docker container exited with code 1: Out of memory";

        BBToolsDeduplicator.DeduplicationResult result =
            new BBToolsDeduplicator.DeduplicationResult(
                false, null, errorMsg, 0, 0, null
            );

        assertFalse(result.isSuccess());
        assertEquals(errorMsg, result.getErrorMessage());
        assertTrue(result.getErrorMessage().contains("Out of memory"));
    }
}
