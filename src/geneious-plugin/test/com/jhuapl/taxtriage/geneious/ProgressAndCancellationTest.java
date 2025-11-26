package com.jhuapl.taxtriage.geneious;

import jebl.util.ProgressListener;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for progress reporting and cancellation functionality (Group 2).
 *
 * Tests:
 * - Progress listener updates
 * - Progress message handling
 * - Cancellation detection
 * - Progress value boundaries
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class ProgressAndCancellationTest {

    private TestProgressListener progressListener;

    @BeforeEach
    void setUp() {
        progressListener = new TestProgressListener();
    }

    /**
     * GROUP 2: Test progress listener creation and initialization
     */
    @Test
    void testProgressListenerCreation() {
        assertNotNull(progressListener, "Progress listener should be created");
        assertEquals(0.0, progressListener.getLastProgress(), 0.001,
            "Initial progress should be 0.0");
        assertEquals("", progressListener.getLastMessage(),
            "Initial message should be empty");
        assertFalse(progressListener.isCanceled(),
            "Should not be cancelled initially");
    }

    /**
     * GROUP 2: Test setting progress values
     */
    @Test
    void testSetProgress() {
        progressListener.setProgress(0.0);
        assertEquals(0.0, progressListener.getLastProgress(), 0.001);

        progressListener.setProgress(0.5);
        assertEquals(0.5, progressListener.getLastProgress(), 0.001);

        progressListener.setProgress(1.0);
        assertEquals(1.0, progressListener.getLastProgress(), 0.001);
    }

    /**
     * GROUP 2: Test progress message updates
     */
    @Test
    void testSetMessage() {
        progressListener.setMessage("Starting analysis...");
        assertEquals("Starting analysis...", progressListener.getLastMessage());

        progressListener.setMessage("Processing samples...");
        assertEquals("Processing samples...", progressListener.getLastMessage());

        progressListener.setMessage("Analysis complete");
        assertEquals("Analysis complete", progressListener.getLastMessage());
    }

    /**
     * GROUP 2: Test progress update counting
     */
    @Test
    void testProgressUpdateCounting() {
        assertEquals(0, progressListener.getProgressUpdateCount(),
            "Should start with 0 updates");

        progressListener.setProgress(0.1);
        assertEquals(1, progressListener.getProgressUpdateCount());

        progressListener.setProgress(0.5);
        assertEquals(2, progressListener.getProgressUpdateCount());

        progressListener.setProgress(1.0);
        assertEquals(3, progressListener.getProgressUpdateCount());
    }

    /**
     * GROUP 2: Test cancellation detection
     */
    @Test
    void testCancellation() {
        assertFalse(progressListener.isCanceled(),
            "Should not be cancelled initially");

        // BasicProgressListener has a cancel() method
        progressListener.cancel();
        assertTrue(progressListener.isCanceled(),
            "Should be cancelled after cancel() called");
    }

    /**
     * GROUP 2: Test progress sequence for typical workflow
     */
    @Test
    void testTypicalWorkflowProgressSequence() {
        // Simulate typical workflow progress
        progressListener.setMessage("Initializing...");
        progressListener.setProgress(0.0);

        progressListener.setMessage("Checking system requirements...");
        progressListener.setProgress(0.05);

        progressListener.setMessage("Setting up workspace...");
        progressListener.setProgress(0.1);

        progressListener.setMessage("Starting TaxTriage workflow...");
        progressListener.setProgress(0.3);

        progressListener.setMessage("Running analysis...");
        progressListener.setProgress(0.6);

        progressListener.setMessage("Processing results...");
        progressListener.setProgress(0.9);

        progressListener.setMessage("Complete");
        progressListener.setProgress(1.0);

        assertEquals(1.0, progressListener.getLastProgress(), 0.001);
        assertEquals("Complete", progressListener.getLastMessage());
        assertEquals(7, progressListener.getProgressUpdateCount());
    }

    /**
     * GROUP 2: Test progress boundaries (0.0 to 1.0)
     */
    @Test
    void testProgressBoundaries() {
        // Test minimum
        progressListener.setProgress(0.0);
        assertEquals(0.0, progressListener.getLastProgress(), 0.001);

        // Test maximum
        progressListener.setProgress(1.0);
        assertEquals(1.0, progressListener.getLastProgress(), 0.001);

        // Test intermediate values
        progressListener.setProgress(0.25);
        assertEquals(0.25, progressListener.getLastProgress(), 0.001);

        progressListener.setProgress(0.75);
        assertEquals(0.75, progressListener.getLastProgress(), 0.001);
    }

    /**
     * GROUP 2: Test empty message handling
     */
    @Test
    void testEmptyMessage() {
        progressListener.setMessage("");
        assertEquals("", progressListener.getLastMessage());

        progressListener.setMessage("Non-empty message");
        assertEquals("Non-empty message", progressListener.getLastMessage());

        progressListener.setMessage("");
        assertEquals("", progressListener.getLastMessage());
    }

    /**
     * GROUP 2: Test null message handling
     */
    @Test
    void testNullMessage() {
        assertDoesNotThrow(() -> {
            progressListener.setMessage(null);
        }, "Should handle null message gracefully");
    }

    /**
     * GROUP 2: Test progress with concurrent updates
     */
    @Test
    void testConcurrentProgressUpdates() {
        progressListener.setProgress(0.1);
        progressListener.setMessage("Message 1");

        progressListener.setProgress(0.2);
        progressListener.setMessage("Message 2");

        progressListener.setProgress(0.3);
        progressListener.setMessage("Message 3");

        assertEquals(0.3, progressListener.getLastProgress(), 0.001);
        assertEquals("Message 3", progressListener.getLastMessage());
        assertEquals(3, progressListener.getProgressUpdateCount());
    }

    /**
     * GROUP 2: Test cancellation during progress
     */
    @Test
    void testCancellationDuringProgress() {
        progressListener.setProgress(0.0);
        progressListener.setMessage("Starting...");
        assertFalse(progressListener.isCanceled());

        progressListener.setProgress(0.5);
        progressListener.setMessage("Halfway...");
        assertFalse(progressListener.isCanceled());

        // Cancel mid-operation
        progressListener.cancel();
        assertTrue(progressListener.isCanceled());

        // Progress can still be updated after cancellation
        progressListener.setProgress(0.6);
        progressListener.setMessage("After cancellation");
        assertEquals(0.6, progressListener.getLastProgress(), 0.001);
        assertTrue(progressListener.isCanceled(),
            "Should remain cancelled");
    }

    /**
     * GROUP 2: Test multiple message updates without progress change
     */
    @Test
    void testMultipleMessagesWithSameProgress() {
        progressListener.setProgress(0.5);
        assertEquals(1, progressListener.getProgressUpdateCount());

        progressListener.setMessage("Message 1");
        progressListener.setMessage("Message 2");
        progressListener.setMessage("Message 3");

        assertEquals(0.5, progressListener.getLastProgress(), 0.001);
        assertEquals("Message 3", progressListener.getLastMessage());
        assertEquals(1, progressListener.getProgressUpdateCount(),
            "Progress count should not change with just message updates");
    }

    /**
     * GROUP 2: Test workflow-specific progress stages
     */
    @Test
    void testWorkflowProgressStages() {
        // Initial validation (0-10%)
        progressListener.setProgress(0.05);
        progressListener.setMessage("Validating inputs...");
        assertTrue(progressListener.getLastProgress() < 0.1);

        // Database check (10-15%)
        progressListener.setProgress(0.12);
        progressListener.setMessage("Checking database...");
        assertTrue(progressListener.getLastProgress() >= 0.1 &&
                  progressListener.getLastProgress() < 0.15);

        // Workspace setup (15-20%)
        progressListener.setProgress(0.18);
        progressListener.setMessage("Setting up workspace...");
        assertTrue(progressListener.getLastProgress() >= 0.15 &&
                  progressListener.getLastProgress() < 0.2);

        // Workflow execution (30-85%)
        progressListener.setProgress(0.6);
        progressListener.setMessage("Running workflow...");
        assertTrue(progressListener.getLastProgress() >= 0.3 &&
                  progressListener.getLastProgress() < 0.85);

        // Result import (90-100%)
        progressListener.setProgress(0.95);
        progressListener.setMessage("Importing results...");
        assertTrue(progressListener.getLastProgress() >= 0.9);
    }

    /**
     * GROUP 2: Test ProgressListener.EMPTY constant
     */
    @Test
    void testEmptyProgressListener() {
        ProgressListener empty = ProgressListener.EMPTY;
        assertNotNull(empty, "EMPTY listener should not be null");

        // Should not throw exceptions
        assertDoesNotThrow(() -> {
            empty.setProgress(0.5);
            empty.setMessage("Test");
            boolean cancelled = empty.isCanceled();
        }, "EMPTY listener should handle calls gracefully");
    }

    // ===== Test Implementation of ProgressListener =====

    /**
     * Test implementation of ProgressListener for testing
     */
    private static class TestProgressListener extends jebl.util.BasicProgressListener {
        private int progressUpdateCount = 0;

        @Override
        protected void _setProgress(double fractionCompleted) {
            progressUpdateCount++;
            super._setProgress(fractionCompleted);
        }

        public int getProgressUpdateCount() {
            return progressUpdateCount;
        }

        public double getLastProgress() {
            return getFractionCompleted();
        }

        public String getLastMessage() {
            return getMessage() != null ? getMessage() : "";
        }
    }
}
