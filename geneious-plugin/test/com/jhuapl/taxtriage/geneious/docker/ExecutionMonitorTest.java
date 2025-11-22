package com.jhuapl.taxtriage.geneious.docker;

import jebl.util.ProgressListener;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for ExecutionMonitor class focusing on enhanced progress tracking (Group 2).
 *
 * Tests:
 * - Process name extraction from various Nextflow output formats
 * - Process name formatting for display
 * - Progress message updates with process details
 * - State tracking across multiple processes
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class ExecutionMonitorTest {

    private ExecutionMonitor monitor;
    private TestProgressListener progressListener;

    @BeforeEach
    void setUp() {
        monitor = new ExecutionMonitor();
        progressListener = new TestProgressListener();
    }

    /**
     * GROUP 2: Test extracting process names from standard Nextflow progress lines
     */
    @Test
    void testExtractProcessNameFromProgressLine() {
        String line = "[12/abc123] process > FASTQC [1 of 3]";
        String processName = monitor.extractProcessName(line);
        assertEquals("FASTQC", processName, "Should extract FASTQC process name");
    }

    /**
     * GROUP 2: Test extracting process names from modern Nextflow format
     */
    @Test
    void testExtractProcessNameFromModernFormat() {
        String line = "[75/679d58] NFCORE_TAXTRIAGE:CLASSIFIER:KRAKEN2_KRAKEN2 (02-AE00001037F538_S2_L002_dedupe) [100%] 1 of 1 ✔";
        String processName = monitor.extractProcessName(line);
        assertEquals("NFCORE_TAXTRIAGE:CLASSIFIER:KRAKEN2_KRAKEN2", processName,
            "Should extract full process path from modern format");
    }

    /**
     * GROUP 2: Test parsing progress from modern Nextflow format
     */
    @Test
    void testParseProgressFromModernFormat() {
        String line = "[bc/98af0e] NFCORE_TAXTRIAGE:TAXTRIAGE:FASTQC (02-AE00001037F538_S2_L002_dedupe) [100%] 1 of 1 ✔";
        double progress = monitor.parseProgress(line);
        assertEquals(1.0, progress, 0.01, "Should parse progress from modern format");
    }

    /**
     * GROUP 2: Test extracting process names from executor submission lines
     */
    @Test
    void testExtractProcessNameFromSubmittedLine() {
        String line = "[12/abc123] Submitted process > KRAKEN2_KRAKEN2 (sample1)";
        String processName = monitor.extractProcessName(line);
        assertEquals("KRAKEN2_KRAKEN2", processName, "Should extract KRAKEN2_KRAKEN2 process name");
    }

    /**
     * GROUP 2: Test extracting process names with colons (module prefixes)
     */
    @Test
    void testExtractProcessNameWithPrefix() {
        String line = "[12/abc123] process > NFCORE_TAXTRIAGE:FASTQC [1 of 2]";
        String processName = monitor.extractProcessName(line);
        assertEquals("NFCORE_TAXTRIAGE:FASTQC", processName, "Should extract full process name with prefix");
    }

    /**
     * GROUP 2: Test extracting process names from lines without process info
     */
    @Test
    void testExtractProcessNameFromNonProcessLine() {
        String line = "Launching workflow...";
        String processName = monitor.extractProcessName(line);
        assertNull(processName, "Should return null for non-process lines");
    }

    /**
     * GROUP 2: Test extracting process names from empty or null lines
     */
    @Test
    void testExtractProcessNameFromEmptyLine() {
        assertNull(monitor.extractProcessName(""), "Should return null for empty line");
        assertNull(monitor.extractProcessName(null), "Should return null for null line");
    }

    /**
     * GROUP 2: Test progress parsing from Nextflow output
     */
    @Test
    void testParseProgressFromNextflowLine() {
        String line = "[12/abc123] process > FASTQC [2 of 4]";
        double progress = monitor.parseProgress(line);
        assertEquals(0.5, progress, 0.01, "Should parse progress as 2/4 = 0.5");
    }

    /**
     * GROUP 2: Test progress parsing with different ratios
     */
    @Test
    void testParseProgressVariousRatios() {
        assertEquals(0.25, monitor.parseProgress("[x] process > TEST [1 of 4]"), 0.01);
        assertEquals(0.33, monitor.parseProgress("[x] process > TEST [1 of 3]"), 0.01);
        assertEquals(1.0, monitor.parseProgress("[x] process > TEST [5 of 5]"), 0.01);
    }

    /**
     * GROUP 2: Test progress parsing from completion line
     */
    @Test
    void testParseProgressFromCompletionLine() {
        String line = "Pipeline completed at: 2024-01-15T10:30:00Z";
        double progress = monitor.parseProgress(line);
        assertEquals(1.0, progress, 0.01, "Should return 1.0 for completion line");
    }

    /**
     * GROUP 2: Test progress parsing from non-progress lines
     */
    @Test
    void testParseProgressFromNonProgressLine() {
        assertEquals(-1.0, monitor.parseProgress("Some random output"), 0.01);
        assertEquals(-1.0, monitor.parseProgress(""), 0.01);
        assertEquals(-1.0, monitor.parseProgress(null), 0.01);
    }

    /**
     * GROUP 2: Test error line detection
     */
    @Test
    void testIsErrorLine() {
        assertTrue(monitor.isErrorLine("ERROR: Something went wrong"));
        assertTrue(monitor.isErrorLine("WARN: Deprecation notice"));
        assertTrue(monitor.isErrorLine("Exception in thread"));
        assertTrue(monitor.isErrorLine("Failed to execute"));
        assertFalse(monitor.isErrorLine("Normal output line"));
        assertFalse(monitor.isErrorLine(""));
        assertFalse(monitor.isErrorLine(null));
    }

    /**
     * GROUP 2: Test completion line detection
     */
    @Test
    void testIsCompletionLine() {
        assertTrue(monitor.isCompletionLine("Pipeline completed at: 2024-01-15"));
        assertTrue(monitor.isCompletionLine("Completed at: 10:30:00"));
        assertTrue(monitor.isCompletionLine("Duration: 1h 30m"));
        assertFalse(monitor.isCompletionLine("Still running..."));
        assertFalse(monitor.isCompletionLine(""));
    }

    /**
     * GROUP 2: Test current process name tracking
     */
    @Test
    void testCurrentProcessTracking() {
        assertNull(monitor.getCurrentProcessName(), "Initially should be null");

        monitor.extractProcessName("[x] process > FASTQC [1 of 1]");
        assertEquals("FASTQC", monitor.getCurrentProcessName(), "Should track current process");

        monitor.extractProcessName("[x] process > KRAKEN2 [1 of 1]");
        assertEquals("KRAKEN2", monitor.getCurrentProcessName(), "Should update to new process");
    }

    /**
     * GROUP 2: Test process tracking reset
     */
    @Test
    void testResetProcessTracking() {
        monitor.extractProcessName("[x] process > FASTQC [1 of 1]");
        assertNotNull(monitor.getCurrentProcessName());

        monitor.resetProcessTracking();
        assertNull(monitor.getCurrentProcessName(), "Should reset to null");
    }

    /**
     * GROUP 2: Test monitor reset includes process tracking
     */
    @Test
    void testMonitorResetIncludesProcessTracking() {
        monitor.extractProcessName("[x] process > FASTQC [1 of 1]");
        assertNotNull(monitor.getCurrentProcessName());

        monitor.reset();
        assertNull(monitor.getCurrentProcessName(), "Reset should clear process tracking");
        assertFalse(monitor.isCancellationRequested(), "Reset should clear cancellation");
    }

    /**
     * GROUP 2: Test multiple process name formats
     */
    @Test
    void testMultipleProcessNameFormats() {
        // Standard format
        assertEquals("FASTQC", monitor.extractProcessName("[x] process > FASTQC [1 of 1]"));

        // With module prefix
        monitor.resetProcessTracking();
        assertEquals("MODULE:PROCESS",
            monitor.extractProcessName("[x] process > MODULE:PROCESS [1 of 1]"));

        // Submitted format
        monitor.resetProcessTracking();
        assertEquals("KRAKEN2",
            monitor.extractProcessName("[x] Submitted process > KRAKEN2 (sample)"));
    }

    /**
     * GROUP 2: Test process name extraction consistency
     */
    @Test
    void testProcessNameExtractionConsistency() {
        String line = "[12/abc] process > TEST_PROCESS [1 of 1]";

        String first = monitor.extractProcessName(line);
        String second = monitor.extractProcessName(line);

        assertEquals(first, second, "Should extract same name consistently");
        assertEquals("TEST_PROCESS", first);
    }

    /**
     * GROUP 2: Test cancellation state
     */
    @Test
    void testCancellationState() {
        assertFalse(monitor.isCancellationRequested(), "Initially not cancelled");

        monitor.requestCancellation();
        assertTrue(monitor.isCancellationRequested(), "Should be cancelled");

        monitor.reset();
        assertFalse(monitor.isCancellationRequested(), "Should reset cancellation");
    }

    /**
     * GROUP 2: Test parsing progress with invalid numbers
     */
    @Test
    void testParseProgressWithInvalidNumbers() {
        assertEquals(-1.0, monitor.parseProgress("[x] process > TEST [x of y]"), 0.01);
        assertEquals(-1.0, monitor.parseProgress("[x] process > TEST [0 of 0]"), 0.01);
    }

    /**
     * GROUP 2: Test process name extraction from complex lines
     */
    @Test
    void testExtractProcessFromComplexLines() {
        // Test with parentheses (sample identifiers) - these don't match the standard pattern
        // but should work with modified pattern if needed
        String line1 = "[12/abc123] process > NFCORE_TAXTRIAGE_QC_FASTQC [1 of 10]";
        String processName1 = monitor.extractProcessName(line1);
        assertNotNull(processName1, "Should extract process from line without parentheses");
        assertEquals("NFCORE_TAXTRIAGE_QC_FASTQC", processName1);

        // Test submitted format with sample name in parentheses
        String line2 = "[12/abc123] Submitted process > FASTQC (sample_R1)";
        String processName2 = monitor.extractProcessName(line2);
        assertNotNull(processName2, "Should extract process from submitted line");
        assertEquals("FASTQC", processName2);
    }

    /**
     * GROUP 2: Test sequential process tracking
     */
    @Test
    void testSequentialProcessTracking() {
        // Simulate workflow progress through multiple processes
        monitor.extractProcessName("[x] process > FASTQC [1 of 1]");
        assertEquals("FASTQC", monitor.getCurrentProcessName());

        monitor.extractProcessName("[x] process > TRIMMING [1 of 1]");
        assertEquals("TRIMMING", monitor.getCurrentProcessName());

        monitor.extractProcessName("[x] process > KRAKEN2 [1 of 1]");
        assertEquals("KRAKEN2", monitor.getCurrentProcessName());
    }

    /**
     * GROUP 2: Test progress parsing edge cases
     */
    @Test
    void testProgressParsingEdgeCases() {
        // Zero total should not cause division by zero
        assertEquals(-1.0, monitor.parseProgress("[x] process > TEST [0 of 0]"), 0.01);

        // Completed before started (should not happen but handle gracefully)
        assertEquals(1.0, monitor.parseProgress("[x] process > TEST [5 of 5]"), 0.01);

        // Single task
        assertEquals(1.0, monitor.parseProgress("[x] process > TEST [1 of 1]"), 0.01);
    }

    // ===== Test Implementation of ProgressListener =====

    /**
     * Test implementation of ProgressListener for testing
     */
    private static class TestProgressListener extends jebl.util.BasicProgressListener {
        private String lastMessage = "";

        @Override
        protected void _setMessage(String message) {
            this.lastMessage = message != null ? message : "";
            super._setMessage(message);
        }

        public String getLastMessage() {
            return lastMessage;
        }

        public double getLastProgress() {
            return getFractionCompleted();
        }
    }
}
