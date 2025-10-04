package com.jhuapl.taxtriage.geneious;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentOperationException;
import com.biomatters.geneious.publicapi.plugin.DocumentSelectionSignature;
import com.biomatters.geneious.publicapi.plugin.Options;
import jebl.util.ProgressListener;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for TaxTriageImportOperation (Group 3).
 *
 * Tests the import operation functionality including:
 * - Operation initialization and properties
 * - Options creation and configuration
 * - Input validation
 * - Error handling for missing directories
 * - Result type selection
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class TaxTriageImportOperationTest {

    private TaxTriageImportOperation operation;

    @TempDir
    Path tempDir;

    @BeforeEach
    void setUp() {
        operation = new TaxTriageImportOperation();
    }

    /**
     * GROUP 3: Test operation initialization
     */
    @Test
    void testOperationInitialization() {
        assertNotNull(operation, "Operation should be initialized");
        assertNotNull(operation.getHelp(), "Help text should not be null");
        assertNotNull(operation.getSelectionSignatures(), "Selection signatures should not be null");
    }

    /**
     * GROUP 3: Test help text content
     */
    @Test
    void testHelpText() {
        String help = operation.getHelp();

        assertNotNull(help, "Help text should not be null");
        assertFalse(help.isEmpty(), "Help text should not be empty");
        assertTrue(help.toLowerCase().contains("import"), "Help should mention importing");
        assertTrue(help.toLowerCase().contains("taxtriage"), "Help should mention TaxTriage");
    }

    /**
     * GROUP 3: Test selection signatures - should not require documents
     */
    @Test
    void testSelectionSignatures() {
        DocumentSelectionSignature[] signatures = operation.getSelectionSignatures();

        assertNotNull(signatures, "Selection signatures should not be null");
        // Import operation doesn't require document selection
        assertEquals(0, signatures.length,
            "Import operation should not require document selection");
    }

    /**
     * GROUP 3: Test options creation
     */
    @Test
    void testOptionsCreation() throws Exception {
        Options options = operation.getOptions(new AnnotatedPluginDocument[0]);

        assertNotNull(options, "Options should be created");

        // Verify required options exist
        assertNotNull(options.getOption("outputPath"),
            "Should have output path option");

        // Verify result type options
        assertNotNull(options.getOption("importReports"),
            "Should have import reports option");
        assertNotNull(options.getOption("importTopHits"),
            "Should have import top hits option");
        assertNotNull(options.getOption("importFilteredReports"),
            "Should have import filtered reports option");
        assertNotNull(options.getOption("importAlignments"),
            "Should have import alignments option");
        assertNotNull(options.getOption("importQualityReports"),
            "Should have import QC reports option");
    }

    /**
     * GROUP 3: Test default option values
     */
    @Test
    void testDefaultOptionValues() throws Exception {
        Options options = operation.getOptions(new AnnotatedPluginDocument[0]);

        // Check default values
        String defaultPath = (String) options.getValue("outputPath");
        assertNotNull(defaultPath, "Default path should be set");
        assertTrue(defaultPath.contains("taxtriage_output"),
            "Default path should contain taxtriage_output");

        // All result types should be enabled by default
        assertTrue((Boolean) options.getValue("importReports"),
            "Import reports should be enabled by default");
        assertTrue((Boolean) options.getValue("importTopHits"),
            "Import top hits should be enabled by default");
        assertTrue((Boolean) options.getValue("importFilteredReports"),
            "Import filtered reports should be enabled by default");
        assertTrue((Boolean) options.getValue("importAlignments"),
            "Import alignments should be enabled by default");
        assertTrue((Boolean) options.getValue("importQualityReports"),
            "Import QC reports should be enabled by default");
    }

    /**
     * GROUP 3: Test operation with null output path
     */
    @Test
    void testPerformOperationWithNullPath() throws Exception {
        Options options = operation.getOptions(new AnnotatedPluginDocument[0]);
        options.setValue("outputPath", null);

        TestProgressListener progressListener = new TestProgressListener();

        assertThrows(DocumentOperationException.class, () -> {
            operation.performOperation(new AnnotatedPluginDocument[0],
                progressListener, options);
        }, "Operation should throw exception with null output path");
    }

    /**
     * GROUP 3: Test operation with empty output path
     */
    @Test
    void testPerformOperationWithEmptyPath() throws Exception {
        Options options = operation.getOptions(new AnnotatedPluginDocument[0]);
        options.setValue("outputPath", "");

        TestProgressListener progressListener = new TestProgressListener();

        assertThrows(DocumentOperationException.class, () -> {
            operation.performOperation(new AnnotatedPluginDocument[0],
                progressListener, options);
        }, "Operation should throw exception with empty output path");
    }

    /**
     * GROUP 3: Test operation with non-existent directory
     */
    @Test
    void testPerformOperationWithNonExistentDirectory() throws Exception {
        Options options = operation.getOptions(new AnnotatedPluginDocument[0]);
        String nonExistentPath = tempDir.resolve("does_not_exist").toString();
        options.setValue("outputPath", nonExistentPath);

        TestProgressListener progressListener = new TestProgressListener();

        DocumentOperationException exception = assertThrows(
            DocumentOperationException.class, () -> {
                operation.performOperation(new AnnotatedPluginDocument[0],
                    progressListener, options);
            }, "Operation should throw exception for non-existent directory");

        assertTrue(exception.getMessage().contains("does not exist") ||
                   exception.getMessage().contains("not exist"),
            "Error message should mention directory doesn't exist");
    }

    /**
     * GROUP 3: Test operation with valid empty directory
     */
    @Test
    void testPerformOperationWithEmptyDirectory() throws Exception {
        // Create empty directory
        Path emptyDir = tempDir.resolve("empty_output");
        Files.createDirectories(emptyDir);

        Options options = operation.getOptions(new AnnotatedPluginDocument[0]);
        options.setValue("outputPath", emptyDir.toString());

        TestProgressListener progressListener = new TestProgressListener();

        // Should not throw exception, but will return empty list
        List<AnnotatedPluginDocument> results = operation.performOperation(
            new AnnotatedPluginDocument[0], progressListener, options);

        assertNotNull(results, "Results should not be null");
        // Empty directory should result in no imported documents
        assertEquals(0, results.size(),
            "Empty directory should produce no results");
    }

    /**
     * GROUP 3: Test progress listener receives updates
     */
    @Test
    void testProgressListenerUpdates() throws Exception {
        Path outputDir = tempDir.resolve("test_output");
        Files.createDirectories(outputDir);

        Options options = operation.getOptions(new AnnotatedPluginDocument[0]);
        options.setValue("outputPath", outputDir.toString());

        TestProgressListener progressListener = new TestProgressListener();

        // Perform operation
        operation.performOperation(new AnnotatedPluginDocument[0],
            progressListener, options);

        // Verify progress was updated
        assertTrue(progressListener.getProgressUpdateCount() > 0,
            "Progress listener should receive updates");
    }

    /**
     * GROUP 3: Test result type filtering options
     */
    @Test
    void testResultTypeFiltering() throws Exception {
        Options options = operation.getOptions(new AnnotatedPluginDocument[0]);

        // Disable all import options
        options.setValue("importReports", false);
        options.setValue("importTopHits", false);
        options.setValue("importFilteredReports", false);
        options.setValue("importAlignments", false);
        options.setValue("importQualityReports", false);

        // Verify values were set
        assertFalse((Boolean) options.getValue("importReports"));
        assertFalse((Boolean) options.getValue("importTopHits"));
        assertFalse((Boolean) options.getValue("importFilteredReports"));
        assertFalse((Boolean) options.getValue("importAlignments"));
        assertFalse((Boolean) options.getValue("importQualityReports"));
    }

    /**
     * GROUP 3: Test with directory containing test files
     */
    @Test
    void testImportWithTestFiles() throws Exception {
        // Create test directory structure
        Path outputDir = tempDir.resolve("test_output");
        Files.createDirectories(outputDir);

        // Create some test files
        Path reportsDir = outputDir.resolve("reports");
        Files.createDirectories(reportsDir);
        Files.writeString(reportsDir.resolve("test_report.txt"),
            "Test report content");

        Options options = operation.getOptions(new AnnotatedPluginDocument[0]);
        options.setValue("outputPath", outputDir.toString());

        TestProgressListener progressListener = new TestProgressListener();

        // Attempt import (may not work fully in unit test environment)
        assertDoesNotThrow(() -> {
            operation.performOperation(new AnnotatedPluginDocument[0],
                progressListener, options);
        }, "Operation should not throw exception with valid directory");
    }

    /**
     * GROUP 3: Test options validation for output path
     */
    @Test
    void testOutputPathValidation() throws Exception {
        Options options = operation.getOptions(new AnnotatedPluginDocument[0]);

        // Test setting various paths
        String testPath = tempDir.toString();
        options.setValue("outputPath", testPath);

        assertEquals(testPath, options.getValue("outputPath"),
            "Output path should be settable");
    }

    /**
     * GROUP 3: Test operation with null progress listener
     */
    @Test
    void testOperationWithNullProgressListener() throws Exception {
        Path outputDir = tempDir.resolve("test_output");
        Files.createDirectories(outputDir);

        Options options = operation.getOptions(new AnnotatedPluginDocument[0]);
        options.setValue("outputPath", outputDir.toString());

        // Should handle null progress listener gracefully
        assertDoesNotThrow(() -> {
            operation.performOperation(new AnnotatedPluginDocument[0],
                null, options);
        }, "Operation should handle null progress listener");
    }

    /**
     * GROUP 3: Test multiple result type combinations
     */
    @Test
    void testMultipleResultTypeCombinations() throws Exception {
        Options options = operation.getOptions(new AnnotatedPluginDocument[0]);

        // Test various combinations
        options.setValue("importReports", true);
        options.setValue("importTopHits", false);
        options.setValue("importFilteredReports", true);
        options.setValue("importAlignments", false);
        options.setValue("importQualityReports", true);

        assertTrue((Boolean) options.getValue("importReports"));
        assertFalse((Boolean) options.getValue("importTopHits"));
        assertTrue((Boolean) options.getValue("importFilteredReports"));
        assertFalse((Boolean) options.getValue("importAlignments"));
        assertTrue((Boolean) options.getValue("importQualityReports"));
    }

    // ===== Helper Classes =====

    /**
     * Test implementation of ProgressListener
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
    }
}
