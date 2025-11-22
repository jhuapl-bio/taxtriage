package com.jhuapl.taxtriage.geneious.execution;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.jhuapl.taxtriage.geneious.TaxTriageOptions;
import com.jhuapl.taxtriage.geneious.config.TaxTriageConfig;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.mock;

/**
 * Unit tests for WorkflowContext class.
 */
public class WorkflowContextTest {

    @TempDir
    private Path tempDir;

    private List<AnnotatedPluginDocument> mockDocuments;
    private TaxTriageOptions mockOptions;
    private TaxTriageConfig mockConfig;
    private WorkflowContext context;

    @BeforeEach
    void setUp() {
        mockDocuments = Arrays.asList(mock(AnnotatedPluginDocument.class));
        mockOptions = mock(TaxTriageOptions.class);
        mockConfig = mock(TaxTriageConfig.class);
        context = new WorkflowContext(mockDocuments, mockOptions, mockConfig, tempDir);
    }

    @Test
    void testInitialization() {
        assertNotNull(context.getWorkflowId());
        assertNotNull(context.getCreatedAt());
        assertNull(context.getStartedAt());
        assertNull(context.getCompletedAt());
        assertEquals(WorkflowContext.ExecutionState.INITIALIZED, context.getState());
        assertEquals(0.0, context.getProgress());
        assertEquals("Workflow initialized", context.getCurrentMessage());
        assertSame(mockDocuments, context.getInputDocuments());
        assertSame(mockOptions, context.getOptions());
        assertSame(mockConfig, context.getConfig());
    }

    @Test
    void testDirectoryStructure() {
        assertEquals(tempDir, context.getWorkingDirectory());
        assertEquals(tempDir.resolve("input"), context.getInputDirectory());
        assertEquals(tempDir.resolve("output"), context.getOutputDirectory());
        assertEquals(tempDir.resolve("config"), context.getConfigDirectory());
    }

    @Test
    void testStateTransitions() {
        assertFalse(context.isTerminal());
        assertFalse(context.isSuccessful());

        // Test state update
        context.updateState(WorkflowContext.ExecutionState.EXPORTING_SEQUENCES, "Exporting...", 0.1);
        assertEquals(WorkflowContext.ExecutionState.EXPORTING_SEQUENCES, context.getState());
        assertEquals("Exporting...", context.getCurrentMessage());
        assertEquals(0.1, context.getProgress());

        // Test completion
        context.markCompleted();
        assertTrue(context.isTerminal());
        assertTrue(context.isSuccessful());
        assertEquals(WorkflowContext.ExecutionState.COMPLETED, context.getState());
        assertEquals(1.0, context.getProgress());
        assertNotNull(context.getCompletedAt());
    }

    @Test
    void testFailureHandling() {
        Exception testError = new RuntimeException("Test error");
        context.markFailed(testError, "Workflow failed");

        assertTrue(context.isTerminal());
        assertFalse(context.isSuccessful());
        assertEquals(WorkflowContext.ExecutionState.FAILED, context.getState());
        assertEquals("Workflow failed", context.getCurrentMessage());
        assertSame(testError, context.getLastError());
        assertNotNull(context.getCompletedAt());
    }

    @Test
    void testCancellation() {
        context.markCancelled();

        assertTrue(context.isTerminal());
        assertFalse(context.isSuccessful());
        assertEquals(WorkflowContext.ExecutionState.CANCELLED, context.getState());
        assertEquals("Workflow cancelled", context.getCurrentMessage());
        assertNotNull(context.getCompletedAt());
    }

    @Test
    void testProgressValidation() {
        // Test valid progress values
        context.setProgress(0.5);
        assertEquals(0.5, context.getProgress());

        // Test clamping to valid range
        context.setProgress(-0.1);
        assertEquals(0.0, context.getProgress());

        context.setProgress(1.5);
        assertEquals(1.0, context.getProgress());
    }

    @Test
    void testTimestamps() {
        LocalDateTime beforeStart = LocalDateTime.now();
        context.setStartedAt(LocalDateTime.now());
        LocalDateTime afterStart = LocalDateTime.now();

        assertTrue(context.getStartedAt().isAfter(beforeStart) ||
                  context.getStartedAt().isEqual(beforeStart));
        assertTrue(context.getStartedAt().isBefore(afterStart) ||
                  context.getStartedAt().isEqual(afterStart));

        LocalDateTime beforeComplete = LocalDateTime.now();
        context.setCompletedAt(LocalDateTime.now());
        LocalDateTime afterComplete = LocalDateTime.now();

        assertTrue(context.getCompletedAt().isAfter(beforeComplete) ||
                  context.getCompletedAt().isEqual(beforeComplete));
        assertTrue(context.getCompletedAt().isBefore(afterComplete) ||
                  context.getCompletedAt().isEqual(afterComplete));
    }

    @Test
    void testFileSetters() {
        File sampleSheet = new File("samplesheet.csv");
        File configFile = new File("nextflow.config");
        File paramsFile = new File("params.json");

        context.setSampleSheetFile(sampleSheet);
        context.setConfigFile(configFile);
        context.setParamsFile(paramsFile);

        assertSame(sampleSheet, context.getSampleSheetFile());
        assertSame(configFile, context.getConfigFile());
        assertSame(paramsFile, context.getParamsFile());
    }

    @Test
    void testResultFiles() {
        File workflowReport = new File("workflow_report.html");
        File taxonomyResults = new File("taxonomy.csv");
        File qualityReport = new File("quality_report.html");

        context.setWorkflowReport(workflowReport);
        context.setTaxonomyResults(taxonomyResults);
        context.setQualityReport(qualityReport);

        assertSame(workflowReport, context.getWorkflowReport());
        assertSame(taxonomyResults, context.getTaxonomyResults());
        assertSame(qualityReport, context.getQualityReport());
    }

    @Test
    void testSummary() {
        context.updateState(WorkflowContext.ExecutionState.EXECUTING_WORKFLOW, "Running...", 0.5);
        context.setStartedAt(LocalDateTime.now());

        String summary = context.getSummary();
        assertNotNull(summary);
        assertTrue(summary.contains(context.getWorkflowId()));
        assertTrue(summary.contains("EXECUTING_WORKFLOW"));
        assertTrue(summary.contains("50.0%"));
        assertTrue(summary.contains("Running..."));
    }

    @Test
    void testSummaryWithError() {
        Exception error = new RuntimeException("Test error");
        context.markFailed(error, "Failed");

        String summary = context.getSummary();
        assertTrue(summary.contains("Test error"));
        assertTrue(summary.contains("FAILED"));
    }

    @Test
    void testUniqueWorkflowIds() {
        WorkflowContext context2 = new WorkflowContext(mockDocuments, mockOptions, mockConfig, tempDir);
        assertNotEquals(context.getWorkflowId(), context2.getWorkflowId());
    }
}