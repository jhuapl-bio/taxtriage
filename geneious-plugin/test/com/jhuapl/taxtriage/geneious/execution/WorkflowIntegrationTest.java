package com.jhuapl.taxtriage.geneious.execution;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentOperationException;
import com.jhuapl.taxtriage.geneious.TaxTriageOptions;
import com.jhuapl.taxtriage.geneious.config.ConfigGenerator;
import com.jhuapl.taxtriage.geneious.config.SampleSheetBuilder;
import com.jhuapl.taxtriage.geneious.config.TaxTriageConfig;
import com.jhuapl.taxtriage.geneious.docker.DockerManager;
import com.jhuapl.taxtriage.geneious.docker.ExecutionResult;
import com.jhuapl.taxtriage.geneious.execution.WorkflowExecutor.SequenceExporter;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Path;
import java.util.concurrent.CompletableFuture;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.*;
import static org.mockito.Mockito.*;

/**
 * Integration tests for the complete TaxTriage workflow execution system.
 *
 * These tests verify that all components work together correctly to execute
 * a complete workflow from start to finish.
 */
public class WorkflowIntegrationTest {

    @TempDir
    private Path tempDir;

    private DockerManager mockDockerManager;
    private ConfigGenerator mockConfigGenerator;
    private SampleSheetBuilder mockSampleSheetBuilder;
    private SequenceExporter mockSequenceExporter;
    private WorkflowExecutor workflowExecutor;
    private TaxTriageWorkflowRunner workflowRunner;
    private MockProgressListener progressListener;

    @BeforeEach
    void setUp() throws Exception {
        // Create mock components
        mockDockerManager = mock(DockerManager.class);
        mockConfigGenerator = mock(ConfigGenerator.class);
        mockSampleSheetBuilder = mock(SampleSheetBuilder.class);
        mockSequenceExporter = mock(SequenceExporter.class);

        // Create real workflow components with mocked dependencies
        workflowExecutor = new WorkflowExecutor(
                mockDockerManager,
                mockConfigGenerator,
                mockSampleSheetBuilder,
                mockSequenceExporter
        );

        workflowRunner = new TaxTriageWorkflowRunner(workflowExecutor, tempDir);
        progressListener = new MockProgressListener();

        // Setup successful Docker execution by default
        ExecutionResult successResult = WorkflowTestUtils.createSuccessResult("Workflow completed successfully");
        when(mockDockerManager.executeNextflowCommand(
                anyString(), any(Path.class), any(Path.class), any(Path.class), any()))
                .thenReturn(successResult);
    }

    @Test
    void testCompleteWorkflowExecution() throws Exception {
        // Create test data
        AnnotatedPluginDocument[] documents = TestDataGenerator.createMockDocuments(3, 200);
        TaxTriageOptions options = TestDataGenerator.createMockOptions(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        // Execute complete workflow
        WorkflowContext result = workflowRunner.runWorkflow(documents, options, progressListener);

        // Verify workflow success
        WorkflowTestUtils.assertWorkflowSuccess(result);

        // Verify all phases were executed
        WorkflowTestUtils.assertAllPhasesExecuted(progressListener);
        WorkflowTestUtils.assertProgressUpdates(progressListener, 5);
        WorkflowTestUtils.assertCompleteProgress(progressListener);
        WorkflowTestUtils.assertMonotonicProgress(progressListener);

        // Verify component interactions
        verify(mockSequenceExporter, times(3)).exportToFastq(any(), any(File.class));
        verify(mockConfigGenerator).generateConfig(any(TaxTriageConfig.class), any(File.class));
        verify(mockConfigGenerator).generateParams(any(TaxTriageConfig.class), any(File.class));
        verify(mockSampleSheetBuilder).buildSampleSheet(anyList(), any(File.class), anyString());
        verify(mockDockerManager).executeNextflowCommand(
                anyString(), any(Path.class), any(Path.class), any(Path.class), any());

        // Verify file structure
        WorkflowTestUtils.assertDirectoriesCreated(result);
        WorkflowTestUtils.assertSequencesExported(result, 3);

        // Verify workflow metadata
        assertNotNull(result.getWorkflowId());
        assertNotNull(result.getStartedAt());
        assertNotNull(result.getCompletedAt());
        assertEquals(TaxTriageOptions.SequencingPreset.ILLUMINA_SE, result.getOptions().getSequencingPreset());
    }

    @Test
    void testWorkflowWithPairedEndReads() throws Exception {
        // Create test data for paired-end reads
        AnnotatedPluginDocument[] documents = TestDataGenerator.createMockDocuments(2, 150);
        TaxTriageOptions options = TestDataGenerator.createMockOptions(TaxTriageOptions.SequencingPreset.ILLUMINA_PE);

        // Execute workflow
        WorkflowContext result = workflowRunner.runWorkflow(documents, options, progressListener);

        // Verify success
        WorkflowTestUtils.assertWorkflowSuccess(result);
        assertEquals(TaxTriageOptions.SequencingPreset.ILLUMINA_PE, result.getOptions().getSequencingPreset());

        // Verify samplesheet was built with correct preset
        verify(mockSampleSheetBuilder).buildSampleSheet(anyList(), any(File.class), eq("ILLUMINA_PE"));
    }

    @Test
    void testWorkflowWithOxfordNanoporeReads() throws Exception {
        // Create test data for Oxford Nanopore reads
        AnnotatedPluginDocument[] documents = TestDataGenerator.createMockDocuments(1, 5000);
        TaxTriageOptions options = TestDataGenerator.createMockOptions(TaxTriageOptions.SequencingPreset.ONT);

        // Execute workflow
        WorkflowContext result = workflowRunner.runWorkflow(documents, options, progressListener);

        // Verify success
        WorkflowTestUtils.assertWorkflowSuccess(result);
        assertEquals(TaxTriageOptions.SequencingPreset.ONT, result.getOptions().getSequencingPreset());

        // Verify samplesheet was built with correct preset
        verify(mockSampleSheetBuilder).buildSampleSheet(anyList(), any(File.class), eq("ONT"));
    }

    @Test
    void testWorkflowFailureDuringSequenceExport() throws Exception {
        // Setup sequence export to fail
        doThrow(new RuntimeException("Export failed")).when(mockSequenceExporter)
                .exportToFastq(any(), any(File.class));

        AnnotatedPluginDocument[] documents = TestDataGenerator.createMockDocuments(1, 100);
        TaxTriageOptions options = TestDataGenerator.createMockOptions(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        // Execute workflow
        WorkflowContext result = workflowRunner.runWorkflow(documents, options, progressListener);

        // Verify failure
        WorkflowTestUtils.assertWorkflowFailure(result);
        assertTrue(result.getLastError().getMessage().contains("Export failed"));

        // Verify that subsequent phases were not executed
        verify(mockConfigGenerator, never()).generateConfig(any(), any());
        verify(mockSampleSheetBuilder, never()).buildSampleSheet(any(), any(), any());
        verify(mockDockerManager, never()).executeNextflowCommand(any(), any(), any(), any(), any());
    }

    @Test
    void testWorkflowFailureDuringConfigGeneration() throws Exception {
        // Setup config generation to fail
        doThrow(new RuntimeException("Config generation failed")).when(mockConfigGenerator)
                .generateConfig(any(TaxTriageConfig.class), any(File.class));

        AnnotatedPluginDocument[] documents = TestDataGenerator.createMockDocuments(1, 100);
        TaxTriageOptions options = TestDataGenerator.createMockOptions(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        // Execute workflow
        WorkflowContext result = workflowRunner.runWorkflow(documents, options, progressListener);

        // Verify failure
        WorkflowTestUtils.assertWorkflowFailure(result);
        assertTrue(result.getLastError().getMessage().contains("Config generation failed"));

        // Verify that sequence export succeeded but subsequent phases failed
        verify(mockSequenceExporter).exportToFastq(any(), any(File.class));
        verify(mockSampleSheetBuilder, never()).buildSampleSheet(any(), any(), any());
        verify(mockDockerManager, never()).executeNextflowCommand(any(), any(), any(), any(), any());
    }

    @Test
    void testWorkflowFailureDuringSampleSheetCreation() throws Exception {
        // Setup samplesheet creation to fail
        doThrow(new RuntimeException("Samplesheet creation failed")).when(mockSampleSheetBuilder)
                .buildSampleSheet(any(), any(), any());

        AnnotatedPluginDocument[] documents = TestDataGenerator.createMockDocuments(1, 100);
        TaxTriageOptions options = TestDataGenerator.createMockOptions(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        // Execute workflow
        WorkflowContext result = workflowRunner.runWorkflow(documents, options, progressListener);

        // Verify failure
        WorkflowTestUtils.assertWorkflowFailure(result);
        assertTrue(result.getLastError().getMessage().contains("Samplesheet creation failed"));

        // Verify that export and config succeeded but execution didn't happen
        verify(mockSequenceExporter).exportToFastq(any(), any(File.class));
        verify(mockConfigGenerator).generateConfig(any(TaxTriageConfig.class), any(File.class));
        verify(mockDockerManager, never()).executeNextflowCommand(any(), any(), any(), any(), any());
    }

    @Test
    void testWorkflowFailureDuringDockerExecution() throws Exception {
        // Setup Docker execution to fail
        ExecutionResult failureResult = WorkflowTestUtils.createFailureResult("Nextflow execution failed");
        when(mockDockerManager.executeNextflowCommand(
                anyString(), any(Path.class), any(Path.class), any(Path.class), any()))
                .thenReturn(failureResult);

        AnnotatedPluginDocument[] documents = TestDataGenerator.createMockDocuments(1, 100);
        TaxTriageOptions options = TestDataGenerator.createMockOptions(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        // Execute workflow
        WorkflowContext result = workflowRunner.runWorkflow(documents, options, progressListener);

        // Verify failure
        WorkflowTestUtils.assertWorkflowFailure(result);
        assertTrue(result.getLastError().getMessage().contains("Nextflow execution failed"));

        // Verify that all preparation phases succeeded
        verify(mockSequenceExporter).exportToFastq(any(), any(File.class));
        verify(mockConfigGenerator).generateConfig(any(TaxTriageConfig.class), any(File.class));
        verify(mockSampleSheetBuilder).buildSampleSheet(any(), any(), any());
        verify(mockDockerManager).executeNextflowCommand(any(), any(), any(), any(), any());
    }

    @Test
    void testAsyncWorkflowExecution() throws Exception {
        AnnotatedPluginDocument[] documents = TestDataGenerator.createMockDocuments(2, 150);
        TaxTriageOptions options = TestDataGenerator.createMockOptions(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        // Execute workflow asynchronously
        CompletableFuture<WorkflowContext> future = workflowRunner.runWorkflowAsync(
                documents, options, progressListener);

        WorkflowContext result = WorkflowTestUtils.waitForCompletion(future, 30);

        // Verify success
        WorkflowTestUtils.assertWorkflowSuccess(result);
        WorkflowTestUtils.assertProgressUpdates(progressListener, 5);
    }

    @Test
    void testWorkflowCancellation() throws Exception {
        // Simulate user cancellation
        progressListener.cancel();

        AnnotatedPluginDocument[] documents = TestDataGenerator.createMockDocuments(1, 100);
        TaxTriageOptions options = TestDataGenerator.createMockOptions(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        // Execute workflow - should detect cancellation
        assertThrows(DocumentOperationException.class, () -> {
            workflowRunner.runWorkflow(documents, options, progressListener);
        });

        assertTrue(progressListener.isCancelled());
    }

    @Test
    void testWorkflowWithLargeDataset() throws Exception {
        // Test with larger dataset
        AnnotatedPluginDocument[] documents = TestDataGenerator.createMockDocuments(10, 300);
        TaxTriageOptions options = TestDataGenerator.createMockOptions(TaxTriageOptions.SequencingPreset.ILLUMINA_PE);

        // Execute workflow
        WorkflowContext result = workflowRunner.runWorkflow(documents, options, progressListener);

        // Verify success
        WorkflowTestUtils.assertWorkflowSuccess(result);
        WorkflowTestUtils.assertSequencesExported(result, 10);

        // Verify all sequences were exported
        verify(mockSequenceExporter, times(10)).exportToFastq(any(), any(File.class));

        // Verify execution time estimate is reasonable
        int estimatedTime = workflowRunner.estimateExecutionTime(documents, options);
        assertTrue(estimatedTime > 20); // Should be higher for more documents
        assertTrue(estimatedTime <= 240); // But capped at maximum
    }

    @Test
    void testWorkflowProgressTracking() throws Exception {
        AnnotatedPluginDocument[] documents = TestDataGenerator.createMockDocuments(2, 150);
        TaxTriageOptions options = TestDataGenerator.createMockOptions(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        // Execute workflow
        WorkflowContext result = workflowRunner.runWorkflow(documents, options, progressListener);

        // Verify success
        WorkflowTestUtils.assertWorkflowSuccess(result);

        // Verify progress milestones were reached
        WorkflowTestUtils.assertProgressMilestones(progressListener,
                0.0,    // Initial
                0.15,   // After export
                0.20,   // After config
                0.25,   // After samplesheet
                0.95,   // After execution
                1.0     // Complete
        );

        // Verify specific messages were sent
        WorkflowTestUtils.assertMessagesContain(progressListener,
                "Starting TaxTriage workflow",
                "Exporting sequences",
                "Generating workflow configuration",
                "Creating samplesheet",
                "Executing TaxTriage workflow",
                "Analysis complete"
        );
    }

    @Test
    void testSystemRequirements() {
        TaxTriageWorkflowRunner.SystemStatus status = workflowRunner.checkSystemRequirements();

        assertNotNull(status);
        assertNotNull(status.getSummary());

        // Basic sanity checks
        assertTrue(status.availableMemoryMB > 0);
        assertTrue(status.availableDiskSpaceGB >= 0);

        // Summary should contain key information
        String summary = status.getSummary();
        assertTrue(summary.contains("System Status"));
        assertTrue(summary.contains("Docker"));
        assertTrue(summary.contains("Disk Space"));
        assertTrue(summary.contains("Memory"));
    }

    @Test
    void testWorkflowDirectoryCleanup() throws Exception {
        AnnotatedPluginDocument[] documents = TestDataGenerator.createMockDocuments(1, 100);
        TaxTriageOptions options = TestDataGenerator.createMockOptions(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        // Execute workflow
        WorkflowContext result = workflowRunner.runWorkflow(documents, options, progressListener);

        // Verify success and directories were created
        WorkflowTestUtils.assertWorkflowSuccess(result);
        WorkflowTestUtils.assertDirectoriesCreated(result);

        // Note: In a real implementation, you might want to test cleanup,
        // but for now we just verify the directories were created properly
        assertTrue(result.getWorkingDirectory().toFile().exists());
    }

    @Test
    void testErrorRecoveryAndReporting() throws Exception {
        // Setup a specific failure scenario
        doThrow(new RuntimeException("Disk space full")).when(mockSequenceExporter)
                .exportToFastq(any(), any(File.class));

        AnnotatedPluginDocument[] documents = TestDataGenerator.createMockDocuments(1, 100);
        TaxTriageOptions options = TestDataGenerator.createMockOptions(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        // Execute workflow
        WorkflowContext result = workflowRunner.runWorkflow(documents, options, progressListener);

        // Verify failure details are captured
        WorkflowTestUtils.assertWorkflowFailure(result);
        assertNotNull(result.getLastError());
        assertTrue(result.getLastError().getMessage().contains("Disk space full"));

        // Verify error state is properly set
        assertEquals(WorkflowContext.ExecutionState.FAILED, result.getState());
        assertNotNull(result.getCompletedAt());

        // Verify progress listener received error message
        String finalMessage = progressListener.getFinalMessage();
        assertTrue(finalMessage.contains("failed") || finalMessage.contains("error"));
    }
}