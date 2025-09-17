package com.jhuapl.taxtriage.geneious.execution;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.sequence.SequenceDocument;
import com.jhuapl.taxtriage.geneious.TaxTriageOptions;
import com.jhuapl.taxtriage.geneious.config.ConfigGenerator;
import com.jhuapl.taxtriage.geneious.config.SampleSheetBuilder;
import com.jhuapl.taxtriage.geneious.config.TaxTriageConfig;
import com.jhuapl.taxtriage.geneious.docker.DockerManager;
import com.jhuapl.taxtriage.geneious.docker.ExecutionResult;
import com.jhuapl.taxtriage.geneious.execution.WorkflowExecutor.SequenceExporter;
import jebl.util.ProgressListener;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.concurrent.CompletableFuture;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.*;
import static org.mockito.Mockito.*;

/**
 * Unit tests for WorkflowExecutor class.
 */
public class WorkflowExecutorTest {

    @TempDir
    private Path tempDir;

    private DockerManager mockDockerManager;
    private ConfigGenerator mockConfigGenerator;
    private SampleSheetBuilder mockSampleSheetBuilder;
    private SequenceExporter mockSequenceExporter;
    private WorkflowExecutor workflowExecutor;
    private WorkflowContext mockContext;
    private ProgressListener mockProgressListener;

    @BeforeEach
    void setUp() throws Exception {
        mockDockerManager = mock(DockerManager.class);
        mockConfigGenerator = mock(ConfigGenerator.class);
        mockSampleSheetBuilder = mock(SampleSheetBuilder.class);
        mockSequenceExporter = mock(SequenceExporter.class);
        mockProgressListener = mock(ProgressListener.class);

        workflowExecutor = new WorkflowExecutor(
                mockDockerManager,
                mockConfigGenerator,
                mockSampleSheetBuilder,
                mockSequenceExporter
        );

        // Setup mock context
        mockContext = createMockContext();
    }

    @Test
    void testSuccessfulWorkflowExecution() throws Exception {
        // Setup successful execution
        ExecutionResult successResult = new ExecutionResult(0, "Success", "");
        when(mockDockerManager.executeNextflowCommand(
                anyString(), any(Path.class), any(Path.class), any(Path.class), any(ProgressListener.class)))
                .thenReturn(successResult);

        // Execute workflow
        CompletableFuture<WorkflowContext> future = workflowExecutor.executeWorkflow(mockContext, mockProgressListener);
        WorkflowContext result = future.get();

        // Verify success
        assertTrue(result.isSuccessful());
        assertEquals(WorkflowContext.ExecutionState.COMPLETED, result.getState());
        assertNotNull(result.getCompletedAt());

        // Verify interactions
        verify(mockConfigGenerator).generateConfig(any(TaxTriageConfig.class), any(File.class));
        verify(mockConfigGenerator).generateParams(any(TaxTriageConfig.class), any(File.class));
        verify(mockSampleSheetBuilder).buildSampleSheet(anyList(), any(File.class), anyString());
        verify(mockDockerManager).executeNextflowCommand(
                anyString(), any(Path.class), any(Path.class), any(Path.class), any(ProgressListener.class));
    }

    @Test
    void testFailedWorkflowExecution() throws Exception {
        // Setup failed execution
        ExecutionResult failureResult = new ExecutionResult(1, "Output", "Error occurred");
        when(mockDockerManager.executeNextflowCommand(
                anyString(), any(Path.class), any(Path.class), any(Path.class), any(ProgressListener.class)))
                .thenReturn(failureResult);

        // Execute workflow
        CompletableFuture<WorkflowContext> future = workflowExecutor.executeWorkflow(mockContext, mockProgressListener);
        WorkflowContext result = future.get();

        // Verify failure
        assertFalse(result.isSuccessful());
        assertEquals(WorkflowContext.ExecutionState.FAILED, result.getState());
        assertNotNull(result.getLastError());
        assertNotNull(result.getCompletedAt());
    }

    @Test
    void testExceptionDuringExecution() throws Exception {
        // Setup exception during Docker execution
        when(mockDockerManager.executeNextflowCommand(
                anyString(), any(Path.class), any(Path.class), any(Path.class), any(ProgressListener.class)))
                .thenThrow(new RuntimeException("Docker execution failed"));

        // Execute workflow
        CompletableFuture<WorkflowContext> future = workflowExecutor.executeWorkflow(mockContext, mockProgressListener);
        WorkflowContext result = future.get();

        // Verify failure handling
        assertFalse(result.isSuccessful());
        assertEquals(WorkflowContext.ExecutionState.FAILED, result.getState());
        assertNotNull(result.getLastError());
        assertTrue(result.getLastError().getMessage().contains("Docker execution failed"));
    }

    @Test
    void testSequenceExportPhase() throws Exception {
        // Create real documents for export testing
        AnnotatedPluginDocument mockDoc = mock(AnnotatedPluginDocument.class);
        SequenceDocument mockSeqDoc = mock(SequenceDocument.class);
        when(mockDoc.getDocument()).thenReturn(mockSeqDoc);
        when(mockDoc.getName()).thenReturn("test_sequence");
        when(mockSeqDoc.getSequenceString()).thenReturn("ATCGATCGATCG");
        when(mockSeqDoc.getSequenceLength()).thenReturn(12);

        TaxTriageOptions mockOptions = mock(TaxTriageOptions.class);
        TaxTriageConfig mockConfig = mock(TaxTriageConfig.class);

        WorkflowContext realContext = new WorkflowContext(
                Arrays.asList(mockDoc), mockOptions, mockConfig, tempDir);

        // Setup successful subsequent phases
        ExecutionResult successResult = new ExecutionResult(0, "Success", "");
        when(mockDockerManager.executeNextflowCommand(
                anyString(), any(Path.class), any(Path.class), any(Path.class), any(ProgressListener.class)))
                .thenReturn(successResult);

        // Execute workflow
        CompletableFuture<WorkflowContext> future = workflowExecutor.executeWorkflow(realContext, mockProgressListener);
        WorkflowContext result = future.get();

        // Verify export occurred
        assertNotNull(result.getExportedFiles());
        assertEquals(1, result.getExportedFiles().size());
        verify(mockSequenceExporter).exportToFastq(eq(mockSeqDoc), any(File.class));
    }

    @Test
    void testProgressUpdates() throws Exception {
        // Setup successful execution
        ExecutionResult successResult = new ExecutionResult(0, "Success", "");
        when(mockDockerManager.executeNextflowCommand(
                anyString(), any(Path.class), any(Path.class), any(Path.class), any(ProgressListener.class)))
                .thenReturn(successResult);

        // Execute workflow
        CompletableFuture<WorkflowContext> future = workflowExecutor.executeWorkflow(mockContext, mockProgressListener);
        WorkflowContext result = future.get();

        // Verify progress updates were made
        verify(mockProgressListener, atLeastOnce()).setMessage(anyString());
        verify(mockProgressListener, atLeastOnce()).setProgress(anyDouble());

        // Verify final progress
        assertEquals(1.0, result.getProgress());
    }

    @Test
    void testWorkflowCancellation() {
        workflowExecutor.cancelWorkflow(mockContext);
        assertEquals(WorkflowContext.ExecutionState.CANCELLED, mockContext.getState());
    }

    @Test
    void testSequenceExporter() throws Exception {
        SequenceExporter exporter = new SequenceExporter();
        SequenceDocument mockSeqDoc = mock(SequenceDocument.class);
        when(mockSeqDoc.getName()).thenReturn("test_sequence");
        when(mockSeqDoc.getSequenceString()).thenReturn("ATCGATCGATCG");
        when(mockSeqDoc.getSequenceLength()).thenReturn(12);

        File outputFile = tempDir.resolve("test_output.fastq").toFile();
        exporter.exportToFastq(mockSeqDoc, outputFile);

        assertTrue(outputFile.exists());
        String content = Files.readString(outputFile.toPath());
        assertTrue(content.contains("@test_sequence"));
        assertTrue(content.contains("ATCGATCGATCG"));
        assertTrue(content.contains("+"));
        // Should contain quality string
        assertTrue(content.split("\n").length >= 4);
    }

    @Test
    void testSequenceExporterCreatesDirectory() throws Exception {
        SequenceExporter exporter = new SequenceExporter();
        SequenceDocument mockSeqDoc = mock(SequenceDocument.class);
        when(mockSeqDoc.getName()).thenReturn("test_sequence");
        when(mockSeqDoc.getSequenceString()).thenReturn("ATCG");
        when(mockSeqDoc.getSequenceLength()).thenReturn(4);

        // Use nested directory that doesn't exist
        File outputFile = tempDir.resolve("nested/dir/test_output.fastq").toFile();
        exporter.exportToFastq(mockSeqDoc, outputFile);

        assertTrue(outputFile.exists());
        assertTrue(outputFile.getParentFile().exists());
    }

    @Test
    void testWorkflowWithNullProgressListener() throws Exception {
        ExecutionResult successResult = new ExecutionResult(0, "Success", "");
        when(mockDockerManager.executeNextflowCommand(
                anyString(), any(Path.class), any(Path.class), any(Path.class), any()))
                .thenReturn(successResult);

        // Execute workflow without progress listener
        CompletableFuture<WorkflowContext> future = workflowExecutor.executeWorkflow(mockContext, null);
        WorkflowContext result = future.get();

        assertTrue(result.isSuccessful());
        assertEquals(WorkflowContext.ExecutionState.COMPLETED, result.getState());
    }

    private WorkflowContext createMockContext() throws Exception {
        AnnotatedPluginDocument mockDoc = mock(AnnotatedPluginDocument.class);
        SequenceDocument mockSeqDoc = mock(SequenceDocument.class);
        when(mockDoc.getDocument()).thenReturn(mockSeqDoc);
        when(mockDoc.getName()).thenReturn("test_doc");
        when(mockSeqDoc.getSequenceString()).thenReturn("ATCG");
        when(mockSeqDoc.getSequenceLength()).thenReturn(4);

        TaxTriageOptions mockOptions = mock(TaxTriageOptions.class);
        when(mockOptions.getSequencingPreset()).thenReturn(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        TaxTriageConfig mockConfig = mock(TaxTriageConfig.class);
        when(mockConfig.validate()).thenReturn(null); // Valid configuration

        WorkflowContext context = new WorkflowContext(
                Arrays.asList(mockDoc), mockOptions, mockConfig, tempDir);

        // Pre-create directories
        Files.createDirectories(context.getInputDirectory());
        Files.createDirectories(context.getOutputDirectory());
        Files.createDirectories(context.getConfigDirectory());

        // Set required files
        File sampleSheet = context.getConfigDirectory().resolve("samplesheet.csv").toFile();
        File configFile = context.getConfigDirectory().resolve("nextflow.config").toFile();
        File paramsFile = context.getConfigDirectory().resolve("params.json").toFile();

        Files.createFile(sampleSheet.toPath());
        Files.createFile(configFile.toPath());
        Files.createFile(paramsFile.toPath());

        context.setSampleSheetFile(sampleSheet);
        context.setConfigFile(configFile);
        context.setParamsFile(paramsFile);
        context.setExportedFiles(Arrays.asList(new File("test.fastq")));

        return context;
    }
}