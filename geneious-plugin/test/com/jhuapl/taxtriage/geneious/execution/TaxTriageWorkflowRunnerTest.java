package com.jhuapl.taxtriage.geneious.execution;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.sequence.SequenceDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentOperationException;
import com.jhuapl.taxtriage.geneious.TaxTriageOptions;
import jebl.util.ProgressListener;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.nio.file.Path;
import java.time.LocalDateTime;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.*;
import static org.mockito.Mockito.*;

/**
 * Unit tests for TaxTriageWorkflowRunner class.
 */
public class TaxTriageWorkflowRunnerTest {

    @TempDir
    private Path tempDir;

    private WorkflowExecutor mockWorkflowExecutor;
    private TaxTriageWorkflowRunner workflowRunner;
    private AnnotatedPluginDocument[] mockDocuments;
    private TaxTriageOptions mockOptions;
    private ProgressListener mockProgressListener;

    @BeforeEach
    void setUp() {
        mockWorkflowExecutor = mock(WorkflowExecutor.class);
        workflowRunner = new TaxTriageWorkflowRunner(mockWorkflowExecutor, tempDir);

        // Setup mock documents
        AnnotatedPluginDocument mockDoc = mock(AnnotatedPluginDocument.class);
        SequenceDocument mockSeqDoc = mock(SequenceDocument.class);
        when(mockDoc.getDocument()).thenReturn(mockSeqDoc);
        when(mockDoc.getName()).thenReturn("test_sequence");
        mockDocuments = new AnnotatedPluginDocument[]{mockDoc};

        mockOptions = mock(TaxTriageOptions.class);
        when(mockOptions.getSequencingPreset()).thenReturn(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        mockProgressListener = mock(ProgressListener.class);
    }

    @Test
    void testSuccessfulWorkflowExecution() throws Exception {
        // Setup successful workflow context
        WorkflowContext successContext = createSuccessfulContext();
        CompletableFuture<WorkflowContext> successFuture = CompletableFuture.completedFuture(successContext);
        when(mockWorkflowExecutor.executeWorkflow(any(WorkflowContext.class), any(ProgressListener.class)))
                .thenReturn(successFuture);

        // Execute workflow
        WorkflowContext result = workflowRunner.runWorkflow(mockDocuments, mockOptions, mockProgressListener);

        // Verify success
        assertNotNull(result);
        assertTrue(result.isSuccessful());
        assertEquals(WorkflowContext.ExecutionState.COMPLETED, result.getState());

        // Verify progress updates
        verify(mockProgressListener, atLeastOnce()).setMessage(anyString());
        verify(mockProgressListener, atLeastOnce()).setProgress(anyDouble());
    }

    @Test
    void testFailedWorkflowExecution() throws Exception {
        // Setup failed workflow context
        WorkflowContext failedContext = createFailedContext();
        CompletableFuture<WorkflowContext> failedFuture = CompletableFuture.completedFuture(failedContext);
        when(mockWorkflowExecutor.executeWorkflow(any(WorkflowContext.class), any(ProgressListener.class)))
                .thenReturn(failedFuture);

        // Execute workflow and expect exception
        assertThrows(DocumentOperationException.class, () -> {
            workflowRunner.runWorkflow(mockDocuments, mockOptions, mockProgressListener);
        });
    }

    @Test
    void testWorkflowTimeout() throws Exception {
        // Setup never-completing future
        CompletableFuture<WorkflowContext> neverCompletingFuture = new CompletableFuture<>();
        when(mockWorkflowExecutor.executeWorkflow(any(WorkflowContext.class), any(ProgressListener.class)))
                .thenReturn(neverCompletingFuture);

        // Execute workflow with very short timeout
        assertThrows(DocumentOperationException.class, () -> {
            workflowRunner.runWorkflow(mockDocuments, mockOptions, mockProgressListener, 0); // 0 minute timeout
        });

        // Verify cancellation was attempted
        verify(mockWorkflowExecutor).cancelWorkflow(any(WorkflowContext.class));
    }

    @Test
    void testWorkflowInterruption() throws Exception {
        // Setup future that throws InterruptedException
        CompletableFuture<WorkflowContext> interruptedFuture = new CompletableFuture<>();
        interruptedFuture.completeExceptionally(new InterruptedException("Thread interrupted"));
        when(mockWorkflowExecutor.executeWorkflow(any(WorkflowContext.class), any(ProgressListener.class)))
                .thenReturn(interruptedFuture);

        // Execute workflow and expect exception
        assertThrows(DocumentOperationException.class, () -> {
            workflowRunner.runWorkflow(mockDocuments, mockOptions, mockProgressListener);
        });

        // Verify cancellation was attempted
        verify(mockWorkflowExecutor).cancelWorkflow(any(WorkflowContext.class));
    }

    @Test
    void testUserCancellation() throws Exception {
        // Setup progress listener that reports cancellation
        when(mockProgressListener.isCancelled()).thenReturn(true);

        WorkflowContext successContext = createSuccessfulContext();
        CompletableFuture<WorkflowContext> successFuture = CompletableFuture.completedFuture(successContext);
        when(mockWorkflowExecutor.executeWorkflow(any(WorkflowContext.class), any(ProgressListener.class)))
                .thenReturn(successFuture);

        // Execute workflow and expect cancellation exception
        assertThrows(DocumentOperationException.class, () -> {
            workflowRunner.runWorkflow(mockDocuments, mockOptions, mockProgressListener);
        });

        // Verify cancellation was attempted
        verify(mockWorkflowExecutor).cancelWorkflow(any(WorkflowContext.class));
    }

    @Test
    void testAsyncWorkflowExecution() throws Exception {
        WorkflowContext successContext = createSuccessfulContext();
        CompletableFuture<WorkflowContext> successFuture = CompletableFuture.completedFuture(successContext);
        when(mockWorkflowExecutor.executeWorkflow(any(WorkflowContext.class), any(ProgressListener.class)))
                .thenReturn(successFuture);

        // Execute workflow asynchronously
        CompletableFuture<WorkflowContext> asyncResult = workflowRunner.runWorkflowAsync(
                mockDocuments, mockOptions, mockProgressListener);

        WorkflowContext result = asyncResult.get();
        assertNotNull(result);
        assertTrue(result.isSuccessful());
    }

    @Test
    void testAsyncWorkflowFailure() throws Exception {
        WorkflowContext failedContext = createFailedContext();
        CompletableFuture<WorkflowContext> failedFuture = CompletableFuture.completedFuture(failedContext);
        when(mockWorkflowExecutor.executeWorkflow(any(WorkflowContext.class), any(ProgressListener.class)))
                .thenReturn(failedFuture);

        // Execute workflow asynchronously
        CompletableFuture<WorkflowContext> asyncResult = workflowRunner.runWorkflowAsync(
                mockDocuments, mockOptions, mockProgressListener);

        // Should complete with RuntimeException wrapping DocumentOperationException
        assertThrows(ExecutionException.class, () -> {
            asyncResult.get();
        });
    }

    @Test
    void testInputValidation() {
        // Test null documents
        assertThrows(DocumentOperationException.class, () -> {
            workflowRunner.runWorkflow(null, mockOptions, mockProgressListener);
        });

        // Test empty documents
        assertThrows(DocumentOperationException.class, () -> {
            workflowRunner.runWorkflow(new AnnotatedPluginDocument[0], mockOptions, mockProgressListener);
        });

        // Test null options
        assertThrows(DocumentOperationException.class, () -> {
            workflowRunner.runWorkflow(mockDocuments, null, mockProgressListener);
        });

        // Test document with null content
        AnnotatedPluginDocument nullContentDoc = mock(AnnotatedPluginDocument.class);
        when(nullContentDoc.getDocument()).thenReturn(null);
        when(nullContentDoc.getName()).thenReturn("null_doc");
        AnnotatedPluginDocument[] invalidDocs = new AnnotatedPluginDocument[]{nullContentDoc};

        assertThrows(DocumentOperationException.class, () -> {
            workflowRunner.runWorkflow(invalidDocs, mockOptions, mockProgressListener);
        });
    }

    @Test
    void testExecutionTimeEstimation() {
        // Test with no documents
        int time1 = workflowRunner.estimateExecutionTime(new AnnotatedPluginDocument[0], mockOptions);
        assertEquals(5, time1); // Minimum time

        // Test with single document
        int time2 = workflowRunner.estimateExecutionTime(mockDocuments, mockOptions);
        assertTrue(time2 > 5);

        // Test with multiple documents
        AnnotatedPluginDocument[] multiDocs = new AnnotatedPluginDocument[5];
        for (int i = 0; i < 5; i++) {
            multiDocs[i] = mock(AnnotatedPluginDocument.class);
        }
        int time3 = workflowRunner.estimateExecutionTime(multiDocs, mockOptions);
        assertTrue(time3 > time2);

        // Test with ONT preset (should take longer)
        when(mockOptions.getSequencingPreset()).thenReturn(TaxTriageOptions.SequencingPreset.ONT);
        int time4 = workflowRunner.estimateExecutionTime(mockDocuments, mockOptions);
        assertTrue(time4 > time2);

        // Test maximum cap
        AnnotatedPluginDocument[] manyDocs = new AnnotatedPluginDocument[1000];
        for (int i = 0; i < 1000; i++) {
            manyDocs[i] = mock(AnnotatedPluginDocument.class);
        }
        int time5 = workflowRunner.estimateExecutionTime(manyDocs, mockOptions);
        assertEquals(240, time5); // Should be capped at 4 hours
    }

    @Test
    void testSystemRequirementsCheck() {
        TaxTriageWorkflowRunner.SystemStatus status = workflowRunner.checkSystemRequirements();
        assertNotNull(status);
        assertNotNull(status.getSummary());

        // Should have some memory available
        assertTrue(status.availableMemoryMB > 0);

        // Should have some disk space
        assertTrue(status.availableDiskSpaceGB >= 0);
    }

    @Test
    void testTempBaseDirectory() {
        assertEquals(tempDir, workflowRunner.getTempBaseDirectory());
    }

    @Test
    void testWorkflowWithNullProgressListener() throws Exception {
        WorkflowContext successContext = createSuccessfulContext();
        CompletableFuture<WorkflowContext> successFuture = CompletableFuture.completedFuture(successContext);
        when(mockWorkflowExecutor.executeWorkflow(any(WorkflowContext.class), any()))
                .thenReturn(successFuture);

        // Execute workflow without progress listener
        WorkflowContext result = workflowRunner.runWorkflow(mockDocuments, mockOptions, null);

        assertNotNull(result);
        assertTrue(result.isSuccessful());
    }

    @Test
    void testSystemStatusSummary() {
        TaxTriageWorkflowRunner.SystemStatus status = new TaxTriageWorkflowRunner.SystemStatus();
        status.overall = true;
        status.dockerAvailable = true;
        status.dockerMessage = "Docker OK";
        status.sufficientDiskSpace = true;
        status.availableDiskSpaceGB = 100;
        status.sufficientMemory = true;
        status.availableMemoryMB = 4096;

        String summary = status.getSummary();
        assertTrue(summary.contains("Ready"));
        assertTrue(summary.contains("Docker OK"));
        assertTrue(summary.contains("100 GB"));
        assertTrue(summary.contains("4096 MB"));
    }

    @Test
    void testSystemStatusWithInsufficientResources() {
        TaxTriageWorkflowRunner.SystemStatus status = new TaxTriageWorkflowRunner.SystemStatus();
        status.overall = false;
        status.dockerAvailable = false;
        status.dockerMessage = "Docker not found";
        status.sufficientDiskSpace = false;
        status.availableDiskSpaceGB = 2;
        status.sufficientMemory = false;
        status.availableMemoryMB = 1024;

        String summary = status.getSummary();
        assertTrue(summary.contains("Not Ready"));
        assertTrue(summary.contains("Not Available"));
        assertTrue(summary.contains("Insufficient"));
    }

    private WorkflowContext createSuccessfulContext() {
        WorkflowContext context = mock(WorkflowContext.class);
        when(context.isSuccessful()).thenReturn(true);
        when(context.getState()).thenReturn(WorkflowContext.ExecutionState.COMPLETED);
        when(context.getWorkflowId()).thenReturn("test-workflow-123");
        when(context.getStartedAt()).thenReturn(LocalDateTime.now().minusMinutes(5));
        when(context.getCompletedAt()).thenReturn(LocalDateTime.now());
        return context;
    }

    private WorkflowContext createFailedContext() {
        WorkflowContext context = mock(WorkflowContext.class);
        when(context.isSuccessful()).thenReturn(false);
        when(context.getState()).thenReturn(WorkflowContext.ExecutionState.FAILED);
        when(context.getLastError()).thenReturn(new RuntimeException("Workflow failed"));
        when(context.getWorkflowId()).thenReturn("test-workflow-456");
        return context;
    }
}