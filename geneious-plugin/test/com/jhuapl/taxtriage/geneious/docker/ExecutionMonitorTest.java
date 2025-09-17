package com.jhuapl.taxtriage.geneious.docker;

import jebl.util.ProgressListener;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.ArgumentCaptor;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import java.io.ByteArrayInputStream;
import java.io.InputStream;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.*;
import static org.mockito.Mockito.*;

/**
 * Comprehensive test suite for ExecutionMonitor functionality.
 *
 * Tests progress parsing, error detection, cancellation handling, and process
 * monitoring capabilities for Nextflow workflow execution.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class ExecutionMonitorTest {

    @Mock
    private ProgressListener mockProgressListener;

    @Mock
    private Process mockProcess;

    private ExecutionMonitor executionMonitor;

    @BeforeEach
    void setUp() {
        MockitoAnnotations.openMocks(this);
        executionMonitor = new ExecutionMonitor();
    }

    @Test
    void testParseProgressFromValidLine() {
        String progressLine = "[task-123] process > PREPROCESSING [5 of 10]";
        double progress = executionMonitor.parseProgress(progressLine);

        assertEquals(0.5, progress, 0.001);
    }

    @Test
    void testParseProgressFromCompletionLine() {
        String completionLine = "Pipeline completed at: 2023-12-01 15:30:45";
        double progress = executionMonitor.parseProgress(completionLine);

        assertEquals(1.0, progress, 0.001);
    }

    @Test
    void testParseProgressFromInvalidLine() {
        String invalidLine = "Some random output line";
        double progress = executionMonitor.parseProgress(invalidLine);

        assertEquals(-1, progress, 0.001);
    }

    @Test
    void testParseProgressFromNullLine() {
        double progress = executionMonitor.parseProgress(null);
        assertEquals(-1, progress, 0.001);
    }

    @Test
    void testParseProgressFromEmptyLine() {
        double progress = executionMonitor.parseProgress("");
        assertEquals(-1, progress, 0.001);
    }

    @Test
    void testParseProgressWithDifferentFormats() {
        // Test various progress line formats
        assertEquals(0.25, executionMonitor.parseProgress("[abc-456] process > ANALYSIS [1 of 4]"), 0.001);
        assertEquals(1.0, executionMonitor.parseProgress("[def-789] process > REPORTING [3 of 3]"), 0.001);
        assertEquals(0.0, executionMonitor.parseProgress("[ghi-012] process > SETUP [0 of 5]"), 0.001);
    }

    @Test
    void testIsErrorLine() {
        assertTrue(executionMonitor.isErrorLine("ERROR: Something went wrong"));
        assertTrue(executionMonitor.isErrorLine("WARN: Warning message"));
        assertTrue(executionMonitor.isErrorLine("Exception in thread main"));
        assertTrue(executionMonitor.isErrorLine("Failed to process file"));
        assertTrue(executionMonitor.isErrorLine("Error occurred during execution"));

        assertFalse(executionMonitor.isErrorLine("INFO: Process completed successfully"));
        assertFalse(executionMonitor.isErrorLine("DEBUG: Debugging information"));
        assertFalse(executionMonitor.isErrorLine("Normal output line"));
        assertFalse(executionMonitor.isErrorLine(null));
        assertFalse(executionMonitor.isErrorLine(""));
    }

    @Test
    void testExtractProcessName() {
        String processLine = "[task-123] process > PREPROCESSING [2 of 5]";
        String processName = executionMonitor.extractProcessName(processLine);

        assertEquals("PREPROCESSING", processName);
    }

    @Test
    void testExtractProcessNameFromInvalidLine() {
        String invalidLine = "Some random output";
        String processName = executionMonitor.extractProcessName(invalidLine);

        assertNull(processName);
    }

    @Test
    void testExtractProcessNameFromNullLine() {
        String processName = executionMonitor.extractProcessName(null);
        assertNull(processName);
    }

    @Test
    void testIsCompletionLine() {
        assertTrue(executionMonitor.isCompletionLine("Pipeline completed at: 2023-12-01 15:30:45"));
        assertTrue(executionMonitor.isCompletionLine("Pipeline completed successfully"));
        assertTrue(executionMonitor.isCompletionLine("Completed at: 2023-12-01"));
        assertTrue(executionMonitor.isCompletionLine("Duration: 5m 30s"));
        assertTrue(executionMonitor.isCompletionLine("Succeeded: 10"));

        assertFalse(executionMonitor.isCompletionLine("Process still running"));
        assertFalse(executionMonitor.isCompletionLine("Starting pipeline"));
        assertFalse(executionMonitor.isCompletionLine(null));
        assertFalse(executionMonitor.isCompletionLine(""));
    }

    @Test
    void testRequestCancellation() {
        assertFalse(executionMonitor.isCancellationRequested());

        executionMonitor.requestCancellation();

        assertTrue(executionMonitor.isCancellationRequested());
    }

    @Test
    void testReset() {
        executionMonitor.requestCancellation();
        assertTrue(executionMonitor.isCancellationRequested());

        executionMonitor.reset();

        assertFalse(executionMonitor.isCancellationRequested());
    }

    @Test
    void testMonitorExecutionSuccess() throws Exception {
        // Mock process streams
        String stdout = "INFO: Starting process\n[task-1] process > ANALYSIS [1 of 2]\n[task-2] process > ANALYSIS [2 of 2]\nPipeline completed at: 2023-12-01";
        String stderr = "";

        InputStream stdoutStream = new ByteArrayInputStream(stdout.getBytes());
        InputStream stderrStream = new ByteArrayInputStream(stderr.getBytes());

        when(mockProcess.getInputStream()).thenReturn(stdoutStream);
        when(mockProcess.getErrorStream()).thenReturn(stderrStream);
        when(mockProcess.isAlive()).thenReturn(false);
        when(mockProcess.waitFor()).thenReturn(0);

        CompletableFuture<ExecutionResult> future = executionMonitor.monitorExecution(mockProcess, mockProgressListener);
        ExecutionResult result = future.get();

        assertNotNull(result);
        assertEquals(0, result.getExitCode());
        assertTrue(result.isSuccess());
        assertTrue(result.getOutput().contains("Starting process"));
        assertTrue(result.getOutput().contains("Pipeline completed"));

        // Verify progress listener was called
        verify(mockProgressListener, atLeastOnce()).setProgress(anyDouble());
        verify(mockProgressListener, atLeastOnce()).setMessage(anyString());
    }

    @Test
    void testMonitorExecutionFailure() throws Exception {
        // Mock process streams with error
        String stdout = "INFO: Starting process\n";
        String stderr = "ERROR: Process failed\n";

        InputStream stdoutStream = new ByteArrayInputStream(stdout.getBytes());
        InputStream stderrStream = new ByteArrayInputStream(stderr.getBytes());

        when(mockProcess.getInputStream()).thenReturn(stdoutStream);
        when(mockProcess.getErrorStream()).thenReturn(stderrStream);
        when(mockProcess.isAlive()).thenReturn(false);
        when(mockProcess.waitFor()).thenReturn(1);

        CompletableFuture<ExecutionResult> future = executionMonitor.monitorExecution(mockProcess, mockProgressListener);
        ExecutionResult result = future.get();

        assertNotNull(result);
        assertEquals(1, result.getExitCode());
        assertFalse(result.isSuccess());
        assertTrue(result.getErrorOutput().contains("ERROR: Process failed"));
    }

    @Test
    void testMonitorExecutionWithCancellation() throws Exception {
        // Mock long-running process
        when(mockProcess.getInputStream()).thenReturn(new ByteArrayInputStream("".getBytes()));
        when(mockProcess.getErrorStream()).thenReturn(new ByteArrayInputStream("".getBytes()));
        when(mockProcess.isAlive()).thenReturn(true); // Process never terminates naturally

        CompletableFuture<ExecutionResult> future = executionMonitor.monitorExecution(mockProcess, mockProgressListener);

        // Request cancellation
        executionMonitor.requestCancellation();

        // Should complete with exception due to cancellation
        assertThrows(ExecutionException.class, () -> {
            future.get();
        });

        // Verify process was destroyed
        verify(mockProcess).destroyForcibly();
    }

    @Test
    void testMonitorExecutionWithoutProgressListener() throws Exception {
        String stdout = "[task-1] process > ANALYSIS [1 of 1]\nPipeline completed";
        InputStream stdoutStream = new ByteArrayInputStream(stdout.getBytes());
        InputStream stderrStream = new ByteArrayInputStream("".getBytes());

        when(mockProcess.getInputStream()).thenReturn(stdoutStream);
        when(mockProcess.getErrorStream()).thenReturn(stderrStream);
        when(mockProcess.isAlive()).thenReturn(false);
        when(mockProcess.waitFor()).thenReturn(0);

        CompletableFuture<ExecutionResult> future = executionMonitor.monitorExecution(mockProcess, null);
        ExecutionResult result = future.get();

        assertNotNull(result);
        assertEquals(0, result.getExitCode());
        assertTrue(result.isSuccess());
    }

    @Test
    void testProgressListenerUpdates() throws Exception {
        String stdout = "Launching workflow...\n" +
                       "[task-1] process > PREPROCESSING [1 of 4]\n" +
                       "[task-2] process > ANALYSIS [2 of 4]\n" +
                       "[task-3] process > FILTERING [3 of 4]\n" +
                       "[task-4] process > REPORTING [4 of 4]\n" +
                       "Pipeline completed at: 2023-12-01";

        InputStream stdoutStream = new ByteArrayInputStream(stdout.getBytes());
        InputStream stderrStream = new ByteArrayInputStream("".getBytes());

        when(mockProcess.getInputStream()).thenReturn(stdoutStream);
        when(mockProcess.getErrorStream()).thenReturn(stderrStream);
        when(mockProcess.isAlive()).thenReturn(false);
        when(mockProcess.waitFor()).thenReturn(0);

        CompletableFuture<ExecutionResult> future = executionMonitor.monitorExecution(mockProcess, mockProgressListener);
        ExecutionResult result = future.get();

        // Capture progress updates
        ArgumentCaptor<Double> progressCaptor = ArgumentCaptor.forClass(Double.class);
        verify(mockProgressListener, atLeastOnce()).setProgress(progressCaptor.capture());

        // Should have progress values: 0.25, 0.5, 0.75, 1.0
        assertTrue(progressCaptor.getAllValues().contains(0.25));
        assertTrue(progressCaptor.getAllValues().contains(0.5));
        assertTrue(progressCaptor.getAllValues().contains(0.75));
        assertTrue(progressCaptor.getAllValues().contains(1.0));

        // Capture message updates
        ArgumentCaptor<String> messageCaptor = ArgumentCaptor.forClass(String.class);
        verify(mockProgressListener, atLeastOnce()).setMessage(messageCaptor.capture());

        // Should have process-specific messages
        assertTrue(messageCaptor.getAllValues().stream()
                .anyMatch(msg -> msg.contains("PREPROCESSING")));
        assertTrue(messageCaptor.getAllValues().stream()
                .anyMatch(msg -> msg.contains("completed")));
    }

    @Test
    void testErrorDetectionInStream() throws Exception {
        String stderr = "WARN: Low memory warning\nERROR: Critical failure\nException: Null pointer";

        InputStream stdoutStream = new ByteArrayInputStream("".getBytes());
        InputStream stderrStream = new ByteArrayInputStream(stderr.getBytes());

        when(mockProcess.getInputStream()).thenReturn(stdoutStream);
        when(mockProcess.getErrorStream()).thenReturn(stderrStream);
        when(mockProcess.isAlive()).thenReturn(false);
        when(mockProcess.waitFor()).thenReturn(1);

        CompletableFuture<ExecutionResult> future = executionMonitor.monitorExecution(mockProcess, mockProgressListener);
        ExecutionResult result = future.get();

        assertFalse(result.isSuccess());
        assertTrue(result.getErrorOutput().contains("WARN: Low memory warning"));
        assertTrue(result.getErrorOutput().contains("ERROR: Critical failure"));
        assertTrue(result.getErrorOutput().contains("Exception: Null pointer"));
    }

    @Test
    void testSpecialNextflowMessages() throws Exception {
        String stdout = "Pulling container image...\n" +
                       "executor > local\n" +
                       "Launching pipeline...\n" +
                       "[task-1] process > MAIN [1 of 1]\n" +
                       "Pipeline completed";

        InputStream stdoutStream = new ByteArrayInputStream(stdout.getBytes());
        InputStream stderrStream = new ByteArrayInputStream("".getBytes());

        when(mockProcess.getInputStream()).thenReturn(stdoutStream);
        when(mockProcess.getErrorStream()).thenReturn(stderrStream);
        when(mockProcess.isAlive()).thenReturn(false);
        when(mockProcess.waitFor()).thenReturn(0);

        CompletableFuture<ExecutionResult> future = executionMonitor.monitorExecution(mockProcess, mockProgressListener);
        ExecutionResult result = future.get();

        ArgumentCaptor<String> messageCaptor = ArgumentCaptor.forClass(String.class);
        verify(mockProgressListener, atLeastOnce()).setMessage(messageCaptor.capture());

        // Check for specific message types
        assertTrue(messageCaptor.getAllValues().stream()
                .anyMatch(msg -> msg.contains("Pulling")));
        assertTrue(messageCaptor.getAllValues().stream()
                .anyMatch(msg -> msg.contains("executor")));
        assertTrue(messageCaptor.getAllValues().stream()
                .anyMatch(msg -> msg.contains("Launching")));
    }
}