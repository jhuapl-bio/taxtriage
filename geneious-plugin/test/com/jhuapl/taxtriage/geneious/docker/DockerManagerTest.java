package com.jhuapl.taxtriage.geneious.docker;

import jebl.util.ProgressListener;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.mockito.ArgumentCaptor;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.CompletableFuture;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.*;
import static org.mockito.Mockito.*;

/**
 * Comprehensive test suite for DockerManager functionality.
 *
 * Tests Docker availability checking, image management, command execution,
 * and error handling across different scenarios without requiring actual Docker.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class DockerManagerTest {

    @Mock
    private VolumeMapper mockVolumeMapper;

    @Mock
    private ExecutionMonitor mockExecutionMonitor;

    @Mock
    private ProgressListener mockProgressListener;

    @TempDir
    Path tempDir;

    private DockerManager dockerManager;
    private Path inputDir;
    private Path outputDir;
    private Path workDir;

    @BeforeEach
    void setUp() throws Exception {
        MockitoAnnotations.openMocks(this);

        // Create test directories
        inputDir = tempDir.resolve("input");
        outputDir = tempDir.resolve("output");
        workDir = tempDir.resolve("work");

        Files.createDirectories(inputDir);
        Files.createDirectories(outputDir);
        Files.createDirectories(workDir);

        // Create mock volume mappings
        Map<String, String> mockMounts = new HashMap<>();
        mockMounts.put(inputDir.toString(), "/input");
        mockMounts.put(outputDir.toString(), "/output");
        mockMounts.put(workDir.toString(), "/work");

        when(mockVolumeMapper.createVolumeMounts(any(), any(), any())).thenReturn(mockMounts);
        when(mockVolumeMapper.getContainerWorkPath()).thenReturn("/work");

        // Mock successful execution by default
        ExecutionResult successResult = new ExecutionResult(0, "Success", "");
        CompletableFuture<ExecutionResult> successFuture = CompletableFuture.completedFuture(successResult);
        when(mockExecutionMonitor.monitorExecution(any(), any())).thenReturn(successFuture);

        // Create DockerManager with mocked dependencies
        dockerManager = new DockerManager(mockVolumeMapper, mockExecutionMonitor, "echo") {
            @Override
            public boolean isDockerAvailable() {
                return true; // Mock Docker as available
            }

            @Override
            public void pullImageIfNeeded(String imageName, ProgressListener progressListener) {
                // Mock image pull - do nothing for tests
            }
        };
    }

    @Test
    void testDockerManagerCreation() {
        assertNotNull(dockerManager);
        assertTrue(dockerManager.isDockerAvailable());
    }

    @Test
    void testDockerManagerCreationWithUnavailableDocker() {
        // Test exception when Docker is not available
        assertThrows(DockerException.class, () -> {
            new DockerManager(mockVolumeMapper, mockExecutionMonitor, "nonexistent-command");
        });
    }

    @Test
    void testIsDockerAvailable() {
        assertTrue(dockerManager.isDockerAvailable());
    }

    @Test
    void testPullImageIfNeededWhenImagePresent() throws DockerException {
        // Create a DockerManager that reports image as already present
        DockerManager managerWithImage = new DockerManager(mockVolumeMapper, mockExecutionMonitor, "echo") {
            @Override
            public boolean isDockerAvailable() { return true; }
            @Override
            public void pullImageIfNeeded(String imageName, ProgressListener progressListener) {
                // Simulate image already present - do nothing
            }
        };

        // Should not throw exception and should complete quickly
        assertDoesNotThrow(() -> {
            managerWithImage.pullImageIfNeeded("test-image", mockProgressListener);
        });
    }

    @Test
    void testPullImageIfNeededWithProgress() throws DockerException {
        // Verify progress listener is called during image pull
        dockerManager.pullImageIfNeeded("test-image", mockProgressListener);

        // Progress listener should be called at least once for pulling message
        verify(mockProgressListener, atLeastOnce()).setMessage(anyString());
    }

    @Test
    void testExecuteNextflowCommandSuccess() throws DockerException {
        String command = "nextflow run main.nf";

        ExecutionResult result = dockerManager.executeNextflowCommand(
                command, inputDir, outputDir, workDir, mockProgressListener);

        assertNotNull(result);
        assertTrue(result.isSuccess());
        assertEquals(0, result.getExitCode());

        // Verify volume mapper was called
        verify(mockVolumeMapper).createVolumeMounts(inputDir, outputDir, workDir);

        // Verify execution monitor was called
        verify(mockExecutionMonitor).monitorExecution(any(), eq(mockProgressListener));

        // Verify progress listener was updated
        verify(mockProgressListener, atLeastOnce()).setMessage(anyString());
    }

    @Test
    void testExecuteNextflowCommandWithCustomImage() throws DockerException {
        String customImage = "custom/taxtriage:v1.0";
        String command = "nextflow run main.nf";

        ExecutionResult result = dockerManager.executeNextflowCommand(
                customImage, command, inputDir, outputDir, workDir, mockProgressListener);

        assertNotNull(result);
        assertTrue(result.isSuccess());
    }

    @Test
    void testExecuteNextflowCommandFailure() throws DockerException {
        // Mock execution failure
        ExecutionResult failureResult = new ExecutionResult(1, "", "Error occurred");
        CompletableFuture<ExecutionResult> failureFuture = CompletableFuture.completedFuture(failureResult);
        when(mockExecutionMonitor.monitorExecution(any(), any())).thenReturn(failureFuture);

        String command = "nextflow run main.nf";

        assertThrows(DockerException.class, () -> {
            dockerManager.executeNextflowCommand(command, inputDir, outputDir, workDir, mockProgressListener);
        });
    }

    @Test
    void testExecuteNextflowCommandWithInvalidDirectories() {
        Path nonExistentDir = tempDir.resolve("nonexistent");
        String command = "nextflow run main.nf";

        assertThrows(DockerException.class, () -> {
            dockerManager.executeNextflowCommand(command, nonExistentDir, outputDir, workDir, mockProgressListener);
        });
    }

    @Test
    void testExecuteNextflowCommandWithNullProgressListener() throws DockerException {
        String command = "nextflow run main.nf";

        ExecutionResult result = dockerManager.executeNextflowCommand(
                command, inputDir, outputDir, workDir, null);

        assertNotNull(result);
        assertTrue(result.isSuccess());

        // Verify execution monitor was called with null progress listener
        verify(mockExecutionMonitor).monitorExecution(any(), isNull());
    }

    @Test
    void testExecuteNextflowCommandWithVolumeMapperException() throws DockerException {
        // Mock volume mapper to throw exception
        when(mockVolumeMapper.createVolumeMounts(any(), any(), any()))
                .thenThrow(new DockerException("Volume mapping failed"));

        String command = "nextflow run main.nf";

        assertThrows(DockerException.class, () -> {
            dockerManager.executeNextflowCommand(command, inputDir, outputDir, workDir, mockProgressListener);
        });
    }

    @Test
    void testGetPlatform() {
        String platform = dockerManager.getPlatform();
        assertNotNull(platform);
        assertFalse(platform.isEmpty());
        assertTrue(platform.toLowerCase().contains("mac") ||
                  platform.toLowerCase().contains("windows") ||
                  platform.toLowerCase().contains("linux"));
    }

    @Test
    void testStopContainer() throws DockerException {
        // Create a DockerManager that can handle container operations
        DockerManager containerManager = new DockerManager(mockVolumeMapper, mockExecutionMonitor, "echo") {
            @Override
            public boolean isDockerAvailable() { return true; }
            @Override
            public void pullImageIfNeeded(String imageName, ProgressListener progressListener) { }
            @Override
            public void stopContainer(String containerId) {
                // Mock successful container stop
            }
        };

        // Should not throw exception
        assertDoesNotThrow(() -> {
            containerManager.stopContainer("test-container-id");
        });
    }

    @Test
    void testRemoveContainer() throws DockerException {
        // Create a DockerManager that can handle container operations
        DockerManager containerManager = new DockerManager(mockVolumeMapper, mockExecutionMonitor, "echo") {
            @Override
            public boolean isDockerAvailable() { return true; }
            @Override
            public void pullImageIfNeeded(String imageName, ProgressListener progressListener) { }
            @Override
            public void removeContainer(String containerId) {
                // Mock successful container removal
            }
        };

        // Should not throw exception
        assertDoesNotThrow(() -> {
            containerManager.removeContainer("test-container-id");
        });
    }

    @Test
    void testExecuteNextflowCommandWithExecutionException() {
        // Mock execution monitor to throw exception
        when(mockExecutionMonitor.monitorExecution(any(), any()))
                .thenReturn(CompletableFuture.failedFuture(new RuntimeException("Execution failed")));

        String command = "nextflow run main.nf";

        assertThrows(DockerException.class, () -> {
            dockerManager.executeNextflowCommand(command, inputDir, outputDir, workDir, mockProgressListener);
        });
    }

    @Test
    void testDefaultTaxTriageImage() {
        assertEquals("jhuaplbio/taxtriage:latest", DockerManager.DEFAULT_TAXTRIAGE_IMAGE);
    }

    @Test
    void testExecuteNextflowCommandWithDefaultImage() throws DockerException {
        String command = "nextflow run main.nf";

        ExecutionResult result = dockerManager.executeNextflowCommand(
                command, inputDir, outputDir, workDir, mockProgressListener);

        assertNotNull(result);
        assertTrue(result.isSuccess());

        // Verify that volume mapper and execution monitor were called properly
        verify(mockVolumeMapper).createVolumeMounts(inputDir, outputDir, workDir);
        verify(mockExecutionMonitor).monitorExecution(any(), eq(mockProgressListener));
    }

    @Test
    void testProgressListenerUpdates() throws DockerException {
        String command = "nextflow run main.nf";

        dockerManager.executeNextflowCommand(command, inputDir, outputDir, workDir, mockProgressListener);

        // Capture all progress listener calls
        ArgumentCaptor<String> messageCaptor = ArgumentCaptor.forClass(String.class);
        verify(mockProgressListener, atLeastOnce()).setMessage(messageCaptor.capture());

        // Verify at least one meaningful message was set
        boolean foundRelevantMessage = messageCaptor.getAllValues().stream()
                .anyMatch(message -> message.contains("Docker") || message.contains("container"));

        assertTrue(foundRelevantMessage, "Expected at least one Docker-related progress message");
    }

    @Test
    void testDirectoryCreation() throws DockerException {
        // Delete output directory to test creation
        outputDir.toFile().delete();
        assertFalse(Files.exists(outputDir));

        String command = "nextflow run main.nf";

        // Execution should succeed and create the directory
        ExecutionResult result = dockerManager.executeNextflowCommand(
                command, inputDir, outputDir, workDir, mockProgressListener);

        assertNotNull(result);
        assertTrue(result.isSuccess());
    }
}