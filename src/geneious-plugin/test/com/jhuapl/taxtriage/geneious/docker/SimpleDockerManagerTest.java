package com.jhuapl.taxtriage.geneious.docker;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Simple test class for DockerManager without external dependencies.
 * Demonstrates basic functionality and can be run manually.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class SimpleDockerManagerTest {

    public static void main(String[] args) {
        System.out.println("Running basic DockerManager tests...");

        try {
            testExecutionResult();
            testDockerException();
            testVolumeMapper();
            testExecutionMonitor();
            // testMockDockerManager(); // Skip mock test for now

            System.out.println("All basic tests passed!");

        } catch (Exception e) {
            System.err.println("Test failed: " + e.getMessage());
            e.printStackTrace();
        }
    }

    private static void testExecutionResult() {
        System.out.println("Testing ExecutionResult...");

        // Test successful result
        ExecutionResult success = ExecutionResult.success("test command", "Success output");
        assert success.isSuccess() : "Should be successful";
        assert success.getExitCode() == 0 : "Exit code should be 0";
        assert "Success output".equals(success.getStandardOutput()) : "Output should match";

        // Test failed result
        ExecutionResult failure = ExecutionResult.failure("test command", 1, "Error occurred");
        assert !failure.isSuccess() : "Should not be successful";
        assert failure.getExitCode() == 1 : "Exit code should be 1";
        assert failure.getErrorOutput().contains("Error") : "Should contain error";

        // Test null handling
        java.time.LocalDateTime now = java.time.LocalDateTime.now();
        ExecutionResult nullResult = new ExecutionResult("test command", 0, null, null, now, now);
        assert "".equals(nullResult.getStandardOutput()) : "Null output should become empty string";
        assert "".equals(nullResult.getErrorOutput()) : "Null error should become empty string";

        System.out.println("ExecutionResult tests passed!");
    }

    private static void testDockerException() {
        System.out.println("Testing DockerException...");

        // Test message constructor
        DockerException ex1 = new DockerException("Test message");
        assert "Test message".equals(ex1.getMessage()) : "Message should match";

        // Test message and cause constructor
        RuntimeException cause = new RuntimeException("Root cause");
        DockerException ex2 = new DockerException("Wrapper message", cause);
        assert "Wrapper message".equals(ex2.getMessage()) : "Message should match";
        assert cause.equals(ex2.getCause()) : "Cause should match";

        // Test cause constructor
        DockerException ex3 = new DockerException(cause);
        assert ex3.getCause().equals(cause) : "Cause should be preserved";

        System.out.println("DockerException tests passed!");
    }

    private static void testVolumeMapper() throws Exception {
        System.out.println("Testing VolumeMapper...");

        VolumeMapper mapper = new VolumeMapper();

        // Test platform detection
        String platform = System.getProperty("os.name").toLowerCase();
        boolean isPlatformDetected = mapper.isWindows() || mapper.isMacOS() || mapper.isLinux();
        assert isPlatformDetected : "Should detect at least one platform";

        // Test path escaping
        String pathWithSpaces = "/path/with spaces/file.txt";
        String escaped = mapper.escapePathForDocker(pathWithSpaces);
        assert escaped.contains("\"") : "Paths with spaces should be quoted";

        String normalPath = "/normal/path";
        String escapedNormal = mapper.escapePathForDocker(normalPath);
        assert normalPath.equals(escapedNormal) : "Normal paths should not be changed";

        // Test container path constants
        assert "/input".equals(mapper.getContainerInputPath()) : "Input path should be /input";
        assert "/output".equals(mapper.getContainerOutputPath()) : "Output path should be /output";
        assert "/work".equals(mapper.getContainerWorkPath()) : "Work path should be /work";

        System.out.println("VolumeMapper tests passed!");
    }

    private static void testExecutionMonitor() {
        System.out.println("Testing ExecutionMonitor...");

        ExecutionMonitor monitor = new ExecutionMonitor();

        // Test progress parsing
        String progressLine = "[task-123] process > PREPROCESSING [2 of 4]";
        double progress = monitor.parseProgress(progressLine);
        assert Math.abs(progress - 0.5) < 0.001 : "Progress should be 0.5 for 2 of 4";

        String completionLine = "Pipeline completed at: 2023-12-01 15:30:45";
        double completionProgress = monitor.parseProgress(completionLine);
        assert Math.abs(completionProgress - 1.0) < 0.001 : "Completion should return 1.0";

        // Test error detection
        assert monitor.isErrorLine("ERROR: Something went wrong") : "Should detect error line";
        assert monitor.isErrorLine("WARN: Warning message") : "Should detect warning line";
        assert !monitor.isErrorLine("INFO: Normal message") : "Should not detect info as error";

        // Test process name extraction
        String processName = monitor.extractProcessName(progressLine);
        assert "PREPROCESSING".equals(processName) : "Should extract process name";

        // Test cancellation
        assert !monitor.isCancellationRequested() : "Should not be cancelled initially";
        monitor.requestCancellation();
        assert monitor.isCancellationRequested() : "Should be cancelled after request";
        monitor.reset();
        assert !monitor.isCancellationRequested() : "Should be reset after reset()";

        System.out.println("ExecutionMonitor tests passed!");
    }

    private static void testMockDockerManager() throws Exception {
        System.out.println("Testing MockDockerManager...");
        // Skipped - requires separate compilation
        System.out.println("MockDockerManager tests skipped!");
    }
}