package com.jhuapl.taxtriage.geneious.execution;

import com.jhuapl.taxtriage.geneious.docker.ExecutionResult;
import org.junit.jupiter.api.Assertions;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.time.temporal.ChronoUnit;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.TimeUnit;

/**
 * Utility methods for workflow testing.
 *
 * This class provides common assertions, validations, and helper methods
 * used across multiple workflow test classes.
 */
public class WorkflowTestUtils {

    /**
     * Asserts that a workflow context represents a successful execution.
     *
     * @param context the workflow context to check
     */
    public static void assertWorkflowSuccess(WorkflowContext context) {
        Assertions.assertNotNull(context, "Workflow context should not be null");
        Assertions.assertTrue(context.isSuccessful(), "Workflow should be successful");
        Assertions.assertEquals(WorkflowContext.ExecutionState.COMPLETED, context.getState(),
                "Workflow state should be COMPLETED");
        Assertions.assertNotNull(context.getStartedAt(), "Start time should be set");
        Assertions.assertNotNull(context.getCompletedAt(), "Completion time should be set");
        Assertions.assertNull(context.getLastError(), "No error should be present");
        Assertions.assertEquals(1.0, context.getProgress(), 0.001, "Progress should be 100%");
    }

    /**
     * Asserts that a workflow context represents a failed execution.
     *
     * @param context the workflow context to check
     */
    public static void assertWorkflowFailure(WorkflowContext context) {
        Assertions.assertNotNull(context, "Workflow context should not be null");
        Assertions.assertFalse(context.isSuccessful(), "Workflow should not be successful");
        Assertions.assertEquals(WorkflowContext.ExecutionState.FAILED, context.getState(),
                "Workflow state should be FAILED");
        Assertions.assertNotNull(context.getLastError(), "Error should be present");
        Assertions.assertNotNull(context.getCompletedAt(), "Completion time should be set");
    }

    /**
     * Asserts that execution time is within reasonable bounds.
     *
     * @param context the workflow context
     * @param maxDurationMinutes maximum expected duration in minutes
     */
    public static void assertExecutionTime(WorkflowContext context, long maxDurationMinutes) {
        if (context.getStartedAt() != null && context.getCompletedAt() != null) {
            long durationMinutes = ChronoUnit.MINUTES.between(context.getStartedAt(), context.getCompletedAt());
            Assertions.assertTrue(durationMinutes <= maxDurationMinutes,
                    String.format("Execution took %d minutes, expected <= %d minutes",
                            durationMinutes, maxDurationMinutes));
        }
    }

    /**
     * Asserts that required directories were created.
     *
     * @param context the workflow context
     */
    public static void assertDirectoriesCreated(WorkflowContext context) {
        Assertions.assertTrue(Files.exists(context.getWorkingDirectory()),
                "Working directory should exist");
        Assertions.assertTrue(Files.exists(context.getInputDirectory()),
                "Input directory should exist");
        Assertions.assertTrue(Files.exists(context.getOutputDirectory()),
                "Output directory should exist");
        Assertions.assertTrue(Files.exists(context.getConfigDirectory()),
                "Config directory should exist");
    }

    /**
     * Asserts that configuration files were generated.
     *
     * @param context the workflow context
     */
    public static void assertConfigFilesGenerated(WorkflowContext context) {
        Assertions.assertNotNull(context.getSampleSheetFile(), "Samplesheet file should be set");
        Assertions.assertNotNull(context.getConfigFile(), "Config file should be set");
        Assertions.assertNotNull(context.getParamsFile(), "Params file should be set");

        Assertions.assertTrue(context.getSampleSheetFile().exists(),
                "Samplesheet file should exist");
        Assertions.assertTrue(context.getConfigFile().exists(),
                "Config file should exist");
        Assertions.assertTrue(context.getParamsFile().exists(),
                "Params file should exist");
    }

    /**
     * Asserts that sequence files were exported.
     *
     * @param context the workflow context
     * @param expectedCount expected number of exported files
     */
    public static void assertSequencesExported(WorkflowContext context, int expectedCount) {
        Assertions.assertNotNull(context.getExportedFiles(), "Exported files should be set");
        Assertions.assertEquals(expectedCount, context.getExportedFiles().size(),
                "Expected number of exported files");

        for (File exportedFile : context.getExportedFiles()) {
            Assertions.assertTrue(exportedFile.exists(),
                    "Exported file should exist: " + exportedFile.getName());
            Assertions.assertTrue(exportedFile.length() > 0,
                    "Exported file should not be empty: " + exportedFile.getName());
        }
    }

    /**
     * Asserts that progress updates were made.
     *
     * @param progressListener the mock progress listener
     * @param minUpdates minimum expected number of updates
     */
    public static void assertProgressUpdates(MockProgressListener progressListener, int minUpdates) {
        Assertions.assertNotNull(progressListener, "Progress listener should not be null");
        Assertions.assertTrue(progressListener.getProgressUpdateCount() >= minUpdates,
                String.format("Expected at least %d progress updates, got %d",
                        minUpdates, progressListener.getProgressUpdateCount()));
        Assertions.assertTrue(progressListener.getMessageUpdateCount() >= minUpdates,
                String.format("Expected at least %d message updates, got %d",
                        minUpdates, progressListener.getMessageUpdateCount()));
    }

    /**
     * Asserts that progress is monotonically increasing.
     *
     * @param progressListener the mock progress listener
     */
    public static void assertMonotonicProgress(MockProgressListener progressListener) {
        Assertions.assertTrue(progressListener.isProgressMonotonic(),
                "Progress should be monotonically increasing");
    }

    /**
     * Asserts that the final progress is 100%.
     *
     * @param progressListener the mock progress listener
     */
    public static void assertCompleteProgress(MockProgressListener progressListener) {
        Assertions.assertEquals(1.0, progressListener.getFinalProgress(), 0.001,
                "Final progress should be 100%");
    }

    /**
     * Asserts that specific progress milestones were reached.
     *
     * @param progressListener the mock progress listener
     * @param milestones the progress milestones to check
     */
    public static void assertProgressMilestones(MockProgressListener progressListener, double... milestones) {
        for (double milestone : milestones) {
            Assertions.assertTrue(progressListener.progressReached(milestone),
                    String.format("Progress milestone %.1f%% was not reached", milestone * 100));
        }
    }

    /**
     * Asserts that specific messages were sent.
     *
     * @param progressListener the mock progress listener
     * @param expectedMessages the messages that should have been sent
     */
    public static void assertMessagesContain(MockProgressListener progressListener, String... expectedMessages) {
        for (String expectedMessage : expectedMessages) {
            Assertions.assertTrue(progressListener.messageContains(expectedMessage),
                    "Expected message not found: " + expectedMessage);
        }
    }

    /**
     * Waits for a CompletableFuture to complete with timeout.
     *
     * @param future the future to wait for
     * @param timeoutSeconds timeout in seconds
     * @param <T> the type of the future result
     * @return the result of the future
     * @throws Exception if the future fails or times out
     */
    public static <T> T waitForCompletion(CompletableFuture<T> future, int timeoutSeconds) throws Exception {
        return future.get(timeoutSeconds, TimeUnit.SECONDS);
    }

    /**
     * Creates a temporary directory and ensures it's cleaned up.
     *
     * @param prefix the directory name prefix
     * @return the created temporary directory
     * @throws IOException if directory creation fails
     */
    public static Path createTempDirectory(String prefix) throws IOException {
        Path tempDir = Files.createTempDirectory(prefix);
        // Add shutdown hook to clean up
        Runtime.getRuntime().addShutdownHook(new Thread(() -> {
            try {
                deleteDirectoryRecursively(tempDir);
            } catch (IOException e) {
                System.err.println("Failed to clean up temp directory: " + e.getMessage());
            }
        }));
        return tempDir;
    }

    /**
     * Recursively deletes a directory and all its contents.
     *
     * @param directory the directory to delete
     * @throws IOException if deletion fails
     */
    public static void deleteDirectoryRecursively(Path directory) throws IOException {
        if (Files.exists(directory)) {
            Files.walk(directory)
                    .sorted((a, b) -> b.compareTo(a)) // Delete files before directories
                    .forEach(path -> {
                        try {
                            Files.deleteIfExists(path);
                        } catch (IOException e) {
                            System.err.println("Failed to delete: " + path + " - " + e.getMessage());
                        }
                    });
        }
    }

    /**
     * Validates the content of a FASTQ file.
     *
     * @param fastqFile the FASTQ file to validate
     * @throws IOException if file reading fails
     */
    public static void validateFastqFile(File fastqFile) throws IOException {
        Assertions.assertTrue(fastqFile.exists(), "FASTQ file should exist");
        Assertions.assertTrue(fastqFile.length() > 0, "FASTQ file should not be empty");

        List<String> lines = Files.readAllLines(fastqFile.toPath());
        Assertions.assertTrue(lines.size() % 4 == 0, "FASTQ file should have multiple of 4 lines");

        for (int i = 0; i < lines.size(); i += 4) {
            Assertions.assertTrue(lines.get(i).startsWith("@"),
                    "FASTQ header line should start with @");
            Assertions.assertFalse(lines.get(i + 1).isEmpty(),
                    "FASTQ sequence line should not be empty");
            Assertions.assertTrue(lines.get(i + 2).startsWith("+"),
                    "FASTQ plus line should start with +");
            Assertions.assertEquals(lines.get(i + 1).length(), lines.get(i + 3).length(),
                    "FASTQ sequence and quality lines should have same length");
        }
    }

    /**
     * Creates a mock successful execution result.
     *
     * @param output the output message
     * @return ExecutionResult representing success
     */
    public static ExecutionResult createSuccessResult(String output) {
        return new ExecutionResult(0, output, "");
    }

    /**
     * Creates a mock failed execution result.
     *
     * @param errorMessage the error message
     * @return ExecutionResult representing failure
     */
    public static ExecutionResult createFailureResult(String errorMessage) {
        return new ExecutionResult(1, "", errorMessage);
    }

    /**
     * Simulates workflow execution time by checking timestamps.
     *
     * @param context the workflow context
     * @param minExecutionMillis minimum expected execution time in milliseconds
     */
    public static void assertMinimumExecutionTime(WorkflowContext context, long minExecutionMillis) {
        if (context.getStartedAt() != null && context.getCompletedAt() != null) {
            long actualMillis = ChronoUnit.MILLIS.between(context.getStartedAt(), context.getCompletedAt());
            Assertions.assertTrue(actualMillis >= minExecutionMillis,
                    String.format("Execution took %d ms, expected at least %d ms",
                            actualMillis, minExecutionMillis));
        }
    }

    /**
     * Validates that all expected workflow phases were executed.
     *
     * @param progressListener the progress listener
     */
    public static void assertAllPhasesExecuted(MockProgressListener progressListener) {
        assertMessagesContain(progressListener,
                "Exporting sequences",
                "Generating workflow configuration",
                "Creating samplesheet",
                "Executing TaxTriage workflow",
                "Processing workflow results"
        );
    }

    /**
     * Creates a workflow context that will succeed.
     *
     * @param tempDir temporary directory for the workflow
     * @return configured successful WorkflowContext
     * @throws IOException if setup fails
     */
    public static WorkflowContext createSuccessfulContext(Path tempDir) throws IOException {
        WorkflowContext context = TestDataGenerator.createRealisticWorkflowContext(
                tempDir, 2, com.jhuapl.taxtriage.geneious.TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        // Simulate successful completion
        context.setStartedAt(LocalDateTime.now().minusMinutes(5));
        context.markCompleted();

        // Add some mock result files
        TestDataGenerator.createMockWorkflowReport(
                context.getOutputDirectory().resolve("workflow_report.html"),
                context.getWorkflowId(), true);
        TestDataGenerator.createMockTaxonomyResults(
                context.getOutputDirectory().resolve("taxonomy_classification.csv"));

        context.setWorkflowReport(context.getOutputDirectory().resolve("workflow_report.html").toFile());
        context.setTaxonomyResults(context.getOutputDirectory().resolve("taxonomy_classification.csv").toFile());

        return context;
    }

    /**
     * Creates a workflow context that will fail.
     *
     * @param tempDir temporary directory for the workflow
     * @return configured failed WorkflowContext
     * @throws IOException if setup fails
     */
    public static WorkflowContext createFailedContext(Path tempDir) throws IOException {
        WorkflowContext context = TestDataGenerator.createRealisticWorkflowContext(
                tempDir, 1, com.jhuapl.taxtriage.geneious.TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

        // Simulate failure
        context.setStartedAt(LocalDateTime.now().minusMinutes(2));
        context.markFailed(new RuntimeException("Simulated workflow failure"), "Workflow execution failed");

        return context;
    }
}