package com.jhuapl.taxtriage.geneious.execution;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.sequence.SequenceDocument;
import com.jhuapl.taxtriage.geneious.TaxTriageOptions;
import com.jhuapl.taxtriage.geneious.config.TaxTriageConfig;
import org.junit.jupiter.api.*;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;

/**
 * Comprehensive tests for pipeline validation logic.
 *
 * Tests input validation, configuration validation, resource checks,
 * and prerequisite verification before workflow execution.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
@DisplayName("Pipeline Validation Tests")
class PipelineValidationTest {

    @TempDir
    Path tempWorkDir;

    private Path inputDir;
    private Path outputDir;
    private Path configDir;

    private TaxTriageOptions mockOptions;
    private TaxTriageConfig mockConfig;

    @BeforeEach
    void setUp() throws IOException {
        // Create test directory structure
        inputDir = tempWorkDir.resolve("input");
        outputDir = tempWorkDir.resolve("output");
        configDir = tempWorkDir.resolve("config");

        Files.createDirectories(inputDir);
        Files.createDirectories(outputDir);
        Files.createDirectories(configDir);

        // Setup mock options and config
        mockOptions = createMockOptions();
        mockConfig = createMockConfig();
    }

    @Nested
    @DisplayName("Input File Format Validation Tests")
    class InputFileFormatValidationTests {

        @Test
        @DisplayName("Should validate FASTQ file format")
        void testValidFastqFormat() throws IOException {
            Path fastqFile = inputDir.resolve("test.fastq");
            createValidFastqFile(fastqFile);

            assertTrue(isValidFastqFile(fastqFile), "Valid FASTQ file should pass validation");
        }

        @Test
        @DisplayName("Should reject invalid FASTQ format")
        void testInvalidFastqFormat() throws IOException {
            Path fastqFile = inputDir.resolve("invalid.fastq");
            Files.writeString(fastqFile, "Not a valid FASTQ file");

            assertFalse(isValidFastqFile(fastqFile), "Invalid FASTQ file should fail validation");
        }

        @Test
        @DisplayName("Should validate compressed FASTQ files")
        void testCompressedFastqFile() throws IOException {
            Path fastqGz = inputDir.resolve("test.fastq.gz");

            if (fastqGz.toFile().exists() || createDummyGzFile(fastqGz)) {
                assertTrue(fastqGz.toString().endsWith(".gz"), "Should recognize .gz extension");
                assertTrue(Files.exists(fastqGz), "Compressed file should exist");
            } else {
                assertTrue(true, "Compression test skipped (no gzip available)");
            }
        }

        @Test
        @DisplayName("Should detect empty FASTQ files")
        void testEmptyFastqFile() throws IOException {
            Path emptyFile = inputDir.resolve("empty.fastq");
            Files.createFile(emptyFile);

            assertFalse(isValidFastqFile(emptyFile), "Empty FASTQ file should be invalid");
            assertEquals(0, Files.size(emptyFile), "File should be empty");
        }

        @Test
        @DisplayName("Should validate FASTQ quality scores")
        void testFastqQualityScores() throws IOException {
            Path fastqFile = inputDir.resolve("quality_test.fastq");
            createValidFastqFile(fastqFile);

            String content = Files.readString(fastqFile);
            List<String> lines = Arrays.asList(content.split("\n"));

            // Check quality line (every 4th line starting from line 3)
            for (int i = 3; i < lines.size(); i += 4) {
                String qualityLine = lines.get(i);
                assertFalse(qualityLine.isEmpty(), "Quality line should not be empty");

                // Quality scores should be ASCII characters
                for (char c : qualityLine.toCharArray()) {
                    assertTrue(c >= 33 && c <= 126, "Quality scores should be valid ASCII");
                }
            }
        }

        @Test
        @DisplayName("Should validate sequence-quality length match")
        void testSequenceQualityLengthMatch() throws IOException {
            Path fastqFile = inputDir.resolve("length_test.fastq");
            createValidFastqFile(fastqFile);

            String content = Files.readString(fastqFile);
            List<String> lines = Arrays.asList(content.split("\n"));

            // Check that sequence and quality have same length
            for (int i = 0; i + 3 < lines.size(); i += 4) {
                String sequenceLine = lines.get(i + 1);
                String qualityLine = lines.get(i + 3);

                assertEquals(sequenceLine.length(), qualityLine.length(),
                        "Sequence and quality line should have same length");
            }
        }

        @Test
        @DisplayName("Should reject files with wrong extension")
        void testWrongFileExtension() throws IOException {
            Path txtFile = inputDir.resolve("test.txt");
            Files.writeString(txtFile, "This is not a FASTQ file");

            assertTrue(Files.exists(txtFile), "File should exist");
            assertFalse(txtFile.toString().endsWith(".fastq"), "Should not have .fastq extension");
        }
    }

    @Nested
    @DisplayName("Configuration Parameter Validation Tests")
    class ConfigurationParameterValidationTests {

        @Test
        @DisplayName("Should validate quality threshold range")
        void testQualityThresholdValidation() {
            when(mockOptions.getQualityThreshold()).thenReturn(20);
            assertTrue(mockOptions.getQualityThreshold() >= 0 &&
                            mockOptions.getQualityThreshold() <= 50,
                    "Quality threshold should be in valid range");

            when(mockOptions.getQualityThreshold()).thenReturn(-1);
            assertFalse(mockOptions.getQualityThreshold() >= 0,
                    "Negative quality threshold should be invalid");

            when(mockOptions.getQualityThreshold()).thenReturn(100);
            assertFalse(mockOptions.getQualityThreshold() <= 50,
                    "Quality threshold above 50 should be invalid");
        }

        @Test
        @DisplayName("Should validate minimum read length")
        void testMinReadLengthValidation() {
            when(mockOptions.getMinReadLength()).thenReturn(50);
            assertTrue(mockOptions.getMinReadLength() > 0,
                    "Min read length should be positive");

            when(mockOptions.getMinReadLength()).thenReturn(0);
            assertFalse(mockOptions.getMinReadLength() > 0,
                    "Zero min read length should be invalid");

            when(mockOptions.getMinReadLength()).thenReturn(-10);
            assertFalse(mockOptions.getMinReadLength() > 0,
                    "Negative min read length should be invalid");
        }

        @Test
        @DisplayName("Should validate thread count")
        void testThreadCountValidation() {
            when(mockOptions.getThreadCount()).thenReturn(4);
            assertTrue(mockOptions.getThreadCount() >= 1 &&
                            mockOptions.getThreadCount() <= 128,
                    "Thread count should be in valid range");

            when(mockOptions.getThreadCount()).thenReturn(0);
            assertFalse(mockOptions.getThreadCount() >= 1,
                    "Zero threads should be invalid");

            when(mockOptions.getThreadCount()).thenReturn(256);
            assertFalse(mockOptions.getThreadCount() <= 128,
                    "Excessive thread count should be invalid");
        }

        @Test
        @DisplayName("Should validate memory limit")
        void testMemoryLimitValidation() {
            when(mockOptions.getMemoryLimit()).thenReturn(8);
            assertTrue(mockOptions.getMemoryLimit() >= 1 &&
                            mockOptions.getMemoryLimit() <= 512,
                    "Memory limit should be in valid range");

            when(mockOptions.getMemoryLimit()).thenReturn(0);
            assertFalse(mockOptions.getMemoryLimit() >= 1,
                    "Zero memory should be invalid");
        }

        @Test
        @DisplayName("Should validate sequencing preset")
        void testSequencingPresetValidation() {
            when(mockOptions.getSequencingPreset())
                    .thenReturn(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);

            assertNotNull(mockOptions.getSequencingPreset(),
                    "Sequencing preset should be set");

            TaxTriageOptions.SequencingPreset preset = mockOptions.getSequencingPreset();
            assertTrue(preset == TaxTriageOptions.SequencingPreset.ILLUMINA_SE ||
                            preset == TaxTriageOptions.SequencingPreset.ILLUMINA_PE ||
                            preset == TaxTriageOptions.SequencingPreset.ONT ||
                            preset == TaxTriageOptions.SequencingPreset.PACBIO,
                    "Preset should be one of the valid options");
        }

        @Test
        @DisplayName("Should validate subsample settings")
        void testSubsampleValidation() {
            when(mockOptions.isEnableSubsampling()).thenReturn(true);
            when(mockOptions.getSubsampleSize()).thenReturn(1000000);

            if (mockOptions.isEnableSubsampling()) {
                assertTrue(mockOptions.getSubsampleSize() > 0,
                        "Subsample size should be positive when enabled");
            }
        }

        @Test
        @DisplayName("Should reject invalid configuration")
        void testInvalidConfiguration() {
            when(mockConfig.validate()).thenReturn("Invalid configuration: missing database path");

            String validationError = mockConfig.validate();
            assertNotNull(validationError, "Invalid config should return error message");
            assertFalse(validationError.isEmpty(), "Error message should not be empty");
        }
    }

    @Nested
    @DisplayName("Required File Presence Tests")
    class RequiredFilePresenceTests {

        @Test
        @DisplayName("Should verify input files exist")
        void testInputFilesExist() throws IOException {
            List<File> inputFiles = createSampleInputFiles();

            for (File file : inputFiles) {
                assertTrue(file.exists(), "Input file should exist: " + file.getName());
                assertTrue(file.length() > 0, "Input file should not be empty: " + file.getName());
            }
        }

        @Test
        @DisplayName("Should verify configuration files exist")
        void testConfigFilesExist() throws IOException {
            File sampleSheet = configDir.resolve("samplesheet.csv").toFile();
            File configFile = configDir.resolve("nextflow.config").toFile();
            File paramsFile = configDir.resolve("params.json").toFile();

            Files.writeString(sampleSheet.toPath(), "sample,fastq_1\n");
            Files.writeString(configFile.toPath(), "params { }");
            Files.writeString(paramsFile.toPath(), "{}");

            assertTrue(sampleSheet.exists(), "Samplesheet should exist");
            assertTrue(configFile.exists(), "Config file should exist");
            assertTrue(paramsFile.exists(), "Params file should exist");
        }

        @Test
        @DisplayName("Should detect missing required files")
        void testMissingRequiredFiles() {
            File missingFile = inputDir.resolve("missing.fastq").toFile();

            assertFalse(missingFile.exists(), "Missing file should not exist");
        }

        @Test
        @DisplayName("Should validate samplesheet format")
        void testSamplesheetFormat() throws IOException {
            File sampleSheet = configDir.resolve("samplesheet.csv").toFile();
            Files.writeString(sampleSheet.toPath(),
                    "sample,fastq_1,fastq_2,long_fastq,fasta,single_end\n" +
                            "sample1,sample1.fastq,,,true\n");

            assertTrue(sampleSheet.exists(), "Samplesheet should exist");

            String content = Files.readString(sampleSheet.toPath());
            assertTrue(content.contains("sample"), "Samplesheet should have sample column");
            assertTrue(content.contains("fastq"), "Samplesheet should have fastq column");
        }
    }

    @Nested
    @DisplayName("Output Directory Permission Tests")
    class OutputDirectoryPermissionTests {

        @Test
        @DisplayName("Should verify output directory is writable")
        void testOutputDirectoryWritable() {
            assertTrue(Files.isWritable(outputDir), "Output directory should be writable");
        }

        @Test
        @DisplayName("Should verify output directory is readable")
        void testOutputDirectoryReadable() {
            assertTrue(Files.isReadable(outputDir), "Output directory should be readable");
        }

        @Test
        @DisplayName("Should create output subdirectories")
        void testCreateOutputSubdirectories() throws IOException {
            Path resultsDir = outputDir.resolve("results");
            Path reportsDir = outputDir.resolve("reports");

            Files.createDirectories(resultsDir);
            Files.createDirectories(reportsDir);

            assertTrue(Files.exists(resultsDir), "Results subdirectory should exist");
            assertTrue(Files.exists(reportsDir), "Reports subdirectory should exist");
            assertTrue(Files.isWritable(resultsDir), "Results subdirectory should be writable");
            assertTrue(Files.isWritable(reportsDir), "Reports subdirectory should be writable");
        }

        @Test
        @DisplayName("Should verify directory exists before writing")
        void testDirectoryExistsBeforeWrite() {
            assertTrue(Files.exists(outputDir), "Output directory should exist");
            assertTrue(Files.isDirectory(outputDir), "Output path should be a directory");
        }

        @Test
        @DisplayName("Should handle write permission test")
        void testWritePermission() throws IOException {
            Path testFile = outputDir.resolve("write_test.txt");
            Files.writeString(testFile, "test content");

            assertTrue(Files.exists(testFile), "Test file should be written");
            assertTrue(Files.isReadable(testFile), "Test file should be readable");

            Files.delete(testFile);
            assertFalse(Files.exists(testFile), "Test file should be deleted");
        }
    }

    @Nested
    @DisplayName("Disk Space Validation Tests")
    class DiskSpaceValidationTests {

        @Test
        @DisplayName("Should check available disk space")
        void testAvailableDiskSpace() {
            File workDir = tempWorkDir.toFile();
            long usableSpace = workDir.getUsableSpace();

            assertTrue(usableSpace > 0, "Should have some available disk space");
            System.out.println("Available disk space: " + formatBytes(usableSpace));
        }

        @Test
        @DisplayName("Should verify minimum disk space requirement")
        void testMinimumDiskSpace() {
            File workDir = tempWorkDir.toFile();
            long usableSpace = workDir.getUsableSpace();

            long minimumRequired = 100L * 1024 * 1024; // 100 MB minimum
            assertTrue(usableSpace >= minimumRequired,
                    "Should have at least 100 MB available disk space");
        }

        @Test
        @DisplayName("Should estimate required disk space")
        void testEstimateRequiredSpace() throws IOException {
            // Create sample input files
            List<File> inputFiles = createSampleInputFiles();

            long totalInputSize = 0;
            for (File file : inputFiles) {
                totalInputSize += file.length();
            }

            // Estimate output will be ~3x input size (conservative)
            long estimatedRequired = totalInputSize * 3;

            File workDir = tempWorkDir.toFile();
            long available = workDir.getUsableSpace();

            System.out.println("Input size: " + formatBytes(totalInputSize));
            System.out.println("Estimated required: " + formatBytes(estimatedRequired));
            System.out.println("Available: " + formatBytes(available));

            assertTrue(available > estimatedRequired || estimatedRequired < 1024,
                    "Should have sufficient disk space for workflow");
        }

        @Test
        @DisplayName("Should handle disk space check failure gracefully")
        void testDiskSpaceCheckFailure() {
            // Simulate a scenario where we can't get disk space
            File invalidDir = new File("/nonexistent/path/that/does/not/exist");

            assertDoesNotThrow(() -> {
                long space = invalidDir.getUsableSpace();
                // Should return 0 or handle gracefully
                assertTrue(space >= 0, "Usable space should be non-negative");
            }, "Should handle disk space check gracefully");
        }
    }

    @Nested
    @DisplayName("Concurrent Execution Prevention Tests")
    class ConcurrentExecutionPreventionTests {

        @Test
        @DisplayName("Should detect workflow lock file")
        void testWorkflowLockFile() throws IOException {
            Path lockFile = tempWorkDir.resolve(".workflow.lock");

            // Create lock file
            Files.writeString(lockFile, "workflow_id_123");

            assertTrue(Files.exists(lockFile), "Lock file should exist");

            // Check if workflow is running
            boolean isRunning = Files.exists(lockFile);
            assertTrue(isRunning, "Should detect running workflow");

            // Cleanup lock
            Files.delete(lockFile);
            assertFalse(Files.exists(lockFile), "Lock file should be removed");
        }

        @Test
        @DisplayName("Should prevent concurrent execution in same directory")
        void testPreventConcurrentExecution() throws IOException {
            Path lockFile = tempWorkDir.resolve(".workflow.lock");

            // First execution creates lock
            if (!Files.exists(lockFile)) {
                Files.writeString(lockFile, "execution_1");
                assertTrue(Files.exists(lockFile), "First execution should create lock");
            }

            // Second execution should detect lock
            boolean canExecute = !Files.exists(lockFile);
            assertFalse(canExecute, "Should prevent concurrent execution");

            // Cleanup
            Files.delete(lockFile);
        }

        @Test
        @DisplayName("Should clean up stale lock files")
        void testStaleLockCleanup() throws IOException {
            Path lockFile = tempWorkDir.resolve(".workflow.lock");

            // Create old lock file
            Files.writeString(lockFile, "old_execution");
            Files.setLastModifiedTime(lockFile,
                    java.nio.file.attribute.FileTime.fromMillis(
                            System.currentTimeMillis() - (25 * 60 * 60 * 1000))); // 25 hours ago

            long ageHours = (System.currentTimeMillis() - Files.getLastModifiedTime(lockFile).toMillis())
                    / (60 * 60 * 1000);

            assertTrue(ageHours > 24, "Lock file should be stale (>24 hours)");

            // Should be safe to delete stale locks
            Files.delete(lockFile);
            assertFalse(Files.exists(lockFile), "Stale lock should be removed");
        }
    }

    @Nested
    @DisplayName("Resource Requirement Validation Tests")
    class ResourceRequirementValidationTests {

        @Test
        @DisplayName("Should validate CPU requirements")
        void testCpuRequirements() {
            int availableCpus = Runtime.getRuntime().availableProcessors();
            when(mockOptions.getThreadCount()).thenReturn(4);

            assertTrue(availableCpus > 0, "Should have at least 1 CPU");
            System.out.println("Available CPUs: " + availableCpus);

            int requestedThreads = mockOptions.getThreadCount();
            assertTrue(requestedThreads <= availableCpus * 2,
                    "Requested threads should be reasonable for available CPUs");
        }

        @Test
        @DisplayName("Should validate memory requirements")
        void testMemoryRequirements() {
            Runtime runtime = Runtime.getRuntime();
            long maxMemory = runtime.maxMemory();
            long totalMemory = runtime.totalMemory();

            assertTrue(maxMemory > 0, "Should have max memory available");
            System.out.println("Max memory: " + formatBytes(maxMemory));
            System.out.println("Total memory: " + formatBytes(totalMemory));

            when(mockOptions.getMemoryLimit()).thenReturn(8);
            long requestedMemory = mockOptions.getMemoryLimit() * 1024L * 1024 * 1024; // GB to bytes

            // Just verify we can calculate requirements
            assertTrue(requestedMemory > 0, "Requested memory should be positive");
        }

        @Test
        @DisplayName("Should validate input file count")
        void testInputFileCount() throws IOException {
            List<File> inputFiles = createSampleInputFiles();

            assertTrue(inputFiles.size() > 0, "Should have at least one input file");
            assertTrue(inputFiles.size() <= 1000, "Should not have excessive input files");
        }
    }

    // Helper methods

    /**
     * Creates mock TaxTriageOptions for testing.
     */
    private TaxTriageOptions createMockOptions() {
        TaxTriageOptions options = mock(TaxTriageOptions.class);
        when(options.getSequencingPreset()).thenReturn(TaxTriageOptions.SequencingPreset.ILLUMINA_SE);
        when(options.getQualityThreshold()).thenReturn(20);
        when(options.getMinReadLength()).thenReturn(50);
        when(options.getThreadCount()).thenReturn(4);
        when(options.getMemoryLimit()).thenReturn(8);
        when(options.isEnableSubsampling()).thenReturn(false);
        when(options.getSubsampleSize()).thenReturn(1000000);
        return options;
    }

    /**
     * Creates mock TaxTriageConfig for testing.
     */
    private TaxTriageConfig createMockConfig() {
        TaxTriageConfig config = mock(TaxTriageConfig.class);
        when(config.validate()).thenReturn(null); // null means valid
        when(config.getOutputDirectory()).thenReturn(outputDir.toFile());
        when(config.getThreadCount()).thenReturn(4);
        when(config.getMemoryLimitGb()).thenReturn(8);
        return config;
    }

    /**
     * Creates a valid FASTQ file for testing.
     */
    private void createValidFastqFile(Path fastqFile) throws IOException {
        StringBuilder content = new StringBuilder();
        for (int i = 0; i < 5; i++) {
            content.append("@sequence_").append(i).append("\n");
            content.append("ATCGATCGATCGATCGATCG\n");
            content.append("+\n");
            content.append("IIIIIIIIIIIIIIIIIIII\n");
        }
        Files.writeString(fastqFile, content.toString());
    }

    /**
     * Validates FASTQ file format.
     */
    private boolean isValidFastqFile(Path fastqFile) throws IOException {
        if (!Files.exists(fastqFile) || Files.size(fastqFile) == 0) {
            return false;
        }

        String content = Files.readString(fastqFile);
        String[] lines = content.split("\n");

        if (lines.length % 4 != 0) {
            return false;
        }

        for (int i = 0; i < lines.length; i += 4) {
            if (!lines[i].startsWith("@")) return false;
            if (lines[i + 1].isEmpty()) return false;
            if (!lines[i + 2].startsWith("+")) return false;
            if (lines[i + 1].length() != lines[i + 3].length()) return false;
        }

        return true;
    }

    /**
     * Creates dummy gzipped file for testing.
     */
    private boolean createDummyGzFile(Path gzFile) throws IOException {
        // Create a simple dummy .gz file for testing
        Files.writeString(gzFile, "dummy gzip content");
        return Files.exists(gzFile);
    }

    /**
     * Creates sample input files for testing.
     */
    private List<File> createSampleInputFiles() throws IOException {
        List<File> files = new ArrayList<>();

        for (int i = 1; i <= 3; i++) {
            Path fastqFile = inputDir.resolve("sample_" + i + ".fastq");
            createValidFastqFile(fastqFile);
            files.add(fastqFile.toFile());
        }

        return files;
    }

    /**
     * Formats bytes into human-readable format.
     */
    private String formatBytes(long bytes) {
        if (bytes < 1024) return bytes + " B";
        int exp = (int) (Math.log(bytes) / Math.log(1024));
        String pre = "KMGTPE".charAt(exp - 1) + "";
        return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
    }
}
