package com.jhuapl.taxtriage.geneious.config;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Comprehensive tests for SampleSheetBuilder class.
 *
 * Tests single-end detection, paired-end pairing, ONT handling,
 * CSV formatting, and file validation.
 */
class SampleSheetBuilderTest {

    @TempDir
    Path tempDir;

    private SampleSheetBuilder sampleSheetBuilder;
    private MockFileValidator mockFileValidator;

    @BeforeEach
    void setUp() {
        mockFileValidator = new MockFileValidator();
        sampleSheetBuilder = new SampleSheetBuilder(mockFileValidator);
    }

    @Nested
    @DisplayName("Single-End Detection")
    class SingleEndDetection {

        @Test
        @DisplayName("Should create single-end entries for Illumina SE preset")
        void shouldCreateSingleEndEntriesForIlluminaSe() throws Exception {
            List<File> inputFiles = Arrays.asList(
                    new File("sample1.fastq.gz"),
                    new File("sample2.fastq"),
                    new File("sample3.fq.gz")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE");

            String content = Files.readString(outputFile.toPath());
            String[] lines = content.split("\n");
            assertEquals(4, lines.length); // Header + 3 samples

            assertTrue(lines[0].contains("sample,fastq_1,fastq_2,long_fastq,fasta,single_end"));
            assertTrue(lines[1].contains("sample1,"));
            assertTrue(lines[1].endsWith(",true"));
            assertTrue(lines[2].contains("sample2,"));
            assertTrue(lines[2].endsWith(",true"));
            assertTrue(lines[3].contains("sample3,"));
            assertTrue(lines[3].endsWith(",true"));
        }

        @Test
        @DisplayName("Should extract sample names correctly from SE files")
        void shouldExtractSampleNamesCorrectlyFromSeFiles() throws Exception {
            List<File> inputFiles = Arrays.asList(
                    new File("sample_01_trimmed.fastq.gz"),
                    new File("experimental.condition.1.fq"),
                    new File("data-2023-sample.fastq")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE");

            String content = Files.readString(outputFile.toPath());
            assertTrue(content.contains("sample_01,"));
            assertTrue(content.contains("experimental.condition.1,"));
            assertTrue(content.contains("data-2023-sample,"));
        }
    }

    @Nested
    @DisplayName("Paired-End Pairing")
    class PairedEndPairing {

        @Test
        @DisplayName("Should pair R1 and R2 files correctly")
        void shouldPairR1AndR2FilesCorrectly() throws Exception {
            List<File> inputFiles = Arrays.asList(
                    new File("sample1_R1.fastq.gz"),
                    new File("sample1_R2.fastq.gz"),
                    new File("sample2_R1.fastq"),
                    new File("sample2_R2.fastq")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_PE");

            String content = Files.readString(outputFile.toPath());
            String[] lines = content.split("\n");
            assertEquals(3, lines.length); // Header + 2 samples

            assertTrue(lines[1].contains("sample1,"));
            assertTrue(lines[1].contains("sample1_R1.fastq.gz"));
            assertTrue(lines[1].contains("sample1_R2.fastq.gz"));
            assertTrue(lines[1].endsWith(",false"));

            assertTrue(lines[2].contains("sample2,"));
            assertTrue(lines[2].contains("sample2_R1.fastq"));
            assertTrue(lines[2].contains("sample2_R2.fastq"));
            assertTrue(lines[2].endsWith(",false"));
        }

        @Test
        @DisplayName("Should handle different R1/R2 naming conventions")
        void shouldHandleDifferentR1R2NamingConventions() throws Exception {
            List<File> inputFiles = Arrays.asList(
                    new File("sample_1.fastq.gz"),
                    new File("sample_2.fastq.gz"),
                    new File("data.r1.fq"),
                    new File("data.r2.fq"),
                    new File("test-R1-001.fastq"),
                    new File("test-R2-001.fastq")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_PE");

            String content = Files.readString(outputFile.toPath());
            String[] lines = content.split("\n");
            assertEquals(4, lines.length); // Header + 3 samples

            // Should create 3 paired samples
            assertTrue(content.contains("sample,"));
            assertTrue(content.contains("data,"));
            assertTrue(content.contains("test,"));
        }

        @Test
        @DisplayName("Should throw exception for incomplete pairs")
        void shouldThrowExceptionForIncompletePairs() {
            List<File> inputFiles = Arrays.asList(
                    new File("sample1_R1.fastq.gz"),
                    new File("sample2_R2.fastq.gz") // Missing R1 for sample2
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            assertThrows(SampleSheetBuilder.SampleSheetException.class, () ->
                    sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_PE"));
        }

        @Test
        @DisplayName("Should throw exception for multiple R1 files")
        void shouldThrowExceptionForMultipleR1Files() {
            List<File> inputFiles = Arrays.asList(
                    new File("sample1_R1.fastq.gz"),
                    new File("sample1_R1_backup.fastq.gz"),
                    new File("sample1_R2.fastq.gz")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            assertThrows(SampleSheetBuilder.SampleSheetException.class, () ->
                    sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_PE"));
        }
    }

    @Nested
    @DisplayName("ONT Handling")
    class OntHandling {

        @Test
        @DisplayName("Should handle ONT FASTQ files")
        void shouldHandleOntFastqFiles() throws Exception {
            List<File> inputFiles = Arrays.asList(
                    new File("long_reads_sample1.fastq.gz"),
                    new File("ont_data.fq"),
                    new File("nanopore_run.fastq")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ONT");

            String content = Files.readString(outputFile.toPath());
            String[] lines = content.split("\n");
            assertEquals(4, lines.length); // Header + 3 samples

            assertTrue(lines[1].contains("long_reads_sample1,,,long_reads_sample1.fastq.gz,,false"));
            assertTrue(lines[2].contains("ont_data,,,ont_data.fq,,false"));
            assertTrue(lines[3].contains("nanopore_run,,,nanopore_run.fastq,,false"));
        }

        @Test
        @DisplayName("Should handle ONT FASTA files")
        void shouldHandleOntFastaFiles() throws Exception {
            List<File> inputFiles = Arrays.asList(
                    new File("consensus.fasta"),
                    new File("assembly.fa.gz"),
                    new File("contigs.fna")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ONT");

            String content = Files.readString(outputFile.toPath());
            String[] lines = content.split("\n");
            assertEquals(4, lines.length);

            assertTrue(lines[1].contains("consensus,,,,consensus.fasta,false"));
            assertTrue(lines[2].contains("assembly,,,,assembly.fa.gz,false"));
            assertTrue(lines[3].contains("contigs,,,,contigs.fna,false"));
        }

        @Test
        @DisplayName("Should handle mixed ONT file types")
        void shouldHandleMixedOntFileTypes() throws Exception {
            List<File> inputFiles = Arrays.asList(
                    new File("reads.fastq.gz"),
                    new File("assembly.fasta")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ONT");

            String content = Files.readString(outputFile.toPath());
            assertTrue(content.contains("reads,,,reads.fastq.gz,,false"));
            assertTrue(content.contains("assembly,,,,assembly.fasta,false"));
        }

        @Test
        @DisplayName("Should reject unsupported file types for ONT")
        void shouldRejectUnsupportedFileTypesForOnt() {
            List<File> inputFiles = Arrays.asList(
                    new File("data.txt"),
                    new File("results.xlsx")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            assertThrows(SampleSheetBuilder.SampleSheetException.class, () ->
                    sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ONT"));
        }
    }

    @Nested
    @DisplayName("CSV Formatting")
    class CsvFormatting {

        @Test
        @DisplayName("Should create correct CSV header")
        void shouldCreateCorrectCsvHeader() throws Exception {
            List<File> inputFiles = Arrays.asList(new File("sample.fastq"));
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE");

            String content = Files.readString(outputFile.toPath());
            String[] lines = content.split("\n");
            assertEquals("sample,fastq_1,fastq_2,long_fastq,fasta,single_end", lines[0]);
        }

        @Test
        @DisplayName("Should escape CSV special characters")
        void shouldEscapeCsvSpecialCharacters() throws Exception {
            List<File> inputFiles = Arrays.asList(
                    new File("sample,with,commas.fastq"),
                    new File("sample\"with\"quotes.fastq"),
                    new File("sample\nwith\nnewlines.fastq")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE");

            String content = Files.readString(outputFile.toPath());
            // Should handle CSV escaping properly
            assertTrue(content.contains("\"sample,with,commas\",\"sample,with,commas.fastq\"") ||
                      content.contains("sample_with_commas"));
        }

        @Test
        @DisplayName("Should handle empty fields correctly")
        void shouldHandleEmptyFieldsCorrectly() throws Exception {
            List<File> inputFiles = Arrays.asList(new File("sample.fastq"));
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE");

            String content = Files.readString(outputFile.toPath());
            String[] lines = content.split("\n");
            String[] fields = lines[1].split(",");
            assertEquals(6, fields.length);
            // fastq_2, long_fastq, fasta should be empty
            assertEquals("", fields[2]); // fastq_2
            assertEquals("", fields[3]); // long_fastq
            assertEquals("", fields[4]); // fasta
        }
    }

    @Nested
    @DisplayName("File Validation")
    class FileValidation {

        @Test
        @DisplayName("Should validate all input files exist")
        void shouldValidateAllInputFilesExist() {
            mockFileValidator.setFileExists("existing.fastq", true);
            mockFileValidator.setFileExists("missing.fastq", false);

            List<File> inputFiles = Arrays.asList(
                    new File("existing.fastq"),
                    new File("missing.fastq")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            SampleSheetBuilder.SampleSheetException exception = assertThrows(
                    SampleSheetBuilder.SampleSheetException.class, () ->
                            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE"));

            assertTrue(exception.getMessage().contains("File does not exist"));
            assertTrue(exception.getMessage().contains("missing.fastq"));
        }

        @Test
        @DisplayName("Should validate files are readable")
        void shouldValidateFilesAreReadable() {
            mockFileValidator.setFileExists("readonly.fastq", true);
            mockFileValidator.setFileReadable("readonly.fastq", false);

            List<File> inputFiles = Arrays.asList(new File("readonly.fastq"));
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            SampleSheetBuilder.SampleSheetException exception = assertThrows(
                    SampleSheetBuilder.SampleSheetException.class, () ->
                            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE"));

            assertTrue(exception.getMessage().contains("File is not readable"));
            assertTrue(exception.getMessage().contains("readonly.fastq"));
        }

        @Test
        @DisplayName("Should validate files are not empty")
        void shouldValidateFilesAreNotEmpty() {
            mockFileValidator.setFileExists("empty.fastq", true);
            mockFileValidator.setFileReadable("empty.fastq", true);
            mockFileValidator.setFileEmpty("empty.fastq", true);

            List<File> inputFiles = Arrays.asList(new File("empty.fastq"));
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            SampleSheetBuilder.SampleSheetException exception = assertThrows(
                    SampleSheetBuilder.SampleSheetException.class, () ->
                            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE"));

            assertTrue(exception.getMessage().contains("File is empty"));
            assertTrue(exception.getMessage().contains("empty.fastq"));
        }

        @Test
        @DisplayName("Should collect multiple validation errors")
        void shouldCollectMultipleValidationErrors() {
            mockFileValidator.setFileExists("missing.fastq", false);
            mockFileValidator.setFileExists("readonly.fastq", true);
            mockFileValidator.setFileReadable("readonly.fastq", false);

            List<File> inputFiles = Arrays.asList(
                    new File("missing.fastq"),
                    new File("readonly.fastq")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            SampleSheetBuilder.SampleSheetException exception = assertThrows(
                    SampleSheetBuilder.SampleSheetException.class, () ->
                            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE"));

            assertTrue(exception.getMessage().contains("missing.fastq"));
            assertTrue(exception.getMessage().contains("readonly.fastq"));
        }
    }

    @Nested
    @DisplayName("Generic File Parsing")
    class GenericFileParsing {

        @Test
        @DisplayName("Should handle generic preset with mixed file types")
        void shouldHandleGenericPresetWithMixedFileTypes() throws Exception {
            List<File> inputFiles = Arrays.asList(
                    new File("reads.fastq.gz"),
                    new File("assembly.fasta"),
                    new File("more_reads.fq")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "CUSTOM");

            String content = Files.readString(outputFile.toPath());
            assertTrue(content.contains("reads,reads.fastq.gz,,,,")); // FASTQ as single-end
            assertTrue(content.contains("assembly,,,,assembly.fasta,")); // FASTA
            assertTrue(content.contains("more_reads,more_reads.fq,,,,")); // FASTQ as single-end
        }

        @Test
        @DisplayName("Should reject unsupported file types for generic preset")
        void shouldRejectUnsupportedFileTypesForGenericPreset() {
            List<File> inputFiles = Arrays.asList(new File("data.txt"));
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            assertThrows(SampleSheetBuilder.SampleSheetException.class, () ->
                    sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "CUSTOM"));
        }
    }

    @Nested
    @DisplayName("Error Conditions")
    class ErrorConditions {

        @Test
        @DisplayName("Should throw exception for empty input file list")
        void shouldThrowExceptionForEmptyInputFileList() {
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            assertThrows(SampleSheetBuilder.SampleSheetException.class, () ->
                    sampleSheetBuilder.buildSampleSheet(Arrays.asList(), outputFile, "ILLUMINA_SE"));
        }

        @Test
        @DisplayName("Should throw exception for null input file list")
        void shouldThrowExceptionForNullInputFileList() {
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            assertThrows(SampleSheetBuilder.SampleSheetException.class, () ->
                    sampleSheetBuilder.buildSampleSheet(null, outputFile, "ILLUMINA_SE"));
        }

        @Test
        @DisplayName("Should reject non-FASTQ files for Illumina presets")
        void shouldRejectNonFastqFilesForIlluminaPresets() {
            List<File> inputFiles = Arrays.asList(new File("assembly.fasta"));
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            assertThrows(SampleSheetBuilder.SampleSheetException.class, () ->
                    sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_PE"));

            assertThrows(SampleSheetBuilder.SampleSheetException.class, () ->
                    sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE"));
        }

        @Test
        @DisplayName("Should create parent directories for output file")
        void shouldCreateParentDirectoriesForOutputFile() throws Exception {
            List<File> inputFiles = Arrays.asList(new File("sample.fastq"));
            File outputFile = tempDir.resolve("subdir/samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE");

            assertTrue(outputFile.exists());
            assertTrue(outputFile.getParentFile().exists());
        }
    }

    @Nested
    @DisplayName("Sample Name Extraction")
    class SampleNameExtraction {

        @Test
        @DisplayName("Should remove common suffixes from sample names")
        void shouldRemoveCommonSuffixesFromSampleNames() throws Exception {
            List<File> inputFiles = Arrays.asList(
                    new File("sample_trimmed.fastq.gz"),
                    new File("data_filtered.fq"),
                    new File("reads_clean.fastq")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE");

            String content = Files.readString(outputFile.toPath());
            assertTrue(content.contains("sample,"));
            assertTrue(content.contains("data,"));
            assertTrue(content.contains("reads,"));
        }

        @Test
        @DisplayName("Should replace invalid characters in sample names")
        void shouldReplaceInvalidCharactersInSampleNames() throws Exception {
            List<File> inputFiles = Arrays.asList(
                    new File("sample@#$.fastq"),
                    new File("data%^&.fq")
            );
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE");

            String content = Files.readString(outputFile.toPath());
            // Invalid characters should be replaced with underscores
            assertTrue(content.contains("sample___,") || content.contains("sample,"));
            assertTrue(content.contains("data___,") || content.contains("data,"));
        }

        @Test
        @DisplayName("Should generate unique sample name for empty filename")
        void shouldGenerateUniqueSampleNameForEmptyFilename() throws Exception {
            List<File> inputFiles = Arrays.asList(new File("_.fastq"));
            File outputFile = tempDir.resolve("samplesheet.csv").toFile();

            sampleSheetBuilder.buildSampleSheet(inputFiles, outputFile, "ILLUMINA_SE");

            String content = Files.readString(outputFile.toPath());
            // Should generate a sample name starting with "sample_"
            assertTrue(content.contains("sample_"));
        }
    }

    /**
     * Mock file validator for testing.
     */
    private static class MockFileValidator implements SampleSheetBuilder.FileValidator {
        private final java.util.Map<String, Boolean> existsMap = new java.util.HashMap<>();
        private final java.util.Map<String, Boolean> readableMap = new java.util.HashMap<>();
        private final java.util.Map<String, Boolean> emptyMap = new java.util.HashMap<>();

        public void setFileExists(String filename, boolean exists) {
            existsMap.put(filename, exists);
        }

        public void setFileReadable(String filename, boolean readable) {
            readableMap.put(filename, readable);
        }

        public void setFileEmpty(String filename, boolean empty) {
            emptyMap.put(filename, empty);
        }

        @Override
        public boolean exists(File file) {
            return existsMap.getOrDefault(file.getName(), true);
        }

        @Override
        public boolean isReadable(File file) {
            return readableMap.getOrDefault(file.getName(), true);
        }

        @Override
        public boolean isEmpty(File file) {
            return emptyMap.getOrDefault(file.getName(), false);
        }
    }
}