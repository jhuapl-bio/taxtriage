package com.jhuapl.taxtriage.geneious.importer;

import com.biomatters.geneious.publicapi.databaseservice.DatabaseServiceException;
import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import jebl.util.ProgressListener;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;

/**
 * Comprehensive test suite for ResultImporter functionality.
 *
 * Tests the import of various TaxTriage output file types, folder structure
 * creation, error handling, and progress tracking.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class ResultImporterTest {

    @Mock
    private WritableDatabaseService mockDatabaseService;

    @Mock
    private ProgressListener mockProgressListener;

    @TempDir
    Path tempDir;

    private ResultImporter resultImporter;
    private File outputDirectory;

    @BeforeEach
    void setUp() throws IOException {
        MockitoAnnotations.openMocks(this);
        resultImporter = new ResultImporter(mockDatabaseService);
        outputDirectory = tempDir.resolve("taxtriage_output").toFile();
        outputDirectory.mkdirs();
    }

    @Test
    void testImportEmptyDirectory() throws DatabaseServiceException {
        ImportResult result = resultImporter.importResults(outputDirectory, "test-workflow", mockProgressListener);

        assertNotNull(result);
        assertTrue(result.isCompleted());
        assertEquals(0, result.getTotalFilesProcessed());
        assertEquals(0, result.getSuccessfulImports());
        assertEquals(0, result.getFailedImports());
    }

    @Test
    void testImportWithValidFiles() throws IOException, DatabaseServiceException {
        // Create test files
        createTestFiles();

        ImportResult result = resultImporter.importResults(outputDirectory, "test-workflow", mockProgressListener);

        assertNotNull(result);
        assertTrue(result.isCompleted());
        assertTrue(result.getTotalFilesProcessed() > 0);
        assertTrue(result.getSuccessfulImports() > 0);
    }

    @Test
    void testImportKrakenReport() throws IOException, DatabaseServiceException {
        // Create a test Kraken report
        File krakenDir = new File(outputDirectory, "kraken2");
        krakenDir.mkdirs();
        File krakenReport = new File(krakenDir, "sample1.kreport");

        String krakenContent = "100.00\t1000\t500\tU\t0\tunclassified\n" +
                              "85.50\t855\t100\tR\t1\troot\n" +
                              "80.00\t800\t50\tD\t2\tBacteria\n" +
                              "75.00\t750\t25\tP\t1224\tProteobacteria\n" +
                              "70.00\t700\t30\tC\t28211\tAlphaproteobacteria\n";
        Files.writeString(krakenReport.toPath(), krakenContent);

        ImportResult result = resultImporter.importResults(outputDirectory, "test-workflow", mockProgressListener);

        assertNotNull(result);
        assertTrue(result.isCompleted());
        assertTrue(result.getSuccessfulImports() > 0);
        assertTrue(result.getFileTypeStats().containsKey("Kraken Report"));
    }

    @Test
    void testImportFastqFile() throws IOException, DatabaseServiceException {
        // Create a test FASTQ file
        File fastpDir = new File(outputDirectory, "fastp");
        fastpDir.mkdirs();
        File fastqFile = new File(fastpDir, "sample1_cleaned.fastq");

        String fastqContent = "@read1\n" +
                             "ATCGATCGATCGATCG\n" +
                             "+\n" +
                             "IIIIIIIIIIIIIIII\n" +
                             "@read2\n" +
                             "GCTAGCTAGCTAGCTA\n" +
                             "+\n" +
                             "JJJJJJJJJJJJJJJJ\n";
        Files.writeString(fastqFile.toPath(), fastqContent);

        ImportResult result = resultImporter.importResults(outputDirectory, "test-workflow", mockProgressListener);

        assertNotNull(result);
        assertTrue(result.isCompleted());
        assertTrue(result.getSuccessfulImports() > 0);
        assertTrue(result.getFileTypeStats().containsValue(1));
    }

    @Test
    void testImportHtmlReport() throws IOException, DatabaseServiceException {
        // Create a test HTML report
        File reportsDir = new File(outputDirectory, "multiqc");
        reportsDir.mkdirs();
        File htmlReport = new File(reportsDir, "multiqc_report.html");

        String htmlContent = "<!DOCTYPE html>\n" +
                           "<html>\n" +
                           "<head><title>MultiQC Report</title></head>\n" +
                           "<body>\n" +
                           "<h1>Quality Control Report</h1>\n" +
                           "<p>Sample: test_sample</p>\n" +
                           "<p>Reads: 1,000,000</p>\n" +
                           "</body>\n" +
                           "</html>";
        Files.writeString(htmlReport.toPath(), htmlContent);

        ImportResult result = resultImporter.importResults(outputDirectory, "test-workflow", mockProgressListener);

        assertNotNull(result);
        assertTrue(result.isCompleted());
        assertTrue(result.getSuccessfulImports() > 0);
        assertTrue(result.getFileTypeStats().containsKey("HTML Report"));
    }

    @Test
    void testImportTextFiles() throws IOException, DatabaseServiceException {
        // Create various text files
        File countDir = new File(outputDirectory, "count");
        countDir.mkdirs();

        // CSV file
        File csvFile = new File(countDir, "counts.csv");
        Files.writeString(csvFile.toPath(), "Sample,Reads,Bases\nSample1,1000,50000\nSample2,2000,100000\n");

        // TSV file
        File tsvFile = new File(countDir, "summary.tsv");
        Files.writeString(tsvFile.toPath(), "Sample\tReads\tBases\nSample1\t1000\t50000\nSample2\t2000\t100000\n");

        // Log file
        File logFile = new File(outputDirectory, "pipeline.log");
        Files.writeString(logFile.toPath(), "INFO: Pipeline started\nINFO: Processing sample1\nINFO: Pipeline completed\n");

        ImportResult result = resultImporter.importResults(outputDirectory, "test-workflow", mockProgressListener);

        assertNotNull(result);
        assertTrue(result.isCompleted());
        assertTrue(result.getSuccessfulImports() >= 3);
        assertTrue(result.getFileTypeStats().size() > 0);
    }

    @Test
    void testSkipUnsupportedFiles() throws IOException, DatabaseServiceException {
        // Create files that should be skipped
        File skipFile1 = new File(outputDirectory, ".DS_Store");
        Files.writeString(skipFile1.toPath(), "binary content");

        File skipFile2 = new File(outputDirectory, "temp.tmp");
        Files.writeString(skipFile2.toPath(), "temporary file");

        File workDir = new File(outputDirectory, "work");
        workDir.mkdirs();
        File workFile = new File(workDir, "task.log");
        Files.writeString(workFile.toPath(), "work file");

        ImportResult result = resultImporter.importResults(outputDirectory, "test-workflow", mockProgressListener);

        assertNotNull(result);
        assertTrue(result.isCompleted());
        assertEquals(0, result.getTotalFilesProcessed()); // All files should be skipped
    }

    @Test
    void testFolderStructureCreation() throws IOException, DatabaseServiceException {
        // Create nested directory structure
        File level1 = new File(outputDirectory, "fastp");
        File level2 = new File(level1, "sample1");
        level2.mkdirs();

        File testFile = new File(level2, "test.txt");
        Files.writeString(testFile.toPath(), "test content");

        ImportResult result = resultImporter.importResults(outputDirectory, "test-workflow", mockProgressListener);

        assertNotNull(result);
        assertTrue(result.isCompleted());
        assertTrue(result.getCreatedFolders().size() > 0);

        FolderStructureBuilder folderBuilder = resultImporter.getFolderBuilder();
        assertTrue(folderBuilder.getFolderCount() > 0);
    }

    @Test
    void testProgressTracking() throws IOException, DatabaseServiceException {
        createTestFiles();

        resultImporter.importResults(outputDirectory, "test-workflow", mockProgressListener);

        // Verify progress listener was called
        verify(mockProgressListener, atLeastOnce()).setMessage(anyString());
        verify(mockProgressListener, atLeastOnce()).setProgress(anyDouble());
    }

    @Test
    void testErrorHandling() throws IOException, DatabaseServiceException {
        // Create a file that will cause an error during import
        File problemFile = new File(outputDirectory, "problem.txt");
        Files.writeString(problemFile.toPath(), "content");

        // Make the file unreadable (this might not work on all systems)
        problemFile.setReadable(false);

        ImportResult result = resultImporter.importResults(outputDirectory, "test-workflow", mockProgressListener);

        assertNotNull(result);
        assertTrue(result.isCompleted());
        // Should still complete even with errors
    }

    @Test
    void testInvalidDirectory() {
        File invalidDir = new File("/nonexistent/directory");

        assertThrows(IllegalArgumentException.class, () -> {
            resultImporter.importResults(invalidDir, "test-workflow", mockProgressListener);
        });
    }

    @Test
    void testImportStatistics() throws IOException, DatabaseServiceException {
        createTestFiles();

        ImportResult result = resultImporter.importResults(outputDirectory, "test-workflow", mockProgressListener);

        String statistics = resultImporter.getImportStatistics();
        assertNotNull(statistics);
        assertTrue(statistics.contains("Import Statistics"));
        assertTrue(statistics.contains("Total Files Processed"));
        assertTrue(statistics.contains("Success Rate"));
    }

    @Test
    void testValidateConfiguration() {
        assertTrue(resultImporter.validateConfiguration());
    }

    @Test
    void testReset() throws IOException, DatabaseServiceException {
        createTestFiles();
        resultImporter.importResults(outputDirectory, "test-workflow", mockProgressListener);

        resultImporter.reset();

        // After reset, folder builder should be reset
        assertEquals(0, resultImporter.getFolderBuilder().getFolderCount());
    }

    /**
     * Creates a variety of test files for comprehensive testing.
     */
    private void createTestFiles() throws IOException {
        // Create directory structure
        File fastpDir = new File(outputDirectory, "fastp");
        File krakenDir = new File(outputDirectory, "kraken2");
        File countDir = new File(outputDirectory, "count");
        File reportsDir = new File(outputDirectory, "multiqc");

        fastpDir.mkdirs();
        krakenDir.mkdirs();
        countDir.mkdirs();
        reportsDir.mkdirs();

        // FASTQ file
        File fastqFile = new File(fastpDir, "sample1.fastq");
        String fastqContent = "@read1\nATCGATCG\n+\nIIIIIIII\n@read2\nGCTAGCTA\n+\nJJJJJJJJ\n";
        Files.writeString(fastqFile.toPath(), fastqContent);

        // Kraken report
        File krakenReport = new File(krakenDir, "sample1.kreport");
        String krakenContent = "100.00\t1000\t500\tU\t0\tunclassified\n85.50\t855\t100\tR\t1\troot\n";
        Files.writeString(krakenReport.toPath(), krakenContent);

        // CSV file
        File csvFile = new File(countDir, "counts.csv");
        Files.writeString(csvFile.toPath(), "Sample,Reads\nSample1,1000\nSample2,2000\n");

        // HTML report
        File htmlFile = new File(reportsDir, "report.html");
        String htmlContent = "<html><head><title>Report</title></head><body><h1>Test Report</h1></body></html>";
        Files.writeString(htmlFile.toPath(), htmlContent);

        // Text file
        File textFile = new File(outputDirectory, "summary.txt");
        Files.writeString(textFile.toPath(), "Analysis Summary\nTotal samples: 2\nTotal reads: 3000\n");
    }
}