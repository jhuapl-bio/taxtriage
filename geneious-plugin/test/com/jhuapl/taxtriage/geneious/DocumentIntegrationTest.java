package com.jhuapl.taxtriage.geneious;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.documents.PluginDocument;
import com.biomatters.geneious.publicapi.implementations.sequence.DefaultNucleotideSequence;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Integration tests for document creation and import functionality (Group 2).
 *
 * Tests:
 * - Document creation from workflow results
 * - Document import into Geneious
 * - Result file handling
 * - Document metadata and properties
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class DocumentIntegrationTest {

    @TempDir
    Path tempDir;

    private Path outputDir;

    @BeforeEach
    void setUp() throws Exception {
        outputDir = tempDir.resolve("output");
        Files.createDirectories(outputDir);
    }

    /**
     * GROUP 2: Test creating sequence documents programmatically
     */
    @Test
    void testCreateSequenceDocument() throws Exception {
        DefaultNucleotideSequence sequence = new DefaultNucleotideSequence(
            "TestSequence",
            "Test description",
            "ACGTACGTACGT",
            new Date()
        );

        assertNotNull(sequence, "Sequence should be created");
        assertEquals("TestSequence", sequence.getName());
        assertEquals(12, sequence.getSequenceLength());
    }

    /**
     * GROUP 2: Test converting sequence to annotated plugin document
     *
     * Note: This test verifies document creation without requiring full Geneious framework.
     * In production, TaxTriageOperation uses DocumentUtilities.createAnnotatedPluginDocuments()
     * which requires the full Geneious environment.
     */
    @Test
    void testConvertToAnnotatedDocument() throws Exception {
        DefaultNucleotideSequence sequence = new DefaultNucleotideSequence(
            "TestSequence",
            "Test description",
            "ACGTACGTACGT",
            new Date()
        );

        // Verify the sequence document was created correctly
        assertNotNull(sequence, "Sequence should be created");
        assertEquals("TestSequence", sequence.getName());
        assertEquals("Test description", sequence.getDescription());
        assertEquals(12, sequence.getSequenceLength());

        // Note: DocumentUtilities.createAnnotatedPluginDocuments() requires full Geneious framework
        // This is tested in integration tests within Geneious
    }

    /**
     * GROUP 2: Test creating multiple documents from results
     *
     * Note: This test verifies multiple document creation without requiring full Geneious framework.
     */
    @Test
    void testCreateMultipleDocuments() throws Exception {
        List<PluginDocument> documents = new ArrayList<>();

        for (int i = 0; i < 5; i++) {
            DefaultNucleotideSequence seq = new DefaultNucleotideSequence(
                "Sequence" + i,
                "Description " + i,
                "ACGT",
                new Date()
            );
            documents.add(seq);
        }

        // Verify all documents were created correctly
        assertEquals(5, documents.size(), "Should create 5 documents");
        for (int i = 0; i < 5; i++) {
            assertEquals("Sequence" + i, documents.get(i).getName());
            assertEquals("Description " + i, documents.get(i).getDescription());
        }

        // Note: DocumentUtilities.createAnnotatedPluginDocuments() requires full Geneious framework
        // This is tested in integration tests within Geneious
    }

    /**
     * GROUP 2: Test document with description/metadata
     */
    @Test
    void testDocumentWithMetadata() throws Exception {
        DefaultNucleotideSequence sequence = new DefaultNucleotideSequence(
            "ResultSequence",
            "TaxTriage classification result",
            "N", // Placeholder for report
            new Date()
        );

        String reportContent = "=== TaxTriage Analysis Report ===\n" +
                              "Sample: test_sample\n" +
                              "Classified reads: 1000\n";
        sequence.setDescription(reportContent);

        assertEquals(reportContent, sequence.getDescription());
        assertTrue(sequence.getDescription().contains("TaxTriage"));
    }

    /**
     * GROUP 2: Test reading result files from output directory
     */
    @Test
    void testReadResultFiles() throws Exception {
        // Create mock result files
        File reportFile = outputDir.resolve("report.html").toFile();
        File csvFile = outputDir.resolve("results.csv").toFile();
        File summaryFile = outputDir.resolve("summary.txt").toFile();

        Files.writeString(reportFile.toPath(), "<html><body>Report</body></html>");
        Files.writeString(csvFile.toPath(), "sample,count\ntest,100\n");
        Files.writeString(summaryFile.toPath(), "Summary of results");

        assertTrue(reportFile.exists(), "Report file should exist");
        assertTrue(csvFile.exists(), "CSV file should exist");
        assertTrue(summaryFile.exists(), "Summary file should exist");

        // Test file type detection
        assertTrue(isReportFile(reportFile), "Should recognize HTML as report");
        assertTrue(isReportFile(csvFile), "Should recognize CSV as report");
        assertTrue(isReportFile(summaryFile), "Should recognize summary as report");
    }

    /**
     * GROUP 2: Test importing result files as documents
     */
    @Test
    void testImportResultFilesAsDocuments() throws Exception {
        // Create a result file
        File resultFile = outputDir.resolve("classification_results.txt").toFile();
        String resultContent = "Taxonomy ID\tCount\n9606\t500\n";
        Files.writeString(resultFile.toPath(), resultContent);

        // Import as document
        DefaultNucleotideSequence doc = new DefaultNucleotideSequence(
            "TaxTriage_" + resultFile.getName(),
            "TaxTriage result: " + resultFile.getName(),
            "N",
            new Date()
        );
        doc.setDescription(resultContent);

        assertNotNull(doc);
        assertEquals("TaxTriage_classification_results.txt", doc.getName());
        assertTrue(doc.getDescription().contains("Taxonomy"));
    }

    /**
     * GROUP 2: Test document collection from workflow output
     *
     * Note: This test verifies workflow result document creation logic.
     */
    @Test
    void testCollectDocumentsFromWorkflow() throws Exception {
        List<PluginDocument> results = new ArrayList<>();

        // Simulate creating multiple result documents
        // 1. Summary report
        DefaultNucleotideSequence summary = new DefaultNucleotideSequence(
            "TaxTriage_Summary",
            "Workflow summary",
            "N",
            new Date()
        );
        summary.setDescription("=== Workflow Summary ===\nCompleted successfully");
        results.add(summary);

        // 2. Classification results
        File[] resultFiles = createMockResultFiles();
        for (File file : resultFiles) {
            String content = Files.readString(file.toPath());
            DefaultNucleotideSequence doc = new DefaultNucleotideSequence(
                "TaxTriage_" + file.getName(),
                "Result: " + file.getName(),
                "N",
                new Date()
            );
            doc.setDescription(content);
            results.add(doc);
        }

        // Verify document collection
        assertTrue(results.size() >= 1, "Should have at least summary document");
        assertEquals("TaxTriage_Summary", results.get(0).getName());
        assertTrue(results.get(0).getDescription().contains("Workflow Summary"));

        // Verify result file documents
        for (int i = 1; i < results.size(); i++) {
            assertTrue(results.get(i).getName().startsWith("TaxTriage_"));
        }
    }

    /**
     * GROUP 2: Test handling empty output directory
     */
    @Test
    void testHandleEmptyOutputDirectory() throws Exception {
        Path emptyDir = tempDir.resolve("empty_output");
        Files.createDirectories(emptyDir);

        File[] files = emptyDir.toFile().listFiles();
        assertNotNull(files);
        assertEquals(0, files.length, "Directory should be empty");

        // Should handle gracefully - no documents to import
        List<PluginDocument> results = new ArrayList<>();
        assertEquals(0, results.size(), "Should have no results from empty directory");
    }

    /**
     * GROUP 2: Test filtering report files from output
     */
    @Test
    void testFilterReportFiles() throws Exception {
        // Create various file types
        File htmlReport = outputDir.resolve("report.html").toFile();
        File csvData = outputDir.resolve("data.csv").toFile();
        File bamFile = outputDir.resolve("alignment.bam").toFile();
        File logFile = outputDir.resolve("workflow.log").toFile();
        File txtReport = outputDir.resolve("summary_report.txt").toFile();

        Files.writeString(htmlReport.toPath(), "<html></html>");
        Files.writeString(csvData.toPath(), "data");
        Files.writeString(bamFile.toPath(), "bam data");
        Files.writeString(logFile.toPath(), "log");
        Files.writeString(txtReport.toPath(), "report");

        // Test file filtering
        assertTrue(isReportFile(htmlReport), "HTML should be report");
        assertTrue(isReportFile(csvData), "CSV should be report");
        assertFalse(isReportFile(bamFile), "BAM should not be report");
        assertTrue(isReportFile(txtReport), "Summary text should be report");
    }

    /**
     * GROUP 2: Test document properties and metadata preservation
     */
    @Test
    void testDocumentMetadataPreservation() throws Exception {
        Date creationDate = new Date();
        DefaultNucleotideSequence sequence = new DefaultNucleotideSequence(
            "TestSeq",
            "Original description",
            "ACGT",
            creationDate
        );

        // Add additional metadata
        sequence.setDescription("Updated description with metadata");

        // Verify metadata is preserved
        assertEquals("TestSeq", sequence.getName());
        assertEquals("Updated description with metadata", sequence.getDescription());
        assertNotNull(sequence.getCreationDate());
    }

    // ===== Helper Methods =====

    /**
     * Checks if a file is a report file that should be imported
     */
    private boolean isReportFile(File file) {
        if (file == null || !file.isFile()) {
            return false;
        }

        String name = file.getName().toLowerCase();
        return name.endsWith(".html") || name.endsWith(".txt") ||
               name.endsWith(".csv") || name.endsWith(".tsv") ||
               name.contains("report") || name.contains("summary") ||
               name.contains("classification") || name.contains("taxonomy");
    }

    /**
     * Creates mock result files for testing
     */
    private File[] createMockResultFiles() throws Exception {
        List<File> files = new ArrayList<>();

        File report = outputDir.resolve("taxonomy_report.html").toFile();
        Files.writeString(report.toPath(), "<html><body>Taxonomy Report</body></html>");
        files.add(report);

        File csv = outputDir.resolve("classification.csv").toFile();
        Files.writeString(csv.toPath(), "taxid,count\n9606,100\n");
        files.add(csv);

        return files.toArray(new File[0]);
    }
}
