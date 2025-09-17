package com.jhuapl.taxtriage.geneious.importer;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Test suite for TextFileImporter functionality.
 *
 * Tests the import and parsing of various text-based files from TaxTriage output.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class TextFileImporterTest {

    @TempDir
    Path tempDir;

    private TextFileImporter importer;

    @BeforeEach
    void setUp() {
        importer = new TextFileImporter();
    }

    @Test
    void testImportCsvFile() throws IOException {
        File csvFile = tempDir.resolve("data.csv").toFile();
        String csvContent = "Sample,Reads,Bases,Quality\n" +
                           "Sample1,1000,50000,35.5\n" +
                           "Sample2,2000,100000,36.2\n" +
                           "Sample3,1500,75000,34.8\n";
        Files.writeString(csvFile.toPath(), csvContent);

        AnnotatedPluginDocument document = importer.importTextFile(csvFile);

        assertNotNull(document);
        assertEquals("data (CSV Table)", document.getName());
        assertEquals("CSV Table", document.getDocument().getDocumentType());

        String description = document.getDocument().getDescription();
        assertTrue(description.contains("4 lines"));
        assertTrue(description.contains("4 columns"));
        assertTrue(description.contains("Contains header row"));

        String content = document.getDocument().toString();
        assertTrue(content.contains("=== CSV Table ==="));
        assertTrue(content.contains("Delimiter: Comma"));
        assertTrue(content.contains("Headers: Sample, Reads, Bases, Quality"));
        assertTrue(content.contains("Sample1\t1000\t50000\t35.5"));
    }

    @Test
    void testImportTsvFile() throws IOException {
        File tsvFile = tempDir.resolve("results.tsv").toFile();
        String tsvContent = "Sample\tReads\tBases\tQuality\n" +
                           "Sample1\t1000\t50000\t35.5\n" +
                           "Sample2\t2000\t100000\t36.2\n";
        Files.writeString(tsvFile.toPath(), tsvContent);

        AnnotatedPluginDocument document = importer.importTextFile(tsvFile);

        assertNotNull(document);
        assertEquals("results (TSV Table)", document.getName());
        assertEquals("TSV Table", document.getDocument().getDocumentType());

        String content = document.getDocument().toString();
        assertTrue(content.contains("Delimiter: Tab"));
        assertTrue(content.contains("Sample\tReads\tBases\tQuality"));
    }

    @Test
    void testImportLogFile() throws IOException {
        File logFile = tempDir.resolve("pipeline.log").toFile();
        String logContent = "INFO: Pipeline started at 2023-01-01 10:00:00\n" +
                           "INFO: Processing sample1\n" +
                           "WARNING: Low quality detected\n" +
                           "INFO: Pipeline completed successfully\n";
        Files.writeString(logFile.toPath(), logContent);

        AnnotatedPluginDocument document = importer.importTextFile(logFile);

        assertNotNull(document);
        assertEquals("pipeline (Log File)", document.getName());
        assertEquals("Log File", document.getDocument().getDocumentType());

        String content = document.getDocument().toString();
        assertTrue(content.contains("=== Log File ==="));
        assertTrue(content.contains("INFO: Pipeline started"));
        assertTrue(content.contains("WARNING: Low quality detected"));
    }

    @Test
    void testImportPlainTextFile() throws IOException {
        File textFile = tempDir.resolve("summary.txt").toFile();
        String textContent = "Analysis Summary Report\n\n" +
                            "Total samples processed: 5\n" +
                            "Total reads: 10,000,000\n" +
                            "Average quality score: 35.2\n" +
                            "\n" +
                            "Quality control passed for all samples.\n";
        Files.writeString(textFile.toPath(), textContent);

        AnnotatedPluginDocument document = importer.importTextFile(textFile);

        assertNotNull(document);
        assertEquals("summary (Text File)", document.getName());
        assertEquals("Text File", document.getDocument().getDocumentType());

        String description = document.getDocument().getDescription();
        assertTrue(description.contains("7 lines"));

        String content = document.getDocument().toString();
        assertTrue(content.contains("=== Text File ==="));
        assertTrue(content.contains("Analysis Summary Report"));
        assertTrue(content.contains("Total samples processed: 5"));
    }

    @Test
    void testImportFileWithoutHeader() throws IOException {
        File csvFile = tempDir.resolve("noheader.csv").toFile();
        String csvContent = "1,1000,50000,35.5\n" +
                           "2,2000,100000,36.2\n" +
                           "3,1500,75000,34.8\n";
        Files.writeString(csvFile.toPath(), csvContent);

        AnnotatedPluginDocument document = importer.importTextFile(csvFile);

        assertNotNull(document);
        String description = document.getDocument().getDescription();
        assertFalse(description.contains("Contains header row"));

        String content = document.getDocument().toString();
        assertFalse(content.contains("Headers:"));
    }

    @Test
    void testImportEmptyFile() throws IOException {
        File emptyFile = tempDir.resolve("empty.txt").toFile();
        Files.writeString(emptyFile.toPath(), "");

        AnnotatedPluginDocument document = importer.importTextFile(emptyFile);

        assertNotNull(document);
        assertEquals("empty (Text File)", document.getName());

        String content = document.getDocument().toString();
        assertTrue(content.contains("Lines: 0"));
    }

    @Test
    void testDetectFileTypes() throws IOException {
        // Test different file types based on name patterns
        File reportFile = tempDir.resolve("analysis_report.txt").toFile();
        Files.writeString(reportFile.toPath(), "Analysis report content");

        AnnotatedPluginDocument document = importer.importTextFile(reportFile);
        assertEquals("analysis_report (Analysis Report)", document.getName());
        assertEquals("Analysis Report", document.getDocument().getDocumentType());

        File statsFile = tempDir.resolve("stats.txt").toFile();
        Files.writeString(statsFile.toPath(), "Statistics content");

        document = importer.importTextFile(statsFile);
        assertEquals("Statistics File", document.getDocument().getDocumentType());

        File classificationFile = tempDir.resolve("classification.txt").toFile();
        Files.writeString(classificationFile.toPath(), "Classification content");

        document = importer.importTextFile(classificationFile);
        assertEquals("Classification Results", document.getDocument().getDocumentType());
    }

    @Test
    void testDelimiterDetection() throws IOException {
        // Test tab delimiter
        File tabFile = tempDir.resolve("tab.txt").toFile();
        Files.writeString(tabFile.toPath(), "col1\tcol2\tcol3\nval1\tval2\tval3\n");

        AnnotatedPluginDocument document = importer.importTextFile(tabFile);
        String content = document.getDocument().toString();
        assertTrue(content.contains("Delimiter: Tab"));
        assertTrue(content.contains("Columns: 3"));

        // Test semicolon delimiter
        File semicolonFile = tempDir.resolve("semicolon.txt").toFile();
        Files.writeString(semicolonFile.toPath(), "col1;col2;col3\nval1;val2;val3\n");

        document = importer.importTextFile(semicolonFile);
        content = document.getDocument().toString();
        assertTrue(content.contains("Columns: 3"));

        // Test pipe delimiter
        File pipeFile = tempDir.resolve("pipe.txt").toFile();
        Files.writeString(pipeFile.toPath(), "col1|col2|col3\nval1|val2|val3\n");

        document = importer.importTextFile(pipeFile);
        content = document.getDocument().toString();
        assertTrue(content.contains("Columns: 3"));
    }

    @Test
    void testHeaderDetection() throws IOException {
        // File with clear header (text columns)
        File headerFile = tempDir.resolve("with_header.csv").toFile();
        Files.writeString(headerFile.toPath(), "Name,Age,Score\nAlice,25,95.5\nBob,30,87.2\n");

        AnnotatedPluginDocument document = importer.importTextFile(headerFile);
        String description = document.getDocument().getDescription();
        assertTrue(description.contains("Contains header row"));

        // File without header (all numeric)
        File noHeaderFile = tempDir.resolve("no_header.csv").toFile();
        Files.writeString(noHeaderFile.toPath(), "1,25,95.5\n2,30,87.2\n3,28,92.1\n");

        document = importer.importTextFile(noHeaderFile);
        description = document.getDocument().getDescription();
        assertFalse(description.contains("Contains header row"));
    }

    @Test
    void testIsTextFile() throws IOException {
        // Test by extension
        File csvFile = tempDir.resolve("test.csv").toFile();
        Files.writeString(csvFile.toPath(), "data");
        assertTrue(importer.isTextFile(csvFile));

        File tsvFile = tempDir.resolve("test.tsv").toFile();
        Files.writeString(tsvFile.toPath(), "data");
        assertTrue(importer.isTextFile(tsvFile));

        File txtFile = tempDir.resolve("test.txt").toFile();
        Files.writeString(txtFile.toPath(), "data");
        assertTrue(importer.isTextFile(txtFile));

        File logFile = tempDir.resolve("test.log").toFile();
        Files.writeString(logFile.toPath(), "data");
        assertTrue(importer.isTextFile(logFile));

        File reportFile = tempDir.resolve("test.report").toFile();
        Files.writeString(reportFile.toPath(), "data");
        assertTrue(importer.isTextFile(reportFile));

        // Test by content
        File unknownFile = tempDir.resolve("unknown").toFile();
        Files.writeString(unknownFile.toPath(), "This is readable text content");
        assertTrue(importer.isTextFile(unknownFile));

        // Test non-text file
        File binaryFile = tempDir.resolve("binary.bin").toFile();
        Files.writeString(binaryFile.toPath(), "\0\1\2\3\4\5");
        assertFalse(importer.isTextFile(binaryFile));
    }

    @Test
    void testImportNullFile() {
        assertThrows(IOException.class, () -> {
            importer.importTextFile(null);
        });
    }

    @Test
    void testImportNonExistentFile() {
        File nonExistent = new File("/nonexistent/file.txt");
        assertThrows(IOException.class, () -> {
            importer.importTextFile(nonExistent);
        });
    }

    @Test
    void testLargeFile() throws IOException {
        File largeFile = tempDir.resolve("large.txt").toFile();
        StringBuilder content = new StringBuilder();
        for (int i = 0; i < 1000; i++) {
            content.append("Line ").append(i + 1).append(" with some content\n");
        }
        Files.writeString(largeFile.toPath(), content.toString());

        AnnotatedPluginDocument document = importer.importTextFile(largeFile);

        assertNotNull(document);
        String description = document.getDocument().getDescription();
        assertTrue(description.contains("1000 lines"));

        String docContent = document.getDocument().toString();
        assertTrue(docContent.contains("Lines: 1000"));
    }

    @Test
    void testFileWithSpecialCharacters() throws IOException {
        File specialFile = tempDir.resolve("special.txt").toFile();
        String specialContent = "Special characters: àáâãäåæçèéêë\n" +
                               "Symbols: !@#$%^&*()_+-=[]{}|;':\",./<>?\n" +
                               "Unicode: 你好世界 🌍 ñoël\n";
        Files.writeString(specialFile.toPath(), specialContent);

        AnnotatedPluginDocument document = importer.importTextFile(specialFile);

        assertNotNull(document);
        String content = document.getDocument().toString();
        assertTrue(content.contains("Special characters: àáâãäåæçèéêë"));
        assertTrue(content.contains("Unicode: 你好世界 🌍 ñoël"));
    }

    @Test
    void testTextFileStructureAnalysis() throws IOException {
        File structuredFile = tempDir.resolve("structured.tsv").toFile();
        String structuredContent = "Column1\tColumn2\tColumn3\tColumn4\n" +
                                  "Value1\tValue2\tValue3\tValue4\n" +
                                  "Data1\tData2\tData3\tData4\n";
        Files.writeString(structuredFile.toPath(), structuredContent);

        AnnotatedPluginDocument document = importer.importTextFile(structuredFile);

        String content = document.getDocument().toString();
        assertTrue(content.contains("Columns: 4"));
        assertTrue(content.contains("Lines: 3"));
        assertTrue(content.contains("Headers: Column1, Column2, Column3, Column4"));
    }

    @Test
    void testIsTextFileEdgeCases() {
        // Test null file
        assertFalse(importer.isTextFile(null));

        // Test non-existent file
        File nonExistent = new File("/nonexistent/file.txt");
        assertFalse(importer.isTextFile(nonExistent));

        // Test directory
        assertFalse(importer.isTextFile(tempDir.toFile()));
    }

    @Test
    void testMixedContentFile() throws IOException {
        File mixedFile = tempDir.resolve("mixed.txt").toFile();
        String mixedContent = "# Comment line\n" +
                             "\n" +
                             "Header1\tHeader2\tHeader3\n" +
                             "Data1\tData2\tData3\n" +
                             "\n" +
                             "# Another comment\n" +
                             "More1\tMore2\tMore3\n";
        Files.writeString(mixedFile.toPath(), mixedContent);

        AnnotatedPluginDocument document = importer.importTextFile(mixedFile);

        assertNotNull(document);
        String content = document.getDocument().toString();
        assertTrue(content.contains("Lines: 7"));
        assertTrue(content.contains("# Comment line"));
        assertTrue(content.contains("Header1\tHeader2\tHeader3"));
    }

    @Test
    void testPreviewLimiting() throws IOException {
        File longFile = tempDir.resolve("long.txt").toFile();
        StringBuilder longContent = new StringBuilder();
        for (int i = 1; i <= 50; i++) {
            longContent.append("Line ").append(i).append("\n");
        }
        Files.writeString(longFile.toPath(), longContent.toString());

        AnnotatedPluginDocument document = importer.importTextFile(longFile);

        String content = document.getDocument().toString();
        // Should show preview of first 10 lines only
        assertTrue(content.contains("Line 1"));
        assertTrue(content.contains("Line 10"));
        assertTrue(content.contains("... (40 more lines)"));
    }
}