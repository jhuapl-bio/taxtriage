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
 * Test suite for HtmlReportImporter functionality.
 *
 * Tests the import and parsing of HTML reports from TaxTriage output.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class HtmlReportImporterTest {

    @TempDir
    Path tempDir;

    private HtmlReportImporter importer;

    @BeforeEach
    void setUp() {
        importer = new HtmlReportImporter();
    }

    @Test
    void testImportMultiQCReport() throws IOException {
        File htmlFile = tempDir.resolve("multiqc_report.html").toFile();
        String htmlContent = "<!DOCTYPE html>\n" +
                           "<html>\n" +
                           "<head>\n" +
                           "    <title>MultiQC Report</title>\n" +
                           "    <meta name=\"description\" content=\"Quality control report for sequencing data\">\n" +
                           "</head>\n" +
                           "<body>\n" +
                           "    <h1>MultiQC Report</h1>\n" +
                           "    <p>Sample: sample_001</p>\n" +
                           "    <p>Total reads: 1,000,000</p>\n" +
                           "    <p>Quality score: 35.2</p>\n" +
                           "    <p>Pass rate: 95.5%</p>\n" +
                           "</body>\n" +
                           "</html>";
        Files.writeString(htmlFile.toPath(), htmlContent);

        AnnotatedPluginDocument document = importer.importHtmlReport(htmlFile);

        assertNotNull(document);
        assertEquals("MultiQC Report", document.getName());
        assertEquals("HTML Report", document.getDocument().getDocumentType());

        String description = document.getDocument().getDescription();
        assertTrue(description.contains("MultiQC Report"));
        assertTrue(description.contains("Quality control report for sequencing data"));
        assertTrue(description.contains("Sample: sample_001"));
        assertTrue(description.contains("Quality: 35.2"));
        assertTrue(description.contains("Pass Rate: 95.5%"));

        String content = document.getDocument().toString();
        assertTrue(content.contains("=== HTML Report ==="));
        assertTrue(content.contains("Type: MultiQC Report"));
        assertTrue(content.contains("Title: MultiQC Report"));
        assertTrue(content.contains("Sample Info: Sample: sample_001"));
        assertTrue(content.contains("Statistics: Reads: 1,000,000, Quality: 35.2, Pass Rate: 95.5%"));
    }

    @Test
    void testImportFastPReport() throws IOException {
        File htmlFile = tempDir.resolve("fastp_report.html").toFile();
        String htmlContent = "<!DOCTYPE html>\n" +
                           "<html>\n" +
                           "<head><title>FastP Quality Control Report</title></head>\n" +
                           "<body>\n" +
                           "    <h1>fastp report</h1>\n" +
                           "    <p>Input file: sample1.fastq</p>\n" +
                           "    <p>Total reads: 500,000</p>\n" +
                           "    <p>Quality filtering: 92.3%</p>\n" +
                           "</body>\n" +
                           "</html>";
        Files.writeString(htmlFile.toPath(), htmlContent);

        AnnotatedPluginDocument document = importer.importHtmlReport(htmlFile);

        assertNotNull(document);
        assertEquals("FastP Quality Control Report", document.getName());

        String description = document.getDocument().getDescription();
        assertTrue(description.contains("FastP Quality Report"));
        assertTrue(description.contains("Input: sample1.fastq"));
        assertTrue(description.contains("Reads: 500,000"));
    }

    @Test
    void testImportBasicHtmlReport() throws IOException {
        File htmlFile = tempDir.resolve("basic_report.html").toFile();
        String htmlContent = "<html>\n" +
                           "<head><title>Analysis Report</title></head>\n" +
                           "<body>\n" +
                           "    <h1>Basic Analysis Report</h1>\n" +
                           "    <p>This is a simple report.</p>\n" +
                           "</body>\n" +
                           "</html>";
        Files.writeString(htmlFile.toPath(), htmlContent);

        AnnotatedPluginDocument document = importer.importHtmlReport(htmlFile);

        assertNotNull(document);
        assertEquals("Analysis Report", document.getName());

        String description = document.getDocument().getDescription();
        assertTrue(description.contains("Analysis Report"));
    }

    @Test
    void testImportHtmlWithoutTitle() throws IOException {
        File htmlFile = tempDir.resolve("no_title.html").toFile();
        String htmlContent = "<!DOCTYPE html>\n" +
                           "<html>\n" +
                           "<body>\n" +
                           "    <h1>Report Content</h1>\n" +
                           "    <p>Some analysis results</p>\n" +
                           "</body>\n" +
                           "</html>";
        Files.writeString(htmlFile.toPath(), htmlContent);

        AnnotatedPluginDocument document = importer.importHtmlReport(htmlFile);

        assertNotNull(document);
        assertEquals("no_title Report", document.getName()); // Uses filename
    }

    @Test
    void testImportClassificationReport() throws IOException {
        File htmlFile = tempDir.resolve("classification.html").toFile();
        String htmlContent = "<!DOCTYPE html>\n" +
                           "<html>\n" +
                           "<head><title>Taxonomic Classification</title></head>\n" +
                           "<body>\n" +
                           "    <h1>Kraken Classification Results</h1>\n" +
                           "    <p>Sample classification completed</p>\n" +
                           "    <p>Quality threshold: 30</p>\n" +
                           "</body>\n" +
                           "</html>";
        Files.writeString(htmlFile.toPath(), htmlContent);

        AnnotatedPluginDocument document = importer.importHtmlReport(htmlFile);

        assertNotNull(document);
        String description = document.getDocument().getDescription();
        assertTrue(description.contains("Taxonomic Classification Report"));
        assertTrue(description.contains("Quality: 30"));
    }

    @Test
    void testImportWorkflowReport() throws IOException {
        File htmlFile = tempDir.resolve("workflow.html").toFile();
        String htmlContent = "<!DOCTYPE html>\n" +
                           "<html>\n" +
                           "<head><title>Workflow Execution Report</title></head>\n" +
                           "<body>\n" +
                           "    <h1>Workflow Summary</h1>\n" +
                           "    <p>Pipeline executed successfully</p>\n" +
                           "    <p>Processing time: 2 hours</p>\n" +
                           "</body>\n" +
                           "</html>";
        Files.writeString(htmlFile.toPath(), htmlContent);

        AnnotatedPluginDocument document = importer.importHtmlReport(htmlFile);

        assertNotNull(document);
        String description = document.getDocument().getDescription();
        assertTrue(description.contains("Workflow Report"));
    }

    @Test
    void testIsHtmlFileByExtension() throws IOException {
        File htmlFile1 = tempDir.resolve("test.html").toFile();
        Files.writeString(htmlFile1.toPath(), "<html><body>content</body></html>");
        assertTrue(importer.isHtmlFile(htmlFile1));

        File htmlFile2 = tempDir.resolve("test.htm").toFile();
        Files.writeString(htmlFile2.toPath(), "<html><body>content</body></html>");
        assertTrue(importer.isHtmlFile(htmlFile2));
    }

    @Test
    void testIsHtmlFileByContent() throws IOException {
        File unknownFile = tempDir.resolve("unknown").toFile();
        Files.writeString(unknownFile.toPath(), "<!DOCTYPE html>\n<html><body>content</body></html>");
        assertTrue(importer.isHtmlFile(unknownFile));

        File htmlFile = tempDir.resolve("no_extension").toFile();
        Files.writeString(htmlFile.toPath(), "<html><head><title>Test</title></head><body>content</body></html>");
        assertTrue(importer.isHtmlFile(htmlFile));
    }

    @Test
    void testIsNotHtmlFile() throws IOException {
        File textFile = tempDir.resolve("text.txt").toFile();
        Files.writeString(textFile.toPath(), "This is just plain text");
        assertFalse(importer.isHtmlFile(textFile));

        File xmlFile = tempDir.resolve("data.xml").toFile();
        Files.writeString(xmlFile.toPath(), "<?xml version=\"1.0\"?><root><data>content</data></root>");
        assertFalse(importer.isHtmlFile(xmlFile));
    }

    @Test
    void testImportNullFile() {
        assertThrows(IOException.class, () -> {
            importer.importHtmlReport(null);
        });
    }

    @Test
    void testImportNonExistentFile() {
        File nonExistent = new File("/nonexistent/file.html");
        assertThrows(IOException.class, () -> {
            importer.importHtmlReport(nonExistent);
        });
    }

    @Test
    void testImportEmptyHtmlFile() throws IOException {
        File emptyFile = tempDir.resolve("empty.html").toFile();
        Files.writeString(emptyFile.toPath(), "");

        AnnotatedPluginDocument document = importer.importHtmlReport(emptyFile);

        assertNotNull(document);
        assertEquals("empty Report", document.getName());
    }

    @Test
    void testImportMalformedHtml() throws IOException {
        File malformedFile = tempDir.resolve("malformed.html").toFile();
        String malformedContent = "<html>\n" +
                                 "<head><title>Malformed HTML</head>\n" + // Missing closing tag
                                 "<body>\n" +
                                 "<p>Some content\n" + // Missing closing tag
                                 "<h1>Header</h1>\n" +
                                 "</body>\n";
        Files.writeString(malformedFile.toPath(), malformedContent);

        AnnotatedPluginDocument document = importer.importHtmlReport(malformedFile);

        assertNotNull(document);
        assertEquals("Malformed HTML", document.getName());
    }

    @Test
    void testExtractMetadataWithSpecialCharacters() throws IOException {
        File htmlFile = tempDir.resolve("special.html").toFile();
        String htmlContent = "<!DOCTYPE html>\n" +
                           "<html>\n" +
                           "<head>\n" +
                           "    <title>Análisis de Calidad - Échantillon #1</title>\n" +
                           "    <meta name=\"description\" content=\"Reporte de análisis con caracteres especiales: ñ, á, é\">\n" +
                           "</head>\n" +
                           "<body>\n" +
                           "    <h1>Reporte Especial</h1>\n" +
                           "    <p>Muestra: échantillon_001</p>\n" +
                           "</body>\n" +
                           "</html>";
        Files.writeString(htmlFile.toPath(), htmlContent);

        AnnotatedPluginDocument document = importer.importHtmlReport(htmlFile);

        assertNotNull(document);
        assertEquals("Análisis de Calidad - Échantillon #1", document.getName());

        String description = document.getDocument().getDescription();
        assertTrue(description.contains("caracteres especiales: ñ, á, é"));
        assertTrue(description.contains("Muestra: échantillon_001"));
    }

    @Test
    void testLargeHtmlFile() throws IOException {
        File largeFile = tempDir.resolve("large.html").toFile();
        StringBuilder largeContent = new StringBuilder();
        largeContent.append("<!DOCTYPE html>\n<html>\n<head><title>Large Report</title></head>\n<body>\n");

        // Add lots of content
        for (int i = 0; i < 1000; i++) {
            largeContent.append("<p>Line ").append(i + 1).append(" with content</p>\n");
        }

        largeContent.append("</body>\n</html>");
        Files.writeString(largeFile.toPath(), largeContent.toString());

        AnnotatedPluginDocument document = importer.importHtmlReport(largeFile);

        assertNotNull(document);
        String description = document.getDocument().getDescription();
        assertTrue(description.contains("File size:")); // Should include file size info
    }

    @Test
    void testHtmlWithJavaScript() throws IOException {
        File htmlFile = tempDir.resolve("interactive.html").toFile();
        String htmlContent = "<!DOCTYPE html>\n" +
                           "<html>\n" +
                           "<head>\n" +
                           "    <title>Interactive Report</title>\n" +
                           "    <script>function showData() { alert('Data'); }</script>\n" +
                           "</head>\n" +
                           "<body>\n" +
                           "    <h1>Interactive MultiQC Report</h1>\n" +
                           "    <button onclick=\"showData()\">Show Data</button>\n" +
                           "    <p>Reads: 2,000,000</p>\n" +
                           "</body>\n" +
                           "</html>";
        Files.writeString(htmlFile.toPath(), htmlContent);

        AnnotatedPluginDocument document = importer.importHtmlReport(htmlFile);

        assertNotNull(document);
        String content = document.getDocument().toString();
        assertTrue(content.contains("Note: This is an HTML report file"));
        assertTrue(content.contains("export this document and open it in a web browser"));
        assertTrue(content.contains("Reads: 2,000,000"));
    }

    @Test
    void testStatisticsExtraction() throws IOException {
        File htmlFile = tempDir.resolve("stats.html").toFile();
        String htmlContent = "<!DOCTYPE html>\n" +
                           "<html>\n" +
                           "<body>\n" +
                           "    <h1>Statistics Report</h1>\n" +
                           "    <p>Processed 1,500,000 reads successfully</p>\n" +
                           "    <p>Average quality score: 34.8</p>\n" +
                           "    <p>Pass rate: 89.2%</p>\n" +
                           "    <p>Adapter trimming: 12.5% of reads</p>\n" +
                           "</body>\n" +
                           "</html>";
        Files.writeString(htmlFile.toPath(), htmlContent);

        AnnotatedPluginDocument document = importer.importHtmlReport(htmlFile);

        assertNotNull(document);
        String description = document.getDocument().getDescription();
        assertTrue(description.contains("Reads: 1,500,000"));
        assertTrue(description.contains("Quality: 34.8"));
        assertTrue(description.contains("Pass Rate: 89.2%"));
    }

    @Test
    void testIsHtmlFileEdgeCases() {
        // Test null file
        assertFalse(importer.isHtmlFile(null));

        // Test non-existent file
        File nonExistent = new File("/nonexistent/file.html");
        assertFalse(importer.isHtmlFile(nonExistent));

        // Test directory
        assertFalse(importer.isHtmlFile(tempDir.toFile()));
    }

    @Test
    void testMinimalHtml() throws IOException {
        File minimalFile = tempDir.resolve("minimal.html").toFile();
        Files.writeString(minimalFile.toPath(), "<html><body>Minimal</body></html>");

        AnnotatedPluginDocument document = importer.importHtmlReport(minimalFile);

        assertNotNull(document);
        assertEquals("minimal Report", document.getName());
    }

    @Test
    void testHtmlWithCaseVariations() throws IOException {
        File htmlFile = tempDir.resolve("case_test.html").toFile();
        String htmlContent = "<!doctype html>\n" +
                           "<HTML>\n" +
                           "<HEAD><TITLE>Case Test</TITLE></HEAD>\n" +
                           "<BODY>\n" +
                           "    <H1>MULTIQC REPORT</H1>\n" +
                           "    <P>Quality: 35</P>\n" +
                           "</BODY>\n" +
                           "</HTML>";
        Files.writeString(htmlFile.toPath(), htmlContent);

        AnnotatedPluginDocument document = importer.importHtmlReport(htmlFile);

        assertNotNull(document);
        assertEquals("Case Test", document.getName());

        String description = document.getDocument().getDescription();
        assertTrue(description.contains("MultiQC Report")); // Should detect despite case
        assertTrue(description.contains("Quality: 35"));
    }
}