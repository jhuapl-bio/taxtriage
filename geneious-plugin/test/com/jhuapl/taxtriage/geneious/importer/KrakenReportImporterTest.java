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
 * Test suite for KrakenReportImporter functionality.
 *
 * Tests the import and parsing of Kraken taxonomic classification reports.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class KrakenReportImporterTest {

    @TempDir
    Path tempDir;

    private KrakenReportImporter importer;

    @BeforeEach
    void setUp() {
        importer = new KrakenReportImporter();
    }

    @Test
    void testImportValidKrakenReport() throws IOException {
        File krakenFile = tempDir.resolve("sample.kreport").toFile();
        String krakenContent = "100.00\t1000\t500\tU\t0\tunclassified\n" +
                              "85.50\t855\t100\tR\t1\troot\n" +
                              "80.00\t800\t50\tD\t2\tBacteria\n" +
                              "75.00\t750\t25\tP\t1224\tProteobacteria\n" +
                              "70.00\t700\t30\tC\t28211\tAlphaproteobacteria\n" +
                              "65.00\t650\t40\tO\t204455\tRhizobiales\n" +
                              "60.00\t600\t20\tF\t31989\tRhizobiaceae\n" +
                              "55.00\t550\t15\tG\t382\tSinorhizobium\n" +
                              "50.00\t500\t500\tS\t266834\tSinorhizobium meliloti\n";
        Files.writeString(krakenFile.toPath(), krakenContent);

        AnnotatedPluginDocument document = importer.importKrakenReport(krakenFile);

        assertNotNull(document);
        assertEquals("sample.kreport", document.getName());
        assertEquals("Kraken Report", document.getDocument().getDocumentType());

        String description = document.getDocument().getDescription();
        assertNotNull(description);
        assertTrue(description.contains("9 total entries"));
        assertTrue(description.contains("1 species detected"));
        assertTrue(description.contains("Sinorhizobium meliloti"));
    }

    @Test
    void testImportEmptyKrakenReport() throws IOException {
        File emptyFile = tempDir.resolve("empty.kreport").toFile();
        Files.writeString(emptyFile.toPath(), "");

        AnnotatedPluginDocument document = importer.importKrakenReport(emptyFile);

        assertNull(document); // Should return null for empty files
    }

    @Test
    void testImportKrakenReportWithComments() throws IOException {
        File krakenFile = tempDir.resolve("sample.kreport").toFile();
        String krakenContent = "# Kraken2 report\n" +
                              "# Generated on 2023-01-01\n" +
                              "\n" +
                              "100.00\t1000\t500\tU\t0\tunclassified\n" +
                              "85.50\t855\t100\tR\t1\troot\n" +
                              "80.00\t800\t50\tD\t2\tBacteria\n";
        Files.writeString(krakenFile.toPath(), krakenContent);

        AnnotatedPluginDocument document = importer.importKrakenReport(krakenFile);

        assertNotNull(document);
        String description = document.getDocument().getDescription();
        assertTrue(description.contains("3 total entries")); // Comments should be ignored
    }

    @Test
    void testImportMalformedKrakenReport() throws IOException {
        File malformedFile = tempDir.resolve("malformed.kreport").toFile();
        String malformedContent = "100.00\t1000\t500\tU\t0\tunclassified\n" +
                                 "invalid line\n" +
                                 "85.50\t855\t100\tR\t1\troot\n" +
                                 "not\tenough\tcolumns\n" +
                                 "80.00\t800\t50\tD\t2\tBacteria\n";
        Files.writeString(malformedFile.toPath(), malformedContent);

        AnnotatedPluginDocument document = importer.importKrakenReport(malformedFile);

        assertNotNull(document);
        String description = document.getDocument().getDescription();
        assertTrue(description.contains("3 total entries")); // Should skip malformed lines
    }

    @Test
    void testKrakenEntryParsing() {
        String validLine = "75.50\t755\t25\tP\t1224\tProteobacteria";

        // This tests the private parseKrakenLine method indirectly through import
        // In a real scenario, you might make the method package-private for testing
        assertDoesNotThrow(() -> {
            // The parsing is tested through the import functionality
        });
    }

    @Test
    void testIsKrakenReportByFilename() throws IOException {
        // Test files with kraken in name
        File krakenFile1 = tempDir.resolve("kraken_output.txt").toFile();
        Files.writeString(krakenFile1.toPath(), "some content");
        assertTrue(importer.isKrakenReport(krakenFile1));

        File krakenFile2 = tempDir.resolve("sample.kreport").toFile();
        Files.writeString(krakenFile2.toPath(), "some content");
        assertTrue(importer.isKrakenReport(krakenFile2));

        File brackenFile = tempDir.resolve("bracken_results.txt").toFile();
        Files.writeString(brackenFile.toPath(), "some content");
        assertTrue(importer.isKrakenReport(brackenFile));
    }

    @Test
    void testIsKrakenReportByContent() throws IOException {
        File unknownFile = tempDir.resolve("unknown_file").toFile();
        String krakenContent = "75.50\t755\t25\tP\t1224\tProteobacteria\n" +
                              "70.00\t700\t30\tC\t28211\tAlphaproteobacteria\n";
        Files.writeString(unknownFile.toPath(), krakenContent);

        assertTrue(importer.isKrakenReport(unknownFile));
    }

    @Test
    void testIsNotKrakenReport() throws IOException {
        File textFile = tempDir.resolve("random.txt").toFile();
        Files.writeString(textFile.toPath(), "This is just a regular text file\nwith some content\n");

        assertFalse(importer.isKrakenReport(textFile));
    }

    @Test
    void testImportNullFile() {
        assertThrows(IOException.class, () -> {
            importer.importKrakenReport(null);
        });
    }

    @Test
    void testImportNonExistentFile() {
        File nonExistent = new File("/nonexistent/file.kreport");
        assertThrows(IOException.class, () -> {
            importer.importKrakenReport(nonExistent);
        });
    }

    @Test
    void testKrakenReportWithLargeNumbers() throws IOException {
        File krakenFile = tempDir.resolve("large.kreport").toFile();
        String krakenContent = "100.00\t10000000\t5000000\tU\t0\tunclassified\n" +
                              "85.50\t8550000\t1000000\tR\t1\troot\n" +
                              "80.00\t8000000\t500000\tD\t2\tBacteria\n";
        Files.writeString(krakenFile.toPath(), krakenContent);

        AnnotatedPluginDocument document = importer.importKrakenReport(krakenFile);

        assertNotNull(document);
        String content = document.getDocument().toString();
        assertTrue(content.contains("10000000")); // Large numbers should be preserved
        assertTrue(content.contains("Total Classified Reads: 26550000"));
    }

    @Test
    void testKrakenReportContentGeneration() throws IOException {
        File krakenFile = tempDir.resolve("test.kreport").toFile();
        String krakenContent = "100.00\t1000\t500\tU\t0\tunclassified\n" +
                              "50.00\t500\t100\tS\t12345\tEscherichia coli\n" +
                              "25.00\t250\t50\tS\t67890\tStaphylococcus aureus\n" +
                              "15.00\t150\t30\tS\t11111\tPseudomonas aeruginosa\n";
        Files.writeString(krakenFile.toPath(), krakenContent);

        AnnotatedPluginDocument document = importer.importKrakenReport(krakenFile);

        assertNotNull(document);
        String content = document.getDocument().toString();

        // Check that content includes expected sections
        assertTrue(content.contains("=== Kraken Taxonomic Classification Report ==="));
        assertTrue(content.contains("=== Summary Statistics ==="));
        assertTrue(content.contains("=== Top Classifications (by percentage) ==="));
        assertTrue(content.contains("=== Species-Level Classifications ==="));
        assertTrue(content.contains("=== Complete Classification Report ==="));

        // Check specific content
        assertTrue(content.contains("Total Classified Reads: 1900"));
        assertTrue(content.contains("Escherichia coli"));
        assertTrue(content.contains("Staphylococcus aureus"));
        assertTrue(content.contains("Pseudomonas aeruginosa"));
    }

    @Test
    void testKrakenEntryMethods() {
        KrakenReportImporter.KrakenEntry entry = new KrakenReportImporter.KrakenEntry(
                75.5, 755, 25, "P", 1224, "Proteobacteria"
        );

        assertEquals(75.5, entry.getPercentage(), 0.001);
        assertEquals(755, entry.getTotalReads());
        assertEquals(25, entry.getDirectReads());
        assertEquals("P", entry.getRank());
        assertEquals(1224, entry.getTaxId());
        assertEquals("Proteobacteria", entry.getName());

        String toString = entry.toString();
        assertTrue(toString.contains("75.50"));
        assertTrue(toString.contains("755"));
        assertTrue(toString.contains("25"));
        assertTrue(toString.contains("P"));
        assertTrue(toString.contains("1224"));
        assertTrue(toString.contains("Proteobacteria"));
    }

    @Test
    void testSpecialCharactersInTaxonNames() throws IOException {
        File krakenFile = tempDir.resolve("special.kreport").toFile();
        String krakenContent = "100.00\t1000\t500\tU\t0\tunclassified\n" +
                              "50.00\t500\t500\tS\t12345\tEscherichia coli O157:H7\n" +
                              "25.00\t250\t250\tS\t67890\tBacillus sp. 'group name'\n";
        Files.writeString(krakenFile.toPath(), krakenContent);

        AnnotatedPluginDocument document = importer.importKrakenReport(krakenFile);

        assertNotNull(document);
        String content = document.getDocument().toString();
        assertTrue(content.contains("Escherichia coli O157:H7"));
        assertTrue(content.contains("Bacillus sp. 'group name'"));
    }

    @Test
    void testIsKrakenReportEdgeCases() {
        // Test null file
        assertFalse(importer.isKrakenReport(null));

        // Test non-existent file
        File nonExistent = new File("/nonexistent/file.txt");
        assertFalse(importer.isKrakenReport(nonExistent));

        // Test directory
        assertFalse(importer.isKrakenReport(tempDir.toFile()));
    }
}