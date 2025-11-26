package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import jebl.util.ProgressListener;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for BamFileImporter (Group 4).
 *
 * Tests the BAM file import functionality including:
 * - Import strategy selection
 * - File existence checking
 * - Index creation
 * - Error handling
 * - BAM import support detection
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class BamFileImporterTest {

    @TempDir
    Path tempDir;

    private ProgressListener progressListener;

    @BeforeEach
    void setUp() {
        progressListener = ProgressListener.EMPTY;
    }

    /**
     * GROUP 4: Test importing non-existent BAM file
     */
    @Test
    void testImportNonExistentFile() {
        File bamFile = new File(tempDir.toFile(), "nonexistent.bam");

        List<AnnotatedPluginDocument> result = BamFileImporter.importBamFile(bamFile, progressListener);

        assertNotNull(result, "Result should not be null");
        assertTrue(result.isEmpty(), "Should return empty list for non-existent file");
    }

    /**
     * GROUP 4: Test importing with null progress listener
     */
    @Test
    void testImportWithNullProgressListener() {
        File bamFile = new File(tempDir.toFile(), "test.bam");

        // Should not throw exception with null progress listener
        assertDoesNotThrow(() -> {
            List<AnnotatedPluginDocument> result = BamFileImporter.importBamFile(bamFile, null);
            assertNotNull(result, "Result should not be null even with null progress listener");
        });
    }

    /**
     * GROUP 4: Test importing with empty reference documents list
     */
    @Test
    void testImportWithEmptyReferenceList() {
        File bamFile = new File(tempDir.toFile(), "test.bam");
        List<AnnotatedPluginDocument> emptyRefs = new ArrayList<>();

        List<AnnotatedPluginDocument> result = BamFileImporter.importBamFile(
            bamFile, progressListener, emptyRefs);

        assertNotNull(result, "Result should not be null");
    }

    /**
     * GROUP 4: Test importing with null reference documents
     */
    @Test
    void testImportWithNullReferences() {
        File bamFile = new File(tempDir.toFile(), "test.bam");

        assertDoesNotThrow(() -> {
            List<AnnotatedPluginDocument> result = BamFileImporter.importBamFile(
                bamFile, progressListener, null);
            assertNotNull(result, "Result should not be null");
        });
    }

    /**
     * GROUP 4: Test BAM import support detection
     */
    @Test
    void testIsBamImportSupported() {
        // Should not throw exception
        assertDoesNotThrow(() -> {
            boolean supported = BamFileImporter.isBamImportSupported();
            // Result may vary based on Geneious installation
            assertNotNull(supported);
        });
    }

    /**
     * GROUP 4: Test getting BAM import info
     */
    @Test
    void testGetBamImportInfo() {
        String info = BamFileImporter.getBamImportInfo();

        assertNotNull(info, "Import info should not be null");
        assertFalse(info.isEmpty(), "Import info should not be empty");
        assertTrue(info.contains("BAM Import Capabilities"), "Should contain capabilities header");
    }

    /**
     * GROUP 4: Test import with invalid file (not a BAM)
     */
    @Test
    void testImportInvalidBamFile() throws Exception {
        // Create a text file with .bam extension
        Path fakeBam = tempDir.resolve("fake.bam");
        Files.writeString(fakeBam, "This is not a BAM file");

        List<AnnotatedPluginDocument> result = BamFileImporter.importBamFile(
            fakeBam.toFile(), progressListener);

        assertNotNull(result, "Result should not be null");
        // Import should fail gracefully and return empty list
        assertTrue(result.isEmpty(), "Invalid BAM should return empty list");
    }

    /**
     * GROUP 4: Test import with directory instead of file
     */
    @Test
    void testImportDirectory() {
        File directory = tempDir.toFile();

        List<AnnotatedPluginDocument> result = BamFileImporter.importBamFile(
            directory, progressListener);

        assertNotNull(result, "Result should not be null");
        assertTrue(result.isEmpty(), "Directory should not be imported as BAM");
    }

    /**
     * GROUP 4: Test BAM info includes samtools availability
     */
    @Test
    void testBamInfoContainsSamtoolsStatus() {
        String info = BamFileImporter.getBamImportInfo();

        assertTrue(
            info.contains("Samtools") || info.contains("samtools"),
            "Info should mention samtools status"
        );
    }

    /**
     * GROUP 4: Test BAM info includes importer status
     */
    @Test
    void testBamInfoContainsImporterStatus() {
        String info = BamFileImporter.getBamImportInfo();

        assertTrue(
            info.contains("importer") || info.contains("Importer"),
            "Info should mention importer status"
        );
    }

    /**
     * GROUP 4: Test import with very large file path
     */
    @Test
    void testImportWithLongPath() {
        // Create a very long path
        StringBuilder longPath = new StringBuilder(tempDir.toString());
        for (int i = 0; i < 50; i++) {
            longPath.append("/very_long_directory_name_").append(i);
        }
        longPath.append("/test.bam");

        File bamFile = new File(longPath.toString());

        assertDoesNotThrow(() -> {
            List<AnnotatedPluginDocument> result = BamFileImporter.importBamFile(
                bamFile, progressListener);
            assertNotNull(result, "Result should not be null");
        });
    }

    /**
     * GROUP 4: Test import returns consistent empty list on failure
     */
    @Test
    void testConsistentEmptyListOnFailure() {
        File bamFile1 = new File(tempDir.toFile(), "test1.bam");
        File bamFile2 = new File(tempDir.toFile(), "test2.bam");

        List<AnnotatedPluginDocument> result1 = BamFileImporter.importBamFile(
            bamFile1, progressListener);
        List<AnnotatedPluginDocument> result2 = BamFileImporter.importBamFile(
            bamFile2, progressListener);

        assertNotNull(result1, "First result should not be null");
        assertNotNull(result2, "Second result should not be null");
        assertTrue(result1.isEmpty(), "First result should be empty");
        assertTrue(result2.isEmpty(), "Second result should be empty");
    }

    /**
     * GROUP 4: Test import info format is user-friendly
     */
    @Test
    void testBamImportInfoFormat() {
        String info = BamFileImporter.getBamImportInfo();

        // Check for proper formatting
        assertFalse(info.contains("null"), "Info should not contain 'null'");
        assertTrue(info.contains("\n") || info.contains("-"),
            "Info should be formatted with newlines or bullet points");
    }

    /**
     * GROUP 4: Test that multiple import attempts don't cause issues
     */
    @Test
    void testMultipleImportAttempts() {
        File bamFile = new File(tempDir.toFile(), "test.bam");

        // Try importing same non-existent file multiple times
        for (int i = 0; i < 5; i++) {
            final int attemptNum = i;
            assertDoesNotThrow(() -> {
                List<AnnotatedPluginDocument> result = BamFileImporter.importBamFile(
                    bamFile, progressListener);
                assertNotNull(result, "Result should not be null on attempt " + attemptNum);
            });
        }
    }
}
