package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for SimpleFastaImporter (Group 4).
 *
 * Tests the FASTA file import functionality including:
 * - Importing valid FASTA files
 * - Handling invalid or malformed FASTA
 * - Empty file handling
 * - Multiple sequence import
 * - Edge cases (whitespace, empty lines, etc.)
 *
 * NOTE: These tests require full Geneious API environment with DocumentUtilities
 * initialized. They are disabled for headless unit testing but can be enabled
 * in integration test environment.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
@Disabled("Requires full Geneious API environment - enable for integration testing")
class SimpleFastaImporterTest {

    @TempDir
    Path tempDir;

    /**
     * GROUP 4: Test importing valid single sequence FASTA
     */
    @Test
    void testImportSingleSequence() throws Exception {
        // Create test FASTA file
        Path fastaFile = tempDir.resolve("test.fasta");
        String fastaContent = ">sequence1\nATGCATGCATGC\n";
        Files.writeString(fastaFile, fastaContent);

        // Import
        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile.toFile());

        assertNotNull(docs, "Result should not be null");
        assertEquals(1, docs.size(), "Should import one sequence");
        assertEquals("sequence1", docs.get(0).getName(), "Sequence name should match");
    }

    /**
     * GROUP 4: Test importing multiple sequences
     */
    @Test
    void testImportMultipleSequences() throws Exception {
        Path fastaFile = tempDir.resolve("multi.fasta");
        String fastaContent = ">seq1\n" +
            "ATGC\n" +
            ">seq2\n" +
            "GCTA\n" +
            ">seq3\n" +
            "TTAA\n";
        Files.writeString(fastaFile, fastaContent);

        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile.toFile());

        assertEquals(3, docs.size(), "Should import three sequences");
        assertEquals("seq1", docs.get(0).getName());
        assertEquals("seq2", docs.get(1).getName());
        assertEquals("seq3", docs.get(2).getName());
    }

    /**
     * GROUP 4: Test importing sequences with wrapped lines
     */
    @Test
    void testImportWrappedSequence() throws Exception {
        Path fastaFile = tempDir.resolve("wrapped.fasta");
        String fastaContent = ">wrapped_sequence\n" +
            "ATGCATGCATGC\n" +
            "GCTAGCTAGCTA\n" +
            "TTAATTAATTAA\n";
        Files.writeString(fastaFile, fastaContent);

        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile.toFile());

        assertEquals(1, docs.size(), "Should import one sequence");
        assertEquals("wrapped_sequence", docs.get(0).getName());
    }

    /**
     * GROUP 4: Test empty file handling
     */
    @Test
    void testImportEmptyFile() throws Exception {
        Path fastaFile = tempDir.resolve("empty.fasta");
        Files.writeString(fastaFile, "");

        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile.toFile());

        assertNotNull(docs, "Result should not be null");
        assertEquals(0, docs.size(), "Empty file should produce no sequences");
    }

    /**
     * GROUP 4: Test file with only whitespace
     */
    @Test
    void testImportWhitespaceOnlyFile() throws Exception {
        Path fastaFile = tempDir.resolve("whitespace.fasta");
        Files.writeString(fastaFile, "   \n\n  \n");

        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile.toFile());

        assertNotNull(docs, "Result should not be null");
        assertEquals(0, docs.size(), "Whitespace-only file should produce no sequences");
    }

    /**
     * GROUP 4: Test non-existent file handling
     */
    @Test
    void testImportNonExistentFile() {
        File nonExistent = new File(tempDir.toFile(), "does_not_exist.fasta");

        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(nonExistent);

        assertNotNull(docs, "Result should not be null");
        assertEquals(0, docs.size(), "Non-existent file should produce empty list");
    }

    /**
     * GROUP 4: Test FASTA with empty lines between sequences
     */
    @Test
    void testImportWithEmptyLines() throws Exception {
        Path fastaFile = tempDir.resolve("empty_lines.fasta");
        String fastaContent = ">seq1\n" +
            "ATGC\n" +
            "\n" +
            ">seq2\n" +
            "\n" +
            "GCTA\n" +
            "\n" +
            ">seq3\n" +
            "TTAA\n";
        Files.writeString(fastaFile, fastaContent);

        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile.toFile());

        assertEquals(3, docs.size(), "Should handle empty lines correctly");
    }

    /**
     * GROUP 4: Test FASTA with header descriptions
     */
    @Test
    void testImportWithHeaderDescriptions() throws Exception {
        Path fastaFile = tempDir.resolve("descriptions.fasta");
        String fastaContent = ">seq1 description with spaces\n" +
            "ATGC\n" +
            ">seq2|accession|description\n" +
            "GCTA\n";
        Files.writeString(fastaFile, fastaContent);

        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile.toFile());

        assertEquals(2, docs.size(), "Should import both sequences");
        assertEquals("seq1 description with spaces", docs.get(0).getName());
        assertEquals("seq2|accession|description", docs.get(1).getName());
    }

    /**
     * GROUP 4: Test sequence with only header, no sequence data
     */
    @Test
    void testImportHeaderOnly() throws Exception {
        Path fastaFile = tempDir.resolve("header_only.fasta");
        String fastaContent = ">seq1\n";
        Files.writeString(fastaFile, fastaContent);

        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile.toFile());

        // Header with no sequence data should not create a document
        assertEquals(0, docs.size(), "Header without sequence should not be imported");
    }

    /**
     * GROUP 4: Test malformed FASTA (sequence before header)
     */
    @Test
    void testImportMalformedFasta() throws Exception {
        Path fastaFile = tempDir.resolve("malformed.fasta");
        String fastaContent = "ATGCATGC\n" +
            ">seq1\n" +
            "GCTAGCTA\n";
        Files.writeString(fastaFile, fastaContent);

        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile.toFile());

        // Should skip lines before first header
        assertEquals(1, docs.size(), "Should import sequence after header");
        assertEquals("seq1", docs.get(0).getName());
    }

    /**
     * GROUP 4: Test FASTA with various nucleotide characters
     */
    @Test
    void testImportVariousNucleotides() throws Exception {
        Path fastaFile = tempDir.resolve("nucleotides.fasta");
        String fastaContent = ">test_seq\n" +
            "ATGCNRYKMSWBDHVatgcnrykmswbdhv\n";
        Files.writeString(fastaFile, fastaContent);

        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile.toFile());

        assertEquals(1, docs.size(), "Should import sequence with various nucleotides");
    }

    /**
     * GROUP 4: Test FASTA with trailing whitespace
     */
    @Test
    void testImportTrailingWhitespace() throws Exception {
        Path fastaFile = tempDir.resolve("trailing.fasta");
        String fastaContent = ">seq1   \nATGC   \nGCTA   \n";
        Files.writeString(fastaFile, fastaContent);

        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile.toFile());

        assertEquals(1, docs.size(), "Should handle trailing whitespace");
        assertEquals("seq1", docs.get(0).getName());
    }

    /**
     * GROUP 4: Test large FASTA file performance
     */
    @Test
    void testImportLargeFastaFile() throws Exception {
        Path fastaFile = tempDir.resolve("large.fasta");
        StringBuilder content = new StringBuilder();

        // Create 100 sequences
        for (int i = 1; i <= 100; i++) {
            content.append(">sequence_").append(i).append("\n");
            content.append("ATGC".repeat(100)).append("\n");
        }

        Files.writeString(fastaFile, content.toString());

        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile.toFile());

        assertEquals(100, docs.size(), "Should import all 100 sequences");
        assertEquals("sequence_1", docs.get(0).getName());
        assertEquals("sequence_100", docs.get(99).getName());
    }

    /**
     * GROUP 4: Test consecutive headers (second header overwrites first)
     */
    @Test
    void testImportConsecutiveHeaders() throws Exception {
        Path fastaFile = tempDir.resolve("consecutive.fasta");
        String fastaContent = ">seq1\n" +
            ">seq2\n" +
            "ATGC\n";
        Files.writeString(fastaFile, fastaContent);

        List<AnnotatedPluginDocument> docs = SimpleFastaImporter.importFasta(fastaFile.toFile());

        // First header has no sequence, only second should be imported
        assertEquals(1, docs.size(), "Should import only sequence with data");
        assertEquals("seq2", docs.get(0).getName());
    }
}
