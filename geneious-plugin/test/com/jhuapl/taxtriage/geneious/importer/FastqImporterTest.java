package com.jhuapl.taxtriage.geneious.importer;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Test suite for FastqImporter functionality.
 *
 * Tests the import and parsing of FASTQ sequence files from TaxTriage output.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class FastqImporterTest {

    @TempDir
    Path tempDir;

    private FastqImporter importer;

    @BeforeEach
    void setUp() {
        importer = new FastqImporter();
    }

    @Test
    void testImportValidFastqFile() throws IOException {
        File fastqFile = tempDir.resolve("test.fastq").toFile();
        String fastqContent = "@read1\n" +
                             "ATCGATCGATCGATCG\n" +
                             "+\n" +
                             "IIIIIIIIIIIIIIII\n" +
                             "@read2\n" +
                             "GCTAGCTAGCTAGCTA\n" +
                             "+\n" +
                             "JJJJJJJJJJJJJJJJ\n";
        Files.writeString(fastqFile.toPath(), fastqContent);

        List<AnnotatedPluginDocument> documents = importer.importFastqFile(fastqFile);

        assertNotNull(documents);
        assertEquals(2, documents.size()); // Two sequences should create two documents

        AnnotatedPluginDocument doc1 = documents.get(0);
        assertEquals("read1", doc1.getName());
        assertTrue(doc1.getDocument().getDescription().contains("Length: 16 bp"));

        AnnotatedPluginDocument doc2 = documents.get(1);
        assertEquals("read2", doc2.getName());
        assertTrue(doc2.getDocument().getDescription().contains("Length: 16 bp"));
    }

    @Test
    void testImportLargeFastqFile() throws IOException {
        File fastqFile = tempDir.resolve("large.fastq").toFile();
        StringBuilder fastqContent = new StringBuilder();

        // Create a file with 200 sequences (should create single multi-sequence document)
        for (int i = 1; i <= 200; i++) {
            fastqContent.append("@read").append(i).append("\n");
            fastqContent.append("ATCGATCGATCGATCG\n");
            fastqContent.append("+\n");
            fastqContent.append("IIIIIIIIIIIIIIII\n");
        }
        Files.writeString(fastqFile.toPath(), fastqContent.toString());

        List<AnnotatedPluginDocument> documents = importer.importFastqFile(fastqFile);

        assertNotNull(documents);
        assertEquals(1, documents.size()); // Large file should create single summary document

        AnnotatedPluginDocument doc = documents.get(0);
        assertEquals("large.fastq", doc.getName());
        assertTrue(doc.getDocument().getDescription().contains("200 sequences"));
    }

    @Test
    void testImportEmptyFastqFile() throws IOException {
        File emptyFile = tempDir.resolve("empty.fastq").toFile();
        Files.writeString(emptyFile.toPath(), "");

        List<AnnotatedPluginDocument> documents = importer.importFastqFile(emptyFile);

        assertNotNull(documents);
        assertTrue(documents.isEmpty()); // Empty file should return empty list
    }

    @Test
    void testImportMalformedFastqFile() throws IOException {
        File malformedFile = tempDir.resolve("malformed.fastq").toFile();
        String malformedContent = "@read1\n" +
                                 "ATCGATCGATCGATCG\n" +
                                 "+\n" +
                                 "IIIIIIIIIIIIIIII\n" +
                                 "@read2\n" +
                                 "GCTAGCTAGCTAGCTA\n" +
                                 "+\n"; // Missing quality line
        Files.writeString(malformedFile.toPath(), malformedContent);

        List<AnnotatedPluginDocument> documents = importer.importFastqFile(malformedFile);

        assertNotNull(documents);
        assertEquals(1, documents.size()); // Should import only valid records
    }

    @Test
    void testImportFastqWithDifferentQualityLengths() throws IOException {
        File fastqFile = tempDir.resolve("bad_quality.fastq").toFile();
        String fastqContent = "@read1\n" +
                             "ATCGATCGATCGATCG\n" +
                             "+\n" +
                             "III\n" + // Quality length doesn't match sequence length
                             "@read2\n" +
                             "GCTAGCTAGCTAGCTA\n" +
                             "+\n" +
                             "JJJJJJJJJJJJJJJJ\n"; // This one is correct
        Files.writeString(fastqFile.toPath(), fastqContent);

        List<AnnotatedPluginDocument> documents = importer.importFastqFile(fastqFile);

        assertNotNull(documents);
        assertEquals(1, documents.size()); // Should only import valid record (read2)
    }

    @Test
    void testFastqRecordMethods() {
        FastqImporter.FastqRecord record = new FastqImporter.FastqRecord(
                "@read1 description", "ATCGATCG", "+", "IIIIIIII"
        );

        assertEquals("@read1 description", record.getHeader());
        assertEquals("ATCGATCG", record.getSequence());
        assertEquals("+", record.getQualityHeader());
        assertEquals("IIIIIIII", record.getQuality());
        assertEquals("read1", record.getSequenceId()); // Should extract ID from header
    }

    @Test
    void testSequenceIdExtraction() {
        FastqImporter.FastqRecord record1 = new FastqImporter.FastqRecord(
                "@read1", "ATCG", "+", "IIII"
        );
        assertEquals("read1", record1.getSequenceId());

        FastqImporter.FastqRecord record2 = new FastqImporter.FastqRecord(
                "@read1 extra information", "ATCG", "+", "IIII"
        );
        assertEquals("read1", record2.getSequenceId());

        FastqImporter.FastqRecord record3 = new FastqImporter.FastqRecord(
                "no_at_symbol", "ATCG", "+", "IIII"
        );
        assertEquals("sequence", record3.getSequenceId()); // Default when no @ symbol
    }

    @Test
    void testIsFastqFileByExtension() throws IOException {
        File fastqFile1 = tempDir.resolve("test.fastq").toFile();
        Files.writeString(fastqFile1.toPath(), "@read1\nATCG\n+\nIIII\n");
        assertTrue(importer.isFastqFile(fastqFile1));

        File fastqFile2 = tempDir.resolve("test.fq").toFile();
        Files.writeString(fastqFile2.toPath(), "@read1\nATCG\n+\nIIII\n");
        assertTrue(importer.isFastqFile(fastqFile2));

        File fastqGzFile = tempDir.resolve("test.fastq.gz").toFile();
        Files.writeString(fastqGzFile.toPath(), "@read1\nATCG\n+\nIIII\n");
        assertTrue(importer.isFastqFile(fastqGzFile));

        File fqGzFile = tempDir.resolve("test.fq.gz").toFile();
        Files.writeString(fqGzFile.toPath(), "@read1\nATCG\n+\nIIII\n");
        assertTrue(importer.isFastqFile(fqGzFile));
    }

    @Test
    void testIsFastqFileByContent() throws IOException {
        File unknownFile = tempDir.resolve("unknown").toFile();
        String fastqContent = "@read1\n" +
                             "ATCGATCGATCGATCG\n" +
                             "+\n" +
                             "IIIIIIIIIIIIIIII\n";
        Files.writeString(unknownFile.toPath(), fastqContent);

        assertTrue(importer.isFastqFile(unknownFile));
    }

    @Test
    void testIsNotFastqFile() throws IOException {
        File textFile = tempDir.resolve("text.txt").toFile();
        Files.writeString(textFile.toPath(), "This is just a regular text file");
        assertFalse(importer.isFastqFile(textFile));

        File fastaFile = tempDir.resolve("sequences.fasta").toFile();
        Files.writeString(fastaFile.toPath(), ">seq1\nATCGATCG\n>seq2\nGCTAGCTA\n");
        assertFalse(importer.isFastqFile(fastaFile));
    }

    @Test
    void testImportNullFile() {
        assertThrows(IOException.class, () -> {
            importer.importFastqFile(null);
        });
    }

    @Test
    void testImportNonExistentFile() {
        File nonExistent = new File("/nonexistent/file.fastq");
        assertThrows(IOException.class, () -> {
            importer.importFastqFile(nonExistent);
        });
    }

    @Test
    void testFastqWithComments() throws IOException {
        File fastqFile = tempDir.resolve("with_comments.fastq").toFile();
        String fastqContent = "@read1 some extra information\n" +
                             "ATCGATCGATCGATCG\n" +
                             "+read1 quality header\n" +
                             "IIIIIIIIIIIIIIII\n";
        Files.writeString(fastqFile.toPath(), fastqContent);

        List<AnnotatedPluginDocument> documents = importer.importFastqFile(fastqFile);

        assertNotNull(documents);
        assertEquals(1, documents.size());

        AnnotatedPluginDocument doc = documents.get(0);
        assertEquals("read1", doc.getName()); // Should extract just the ID
        assertTrue(doc.getDocument().getDescription().contains("@read1 some extra information"));
        assertTrue(doc.getDocument().getDescription().contains("+read1 quality header"));
    }

    @Test
    void testFastqWithEmptyLines() throws IOException {
        File fastqFile = tempDir.resolve("with_empty_lines.fastq").toFile();
        String fastqContent = "\n" +
                             "@read1\n" +
                             "ATCGATCGATCGATCG\n" +
                             "+\n" +
                             "IIIIIIIIIIIIIIII\n" +
                             "\n" +
                             "@read2\n" +
                             "GCTAGCTAGCTAGCTA\n" +
                             "+\n" +
                             "JJJJJJJJJJJJJJJJ\n" +
                             "\n";
        Files.writeString(fastqFile.toPath(), fastqContent);

        List<AnnotatedPluginDocument> documents = importer.importFastqFile(fastqFile);

        assertNotNull(documents);
        assertEquals(2, documents.size()); // Should handle empty lines gracefully
    }

    @Test
    void testMultiSequenceDocumentGeneration() throws IOException {
        File fastqFile = tempDir.resolve("multi.fastq").toFile();
        StringBuilder fastqContent = new StringBuilder();

        // Create 150 sequences to trigger multi-sequence document
        for (int i = 1; i <= 150; i++) {
            fastqContent.append("@read").append(i).append("\n");
            fastqContent.append("ATCGATCGATCGATCGATCGATCGATCGATCG\n"); // 32 bp
            fastqContent.append("+\n");
            fastqContent.append("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");
        }
        Files.writeString(fastqFile.toPath(), fastqContent.toString());

        List<AnnotatedPluginDocument> documents = importer.importFastqFile(fastqFile);

        assertNotNull(documents);
        assertEquals(1, documents.size()); // Should create single multi-sequence document

        AnnotatedPluginDocument doc = documents.get(0);
        String description = doc.getDocument().getDescription();
        assertTrue(description.contains("150 sequences"));
        assertTrue(description.contains("Average Length: 32.0"));
        assertTrue(description.contains("Total Bases: 4800"));

        String content = doc.getDocument().toString();
        assertTrue(content.contains("=== FASTQ File Summary ==="));
        assertTrue(content.contains("=== Sequence Statistics ==="));
        assertTrue(content.contains("=== Sample Sequences (first 5) ==="));
    }

    @Test
    void testIsFastqFileEdgeCases() {
        // Test null file
        assertFalse(importer.isFastqFile(null));

        // Test non-existent file
        File nonExistent = new File("/nonexistent/file.fastq");
        assertFalse(importer.isFastqFile(nonExistent));

        // Test directory
        assertFalse(importer.isFastqFile(tempDir.toFile()));
    }

    @Test
    void testFastqWithVaryingSequenceLengths() throws IOException {
        File fastqFile = tempDir.resolve("varying_lengths.fastq").toFile();
        String fastqContent = "@read1\n" +
                             "ATCG\n" +
                             "+\n" +
                             "IIII\n" +
                             "@read2\n" +
                             "ATCGATCGATCGATCGATCGATCGATCGATCG\n" +
                             "+\n" +
                             "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n" +
                             "@read3\n" +
                             "ATCGATCGATCG\n" +
                             "+\n" +
                             "IIIIIIIIIIII\n";
        Files.writeString(fastqFile.toPath(), fastqContent.toString());

        List<AnnotatedPluginDocument> documents = importer.importFastqFile(fastqFile);

        assertNotNull(documents);
        assertEquals(3, documents.size());

        // Check that all different lengths are handled correctly
        assertEquals("read1", documents.get(0).getName());
        assertEquals("read2", documents.get(1).getName());
        assertEquals("read3", documents.get(2).getName());
    }
}