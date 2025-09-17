package com.jhuapl.taxtriage.geneious.importer;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Test suite for FileFormatDetector functionality.
 *
 * Tests file format detection based on extensions and content analysis
 * for various TaxTriage output file types.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class FileFormatDetectorTest {

    @TempDir
    Path tempDir;

    private FileFormatDetector detector;

    @BeforeEach
    void setUp() {
        detector = new FileFormatDetector();
    }

    @Test
    void testDetectFastaByExtension() throws IOException {
        File fastaFile = tempDir.resolve("test.fasta").toFile();
        Files.writeString(fastaFile.toPath(), ">seq1\nATCGATCG\n>seq2\nGCTAGCTA\n");

        FileFormatDetector.FileFormat format = detector.detectFormat(fastaFile);

        assertEquals(FileFormatDetector.FileType.SEQUENCE, format.getType());
        assertEquals("fasta", format.getSubtype());
        assertFalse(format.isCompressed());
        assertTrue(format.isImportable());
    }

    @Test
    void testDetectFastqByExtension() throws IOException {
        File fastqFile = tempDir.resolve("test.fastq").toFile();
        Files.writeString(fastqFile.toPath(), "@read1\nATCGATCG\n+\nIIIIIIII\n");

        FileFormatDetector.FileFormat format = detector.detectFormat(fastqFile);

        assertEquals(FileFormatDetector.FileType.SEQUENCE, format.getType());
        assertEquals("fastq", format.getSubtype());
        assertFalse(format.isCompressed());
        assertTrue(format.isImportable());
    }

    @Test
    void testDetectCompressedFastq() throws IOException {
        File fastqGzFile = tempDir.resolve("test.fastq.gz").toFile();
        Files.writeString(fastqGzFile.toPath(), "@read1\nATCGATCG\n+\nIIIIIIII\n");

        FileFormatDetector.FileFormat format = detector.detectFormat(fastqGzFile);

        assertEquals(FileFormatDetector.FileType.SEQUENCE, format.getType());
        assertEquals("fastq", format.getSubtype());
        assertTrue(format.isCompressed());
        assertEquals("gz", format.getCompressionFormat());
        assertTrue(format.isImportable());
    }

    @Test
    void testDetectKrakenReport() throws IOException {
        File krakenFile = tempDir.resolve("sample.kreport").toFile();
        String krakenContent = "100.00\t1000\t500\tU\t0\tunclassified\n" +
                              "85.50\t855\t100\tR\t1\troot\n" +
                              "80.00\t800\t50\tD\t2\tBacteria\n";
        Files.writeString(krakenFile.toPath(), krakenContent);

        FileFormatDetector.FileFormat format = detector.detectFormat(krakenFile);

        assertEquals(FileFormatDetector.FileType.KRAKEN_REPORT, format.getType());
        assertEquals("kreport", format.getSubtype());
        assertTrue(format.isImportable());
    }

    @Test
    void testDetectHtmlReport() throws IOException {
        File htmlFile = tempDir.resolve("report.html").toFile();
        String htmlContent = "<!DOCTYPE html>\n<html><head><title>Report</title></head>" +
                           "<body><h1>MultiQC Report</h1></body></html>";
        Files.writeString(htmlFile.toPath(), htmlContent);

        FileFormatDetector.FileFormat format = detector.detectFormat(htmlFile);

        assertEquals(FileFormatDetector.FileType.REPORT, format.getType());
        assertEquals("html", format.getSubtype());
        assertTrue(format.isImportable());
    }

    @Test
    void testDetectCsvFile() throws IOException {
        File csvFile = tempDir.resolve("data.csv").toFile();
        Files.writeString(csvFile.toPath(), "Sample,Reads,Bases\nSample1,1000,50000\nSample2,2000,100000\n");

        FileFormatDetector.FileFormat format = detector.detectFormat(csvFile);

        assertEquals(FileFormatDetector.FileType.TEXT, format.getType());
        assertEquals("csv", format.getSubtype());
        assertTrue(format.isImportable());
    }

    @Test
    void testDetectTsvFile() throws IOException {
        File tsvFile = tempDir.resolve("data.tsv").toFile();
        Files.writeString(tsvFile.toPath(), "Sample\tReads\tBases\nSample1\t1000\t50000\nSample2\t2000\t100000\n");

        FileFormatDetector.FileFormat format = detector.detectFormat(tsvFile);

        assertEquals(FileFormatDetector.FileType.TEXT, format.getType());
        assertEquals("tsv", format.getSubtype());
        assertTrue(format.isImportable());
    }

    @Test
    void testDetectSamFile() throws IOException {
        File samFile = tempDir.resolve("alignment.sam").toFile();
        String samContent = "@HD\tVN:1.0\tSO:coordinate\n" +
                           "@SQ\tSN:chr1\tLN:249250621\n" +
                           "read1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tATCGATCG\tIIIIIIII\n";
        Files.writeString(samFile.toPath(), samContent);

        FileFormatDetector.FileFormat format = detector.detectFormat(samFile);

        assertEquals(FileFormatDetector.FileType.ALIGNMENT, format.getType());
        assertEquals("sam", format.getSubtype());
        assertTrue(format.isImportable());
    }

    @Test
    void testDetectVcfFile() throws IOException {
        File vcfFile = tempDir.resolve("variants.vcf").toFile();
        String vcfContent = "##fileformat=VCFv4.2\n" +
                           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                           "chr1\t100\t.\tA\tT\t60\tPASS\t.\n";
        Files.writeString(vcfFile.toPath(), vcfContent);

        FileFormatDetector.FileFormat format = detector.detectFormat(vcfFile);

        assertEquals(FileFormatDetector.FileType.VARIANT, format.getType());
        assertEquals("vcf", format.getSubtype());
        assertTrue(format.isImportable());
    }

    @Test
    void testDetectByContent() throws IOException {
        // Test FASTA detection by content (no extension)
        File unknownFile = tempDir.resolve("unknown").toFile();
        Files.writeString(unknownFile.toPath(), ">sequence1\nATCGATCGATCG\n>sequence2\nGCTAGCTAGCTA\n");

        FileFormatDetector.FileFormat format = detector.detectFormat(unknownFile);

        assertEquals(FileFormatDetector.FileType.SEQUENCE, format.getType());
        assertEquals("fasta", format.getSubtype());
        assertTrue(format.isImportable());
    }

    @Test
    void testDetectFastqByContent() throws IOException {
        File unknownFile = tempDir.resolve("unknown").toFile();
        Files.writeString(unknownFile.toPath(), "@read1\nATCGATCG\n+\nIIIIIIII\n@read2\nGCTAGCTA\n+\nJJJJJJJJ\n");

        FileFormatDetector.FileFormat format = detector.detectFormat(unknownFile);

        assertEquals(FileFormatDetector.FileType.SEQUENCE, format.getType());
        assertEquals("fastq", format.getSubtype());
        assertTrue(format.isImportable());
    }

    @Test
    void testDetectTaxTriagePattern() throws IOException {
        File krakenFile = tempDir.resolve("kraken_classification.txt").toFile();
        Files.writeString(krakenFile.toPath(), "Some kraken output content");

        FileFormatDetector.FileFormat format = detector.detectFormat(krakenFile);

        assertEquals(FileFormatDetector.FileType.TEXT, format.getType());
        assertEquals("taxtriage", format.getSubtype());
        assertTrue(format.isImportable());
    }

    @Test
    void testDetectUnknownFile() throws IOException {
        File unknownFile = tempDir.resolve("unknown.xyz").toFile();
        Files.writeString(unknownFile.toPath(), "\0\1\2\3binary content");

        FileFormatDetector.FileFormat format = detector.detectFormat(unknownFile);

        assertEquals(FileFormatDetector.FileType.UNKNOWN, format.getType());
        assertFalse(format.isImportable());
    }

    @Test
    void testShouldImportValidFiles() throws IOException {
        File fastqFile = tempDir.resolve("test.fastq").toFile();
        Files.writeString(fastqFile.toPath(), "@read1\nATCG\n+\nIIII\n");

        assertTrue(detector.shouldImport(fastqFile));
    }

    @Test
    void testShouldNotImportInvalidFiles() throws IOException {
        File unknownFile = tempDir.resolve("unknown.xyz").toFile();
        Files.writeString(unknownFile.toPath(), "\0\1\2\3binary content");

        assertFalse(detector.shouldImport(unknownFile));
    }

    @Test
    void testNullFile() {
        FileFormatDetector.FileFormat format = detector.detectFormat(null);

        assertEquals(FileFormatDetector.FileType.UNKNOWN, format.getType());
        assertFalse(format.isImportable());
    }

    @Test
    void testNonExistentFile() {
        File nonExistent = new File("/nonexistent/file.txt");
        FileFormatDetector.FileFormat format = detector.detectFormat(nonExistent);

        assertEquals(FileFormatDetector.FileType.UNKNOWN, format.getType());
        assertFalse(format.isImportable());
    }

    @Test
    void testEmptyFile() throws IOException {
        File emptyFile = tempDir.resolve("empty.txt").toFile();
        Files.writeString(emptyFile.toPath(), "");

        FileFormatDetector.FileFormat format = detector.detectFormat(emptyFile);

        assertEquals(FileFormatDetector.FileType.UNKNOWN, format.getType());
        assertFalse(format.isImportable());
    }

    @Test
    void testGetFormatDescription() throws IOException {
        File fastqFile = tempDir.resolve("test.fastq").toFile();
        Files.writeString(fastqFile.toPath(), "@read1\nATCG\n+\nIIII\n");

        String description = detector.getFormatDescription(fastqFile);

        assertNotNull(description);
        assertTrue(description.contains("Sequence"));
    }

    @Test
    void testMultipleExtensions() throws IOException {
        File compressedVcf = tempDir.resolve("variants.vcf.gz").toFile();
        Files.writeString(compressedVcf.toPath(), "##fileformat=VCFv4.2\n");

        FileFormatDetector.FileFormat format = detector.detectFormat(compressedVcf);

        assertEquals(FileFormatDetector.FileType.VARIANT, format.getType());
        assertEquals("vcf", format.getSubtype());
        assertTrue(format.isCompressed());
        assertEquals("gz", format.getCompressionFormat());
        assertTrue(format.isImportable());
    }

    @Test
    void testLogFile() throws IOException {
        File logFile = tempDir.resolve("pipeline.log").toFile();
        Files.writeString(logFile.toPath(), "INFO: Starting pipeline\nINFO: Processing sample\nINFO: Complete\n");

        FileFormatDetector.FileFormat format = detector.detectFormat(logFile);

        assertEquals(FileFormatDetector.FileType.TEXT, format.getType());
        assertEquals("log", format.getSubtype());
        assertTrue(format.isImportable());
    }

    @Test
    void testFormatToString() throws IOException {
        File fastqGzFile = tempDir.resolve("test.fastq.gz").toFile();
        Files.writeString(fastqGzFile.toPath(), "@read1\nATCG\n+\nIIII\n");

        FileFormatDetector.FileFormat format = detector.detectFormat(fastqGzFile);
        String toString = format.toString();

        assertNotNull(toString);
        assertTrue(toString.contains("Sequence"));
        assertTrue(toString.contains("fastq"));
        assertTrue(toString.contains("gz"));
    }
}