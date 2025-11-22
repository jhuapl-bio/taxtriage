package com.jhuapl.taxtriage.geneious.utils;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for FileTypeUtil (Group 3).
 *
 * Tests file type detection and utility functionality including:
 * - File extension detection
 * - File validation
 * - Format checking
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class FileTypeUtilTest {

    @TempDir
    Path tempDir;

    /**
     * GROUP 3: Test FASTQ file detection
     */
    @Test
    void testIsFastqFile() throws Exception {
        File fastqFile = createTestFile("test.fastq", "@read1\nACGT\n+\nIIII\n");
        assertTrue(FileTypeUtil.isFastqFile(fastqFile),
            "Should recognize .fastq files");

        File fqFile = createTestFile("test.fq", "@read1\nACGT\n+\nIIII\n");
        assertTrue(FileTypeUtil.isFastqFile(fqFile),
            "Should recognize .fq files");

        File gzFile = createTestFile("test.fastq.gz", "dummy");
        assertTrue(FileTypeUtil.isFastqFile(gzFile),
            "Should recognize .fastq.gz files");
    }

    /**
     * GROUP 3: Test FASTA file detection
     */
    @Test
    void testIsFastaFile() throws Exception {
        File fastaFile = createTestFile("test.fasta", ">seq1\nACGT\n");
        assertTrue(FileTypeUtil.isFastaFile(fastaFile),
            "Should recognize .fasta files");

        File faFile = createTestFile("test.fa", ">seq1\nACGT\n");
        assertTrue(FileTypeUtil.isFastaFile(faFile),
            "Should recognize .fa files");

        File fnaFile = createTestFile("test.fna", ">seq1\nACGT\n");
        assertTrue(FileTypeUtil.isFastaFile(fnaFile),
            "Should recognize .fna files");
    }

    /**
     * GROUP 3: Test BAM file detection
     */
    @Test
    void testIsBamFile() throws Exception {
        File bamFile = createTestFile("test.bam", "dummy");
        assertTrue(FileTypeUtil.isBamFile(bamFile),
            "Should recognize .bam files");

        File txtFile = createTestFile("test.txt", "text");
        assertFalse(FileTypeUtil.isBamFile(txtFile),
            "Should not recognize .txt as BAM");
    }

    /**
     * GROUP 3: Test GenBank file detection
     */
    @Test
    void testIsGenBankFile() throws Exception {
        File gbFile = createTestFile("test.gb", "LOCUS");
        assertTrue(FileTypeUtil.isGenBankFile(gbFile),
            "Should recognize .gb files");

        File gbkFile = createTestFile("test.gbk", "LOCUS");
        assertTrue(FileTypeUtil.isGenBankFile(gbkFile),
            "Should recognize .gbk files");

        File genbankFile = createTestFile("test.genbank", "LOCUS");
        assertTrue(FileTypeUtil.isGenBankFile(genbankFile),
            "Should recognize .genbank files");
    }

    /**
     * GROUP 3: Test alignment file detection
     */
    @Test
    void testIsAlignmentFile() throws Exception {
        File bamFile = createTestFile("test.bam", "dummy");
        assertTrue(FileTypeUtil.isAlignmentFile(bamFile),
            "BAM should be alignment file");

        File samFile = createTestFile("test.sam", "dummy");
        assertTrue(FileTypeUtil.isAlignmentFile(samFile),
            "SAM should be alignment file");

        File fastaFile = createTestFile("test.fasta", ">seq");
        assertFalse(FileTypeUtil.isAlignmentFile(fastaFile),
            "FASTA should not be alignment file");
    }

    /**
     * GROUP 3: Test valid sequence file detection
     */
    @Test
    void testIsValidSequenceFile() throws Exception {
        File fastqFile = createTestFile("test.fastq", "@read\nACGT\n+\nIIII");
        assertTrue(FileTypeUtil.isValidSequenceFile(fastqFile),
            "FASTQ should be valid sequence file");

        File excludedFile = createTestFile("test.fastp.fastq", "@read\n");
        assertFalse(FileTypeUtil.isValidSequenceFile(excludedFile),
            "Files with .fastp. should be excluded");

        File dwnldFile = createTestFile("test.dwnld.fastq", "@read\n");
        assertFalse(FileTypeUtil.isValidSequenceFile(dwnldFile),
            "Files with .dwnld. should be excluded");
    }

    /**
     * GROUP 3: Test null file handling
     */
    @Test
    void testNullFileHandling() {
        assertFalse(FileTypeUtil.isFastqFile((File) null),
            "Should handle null file gracefully");
        assertFalse(FileTypeUtil.isFastaFile(null),
            "Should handle null file gracefully");
        assertFalse(FileTypeUtil.isBamFile(null),
            "Should handle null file gracefully");
        assertFalse(FileTypeUtil.isValidSequenceFile(null),
            "Should handle null file gracefully");
    }

    /**
     * GROUP 3: Test non-existent file handling
     */
    @Test
    void testNonExistentFile() {
        File nonExistent = new File(tempDir.toFile(), "does_not_exist.fastq");
        assertFalse(FileTypeUtil.isFastqFile(nonExistent),
            "Non-existent files should return false");
    }

    /**
     * GROUP 3: Test case insensitive detection
     */
    @Test
    void testCaseInsensitiveDetection() throws Exception {
        assertTrue(FileTypeUtil.isFastqFile(createTestFile("test.FASTQ", "data")),
            "Should be case insensitive for FASTQ");
        assertTrue(FileTypeUtil.isFastaFile(createTestFile("test.FASTA", "data")),
            "Should be case insensitive for FASTA");
        assertTrue(FileTypeUtil.isBamFile(createTestFile("test.BAM", "data")),
            "Should be case insensitive for BAM");
    }

    /**
     * GROUP 3: Test multiple extension handling
     */
    @Test
    void testMultipleExtensions() throws Exception {
        File fastqGz = createTestFile("sample.fastq.gz", "data");
        assertTrue(FileTypeUtil.isFastqFile(fastqGz),
            "Should handle .fastq.gz");

        File fastaGz = createTestFile("sample.fasta.gz", "data");
        assertTrue(FileTypeUtil.isFastaFile(fastaGz),
            "Should handle .fasta.gz");
    }

    /**
     * GROUP 3: Test is text file
     */
    @Test
    void testIsTextFile() throws Exception {
        assertTrue(FileTypeUtil.isTextFile(createTestFile("test.txt", "text")));
        assertTrue(FileTypeUtil.isTextFile(createTestFile("test.tsv", "tab\tsep")));
        assertTrue(FileTypeUtil.isTextFile(createTestFile("test.csv", "comma,sep")));
        assertFalse(FileTypeUtil.isTextFile(createTestFile("test.bam", "data")));
    }

    /**
     * GROUP 3: Test Kraken database file detection
     */
    @Test
    void testIsKrakenDatabaseFile() throws Exception {
        File k2dFile = createTestFile("test.k2d", "data");
        assertTrue(FileTypeUtil.isKrakenDatabaseFile(k2dFile),
            "Should recognize .k2d files");

        File kmerFile = createTestFile("test.kmer_distrib", "data");
        assertTrue(FileTypeUtil.isKrakenDatabaseFile(kmerFile),
            "Should recognize .kmer_distrib files");
    }

    /**
     * GROUP 3: Test deduplicated BAM detection
     */
    @Test
    void testIsDeduplicatedBamFile() throws Exception {
        File dedupBam = createTestFile("sample.dedup.bam", "data");
        assertTrue(FileTypeUtil.isDeduplicatedBamFile(dedupBam),
            "Should recognize deduplicated BAM files");

        File normalBam = createTestFile("sample.bam", "data");
        assertFalse(FileTypeUtil.isDeduplicatedBamFile(normalBam),
            "Should not recognize normal BAM as deduplicated");
    }

    /**
     * GROUP 3: Test original BAM detection
     */
    @Test
    void testIsOriginalBamFile() throws Exception {
        File normalBam = createTestFile("sample.bam", "data");
        assertTrue(FileTypeUtil.isOriginalBamFile(normalBam),
            "Should recognize original BAM files");

        File dedupBam = createTestFile("sample.dedup.bam", "data");
        assertFalse(FileTypeUtil.isOriginalBamFile(dedupBam),
            "Should not recognize deduplicated BAM as original");
    }

    /**
     * GROUP 3: Test get base name
     */
    @Test
    void testGetBaseName() {
        assertEquals("test", FileTypeUtil.getBaseName("test.txt"));
        assertEquals("sample", FileTypeUtil.getBaseName("sample.fastq"));
        assertEquals("file", FileTypeUtil.getBaseName("file.fastq.gz"));
        assertEquals("", FileTypeUtil.getBaseName(""));
        assertEquals("", FileTypeUtil.getBaseName(null));
        assertEquals("noextension", FileTypeUtil.getBaseName("noextension"));
    }

    /**
     * GROUP 3: Test path variants
     */
    @Test
    void testPathVariants() throws Exception {
        Path fastqPath = tempDir.resolve("test.fastq");
        Files.writeString(fastqPath, "@read\nACGT\n+\nIIII");

        assertTrue(FileTypeUtil.isFastqFile(fastqPath),
            "Should handle Path objects");
    }

    // ===== Helper Methods =====

    private File createTestFile(String name, String content) throws Exception {
        File file = tempDir.resolve(name).toFile();
        Files.writeString(file.toPath(), content);
        return file;
    }
}
