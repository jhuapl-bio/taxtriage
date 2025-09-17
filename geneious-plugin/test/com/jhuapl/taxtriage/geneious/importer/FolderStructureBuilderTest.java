package com.jhuapl.taxtriage.geneious.importer;

import com.biomatters.geneious.publicapi.databaseservice.DatabaseServiceException;
import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
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
 * Test suite for FolderStructureBuilder functionality.
 *
 * Tests the creation of folder structures in Geneious that mirror
 * TaxTriage output directories.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class FolderStructureBuilderTest {

    @Mock
    private WritableDatabaseService mockDatabaseService;

    @TempDir
    Path tempDir;

    private FolderStructureBuilder folderBuilder;
    private File outputDirectory;

    @BeforeEach
    void setUp() throws IOException {
        MockitoAnnotations.openMocks(this);
        folderBuilder = new FolderStructureBuilder(mockDatabaseService);
        outputDirectory = tempDir.resolve("taxtriage_output").toFile();
        outputDirectory.mkdirs();
    }

    @Test
    void testCreateRootFolder() throws DatabaseServiceException {
        String workflowId = "test-workflow-123";

        AnnotatedPluginDocument rootFolder = folderBuilder.createRootFolder(workflowId);

        assertNotNull(rootFolder);
        assertEquals("TaxTriage_Results_" + workflowId, rootFolder.getName());
        assertEquals(1, folderBuilder.getFolderCount()); // Root folder should be counted
    }

    @Test
    void testCreateRootFolderTwice() throws DatabaseServiceException {
        String workflowId = "test-workflow-123";

        AnnotatedPluginDocument rootFolder1 = folderBuilder.createRootFolder(workflowId);
        AnnotatedPluginDocument rootFolder2 = folderBuilder.createRootFolder(workflowId);

        assertNotNull(rootFolder1);
        assertNotNull(rootFolder2);
        assertSame(rootFolder1, rootFolder2); // Should return the same instance
    }

    @Test
    void testGetOrCreateFolderForRootFile() throws IOException, DatabaseServiceException {
        folderBuilder.createRootFolder("test-workflow");

        File rootFile = new File(outputDirectory, "summary.txt");
        Files.writeString(rootFile.toPath(), "test content");

        AnnotatedPluginDocument folder = folderBuilder.getOrCreateFolder(outputDirectory, rootFile);

        assertNotNull(folder);
        assertEquals("TaxTriage_Results_test-workflow", folder.getName());
    }

    @Test
    void testGetOrCreateFolderForNestedFile() throws IOException, DatabaseServiceException {
        folderBuilder.createRootFolder("test-workflow");

        // Create nested directory structure
        File krakenDir = new File(outputDirectory, "kraken2");
        File sampleDir = new File(krakenDir, "sample1");
        sampleDir.mkdirs();

        File nestedFile = new File(sampleDir, "classification.txt");
        Files.writeString(nestedFile.toPath(), "test content");

        AnnotatedPluginDocument folder = folderBuilder.getOrCreateFolder(outputDirectory, nestedFile);

        assertNotNull(folder);
        assertTrue(folderBuilder.getFolderCount() > 1); // Should have created multiple folders
    }

    @Test
    void testCachingFolders() throws IOException, DatabaseServiceException {
        folderBuilder.createRootFolder("test-workflow");

        File krakenDir = new File(outputDirectory, "kraken2");
        krakenDir.mkdirs();

        File file1 = new File(krakenDir, "file1.txt");
        File file2 = new File(krakenDir, "file2.txt");
        Files.writeString(file1.toPath(), "content1");
        Files.writeString(file2.toPath(), "content2");

        AnnotatedPluginDocument folder1 = folderBuilder.getOrCreateFolder(outputDirectory, file1);
        AnnotatedPluginDocument folder2 = folderBuilder.getOrCreateFolder(outputDirectory, file2);

        assertNotNull(folder1);
        assertNotNull(folder2);
        assertSame(folder1, folder2); // Should return the same cached folder
    }

    @Test
    void testGetCreatedFolderPaths() throws IOException, DatabaseServiceException {
        folderBuilder.createRootFolder("test-workflow");

        // Create multi-level structure
        File level1 = new File(outputDirectory, "fastp");
        File level2 = new File(level1, "sample1");
        File level3 = new File(level2, "qc");
        level3.mkdirs();

        File deepFile = new File(level3, "report.html");
        Files.writeString(deepFile.toPath(), "test content");

        folderBuilder.getOrCreateFolder(outputDirectory, deepFile);

        String[] paths = folderBuilder.getCreatedFolderPaths();
        assertTrue(paths.length > 1);

        // Should contain the nested path
        boolean foundNestedPath = false;
        for (String path : paths) {
            if (path.contains("fastp") && path.contains("sample1") && path.contains("qc")) {
                foundNestedPath = true;
                break;
            }
        }
        assertTrue(foundNestedPath);
    }

    @Test
    void testHasCachedFolder() throws IOException, DatabaseServiceException {
        folderBuilder.createRootFolder("test-workflow");

        File krakenDir = new File(outputDirectory, "kraken2");
        krakenDir.mkdirs();
        File testFile = new File(krakenDir, "test.txt");
        Files.writeString(testFile.toPath(), "content");

        assertFalse(folderBuilder.hasCachedFolder("kraken2"));

        folderBuilder.getOrCreateFolder(outputDirectory, testFile);

        assertTrue(folderBuilder.hasCachedFolder("kraken2"));
    }

    @Test
    void testGetCachedFolder() throws IOException, DatabaseServiceException {
        folderBuilder.createRootFolder("test-workflow");

        File krakenDir = new File(outputDirectory, "kraken2");
        krakenDir.mkdirs();
        File testFile = new File(krakenDir, "test.txt");
        Files.writeString(testFile.toPath(), "content");

        assertNull(folderBuilder.getCachedFolder("kraken2"));

        AnnotatedPluginDocument folder = folderBuilder.getOrCreateFolder(outputDirectory, testFile);
        AnnotatedPluginDocument cachedFolder = folderBuilder.getCachedFolder("kraken2");

        assertNotNull(cachedFolder);
        assertSame(folder, cachedFolder);
    }

    @Test
    void testValidateStructure() throws DatabaseServiceException {
        // Before creating root folder
        assertFalse(folderBuilder.validateStructure());

        folderBuilder.createRootFolder("test-workflow");

        // After creating root folder
        assertTrue(folderBuilder.validateStructure());
    }

    @Test
    void testReset() throws IOException, DatabaseServiceException {
        folderBuilder.createRootFolder("test-workflow");

        File krakenDir = new File(outputDirectory, "kraken2");
        krakenDir.mkdirs();
        File testFile = new File(krakenDir, "test.txt");
        Files.writeString(testFile.toPath(), "content");

        folderBuilder.getOrCreateFolder(outputDirectory, testFile);

        assertTrue(folderBuilder.getFolderCount() > 0);
        assertTrue(folderBuilder.hasCachedFolder("kraken2"));

        folderBuilder.reset();

        assertEquals(0, folderBuilder.getFolderCount());
        assertFalse(folderBuilder.hasCachedFolder("kraken2"));
        assertFalse(folderBuilder.validateStructure()); // No root folder after reset
    }

    @Test
    void testGetFolderStructureSummary() throws IOException, DatabaseServiceException {
        folderBuilder.createRootFolder("test-workflow");

        // Create some folders
        File fastpDir = new File(outputDirectory, "fastp");
        File krakenDir = new File(outputDirectory, "kraken2");
        File sampleDir = new File(fastpDir, "sample1");
        sampleDir.mkdirs();
        krakenDir.mkdirs();

        File file1 = new File(sampleDir, "file1.txt");
        File file2 = new File(krakenDir, "file2.txt");
        Files.writeString(file1.toPath(), "content1");
        Files.writeString(file2.toPath(), "content2");

        folderBuilder.getOrCreateFolder(outputDirectory, file1);
        folderBuilder.getOrCreateFolder(outputDirectory, file2);

        String summary = folderBuilder.getFolderStructureSummary();

        assertNotNull(summary);
        assertTrue(summary.contains("Created Folder Structure"));
        assertTrue(summary.contains("TaxTriage_Results_test-workflow"));
    }

    @Test
    void testKnownTaxTriageFolders() throws IOException, DatabaseServiceException {
        folderBuilder.createRootFolder("test-workflow");

        // Test known TaxTriage folder names
        String[] knownFolders = {"fastp", "kraken2", "bracken", "count", "pipeline_info", "multiqc"};

        for (String folderName : knownFolders) {
            File folder = new File(outputDirectory, folderName);
            folder.mkdirs();
            File testFile = new File(folder, "test.txt");
            Files.writeString(testFile.toPath(), "content");

            AnnotatedPluginDocument folderDoc = folderBuilder.getOrCreateFolder(outputDirectory, testFile);
            assertNotNull(folderDoc);
        }

        assertTrue(folderBuilder.getFolderCount() > knownFolders.length); // Including root folder
    }

    @Test
    void testExceptionWhenNoRootFolder() throws IOException {
        // Don't create root folder

        File testFile = new File(outputDirectory, "test.txt");
        Files.writeString(testFile.toPath(), "content");

        assertThrows(IllegalStateException.class, () -> {
            folderBuilder.getOrCreateFolder(outputDirectory, testFile);
        });
    }

    @Test
    void testComplexFolderStructure() throws IOException, DatabaseServiceException {
        folderBuilder.createRootFolder("test-workflow");

        // Create a complex nested structure
        File[] nestedDirs = {
            new File(outputDirectory, "fastp/sample1/qc"),
            new File(outputDirectory, "fastp/sample2/qc"),
            new File(outputDirectory, "kraken2/sample1"),
            new File(outputDirectory, "kraken2/sample2"),
            new File(outputDirectory, "count/summary"),
            new File(outputDirectory, "multiqc/data")
        };

        for (File dir : nestedDirs) {
            dir.mkdirs();
            File testFile = new File(dir, "test.txt");
            Files.writeString(testFile.toPath(), "content");
            folderBuilder.getOrCreateFolder(outputDirectory, testFile);
        }

        // Should have created many folders
        assertTrue(folderBuilder.getFolderCount() > nestedDirs.length);

        // Validate structure
        assertTrue(folderBuilder.validateStructure());

        // Check that all paths are tracked
        String[] paths = folderBuilder.getCreatedFolderPaths();
        assertTrue(paths.length > 0);
    }

    @Test
    void testEdgeCaseFilePaths() throws IOException, DatabaseServiceException {
        folderBuilder.createRootFolder("test-workflow");

        // Test file with special characters in path (if supported by filesystem)
        File specialDir = new File(outputDirectory, "test_dir-123");
        specialDir.mkdirs();
        File specialFile = new File(specialDir, "test_file-456.txt");
        Files.writeString(specialFile.toPath(), "content");

        AnnotatedPluginDocument folder = folderBuilder.getOrCreateFolder(outputDirectory, specialFile);
        assertNotNull(folder);
    }
}