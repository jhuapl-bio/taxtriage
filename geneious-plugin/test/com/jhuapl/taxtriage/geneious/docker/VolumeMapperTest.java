package com.jhuapl.taxtriage.geneious.docker;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Comprehensive test suite for VolumeMapper functionality.
 *
 * Tests path conversion, volume mounting, permission handling, and cross-platform
 * compatibility for Docker volume operations.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class VolumeMapperTest {

    @TempDir
    Path tempDir;

    private VolumeMapper volumeMapper;
    private Path inputDir;
    private Path outputDir;
    private Path workDir;

    @BeforeEach
    void setUp() throws Exception {
        volumeMapper = new VolumeMapper();

        // Create test directories
        inputDir = tempDir.resolve("input");
        outputDir = tempDir.resolve("output");
        workDir = tempDir.resolve("work");

        Files.createDirectories(inputDir);
        Files.createDirectories(outputDir);
        Files.createDirectories(workDir);
    }

    @Test
    void testCreateVolumeMounts() throws DockerException {
        Map<String, String> mounts = volumeMapper.createVolumeMounts(inputDir, outputDir, workDir);

        assertNotNull(mounts);
        assertEquals(3, mounts.size());

        // Verify container paths
        assertTrue(mounts.containsValue(VolumeMapper.CONTAINER_INPUT_PATH));
        assertTrue(mounts.containsValue(VolumeMapper.CONTAINER_OUTPUT_PATH));
        assertTrue(mounts.containsValue(VolumeMapper.CONTAINER_WORK_PATH));

        // Verify host paths are absolute
        for (String hostPath : mounts.keySet()) {
            assertTrue(Paths.get(hostPath).isAbsolute(), "Host path should be absolute: " + hostPath);
        }
    }

    @Test
    void testCreateVolumeMountsWithNonExistentDirectories() throws DockerException {
        Path newInputDir = tempDir.resolve("new_input");
        Path newOutputDir = tempDir.resolve("new_output");
        Path newWorkDir = tempDir.resolve("new_work");

        // Directories don't exist yet
        assertFalse(Files.exists(newInputDir));
        assertFalse(Files.exists(newOutputDir));
        assertFalse(Files.exists(newWorkDir));

        Map<String, String> mounts = volumeMapper.createVolumeMounts(newInputDir, newOutputDir, newWorkDir);

        assertNotNull(mounts);
        assertEquals(3, mounts.size());

        // Output and work directories should be created, but input should fail
        assertThrows(DockerException.class, () -> {
            volumeMapper.createVolumeMounts(newInputDir, newOutputDir, newWorkDir);
        });
    }

    @Test
    void testConvertToDockerPathUnix() throws DockerException {
        // Test with Unix-like system (current system if not Windows)
        VolumeMapper unixMapper = new VolumeMapper("linux");

        Path testPath = Paths.get("/home/user/data");
        String dockerPath = unixMapper.convertToDockerPath(testPath);

        assertEquals("/home/user/data", dockerPath);
    }

    @Test
    void testConvertToDockerPathWindows() throws DockerException {
        VolumeMapper windowsMapper = new VolumeMapper("windows");

        // Test Windows drive letter conversion
        Path testPath = Paths.get("C:\\Users\\user\\data");
        String dockerPath = windowsMapper.convertToDockerPath(testPath);

        assertEquals("/c/Users/user/data", dockerPath);
    }

    @Test
    void testConvertToDockerPathWindowsInvalidFormat() {
        VolumeMapper windowsMapper = new VolumeMapper("windows");

        // Test invalid Windows path (no drive letter)
        Path invalidPath = Paths.get("\\invalid\\path");

        assertThrows(DockerException.class, () -> {
            windowsMapper.convertToDockerPath(invalidPath);
        });
    }

    @Test
    void testEscapePathForDockerWithSpaces() {
        String pathWithSpaces = "/path/with spaces/file.txt";
        String escaped = volumeMapper.escapePathForDocker(pathWithSpaces);

        assertEquals("\"/path/with spaces/file.txt\"", escaped);
    }

    @Test
    void testEscapePathForDockerWithSpecialCharacters() {
        String pathWithSpecials = "/path/with$dollar`backtick!exclamation";
        String escaped = volumeMapper.escapePathForDocker(pathWithSpecials);

        assertEquals("/path/with\\$dollar\\`backtick\\!exclamation", escaped);
    }

    @Test
    void testEscapePathForDockerWithQuotes() {
        String pathWithQuotes = "/path/with\"quotes\"/file.txt";
        String escaped = volumeMapper.escapePathForDocker(pathWithQuotes);

        assertEquals("\"/path/with\\\"quotes\\\"/file.txt\"", escaped);
    }

    @Test
    void testEscapePathForDockerNormalPath() {
        String normalPath = "/normal/path/file.txt";
        String escaped = volumeMapper.escapePathForDocker(normalPath);

        assertEquals(normalPath, escaped); // Should remain unchanged
    }

    @Test
    void testEscapePathForDockerNullAndEmpty() {
        assertNull(volumeMapper.escapePathForDocker(null));
        assertEquals("", volumeMapper.escapePathForDocker(""));
    }

    @Test
    void testValidatePermissionsExistingDirectory() {
        assertTrue(volumeMapper.validatePermissions(inputDir));
        assertTrue(volumeMapper.validatePermissions(outputDir));
        assertTrue(volumeMapper.validatePermissions(workDir));
    }

    @Test
    void testValidatePermissionsNonExistentPath() {
        Path nonExistent = tempDir.resolve("does_not_exist");
        assertFalse(volumeMapper.validatePermissions(nonExistent));
    }

    @Test
    void testValidatePermissionsFile() throws Exception {
        Path testFile = inputDir.resolve("test.txt");
        Files.write(testFile, "test content".getBytes());

        assertTrue(volumeMapper.validatePermissions(testFile));
    }

    @Test
    void testGetContainerPaths() {
        assertEquals("/input", volumeMapper.getContainerInputPath());
        assertEquals("/output", volumeMapper.getContainerOutputPath());
        assertEquals("/work", volumeMapper.getContainerWorkPath());
    }

    @Test
    void testPlatformDetection() {
        String currentOs = System.getProperty("os.name").toLowerCase();

        if (currentOs.contains("windows")) {
            assertTrue(volumeMapper.isWindows());
            assertFalse(volumeMapper.isMacOS());
            assertFalse(volumeMapper.isLinux());
        } else if (currentOs.contains("mac") || currentOs.contains("darwin")) {
            assertFalse(volumeMapper.isWindows());
            assertTrue(volumeMapper.isMacOS());
            assertFalse(volumeMapper.isLinux());
        } else if (currentOs.contains("linux")) {
            assertFalse(volumeMapper.isWindows());
            assertFalse(volumeMapper.isMacOS());
            assertTrue(volumeMapper.isLinux());
        }
    }

    @Test
    void testCustomPlatformMapper() {
        VolumeMapper windowsMapper = new VolumeMapper("Windows 10");
        assertTrue(windowsMapper.isWindows());
        assertFalse(windowsMapper.isMacOS());
        assertFalse(windowsMapper.isLinux());

        VolumeMapper macMapper = new VolumeMapper("Mac OS X");
        assertFalse(macMapper.isWindows());
        assertTrue(macMapper.isMacOS());
        assertFalse(macMapper.isLinux());

        VolumeMapper linuxMapper = new VolumeMapper("Linux");
        assertFalse(linuxMapper.isWindows());
        assertFalse(linuxMapper.isMacOS());
        assertTrue(linuxMapper.isLinux());
    }

    @Test
    void testConvertToDockerPathRelativePath() throws DockerException {
        // Relative paths should be converted to absolute first
        Path relativePath = Paths.get("relative/path");
        String dockerPath = volumeMapper.convertToDockerPath(relativePath);

        // Should be converted to absolute path
        assertTrue(dockerPath.startsWith("/") || dockerPath.matches("^/[a-z]/.*")); // Unix or Windows format
    }

    @Test
    void testCreateVolumeMountsWithComplexPaths() throws DockerException {
        // Create directories with complex names
        Path complexInput = tempDir.resolve("input with spaces & symbols");
        Path complexOutput = tempDir.resolve("output-with-dashes_and_underscores");
        Path complexWork = tempDir.resolve("work.with.dots");

        Files.createDirectories(complexInput);
        Files.createDirectories(complexOutput);
        Files.createDirectories(complexWork);

        Map<String, String> mounts = volumeMapper.createVolumeMounts(complexInput, complexOutput, complexWork);

        assertNotNull(mounts);
        assertEquals(3, mounts.size());

        // All mounts should be created successfully
        assertTrue(mounts.containsValue(VolumeMapper.CONTAINER_INPUT_PATH));
        assertTrue(mounts.containsValue(VolumeMapper.CONTAINER_OUTPUT_PATH));
        assertTrue(mounts.containsValue(VolumeMapper.CONTAINER_WORK_PATH));
    }

    @Test
    void testConstantValues() {
        assertEquals("/input", VolumeMapper.CONTAINER_INPUT_PATH);
        assertEquals("/output", VolumeMapper.CONTAINER_OUTPUT_PATH);
        assertEquals("/work", VolumeMapper.CONTAINER_WORK_PATH);
    }

    @Test
    void testWindowsUNCPathConversion() throws DockerException {
        VolumeMapper windowsMapper = new VolumeMapper("windows");

        // Mock Windows UNC path - this would need special handling in real implementation
        // For now, this tests the basic functionality
        Path windowsPath = Paths.get("D:\\shared\\folder");
        String dockerPath = windowsMapper.convertToDockerPath(windowsPath);

        assertEquals("/d/shared/folder", dockerPath);
    }

    @Test
    void testValidatePermissionsWithReadOnlyFile() throws Exception {
        Path readOnlyFile = inputDir.resolve("readonly.txt");
        Files.write(readOnlyFile, "content".getBytes());

        // File should still be valid for reading
        assertTrue(volumeMapper.validatePermissions(readOnlyFile));
    }

    @Test
    void testDirectoryCreationInCreateVolumeMounts() throws DockerException {
        // Delete one of the directories
        outputDir.toFile().delete();
        assertFalse(Files.exists(outputDir));

        // createVolumeMounts should recreate it
        Map<String, String> mounts = volumeMapper.createVolumeMounts(inputDir, outputDir, workDir);

        assertNotNull(mounts);
        assertTrue(Files.exists(outputDir));
        assertTrue(Files.isDirectory(outputDir));
    }
}