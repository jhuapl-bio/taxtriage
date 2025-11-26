package com.jhuapl.taxtriage.geneious.docker;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.attribute.PosixFilePermission;
import java.nio.file.attribute.PosixFilePermissions;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Handles path mapping and volume mounting between the host system and Docker containers.
 *
 * This class provides cross-platform support for converting local file paths to Docker
 * volume mounts, handling special characters, spaces in paths, and ensuring proper
 * permissions for container access.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class VolumeMapper {

    private static final Logger logger = Logger.getLogger(VolumeMapper.class.getName());

    /** Container mount points */
    public static final String CONTAINER_INPUT_PATH = "/input";
    public static final String CONTAINER_OUTPUT_PATH = "/output";
    public static final String CONTAINER_WORK_PATH = "/work";

    /** Windows drive letter pattern */
    private static final String WINDOWS_DRIVE_PATTERN = "^[A-Za-z]:.*";

    private final String platform;

    /**
     * Creates a new VolumeMapper for the current platform.
     */
    public VolumeMapper() {
        this.platform = System.getProperty("os.name").toLowerCase();
    }

    /**
     * Creates a new VolumeMapper for a specific platform (for testing).
     *
     * @param platform the target platform (windows, linux, darwin)
     */
    public VolumeMapper(String platform) {
        this.platform = platform.toLowerCase();
    }

    /**
     * Creates volume mount mappings for Docker containers.
     *
     * @param inputDirectory the local input directory
     * @param outputDirectory the local output directory
     * @param workDirectory the local work directory
     * @return map of host paths to container paths for volume mounting
     * @throws DockerException if path conversion fails
     */
    public Map<String, String> createVolumeMounts(Path inputDirectory, Path outputDirectory, Path workDirectory)
            throws DockerException {
        Map<String, String> mounts = new HashMap<>();

        try {
            // Convert and validate each path
            String inputMount = convertToDockerPath(inputDirectory.toAbsolutePath());
            String outputMount = convertToDockerPath(outputDirectory.toAbsolutePath());
            String workMount = convertToDockerPath(workDirectory.toAbsolutePath());

            // Ensure directories exist and have proper permissions
            ensureDirectoryAccess(inputDirectory);
            ensureDirectoryAccess(outputDirectory);
            ensureDirectoryAccess(workDirectory);

            mounts.put(inputMount, CONTAINER_INPUT_PATH);
            mounts.put(outputMount, CONTAINER_OUTPUT_PATH);
            mounts.put(workMount, CONTAINER_WORK_PATH);

            logger.info("Created volume mounts: " + mounts);
            return mounts;

        } catch (Exception e) {
            throw new DockerException("Failed to create volume mounts", e);
        }
    }

    /**
     * Converts a local file path to a Docker-compatible path format.
     *
     * @param localPath the local file system path
     * @return Docker-compatible path string
     * @throws DockerException if path conversion fails
     */
    public String convertToDockerPath(Path localPath) throws DockerException {
        String pathStr = localPath.toAbsolutePath().toString();

        if (isWindows()) {
            return convertWindowsPath(pathStr);
        } else {
            return convertUnixPath(pathStr);
        }
    }

    /**
     * Handles special characters and spaces in paths for Docker compatibility.
     *
     * @param path the path to escape
     * @return escaped path safe for Docker commands
     */
    public String escapePathForDocker(String path) {
        if (path == null || path.isEmpty()) {
            return path;
        }

        // For paths with spaces, wrap in quotes
        if (path.contains(" ")) {
            return "\"" + path.replace("\"", "\\\"") + "\"";
        }

        // Escape other special characters if needed
        return path.replace("$", "\\$")
                  .replace("`", "\\`")
                  .replace("!", "\\!");
    }

    /**
     * Validates that a path has proper permissions for Docker access.
     *
     * @param path the path to validate
     * @return true if permissions are valid, false otherwise
     */
    public boolean validatePermissions(Path path) {
        try {
            if (!Files.exists(path)) {
                logger.warning("Path does not exist: " + path);
                return false;
            }

            if (!Files.isReadable(path)) {
                logger.warning("Path is not readable: " + path);
                return false;
            }

            if (Files.isDirectory(path) && !Files.isWritable(path)) {
                logger.warning("Directory is not writable: " + path);
                return false;
            }

            // On Unix-like systems, check POSIX permissions
            if (!isWindows()) {
                return validatePosixPermissions(path);
            }

            return true;

        } catch (Exception e) {
            logger.log(Level.WARNING, "Error validating permissions for path: " + path, e);
            return false;
        }
    }

    /**
     * Gets the container work directory path.
     *
     * @return the container work directory path
     */
    public String getContainerWorkPath() {
        return CONTAINER_WORK_PATH;
    }

    /**
     * Gets the container input directory path.
     *
     * @return the container input directory path
     */
    public String getContainerInputPath() {
        return CONTAINER_INPUT_PATH;
    }

    /**
     * Gets the container output directory path.
     *
     * @return the container output directory path
     */
    public String getContainerOutputPath() {
        return CONTAINER_OUTPUT_PATH;
    }

    /**
     * Checks if the current platform is Windows.
     *
     * @return true if running on Windows, false otherwise
     */
    public boolean isWindows() {
        return platform.contains("windows");
    }

    /**
     * Checks if the current platform is macOS.
     *
     * @return true if running on macOS, false otherwise
     */
    public boolean isMacOS() {
        return platform.contains("mac") || platform.contains("darwin");
    }

    /**
     * Checks if the current platform is Linux.
     *
     * @return true if running on Linux, false otherwise
     */
    public boolean isLinux() {
        return platform.contains("linux");
    }

    /**
     * Converts a Windows path to Docker format.
     *
     * @param windowsPath the Windows file path
     * @return Docker-compatible path
     * @throws DockerException if conversion fails
     */
    private String convertWindowsPath(String windowsPath) throws DockerException {
        if (!windowsPath.matches(WINDOWS_DRIVE_PATTERN)) {
            throw new DockerException("Invalid Windows path format: " + windowsPath);
        }

        // Convert C:\path\to\file to /c/path/to/file
        String driveLetter = windowsPath.substring(0, 1).toLowerCase();
        String pathWithoutDrive = windowsPath.substring(2); // Remove "C:"

        String dockerPath = "/" + driveLetter + pathWithoutDrive.replace("\\", "/");

        // Handle UNC paths (\\server\share)
        if (windowsPath.startsWith("\\\\")) {
            // Convert \\server\share\path to //server/share/path
            dockerPath = windowsPath.replace("\\", "/");
        }

        logger.fine("Converted Windows path '" + windowsPath + "' to Docker path '" + dockerPath + "'");
        return dockerPath;
    }

    /**
     * Converts a Unix path to Docker format (typically no conversion needed).
     *
     * @param unixPath the Unix file path
     * @return Docker-compatible path
     */
    private String convertUnixPath(String unixPath) {
        // Unix paths are already compatible with Docker
        logger.fine("Unix path '" + unixPath + "' is Docker-compatible");
        return unixPath;
    }

    /**
     * Ensures a directory exists and has proper access permissions.
     *
     * @param directory the directory to check/create
     * @throws IOException if directory operations fail
     */
    private void ensureDirectoryAccess(Path directory) throws IOException {
        if (!Files.exists(directory)) {
            Files.createDirectories(directory);
            logger.info("Created directory: " + directory);
        }

        if (!Files.isDirectory(directory)) {
            throw new IOException("Path is not a directory: " + directory);
        }

        // Try to set appropriate permissions for Docker access
        if (!isWindows()) {
            try {
                Set<PosixFilePermission> perms = PosixFilePermissions.fromString("rwxr-xr-x");
                Files.setPosixFilePermissions(directory, perms);
                logger.fine("Set POSIX permissions for directory: " + directory);
            } catch (Exception e) {
                logger.log(Level.WARNING, "Could not set POSIX permissions for: " + directory, e);
            }
        }
    }

    /**
     * Validates POSIX permissions for Unix-like systems.
     *
     * @param path the path to validate
     * @return true if permissions are adequate, false otherwise
     */
    private boolean validatePosixPermissions(Path path) {
        try {
            Set<PosixFilePermission> perms = Files.getPosixFilePermissions(path);

            // Check for read permission
            if (!perms.contains(PosixFilePermission.OWNER_READ) &&
                !perms.contains(PosixFilePermission.GROUP_READ) &&
                !perms.contains(PosixFilePermission.OTHERS_READ)) {
                logger.warning("Path lacks read permissions: " + path);
                return false;
            }

            // For directories, check execute permission (needed to traverse)
            if (Files.isDirectory(path)) {
                if (!perms.contains(PosixFilePermission.OWNER_EXECUTE) &&
                    !perms.contains(PosixFilePermission.GROUP_EXECUTE) &&
                    !perms.contains(PosixFilePermission.OTHERS_EXECUTE)) {
                    logger.warning("Directory lacks execute permissions: " + path);
                    return false;
                }
            }

            return true;

        } catch (Exception e) {
            logger.log(Level.WARNING, "Error checking POSIX permissions for: " + path, e);
            return false;
        }
    }
}