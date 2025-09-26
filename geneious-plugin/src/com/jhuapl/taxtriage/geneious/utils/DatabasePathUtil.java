package com.jhuapl.taxtriage.geneious.utils;

import com.jhuapl.taxtriage.geneious.database.DatabaseManager;
import com.jhuapl.taxtriage.geneious.database.DatabaseManager.DatabaseType;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.logging.Logger;

/**
 * Utility class for database path management and validation in TaxTriage workflows.
 *
 * <p>This utility centralizes database-related path operations, providing consistent
 * database location resolution, validation, and caching logic across the TaxTriage
 * plugin components.</p>
 *
 * <h3>Database Management Features:</h3>
 * <ul>
 *   <li><strong>Path Resolution:</strong> Consistent database path calculation</li>
 *   <li><strong>Validation:</strong> Database integrity and completeness checks</li>
 *   <li><strong>Type Mapping:</strong> Database name to type conversions</li>
 *   <li><strong>Cache Management:</strong> Database caching and retrieval</li>
 * </ul>
 *
 * <h3>Supported Database Types:</h3>
 * <ul>
 *   <li><strong>Viral:</strong> Kraken2 viral genome database</li>
 *   <li><strong>Standard:</strong> Kraken2 standard database</li>
 *   <li><strong>MiniKraken:</strong> Reduced-size Kraken2 database</li>
 * </ul>
 *
 * @author TaxTriage Development Team
 * @version 2.0
 * @since 2.0
 */
public final class DatabasePathUtil {

    private static final Logger logger = Logger.getLogger(DatabasePathUtil.class.getName());

    // Database cache constants
    private static final String CACHE_DIR_NAME = ".taxtriage-geneious";
    private static final String STANDARD_DB_NAME = "standard";

    // Required database files for validation
    private static final String HASH_FILE = "hash.k2d";
    private static final String TAXO_FILE = "taxo.k2d";
    private static final String OPTS_FILE = "opts.k2d";

    /**
     * Private constructor to prevent instantiation.
     */
    private DatabasePathUtil() {
        throw new UnsupportedOperationException("Utility class cannot be instantiated");
    }

    /**
     * Gets the path to the cached database directory for the specified database type.
     *
     * @param databaseName the database name/type
     * @return absolute path to the database directory
     */
    public static String getDatabasePath(String databaseName) {
        String userHome = System.getProperty("user.home");
        Path cacheDir = Paths.get(userHome, CACHE_DIR_NAME);

        String databaseType = databaseName;
        if (databaseType == null || databaseType.isEmpty()) {
            databaseType = STANDARD_DB_NAME;
        }

        Path databasePath = cacheDir.resolve(databaseType);
        return databasePath.toAbsolutePath().toString();
    }

    /**
     * Maps a database name string to a DatabaseType enum.
     *
     * @param databaseName the database name to map
     * @return corresponding DatabaseType or null if unknown
     */
    public static DatabaseType mapDatabaseNameToType(String databaseName) {
        if (databaseName == null) {
            return null;
        }

        String normalized = databaseName.toLowerCase().trim();

        if (normalized.equals("viral") || normalized.contains("viral")) {
            return DatabaseType.VIRAL;
        } else if (normalized.equals("standard") || normalized.contains("standard")) {
            return DatabaseType.STANDARD;
        } else if (normalized.equals("minikraken") || normalized.contains("mini")) {
            return DatabaseType.MINIKRAKEN;
        }

        return null; // Unknown database type
    }

    /**
     * Checks if a directory contains valid Kraken2 database files.
     *
     * @param dir the directory to check
     * @return true if the directory contains a valid Kraken2 database
     */
    public static boolean isValidDatabaseDirectory(Path dir) {
        if (dir == null || !Files.exists(dir) || !Files.isDirectory(dir)) {
            return false;
        }

        try {
            // Check for required Kraken2 database files
            boolean hasHashFile = Files.exists(dir.resolve(HASH_FILE));
            boolean hasTaxoFile = Files.exists(dir.resolve(TAXO_FILE));

            // opts.k2d is optional but common
            boolean hasOptsFile = Files.exists(dir.resolve(OPTS_FILE));

            // Valid if has at least hash and taxo files
            boolean isValid = hasHashFile && hasTaxoFile;

            if (isValid) {
                logger.fine("Valid database directory found: " + dir +
                           " (hash=" + hasHashFile + ", taxo=" + hasTaxoFile + ", opts=" + hasOptsFile + ")");
            } else {
                logger.fine("Invalid database directory: " + dir +
                           " (hash=" + hasHashFile + ", taxo=" + hasTaxoFile + ")");
            }

            return isValid;
        } catch (Exception e) {
            logger.warning("Error validating database directory " + dir + ": " + e.getMessage());
            return false;
        }
    }

    /**
     * Checks if a database exists and is valid at the specified path.
     *
     * @param databasePath the path to check
     * @return true if a valid database exists at the path
     */
    public static boolean databaseExists(String databasePath) {
        if (databasePath == null || databasePath.isEmpty()) {
            return false;
        }

        Path dbPath = Paths.get(databasePath);
        return isValidDatabaseDirectory(dbPath);
    }

    /**
     * Checks if a database has content by looking for database files.
     *
     * @param dbPath the database directory path
     * @return true if the database has content
     */
    public static boolean databaseHasContent(Path dbPath) {
        if (!Files.exists(dbPath) || !Files.isDirectory(dbPath)) {
            return false;
        }

        // Check if directory has database files (any .k2d files indicate a valid Kraken database)
        try {
            return Files.walk(dbPath)
                .anyMatch(path -> FileTypeUtil.isKrakenDatabaseFile(path.toFile()) ||
                                path.toString().contains("taxonomy"));
        } catch (Exception e) {
            logger.warning("Error checking database content: " + e.getMessage());
            return false;
        }
    }

    /**
     * Caches a downloaded database by copying it to the cache directory.
     *
     * @param dbType the database type
     * @param sourcePath the path to the downloaded database
     * @return true if caching was successful
     */
    public static boolean cacheDatabase(DatabaseType dbType, Path sourcePath) {
        if (dbType == null || sourcePath == null || !Files.exists(sourcePath)) {
            return false;
        }

        try {
            String cachePath = getDatabasePath(dbType.getId());
            Path targetPath = Paths.get(cachePath);

            // Create cache directory if needed
            Files.createDirectories(targetPath.getParent());

            // Copy database files
            if (Files.isDirectory(sourcePath)) {
                // Copy entire directory
                copyDatabaseDirectory(sourcePath, targetPath);
            } else {
                // Single file - copy to target directory
                Files.createDirectories(targetPath);
                Files.copy(sourcePath, targetPath.resolve(sourcePath.getFileName()));
            }

            // Use DatabaseManager to register the cached database
            DatabaseManager dbManager = DatabaseManager.getInstance();
            String version = "Downloaded " + LocalDateTime.now().format(DateTimeFormatter.ISO_LOCAL_DATE);
            dbManager.cacheDownloadedDatabase(dbType, targetPath, version);

            logger.info("Successfully cached database: " + dbType.getId() + " at " + targetPath);
            return true;

        } catch (Exception e) {
            logger.warning("Failed to cache database " + dbType.getId() + ": " + e.getMessage());
            return false;
        }
    }

    /**
     * Copies a database directory recursively.
     *
     * @param source the source directory
     * @param target the target directory
     * @throws IOException if copying fails
     */
    private static void copyDatabaseDirectory(Path source, Path target) throws IOException {
        Files.createDirectories(target);

        Files.walk(source)
            .forEach(sourcePath -> {
                try {
                    Path targetPath = target.resolve(source.relativize(sourcePath));
                    if (Files.isDirectory(sourcePath)) {
                        Files.createDirectories(targetPath);
                    } else {
                        Files.copy(sourcePath, targetPath);
                    }
                } catch (IOException e) {
                    logger.warning("Failed to copy database file: " + sourcePath + " -> " + e.getMessage());
                }
            });
    }

    /**
     * Gets the cache directory path.
     *
     * @return the cache directory path
     */
    public static Path getCacheDirectory() {
        String userHome = System.getProperty("user.home");
        return Paths.get(userHome, CACHE_DIR_NAME);
    }

    /**
     * Ensures the cache directory exists.
     *
     * @return true if the cache directory exists or was created successfully
     */
    public static boolean ensureCacheDirectoryExists() {
        try {
            Path cacheDir = getCacheDirectory();
            Files.createDirectories(cacheDir);
            return true;
        } catch (IOException e) {
            logger.warning("Failed to create cache directory: " + e.getMessage());
            return false;
        }
    }
}