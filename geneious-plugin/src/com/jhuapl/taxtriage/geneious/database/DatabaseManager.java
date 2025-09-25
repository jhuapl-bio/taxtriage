package com.jhuapl.taxtriage.geneious.database;

import com.biomatters.geneious.publicapi.plugin.PluginUtilities;

import java.io.*;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.text.SimpleDateFormat;
import java.time.LocalDateTime;
import java.time.ZoneId;
import java.time.format.DateTimeFormatter;
import java.time.temporal.ChronoUnit;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Manages TaxTriage databases including bundled databases and cached downloads.
 * Provides intelligent database selection with fallback mechanisms and versioning.
 */
public class DatabaseManager {

    private static final Logger logger = Logger.getLogger(DatabaseManager.class.getName());

    // Cache configuration
    private static final String CACHE_BASE_DIR = System.getProperty("user.home") +
                                                  File.separator + ".geneious" +
                                                  File.separator + "plugin-data" +
                                                  File.separator + "taxtriage" +
                                                  File.separator + "databases";

    private static final String METADATA_FILE = "database-metadata.json";
    private static final long CACHE_EXPIRY_DAYS = 30; // Re-validate cached databases after 30 days
    private static final long MIN_DB_SIZE_BYTES = 1024 * 1024; // Minimum 1MB for valid database

    // Database types
    public enum DatabaseType {
        VIRAL("viral", "Viral reference database", "v1.0"),
        STANDARD("standard", "Standard reference database", "v1.0"),
        MINIKRAKEN("minikraken", "MiniKraken reference database", "v1.0");

        private final String id;
        private final String description;
        private final String bundledVersion;

        DatabaseType(String id, String description, String bundledVersion) {
            this.id = id;
            this.description = description;
            this.bundledVersion = bundledVersion;
        }

        public String getId() { return id; }
        public String getDescription() { return description; }
        public String getBundledVersion() { return bundledVersion; }
    }

    /**
     * Database metadata for tracking versions and validity
     */
    public static class DatabaseMetadata {
        public String databaseType;
        public String version;
        public long downloadTimestamp;
        public long sizeBytes;
        public String checksum;
        public String source;

        public DatabaseMetadata() {}

        public DatabaseMetadata(String type, String version, String source) {
            this.databaseType = type;
            this.version = version;
            this.source = source;
            this.downloadTimestamp = System.currentTimeMillis();
        }

        public boolean isExpired(long expiryDays) {
            long ageInDays = ChronoUnit.DAYS.between(
                LocalDateTime.ofInstant(new Date(downloadTimestamp).toInstant(), ZoneId.systemDefault()),
                LocalDateTime.now()
            );
            return ageInDays > expiryDays;
        }

        public String getFormattedDate() {
            SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
            return sdf.format(new Date(downloadTimestamp));
        }
    }

    private static DatabaseManager instance;
    private final Map<String, DatabaseMetadata> metadataCache = new HashMap<>();

    private DatabaseManager() {
        loadMetadata();
        ensureCacheDirectoryExists();
    }

    /**
     * Get singleton instance
     */
    public static synchronized DatabaseManager getInstance() {
        if (instance == null) {
            instance = new DatabaseManager();
        }
        return instance;
    }

    /**
     * Get the best available database path for a given type.
     * Priority: 1. Valid cached database, 2. Bundled database, 3. Download new
     *
     * @param type Database type
     * @param allowDownload Whether to download if not available locally
     * @return Path to database or null if not available
     */
    public String getDatabasePath(DatabaseType type, boolean allowDownload) {
        logger.info("==========================================");
        logger.info("DATABASE MANAGER: Requesting " + type.getDescription());
        logger.info("==========================================");

        // Check for cached database first
        DatabaseInfo cachedDb = getCachedDatabase(type);
        if (cachedDb != null && cachedDb.isValid()) {
            logger.info("✓ USING CACHED DATABASE:");
            logger.info("  Type: " + type.getId());
            logger.info("  Version: " + cachedDb.version);
            logger.info("  Location: " + cachedDb.path);
            logger.info("  Downloaded: " + cachedDb.metadata.getFormattedDate());
            logger.info("  Size: " + formatFileSize(cachedDb.metadata.sizeBytes));
            logger.info("  Source: " + cachedDb.metadata.source);
            logger.info("  Cache age: " + getCacheAge(cachedDb.metadata) + " days");
            logger.info("==========================================");
            return cachedDb.path;
        }

        // Check for bundled database
        DatabaseInfo bundledDb = getBundledDatabase(type);
        if (bundledDb != null && bundledDb.isValid()) {
            logger.info("✓ USING BUNDLED DATABASE:");
            logger.info("  Type: " + type.getId());
            logger.info("  Version: " + bundledDb.version + " (bundled)");
            logger.info("  Location: " + bundledDb.path);
            logger.info("  Source: Plugin bundle");
            logger.info("==========================================");

            // Copy bundled database to cache for future use
            copyBundledToCache(type, bundledDb);
            return bundledDb.path;
        }

        // Download if allowed
        if (allowDownload) {
            logger.info("⚠ No cached or bundled database found");
            logger.info("  Will attempt to download " + type.getId() + " database");
            logger.info("  Note: Download will be handled by Nextflow workflow");
            logger.info("==========================================");

            // Return special marker that tells workflow to download
            return "DOWNLOAD_REQUIRED";
        }

        logger.warning("✗ NO DATABASE AVAILABLE for " + type.getId());
        logger.warning("  No cached version found");
        logger.warning("  No bundled version found");
        logger.warning("  Download not permitted");
        logger.info("==========================================");
        return null;
    }

    /**
     * Store a downloaded database in the cache
     */
    public void cacheDownloadedDatabase(DatabaseType type, Path downloadedPath, String version) {
        logger.info("==========================================");
        logger.info("DATABASE MANAGER: Caching downloaded database");
        logger.info("  Type: " + type.getId());
        logger.info("  Version: " + version);
        logger.info("  Source path: " + downloadedPath);

        try {
            Path cacheDir = Paths.get(CACHE_BASE_DIR, type.getId());
            Files.createDirectories(cacheDir);

            // Copy database files to cache
            if (Files.isDirectory(downloadedPath)) {
                copyDirectory(downloadedPath, cacheDir);
            } else {
                Files.copy(downloadedPath, cacheDir.resolve(downloadedPath.getFileName()),
                          StandardCopyOption.REPLACE_EXISTING);
            }

            // Update metadata
            DatabaseMetadata metadata = new DatabaseMetadata(type.getId(), version, "Downloaded by workflow");
            metadata.sizeBytes = calculateDirectorySize(cacheDir);
            metadata.checksum = calculateChecksum(cacheDir);

            metadataCache.put(type.getId(), metadata);
            saveMetadata();

            logger.info("✓ Successfully cached database");
            logger.info("  Cache location: " + cacheDir);
            logger.info("  Size: " + formatFileSize(metadata.sizeBytes));
            logger.info("==========================================");

        } catch (IOException e) {
            logger.log(Level.WARNING, "Failed to cache database", e);
        }
    }

    /**
     * Clear cached databases
     */
    public void clearCache(DatabaseType type) {
        logger.info("Clearing cache for database: " + type.getId());

        Path cacheDir = Paths.get(CACHE_BASE_DIR, type.getId());
        if (Files.exists(cacheDir)) {
            try {
                deleteDirectory(cacheDir);
                metadataCache.remove(type.getId());
                saveMetadata();
                logger.info("Cache cleared for " + type.getId());
            } catch (IOException e) {
                logger.log(Level.WARNING, "Failed to clear cache", e);
            }
        }
    }

    /**
     * Get information about all available databases
     */
    public Map<DatabaseType, DatabaseInfo> getAllDatabaseInfo() {
        Map<DatabaseType, DatabaseInfo> info = new HashMap<>();

        for (DatabaseType type : DatabaseType.values()) {
            DatabaseInfo cached = getCachedDatabase(type);
            DatabaseInfo bundled = getBundledDatabase(type);

            if (cached != null && cached.isValid()) {
                info.put(type, cached);
            } else if (bundled != null && bundled.isValid()) {
                info.put(type, bundled);
            }
        }

        return info;
    }

    // Private helper methods

    private DatabaseInfo getCachedDatabase(DatabaseType type) {
        Path cacheDir = Paths.get(CACHE_BASE_DIR, type.getId());

        // Resolve symlinks and handle macOS /var -> /private/var mapping
        try {
            if (Files.exists(cacheDir)) {
                cacheDir = cacheDir.toRealPath();
            }
        } catch (IOException e) {
            // Ignore and use original path
        }

        if (!Files.exists(cacheDir)) {
            return null;
        }

        DatabaseMetadata metadata = metadataCache.get(type.getId());
        if (metadata == null) {
            return null;
        }

        // Check if cache is expired
        if (metadata.isExpired(CACHE_EXPIRY_DAYS)) {
            logger.info("Cached database " + type.getId() + " is expired (>" + CACHE_EXPIRY_DAYS + " days old)");
        }

        return new DatabaseInfo(cacheDir.toString(), metadata.version, metadata, true);
    }

    private DatabaseInfo getBundledDatabase(DatabaseType type) {
        // Check if bundled database exists in plugin resources
        try {
            String resourcePath = "/databases/" + type.getId();
            InputStream resourceStream = getClass().getResourceAsStream(resourcePath + "/hash.k2d");

            if (resourceStream != null) {
                resourceStream.close();

                // Extract bundled database to temp location if needed
                Path extractedPath = extractBundledDatabase(type);
                if (extractedPath != null) {
                    DatabaseMetadata metadata = new DatabaseMetadata(type.getId(),
                                                                    type.getBundledVersion(),
                                                                    "Bundled with plugin");
                    return new DatabaseInfo(extractedPath.toString(), type.getBundledVersion(), metadata, false);
                }
            }
        } catch (Exception e) {
            logger.log(Level.FINE, "No bundled database found for " + type.getId(), e);
        }

        return null;
    }

    private Path extractBundledDatabase(DatabaseType type) {
        try {
            Path tempDir = Files.createTempDirectory("taxtriage_bundled_" + type.getId());
            String resourceBase = "/databases/" + type.getId() + "/";

            // List of expected database files for Kraken2
            String[] dbFiles = {"hash.k2d", "opts.k2d", "taxo.k2d", "seqid2taxid.map"};

            for (String file : dbFiles) {
                InputStream is = getClass().getResourceAsStream(resourceBase + file);
                if (is != null) {
                    Path targetFile = tempDir.resolve(file);
                    Files.copy(is, targetFile, StandardCopyOption.REPLACE_EXISTING);
                    is.close();
                }
            }

            return tempDir;
        } catch (IOException e) {
            logger.log(Level.WARNING, "Failed to extract bundled database", e);
            return null;
        }
    }

    private void copyBundledToCache(DatabaseType type, DatabaseInfo bundledDb) {
        try {
            Path cacheDir = Paths.get(CACHE_BASE_DIR, type.getId());
            Files.createDirectories(cacheDir);

            Path sourcePath = Paths.get(bundledDb.path);
            if (Files.isDirectory(sourcePath)) {
                copyDirectory(sourcePath, cacheDir);
            }

            DatabaseMetadata metadata = new DatabaseMetadata(type.getId(),
                                                            bundledDb.version,
                                                            "Copied from bundle");
            metadata.sizeBytes = calculateDirectorySize(cacheDir);

            metadataCache.put(type.getId(), metadata);
            saveMetadata();

            logger.info("Copied bundled database to cache: " + type.getId());
        } catch (IOException e) {
            logger.log(Level.WARNING, "Failed to copy bundled database to cache", e);
        }
    }

    private void loadMetadata() {
        Path metadataFile = Paths.get(CACHE_BASE_DIR, METADATA_FILE);
        if (Files.exists(metadataFile)) {
            try {
                String json = new String(Files.readAllBytes(metadataFile));
                // Simple JSON parsing (could use a JSON library for production)
                parseMetadataJson(json);
            } catch (IOException e) {
                logger.log(Level.WARNING, "Failed to load metadata", e);
            }
        }
    }

    private void saveMetadata() {
        Path metadataFile = Paths.get(CACHE_BASE_DIR, METADATA_FILE);
        try {
            Files.createDirectories(metadataFile.getParent());
            String json = generateMetadataJson();
            Files.write(metadataFile, json.getBytes());
        } catch (IOException e) {
            logger.log(Level.WARNING, "Failed to save metadata", e);
        }
    }

    private void parseMetadataJson(String json) {
        // Simple JSON parsing without external dependencies
        // In production, use a proper JSON library
        metadataCache.clear();

        if (json == null || json.trim().isEmpty()) {
            return;
        }

        // Basic parsing logic - would be better with proper JSON library
        String[] entries = json.replaceAll("[{}\\[\\]]", "").split("},");
        for (String entry : entries) {
            DatabaseMetadata metadata = new DatabaseMetadata();
            String[] fields = entry.split(",");
            for (String field : fields) {
                String[] kv = field.split(":");
                if (kv.length == 2) {
                    String key = kv[0].trim().replaceAll("\"", "");
                    String value = kv[1].trim().replaceAll("\"", "");

                    switch (key) {
                        case "databaseType":
                            metadata.databaseType = value;
                            break;
                        case "version":
                            metadata.version = value;
                            break;
                        case "downloadTimestamp":
                            metadata.downloadTimestamp = Long.parseLong(value);
                            break;
                        case "sizeBytes":
                            metadata.sizeBytes = Long.parseLong(value);
                            break;
                        case "source":
                            metadata.source = value;
                            break;
                    }
                }
            }
            if (metadata.databaseType != null) {
                metadataCache.put(metadata.databaseType, metadata);
            }
        }
    }

    private String generateMetadataJson() {
        StringBuilder json = new StringBuilder("{\n");
        int count = 0;
        for (Map.Entry<String, DatabaseMetadata> entry : metadataCache.entrySet()) {
            if (count > 0) json.append(",\n");
            DatabaseMetadata m = entry.getValue();
            json.append("  \"").append(entry.getKey()).append("\": {\n");
            json.append("    \"databaseType\": \"").append(m.databaseType).append("\",\n");
            json.append("    \"version\": \"").append(m.version).append("\",\n");
            json.append("    \"downloadTimestamp\": ").append(m.downloadTimestamp).append(",\n");
            json.append("    \"sizeBytes\": ").append(m.sizeBytes).append(",\n");
            json.append("    \"source\": \"").append(m.source).append("\"\n");
            json.append("  }");
            count++;
        }
        json.append("\n}");
        return json.toString();
    }

    private void ensureCacheDirectoryExists() {
        try {
            Files.createDirectories(Paths.get(CACHE_BASE_DIR));
        } catch (IOException e) {
            logger.log(Level.WARNING, "Failed to create cache directory", e);
        }
    }

    private void copyDirectory(Path source, Path target) throws IOException {
        // Create target directory if it doesn't exist
        Files.createDirectories(target);

        Files.walk(source).forEach(sourcePath -> {
            try {
                Path targetPath = target.resolve(source.relativize(sourcePath));
                if (Files.isDirectory(sourcePath)) {
                    Files.createDirectories(targetPath);
                } else {
                    Files.copy(sourcePath, targetPath, StandardCopyOption.REPLACE_EXISTING);
                }
            } catch (IOException e) {
                if (!e.getMessage().contains("already exists")) {
                    logger.log(Level.WARNING, "Failed to copy file: " + sourcePath, e);
                }
            }
        });
    }

    private void deleteDirectory(Path dir) throws IOException {
        if (Files.exists(dir)) {
            Files.walk(dir)
                .sorted(Comparator.reverseOrder())
                .map(Path::toFile)
                .forEach(File::delete);
        }
    }

    private long calculateDirectorySize(Path dir) {
        try {
            return Files.walk(dir)
                .filter(Files::isRegularFile)
                .mapToLong(p -> {
                    try {
                        return Files.size(p);
                    } catch (IOException e) {
                        return 0;
                    }
                }).sum();
        } catch (IOException e) {
            return 0;
        }
    }

    private String calculateChecksum(Path path) {
        // Simple checksum based on file sizes and names
        // In production, use proper checksumming
        try {
            MessageDigest md = MessageDigest.getInstance("MD5");
            Files.walk(path)
                .filter(Files::isRegularFile)
                .sorted()
                .forEach(p -> {
                    md.update(p.getFileName().toString().getBytes());
                    try {
                        md.update(Long.toString(Files.size(p)).getBytes());
                    } catch (IOException e) {
                        // Ignore
                    }
                });

            byte[] digest = md.digest();
            StringBuilder sb = new StringBuilder();
            for (byte b : digest) {
                sb.append(String.format("%02x", b));
            }
            return sb.toString();
        } catch (Exception e) {
            return "unknown";
        }
    }

    private String formatFileSize(long bytes) {
        if (bytes < 1024) return bytes + " B";
        int exp = (int) (Math.log(bytes) / Math.log(1024));
        String pre = "KMGTPE".charAt(exp - 1) + "";
        return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
    }

    private long getCacheAge(DatabaseMetadata metadata) {
        return ChronoUnit.DAYS.between(
            LocalDateTime.ofInstant(new Date(metadata.downloadTimestamp).toInstant(), ZoneId.systemDefault()),
            LocalDateTime.now()
        );
    }

    /**
     * Database information holder
     */
    public static class DatabaseInfo {
        public final String path;
        public final String version;
        public final DatabaseMetadata metadata;
        public final boolean isCached;

        public DatabaseInfo(String path, String version, DatabaseMetadata metadata, boolean isCached) {
            this.path = path;
            this.version = version;
            this.metadata = metadata;
            this.isCached = isCached;
        }

        public boolean isValid() {
            if (path == null) return false;

            Path dbPath = Paths.get(path);
            if (!Files.exists(dbPath)) return false;

            // Check for minimum expected files for Kraken2 database
            Path hashFile = dbPath.resolve("hash.k2d");
            Path taxoFile = dbPath.resolve("taxo.k2d");

            return Files.exists(hashFile) || Files.exists(taxoFile);
        }
    }
}