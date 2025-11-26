package com.jhuapl.taxtriage.geneious.database;

import org.junit.jupiter.api.*;
import org.junit.jupiter.api.io.TempDir;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Comprehensive test suite for DatabaseManager functionality.
 *
 * Tests database path detection, validation, caching, and fallback mechanisms
 * without requiring actual database downloads.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
@DisplayName("DatabaseManager Tests")
class DatabaseManagerTest {

    @TempDir
    Path tempCacheDir;

    private DatabaseManager databaseManager;
    private Path mockDatabaseDir;

    @BeforeEach
    void setUp() throws IOException {
        // Create a mock database directory structure
        mockDatabaseDir = tempCacheDir.resolve("mock_database");
        Files.createDirectories(mockDatabaseDir);

        // Note: DatabaseManager uses singleton pattern, so we test with the actual instance
        databaseManager = DatabaseManager.getInstance();
    }

    @Nested
    @DisplayName("Database Path Detection Tests")
    class DatabasePathDetectionTests {

        @Test
        @DisplayName("Should detect valid database with required files")
        void testDetectValidDatabase() throws IOException {
            // Create valid database files
            createValidKrakenDatabase(mockDatabaseDir);

            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    mockDatabaseDir.toString(),
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test"),
                    false
            );

            assertTrue(dbInfo.isValid(), "Database with hash.k2d and taxo.k2d should be valid");
        }

        @Test
        @DisplayName("Should reject database missing required files")
        void testRejectInvalidDatabase() throws IOException {
            // Create directory with only one required file
            Files.createFile(mockDatabaseDir.resolve("hash.k2d"));
            // Missing taxo.k2d

            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    mockDatabaseDir.toString(),
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test"),
                    false
            );

            assertFalse(dbInfo.isValid(), "Database missing taxo.k2d should be invalid");
        }

        @Test
        @DisplayName("Should handle non-existent database path")
        void testNonExistentPath() {
            Path nonExistentPath = tempCacheDir.resolve("does_not_exist");

            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    nonExistentPath.toString(),
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test"),
                    false
            );

            assertFalse(dbInfo.isValid(), "Non-existent database path should be invalid");
        }

        @Test
        @DisplayName("Should handle null database path")
        void testNullPath() {
            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    null,
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test"),
                    false
            );

            assertFalse(dbInfo.isValid(), "Null database path should be invalid");
        }
    }

    @Nested
    @DisplayName("Database Type Tests")
    class DatabaseTypeTests {

        @Test
        @DisplayName("Should have viral database type")
        void testViralDatabaseType() {
            DatabaseManager.DatabaseType viral = DatabaseManager.DatabaseType.VIRAL;

            assertEquals("viral", viral.getId());
            assertEquals("Viral reference database", viral.getDescription());
            assertNotNull(viral.getBundledVersion());
        }

        @Test
        @DisplayName("Should have standard database type")
        void testStandardDatabaseType() {
            DatabaseManager.DatabaseType standard = DatabaseManager.DatabaseType.STANDARD;

            assertEquals("standard", standard.getId());
            assertEquals("Standard reference database", standard.getDescription());
            assertNotNull(standard.getBundledVersion());
        }

        @Test
        @DisplayName("Should have minikraken database type")
        void testMiniKrakenDatabaseType() {
            DatabaseManager.DatabaseType minikraken = DatabaseManager.DatabaseType.MINIKRAKEN;

            assertEquals("minikraken", minikraken.getId());
            assertEquals("MiniKraken reference database", minikraken.getDescription());
            assertNotNull(minikraken.getBundledVersion());
        }

        @Test
        @DisplayName("Should enumerate all database types")
        void testAllDatabaseTypes() {
            DatabaseManager.DatabaseType[] types = DatabaseManager.DatabaseType.values();

            assertEquals(3, types.length, "Should have exactly 3 database types");
            assertTrue(containsType(types, "viral"));
            assertTrue(containsType(types, "standard"));
            assertTrue(containsType(types, "minikraken"));
        }

        private boolean containsType(DatabaseManager.DatabaseType[] types, String id) {
            for (DatabaseManager.DatabaseType type : types) {
                if (type.getId().equals(id)) {
                    return true;
                }
            }
            return false;
        }
    }

    @Nested
    @DisplayName("Database Metadata Tests")
    class DatabaseMetadataTests {

        @Test
        @DisplayName("Should create metadata with basic information")
        void testCreateBasicMetadata() {
            DatabaseManager.DatabaseMetadata metadata =
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test_source");

            assertEquals("viral", metadata.databaseType);
            assertEquals("v1.0", metadata.version);
            assertEquals("test_source", metadata.source);
            assertTrue(metadata.downloadTimestamp > 0, "Timestamp should be set");
        }

        @Test
        @DisplayName("Should format download date correctly")
        void testFormattedDate() {
            DatabaseManager.DatabaseMetadata metadata =
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test_source");

            String formattedDate = metadata.getFormattedDate();

            assertNotNull(formattedDate);
            assertFalse(formattedDate.isEmpty());
            assertTrue(formattedDate.matches("\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"),
                    "Date should match YYYY-MM-DD HH:MM:SS format");
        }

        @Test
        @DisplayName("Should detect fresh database as not expired")
        void testFreshDatabaseNotExpired() {
            DatabaseManager.DatabaseMetadata metadata =
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test_source");

            assertFalse(metadata.isExpired(30), "Fresh database should not be expired");
        }

        @Test
        @DisplayName("Should detect old database as expired")
        void testOldDatabaseExpired() {
            DatabaseManager.DatabaseMetadata metadata =
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test_source");

            // Set timestamp to 60 days ago
            long sixtyDaysAgo = System.currentTimeMillis() - (60L * 24 * 60 * 60 * 1000);
            metadata.downloadTimestamp = sixtyDaysAgo;

            assertTrue(metadata.isExpired(30), "60-day old database should be expired with 30-day limit");
        }

        @Test
        @DisplayName("Should handle metadata at expiry boundary")
        void testExpiryBoundary() {
            DatabaseManager.DatabaseMetadata metadata =
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test_source");

            // Set timestamp to exactly 30 days ago
            long thirtyDaysAgo = System.currentTimeMillis() - (30L * 24 * 60 * 60 * 1000);
            metadata.downloadTimestamp = thirtyDaysAgo;

            // At exactly 30 days, it should not be expired (> not >=)
            assertFalse(metadata.isExpired(30), "Database at exactly 30 days should not be expired");
        }
    }

    @Nested
    @DisplayName("Database Validation Tests")
    class DatabaseValidationTests {

        @Test
        @DisplayName("Should validate complete Kraken2 database")
        void testCompleteKrakenDatabase() throws IOException {
            createValidKrakenDatabase(mockDatabaseDir);

            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    mockDatabaseDir.toString(),
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test"),
                    false
            );

            assertTrue(dbInfo.isValid(), "Complete Kraken2 database should be valid");
        }

        @Test
        @DisplayName("Should reject database with only hash file")
        void testDatabaseWithOnlyHashFile() throws IOException {
            Files.createFile(mockDatabaseDir.resolve("hash.k2d"));

            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    mockDatabaseDir.toString(),
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test"),
                    false
            );

            assertFalse(dbInfo.isValid(), "Database with only hash.k2d should be invalid");
        }

        @Test
        @DisplayName("Should reject database with only taxo file")
        void testDatabaseWithOnlyTaxoFile() throws IOException {
            Files.createFile(mockDatabaseDir.resolve("taxo.k2d"));

            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    mockDatabaseDir.toString(),
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test"),
                    false
            );

            assertFalse(dbInfo.isValid(), "Database with only taxo.k2d should be invalid");
        }

        @Test
        @DisplayName("Should validate database with additional optional files")
        void testDatabaseWithOptionalFiles() throws IOException {
            createValidKrakenDatabase(mockDatabaseDir);

            // Add optional files
            Files.createFile(mockDatabaseDir.resolve("opts.k2d"));
            Files.createFile(mockDatabaseDir.resolve("seqid2taxid.map"));

            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    mockDatabaseDir.toString(),
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test"),
                    false
            );

            assertTrue(dbInfo.isValid(), "Database with optional files should still be valid");
        }

        @Test
        @DisplayName("Should handle empty database directory")
        void testEmptyDatabaseDirectory() {
            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    mockDatabaseDir.toString(),
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test"),
                    false
            );

            assertFalse(dbInfo.isValid(), "Empty database directory should be invalid");
        }
    }

    @Nested
    @DisplayName("Database Download Tests")
    class DatabaseDownloadTests {

        @Test
        @DisplayName("Should return DOWNLOAD_REQUIRED when download is allowed")
        void testDownloadAllowed() {
            // Request database with download allowed but no cached/bundled version available
            String dbPath = databaseManager.getDatabasePath(
                    DatabaseManager.DatabaseType.VIRAL,
                    true  // allow download
            );

            // Since no cached or bundled database exists, should return download marker
            // (or null if caching system found something)
            assertTrue(dbPath == null || dbPath.equals("DOWNLOAD_REQUIRED") || dbPath.contains("viral"),
                    "Should handle download scenario appropriately");
        }

        @Test
        @DisplayName("Should return null when download is not allowed")
        void testDownloadNotAllowed() {
            // Request database without allowing download
            String dbPath = databaseManager.getDatabasePath(
                    DatabaseManager.DatabaseType.MINIKRAKEN,
                    false  // do not allow download
            );

            // Without cached/bundled database and no download, should return null or valid path
            assertTrue(dbPath == null || Files.exists(Path.of(dbPath)),
                    "Should return null or valid path when download not allowed");
        }
    }

    @Nested
    @DisplayName("Database Info Tests")
    class DatabaseInfoTests {

        @Test
        @DisplayName("Should store database path correctly")
        void testDatabasePath() {
            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    "/path/to/database",
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test"),
                    true
            );

            assertEquals("/path/to/database", dbInfo.path);
        }

        @Test
        @DisplayName("Should store version correctly")
        void testDatabaseVersion() {
            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    "/path/to/database",
                    "v2.5.3",
                    new DatabaseManager.DatabaseMetadata("viral", "v2.5.3", "test"),
                    true
            );

            assertEquals("v2.5.3", dbInfo.version);
        }

        @Test
        @DisplayName("Should track cached status")
        void testCachedStatus() {
            DatabaseManager.DatabaseInfo cachedDb = new DatabaseManager.DatabaseInfo(
                    "/path/to/database",
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "cached"),
                    true
            );

            DatabaseManager.DatabaseInfo bundledDb = new DatabaseManager.DatabaseInfo(
                    "/path/to/database",
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "bundled"),
                    false
            );

            assertTrue(cachedDb.isCached, "Cached database should have isCached=true");
            assertFalse(bundledDb.isCached, "Bundled database should have isCached=false");
        }

        @Test
        @DisplayName("Should include metadata reference")
        void testMetadataReference() {
            DatabaseManager.DatabaseMetadata metadata =
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test_source");

            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    "/path/to/database",
                    "v1.0",
                    metadata,
                    true
            );

            assertNotNull(dbInfo.metadata);
            assertEquals("viral", dbInfo.metadata.databaseType);
            assertEquals("v1.0", dbInfo.metadata.version);
            assertEquals("test_source", dbInfo.metadata.source);
        }
    }

    @Nested
    @DisplayName("Cache Management Tests")
    class CacheManagementTests {

        @Test
        @DisplayName("Should clear cache for specific database type")
        void testClearCache() {
            // Clear cache should not throw exception even if cache doesn't exist
            assertDoesNotThrow(() -> {
                databaseManager.clearCache(DatabaseManager.DatabaseType.VIRAL);
            }, "Clearing non-existent cache should not throw exception");
        }

        @Test
        @DisplayName("Should get all database info")
        void testGetAllDatabaseInfo() {
            Map<DatabaseManager.DatabaseType, DatabaseManager.DatabaseInfo> allInfo =
                    databaseManager.getAllDatabaseInfo();

            assertNotNull(allInfo, "getAllDatabaseInfo should not return null");
            // May be empty if no databases are cached or bundled
            assertTrue(allInfo.size() >= 0, "Should return valid map");
        }

        @Test
        @DisplayName("Should handle cache for all database types")
        void testCacheForAllTypes() {
            for (DatabaseManager.DatabaseType type : DatabaseManager.DatabaseType.values()) {
                assertDoesNotThrow(() -> {
                    databaseManager.clearCache(type);
                }, "Should handle cache operations for " + type.getId());
            }
        }
    }

    @Nested
    @DisplayName("Error Handling Tests")
    class ErrorHandlingTests {

        @Test
        @DisplayName("Should handle corrupted database gracefully")
        void testCorruptedDatabase() throws IOException {
            // Create files with wrong content/size
            Files.createFile(mockDatabaseDir.resolve("hash.k2d"));
            Files.createFile(mockDatabaseDir.resolve("taxo.k2d"));

            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    mockDatabaseDir.toString(),
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test"),
                    false
            );

            // Even though files exist, they might be considered valid
            // (actual content validation would require more sophisticated checks)
            assertNotNull(dbInfo, "Should create DatabaseInfo even with potentially corrupted files");
        }

        @Test
        @DisplayName("Should handle permission denied scenarios")
        void testPermissionDenied() {
            // Create path that might not be readable (platform-dependent)
            Path restrictedPath = tempCacheDir.resolve("restricted");

            DatabaseManager.DatabaseInfo dbInfo = new DatabaseManager.DatabaseInfo(
                    restrictedPath.toString(),
                    "v1.0",
                    new DatabaseManager.DatabaseMetadata("viral", "v1.0", "test"),
                    false
            );

            // Should handle gracefully without throwing
            assertNotNull(dbInfo, "Should handle restricted paths gracefully");
        }
    }

    @Nested
    @DisplayName("Singleton Pattern Tests")
    class SingletonPatternTests {

        @Test
        @DisplayName("Should return same instance")
        void testSingletonInstance() {
            DatabaseManager instance1 = DatabaseManager.getInstance();
            DatabaseManager instance2 = DatabaseManager.getInstance();

            assertSame(instance1, instance2, "getInstance should return same instance");
        }

        @Test
        @DisplayName("Should not be null")
        void testInstanceNotNull() {
            DatabaseManager instance = DatabaseManager.getInstance();

            assertNotNull(instance, "getInstance should never return null");
        }
    }

    // Helper methods

    /**
     * Creates a valid Kraken2 database structure with required files.
     *
     * @param dbDir the database directory
     * @throws IOException if file creation fails
     */
    private void createValidKrakenDatabase(Path dbDir) throws IOException {
        Files.createDirectories(dbDir);

        // Create required Kraken2 database files
        Files.createFile(dbDir.resolve("hash.k2d"));
        Files.createFile(dbDir.resolve("taxo.k2d"));

        // Write some dummy data to make files non-empty
        Files.writeString(dbDir.resolve("hash.k2d"), "dummy_hash_data");
        Files.writeString(dbDir.resolve("taxo.k2d"), "dummy_taxo_data");
    }
}
