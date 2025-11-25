package com.jhuapl.taxtriage.geneious.utils;

import com.jhuapl.taxtriage.geneious.database.DatabaseManager.DatabaseType;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for DatabasePathUtil (Group 3).
 *
 * Tests database path utility functionality including:
 * - Database path resolution
 * - Database type mapping
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class DatabasePathUtilTest {

    /**
     * GROUP 3: Test get database path
     */
    @Test
    void testGetDatabasePath() {
        String path = DatabasePathUtil.getDatabasePath("standard");

        assertNotNull(path);
        assertTrue(path.contains(".taxtriage-geneious"));
        assertTrue(path.contains("standard"));
    }

    /**
     * GROUP 3: Test get database path with null name
     */
    @Test
    void testGetDatabasePathWithNullName() {
        String path = DatabasePathUtil.getDatabasePath(null);

        assertNotNull(path);
        assertTrue(path.contains(".taxtriage-geneious"));
        assertTrue(path.contains("standard"), "Should default to standard");
    }

    /**
     * GROUP 3: Test get database path with empty name
     */
    @Test
    void testGetDatabasePathWithEmptyName() {
        String path = DatabasePathUtil.getDatabasePath("");

        assertNotNull(path);
        assertTrue(path.contains("standard"), "Should default to standard");
    }

    /**
     * GROUP 3: Test map database name to type - viral
     */
    @Test
    void testMapDatabaseNameToTypeViral() {
        assertEquals(DatabaseType.VIRAL,
            DatabasePathUtil.mapDatabaseNameToType("viral"));
        assertEquals(DatabaseType.VIRAL,
            DatabasePathUtil.mapDatabaseNameToType("Viral"));
        assertEquals(DatabaseType.VIRAL,
            DatabasePathUtil.mapDatabaseNameToType("kraken2-viral"));
    }

    /**
     * GROUP 3: Test map database name to type - standard
     */
    @Test
    void testMapDatabaseNameToTypeStandard() {
        assertEquals(DatabaseType.STANDARD,
            DatabasePathUtil.mapDatabaseNameToType("standard"));
        assertEquals(DatabaseType.STANDARD,
            DatabasePathUtil.mapDatabaseNameToType("Standard"));
        assertEquals(DatabaseType.STANDARD,
            DatabasePathUtil.mapDatabaseNameToType("kraken2-standard"));
    }

    /**
     * GROUP 3: Test map database name to type - minikraken
     */
    @Test
    void testMapDatabaseNameToTypeMiniKraken() {
        assertEquals(DatabaseType.MINIKRAKEN,
            DatabasePathUtil.mapDatabaseNameToType("minikraken"));
        assertEquals(DatabaseType.MINIKRAKEN,
            DatabasePathUtil.mapDatabaseNameToType("MiniKraken"));
        assertEquals(DatabaseType.MINIKRAKEN,
            DatabasePathUtil.mapDatabaseNameToType("mini"));
    }

    /**
     * GROUP 3: Test map database name to type - null
     */
    @Test
    void testMapDatabaseNameToTypeNull() {
        assertNull(DatabasePathUtil.mapDatabaseNameToType(null));
    }

    /**
     * GROUP 3: Test map database name to type - unknown
     */
    @Test
    void testMapDatabaseNameToTypeUnknown() {
        assertNull(DatabasePathUtil.mapDatabaseNameToType("unknown"));
        assertNull(DatabasePathUtil.mapDatabaseNameToType("custom"));
        assertNull(DatabasePathUtil.mapDatabaseNameToType("xyz"));
    }

    /**
     * GROUP 3: Test case insensitive mapping
     */
    @Test
    void testCaseInsensitiveMapping() {
        assertEquals(DatabaseType.VIRAL,
            DatabasePathUtil.mapDatabaseNameToType("VIRAL"));
        assertEquals(DatabaseType.STANDARD,
            DatabasePathUtil.mapDatabaseNameToType("STANDARD"));
        assertEquals(DatabaseType.MINIKRAKEN,
            DatabasePathUtil.mapDatabaseNameToType("MINIKRAKEN"));
    }

    /**
     * GROUP 3: Test trimming in mapping
     */
    @Test
    void testTrimmingInMapping() {
        assertEquals(DatabaseType.VIRAL,
            DatabasePathUtil.mapDatabaseNameToType("  viral  "));
        assertEquals(DatabaseType.STANDARD,
            DatabasePathUtil.mapDatabaseNameToType("  standard  "));
    }

    /**
     * GROUP 3: Test different database names
     */
    @Test
    void testDifferentDatabaseNames() {
        String viralPath = DatabasePathUtil.getDatabasePath("viral");
        String standardPath = DatabasePathUtil.getDatabasePath("standard");
        String customPath = DatabasePathUtil.getDatabasePath("custom");

        assertNotNull(viralPath);
        assertNotNull(standardPath);
        assertNotNull(customPath);

        assertTrue(viralPath.contains("viral"));
        assertTrue(standardPath.contains("standard"));
        assertTrue(customPath.contains("custom"));

        assertNotEquals(viralPath, standardPath);
        assertNotEquals(standardPath, customPath);
    }
}
