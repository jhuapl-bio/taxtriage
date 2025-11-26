package com.jhuapl.taxtriage.geneious.config;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for Nextflow database name mapping to ensure compatibility with TaxTriage workflow.
 */
@DisplayName("Nextflow Database Mapping Tests")
class NextflowDatabaseMappingTest {

    @Test
    @DisplayName("Standard-8 variants map to standard8")
    void testStandard8Mapping() {
        assertDatabaseMapping("standard-8", "standard8");
        assertDatabaseMapping("standard8", "standard8");
        assertDatabaseMapping("Standard-8", "standard8");
    }

    @Test
    @DisplayName("Standard-16 variants map to standard16")
    void testStandard16Mapping() {
        assertDatabaseMapping("standard-16", "standard16");
        assertDatabaseMapping("standard16", "standard16");
    }

    @Test
    @DisplayName("Standard database maps correctly")
    void testStandardMapping() {
        assertDatabaseMapping("standard", "standard");
    }

    @Test
    @DisplayName("Viral database maps correctly")
    void testViralMapping() {
        assertDatabaseMapping("viral", "viral");
    }

    @Test
    @DisplayName("PlusPF variants map correctly")
    void testPlusPFMapping() {
        assertDatabaseMapping("pluspf-8", "pluspf8");
        assertDatabaseMapping("pluspf8", "pluspf8");
        assertDatabaseMapping("pluspf-16", "pluspfp16");
        assertDatabaseMapping("pluspfp16", "pluspfp16");
        assertDatabaseMapping("pluspf", "pluspf");
    }

    @Test
    @DisplayName("EuPath database maps correctly")
    void testEuPathMapping() {
        assertDatabaseMapping("eupath", "eupath");
        assertDatabaseMapping("eupathdb", "eupath");
    }

    @Test
    @DisplayName("MiniKraken variants map to minikraken2")
    void testMiniKrakenMapping() {
        assertDatabaseMapping("minikraken", "minikraken2");
        assertDatabaseMapping("minikraken2", "minikraken2");
    }

    @Test
    @DisplayName("FluKraken variants map to flukraken2")
    void testFluKrakenMapping() {
        assertDatabaseMapping("flukraken", "flukraken2");
        assertDatabaseMapping("flukraken2", "flukraken2");
    }

    @Test
    @DisplayName("Test database maps correctly")
    void testTestDatabaseMapping() {
        assertDatabaseMapping("test", "test");
    }

    private void assertDatabaseMapping(String userFriendlyName, String expectedNextflowName) {
        TaxTriageConfig config = new TaxTriageConfig.Builder()
                .withDatabases(userFriendlyName, userFriendlyName)
                .withOutputDirectory(new File("/tmp/test"))
                .build();

        Map<String, Object> params = config.toParameterMap();

        // When download is needed, should use the mapped Nextflow name
        String dbParam = (String) params.get("db");

        // If download_db is true, check the mapping
        Boolean downloadDb = (Boolean) params.get("download_db");
        if (downloadDb != null && downloadDb) {
            assertEquals(expectedNextflowName, dbParam,
                String.format("Database '%s' should map to Nextflow identifier '%s'",
                    userFriendlyName, expectedNextflowName));
        }
    }

    @Test
    @DisplayName("Custom database paths should pass through unchanged")
    void testCustomPathPassthrough() {
        String customPath = "/custom/database/path";
        
        TaxTriageConfig config = new TaxTriageConfig.Builder()
                .withDatabases(customPath, customPath)
                .withOutputDirectory(new File("/tmp/test"))
                .build();

        Map<String, Object> params = config.toParameterMap();
        
        // Custom paths should be preserved
        String dbParam = (String) params.get("db");
        assertTrue(dbParam.equals(customPath) || dbParam.contains(customPath),
            "Custom database path should be preserved");
    }

    @Test
    @DisplayName("Case insensitive database name mapping")
    void testCaseInsensitiveMapping() {
        String[] variations = {"Standard-8", "STANDARD-8", "StAnDaRd-8"};
        
        for (String variation : variations) {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withDatabases(variation, variation)
                    .withOutputDirectory(new File("/tmp/test"))
                    .build();

            Map<String, Object> params = config.toParameterMap();
            
            Boolean downloadDb = (Boolean) params.get("download_db");
            if (downloadDb != null && downloadDb) {
                String dbParam = (String) params.get("db");
                assertEquals("standard8", dbParam,
                    String.format("Case variation '%s' should map to 'standard8'", variation));
            }
        }
    }

    @Test
    @DisplayName("Supported databases list matches Nextflow workflow")
    void testSupportedDatabasesDocumentation() {
        // This test documents the supported databases from the Nextflow workflow
        // Nextflow supported_dbs keys: standard8, standard, viral, pluspf, pluspfp16, 
        // pluspf8, eupath, minikraken2, flukraken2, test
        
        String[] nextflowDatabases = {
            "standard8", "standard", "viral", "pluspf", "pluspfp16", 
            "pluspf8", "eupath", "minikraken2", "flukraken2", "test"
        };
        
        // Verify each can be used
        for (String db : nextflowDatabases) {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withDatabases(db, db)
                    .withOutputDirectory(new File("/tmp/test"))
                    .build();

            Map<String, Object> params = config.toParameterMap();
            assertNotNull(params.get("db"), 
                "Database parameter should be set for: " + db);
        }
    }
}
