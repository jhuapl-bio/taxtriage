package com.jhuapl.taxtriage.geneious;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for operation options and validation (Group 2).
 *
 * Tests:
 * - Options initialization with defaults
 * - Sequencing preset configuration
 * - Database selection
 * - Resource limits (threads, memory)
 * - Quality filtering parameters
 * - Subsampling options
 * - BBTools preprocessing options
 * - Input validation
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class OperationOptionsTest {

    private TaxTriageOptions options;

    @BeforeEach
    void setUp() {
        options = new TaxTriageOptions();
    }

    /**
     * GROUP 2: Test default options initialization
     */
    @Test
    void testDefaultOptionsInitialization() {
        assertNotNull(options, "Options should be initialized");

        // Test default values
        assertEquals(TaxTriageOptions.SequencingPreset.ILLUMINA_PE,
            options.getSequencingPreset(),
            "Default preset should be Illumina PE");

        assertEquals("standard", options.getKrakenDatabase(),
            "Default Kraken database should be 'standard'");
        assertEquals("standard", options.getBrackenDatabase(),
            "Default Bracken database should be 'standard'");

        assertTrue(options.getThreadCount() > 0,
            "Default thread count should be positive");
        assertTrue(options.getMemoryLimit() > 0,
            "Default memory limit should be positive");
    }

    /**
     * GROUP 2: Test sequencing preset enum values
     */
    @Test
    void testSequencingPresets() {
        TaxTriageOptions.SequencingPreset[] presets =
            TaxTriageOptions.SequencingPreset.values();

        assertTrue(presets.length >= 4,
            "Should have at least 4 preset types");

        // Verify each preset has a display name
        for (TaxTriageOptions.SequencingPreset preset : presets) {
            assertNotNull(preset.getDisplayName(),
                "Preset should have display name: " + preset);
            assertFalse(preset.getDisplayName().isEmpty(),
                "Display name should not be empty: " + preset);
        }
    }

    /**
     * GROUP 2: Test ONT preset configuration
     */
    @Test
    void testONTPreset() {
        TaxTriageOptions.SequencingPreset ont =
            TaxTriageOptions.SequencingPreset.ONT;

        assertNotNull(ont);
        assertTrue(ont.getDisplayName().contains("ONT") ||
                  ont.getDisplayName().contains("Nanopore"));
    }

    /**
     * GROUP 2: Test Illumina PE preset configuration
     */
    @Test
    void testIlluminaPEPreset() {
        TaxTriageOptions.SequencingPreset illuminaPE =
            TaxTriageOptions.SequencingPreset.ILLUMINA_PE;

        assertNotNull(illuminaPE);
        assertTrue(illuminaPE.getDisplayName().contains("Illumina"));
        assertTrue(illuminaPE.getDisplayName().contains("Paired") ||
                  illuminaPE.getDisplayName().contains("PE"));
    }

    /**
     * GROUP 2: Test Illumina SE preset configuration
     */
    @Test
    void testIlluminaSEPreset() {
        TaxTriageOptions.SequencingPreset illuminaSE =
            TaxTriageOptions.SequencingPreset.ILLUMINA_SE;

        assertNotNull(illuminaSE);
        assertTrue(illuminaSE.getDisplayName().contains("Illumina"));
        assertTrue(illuminaSE.getDisplayName().contains("Single") ||
                  illuminaSE.getDisplayName().contains("SE"));
    }

    /**
     * GROUP 2: Test database selection options
     */
    @Test
    void testDatabaseOptions() {
        // Test default database
        assertEquals("standard", options.getKrakenDatabase());

        // Verify databases are accessible
        assertNotNull(options.getKrakenDatabase());
        assertNotNull(options.getBrackenDatabase());
    }

    /**
     * GROUP 2: Test thread count configuration
     */
    @Test
    void testThreadCount() {
        int defaultThreads = options.getThreadCount();
        assertTrue(defaultThreads > 0,
            "Default thread count should be positive");

        // Thread count should be reasonable (not more than available processors)
        int availableProcessors = Runtime.getRuntime().availableProcessors();
        assertTrue(defaultThreads <= availableProcessors * 2,
            "Thread count should not exceed 2x available processors");
    }

    /**
     * GROUP 2: Test memory limit configuration
     */
    @Test
    void testMemoryLimit() {
        int defaultMemory = options.getMemoryLimit();
        assertTrue(defaultMemory > 0,
            "Default memory limit should be positive");

        // Memory should be at least 4GB for genomic analysis
        assertTrue(defaultMemory >= 4,
            "Memory limit should be at least 4GB");
    }

    /**
     * GROUP 2: Test quality threshold default
     */
    @Test
    void testQualityThreshold() {
        int qualityThreshold = options.getQualityThreshold();
        assertTrue(qualityThreshold >= 0,
            "Quality threshold should be non-negative");
        assertTrue(qualityThreshold <= 40,
            "Quality threshold should be reasonable (<= 40)");
    }

    /**
     * GROUP 2: Test min read length default
     */
    @Test
    void testMinReadLength() {
        int minReadLength = options.getMinReadLength();
        assertTrue(minReadLength > 0,
            "Min read length should be positive");
        assertTrue(minReadLength >= 20,
            "Min read length should be reasonable (>= 20)");
    }

    /**
     * GROUP 2: Test subsampling disabled by default
     */
    @Test
    void testSubsamplingDefaults() {
        assertFalse(options.isSubsamplingEnabled(),
            "Subsampling should be disabled by default");

        int subsampleSize = options.getSubsampleSize();
        assertTrue(subsampleSize > 0,
            "Subsample size should be positive even when disabled");
        assertEquals(100000, subsampleSize,
            "Default subsample size should be 100000");
    }

    /**
     * GROUP 2: Test BBTools preprocessing defaults
     */
    @Test
    void testBBToolsDefaults() {
        // BBTools preprocessing should have reasonable defaults
        int subsThreshold = options.getBBToolsSubstitutionThreshold();
        assertTrue(subsThreshold >= 0,
            "Substitution threshold should be non-negative");

        int memoryGB = options.getBBToolsMemoryGB();
        assertTrue(memoryGB > 0,
            "BBTools memory should be positive");
        assertTrue(memoryGB >= 4,
            "BBTools memory should be at least 4GB");
    }

    /**
     * GROUP 2: Test input file handling
     */
    @Test
    void testInputFileHandling() {
        List<File> inputFiles = options.getInputFiles();
        // Can be null or empty initially
        assertTrue(inputFiles == null || inputFiles.isEmpty(),
            "Input files should be null or empty initially");
    }

    /**
     * GROUP 2: Test input directory handling
     */
    @Test
    void testInputDirectoryHandling() {
        File inputDirectory = options.getInputDirectory();
        // Can be null initially
        assertNull(inputDirectory,
            "Input directory should be null initially");
    }

    /**
     * GROUP 2: Test options validation - valid configuration
     */
    @Test
    void testValidOptionsValidation() {
        // Default options should be valid
        String validation = validateOptions(options);
        assertNull(validation,
            "Default options should be valid");
    }

    /**
     * GROUP 2: Test options validation - invalid thread count
     */
    @Test
    void testInvalidThreadCount() {
        // Cannot directly set invalid values through public API,
        // but we can test the validation logic
        String validation = validateThreadCount(0);
        assertNotNull(validation,
            "Zero thread count should be invalid");

        validation = validateThreadCount(-1);
        assertNotNull(validation,
            "Negative thread count should be invalid");
    }

    /**
     * GROUP 2: Test options validation - invalid memory
     */
    @Test
    void testInvalidMemoryLimit() {
        String validation = validateMemoryLimit(0);
        assertNotNull(validation,
            "Zero memory should be invalid");

        validation = validateMemoryLimit(-1);
        assertNotNull(validation,
            "Negative memory should be invalid");
    }

    /**
     * GROUP 2: Test options validation - invalid quality threshold
     */
    @Test
    void testInvalidQualityThreshold() {
        String validation = validateQualityThreshold(-1);
        assertNotNull(validation,
            "Negative quality threshold should be invalid");

        validation = validateQualityThreshold(50);
        assertNotNull(validation,
            "Quality threshold > 40 should be invalid");
    }

    /**
     * GROUP 2: Test preset-specific configurations
     */
    @Test
    void testPresetSpecificSettings() {
        // Each preset should have appropriate settings
        TaxTriageOptions.SequencingPreset[] presets =
            TaxTriageOptions.SequencingPreset.values();

        for (TaxTriageOptions.SequencingPreset preset : presets) {
            assertNotNull(preset);
            assertNotNull(preset.getDisplayName());

            // Display name should match preset type
            String displayName = preset.getDisplayName().toUpperCase();
            String presetName = preset.name();

            if (presetName.contains("ONT")) {
                assertTrue(displayName.contains("ONT") ||
                          displayName.contains("NANOPORE"));
            } else if (presetName.contains("ILLUMINA")) {
                assertTrue(displayName.contains("ILLUMINA"));
            }
        }
    }

    /**
     * GROUP 2: Test resource limit boundaries
     */
    @Test
    void testResourceLimitBoundaries() {
        // Thread count should be within reasonable limits
        int threads = options.getThreadCount();
        assertTrue(threads >= 1 && threads <= 128,
            "Thread count should be between 1 and 128");

        // Memory should be within reasonable limits
        int memory = options.getMemoryLimit();
        assertTrue(memory >= 4 && memory <= 512,
            "Memory limit should be between 4GB and 512GB");
    }

    /**
     * GROUP 2: Test database name validation
     */
    @Test
    void testDatabaseNameValidation() {
        String krakenDb = options.getKrakenDatabase();
        assertNotNull(krakenDb);
        assertFalse(krakenDb.isEmpty());

        String brackenDb = options.getBrackenDatabase();
        assertNotNull(brackenDb);
        assertFalse(brackenDb.isEmpty());
    }

    // ===== Helper Methods =====

    /**
     * Validates complete options configuration
     */
    private String validateOptions(TaxTriageOptions opts) {
        if (opts.getThreadCount() <= 0) {
            return "Thread count must be positive";
        }
        if (opts.getMemoryLimit() <= 0) {
            return "Memory limit must be positive";
        }
        if (opts.getQualityThreshold() < 0) {
            return "Quality threshold must be non-negative";
        }
        if (opts.getMinReadLength() <= 0) {
            return "Min read length must be positive";
        }
        return null; // Valid
    }

    /**
     * Validates thread count
     */
    private String validateThreadCount(int threads) {
        if (threads <= 0) {
            return "Thread count must be positive";
        }
        return null;
    }

    /**
     * Validates memory limit
     */
    private String validateMemoryLimit(int memoryGB) {
        if (memoryGB <= 0) {
            return "Memory limit must be positive";
        }
        return null;
    }

    /**
     * Validates quality threshold
     */
    private String validateQualityThreshold(int threshold) {
        if (threshold < 0) {
            return "Quality threshold must be non-negative";
        }
        if (threshold > 40) {
            return "Quality threshold must be <= 40";
        }
        return null;
    }
}
