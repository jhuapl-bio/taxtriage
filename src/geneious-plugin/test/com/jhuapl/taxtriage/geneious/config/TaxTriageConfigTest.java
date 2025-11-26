package com.jhuapl.taxtriage.geneious.config;

import com.jhuapl.taxtriage.geneious.TaxTriageOptions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;

import java.io.File;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Comprehensive tests for TaxTriageConfig class.
 *
 * Tests validation rules, default values, parameter mapping,
 * and configuration building from various sources.
 */
class TaxTriageConfigTest {

    private File tempOutputDir;

    @BeforeEach
    void setUp() {
        tempOutputDir = new File(System.getProperty("java.io.tmpdir"), "taxtriage_test");
    }

    @Nested
    @DisplayName("Configuration Building")
    class ConfigurationBuilding {

        @Test
        @DisplayName("Should create default configuration successfully")
        void shouldCreateDefaultConfiguration() {
            TaxTriageConfig config = TaxTriageConfig.createDefault();

            assertNotNull(config);
            assertEquals("ILLUMINA_PE", config.getPreset());
            assertEquals(20, config.getQualityThreshold());
            assertEquals(50, config.getMinReadLength());
            assertEquals(4, config.getThreadCount());
            assertEquals(8, config.getMemoryLimitGb());
            assertNotNull(config.getOutputDirectory());
        }

        @Test
        @DisplayName("Should build configuration from TaxTriageOptions")
        void shouldBuildFromTaxTriageOptions() {
            TaxTriageOptions options = new TaxTriageOptions();

            TaxTriageConfig config = TaxTriageConfig.fromOptions(options);

            assertNotNull(config);
            assertEquals(options.getSequencingPreset().name(), config.getPreset());
            assertEquals(options.getQualityThreshold(), config.getQualityThreshold());
            assertEquals(options.getMinReadLength(), config.getMinReadLength());
            assertEquals(options.getThreadCount(), config.getThreadCount());
            assertEquals(options.getMemoryLimit(), config.getMemoryLimitGb());
        }

        @Test
        @DisplayName("Should build configuration with custom parameters")
        void shouldBuildWithCustomParameters() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ONT")
                    .withQualityThreshold(15)
                    .withMinReadLength(1000)
                    .withSubsampling(true, 50000)
                    .withComputeResources(8, 32)
                    .withDatabases("/custom/kraken", "/custom/bracken")
                    .withOutputDirectory(tempOutputDir)
                    .withTopTaxa(30)
                    .withConfidenceThreshold(0.05)
                    .withReporting(true, false)
                    .withDockerProfile("singularity")
                    .build();

            assertEquals("ONT", config.getPreset());
            assertEquals(15, config.getQualityThreshold());
            assertEquals(1000, config.getMinReadLength());
            assertTrue(config.isEnableSubsampling());
            assertEquals(50000, config.getSubsampleSize());
            assertEquals(8, config.getThreadCount());
            assertEquals(32, config.getMemoryLimitGb());
            assertEquals("/custom/kraken", config.getKrakenDatabase());
            assertEquals("/custom/bracken", config.getBrackenDatabase());
            assertEquals(tempOutputDir, config.getOutputDirectory());
            assertEquals(30, config.getTopTaxa());
            assertEquals(0.05, config.getConfidenceThreshold(), 0.001);
            assertTrue(config.isEnableKrona());
            assertFalse(config.isEnableMultiQC());
            assertEquals("singularity", config.getDockerProfile());
        }
    }

    @Nested
    @DisplayName("Validation")
    class Validation {

        @Test
        @DisplayName("Should validate valid configuration")
        void shouldValidateValidConfiguration() {
            TaxTriageConfig config = TaxTriageConfig.createDefault();

            String result = config.validate();

            assertNull(result);
        }

        @Test
        @DisplayName("Should reject null preset")
        void shouldRejectNullPreset() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset(null)
                    .withDefaults()
                    .build();

            String result = config.validate();

            assertNotNull(result);
            assertTrue(result.contains("Preset must be specified"));
        }

        @Test
        @DisplayName("Should reject empty preset")
        void shouldRejectEmptyPreset() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("")
                    .withDefaults()
                    .build();

            String result = config.validate();

            assertNotNull(result);
            assertTrue(result.contains("Preset must be specified"));
        }

        @Test
        @DisplayName("Should reject invalid quality threshold")
        void shouldRejectInvalidQualityThreshold() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withQualityThreshold(0)
                    .withDefaults()
                    .build();

            String result = config.validate();

            assertNotNull(result);
            assertTrue(result.contains("Quality threshold must be between 1 and 40"));
        }

        @Test
        @DisplayName("Should reject quality threshold too high")
        void shouldRejectQualityThresholdTooHigh() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withQualityThreshold(50)
                    .withDefaults()
                    .build();

            String result = config.validate();

            assertNotNull(result);
            assertTrue(result.contains("Quality threshold must be between 1 and 40"));
        }

        @Test
        @DisplayName("Should reject invalid min read length")
        void shouldRejectInvalidMinReadLength() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withMinReadLength(0)
                    .withDefaults()
                    .build();

            String result = config.validate();

            assertNotNull(result);
            assertTrue(result.contains("Minimum read length must be at least 1"));
        }

        @Test
        @DisplayName("Should reject small subsample size when subsampling enabled")
        void shouldRejectSmallSubsampleSize() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withSubsampling(true, 500)
                    .withDefaults()
                    .build();

            String result = config.validate();

            assertNotNull(result);
            assertTrue(result.contains("Subsample size must be at least 1000 reads"));
        }

        @Test
        @DisplayName("Should accept small subsample size when subsampling disabled")
        void shouldAcceptSmallSubsampleSizeWhenDisabled() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withSubsampling(false, 500)
                    .withDefaults()
                    .build();

            String result = config.validate();

            assertNull(result);
        }

        @Test
        @DisplayName("Should reject invalid thread count")
        void shouldRejectInvalidThreadCount() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withComputeResources(0, 8)
                    .withDefaults()
                    .build();

            String result = config.validate();

            assertNotNull(result);
            assertTrue(result.contains("Thread count must be between 1 and 128"));
        }

        @Test
        @DisplayName("Should reject thread count too high")
        void shouldRejectThreadCountTooHigh() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withComputeResources(200, 8)
                    .withDefaults()
                    .build();

            String result = config.validate();

            assertNotNull(result);
            assertTrue(result.contains("Thread count must be between 1 and 128"));
        }

        @Test
        @DisplayName("Should reject invalid memory limit")
        void shouldRejectInvalidMemoryLimit() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withComputeResources(4, 0)
                    .withDefaults()
                    .build();

            String result = config.validate();

            assertNotNull(result);
            assertTrue(result.contains("Memory limit must be between 1 and 1024 GB"));
        }

        @Test
        @DisplayName("Should reject null databases")
        void shouldRejectNullDatabases() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withDatabases(null, "standard")
                    .withDefaults()
                    .build();

            String result = config.validate();

            assertNotNull(result);
            assertTrue(result.contains("Kraken database must be specified"));
        }

        @Test
        @DisplayName("Should reject null output directory")
        void shouldRejectNullOutputDirectory() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withOutputDirectory(null)
                    .build();

            String result = config.validate();

            assertNotNull(result);
            assertTrue(result.contains("Output directory must be specified"));
        }

        @Test
        @DisplayName("Should reject invalid top taxa count")
        void shouldRejectInvalidTopTaxaCount() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withTopTaxa(0)
                    .withDefaults()
                    .build();

            String result = config.validate();

            assertNotNull(result);
            assertTrue(result.contains("Top taxa count must be between 1 and 1000"));
        }

        @Test
        @DisplayName("Should reject invalid confidence threshold")
        void shouldRejectInvalidConfidenceThreshold() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withConfidenceThreshold(-0.1)
                    .withDefaults()
                    .build();

            String result = config.validate();

            assertNotNull(result);
            assertTrue(result.contains("Confidence threshold must be between 0.0 and 1.0"));
        }
    }

    @Nested
    @DisplayName("Parameter Mapping")
    class ParameterMapping {

        @Test
        @DisplayName("Should convert to parameter map correctly")
        void shouldConvertToParameterMap() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ONT")
                    .withQualityThreshold(15)
                    .withMinReadLength(1000)
                    .withSubsampling(true, 50000)
                    .withComputeResources(8, 32)
                    .withDatabases("/kraken/db", "/bracken/db")
                    .withOutputDirectory(tempOutputDir)
                    .withTopTaxa(30)
                    .withConfidenceThreshold(0.05)
                    .withReporting(true, false)
                    .build();

            Map<String, Object> params = config.toParameterMap();

            assertEquals("ONT", params.get("preset"));
            assertEquals(15, params.get("quality_threshold"));
            assertEquals(1000, params.get("min_read_length"));
            assertEquals(true, params.get("enable_subsampling"));
            assertEquals(50000, params.get("subsample_size"));
            assertEquals(8, params.get("max_cpus"));
            assertEquals("32.GB", params.get("max_memory"));
            assertEquals("/kraken/db", params.get("kraken_db"));
            assertEquals("/bracken/db", params.get("bracken_db"));
            assertEquals(tempOutputDir.getAbsolutePath(), params.get("outdir"));
            assertEquals(30, params.get("top_taxa"));
            assertEquals(0.05, params.get("confidence_threshold"));
            assertEquals(true, params.get("enable_krona"));
            assertEquals(false, params.get("enable_multiqc"));
        }
    }

    @Nested
    @DisplayName("Preset Detection")
    class PresetDetection {

        @Test
        @DisplayName("Should detect ONT as long read preset")
        void shouldDetectOntAsLongRead() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ONT")
                    .withDefaults()
                    .build();

            assertTrue(config.isLongReadPreset());
            assertFalse(config.isPairedEndPreset());
        }

        @Test
        @DisplayName("Should detect Illumina PE as paired end preset")
        void shouldDetectIlluminaPeAsPairedEnd() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withDefaults()
                    .build();

            assertFalse(config.isLongReadPreset());
            assertTrue(config.isPairedEndPreset());
        }

        @Test
        @DisplayName("Should detect Illumina SE as neither long read nor paired end")
        void shouldDetectIlluminaSeAsNeither() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_SE")
                    .withDefaults()
                    .build();

            assertFalse(config.isLongReadPreset());
            assertFalse(config.isPairedEndPreset());
        }
    }

    @Nested
    @DisplayName("Docker Profile")
    class DockerProfile {

        @Test
        @DisplayName("Should return default docker profile when none specified")
        void shouldReturnDefaultDockerProfile() {
            TaxTriageConfig config = TaxTriageConfig.createDefault();

            assertEquals("docker", config.getDockerProfile());
        }

        @Test
        @DisplayName("Should return custom docker profile when specified")
        void shouldReturnCustomDockerProfile() {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withDockerProfile("singularity")
                    .withDefaults()
                    .build();

            assertEquals("singularity", config.getDockerProfile());
        }
    }

    @Nested
    @DisplayName("Equality and HashCode")
    class EqualityAndHashCode {

        @Test
        @DisplayName("Should be equal when all parameters match")
        void shouldBeEqualWhenAllParametersMatch() {
            TaxTriageConfig config1 = TaxTriageConfig.createDefault();
            TaxTriageConfig config2 = TaxTriageConfig.createDefault();

            assertEquals(config1, config2);
            assertEquals(config1.hashCode(), config2.hashCode());
        }

        @Test
        @DisplayName("Should not be equal when parameters differ")
        void shouldNotBeEqualWhenParametersDiffer() {
            TaxTriageConfig config1 = TaxTriageConfig.createDefault();
            TaxTriageConfig config2 = new TaxTriageConfig.Builder()
                    .withPreset("ONT")
                    .withDefaults()
                    .build();

            assertNotEquals(config1, config2);
        }

        @Test
        @DisplayName("Should handle null comparison")
        void shouldHandleNullComparison() {
            TaxTriageConfig config = TaxTriageConfig.createDefault();

            assertNotEquals(config, null);
        }

        @Test
        @DisplayName("Should handle different class comparison")
        void shouldHandleDifferentClassComparison() {
            TaxTriageConfig config = TaxTriageConfig.createDefault();

            assertNotEquals(config, "not a config");
        }
    }

    @Nested
    @DisplayName("toString")
    class ToStringTest {

        @Test
        @DisplayName("Should provide meaningful string representation")
        void shouldProvideMeaningfulStringRepresentation() {
            TaxTriageConfig config = TaxTriageConfig.createDefault();

            String str = config.toString();

            assertNotNull(str);
            assertTrue(str.contains("TaxTriageConfig"));
            assertTrue(str.contains("preset"));
            assertTrue(str.contains("qualityThreshold"));
            assertTrue(str.contains("ILLUMINA_PE"));
        }
    }
}