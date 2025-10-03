package com.jhuapl.taxtriage.geneious.config;

import com.jhuapl.taxtriage.geneious.TaxTriageOptions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Comprehensive tests for ConfigGenerator class.
 *
 * Tests configuration generation for different presets,
 * parameter substitution, path handling, and error conditions.
 */
class ConfigGeneratorTest {

    @TempDir
    Path tempDir;

    private ConfigGenerator configGenerator;
    private MockTemplateLoader mockTemplateLoader;

    @BeforeEach
    void setUp() {
        mockTemplateLoader = new MockTemplateLoader();
        configGenerator = new ConfigGenerator(mockTemplateLoader);
    }

    @Nested
    @DisplayName("Configuration Generation")
    class ConfigurationGeneration {

        @Test
        @DisplayName("Should generate config from TaxTriageOptions")
        void shouldGenerateConfigFromOptions() throws Exception {
            TaxTriageOptions options = new TaxTriageOptions();
            File outputFile = tempDir.resolve("nextflow.config").toFile();

            configGenerator.generateConfig(options, outputFile);

            assertTrue(outputFile.exists());
            String content = Files.readString(outputFile.toPath());
            assertNotNull(content);
            assertTrue(content.contains("params {"));
            assertTrue(content.contains("preset"));
        }

        @Test
        @DisplayName("Should generate config from TaxTriageConfig")
        void shouldGenerateConfigFromConfig() throws Exception {
            TaxTriageConfig config = TaxTriageConfig.createDefault();
            File outputFile = tempDir.resolve("nextflow.config").toFile();

            configGenerator.generateConfig(config, outputFile);

            assertTrue(outputFile.exists());
            String content = Files.readString(outputFile.toPath());
            assertNotNull(content);
            assertTrue(content.contains("params {"));
        }

        @Test
        @DisplayName("Should generate params.json from config")
        void shouldGenerateParamsFromConfig() throws Exception {
            TaxTriageConfig config = TaxTriageConfig.createDefault();
            File outputFile = tempDir.resolve("params.json").toFile();

            configGenerator.generateParams(config, outputFile);

            assertTrue(outputFile.exists());
            String content = Files.readString(outputFile.toPath());
            assertNotNull(content);
            assertTrue(content.contains("{"));
            assertTrue(content.contains("preset"));
        }

        @Test
        @DisplayName("Should create parent directories if they don't exist")
        void shouldCreateParentDirectories() throws Exception {
            TaxTriageConfig config = TaxTriageConfig.createDefault();
            File outputFile = tempDir.resolve("subdir/nextflow.config").toFile();

            configGenerator.generateConfig(config, outputFile);

            assertTrue(outputFile.exists());
            assertTrue(outputFile.getParentFile().exists());
        }
    }

    @Nested
    @DisplayName("Preset-Specific Generation")
    class PresetSpecificGeneration {

        @Test
        @DisplayName("Should generate ONT preset configuration")
        void shouldGenerateOntPresetConfiguration() throws Exception {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ONT")
                    .withQualityThreshold(10)
                    .withMinReadLength(500)
                    .withComputeResources(8, 32)
                    .withDefaults()
                    .build();
            File outputFile = tempDir.resolve("nextflow.config").toFile();

            configGenerator.generateConfig(config, outputFile);

            String content = Files.readString(outputFile.toPath());
            assertTrue(content.contains("preset = 'ONT'"));
            assertTrue(content.contains("quality_threshold = 10"));
            assertTrue(content.contains("min_read_length = 500"));
            assertTrue(content.contains("max_cpus = 8"));
            assertTrue(content.contains("max_memory = '32.GB'"));
        }

        @Test
        @DisplayName("Should generate Illumina PE preset configuration")
        void shouldGenerateIlluminaPePresetConfiguration() throws Exception {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_PE")
                    .withQualityThreshold(20)
                    .withMinReadLength(50)
                    .withComputeResources(4, 16)
                    .withDefaults()
                    .build();
            File outputFile = tempDir.resolve("nextflow.config").toFile();

            configGenerator.generateConfig(config, outputFile);

            String content = Files.readString(outputFile.toPath());
            assertTrue(content.contains("preset = 'ILLUMINA_PE'"));
            assertTrue(content.contains("quality_threshold = 20"));
            assertTrue(content.contains("min_read_length = 50"));
            assertTrue(content.contains("max_cpus = 4"));
            assertTrue(content.contains("max_memory = '16.GB'"));
        }

        @Test
        @DisplayName("Should generate Illumina SE preset configuration")
        void shouldGenerateIlluminaSePresetConfiguration() throws Exception {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("ILLUMINA_SE")
                    .withQualityThreshold(20)
                    .withMinReadLength(50)
                    .withComputeResources(4, 16)
                    .withDefaults()
                    .build();
            File outputFile = tempDir.resolve("nextflow.config").toFile();

            configGenerator.generateConfig(config, outputFile);

            String content = Files.readString(outputFile.toPath());
            assertTrue(content.contains("preset = 'ILLUMINA_SE'"));
            assertTrue(content.contains("quality_threshold = 20"));
            assertTrue(content.contains("min_read_length = 50"));
        }
    }

    @Nested
    @DisplayName("Parameter Substitution")
    class ParameterSubstitution {

        @Test
        @DisplayName("Should substitute all basic parameters")
        void shouldSubstituteAllBasicParameters() throws Exception {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("TEST")
                    .withQualityThreshold(25)
                    .withMinReadLength(100)
                    .withSubsampling(true, 75000)
                    .withComputeResources(6, 24)
                    .withDatabases("/test/kraken", "/test/bracken")
                    .withOutputDirectory(tempDir.toFile())
                    .withTopTaxa(15)
                    .withConfidenceThreshold(0.15)
                    .withReporting(false, true)
                    .withDockerProfile("singularity")
                    .build();
            File outputFile = tempDir.resolve("nextflow.config").toFile();

            configGenerator.generateConfig(config, outputFile);

            String content = Files.readString(outputFile.toPath());
            assertTrue(content.contains("preset = 'TEST'"));
            assertTrue(content.contains("quality_threshold = 25"));
            assertTrue(content.contains("min_read_length = 100"));
            assertTrue(content.contains("enable_subsampling = true"));
            assertTrue(content.contains("subsample_size = 75000"));
            assertTrue(content.contains("max_cpus = 6"));
            assertTrue(content.contains("max_memory = '24.GB'"));
            assertTrue(content.contains("kraken_db = '/test/kraken'"));
            assertTrue(content.contains("bracken_db = '/test/bracken'"));
            assertTrue(content.contains("top_taxa = 15"));
            assertTrue(content.contains("confidence_threshold = 0.15"));
            assertTrue(content.contains("enable_krona = false"));
            assertTrue(content.contains("enable_multiqc = true"));
            assertTrue(content.contains("docker_profile = singularity"));
        }

        @Test
        @DisplayName("Should substitute parameters in JSON format")
        void shouldSubstituteParametersInJsonFormat() throws Exception {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("TEST")
                    .withQualityThreshold(25)
                    .withSubsampling(true, 75000)
                    .withDefaults()
                    .build();
            File outputFile = tempDir.resolve("params.json").toFile();

            configGenerator.generateParams(config, outputFile);

            String content = Files.readString(outputFile.toPath());
            assertTrue(content.contains("\"preset\": \"TEST\""));
            assertTrue(content.contains("\"quality_threshold\": 25"));
            assertTrue(content.contains("\"enable_subsampling\": true"));
            assertTrue(content.contains("\"subsample_size\": 75000"));
        }
    }

    @Nested
    @DisplayName("Path Handling")
    class PathHandling {

        @Test
        @DisplayName("Should escape paths with spaces for Docker")
        void shouldEscapePathsWithSpacesForDocker() throws Exception {
            File outputDirWithSpaces = new File(tempDir.toFile(), "output with spaces");
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("TEST")
                    .withOutputDirectory(outputDirWithSpaces)
                    .withDatabases("/path with spaces/kraken", "/another path/bracken")
                    .withDefaults()
                    .build();
            File outputFile = tempDir.resolve("nextflow.config").toFile();

            configGenerator.generateConfig(config, outputFile);

            String content = Files.readString(outputFile.toPath());
            assertTrue(content.contains("'/path with spaces/kraken'") ||
                      content.contains("\"/path with spaces/kraken\""));
            assertTrue(content.contains("'/another path/bracken'") ||
                      content.contains("\"/another path/bracken\""));
        }

        @Test
        @DisplayName("Should handle Windows-style paths")
        void shouldHandleWindowsStylePaths() throws Exception {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("TEST")
                    .withDatabases("C:\\Users\\test\\kraken", "D:\\data\\bracken")
                    .withDefaults()
                    .build();
            File outputFile = tempDir.resolve("nextflow.config").toFile();

            configGenerator.generateConfig(config, outputFile);

            String content = Files.readString(outputFile.toPath());
            // Windows paths should be converted to Unix-style for Docker
            assertTrue(content.contains("C:/Users/test/kraken") ||
                      content.contains("C:\\\\Users\\\\test\\\\kraken"));
            assertTrue(content.contains("D:/data/bracken") ||
                      content.contains("D:\\\\data\\\\bracken"));
        }

        @Test
        @DisplayName("Should generate Docker volumes for file paths")
        void shouldGenerateDockerVolumesForFilePaths() throws Exception {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("TEST")
                    .withDatabases("/data/kraken/db.k2d", "/data/bracken/db.bracken")
                    .withOutputDirectory(tempDir.toFile())
                    .withDefaults()
                    .build();
            File outputFile = tempDir.resolve("nextflow.config").toFile();

            configGenerator.generateConfig(config, outputFile);

            String content = Files.readString(outputFile.toPath());
            // Should include volume mounts for database directories
            assertTrue(content.contains(tempDir.toString()));
        }

        @Test
        @DisplayName("Should not generate volumes for standard database names")
        void shouldNotGenerateVolumesForStandardDatabaseNames() throws Exception {
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("TEST")
                    .withDatabases("standard", "standard")
                    .withDefaults()
                    .build();
            File outputFile = tempDir.resolve("nextflow.config").toFile();

            configGenerator.generateConfig(config, outputFile);

            String content = Files.readString(outputFile.toPath());
            assertTrue(content.contains("kraken_db = 'standard'"));
            assertTrue(content.contains("bracken_db = 'standard'"));
        }
    }

    @Nested
    @DisplayName("Error Handling")
    class ErrorHandling {

        @Test
        @DisplayName("Should throw ConfigurationException for invalid config")
        void shouldThrowConfigurationExceptionForInvalidConfig() {
            TaxTriageConfig invalidConfig = new TaxTriageConfig.Builder()
                    .withPreset(null)
                    .build();
            File outputFile = tempDir.resolve("nextflow.config").toFile();

            assertThrows(ConfigGenerator.ConfigurationException.class, () ->
                    configGenerator.generateConfig(invalidConfig, outputFile));
        }

        @Test
        @DisplayName("Should throw ConfigurationException for invalid params config")
        void shouldThrowConfigurationExceptionForInvalidParamsConfig() {
            TaxTriageConfig invalidConfig = new TaxTriageConfig.Builder()
                    .withPreset("TEST")
                    .withQualityThreshold(0)
                    .build();
            File outputFile = tempDir.resolve("params.json").toFile();

            assertThrows(ConfigGenerator.ConfigurationException.class, () ->
                    configGenerator.generateParams(invalidConfig, outputFile));
        }

        @Test
        @DisplayName("Should handle IOException when template loading fails")
        void shouldHandleIOExceptionWhenTemplateLoadingFails() {
            MockTemplateLoader failingLoader = new MockTemplateLoader();
            failingLoader.setShouldFail(true);
            ConfigGenerator failingGenerator = new ConfigGenerator(failingLoader);
            TaxTriageConfig config = TaxTriageConfig.createDefault();
            File outputFile = tempDir.resolve("nextflow.config").toFile();

            assertThrows(IOException.class, () ->
                    failingGenerator.generateConfig(config, outputFile));
        }

        @Test
        @DisplayName("Should handle read-only output directory")
        void shouldHandleReadOnlyOutputDirectory() throws Exception {
            File readOnlyDir = tempDir.resolve("readonly").toFile();
            readOnlyDir.mkdirs();
            readOnlyDir.setWritable(false);

            TaxTriageConfig config = TaxTriageConfig.createDefault();
            File outputFile = new File(readOnlyDir, "nextflow.config");

            try {
                assertThrows(IOException.class, () ->
                        configGenerator.generateConfig(config, outputFile));
            } finally {
                readOnlyDir.setWritable(true); // Cleanup
            }
        }
    }

    @Nested
    @DisplayName("Template Loading")
    class TemplateLoading {

        @Test
        @DisplayName("Should use template content for generation")
        void shouldUseTemplateContentForGeneration() throws Exception {
            mockTemplateLoader.setConfigTemplate("Custom config: preset = '${preset}'");
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("CUSTOM")
                    .withDefaults()
                    .build();
            File outputFile = tempDir.resolve("nextflow.config").toFile();

            configGenerator.generateConfig(config, outputFile);

            String content = Files.readString(outputFile.toPath());
            assertTrue(content.contains("Custom config: preset = 'CUSTOM'"));
        }

        @Test
        @DisplayName("Should use params template for JSON generation")
        void shouldUseParamsTemplateForJsonGeneration() throws Exception {
            mockTemplateLoader.setParamsTemplate("{ \"custom_preset\": ${preset} }");
            TaxTriageConfig config = new TaxTriageConfig.Builder()
                    .withPreset("CUSTOM")
                    .withDefaults()
                    .build();
            File outputFile = tempDir.resolve("params.json").toFile();

            configGenerator.generateParams(config, outputFile);

            String content = Files.readString(outputFile.toPath());
            assertTrue(content.contains("\"custom_preset\": \"CUSTOM\""));
        }
    }

    /**
     * Mock template loader for testing.
     */
    private static class MockTemplateLoader implements ConfigGenerator.TemplateLoader {
        private String configTemplate =
                "params {\n" +
                "    preset = '${preset}'\n" +
                "    quality_threshold = ${quality_threshold}\n" +
                "    min_read_length = ${min_read_length}\n" +
                "    enable_subsampling = ${enable_subsampling}\n" +
                "    subsample_size = ${subsample_size}\n" +
                "    max_cpus = ${max_cpus}\n" +
                "    max_memory = '${max_memory}'\n" +
                "    kraken_db = '${kraken_db}'\n" +
                "    bracken_db = '${bracken_db}'\n" +
                "    outdir = '${outdir}'\n" +
                "    top_taxa = ${top_taxa}\n" +
                "    confidence_threshold = ${confidence_threshold}\n" +
                "    enable_krona = ${enable_krona}\n" +
                "    enable_multiqc = ${enable_multiqc}\n" +
                "}\n\n" +
                "profiles {\n" +
                "    ${docker_profile} {\n" +
                "        docker.enabled = true\n" +
                "        docker.runOptions = \"-v ${docker_volumes}\"\n" +
                "    }\n" +
                "}\n";

        private String paramsTemplate =
                "{\n" +
                "    \"preset\": ${preset},\n" +
                "    \"quality_threshold\": ${quality_threshold},\n" +
                "    \"min_read_length\": ${min_read_length},\n" +
                "    \"enable_subsampling\": ${enable_subsampling},\n" +
                "    \"subsample_size\": ${subsample_size},\n" +
                "    \"max_cpus\": ${max_cpus},\n" +
                "    \"max_memory\": ${max_memory},\n" +
                "    \"kraken_db\": ${kraken_db},\n" +
                "    \"bracken_db\": ${bracken_db},\n" +
                "    \"outdir\": ${outdir},\n" +
                "    \"top_taxa\": ${top_taxa},\n" +
                "    \"confidence_threshold\": ${confidence_threshold},\n" +
                "    \"enable_krona\": ${enable_krona},\n" +
                "    \"enable_multiqc\": ${enable_multiqc}\n" +
                "}\n";

        private boolean shouldFail = false;

        public void setConfigTemplate(String template) {
            this.configTemplate = template;
        }

        public void setParamsTemplate(String template) {
            this.paramsTemplate = template;
        }

        public void setShouldFail(boolean shouldFail) {
            this.shouldFail = shouldFail;
        }

        @Override
        public String loadTemplate(String resourcePath) throws IOException {
            if (shouldFail) {
                throw new IOException("Template loading failed");
            }

            if (resourcePath.contains("nextflow.config")) {
                return configTemplate;
            } else if (resourcePath.contains("params.json")) {
                return paramsTemplate;
            } else {
                throw new IOException("Template not found: " + resourcePath);
            }
        }
    }
}