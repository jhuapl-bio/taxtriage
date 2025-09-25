package com.jhuapl.taxtriage.geneious.config;

import com.jhuapl.taxtriage.geneious.TaxTriageOptions;
import com.jhuapl.taxtriage.geneious.database.DatabaseManager;
import com.jhuapl.taxtriage.geneious.database.DatabaseManager.DatabaseType;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

/**
 * Data class representing TaxTriage Nextflow configuration parameters.
 *
 * This class encapsulates all configuration parameters needed to generate
 * a Nextflow configuration file for TaxTriage workflows. It provides
 * validation methods and conversion utilities for different output formats.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class TaxTriageConfig {

    // Core workflow parameters
    private final String preset;
    private final int qualityThreshold;
    private final int minReadLength;
    private final boolean enableSubsampling;
    private final int subsampleSize;
    private final int threadCount;
    private final int memoryLimitGb;
    private final String krakenDatabase;
    private final String brackenDatabase;
    private final File outputDirectory;

    // Advanced parameters
    private final int topTaxa;
    private final double confidenceThreshold;
    private final boolean enableKrona;
    private final boolean enableMultiQC;
    private final String dockerProfile;

    /**
     * Private constructor for builder pattern.
     */
    private TaxTriageConfig(Builder builder) {
        this.preset = builder.preset;
        this.qualityThreshold = builder.qualityThreshold;
        this.minReadLength = builder.minReadLength;
        this.enableSubsampling = builder.enableSubsampling;
        this.subsampleSize = builder.subsampleSize;
        this.threadCount = builder.threadCount;
        this.memoryLimitGb = builder.memoryLimitGb;
        this.krakenDatabase = builder.krakenDatabase;
        this.brackenDatabase = builder.brackenDatabase;
        this.outputDirectory = builder.outputDirectory;
        this.topTaxa = builder.topTaxa;
        this.confidenceThreshold = builder.confidenceThreshold;
        this.enableKrona = builder.enableKrona;
        this.enableMultiQC = builder.enableMultiQC;
        this.dockerProfile = builder.dockerProfile;
    }

    /**
     * Creates a TaxTriageConfig from TaxTriageOptions.
     *
     * @param options the TaxTriage options from the Geneious plugin
     * @return configured TaxTriageConfig instance
     */
    public static TaxTriageConfig fromOptions(TaxTriageOptions options) {
        return new Builder()
                .withPreset(options.getSequencingPreset().name())
                .withQualityThreshold(options.getQualityThreshold())
                .withMinReadLength(options.getMinReadLength())
                .withSubsampling(options.isSubsamplingEnabled(), options.getSubsampleSize())
                .withComputeResources(options.getThreadCount(), options.getMemoryLimit())
                .withDatabases(options.getKrakenDatabase(), options.getBrackenDatabase())
                .withOutputDirectory(null) // Will use temp directory
                .withDefaults()
                .build();
    }

    /**
     * Creates a default configuration for testing.
     *
     * @return default TaxTriageConfig instance
     */
    public static TaxTriageConfig createDefault() {
        return new Builder()
                .withPreset("ILLUMINA_PE")
                .withQualityThreshold(20)
                .withMinReadLength(50)
                .withSubsampling(false, 100000)
                .withComputeResources(4, 8)
                .withDatabases("standard", "standard")
                .withOutputDirectory(new File(System.getProperty("user.home"), "TaxTriage_Results"))
                .withDefaults()
                .build();
    }

    /**
     * Validates all configuration parameters.
     *
     * @return validation error message or null if valid
     */
    public String validate() {
        if (preset == null || preset.trim().isEmpty()) {
            return "Preset must be specified";
        }

        if (qualityThreshold < 1 || qualityThreshold > 40) {
            return "Quality threshold must be between 1 and 40";
        }

        if (minReadLength < 1) {
            return "Minimum read length must be at least 1";
        }

        if (enableSubsampling && subsampleSize < 1000) {
            return "Subsample size must be at least 1000 reads";
        }

        if (threadCount < 1 || threadCount > 128) {
            return "Thread count must be between 1 and 128";
        }

        if (memoryLimitGb < 1 || memoryLimitGb > 1024) {
            return "Memory limit must be between 1 and 1024 GB";
        }

        if (krakenDatabase == null || krakenDatabase.trim().isEmpty()) {
            return "Kraken database must be specified";
        }

        if (brackenDatabase == null || brackenDatabase.trim().isEmpty()) {
            return "Bracken database must be specified";
        }

        if (outputDirectory == null) {
            return "Output directory must be specified";
        }

        if (topTaxa < 1 || topTaxa > 1000) {
            return "Top taxa count must be between 1 and 1000";
        }

        if (confidenceThreshold < 0.0 || confidenceThreshold > 1.0) {
            return "Confidence threshold must be between 0.0 and 1.0";
        }

        return null; // All validations passed
    }

    /**
     * Converts configuration to a map suitable for Nextflow parameters.
     *
     * @return map of parameter names to values
     */
    public Map<String, Object> toParameterMap() {
        Map<String, Object> params = new HashMap<>();

        // Core parameters
        params.put("preset", preset);
        params.put("quality_threshold", qualityThreshold);
        params.put("min_read_length", minReadLength);
        params.put("enable_subsampling", enableSubsampling);
        params.put("subsample_size", subsampleSize);
        params.put("max_cpus", threadCount);
        params.put("max_memory", memoryLimitGb + ".GB");

        // Use DatabaseManager to resolve database paths
        DatabaseManager dbManager = DatabaseManager.getInstance();
        String resolvedKrakenDb = resolveDatabasePath(krakenDatabase, true);
        String resolvedBrackenDb = resolveDatabasePath(brackenDatabase, false);

        params.put("kraken_db", resolvedKrakenDb);
        params.put("bracken_db", resolvedBrackenDb);

        // Add download_db flag based on whether databases need downloading
        boolean needsDownload = "DOWNLOAD_REQUIRED".equals(resolvedKrakenDb) ||
                               "DOWNLOAD_REQUIRED".equals(resolvedBrackenDb);
        params.put("download_db", needsDownload);

        params.put("outdir", outputDirectory.getAbsolutePath());

        // Advanced parameters
        params.put("top_taxa", topTaxa);
        params.put("confidence_threshold", confidenceThreshold);
        params.put("enable_krona", enableKrona);
        params.put("enable_multiqc", enableMultiQC);

        return params;
    }

    /**
     * Resolves a database name to a path using DatabaseManager.
     *
     * @param databaseName The database name (e.g., "viral", "standard")
     * @param allowDownload Whether to allow downloading if not cached
     * @return The resolved database path or "DOWNLOAD_REQUIRED" if download is needed
     */
    private String resolveDatabasePath(String databaseName, boolean allowDownload) {
        DatabaseManager dbManager = DatabaseManager.getInstance();

        // Map database name to DatabaseType
        DatabaseType dbType = mapDatabaseNameToType(databaseName);
        if (dbType == null) {
            // If not a recognized type, return as-is (might be a custom path)
            return databaseName;
        }

        // Get the best available path for this database
        String dbPath = dbManager.getDatabasePath(dbType, allowDownload);

        // If no path available and download not allowed, use the name as-is
        // (will trigger download via workflow)
        if (dbPath == null) {
            return databaseName;
        }

        return dbPath;
    }

    /**
     * Maps a database name string to a DatabaseType enum.
     */
    private DatabaseType mapDatabaseNameToType(String databaseName) {
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
     * Gets the Docker profile configuration based on preset.
     *
     * @return Docker profile name
     */
    public String getDockerProfile() {
        return dockerProfile != null ? dockerProfile : "docker";
    }

    /**
     * Checks if this is a long-read preset (ONT).
     *
     * @return true if ONT preset
     */
    public boolean isLongReadPreset() {
        return "ONT".equals(preset);
    }

    /**
     * Checks if this is a paired-end preset.
     *
     * @return true if paired-end preset
     */
    public boolean isPairedEndPreset() {
        return "ILLUMINA_PE".equals(preset);
    }

    // Getters
    public String getPreset() { return preset; }
    public int getQualityThreshold() { return qualityThreshold; }
    public int getMinReadLength() { return minReadLength; }
    public boolean isEnableSubsampling() { return enableSubsampling; }
    public int getSubsampleSize() { return subsampleSize; }
    public int getThreadCount() { return threadCount; }
    public int getMemoryLimitGb() { return memoryLimitGb; }
    public String getKrakenDatabase() { return krakenDatabase; }
    public String getBrackenDatabase() { return brackenDatabase; }
    public File getOutputDirectory() { return outputDirectory; }
    public int getTopTaxa() { return topTaxa; }
    public double getConfidenceThreshold() { return confidenceThreshold; }
    public boolean isEnableKrona() { return enableKrona; }
    public boolean isEnableMultiQC() { return enableMultiQC; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        TaxTriageConfig that = (TaxTriageConfig) o;
        return qualityThreshold == that.qualityThreshold &&
               minReadLength == that.minReadLength &&
               enableSubsampling == that.enableSubsampling &&
               subsampleSize == that.subsampleSize &&
               threadCount == that.threadCount &&
               memoryLimitGb == that.memoryLimitGb &&
               topTaxa == that.topTaxa &&
               Double.compare(that.confidenceThreshold, confidenceThreshold) == 0 &&
               enableKrona == that.enableKrona &&
               enableMultiQC == that.enableMultiQC &&
               Objects.equals(preset, that.preset) &&
               Objects.equals(krakenDatabase, that.krakenDatabase) &&
               Objects.equals(brackenDatabase, that.brackenDatabase) &&
               Objects.equals(outputDirectory, that.outputDirectory) &&
               Objects.equals(dockerProfile, that.dockerProfile);
    }

    @Override
    public int hashCode() {
        return Objects.hash(preset, qualityThreshold, minReadLength, enableSubsampling,
                          subsampleSize, threadCount, memoryLimitGb, krakenDatabase,
                          brackenDatabase, outputDirectory, topTaxa, confidenceThreshold,
                          enableKrona, enableMultiQC, dockerProfile);
    }

    @Override
    public String toString() {
        return "TaxTriageConfig{" +
               "preset='" + preset + '\'' +
               ", qualityThreshold=" + qualityThreshold +
               ", minReadLength=" + minReadLength +
               ", enableSubsampling=" + enableSubsampling +
               ", subsampleSize=" + subsampleSize +
               ", threadCount=" + threadCount +
               ", memoryLimitGb=" + memoryLimitGb +
               ", krakenDatabase='" + krakenDatabase + '\'' +
               ", brackenDatabase='" + brackenDatabase + '\'' +
               ", outputDirectory=" + outputDirectory +
               ", topTaxa=" + topTaxa +
               ", confidenceThreshold=" + confidenceThreshold +
               ", enableKrona=" + enableKrona +
               ", enableMultiQC=" + enableMultiQC +
               ", dockerProfile='" + dockerProfile + '\'' +
               '}';
    }

    /**
     * Builder class for TaxTriageConfig.
     */
    public static class Builder {
        private String preset;
        private int qualityThreshold = 20;
        private int minReadLength = 50;
        private boolean enableSubsampling = false;
        private int subsampleSize = 100000;
        private int threadCount = 4;
        private int memoryLimitGb = 8;
        private String krakenDatabase = "standard";
        private String brackenDatabase = "standard";
        private File outputDirectory;
        private int topTaxa = 20;
        private double confidenceThreshold = 0.1;
        private boolean enableKrona = true;
        private boolean enableMultiQC = true;
        private String dockerProfile = "docker";

        public Builder withPreset(String preset) {
            this.preset = preset;
            return this;
        }

        public Builder withQualityThreshold(int qualityThreshold) {
            this.qualityThreshold = qualityThreshold;
            return this;
        }

        public Builder withMinReadLength(int minReadLength) {
            this.minReadLength = minReadLength;
            return this;
        }

        public Builder withSubsampling(boolean enable, int size) {
            this.enableSubsampling = enable;
            this.subsampleSize = size;
            return this;
        }

        public Builder withComputeResources(int threads, int memoryGb) {
            this.threadCount = threads;
            this.memoryLimitGb = memoryGb;
            return this;
        }

        public Builder withDatabases(String kraken, String bracken) {
            this.krakenDatabase = kraken;
            this.brackenDatabase = bracken;
            return this;
        }

        public Builder withOutputDirectory(File outputDirectory) {
            this.outputDirectory = outputDirectory;
            return this;
        }

        public Builder withTopTaxa(int topTaxa) {
            this.topTaxa = topTaxa;
            return this;
        }

        public Builder withConfidenceThreshold(double confidenceThreshold) {
            this.confidenceThreshold = confidenceThreshold;
            return this;
        }

        public Builder withReporting(boolean enableKrona, boolean enableMultiQC) {
            this.enableKrona = enableKrona;
            this.enableMultiQC = enableMultiQC;
            return this;
        }

        public Builder withDockerProfile(String dockerProfile) {
            this.dockerProfile = dockerProfile;
            return this;
        }

        public Builder withDefaults() {
            if (outputDirectory == null) {
                outputDirectory = new File(System.getProperty("user.home"), "TaxTriage_Results");
            }
            return this;
        }

        public TaxTriageConfig build() {
            return new TaxTriageConfig(this);
        }
    }
}