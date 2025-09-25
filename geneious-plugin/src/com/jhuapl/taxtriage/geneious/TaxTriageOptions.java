package com.jhuapl.taxtriage.geneious;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.plugin.Options;

import javax.swing.JFileChooser;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Configuration options for TaxTriage analysis operation.
 *
 * This class defines all user-configurable parameters for running TaxTriage
 * workflows within Geneious. It provides preset configurations for common
 * sequencing platforms and allows advanced customization of analysis parameters.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class TaxTriageOptions extends Options {

    // Option keys for serialization and retrieval
    private static final String INPUT_FILES_KEY = "inputFiles";
    private static final String INPUT_DIRECTORY_KEY = "inputDirectory";
    private static final String PRESET_KEY = "preset";
    private static final String DEDUPLICATE_READS_KEY = "deduplicateReads";
    private static final String QUALITY_THRESHOLD_KEY = "qualityThreshold";
    private static final String MIN_READ_LENGTH_KEY = "minReadLength";
    private static final String SUBSAMPLE_SIZE_KEY = "subsampleSize";
    private static final String ENABLE_SUBSAMPLING_KEY = "enableSubsampling";
    private static final String THREAD_COUNT_KEY = "threadCount";
    private static final String MEMORY_LIMIT_KEY = "memoryLimit";
    private static final String KRAKEN_DATABASE_KEY = "krakenDatabase";
    private static final String BRACKEN_DATABASE_KEY = "brackenDatabase";

    // Preset options
    public enum SequencingPreset {
        ONT("ONT (Oxford Nanopore)", "Optimized for long-read Oxford Nanopore sequencing"),
        ILLUMINA_PE("Illumina PE", "Optimized for Illumina paired-end sequencing"),
        ILLUMINA_SE("Illumina SE", "Optimized for Illumina single-end sequencing"),
        CUSTOM("Custom", "Custom configuration with manual parameter settings");

        private final String displayName;
        private final String description;

        SequencingPreset(String displayName, String description) {
            this.displayName = displayName;
            this.description = description;
        }

        public String getDisplayName() {
            return displayName;
        }

        public String getDescription() {
            return description;
        }

        @Override
        public String toString() {
            return displayName;
        }
    }

    /**
     * Constructs TaxTriage options for the given input documents.
     *
     * @param inputDocuments the sequence documents to be analyzed
     */
    public TaxTriageOptions(List<AnnotatedPluginDocument> inputDocuments) {
        super(TaxTriageOptions.class);
        initializeOptions();
    }

    /**
     * Default constructor for options deserialization.
     */
    public TaxTriageOptions() {
        super(TaxTriageOptions.class);
        initializeOptions();
    }

    /**
     * Initializes all option definitions and default values.
     */
    private void initializeOptions() {
        // === INPUT SELECTION SECTION ===
        addLabel("Input Selection (choose one of the following):");

        // Individual input files selection with browse button
        FileSelectionOption filesOption = addFileSelectionOption(INPUT_FILES_KEY, "Input Files:", "");
        filesOption.setDescription("Browse to select FASTQ or FASTA files for analysis");
        filesOption.setSelectionType(JFileChooser.FILES_ONLY);
        filesOption.setAllowMultipleSelection(true);

        // Input directory selection with browse button
        FileSelectionOption dirOption = addFileSelectionOption(INPUT_DIRECTORY_KEY, "Input Directory:", "");
        dirOption.setDescription("Browse to select directory containing FASTQ or FASTA files");
        dirOption.setSelectionType(JFileChooser.DIRECTORIES_ONLY);

        addLabel(" ");

        // === ANALYSIS CONFIGURATION ===
        addLabel("Analysis Configuration:");

        // Sequencing preset selection using combo box
        OptionValue[] presetValues = new OptionValue[SequencingPreset.values().length];
        for (int i = 0; i < SequencingPreset.values().length; i++) {
            SequencingPreset preset = SequencingPreset.values()[i];
            presetValues[i] = new OptionValue(preset.name(), preset.getDisplayName(), preset.getDescription());
        }
        addComboBoxOption(PRESET_KEY, "Sequencing Preset:", presetValues, presetValues[1]); // Default to Illumina PE

        // Basic parameters
        addIntegerOption(QUALITY_THRESHOLD_KEY, "Quality Threshold:", 20, 1, 40);
        getOption(QUALITY_THRESHOLD_KEY).setDescription("Minimum Phred quality score for sequence filtering");

        addIntegerOption(MIN_READ_LENGTH_KEY, "Min Read Length:", 50, 1, 10000);
        getOption(MIN_READ_LENGTH_KEY).setDescription("Minimum read length to retain after filtering");

        // Subsampling options
        addBooleanOption(ENABLE_SUBSAMPLING_KEY, "Enable Subsampling:", false);
        getOption(ENABLE_SUBSAMPLING_KEY).setDescription("Randomly subsample reads to reduce computational requirements");

        addIntegerOption(SUBSAMPLE_SIZE_KEY, "Subsample Size:", 100000, 1000, 10000000);
        getOption(SUBSAMPLE_SIZE_KEY).setDescription("Number of reads to retain when subsampling");

        addLabel(" ");

        // === DATABASE CONFIGURATION ===
        addLabel("Database Configuration:");

        // Database options using combo boxes
        OptionValue[] krakenDatabases = {
            new OptionValue("standard", "Standard", "Standard Kraken2 database with common organisms"),
            new OptionValue("standard-8", "Standard-8", "Standard database with 8-mer minimizers"),
            new OptionValue("standard-16", "Standard-16", "Standard database with 16-mer minimizers"),
            new OptionValue("minikraken", "MiniKraken", "Smaller database for faster analysis"),
            new OptionValue("viral", "Viral", "Viral-specific database"),
            new OptionValue("bacteria", "Bacteria", "Bacterial-specific database"),
            new OptionValue("custom", "Custom", "Custom database path (specify in advanced options)")
        };
        addComboBoxOption(KRAKEN_DATABASE_KEY, "Kraken Database:", krakenDatabases, krakenDatabases[0]);

        OptionValue[] brackenDatabases = {
            new OptionValue("standard", "Standard", "Standard Bracken database"),
            new OptionValue("standard-8", "Standard-8", "Standard database with 8-mer minimizers"),
            new OptionValue("standard-16", "Standard-16", "Standard database with 16-mer minimizers"),
            new OptionValue("minikraken", "MiniKraken", "Smaller database for faster analysis"),
            new OptionValue("viral", "Viral", "Viral-specific database"),
            new OptionValue("bacteria", "Bacteria", "Bacterial-specific database"),
            new OptionValue("custom", "Custom", "Custom database path (specify in advanced options)")
        };
        addComboBoxOption(BRACKEN_DATABASE_KEY, "Bracken Database:", brackenDatabases, brackenDatabases[0]);

        addLabel(" ");

        // === ADVANCED PARAMETERS ===
        addLabel("Advanced Parameters:");

        // Advanced parameters
        addIntegerOption(THREAD_COUNT_KEY, "Thread Count:", Runtime.getRuntime().availableProcessors(), 1, 32);
        getOption(THREAD_COUNT_KEY).setDescription("Number of CPU threads to use for analysis");

        addIntegerOption(MEMORY_LIMIT_KEY, "Memory Limit (GB):", 8, 1, 64);
        getOption(MEMORY_LIMIT_KEY).setDescription("Maximum memory allocation for TaxTriage processes");

        addLabel(" ");

        // === POST-PROCESSING OPTIONS ===
        addLabel("Post-Processing Options:");

        // Deduplication option
        addBooleanOption(DEDUPLICATE_READS_KEY, "Deduplicate Mapped Reads:", false);
        getOption(DEDUPLICATE_READS_KEY).setDescription("Remove PCR duplicates from BAM files using samtools markdup (requires samtools installed)");
    }

    /**
     * Gets the default output directory based on user preferences.
     *
     * @return default output directory
     */
    private File getDefaultOutputDirectory() {
        String userHome = System.getProperty("user.home");
        return new File(userHome, "TaxTriage_Results");
    }

    // Getter methods for accessing option values

    /**
     * Gets the input files from the file selection option.
     *
     * @return list of input files
     */
    public List<File> getInputFiles() {
        String value = (String) getValue(INPUT_FILES_KEY);
        List<File> files = new ArrayList<>();

        if (value != null && !value.trim().isEmpty()) {
            // File selection returns comma-separated paths for multiple files
            String[] filePaths = value.split(",");
            for (String path : filePaths) {
                String trimmedPath = path.trim();
                if (!trimmedPath.isEmpty()) {
                    files.add(new File(trimmedPath));
                }
            }
        }
        return files;
    }

    /**
     * Gets the input directory from the file selection option.
     *
     * @return input directory or null if not specified
     */
    public File getInputDirectory() {
        String value = (String) getValue(INPUT_DIRECTORY_KEY);
        if (value != null && !value.trim().isEmpty()) {
            return new File(value.trim());
        }
        return null;
    }

    /**
     * Gets whether to deduplicate mapped reads.
     *
     * @return true if deduplication is enabled
     */
    public boolean isDeduplicationEnabled() {
        Boolean value = (Boolean) getValue(DEDUPLICATE_READS_KEY);
        return value != null ? value : false;
    }

    /**
     * Gets the selected sequencing preset.
     *
     * @return sequencing preset enum value
     */
    public SequencingPreset getSequencingPreset() {
        OptionValue presetValue = (OptionValue) getValue(PRESET_KEY);
        if (presetValue != null) {
            return SequencingPreset.valueOf(presetValue.getName());
        }
        return SequencingPreset.ILLUMINA_PE; // default fallback
    }

    /**
     * Gets the quality threshold value.
     *
     * @return quality threshold
     */
    public int getQualityThreshold() {
        return (Integer) getValue(QUALITY_THRESHOLD_KEY);
    }

    /**
     * Gets the minimum read length.
     *
     * @return minimum read length
     */
    public int getMinReadLength() {
        return (Integer) getValue(MIN_READ_LENGTH_KEY);
    }

    /**
     * Gets whether subsampling is enabled.
     *
     * @return true if subsampling is enabled
     */
    public boolean isSubsamplingEnabled() {
        return (Boolean) getValue(ENABLE_SUBSAMPLING_KEY);
    }

    /**
     * Gets the subsample size.
     *
     * @return number of reads to subsample
     */
    public int getSubsampleSize() {
        return (Integer) getValue(SUBSAMPLE_SIZE_KEY);
    }

    /**
     * Gets the thread count for analysis.
     *
     * @return number of threads
     */
    public int getThreadCount() {
        return (Integer) getValue(THREAD_COUNT_KEY);
    }

    /**
     * Gets the memory limit in GB.
     *
     * @return memory limit
     */
    public int getMemoryLimit() {
        return (Integer) getValue(MEMORY_LIMIT_KEY);
    }

    /**
     * Gets the Kraken database configuration.
     *
     * @return Kraken database path or preset name
     */
    public String getKrakenDatabase() {
        OptionValue krakenValue = (OptionValue) getValue(KRAKEN_DATABASE_KEY);
        return krakenValue != null ? krakenValue.getName() : "standard";
    }

    /**
     * Gets the Bracken database configuration.
     *
     * @return Bracken database path or preset name
     */
    public String getBrackenDatabase() {
        OptionValue brackenValue = (OptionValue) getValue(BRACKEN_DATABASE_KEY);
        return brackenValue != null ? brackenValue.getName() : "standard";
    }

    /**
     * Validates all option values and returns any validation errors.
     *
     * @return validation error message or null if valid
     */
    @Override
    public String verifyOptionsAreValid() {
        // Validate input selection
        List<File> inputFiles = getInputFiles();
        File inputDirectory = getInputDirectory();

        if ((inputFiles == null || inputFiles.isEmpty()) && inputDirectory == null) {
            return "Please select either input files or an input directory";
        }

        // Validate input files if provided
        if (inputFiles != null && !inputFiles.isEmpty()) {
            for (File file : inputFiles) {
                if (!file.exists()) {
                    return "Input file does not exist: " + file.getAbsolutePath();
                }
                if (!file.isFile()) {
                    return "Input path is not a file: " + file.getAbsolutePath();
                }
            }
        }

        // Validate input directory if provided
        if (inputDirectory != null) {
            if (!inputDirectory.exists()) {
                return "Input directory does not exist: " + inputDirectory.getAbsolutePath();
            }
            if (!inputDirectory.isDirectory()) {
                return "Input path is not a directory: " + inputDirectory.getAbsolutePath();
            }
        }

        // Output directory no longer needed - results are imported to Geneious

        // Validate quality threshold
        int quality = getQualityThreshold();
        if (quality < 1 || quality > 40) {
            return "Quality threshold must be between 1 and 40";
        }

        // Validate read length
        int minLength = getMinReadLength();
        if (minLength < 1) {
            return "Minimum read length must be at least 1";
        }

        // Validate subsampling parameters
        if (isSubsamplingEnabled()) {
            int subsampleSize = getSubsampleSize();
            if (subsampleSize < 1000) {
                return "Subsample size must be at least 1000 reads";
            }
        }

        // Validate thread count
        int threads = getThreadCount();
        if (threads < 1) {
            return "Thread count must be at least 1";
        }

        // Validate memory limit
        int memory = getMemoryLimit();
        if (memory < 1) {
            return "Memory limit must be at least 1 GB";
        }

        return null; // All validations passed
    }
}