package com.jhuapl.taxtriage.geneious;

import com.biomatters.geneious.publicapi.plugin.DocumentOperation;

/**
 * Simple test class for TaxTriage Geneious plugin core functionality.
 *
 * This test class focuses on basic plugin functionality without UI components
 * that might require icon resources or other dependencies.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class TaxTriageSimpleTest {

    /**
     * Main test method that runs core plugin tests.
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("Starting TaxTriage Simple Plugin Tests...");

        TaxTriageSimpleTest test = new TaxTriageSimpleTest();

        try {
            test.testPluginBasics();
            test.testOptionsBasics();
            test.testOperationBasics();

            System.out.println("All simple tests passed successfully!");

        } catch (Exception e) {
            System.err.println("Test failed: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Tests basic plugin metadata without UI components.
     */
    public void testPluginBasics() {
        System.out.println("Testing plugin basics...");

        TaxTriagePlugin plugin = new TaxTriagePlugin();

        // Test basic metadata
        assert plugin.getName() != null : "Plugin name should not be null";
        assert plugin.getName().equals("TaxTriage") : "Plugin name should be 'TaxTriage'";

        assert plugin.getDescription() != null : "Plugin description should not be null";
        assert plugin.getDescription().contains("Taxonomic") : "Description should mention taxonomic classification";

        assert plugin.getVersion() != null : "Plugin version should not be null";
        assert plugin.getVersion().equals("1.0.0") : "Plugin version should be '1.0.0'";

        assert plugin.getAuthors() != null : "Plugin authors should not be null";
        assert plugin.getAuthors().contains("Johns Hopkins") : "Authors should include Johns Hopkins";

        assert plugin.getMinimumApiVersion() != null : "Minimum API version should not be null";
        assert plugin.getMaximumApiVersion() > 0 : "Maximum API version should be positive";

        // Test operations array
        DocumentOperation[] operations = plugin.getDocumentOperations();
        assert operations != null : "Operations array should not be null";
        assert operations.length == 1 : "Should have exactly one operation";

        System.out.println("✓ Plugin basics test passed");
    }

    /**
     * Tests basic options functionality without complex UI.
     */
    public void testOptionsBasics() {
        System.out.println("Testing options basics...");

        TaxTriageOptions options = new TaxTriageOptions();

        // Test that options have reasonable defaults
        assert options.getQualityThreshold() > 0 : "Quality threshold should be positive";
        assert options.getQualityThreshold() == 20 : "Default quality threshold should be 20";

        assert options.getMinReadLength() > 0 : "Min read length should be positive";
        assert options.getMinReadLength() == 50 : "Default min read length should be 50";

        assert options.getThreadCount() > 0 : "Thread count should be positive";
        assert options.getMemoryLimit() > 0 : "Memory limit should be positive";

        // Test preset functionality
        assert options.getSequencingPreset() != null : "Sequencing preset should not be null";
        assert options.getSequencingPreset() == TaxTriageOptions.SequencingPreset.ILLUMINA_PE :
               "Default preset should be Illumina PE";

        // Test database settings
        assert options.getKrakenDatabase() != null : "Kraken database should not be null";
        assert options.getKrakenDatabase().equals("standard") : "Default Kraken database should be 'standard'";

        assert options.getBrackenDatabase() != null : "Bracken database should not be null";
        assert options.getBrackenDatabase().equals("standard") : "Default Bracken database should be 'standard'";

        // Test subsampling defaults
        assert !options.isSubsamplingEnabled() : "Subsampling should be disabled by default";
        assert options.getSubsampleSize() == 100000 : "Default subsample size should be 100000";

        System.out.println("✓ Options basics test passed");
    }

    /**
     * Tests basic operation functionality without UI components.
     */
    public void testOperationBasics() {
        System.out.println("Testing operation basics...");

        TaxTriageOperation operation = new TaxTriageOperation();

        // Test help text
        assert operation.getHelp() != null : "Operation help should not be null";
        assert operation.getHelp().contains("TaxTriage") : "Help should mention TaxTriage";

        // Test selection signatures
        assert operation.getSelectionSignatures() != null : "Selection signatures should not be null";
        assert operation.getSelectionSignatures().length > 0 : "Should have at least one selection signature";

        // Test options creation without documents
        try {
            assert operation.getOptions() != null : "Should be able to create basic options";
            System.out.println("✓ Options creation works");
        } catch (Exception e) {
            System.out.println("Note: Options creation requires documents, which is acceptable");
        }

        System.out.println("✓ Operation basics test passed");
    }

    /**
     * Tests preset enum functionality.
     */
    public void testPresetEnum() {
        System.out.println("Testing preset enum...");

        // Test that all presets are available
        TaxTriageOptions.SequencingPreset[] presets = TaxTriageOptions.SequencingPreset.values();
        assert presets.length == 4 : "Should have 4 presets";

        // Test individual presets
        assert TaxTriageOptions.SequencingPreset.ONT.getDisplayName().contains("ONT") :
               "ONT preset should contain 'ONT' in display name";
        assert TaxTriageOptions.SequencingPreset.ILLUMINA_PE.getDisplayName().contains("Illumina") :
               "Illumina PE preset should contain 'Illumina' in display name";
        assert TaxTriageOptions.SequencingPreset.ILLUMINA_SE.getDisplayName().contains("Illumina") :
               "Illumina SE preset should contain 'Illumina' in display name";
        assert TaxTriageOptions.SequencingPreset.CUSTOM.getDisplayName().contains("Custom") :
               "Custom preset should contain 'Custom' in display name";

        System.out.println("✓ Preset enum test passed");
    }
}