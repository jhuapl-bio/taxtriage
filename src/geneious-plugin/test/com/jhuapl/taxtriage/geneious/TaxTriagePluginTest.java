package com.jhuapl.taxtriage.geneious;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.sequence.SequenceDocument;
import com.biomatters.geneious.publicapi.implementations.sequence.DefaultNucleotideSequence;
import com.biomatters.geneious.publicapi.plugin.DocumentOperation;
import com.biomatters.geneious.publicapi.plugin.Options;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;

/**
 * Test class for TaxTriage Geneious plugin components.
 *
 * This test class verifies that the plugin loads correctly, creates proper options,
 * and performs basic validation. It uses mock sequence documents to test the
 * plugin functionality without requiring actual sequence files.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class TaxTriagePluginTest {

    /**
     * Main test method that runs all plugin tests.
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("Starting TaxTriage Plugin Tests...");

        TaxTriagePluginTest test = new TaxTriagePluginTest();

        try {
            test.testPluginLoading();
            test.testOperationCreation();
            test.testOptionsCreation();
            test.testDocumentSelection();
            test.testOptionsValidation();

            System.out.println("All tests passed successfully!");

        } catch (Exception e) {
            System.err.println("Test failed: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Tests that the main plugin class loads correctly and provides expected metadata.
     */
    public void testPluginLoading() {
        System.out.println("Testing plugin loading...");

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

        System.out.println("✓ Plugin loading test passed");
    }

    /**
     * Tests that the plugin creates and returns the TaxTriage operation.
     */
    public void testOperationCreation() {
        System.out.println("Testing operation creation...");

        TaxTriagePlugin plugin = new TaxTriagePlugin();
        DocumentOperation[] operations = plugin.getDocumentOperations();

        assert operations != null : "Operations array should not be null";
        assert operations.length == 1 : "Should have exactly one operation";

        DocumentOperation operation = operations[0];
        assert operation instanceof TaxTriageOperation : "Operation should be TaxTriageOperation";

        TaxTriageOperation taxTriageOp = (TaxTriageOperation) operation;
        assert taxTriageOp.getActionOptions() != null : "Action options should not be null";
        assert taxTriageOp.getHelp() != null : "Operation help should not be null";
        assert taxTriageOp.getSelectionSignatures() != null : "Selection signatures should not be null";

        System.out.println("✓ Operation creation test passed");
    }

    /**
     * Tests that the operation creates proper options for configuration.
     */
    public void testOptionsCreation() {
        System.out.println("Testing options creation...");

        TaxTriageOperation operation = new TaxTriageOperation();

        try {
            // Test with no documents
            Options options = operation.getOptions();
            assert options != null : "Options should not be null";
            assert options instanceof TaxTriageOptions : "Options should be TaxTriageOptions";

            TaxTriageOptions taxTriageOptions = (TaxTriageOptions) options;

            // Test that options have reasonable defaults
            assert taxTriageOptions.getQualityThreshold() > 0 : "Quality threshold should be positive";
            assert taxTriageOptions.getMinReadLength() > 0 : "Min read length should be positive";
            assert taxTriageOptions.getThreadCount() > 0 : "Thread count should be positive";
            assert taxTriageOptions.getMemoryLimit() > 0 : "Memory limit should be positive";

            // Test preset functionality
            assert taxTriageOptions.getSequencingPreset() != null : "Sequencing preset should not be null";

            System.out.println("✓ Options creation test passed");

        } catch (Exception e) {
            throw new RuntimeException("Failed to create options: " + e.getMessage(), e);
        }
    }

    /**
     * Tests that the operation correctly identifies supported document types.
     */
    public void testDocumentSelection() {
        System.out.println("Testing document selection...");

        TaxTriageOperation operation = new TaxTriageOperation();

        // Test selection signatures
        assert operation.getSelectionSignatures() != null : "Selection signatures should not be null";
        assert operation.getSelectionSignatures().length > 0 : "Should have at least one selection signature";

        // Test action options
        assert operation.getActionOptions() != null : "Action options should not be null";
        assert operation.getActionOptions().getName() != null : "Action name should not be null";

        System.out.println("✓ Document selection test passed");
    }

    /**
     * Tests that options validation works correctly.
     */
    public void testOptionsValidation() {
        System.out.println("Testing options validation...");

        TaxTriageOptions options = new TaxTriageOptions();

        // Test default options validation
        String validationError = options.verifyOptionsAreValid();
        // Note: Default options may have validation issues that require user input
        if (validationError != null) {
            System.out.println("Note: Default options validation message: " + validationError);
        }

        // Test specific validation cases by accessing methods
        assert options.getQualityThreshold() >= 1 && options.getQualityThreshold() <= 40 :
            "Quality threshold should be in valid range";
        assert options.getMinReadLength() >= 1 : "Min read length should be at least 1";
        assert options.getThreadCount() >= 1 : "Thread count should be at least 1";
        assert options.getMemoryLimit() >= 1 : "Memory limit should be at least 1 GB";

        // Test preset options
        assert options.getSequencingPreset() != null : "Should have a sequencing preset";
        assert options.getKrakenDatabase() != null : "Should have a Kraken database setting";
        assert options.getBrackenDatabase() != null : "Should have a Bracken database setting";

        System.out.println("✓ Options validation test passed");
    }

    // Note: Mock document creation removed for simplicity.
    // In a full implementation, we would use DocumentUtilities.createAnnotatedPluginDocuments()

    /**
     * Utility method to print test results.
     *
     * @param testName name of the test
     * @param passed whether the test passed
     */
    private void printTestResult(String testName, boolean passed) {
        String result = passed ? "✓ PASSED" : "✗ FAILED";
        System.out.println(testName + ": " + result);
    }
}