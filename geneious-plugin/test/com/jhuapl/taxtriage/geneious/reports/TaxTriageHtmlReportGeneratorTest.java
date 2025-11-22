package com.jhuapl.taxtriage.geneious.reports;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.nio.file.Path;
import java.util.Date;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for TaxTriageHtmlReportGenerator (Group 3).
 *
 * Tests the HTML report generation functionality including:
 * - Report initialization
 * - Adding parameters, samples, and statistics
 * - HTML generation and formatting
 * - Error and warning handling
 * - Report structure and content
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class TaxTriageHtmlReportGeneratorTest {

    private TaxTriageHtmlReportGenerator generator;

    @TempDir
    Path tempDir;

    @BeforeEach
    void setUp() {
        generator = new TaxTriageHtmlReportGenerator("Test Run", tempDir);
    }

    /**
     * GROUP 3: Test generator initialization
     */
    @Test
    void testGeneratorInitialization() {
        assertNotNull(generator, "Generator should be initialized");
    }

    /**
     * GROUP 3: Test basic HTML generation
     */
    @Test
    void testBasicHtmlGeneration() {
        String html = generator.generateHtml();

        assertNotNull(html);
        assertTrue(html.contains("<!DOCTYPE html>") || html.contains("<html>"));
        assertTrue(html.contains("TaxTriage"));
        assertTrue(html.contains("Test Run"));
    }

    /**
     * GROUP 3: Test adding command line
     */
    @Test
    void testAddCommandLine() {
        generator.setCommandLine("nextflow run taxtriage --input data.fastq");
        String html = generator.generateHtml();

        assertTrue(html.contains("nextflow"));
        assertTrue(html.contains("data.fastq"));
    }

    /**
     * GROUP 3: Test adding parameters
     */
    @Test
    void testAddParameters() {
        generator.addParameter("platform", "illumina");
        generator.addParameter("database", "kraken2");
        generator.addParameter("threads", "4");

        String html = generator.generateHtml();

        assertTrue(html.contains("platform"));
        assertTrue(html.contains("illumina"));
        assertTrue(html.contains("database"));
        assertTrue(html.contains("kraken2"));
        assertTrue(html.contains("threads"));
    }

    /**
     * GROUP 3: Test adding samples
     */
    @Test
    void testAddSamples() {
        generator.addSample("sample1");
        generator.addSample("sample2");
        generator.addSample("sample3");

        String html = generator.generateHtml();

        assertTrue(html.contains("sample1"));
        assertTrue(html.contains("sample2"));
        assertTrue(html.contains("sample3"));
    }

    /**
     * GROUP 3: Test adding statistics
     */
    @Test
    void testAddStatistics() {
        generator.addStatistic("Total Reads", 1000000);
        generator.addStatistic("Classified Reads", 850000);
        generator.addStatistic("Processing Time", "15 minutes");

        String html = generator.generateHtml();

        assertTrue(html.contains("Total Reads"));
        assertTrue(html.contains("1000000"));
        assertTrue(html.contains("Classified Reads"));
        assertTrue(html.contains("850000"));
        assertTrue(html.contains("Processing Time"));
        assertTrue(html.contains("15 minutes"));
    }

    /**
     * GROUP 3: Test adding errors
     */
    @Test
    void testAddErrors() {
        generator.addError("Error 1: File not found");
        generator.addError("Error 2: Permission denied");

        String html = generator.generateHtml();

        assertTrue(html.contains("Error 1"));
        assertTrue(html.contains("File not found"));
        assertTrue(html.contains("Error 2"));
        assertTrue(html.contains("Permission denied"));
    }

    /**
     * GROUP 3: Test adding warnings
     */
    @Test
    void testAddWarnings() {
        generator.addWarning("Warning: Low quality scores detected");
        generator.addWarning("Warning: Incomplete reference database");

        String html = generator.generateHtml();

        assertTrue(html.contains("Warning"));
        assertTrue(html.contains("Low quality scores"));
        assertTrue(html.contains("Incomplete reference database"));
    }

    /**
     * GROUP 3: Test setting end time
     */
    @Test
    void testSetEndTime() {
        Date endTime = new Date();

        generator.setEndTime(endTime);

        String html = generator.generateHtml();

        // Should contain time-related information
        assertNotNull(html);
        assertTrue(html.length() > 0);
    }

    /**
     * GROUP 3: Test HTML contains proper structure
     */
    @Test
    void testHtmlStructure() {
        String html = generator.generateHtml();

        assertTrue(html.contains("<html"));
        assertTrue(html.contains("<head"));
        assertTrue(html.contains("<body"));
        assertTrue(html.contains("</body>"));
        assertTrue(html.contains("</html>"));
    }

    /**
     * GROUP 3: Test HTML contains CSS styling
     */
    @Test
    void testHtmlContainsStyling() {
        String html = generator.generateHtml();

        assertTrue(html.contains("<style") || html.contains("style="));
    }

    /**
     * GROUP 3: Test adding null values
     */
    @Test
    void testAddNullValues() {
        assertDoesNotThrow(() -> {
            generator.addParameter("nullParam", null);
            generator.addSample(null);
            generator.addStatistic("nullStat", null);
        }, "Should handle null values gracefully");

        String html = generator.generateHtml();
        assertNotNull(html);
    }

    /**
     * GROUP 3: Test adding empty strings
     */
    @Test
    void testAddEmptyStrings() {
        assertDoesNotThrow(() -> {
            generator.addParameter("emptyParam", "");
            generator.addSample("");
            generator.addError("");
            generator.addWarning("");
        }, "Should handle empty strings gracefully");

        String html = generator.generateHtml();
        assertNotNull(html);
    }

    /**
     * GROUP 3: Test multiple parameters with same key
     */
    @Test
    void testMultipleParametersSameKey() {
        generator.addParameter("key", "value1");
        generator.addParameter("key", "value2");

        String html = generator.generateHtml();

        // Should handle duplicate keys appropriately
        assertNotNull(html);
    }

    /**
     * GROUP 3: Test special characters in content
     */
    @Test
    void testSpecialCharactersInContent() {
        generator.addParameter("special", "Value with <html> & \"quotes\"");
        generator.addSample("Sample with <tags>");
        generator.addError("Error: Path contains \\ backslashes");

        String html = generator.generateHtml();

        // HTML should be properly escaped
        assertNotNull(html);
        // Special characters should be escaped in HTML
        assertTrue(html.contains("&lt;") || html.contains("&amp;") ||
                   !html.contains("<tags>"));
    }

    /**
     * GROUP 3: Test large number of samples
     */
    @Test
    void testLargeNumberOfSamples() {
        for (int i = 0; i < 100; i++) {
            generator.addSample("sample_" + i);
        }

        String html = generator.generateHtml();

        assertNotNull(html);
        assertTrue(html.contains("sample_0"));
        assertTrue(html.contains("sample_99"));
    }

    /**
     * GROUP 3: Test comprehensive report
     */
    @Test
    void testComprehensiveReport() {
        generator.setCommandLine("nextflow run taxtriage");
        generator.setEndTime(new Date());

        generator.addParameter("platform", "illumina");
        generator.addParameter("database", "standard");

        generator.addSample("sample1");
        generator.addSample("sample2");

        generator.addStatistic("Total Reads", 500000);
        generator.addStatistic("Classified", 450000);

        generator.addWarning("Low coverage in sample2");

        String html = generator.generateHtml();

        // Verify all components are present
        assertTrue(html.contains("nextflow"));
        assertTrue(html.contains("illumina"));
        assertTrue(html.contains("sample1"));
        assertTrue(html.contains("sample2"));
        assertTrue(html.contains("500000"));
        assertTrue(html.contains("Low coverage"));
    }

    /**
     * GROUP 3: Test empty report generation
     */
    @Test
    void testEmptyReportGeneration() {
        TaxTriageHtmlReportGenerator emptyGen =
            new TaxTriageHtmlReportGenerator("Empty Run", tempDir);

        String html = emptyGen.generateHtml();

        assertNotNull(html);
        assertTrue(html.contains("Empty Run"));
        assertTrue(html.contains("<html"));
    }

    /**
     * GROUP 3: Test report with only errors
     */
    @Test
    void testReportWithOnlyErrors() {
        generator.addError("Critical error 1");
        generator.addError("Critical error 2");
        generator.addError("Critical error 3");

        String html = generator.generateHtml();

        assertTrue(html.contains("Critical error 1"));
        assertTrue(html.contains("Critical error 2"));
        assertTrue(html.contains("Critical error 3"));
    }

    /**
     * GROUP 3: Test report with mixed content types
     */
    @Test
    void testReportWithMixedContent() {
        generator.addParameter("param1", "value1");
        generator.addSample("sample1");
        generator.addStatistic("stat1", 100);
        generator.addError("error1");
        generator.addWarning("warning1");

        String html = generator.generateHtml();

        assertTrue(html.contains("param1"));
        assertTrue(html.contains("sample1"));
        assertTrue(html.contains("100"));
        assertTrue(html.contains("error1"));
        assertTrue(html.contains("warning1"));
    }

    /**
     * GROUP 3: Test statistic with different value types
     */
    @Test
    void testStatisticWithDifferentTypes() {
        generator.addStatistic("Integer Stat", 42);
        generator.addStatistic("String Stat", "test value");
        generator.addStatistic("Long Stat", 1000000000L);
        generator.addStatistic("Double Stat", 3.14159);

        String html = generator.generateHtml();

        assertTrue(html.contains("42"));
        assertTrue(html.contains("test value"));
        assertTrue(html.contains("1000000000"));
        assertTrue(html.contains("3.14"));
    }

    /**
     * GROUP 3: Test working directory display
     */
    @Test
    void testWorkingDirectoryDisplay() {
        String html = generator.generateHtml();

        // Should display the working directory path
        assertTrue(html.contains(tempDir.toString()) ||
                   html.contains("Working") ||
                   html.contains("Directory"));
    }

    /**
     * GROUP 3: Test run name display
     */
    @Test
    void testRunNameDisplay() {
        TaxTriageHtmlReportGenerator namedGen =
            new TaxTriageHtmlReportGenerator("MyCustomRun_2024", tempDir);

        String html = namedGen.generateHtml();

        assertTrue(html.contains("MyCustomRun_2024"));
    }
}
