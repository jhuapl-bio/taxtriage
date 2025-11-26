package com.jhuapl.taxtriage.geneious.documents;

import com.biomatters.geneious.publicapi.documents.DocumentField;
import org.jdom.Element;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Date;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for TaxTriageResultDocument (Group 3).
 *
 * Tests the result document functionality including:
 * - Document creation and initialization
 * - HTML conversion for different formats
 * - XML serialization and deserialization
 * - Field values and display
 * - Content handling
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class TaxTriageResultDocumentTest {

    private TaxTriageResultDocument document;

    @BeforeEach
    void setUp() {
        document = new TaxTriageResultDocument(
            "Test Document",
            "Sample content",
            "Test Type",
            "/path/to/file.txt",
            "TXT"
        );
    }

    /**
     * GROUP 3: Test document initialization with parameters
     */
    @Test
    void testDocumentInitialization() {
        assertNotNull(document, "Document should be initialized");
        assertEquals("Test Document", document.getName());
        assertEquals("Sample content", document.getContent());
        assertEquals("Test Type", document.getResultType());
        assertEquals("TXT", document.getFileFormat());
        assertNotNull(document.getCreationDate());
    }

    /**
     * GROUP 3: Test empty constructor
     */
    @Test
    void testEmptyConstructor() {
        TaxTriageResultDocument emptyDoc = new TaxTriageResultDocument();
        assertNotNull(emptyDoc, "Document should be created with empty constructor");
    }

    /**
     * GROUP 3: Test getName
     */
    @Test
    void testGetName() {
        assertEquals("Test Document", document.getName());
    }

    /**
     * GROUP 3: Test getDescription
     */
    @Test
    void testGetDescription() {
        String description = document.getDescription();
        assertNotNull(description);
        assertTrue(description.contains("TaxTriage"));
        assertTrue(description.contains("Test Type"));
    }

    /**
     * GROUP 3: Test getCreationDate
     */
    @Test
    void testGetCreationDate() {
        Date creationDate = document.getCreationDate();
        assertNotNull(creationDate);
        assertTrue(creationDate.getTime() <= System.currentTimeMillis());
    }

    /**
     * GROUP 3: Test getURN returns null (not persisted)
     */
    @Test
    void testGetURN() {
        assertNull(document.getURN(), "URN should be null for unpersisted document");
    }

    /**
     * GROUP 3: Test toHTML for text format
     */
    @Test
    void testToHtmlForTextFormat() {
        String html = document.toHTML();

        assertNotNull(html);
        assertTrue(html.contains("<pre>"));
        assertTrue(html.contains("Sample content"));
        assertTrue(html.contains("</pre>"));
    }

    /**
     * GROUP 3: Test toHTML for HTML format
     */
    @Test
    void testToHtmlForHtmlFormat() {
        TaxTriageResultDocument htmlDoc = new TaxTriageResultDocument(
            "HTML Doc",
            "<h1>Title</h1><p>Content</p>",
            "HTML Report",
            "/path/to/report.html",
            "HTML"
        );

        String html = htmlDoc.toHTML();
        assertNotNull(html);
        assertTrue(html.contains("Title"));
        assertTrue(html.contains("Content"));
        // Script tags should be removed for security
        assertFalse(html.contains("<script"));
    }

    /**
     * GROUP 3: Test toHTML removes script tags
     */
    @Test
    void testToHtmlRemovesScriptTags() {
        TaxTriageResultDocument htmlDoc = new TaxTriageResultDocument(
            "HTML Doc",
            "<h1>Title</h1><script>alert('test');</script><p>Content</p>",
            "HTML Report",
            "/path/to/report.html",
            "HTML"
        );

        String html = htmlDoc.toHTML();
        assertNotNull(html);
        assertFalse(html.toLowerCase().contains("<script"),
            "Script tags should be removed");
    }

    /**
     * GROUP 3: Test toHTML for TSV format
     */
    @Test
    void testToHtmlForTsvFormat() {
        String tsvContent = "Header1\tHeader2\tHeader3\n" +
                           "Value1\tValue2\tValue3\n" +
                           "Value4\tValue5\tValue6";

        TaxTriageResultDocument tsvDoc = new TaxTriageResultDocument(
            "TSV Doc",
            tsvContent,
            "Data Table",
            "/path/to/data.tsv",
            "TSV"
        );

        String html = tsvDoc.toHTML();
        assertNotNull(html);
        assertTrue(html.contains("<table>"));
        assertTrue(html.contains("<th>"));
        assertTrue(html.contains("<td>"));
        assertTrue(html.contains("Header1"));
        assertTrue(html.contains("Value1"));
    }

    /**
     * GROUP 3: Test toHTML for CSV format
     */
    @Test
    void testToHtmlForCsvFormat() {
        String csvContent = "Name,Age,City\n" +
                           "Alice,30,NYC\n" +
                           "Bob,25,LA";

        TaxTriageResultDocument csvDoc = new TaxTriageResultDocument(
            "CSV Doc",
            csvContent,
            "Data Table",
            "/path/to/data.csv",
            "CSV"
        );

        String html = csvDoc.toHTML();
        assertNotNull(html);
        assertTrue(html.contains("<table>"));
        assertTrue(html.contains("Name"));
        assertTrue(html.contains("Alice"));
    }

    /**
     * GROUP 3: Test HTML escaping
     */
    @Test
    void testHtmlEscaping() {
        TaxTriageResultDocument doc = new TaxTriageResultDocument(
            "Test",
            "Content with <html> & \"special\" characters",
            "Test",
            "/path",
            "TXT"
        );

        String html = doc.toHTML();
        assertTrue(html.contains("&lt;html&gt;"));
        assertTrue(html.contains("&amp;"));
        assertTrue(html.contains("&quot;"));
    }

    /**
     * GROUP 3: Test toXML
     */
    @Test
    void testToXml() {
        Element xml = document.toXML();

        assertNotNull(xml);
        assertEquals("TaxTriageResultDocument", xml.getName());
        assertEquals("Test Document", xml.getChildText("name"));
        assertEquals("Test Type", xml.getChildText("resultType"));
        assertEquals("TXT", xml.getChildText("fileFormat"));
        assertEquals("Sample content", xml.getChildText("content"));
        assertNotNull(xml.getChildText("creationDate"));
    }

    /**
     * GROUP 3: Test fromXML
     */
    @Test
    void testFromXml() {
        // Create XML element
        Element xml = new Element("TaxTriageResultDocument");
        xml.addContent(new Element("name").setText("Restored Doc"));
        xml.addContent(new Element("resultType").setText("Restored Type"));
        xml.addContent(new Element("fileFormat").setText("HTML"));
        xml.addContent(new Element("filePath").setText("/restored/path"));
        xml.addContent(new Element("creationDate").setText("1000000000000"));
        xml.addContent(new Element("content").setText("Restored content"));

        TaxTriageResultDocument restoredDoc = new TaxTriageResultDocument();
        restoredDoc.fromXML(xml);

        assertEquals("Restored Doc", restoredDoc.getName());
        assertEquals("Restored Type", restoredDoc.getResultType());
        assertEquals("HTML", restoredDoc.getFileFormat());
        assertEquals("Restored content", restoredDoc.getContent());
        assertNotNull(restoredDoc.getCreationDate());
    }

    /**
     * GROUP 3: Test fromXML with invalid date
     */
    @Test
    void testFromXmlWithInvalidDate() {
        Element xml = new Element("TaxTriageResultDocument");
        xml.addContent(new Element("name").setText("Test"));
        xml.addContent(new Element("resultType").setText("Type"));
        xml.addContent(new Element("fileFormat").setText("TXT"));
        xml.addContent(new Element("filePath").setText("/path"));
        xml.addContent(new Element("creationDate").setText("invalid"));
        xml.addContent(new Element("content").setText("Content"));

        TaxTriageResultDocument doc = new TaxTriageResultDocument();
        doc.fromXML(xml);

        assertNotNull(doc.getCreationDate(),
            "Should set current date if parsing fails");
    }

    /**
     * GROUP 3: Test getDisplayableFields
     */
    @Test
    void testGetDisplayableFields() {
        List<DocumentField> fields = document.getDisplayableFields();

        assertNotNull(fields);
        assertTrue(fields.size() >= 3, "Should have at least 3 displayable fields");

        // Verify field types exist
        boolean hasResultType = false;
        boolean hasFileFormat = false;
        boolean hasLineCount = false;

        for (DocumentField field : fields) {
            String name = field.getName();
            if (name.contains("Result Type")) hasResultType = true;
            if (name.contains("File Format")) hasFileFormat = true;
            if (name.contains("Line Count")) hasLineCount = true;
        }

        assertTrue(hasResultType, "Should have Result Type field");
        assertTrue(hasFileFormat, "Should have File Format field");
        assertTrue(hasLineCount, "Should have Line Count field");
    }

    /**
     * GROUP 3: Test getFieldValue for result type
     */
    @Test
    void testGetFieldValueResultType() {
        Object value = document.getFieldValue("RESULT_TYPE");
        assertEquals("Test Type", value);
    }

    /**
     * GROUP 3: Test getFieldValue for file format
     */
    @Test
    void testGetFieldValueFileFormat() {
        Object value = document.getFieldValue("FILE_FORMAT");
        assertEquals("TXT", value);
    }

    /**
     * GROUP 3: Test getFieldValue for line count
     */
    @Test
    void testGetFieldValueLineCount() {
        TaxTriageResultDocument multiLineDoc = new TaxTriageResultDocument(
            "Multi Line",
            "Line 1\nLine 2\nLine 3",
            "Test",
            "/path",
            "TXT"
        );

        Object value = multiLineDoc.getFieldValue("LINE_COUNT");
        assertEquals(3, value);
    }

    /**
     * GROUP 3: Test getFieldValue for unknown field
     */
    @Test
    void testGetFieldValueUnknownField() {
        Object value = document.getFieldValue("UNKNOWN_FIELD");
        assertNull(value, "Unknown field should return null");
    }

    /**
     * GROUP 3: Test line count with empty content
     */
    @Test
    void testLineCountWithEmptyContent() {
        TaxTriageResultDocument emptyDoc = new TaxTriageResultDocument(
            "Empty",
            "",
            "Test",
            "/path",
            "TXT"
        );

        Object lineCount = emptyDoc.getFieldValue("LINE_COUNT");
        assertEquals(0, lineCount);
    }

    /**
     * GROUP 3: Test line count with null content
     */
    @Test
    void testLineCountWithNullContent() {
        TaxTriageResultDocument nullDoc = new TaxTriageResultDocument(
            "Null",
            null,
            "Test",
            "/path",
            "TXT"
        );

        Object lineCount = nullDoc.getFieldValue("LINE_COUNT");
        assertEquals(0, lineCount);
    }

    /**
     * GROUP 3: Test table HTML with empty lines
     */
    @Test
    void testTableHtmlWithEmptyLines() {
        String tsvContent = "Header1\tHeader2\n" +
                           "\n" +  // Empty line
                           "Value1\tValue2\n" +
                           "   \n" +  // Whitespace line
                           "Value3\tValue4";

        TaxTriageResultDocument doc = new TaxTriageResultDocument(
            "TSV",
            tsvContent,
            "Table",
            "/path",
            "TSV"
        );

        String html = doc.toHTML();
        assertNotNull(html);
        // Should skip empty lines
        assertTrue(html.contains("Value1"));
        assertTrue(html.contains("Value3"));
    }

    /**
     * GROUP 3: Test document with all supported formats
     */
    @Test
    void testAllSupportedFormats() {
        String[] formats = {"TXT", "TSV", "CSV", "HTML", "JSON"};

        for (String format : formats) {
            TaxTriageResultDocument doc = new TaxTriageResultDocument(
                "Test " + format,
                "Content for " + format,
                "Type",
                "/path",
                format
            );

            assertNotNull(doc.toHTML(),
                "Should generate HTML for " + format);
            assertEquals(format, doc.getFileFormat());
        }
    }

    /**
     * GROUP 3: Test getContent
     */
    @Test
    void testGetContent() {
        assertEquals("Sample content", document.getContent());
    }

    /**
     * GROUP 3: Test getResultType
     */
    @Test
    void testGetResultType() {
        assertEquals("Test Type", document.getResultType());
    }

    /**
     * GROUP 3: Test getFileFormat
     */
    @Test
    void testGetFileFormat() {
        assertEquals("TXT", document.getFileFormat());
    }

    /**
     * GROUP 3: Test XML round-trip
     */
    @Test
    void testXmlRoundTrip() {
        // Convert to XML
        Element xml = document.toXML();

        // Create new document from XML
        TaxTriageResultDocument restoredDoc = new TaxTriageResultDocument();
        restoredDoc.fromXML(xml);

        // Verify all fields match
        assertEquals(document.getName(), restoredDoc.getName());
        assertEquals(document.getContent(), restoredDoc.getContent());
        assertEquals(document.getResultType(), restoredDoc.getResultType());
        assertEquals(document.getFileFormat(), restoredDoc.getFileFormat());
    }

    /**
     * GROUP 3: Test table styling in HTML output
     */
    @Test
    void testTableStyling() {
        TaxTriageResultDocument tsvDoc = new TaxTriageResultDocument(
            "TSV",
            "H1\tH2\nV1\tV2",
            "Table",
            "/path",
            "TSV"
        );

        String html = tsvDoc.toHTML();
        assertTrue(html.contains("<style>"));
        assertTrue(html.contains("border-collapse"));
        assertTrue(html.contains("background-color"));
    }

    /**
     * GROUP 3: Test field value for file path
     */
    @Test
    void testGetFieldValueFilePath() {
        Object value = document.getFieldValue("FILE_PATH");
        assertEquals("/path/to/file.txt", value);
    }
}
