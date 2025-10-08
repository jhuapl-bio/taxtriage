package com.jhuapl.taxtriage.geneious.documents;

import com.biomatters.geneious.publicapi.documents.DocumentField;
import com.biomatters.geneious.publicapi.documents.PluginDocument;
import com.biomatters.geneious.publicapi.documents.URN;
import org.jdom.Element;
import org.jdom.CDATA;

import java.util.Arrays;
import java.util.Date;
import java.util.List;

/**
 * A document representing TaxTriage analysis results.
 * This can contain various types of results including reports, tables, and HTML content.
 */
public class TaxTriageResultDocument implements PluginDocument {

    private String name;
    private String content;
    private String resultType;
    private Date creationDate;
    private String filePath;
    private String fileFormat; // e.g., "TSV", "CSV", "HTML", "JSON", "TXT"

    // Document field codes
    private static final String KEY_RESULT_TYPE = "RESULT_TYPE";
    private static final String KEY_FILE_FORMAT = "FILE_FORMAT";
    private static final String KEY_FILE_PATH = "FILE_PATH";
    private static final String KEY_LINE_COUNT = "LINE_COUNT";

    /**
     * Constructor for creating a new result document.
     */
    public TaxTriageResultDocument(String name, String content, String resultType,
                                   String filePath, String fileFormat) {
        this.name = name;
        this.content = content;
        this.resultType = resultType;
        this.filePath = filePath;
        this.fileFormat = fileFormat;
        this.creationDate = new Date();
    }

    /**
     * Empty constructor required by PluginDocument interface.
     */
    public TaxTriageResultDocument() {
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public URN getURN() {
        return null;
    }

    @Override
    public Date getCreationDate() {
        return creationDate;
    }

    @Override
    public String getDescription() {
        return "TaxTriage " + resultType;
    }

    @Override
    public String toHTML() {
        if ("HTML".equals(fileFormat)) {
            // If it's already HTML, return as-is (with some safety escaping for script tags)
            return content.replaceAll("(?i)<script[^>]*>.*?</script>", "");
        } else if ("TSV".equals(fileFormat) || "CSV".equals(fileFormat)) {
            // Convert table to HTML
            return convertTableToHtml();
        } else {
            // For text files, wrap in pre tags
            return "<pre>" + escapeHtml(content) + "</pre>";
        }
    }

    @Override
    public Element toXML() {
        Element root = new Element("TaxTriageResultDocument");
        root.addContent(new Element("name").setText(name));
        root.addContent(new Element("resultType").setText(resultType));
        root.addContent(new Element("fileFormat").setText(fileFormat));
        root.addContent(new Element("filePath").setText(filePath != null ? filePath : ""));
        root.addContent(new Element("creationDate").setText(String.valueOf(creationDate.getTime())));
        root.addContent(new Element("content").setContent(new CDATA(content)));
        return root;
    }

    @Override
    public void fromXML(Element element) {
        this.name = element.getChildText("name");
        this.resultType = element.getChildText("resultType");
        this.fileFormat = element.getChildText("fileFormat");
        this.filePath = element.getChildText("filePath");

        String dateText = element.getChildText("creationDate");
        if (dateText != null) {
            try {
                this.creationDate = new Date(Long.parseLong(dateText));
            } catch (NumberFormatException e) {
                this.creationDate = new Date();
            }
        }

        this.content = element.getChildText("content");
    }

    @Override
    public List<DocumentField> getDisplayableFields() {
        DocumentField resultTypeField = DocumentField.createStringField("Result Type",
                "Type of TaxTriage result", KEY_RESULT_TYPE);
        DocumentField fileFormatField = DocumentField.createStringField("File Format",
                "Format of the result file", KEY_FILE_FORMAT);
        DocumentField filePathField = DocumentField.createStringField("Source File",
                "Original file path", KEY_FILE_PATH);
        DocumentField lineCountField = DocumentField.createIntegerField("Line Count",
                "Number of lines in the document", KEY_LINE_COUNT);

        return Arrays.asList(resultTypeField, fileFormatField, lineCountField);
    }

    @Override
    public Object getFieldValue(String fieldCode) {
        switch (fieldCode) {
            case KEY_RESULT_TYPE:
                return resultType;
            case KEY_FILE_FORMAT:
                return fileFormat;
            case KEY_FILE_PATH:
                return filePath;
            case KEY_LINE_COUNT:
                return getLineCount();
            default:
                return null;
        }
    }

    /**
     * Converts table content (TSV/CSV) to HTML table.
     */
    private String convertTableToHtml() {
        StringBuilder html = new StringBuilder();
        html.append("<style>")
            .append("table { border-collapse: collapse; width: 100%; }")
            .append("th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }")
            .append("th { background-color: #f2f2f2; font-weight: bold; }")
            .append("tr:nth-child(even) { background-color: #f9f9f9; }")
            .append("</style>");

        html.append("<table>");

        String delimiter = "TSV".equals(fileFormat) ? "\t" : ",";
        String[] lines = content.split("\n");

        boolean firstLine = true;
        for (String line : lines) {
            if (line.trim().isEmpty()) continue;

            html.append("<tr>");
            String[] cells = line.split(delimiter);

            for (String cell : cells) {
                if (firstLine) {
                    html.append("<th>").append(escapeHtml(cell)).append("</th>");
                } else {
                    html.append("<td>").append(escapeHtml(cell)).append("</td>");
                }
            }
            html.append("</tr>");

            firstLine = false;
        }

        html.append("</table>");
        return html.toString();
    }

    /**
     * Escapes HTML special characters.
     */
    private String escapeHtml(String text) {
        if (text == null) return "";
        return text.replace("&", "&amp;")
                  .replace("<", "&lt;")
                  .replace(">", "&gt;")
                  .replace("\"", "&quot;")
                  .replace("'", "&#39;");
    }

    /**
     * Gets the number of lines in the content.
     */
    private int getLineCount() {
        if (content == null || content.isEmpty()) return 0;
        return content.split("\n").length;
    }

    // Getters for internal use
    public String getContent() {
        return content;
    }

    public String getResultType() {
        return resultType;
    }

    public String getFileFormat() {
        return fileFormat;
    }
}