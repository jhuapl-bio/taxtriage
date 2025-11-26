package com.jhuapl.taxtriage.geneious.reports;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.logging.Logger;

public class TaxTriageHtmlReportGenerator {

    private static final Logger logger = Logger.getLogger(TaxTriageHtmlReportGenerator.class.getName());

    private String runName;
    private Path workingDirectory;
    private Date startTime;
    private Date endTime;
    private List<String> samples;
    private String commandLine;
    private List<String> errors;
    private List<String> warnings;
    private Map<String, String> parameters;
    private Map<String, Object> statistics;
    private Map<String, String> toolVersions;
    private String nextflowCommand;
    private Path outputDirectory;

    public TaxTriageHtmlReportGenerator(String runName, Path workingDirectory) {
        this.runName = runName;
        this.workingDirectory = workingDirectory;
        this.startTime = new Date();
        this.samples = new ArrayList<>();
        this.errors = new ArrayList<>();
        this.warnings = new ArrayList<>();
        this.parameters = new LinkedHashMap<>();
        this.statistics = new LinkedHashMap<>();
        this.toolVersions = new LinkedHashMap<>();
    }

    public void setEndTime(Date endTime) {
        this.endTime = endTime;
    }

    public void addSample(String sample) {
        this.samples.add(sample);
    }

    public void setCommandLine(String commandLine) {
        this.commandLine = commandLine;
    }

    public void addError(String error) {
        this.errors.add(error);
    }

    public void addWarning(String warning) {
        this.warnings.add(warning);
    }

    public void addParameter(String name, String value) {
        this.parameters.put(name, value);
    }

    public void addStatistic(String name, Object value) {
        this.statistics.put(name, value);
    }

    public void addToolVersion(String tool, String version) {
        this.toolVersions.put(tool, version);
    }

    public void setNextflowCommand(String command) {
        this.nextflowCommand = command;
    }

    public void setOutputDirectory(Path outputDir) {
        this.outputDirectory = outputDir;
    }

    public String generateHtml() {
        SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

        StringBuilder html = new StringBuilder();
        html.append("<!DOCTYPE html>\n");
        html.append("<html>\n");
        html.append("<head>\n");
        html.append("    <meta charset=\"UTF-8\">\n");
        html.append("    <title>TaxTriage Analysis Overview - ").append(escapeHtml(runName)).append("</title>\n");
        html.append("    <style>\n");
        html.append(getStyleSheet());
        html.append("    </style>\n");
        html.append("</head>\n");
        html.append("<body>\n");

        html.append("    <div class=\"header\">\n");
        html.append("        <h1>TaxTriage Analysis Overview</h1>\n");
        html.append("        <div class=\"run-name\">").append(escapeHtml(runName)).append("</div>\n");
        html.append("    </div>\n");

        html.append("    <div class=\"container\">\n");

        html.append("        <section class=\"card\">\n");
        html.append("            <h2>Analysis Summary</h2>\n");
        html.append("            <table>\n");
        html.append("                <tr><td class=\"label\">Start Time:</td><td>").append(dateFormat.format(startTime)).append("</td></tr>\n");
        if (endTime != null) {
            html.append("                <tr><td class=\"label\">End Time:</td><td>").append(dateFormat.format(endTime)).append("</td></tr>\n");
            long durationMs = endTime.getTime() - startTime.getTime();
            html.append("                <tr><td class=\"label\">Duration:</td><td>").append(formatDuration(durationMs)).append("</td></tr>\n");
        }
        html.append("                <tr><td class=\"label\">Working Directory:</td><td class=\"path\">").append(escapeHtml(workingDirectory.toString())).append("</td></tr>\n");
        html.append("                <tr><td class=\"label\">Samples Processed:</td><td>").append(samples.size()).append("</td></tr>\n");
        html.append("            </table>\n");
        html.append("        </section>\n");

        if (!samples.isEmpty()) {
            html.append("        <section class=\"card\">\n");
            html.append("            <h2>Samples</h2>\n");
            html.append("            <ul class=\"sample-list\">\n");
            for (String sample : samples) {
                html.append("                <li>").append(escapeHtml(sample)).append("</li>\n");
            }
            html.append("            </ul>\n");
            html.append("        </section>\n");
        }

        if (!parameters.isEmpty()) {
            html.append("        <section class=\"card\">\n");
            html.append("            <h2>Analysis Parameters</h2>\n");
            html.append("            <table>\n");
            for (Map.Entry<String, String> entry : parameters.entrySet()) {
                html.append("                <tr><td class=\"label\">").append(escapeHtml(entry.getKey())).append(":</td>");
                html.append("<td>").append(escapeHtml(entry.getValue())).append("</td></tr>\n");
            }
            html.append("            </table>\n");
            html.append("        </section>\n");
        }

        if (!toolVersions.isEmpty()) {
            html.append("        <section class=\"card\">\n");
            html.append("            <h2>Tool Versions</h2>\n");
            html.append("            <table>\n");
            for (Map.Entry<String, String> entry : toolVersions.entrySet()) {
                html.append("                <tr><td class=\"label\">").append(escapeHtml(entry.getKey())).append(":</td>");
                html.append("<td>").append(escapeHtml(entry.getValue())).append("</td></tr>\n");
            }
            html.append("            </table>\n");
            html.append("        </section>\n");
        }

        if ((commandLine != null && !commandLine.isEmpty()) || (nextflowCommand != null && !nextflowCommand.isEmpty())) {
            html.append("        <section class=\"card\">\n");
            html.append("            <h2>Command Line Replication</h2>\n");

            if (workingDirectory != null) {
                html.append("            <div class=\"info-section\">\n");
                html.append("                <p><strong>Working Directory:</strong> <code>").append(escapeHtml(workingDirectory.toString())).append("</code></p>\n");
                html.append("            </div>\n");
            }

            if (nextflowCommand != null && !nextflowCommand.isEmpty()) {
                html.append("            <div class=\"command-box\">\n");
                html.append("                <h3>Nextflow Command</h3>\n");
                html.append("                <pre>").append(escapeHtml(nextflowCommand)).append("</pre>\n");
                html.append("            </div>\n");
            }

            if (commandLine != null && !commandLine.isEmpty() && !commandLine.equals(nextflowCommand)) {
                html.append("            <div class=\"command-box\">\n");
                html.append("                <h3>Full Command Line</h3>\n");
                html.append("                <pre>").append(escapeHtml(commandLine)).append("</pre>\n");
                html.append("            </div>\n");
            }

            if (outputDirectory != null) {
                html.append("            <div class=\"info-section\">\n");
                html.append("                <p><strong>Output Directory:</strong> <code>").append(escapeHtml(outputDirectory.toString())).append("</code></p>\n");
                html.append("            </div>\n");
            }

            html.append("            <p class=\"help-text\">Copy and execute the Nextflow command in the working directory to replicate the analysis.</p>\n");
            html.append("        </section>\n");
        }

        if (!statistics.isEmpty()) {
            html.append("        <section class=\"card\">\n");
            html.append("            <h2>Statistics</h2>\n");
            html.append("            <table>\n");
            for (Map.Entry<String, Object> entry : statistics.entrySet()) {
                html.append("                <tr><td class=\"label\">").append(escapeHtml(entry.getKey())).append(":</td>");
                html.append("<td>").append(escapeHtml(String.valueOf(entry.getValue()))).append("</td></tr>\n");
            }
            html.append("            </table>\n");
            html.append("        </section>\n");
        }

        if (!warnings.isEmpty()) {
            html.append("        <section class=\"card warning-card\">\n");
            html.append("            <h2>⚠️ Warnings (").append(warnings.size()).append(")</h2>\n");
            html.append("            <ul class=\"message-list\">\n");
            for (String warning : warnings) {
                html.append("                <li>").append(escapeHtml(warning)).append("</li>\n");
            }
            html.append("            </ul>\n");
            html.append("        </section>\n");
        }

        if (!errors.isEmpty()) {
            html.append("        <section class=\"card error-card\">\n");
            html.append("            <h2>❌ Errors (").append(errors.size()).append(")</h2>\n");
            html.append("            <ul class=\"message-list\">\n");
            for (String error : errors) {
                html.append("                <li>").append(escapeHtml(error)).append("</li>\n");
            }
            html.append("            </ul>\n");
            html.append("        </section>\n");
        }

        html.append("    </div>\n");

        html.append("    <div class=\"footer\">\n");
        html.append("        <p>Generated by TaxTriage Geneious Plugin</p>\n");
        html.append("        <p>Report generated on: ").append(dateFormat.format(new Date())).append("</p>\n");
        html.append("    </div>\n");

        html.append("</body>\n");
        html.append("</html>\n");

        return html.toString();
    }

    public void saveToFile(Path outputPath) throws IOException {
        String htmlContent = generateHtml();
        Files.writeString(outputPath, htmlContent);
        logger.info("HTML report saved to: " + outputPath);
    }

    private String getStyleSheet() {
        return "body {" + "\n" +
                "    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;" + "\n" +
                "    line-height: 1.6;" + "\n" +
                "    color: #333;" + "\n" +
                "    margin: 0;" + "\n" +
                "    padding: 0;" + "\n" +
                "    background: #f5f5f5;" + "\n" +
                "}" + "\n" +
                ".header {" + "\n" +
                "    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);" + "\n" +
                "    color: white;" + "\n" +
                "    padding: 2rem;" + "\n" +
                "    box-shadow: 0 2px 4px rgba(0,0,0,0.1);" + "\n" +
                "}" + "\n" +
                ".header h1 {" + "\n" +
                "    margin: 0;" + "\n" +
                "    font-size: 2rem;" + "\n" +
                "}" + "\n" +
                ".run-name {" + "\n" +
                "    font-size: 1.2rem;" + "\n" +
                "    opacity: 0.9;" + "\n" +
                "    margin-top: 0.5rem;" + "\n" +
                "}" + "\n" +
                ".container {" + "\n" +
                "    max-width: 1200px;" + "\n" +
                "    margin: 2rem auto;" + "\n" +
                "    padding: 0 1rem;" + "\n" +
                "}" + "\n" +
                ".card {" + "\n" +
                "    background: white;" + "\n" +
                "    border-radius: 8px;" + "\n" +
                "    padding: 1.5rem;" + "\n" +
                "    margin-bottom: 1.5rem;" + "\n" +
                "    box-shadow: 0 2px 4px rgba(0,0,0,0.1);" + "\n" +
                "}" + "\n" +
                ".card h2 {" + "\n" +
                "    margin-top: 0;" + "\n" +
                "    color: #667eea;" + "\n" +
                "    border-bottom: 2px solid #f0f0f0;" + "\n" +
                "    padding-bottom: 0.5rem;" + "\n" +
                "}" + "\n" +
                "table {" + "\n" +
                "    width: 100%;" + "\n" +
                "    border-collapse: collapse;" + "\n" +
                "}" + "\n" +
                "td {" + "\n" +
                "    padding: 0.5rem;" + "\n" +
                "    border-bottom: 1px solid #f0f0f0;" + "\n" +
                "}" + "\n" +
                ".label {" + "\n" +
                "    font-weight: 600;" + "\n" +
                "    width: 200px;" + "\n" +
                "    color: #666;" + "\n" +
                "}" + "\n" +
                ".path {" + "\n" +
                "    font-family: 'Courier New', monospace;" + "\n" +
                "    font-size: 0.9rem;" + "\n" +
                "    background: #f8f8f8;" + "\n" +
                "    padding: 0.25rem 0.5rem;" + "\n" +
                "    border-radius: 4px;" + "\n" +
                "}" + "\n" +
                ".sample-list {" + "\n" +
                "    list-style: none;" + "\n" +
                "    padding: 0;" + "\n" +
                "    display: grid;" + "\n" +
                "    grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));" + "\n" +
                "    gap: 0.5rem;" + "\n" +
                "}" + "\n" +
                ".sample-list li {" + "\n" +
                "    background: #f8f9fa;" + "\n" +
                "    padding: 0.5rem 1rem;" + "\n" +
                "    border-radius: 4px;" + "\n" +
                "    border-left: 3px solid #667eea;" + "\n" +
                "}" + "\n" +
                ".command-box {" + "\n" +
                "    background: #2d2d2d;" + "\n" +
                "    color: #f8f8f2;" + "\n" +
                "    padding: 1rem;" + "\n" +
                "    border-radius: 4px;" + "\n" +
                "    overflow-x: auto;" + "\n" +
                "}" + "\n" +
                ".command-box pre {" + "\n" +
                "    margin: 0;" + "\n" +
                "    font-family: 'Courier New', monospace;" + "\n" +
                "    font-size: 0.9rem;" + "\n" +
                "    white-space: pre-wrap;" + "\n" +
                "    word-wrap: break-word;" + "\n" +
                "}" + "\n" +
                ".help-text {" + "\n" +
                "    color: #666;" + "\n" +
                "    font-size: 0.9rem;" + "\n" +
                "    margin-top: 0.5rem;" + "\n" +
                "}" + "\n" +
                ".info-section {" + "\n" +
                "    background: #f8f9fa;" + "\n" +
                "    padding: 1rem;" + "\n" +
                "    border-radius: 4px;" + "\n" +
                "    margin: 1rem 0;" + "\n" +
                "}" + "\n" +
                ".command-box h3 {" + "\n" +
                "    color: #f8f8f2;" + "\n" +
                "    margin: 0 0 0.5rem 0;" + "\n" +
                "    font-size: 1rem;" + "\n" +
                "}" + "\n" +
                ".message-list {" + "\n" +
                "    list-style: none;" + "\n" +
                "    padding: 0;" + "\n" +
                "}" + "\n" +
                ".message-list li {" + "\n" +
                "    padding: 0.5rem;" + "\n" +
                "    margin-bottom: 0.5rem;" + "\n" +
                "    border-radius: 4px;" + "\n" +
                "    background: #f8f9fa;" + "\n" +
                "}" + "\n" +
                ".warning-card {" + "\n" +
                "    border-left: 4px solid #ff9800;" + "\n" +
                "}" + "\n" +
                ".warning-card h2 {" + "\n" +
                "    color: #ff9800;" + "\n" +
                "}" + "\n" +
                ".error-card {" + "\n" +
                "    border-left: 4px solid #f44336;" + "\n" +
                "}" + "\n" +
                ".error-card h2 {" + "\n" +
                "    color: #f44336;" + "\n" +
                "}" + "\n" +
                ".footer {" + "\n" +
                "    text-align: center;" + "\n" +
                "    padding: 2rem;" + "\n" +
                "    color: #666;" + "\n" +
                "    font-size: 0.9rem;" + "\n" +
                "}";
    }

    private String escapeHtml(String text) {
        if (text == null) return "";
        return text.replace("&", "&amp;")
                   .replace("<", "&lt;")
                   .replace(">", "&gt;")
                   .replace("\"", "&quot;")
                   .replace("'", "&#39;");
    }

    private String formatDuration(long milliseconds) {
        long seconds = milliseconds / 1000;
        long minutes = seconds / 60;
        long hours = minutes / 60;

        if (hours > 0) {
            return String.format("%d hours, %d minutes", hours, minutes % 60);
        } else if (minutes > 0) {
            return String.format("%d minutes, %d seconds", minutes, seconds % 60);
        } else {
            return String.format("%d seconds", seconds);
        }
    }
}