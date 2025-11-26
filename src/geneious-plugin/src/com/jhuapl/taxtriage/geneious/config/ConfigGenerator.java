package com.jhuapl.taxtriage.geneious.config;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

/**
 * Generates Nextflow configuration and parameter files for TaxTriage workflows.
 *
 * This class is responsible for creating the necessary configuration files
 * that Nextflow needs to execute TaxTriage workflows with the correct
 * parameters and Docker settings.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class ConfigGenerator {

    /**
     * Generates a Nextflow configuration file.
     *
     * @param config the TaxTriage configuration
     * @param outputFile the output configuration file
     * @throws IOException if file writing fails
     */
    public void generateConfig(TaxTriageConfig config, File outputFile) throws IOException {
        try (FileWriter writer = new FileWriter(outputFile)) {
            writer.write("// TaxTriage Nextflow Configuration\n");
            writer.write("// Generated automatically - do not edit manually\n\n");

            // Docker configuration
            writer.write("docker {\n");
            writer.write("    enabled = true\n");
            writer.write("    runOptions = '-u $(id -u):$(id -g)'\n");
            writer.write("}\n\n");

            // Process configuration
            writer.write("process {\n");
            writer.write("    cpus = " + config.getThreadCount() + "\n");
            writer.write("    memory = '" + config.getMemoryLimitGb() + ".GB'\n");
            writer.write("    container = 'taxtriage:latest'\n");
            writer.write("}\n\n");

            // Executor configuration
            writer.write("executor {\n");
            writer.write("    name = 'local'\n");
            writer.write("    cpus = " + config.getThreadCount() + "\n");
            writer.write("}\n\n");

            // Singularity configuration (alternative to Docker)
            writer.write("singularity {\n");
            writer.write("    enabled = false\n");
            writer.write("    autoMounts = true\n");
            writer.write("}\n\n");

            // Timeline and report configuration
            writer.write("timeline {\n");
            writer.write("    enabled = true\n");
            writer.write("    file = 'timeline.html'\n");
            writer.write("}\n\n");

            writer.write("report {\n");
            writer.write("    enabled = true\n");
            writer.write("    file = 'report.html'\n");
            writer.write("}\n\n");

            writer.write("trace {\n");
            writer.write("    enabled = true\n");
            writer.write("    file = 'trace.txt'\n");
            writer.write("}\n");
        }
    }

    /**
     * Generates a Nextflow parameters JSON file.
     *
     * @param config the TaxTriage configuration
     * @param outputFile the output parameters file
     * @throws IOException if file writing fails
     */
    public void generateParams(TaxTriageConfig config, File outputFile) throws IOException {
        Map<String, Object> params = config.toParameterMap();

        // Use a simple JSON generation approach to avoid external dependencies
        try (FileWriter writer = new FileWriter(outputFile)) {
            writer.write("{\n");

            String[] entries = new String[params.size()];
            int i = 0;
            for (Map.Entry<String, Object> entry : params.entrySet()) {
                String key = entry.getKey();
                Object value = entry.getValue();
                String jsonValue;

                if (value instanceof String) {
                    jsonValue = "\"" + escapeJson((String) value) + "\"";
                } else if (value instanceof Boolean) {
                    jsonValue = value.toString();
                } else if (value instanceof Number) {
                    jsonValue = value.toString();
                } else {
                    jsonValue = "\"" + value.toString() + "\"";
                }

                entries[i++] = "  \"" + key + "\": " + jsonValue;
            }

            for (int j = 0; j < entries.length; j++) {
                writer.write(entries[j]);
                if (j < entries.length - 1) {
                    writer.write(",");
                }
                writer.write("\n");
            }

            writer.write("}\n");
        }
    }

    /**
     * Escapes special characters in JSON strings.
     *
     * @param input the input string
     * @return escaped string
     */
    private String escapeJson(String input) {
        if (input == null) {
            return "";
        }
        return input.replace("\\", "\\\\")
                   .replace("\"", "\\\"")
                   .replace("\n", "\\n")
                   .replace("\r", "\\r")
                   .replace("\t", "\\t");
    }
}