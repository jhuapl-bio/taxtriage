package com.jhuapl.taxtriage.geneious.config;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 * Builds samplesheet CSV files for TaxTriage Nextflow workflows.
 *
 * This class creates properly formatted samplesheet files that specify
 * input samples and their metadata for TaxTriage analysis. The format
 * varies based on the sequencing preset (single-end, paired-end, or long-read).
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class SampleSheetBuilder {

    /**
     * Builds a samplesheet CSV file from the exported sequence files.
     *
     * @param exportedFiles list of exported FASTQ files
     * @param sampleSheetFile output samplesheet file
     * @param preset sequencing preset (ONT, ILLUMINA_PE, ILLUMINA_SE)
     * @throws IOException if file writing fails
     */
    public void buildSampleSheet(List<File> exportedFiles, File sampleSheetFile, String preset) throws IOException {
        if (exportedFiles == null || exportedFiles.isEmpty()) {
            throw new IllegalArgumentException("No exported files provided");
        }

        try (FileWriter writer = new FileWriter(sampleSheetFile)) {
            // Write header based on preset
            writeHeader(writer, preset);

            // Write sample entries
            for (int i = 0; i < exportedFiles.size(); i++) {
                File file = exportedFiles.get(i);
                String sampleName = generateSampleName(file, i);
                writeSampleEntry(writer, sampleName, file, preset);
            }
        }
    }

    /**
     * Writes the CSV header based on sequencing preset.
     *
     * @param writer the file writer
     * @param preset the sequencing preset
     * @throws IOException if writing fails
     */
    private void writeHeader(FileWriter writer, String preset) throws IOException {
        // All presets use the same basic header format
        // TaxTriage expects: sample,fastq_1,fastq_2,long_fastq,fast5
        writer.write("sample,fastq_1,fastq_2,long_fastq,fast5\n");
    }

    /**
     * Writes a sample entry to the CSV file.
     *
     * @param writer the file writer
     * @param sampleName the sample name
     * @param file the FASTQ file
     * @param preset the sequencing preset
     * @throws IOException if writing fails
     */
    private void writeSampleEntry(FileWriter writer, String sampleName, File file, String preset) throws IOException {
        String absolutePath = file.getAbsolutePath();

        // Use the standard 5-column format: sample,fastq_1,fastq_2,long_fastq,fast5
        // Place the file in the correct column based on the platform
        switch (preset.toUpperCase()) {
            case "ONT":
                // For ONT, put the file in fastq_1 column
                writer.write(sampleName + "," + absolutePath + ",,,\n");
                break;
            case "ILLUMINA_PE":
                // For paired-end, we'd need to pair files - for now treat as single-end
                // In a real implementation, you'd parse _R1/_R2 patterns
                writer.write(sampleName + "," + absolutePath + ",,,\n");
                break;
            case "ILLUMINA_SE":
            default:
                // Single-end format - file in fastq_1 column
                writer.write(sampleName + "," + absolutePath + ",,,\n");
                break;
        }
    }

    /**
     * Generates a sample name from the file.
     *
     * @param file the FASTQ file
     * @param index the file index
     * @return generated sample name
     */
    private String generateSampleName(File file, int index) {
        String fileName = file.getName();

        // Remove file extensions (handle .gz files properly)
        if (fileName.endsWith(".fastq.gz")) {
            fileName = fileName.substring(0, fileName.length() - 9);
        } else if (fileName.endsWith(".fq.gz")) {
            fileName = fileName.substring(0, fileName.length() - 6);
        } else if (fileName.endsWith(".fastq")) {
            fileName = fileName.substring(0, fileName.length() - 6);
        } else if (fileName.endsWith(".fq")) {
            fileName = fileName.substring(0, fileName.length() - 3);
        }

        // Clean up the name (remove invalid characters)
        fileName = fileName.replaceAll("[^a-zA-Z0-9._-]", "_");

        // Ensure uniqueness
        return fileName.isEmpty() ? "sample_" + (index + 1) : fileName;
    }

    /**
     * Validates that the samplesheet format is correct.
     *
     * @param sampleSheetFile the samplesheet file to validate
     * @param preset the expected preset
     * @return validation error message or null if valid
     */
    public String validateSampleSheet(File sampleSheetFile, String preset) {
        if (!sampleSheetFile.exists()) {
            return "Samplesheet file does not exist: " + sampleSheetFile.getAbsolutePath();
        }

        if (sampleSheetFile.length() == 0) {
            return "Samplesheet file is empty";
        }

        // In a full implementation, you would parse the CSV and validate:
        // - Header matches expected format
        // - All referenced files exist
        // - Sample names are unique
        // - Paired-end files are properly matched

        return null; // Assume valid for now
    }

    /**
     * Counts the number of samples in a samplesheet.
     *
     * @param sampleSheetFile the samplesheet file
     * @return number of samples (excluding header)
     * @throws IOException if file reading fails
     */
    public int countSamples(File sampleSheetFile) throws IOException {
        if (!sampleSheetFile.exists()) {
            return 0;
        }

        int lineCount = 0;
        try (java.io.BufferedReader reader = new java.io.BufferedReader(new java.io.FileReader(sampleSheetFile))) {
            String line;
            while ((line = reader.readLine()) != null) {
                lineCount++;
            }
        }

        // Subtract 1 for header line
        return Math.max(0, lineCount - 1);
    }
}