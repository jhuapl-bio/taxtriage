package com.jhuapl.taxtriage.geneious.docker;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.sequence.SequenceDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentOperationException;
import jebl.util.ProgressListener;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.List;

/**
 * Integration helper class for using DockerManager with TaxTriage workflows.
 *
 * This class demonstrates how to integrate the Docker execution manager
 * with the existing TaxTriageOperation to run actual workflows.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class TaxTriageDockerIntegration {

    private final DockerManager dockerManager;

    /**
     * Creates a new integration helper.
     *
     * @throws DockerException if Docker is not available
     */
    public TaxTriageDockerIntegration() throws DockerException {
        this.dockerManager = new DockerManager();
    }

    /**
     * Creates a new integration helper with a custom DockerManager.
     *
     * @param dockerManager the Docker manager to use
     */
    public TaxTriageDockerIntegration(DockerManager dockerManager) {
        this.dockerManager = dockerManager;
    }

    /**
     * Executes a TaxTriage workflow on the provided sequence documents.
     *
     * @param documents the input sequence documents
     * @param workspaceDir the workspace directory for temporary files
     * @param progressListener optional progress listener
     * @return the execution result
     * @throws DockerException if the workflow execution fails
     */
    public ExecutionResult executeTaxTriageWorkflow(
            AnnotatedPluginDocument[] documents,
            Path workspaceDir,
            ProgressListener progressListener) throws DockerException {

        try {
            // Create required directories
            Path inputDir = workspaceDir.resolve("input");
            Path outputDir = workspaceDir.resolve("output");
            Path workDir = workspaceDir.resolve("work");

            Files.createDirectories(inputDir);
            Files.createDirectories(outputDir);
            Files.createDirectories(workDir);

            if (progressListener != null) {
                progressListener.setMessage("Preparing input files...");
                progressListener.setProgress(0.1);
            }

            // Export sequence documents to FASTQ files
            List<Path> inputFiles = exportSequencesToFiles(documents, inputDir);

            if (progressListener != null) {
                progressListener.setMessage("Configuring TaxTriage workflow...");
                progressListener.setProgress(0.2);
            }

            // Build Nextflow command
            String nextflowCommand = buildNextflowCommand(inputFiles);

            if (progressListener != null) {
                progressListener.setMessage("Starting TaxTriage execution...");
                progressListener.setProgress(0.3);
            }

            // Execute the workflow
            ExecutionResult result = dockerManager.executeNextflowCommand(
                    nextflowCommand,
                    inputDir,
                    outputDir,
                    workDir,
                    progressListener
            );

            if (progressListener != null) {
                if (result.isSuccess()) {
                    progressListener.setMessage("TaxTriage workflow completed successfully");
                    progressListener.setProgress(1.0);
                } else {
                    progressListener.setMessage("TaxTriage workflow failed");
                }
            }

            return result;

        } catch (IOException e) {
            throw new DockerException("Failed to prepare workflow files", e);
        }
    }

    /**
     * Executes a TaxTriage workflow with custom parameters.
     *
     * @param documents the input sequence documents
     * @param workspaceDir the workspace directory
     * @param customParams custom Nextflow parameters
     * @param progressListener optional progress listener
     * @return the execution result
     * @throws DockerException if the workflow execution fails
     */
    public ExecutionResult executeTaxTriageWorkflow(
            AnnotatedPluginDocument[] documents,
            Path workspaceDir,
            String customParams,
            ProgressListener progressListener) throws DockerException {

        try {
            Path inputDir = workspaceDir.resolve("input");
            Path outputDir = workspaceDir.resolve("output");
            Path workDir = workspaceDir.resolve("work");

            Files.createDirectories(inputDir);
            Files.createDirectories(outputDir);
            Files.createDirectories(workDir);

            // Export sequences and build custom command
            List<Path> inputFiles = exportSequencesToFiles(documents, inputDir);
            String nextflowCommand = buildCustomNextflowCommand(inputFiles, customParams);

            return dockerManager.executeNextflowCommand(
                    nextflowCommand,
                    inputDir,
                    outputDir,
                    workDir,
                    progressListener
            );

        } catch (IOException e) {
            throw new DockerException("Failed to prepare workflow files", e);
        }
    }

    /**
     * Checks if Docker and TaxTriage are available for execution.
     *
     * @return true if everything is ready, false otherwise
     */
    public boolean isReadyForExecution() {
        try {
            return dockerManager.isDockerAvailable();
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Exports sequence documents to FASTQ files in the input directory.
     *
     * @param documents the sequence documents to export
     * @param inputDir the target input directory
     * @return list of created file paths
     * @throws IOException if file writing fails
     */
    private List<Path> exportSequencesToFiles(AnnotatedPluginDocument[] documents, Path inputDir)
            throws IOException {

        List<Path> files = new ArrayList<>();

        for (int i = 0; i < documents.length; i++) {
            AnnotatedPluginDocument doc = documents[i];
            try {
                if (!(doc.getDocument() instanceof SequenceDocument)) {
                    continue;
                }

                SequenceDocument seqDoc = (SequenceDocument) doc.getDocument();
            String fileName = sanitizeFileName(doc.getName()) + "_" + i + ".fastq";
            Path filePath = inputDir.resolve(fileName);

            // Create a simple FASTQ representation
            String fastqContent = createFastqContent(seqDoc);
                Files.write(filePath, fastqContent.getBytes(), StandardOpenOption.CREATE);

                files.add(filePath);
            } catch (DocumentOperationException e) {
                throw new IOException("Failed to access document: " + doc.getName(), e);
            }
        }

        return files;
    }

    /**
     * Builds the standard Nextflow command for TaxTriage.
     *
     * @param inputFiles the input FASTQ files
     * @return the Nextflow command string
     */
    private String buildNextflowCommand(List<Path> inputFiles) {
        StringBuilder cmd = new StringBuilder();
        cmd.append("nextflow run jhuapl-bio/taxtriage ");
        cmd.append("--input /input ");
        cmd.append("--outdir /output ");
        cmd.append("--workdir /work ");

        // Add input files as a glob pattern
        if (!inputFiles.isEmpty()) {
            cmd.append("--input_pattern '*.fastq' ");
        }

        // Add default parameters
        cmd.append("--skip_kraken2 false ");
        cmd.append("--skip_metaphlan false ");
        cmd.append("--generate_krona true ");

        return cmd.toString();
    }

    /**
     * Builds a custom Nextflow command with user-specified parameters.
     *
     * @param inputFiles the input FASTQ files
     * @param customParams custom parameters string
     * @return the Nextflow command string
     */
    private String buildCustomNextflowCommand(List<Path> inputFiles, String customParams) {
        StringBuilder cmd = new StringBuilder();
        cmd.append("nextflow run jhuapl-bio/taxtriage ");
        cmd.append("--input /input ");
        cmd.append("--outdir /output ");
        cmd.append("--workdir /work ");

        // Add custom parameters
        if (customParams != null && !customParams.trim().isEmpty()) {
            cmd.append(customParams.trim()).append(" ");
        }

        return cmd.toString();
    }

    /**
     * Creates FASTQ content from a sequence document.
     *
     * @param seqDoc the sequence document
     * @return FASTQ formatted string
     */
    private String createFastqContent(SequenceDocument seqDoc) {
        StringBuilder fastq = new StringBuilder();
        fastq.append("@").append(seqDoc.getName()).append("\n");
        fastq.append(seqDoc.getSequenceString()).append("\n");
        fastq.append("+\n");

        // Create dummy quality scores (all high quality 'I' = 40)
        int length = seqDoc.getSequenceLength();
        for (int i = 0; i < length; i++) {
            fastq.append("I");
        }
        fastq.append("\n");

        return fastq.toString();
    }

    /**
     * Sanitizes a file name by removing invalid characters.
     *
     * @param fileName the original file name
     * @return sanitized file name
     */
    private String sanitizeFileName(String fileName) {
        if (fileName == null) {
            return "sequence";
        }

        return fileName.replaceAll("[^a-zA-Z0-9._-]", "_");
    }

    /**
     * Gets the underlying DockerManager instance.
     *
     * @return the Docker manager
     */
    public DockerManager getDockerManager() {
        return dockerManager;
    }
}