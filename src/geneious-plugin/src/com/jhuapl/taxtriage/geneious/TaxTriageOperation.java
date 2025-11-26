package com.jhuapl.taxtriage.geneious;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.documents.PluginDocument;
import com.biomatters.geneious.publicapi.documents.sequence.SequenceDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentOperation;
import com.biomatters.geneious.publicapi.plugin.DocumentOperationException;
import com.biomatters.geneious.publicapi.plugin.DocumentSelectionSignature;
import com.biomatters.geneious.publicapi.plugin.GeneiousActionOptions;
import com.biomatters.geneious.publicapi.plugin.Options;
import com.biomatters.geneious.publicapi.implementations.sequence.DefaultNucleotideSequence;
import com.jhuapl.taxtriage.geneious.execution.TaxTriageWorkflowRunner;
import com.jhuapl.taxtriage.geneious.execution.WorkflowContext;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * TaxTriage document operation for taxonomic classification and pathogen detection.
 *
 * This operation processes sequence documents (FASTQ, FASTA) through the TaxTriage
 * workflow to identify pathogens and perform taxonomic classification. It supports
 * both single-end and paired-end sequencing data with configurable parameters
 * for different sequencing platforms.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
public class TaxTriageOperation extends DocumentOperation {

    private static final Logger logger = Logger.getLogger(TaxTriageOperation.class.getName());

    /**
     * Performs the TaxTriage analysis operation.
     *
     * @param documents input sequence documents
     * @param progressListener progress tracking for the operation
     * @param options configuration options for the analysis
     * @return list of result documents (reports, classifications, etc.)
     * @throws DocumentOperationException if the operation fails
     */
    @Override
    public List<AnnotatedPluginDocument> performOperation(AnnotatedPluginDocument[] documents,
                                                          ProgressListener progressListener,
                                                          Options options) throws DocumentOperationException {
        if (progressListener != null) {
            progressListener.setMessage("Initializing TaxTriage analysis...");
            progressListener.setProgress(0.0);
        }

        try {
            // Validate inputs
            if (documents == null || documents.length == 0) {
                throw new DocumentOperationException("No input documents provided");
            }

            // Cast to our specific options class
            TaxTriageOptions taxTriageOptions;
            if (options instanceof TaxTriageOptions) {
                taxTriageOptions = (TaxTriageOptions) options;
            } else {
                throw new DocumentOperationException("Invalid options type provided");
            }

            // Check if we have input files from browsers
            List<File> inputFiles = taxTriageOptions.getInputFiles();
            File inputDirectory = taxTriageOptions.getInputDirectory();
            boolean hasFileInputs = (inputFiles != null && !inputFiles.isEmpty()) || inputDirectory != null;

            // Log input sources
            logger.info("[TaxTriage] Input files from browser: " + (inputFiles != null ? inputFiles.size() : 0));
            logger.info("[TaxTriage] Input directory: " + (inputDirectory != null ? inputDirectory.getAbsolutePath() : "null"));
            logger.info("[TaxTriage] Number of selected documents: " + (documents != null ? documents.length : 0));

            if (!hasFileInputs) {
                // Only check selected documents if no file browser input was provided
                boolean hasSelectedDocuments = documents != null && documents.length > 0;

                if (!hasSelectedDocuments) {
                    throw new DocumentOperationException("No input provided. Please use the file browsers to specify input files or directory.");
                }

                // If using selected documents (backward compatibility), validate they are sequence documents
                for (AnnotatedPluginDocument doc : documents) {
                    if (!(doc.getDocument() instanceof SequenceDocument)) {
                        throw new DocumentOperationException("Document '" + doc.getName() + "' is not a sequence document. Please use the file browsers instead.");
                    }
                }
            } else {
                // When using file browsers, ignore selected documents completely
                logger.info("[TaxTriage] Using file browser inputs, ignoring any selected documents");
            }

            if (progressListener != null) {
                progressListener.setMessage("Checking system requirements...");
                progressListener.setProgress(0.05);
            }

            // Check system requirements
            TaxTriageWorkflowRunner workflowRunner;
            try {
                workflowRunner = new TaxTriageWorkflowRunner();
                TaxTriageWorkflowRunner.SystemStatus systemStatus = workflowRunner.checkSystemRequirements();

                if (!systemStatus.overall) {
                    throw new DocumentOperationException("System requirements not met:\n" + systemStatus.getSummary());
                }

                logger.info("System requirements check passed");
            } catch (Exception e) {
                throw new DocumentOperationException("Failed to initialize workflow runner: " + e.getMessage(), e);
            }

            if (progressListener != null) {
                progressListener.setMessage("Starting TaxTriage workflow...");
                progressListener.setProgress(0.1);
            }

            // Execute the workflow
            WorkflowContext workflowContext;
            try {
                workflowContext = workflowRunner.runWorkflow(documents, taxTriageOptions, progressListener);
            } catch (DocumentOperationException e) {
                throw e;
            } catch (Exception e) {
                throw new DocumentOperationException("Workflow execution failed: " + e.getMessage(), e);
            }

            // Check if workflow was successful
            if (!workflowContext.isSuccessful()) {
                String errorMessage = "Workflow failed";
                if (workflowContext.getLastError() != null) {
                    errorMessage = workflowContext.getLastError().getMessage();
                }
                throw new DocumentOperationException(errorMessage, workflowContext.getLastError());
            }

            if (progressListener != null) {
                progressListener.setMessage("Processing workflow results...");
                progressListener.setProgress(0.95);
            }

            // Convert workflow results to Geneious documents
            List<PluginDocument> results = convertWorkflowResults(workflowContext, documents);

            if (progressListener != null) {
                progressListener.setMessage("Analysis complete");
                progressListener.setProgress(1.0);
            }

            logger.info("TaxTriage analysis completed successfully for " + documents.length + " documents");
            return DocumentUtilities.createAnnotatedPluginDocuments(results);

        } catch (Exception e) {
            logger.log(Level.SEVERE, "TaxTriage analysis failed", e);
            if (e instanceof DocumentOperationException) {
                throw e;
            }
            throw new DocumentOperationException("TaxTriage analysis failed: " + e.getMessage(), e);
        }
    }

    /**
     * Converts workflow results to Geneious plugin documents.
     *
     * @param workflowContext the completed workflow context
     * @param originalDocuments the original input documents
     * @return list of result documents
     * @throws DocumentOperationException if conversion fails
     */
    private List<PluginDocument> convertWorkflowResults(WorkflowContext workflowContext,
                                                        AnnotatedPluginDocument[] originalDocuments)
            throws DocumentOperationException {
        List<PluginDocument> results = new ArrayList<>();

        try {
            // Create summary report document
            String summaryReport = generateSummaryReport(workflowContext, originalDocuments);
            DefaultNucleotideSequence summaryDoc = new DefaultNucleotideSequence(
                    "TaxTriage_Summary_Report",
                    "TaxTriage workflow summary and results",
                    "N", // Placeholder sequence for report
                    new Date()
            );
            summaryDoc.setDescription(summaryReport);
            results.add(summaryDoc);

            // Add individual result files as documents if available
            if (workflowContext.getOutputFiles() != null) {
                for (File outputFile : workflowContext.getOutputFiles()) {
                    if (isReportFile(outputFile)) {
                        try {
                            String content = Files.readString(outputFile.toPath());
                            DefaultNucleotideSequence reportDoc = new DefaultNucleotideSequence(
                                    "TaxTriage_" + outputFile.getName(),
                                    "TaxTriage analysis result: " + outputFile.getName(),
                                    "N", // Placeholder sequence for report
                                    new Date()
                            );
                            reportDoc.setDescription(content);
                            results.add(reportDoc);
                        } catch (IOException e) {
                            logger.log(Level.WARNING, "Failed to read result file: " + outputFile.getName(), e);
                        }
                    }
                }
            }

            logger.info("Converted " + results.size() + " workflow results to Geneious documents");
            return results;

        } catch (Exception e) {
            throw new DocumentOperationException("Failed to convert workflow results", e);
        }
    }

    /**
     * Generates a summary report of the workflow execution and results.
     *
     * @param workflowContext the workflow context
     * @param originalDocuments the original input documents
     * @return formatted summary report
     */
    private String generateSummaryReport(WorkflowContext workflowContext,
                                        AnnotatedPluginDocument[] originalDocuments) {
        StringBuilder report = new StringBuilder();
        report.append("=== TaxTriage Analysis Report ===\n\n");

        // Workflow information
        report.append("Workflow ID: ").append(workflowContext.getWorkflowId()).append("\n");
        report.append("Execution Status: ").append(workflowContext.getState()).append("\n");
        report.append("Start Time: ").append(workflowContext.getStartedAt()).append("\n");
        report.append("Completion Time: ").append(workflowContext.getCompletedAt()).append("\n");

        // Input information
        report.append("\nInput Documents:\n");
        for (AnnotatedPluginDocument doc : originalDocuments) {
            try {
                SequenceDocument seqDoc = (SequenceDocument) doc.getDocument();
                report.append("- ").append(doc.getName())
                      .append(" (").append(seqDoc.getSequenceLength()).append(" bp)\n");
            } catch (DocumentOperationException e) {
                report.append("- ").append(doc.getName())
                      .append(" (length unknown)\n");
            }
        }

        // Configuration
        TaxTriageOptions options = workflowContext.getOptions();
        if (options != null) {
            report.append("\nAnalysis Parameters:\n");
            report.append("- Sequencing Preset: ").append(options.getSequencingPreset().getDisplayName()).append("\n");
            report.append("- Quality Threshold: ").append(options.getQualityThreshold()).append("\n");
            report.append("- Min Read Length: ").append(options.getMinReadLength()).append("\n");
            report.append("- Thread Count: ").append(options.getThreadCount()).append("\n");
            report.append("- Memory Limit: ").append(options.getMemoryLimit()).append(" GB\n");
        }

        // Results
        report.append("\nResults:\n");
        if (workflowContext.getOutputFiles() != null) {
            report.append("- Output Files Generated: ").append(workflowContext.getOutputFiles().size()).append("\n");

            // List key result files
            if (workflowContext.getWorkflowReport() != null) {
                report.append("- Workflow Report: ").append(workflowContext.getWorkflowReport().getName()).append("\n");
            }
            if (workflowContext.getTaxonomyResults() != null) {
                report.append("- Taxonomy Results: ").append(workflowContext.getTaxonomyResults().getName()).append("\n");
            }
            if (workflowContext.getQualityReport() != null) {
                report.append("- Quality Report: ").append(workflowContext.getQualityReport().getName()).append("\n");
            }
        } else {
            report.append("- No output files were generated\n");
        }

        // Workflow summary
        report.append("\nWorkflow Summary:\n");
        report.append(workflowContext.getSummary());

        return report.toString();
    }

    /**
     * Checks if a file is a report file that should be imported as a document.
     *
     * @param file the file to check
     * @return true if it's a report file
     */
    private boolean isReportFile(File file) {
        if (file == null || !file.isFile()) {
            return false;
        }

        String name = file.getName().toLowerCase();
        return name.endsWith(".html") || name.endsWith(".txt") ||
               name.endsWith(".csv") || name.endsWith(".tsv") ||
               name.contains("report") || name.contains("summary") ||
               name.contains("classification") || name.contains("taxonomy");
    }

    /**
     * Gets the action options for this operation in the Geneious UI.
     *
     * @return action options configuration
     */
    @Override
    public GeneiousActionOptions getActionOptions() {
        return new GeneiousActionOptions("TaxTriage Analysis",
                "Perform taxonomic classification and pathogen detection using TaxTriage workflows")
                .setMainMenuLocation(GeneiousActionOptions.MainMenu.Tools)
                .setInMainToolbar(true)
                .setInPopupMenu(true)
                .setAvailableToWorkflows(true);
    }

    /**
     * Gets the help text for this operation.
     *
     * @return help text
     */
    @Override
    public String getHelp() {
        return "TaxTriage Analysis performs comprehensive taxonomic classification to identify " +
               "pathogens and characterize microbial communities from sequencing data.\n\n" +
               "Input Options:\n" +
               "• Select documents from Geneious (sequences, FASTQ, FASTA)\n" +
               "• Use file browser to select individual files from filesystem\n" +
               "• Use directory browser to select a folder containing input files\n\n" +
               "Supported File Types:\n" +
               "• FASTQ files (single-end or paired-end)\n" +
               "• FASTA sequences\n" +
               "• Multiple sequence alignments\n\n" +
               "Workflow Features:\n" +
               "• Quality filtering and preprocessing\n" +
               "• Taxonomic classification using multiple databases\n" +
               "• Pathogen detection and identification\n" +
               "• Comprehensive reporting with visualizations\n" +
               "• Support for ONT, Illumina PE, and Illumina SE presets\n\n" +
               "Output includes classification reports, quality metrics, and taxonomic summaries.";
    }

    /**
     * Defines which document types and combinations this operation can process.
     *
     * @return array of document selection signatures
     */
    @Override
    public DocumentSelectionSignature[] getSelectionSignatures() {
        // Accept any document type or no selection (0 to MAX)
        // This allows the plugin to be launched without requiring sequence selection
        DocumentSelectionSignature anySignature = new DocumentSelectionSignature(
            PluginDocument.class, 0, Integer.MAX_VALUE
        );

        return new DocumentSelectionSignature[]{anySignature};
    }

    /**
     * Creates and returns the options for this operation.
     *
     * @param documents documents that will be processed
     * @return configured options object
     */
    @Override
    public Options getOptions(AnnotatedPluginDocument... documents) {
        try {
            List<AnnotatedPluginDocument> docList = new ArrayList<AnnotatedPluginDocument>();
            if (documents != null) {
                for (AnnotatedPluginDocument doc : documents) {
                    docList.add(doc);
                }
            }
            return new TaxTriageOptions(docList);
        } catch (Exception e) {
            // Return basic options if there's an error creating advanced options
            return new Options(TaxTriageOptions.class);
        }
    }
}
