package com.jhuapl.taxtriage.geneious;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.plugin.*;
import com.biomatters.geneious.publicapi.utilities.IconUtilities;
import com.jhuapl.taxtriage.geneious.results.TaxTriageResultImporter;
import jebl.util.ProgressListener;

import javax.swing.*;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * Test operation to import TaxTriage results from a specific output directory.
 * This is a development tool to test importing results without running the full workflow.
 */
public class TaxTriageImportOperation extends DocumentOperation {

    private static final Logger logger = Logger.getLogger(TaxTriageImportOperation.class.getName());

    // Test path - you can modify this to point to your output directory
    private static final String DEFAULT_OUTPUT_PATH = "/private/var/folders/pw/r9x37tvs2q5d6c044x2fhp200000gn/T/taxtriage_20250916_202306/output";

    @Override
    public String getHelp() {
        return "Imports TaxTriage analysis results from a specified output directory (for testing).";
    }

    @Override
    public DocumentSelectionSignature[] getSelectionSignatures() {
        // This operation doesn't require any document selection
        return new DocumentSelectionSignature[0];
    }

    @Override
    public GeneiousActionOptions getActionOptions() {
        return new GeneiousActionOptions("Import TaxTriage Results (Test)",
                "Import completed TaxTriage results for testing")
                .setMainMenuLocation(GeneiousActionOptions.MainMenu.Tools)
                .setInMainToolbar(false);
    }

    @Override
    public Options getOptions(AnnotatedPluginDocument... documents) throws DocumentOperationException {
        Options options = new Options(getClass());

        options.addLabel("Import TaxTriage Results from Output Directory");

        // Add a string option for the output directory path
        options.addStringOption("outputPath", "Output Directory:", DEFAULT_OUTPUT_PATH);

        options.addLabel("<html><b>Select Result Types to Import:</b></html>");

        // Add checkboxes for different result types
        options.addBooleanOption("importReports", "Import Kraken Reports", true);
        options.addBooleanOption("importTopHits", "Import Top Hits", true);
        options.addBooleanOption("importFilteredReports", "Import Filtered Reports", true);
        options.addBooleanOption("importAlignments", "Import Alignment Results", true);
        options.addBooleanOption("importQualityReports", "Import QC Reports", true);

        return options;
    }

    @Override
    public List<AnnotatedPluginDocument> performOperation(AnnotatedPluginDocument[] documents,
                                                         ProgressListener progressListener,
                                                         Options options) throws DocumentOperationException {

        List<AnnotatedPluginDocument> results = new ArrayList<>();

        String outputPath = (String) options.getValue("outputPath");
        if (outputPath == null || outputPath.trim().isEmpty()) {
            throw new DocumentOperationException("Please specify an output directory path");
        }

        Path outputDir = Paths.get(outputPath);
        if (!Files.exists(outputDir)) {
            throw new DocumentOperationException("Output directory does not exist: " + outputPath);
        }

        if (progressListener != null) {
            progressListener.setMessage("Importing TaxTriage results from: " + outputPath);
            progressListener.setProgress(0.1);
        }

        try {
            // Use the TaxTriageResultImporter
            TaxTriageResultImporter importer = new TaxTriageResultImporter();
            results = importer.importResults(outputDir, progressListener);

        } catch (Exception e) {
            logger.log(Level.SEVERE, "Error importing TaxTriage results", e);
            throw new DocumentOperationException("Error importing results: " + e.getMessage(), e);
        }

        if (progressListener != null) {
            progressListener.setMessage("Successfully imported " + results.size() + " documents");
            progressListener.setProgress(1.0);
        }

        logger.info("Imported " + results.size() + " TaxTriage result documents");
        return results;
    }
}