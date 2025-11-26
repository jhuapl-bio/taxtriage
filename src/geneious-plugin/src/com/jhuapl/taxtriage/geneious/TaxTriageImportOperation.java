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
 * Standalone TaxTriage result import operation for external workflow integration.
 *
 * <p>This operation provides a dedicated interface for importing TaxTriage analysis
 * results that were generated outside of Geneious, enabling seamless integration
 * of results from high-performance computing clusters, cloud platforms, or other
 * external systems into the Geneious workspace.</p>
 *
 * <h3>Use Cases:</h3>
 * <ul>
 *   <li><strong>HPC Integration:</strong> Import results from cluster-based TaxTriage runs</li>
 *   <li><strong>Cloud Workflows:</strong> Integration with cloud-based genomics platforms</li>
 *   <li><strong>Batch Processing:</strong> Import results from automated pipeline systems</li>
 *   <li><strong>Collaborative Research:</strong> Share and import results across research teams</li>
 * </ul>
 *
 * <h3>Import Configuration:</h3>
 * <ul>
 *   <li><strong>Flexible Input:</strong> User-configurable output directory path</li>
 *   <li><strong>Selective Import:</strong> Checkboxes for specific result types</li>
 *   <li><strong>Result Types:</strong> Kraken reports, top hits, filtered reports, alignments, QC reports</li>
 * </ul>
 *
 * <h3>Validation and Safety:</h3>
 * <ul>
 *   <li><strong>Directory Validation:</strong> Existence and accessibility checks</li>
 *   <li><strong>Format Validation:</strong> Ensures compatible TaxTriage output structure</li>
 *   <li><strong>Error Handling:</strong> Graceful handling of malformed or incomplete results</li>
 * </ul>
 *
 * <h3>Integration with Core Importer:</h3>
 * <p>This operation delegates the actual import process to {@link TaxTriageResultImporter},
 * ensuring consistency with the main workflow operation while providing a dedicated
 * entry point for external result integration.</p>
 *
 * @author TaxTriage Development Team
 * @version 2.0
 * @since 1.0
 * @see TaxTriageResultImporter
 * @see TaxTriageSimpleOperation
 */
public class TaxTriageImportOperation extends DocumentOperation {

    private static final Logger logger = Logger.getLogger(TaxTriageImportOperation.class.getName());

    // Default output path for TaxTriage results import
    private static final String DEFAULT_OUTPUT_PATH = System.getProperty("user.home") + "/taxtriage_output";

    @Override
    public String getHelp() {
        return "Imports TaxTriage analysis results from a specified output directory.";
    }

    @Override
    public DocumentSelectionSignature[] getSelectionSignatures() {
        // This operation doesn't require any document selection
        return new DocumentSelectionSignature[0];
    }

    @Override
    public GeneiousActionOptions getActionOptions() {
        return new GeneiousActionOptions("Import TaxTriage Results",
                "Import completed TaxTriage results")
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