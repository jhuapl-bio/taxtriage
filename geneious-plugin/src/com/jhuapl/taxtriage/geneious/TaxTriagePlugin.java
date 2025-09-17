package com.jhuapl.taxtriage.geneious;

import com.biomatters.geneious.publicapi.plugin.GeneiousPlugin;
import com.biomatters.geneious.publicapi.plugin.DocumentOperation;

/**
 * Main plugin class for TaxTriage Geneious integration.
 *
 * This plugin provides taxonomic classification and analysis capabilities
 * within the Geneious bioinformatics platform, integrating TaxTriage
 * workflows with Geneious sequence data management.
 */
public class TaxTriagePlugin extends GeneiousPlugin {

    @Override
    public String getName() {
        return "TaxTriage";
    }

    @Override
    public String getDescription() {
        return "Taxonomic classification and analysis plugin for pathogen detection and identification";
    }

    @Override
    public String getHelp() {
        return "TaxTriage provides advanced taxonomic classification capabilities for " +
               "identifying pathogens and analyzing microbial communities from sequencing data. " +
               "The plugin integrates TaxTriage workflows directly into Geneious for seamless " +
               "analysis of FASTQ files and sequence assemblies.";
    }

    @Override
    public String getAuthors() {
        return "Johns Hopkins University Applied Physics Laboratory";
    }

    @Override
    public String getVersion() {
        return "1.0.0";
    }

    @Override
    public String getMinimumApiVersion() {
        return "4.1";
    }

    @Override
    public int getMaximumApiVersion() {
        return 4;
    }

    @Override
    public DocumentOperation[] getDocumentOperations() {
        return new DocumentOperation[] {
            new TaxTriageSimpleOperation(),
            new TaxTriageImportOperation()  // Test import operation
        };
    }
}