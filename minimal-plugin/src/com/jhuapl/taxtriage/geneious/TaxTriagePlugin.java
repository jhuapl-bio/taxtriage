package com.jhuapl.taxtriage.geneious;

import com.biomatters.geneious.publicapi.plugin.GeneiousPlugin;
import com.biomatters.geneious.publicapi.plugin.DocumentOperation;

public class TaxTriagePlugin extends GeneiousPlugin {

    @Override
    public String getName() {
        return "TaxTriage";
    }

    @Override
    public String getDescription() {
        return "Taxonomic classification using TaxTriage";
    }

    @Override
    public String getHelp() {
        return "TaxTriage provides taxonomic classification for sequence data.";
    }

    @Override
    public String getAuthors() {
        return "JHU APL";
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
            new TaxTriageOperation()
        };
    }
}