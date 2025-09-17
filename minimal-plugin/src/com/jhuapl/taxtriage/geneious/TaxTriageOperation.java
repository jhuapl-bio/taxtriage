package com.jhuapl.taxtriage.geneious;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.plugin.DocumentOperation;
import com.biomatters.geneious.publicapi.plugin.DocumentOperationException;
import com.biomatters.geneious.publicapi.plugin.DocumentSelectionSignature;
import com.biomatters.geneious.publicapi.plugin.GeneiousActionOptions;
import com.biomatters.geneious.publicapi.plugin.Options;
import com.biomatters.geneious.publicapi.documents.sequence.SequenceDocument;
import jebl.util.ProgressListener;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class TaxTriageOperation extends DocumentOperation {

    @Override
    public GeneiousActionOptions getActionOptions() {
        return new GeneiousActionOptions("TaxTriage",
                "Run TaxTriage taxonomic classification")
                .setMainMenuLocation(GeneiousActionOptions.MainMenu.Tools)
                .setInMainToolbar(true);
    }

    @Override
    public String getHelp() {
        return "TaxTriage performs taxonomic classification on sequence data.";
    }

    @Override
    public DocumentSelectionSignature[] getSelectionSignatures() {
        return new DocumentSelectionSignature[] {
            new DocumentSelectionSignature(SequenceDocument.class, 1, Integer.MAX_VALUE)
        };
    }

    @Override
    public Options getOptions(AnnotatedPluginDocument... documents) throws DocumentOperationException {
        return new TaxTriageOptions();
    }

    @Override
    public List<AnnotatedPluginDocument> performOperation(AnnotatedPluginDocument[] documents,
                                                          ProgressListener progressListener,
                                                          Options options) throws DocumentOperationException {

        if (progressListener != null) {
            progressListener.setMessage("TaxTriage Analysis Started");
        }

        // Show a dialog to confirm the plugin is working
        JOptionPane.showMessageDialog(null,
            "TaxTriage plugin is installed!\n\n" +
            "Selected " + documents.length + " sequence(s) for analysis.\n\n" +
            "Full Docker integration is pending.\n" +
            "This will run the TaxTriage workflow on your sequences.",
            "TaxTriage",
            JOptionPane.INFORMATION_MESSAGE);

        return new ArrayList<>();
    }
}