package com.jhuapl.taxtriage.geneious;

import com.biomatters.geneious.publicapi.plugin.Options;

public class TaxTriageOptions extends Options {

    public TaxTriageOptions() {
        addLabel("TaxTriage Settings");
        addDivider("Settings");
        addStringOption("preset", "Preset", "ONT");
    }
}