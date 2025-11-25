#!/bin/bash

cd /tmp/test-plugin-nextflow-20250930_174457

export NXF_DEBUG=x
export NXF_ANSI_LOG=false

"/Users/dho/Library/Application Support/Geneious/plugins/com.jhuapl.taxtriage.geneious.TaxTriagePlugin/bin/nextflow" run https://github.com/jhuapl-bio/taxtriage -r main
