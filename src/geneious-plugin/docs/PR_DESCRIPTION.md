# Geneious Plugin for TaxTriage

## Summary

This is a Geneious Prime plugin that integrates TaxTriage taxonomic classification into the Geneious bioinformatics platform. This plugin enables researchers to perform metagenomic analysis directly from within Geneious using Docker-based workflows.

### Package Structure

```
com.jhuapl.taxtriage.geneious/
├── TaxTriagePlugin.java          # Plugin entry point
├── TaxTriageOperation.java       # Operation interface
├── TaxTriageSimpleOperation.java # Main implementation
├── TaxTriageOptions.java         # Configuration UI
├── config/                       # Configuration generation
│   ├── TaxTriageConfig.java
│   ├── ConfigGenerator.java
│   └── SampleSheetBuilder.java
├── database/                     # Database management
│   └── DatabaseManager.java
├── docker/                       # Docker integration
│   ├── DockerManager.java
│   ├── ExecutionMonitor.java
│   ├── NextflowTraceMonitor.java
│   └── VolumeMapper.java
├── execution/                    # Workflow execution
│   ├── TaxTriageWorkflowRunner.java
│   ├── WorkflowExecutor.java
│   └── WorkflowContext.java
├── importers/                    # File import system
│   ├── BamImportHelper.java
│   ├── SampleBasedImporter.java
│   └── NCBIReferenceDownloader.java
└── results/                      # Result processing
    └── TaxTriageResultImporter.java
```

## Installation

### Building the Plugin

```bash
cd geneious-plugin
ant clean build distribute
```

This creates `build/TaxTriage.gplugin` for installation.

### Installing in Geneious

Copy `TaxTriage.gplugin` to the Geneious plugins directory:

- **macOS:** `~/.geneious/plugins/`
- **Windows:** `%USERPROFILE%\.geneious\plugins\`
- **Linux:** `~/.geneious/plugins/`

Restart Geneious to load the plugin.

Note that this has only been tested on macOS.

## Requirements

### System Requirements

- Tested on Geneious Prime 2025.2.2
- Docker Desktop installed and running
- Minimum 8GB RAM (16GB+ recommended)
- Disk space for TaxTriage databases and Nextflow work directories. Maybe 50GB?

---

**Built by dave o'connor and claude code** ✅
