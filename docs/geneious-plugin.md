# Geneious Plugin

Dave O'Connor's Laboratory at the University of Wisconsin-Madison developed a custom **Geneious Prime plugin** that allows TaxTriage analyses to be run directly from within the Geneious bioinformatics platform — no command-line knowledge required.

The plugin integrates with Docker and Nextflow to execute the full TaxTriage workflow, accepting file inputs through Geneious's GUI and importing results back into the project.

---

## Prerequisites

### Required Software

- **Geneious Prime** version 2024.0.2 or later
- **Docker** — must be installed and running on your system
- **Nextflow** — pulled automatically via Docker on first run
- **TaxTriage databases** — downloaded automatically by the workflow on first run

---

## Installation

1. Locate the `TaxTriage.gplugin` file in the repository:

   ```
   src/geneious-plugin/TaxTriage.gplugin
   ```

   > This file is built for macOS ARM (Apple Silicon). It may work on Windows/Linux but has not been officially tested on those platforms.

2. Copy it to your Geneious plugins directory:

   | OS      | Path                               |
   | ------- | ---------------------------------- |
   | macOS   | `~/.geneious/plugins/`             |
   | Windows | `%USERPROFILE%\.geneious\plugins\` |
   | Linux   | `~/.geneious/plugins/`             |

3. Restart Geneious.

---

## Usage

### Launch TaxTriage

Go to **Tools > TaxTriage Analysis** in Geneious.

### Select Input Files

Choose one of:

- **Input Files** — Click "Browse" to select individual FASTQ/FASTA files
- **Input Directory** — Click "Browse" to select a folder of sequence files

### Configure Sequencing Preset

| Preset                | Use For                   |
| --------------------- | ------------------------- |
| ONT (Oxford Nanopore) | Long-read sequencing data |
| Illumina PE           | Paired-end Illumina reads |
| Illumina SE           | Single-end Illumina reads |

### Quality Parameters

Set minimum quality thresholds for filtering (equivalent to `--minq`).

### Run the Analysis

Click **Run** to start the pipeline. The plugin generates a configuration file, launches Docker + Nextflow, and runs the TaxTriage workflow.

Results are imported back into Geneious upon completion.

---

## Building the Plugin from Source

### Prerequisites

- Java 11 or higher
- Apache Ant
- Geneious devkit (for API dependencies)

### Project Structure

```
src/geneious-plugin/
├── src/                  # Java source code
│   └── com/jhuapl/taxtriage/geneious/
├── test/                 # Unit tests
├── lib/                  # Required JARs
│   ├── GeneiousPublicAPI.jar
│   ├── jdom.jar
│   └── jebl.jar
├── resources/
│   ├── presets/          # Analysis presets
│   └── templates/        # Configuration templates
├── plugin.properties     # Plugin metadata
├── build.xml             # Ant build script
└── README.md
```

### Dependencies

| Dependency              | Purpose                           |
| ----------------------- | --------------------------------- |
| `GeneiousPublicAPI.jar` | Core Geneious plugin API          |
| `jdom.jar`              | XML document processing           |
| `jebl.jar`              | Java Evolutionary Biology Library |

These JARs are automatically copied from the Geneious devkit during project setup.

### Build Commands

From the `src/geneious-plugin/` directory:

```bash
# Set up the project (copies devkit JARs)
ant setup

# Build the plugin
ant build

# Run unit tests
ant test

# Package the .gplugin file
ant package
```

The resulting `TaxTriage.gplugin` will be in `build/`.

---

## Docker Integration

The plugin handles Docker automatically. On first run:

1. The Nextflow Docker image is pulled
2. TaxTriage container images are pulled
3. Databases are downloaded to the working directory

For details on the Docker integration architecture, see [`src/geneious-plugin/docs/DOCKER_INTEGRATION.md`](https://github.com/jhuapl-bio/taxtriage/blob/main/src/geneious-plugin/docs/DOCKER_INTEGRATION.md) in the repository.

---

## Support and Compatibility

- **Tested on:** macOS ARM (Apple Silicon) with Geneious Prime 2024.0.2+
- **Untested on:** Windows, Linux (may work but not officially supported)
- For issues with the plugin, open a GitHub issue or contact the development team.
