# TaxTriage Geneious Plugin - Usage Instructions

## Overview

The TaxTriage Geneious plugin now includes full Docker/Nextflow workflow execution capabilities. The plugin can:

1. **Accept file inputs** via file browsers (individual files or directories)
2. **Execute TaxTriage workflows** using Docker and Nextflow
3. **Generate configuration files** automatically based on user parameters
4. **Import results** back to Geneious (basic implementation)

## Prerequisites

### Required Software
- **Docker** must be installed and running on your system
- **Nextflow** Docker image will be automatically pulled when first run
- **TaxTriage databases** will be downloaded automatically by the workflow

### Installation
1. Copy `TaxTriage.gplugin` to your Geneious plugins directory:
   - macOS: `~/.geneious/plugins/`
   - Windows: `%USERPROFILE%\.geneious\plugins\`
   - Linux: `~/.geneious/plugins/`
2. Restart Geneious

## Usage

### 1. Launch TaxTriage
- Go to **Tools > TaxTriage Analysis** in Geneious
- The plugin will open with configuration options

### 2. Select Input Files
Choose one of these input methods:
- **Input Files**: Click "Browse" to select individual FASTQ/FASTA files
- **Input Directory**: Click "Browse" to select a directory containing sequence files

### 3. Configure Analysis Parameters

#### Sequencing Preset
- **ONT (Oxford Nanopore)**: For long-read sequencing data
- **Illumina PE**: For paired-end Illumina sequencing
- **Illumina SE**: For single-end Illumina sequencing

#### Quality Parameters
- **Quality Threshold**: Minimum Phred quality score (default: 20)
- **Min Read Length**: Minimum read length to retain (default: 50)

#### Database Configuration
- **Kraken Database**: Taxonomic classification database
- **Bracken Database**: Species-level abundance estimation database

#### Resource Settings
- **Thread Count**: Number of CPU cores to use
- **Memory Limit**: Maximum memory allocation in GB

### 4. Output Configuration
- **Output Directory**: Where TaxTriage results will be saved
- Default: `~/TaxTriage_Results`

### 5. Execute Workflow
- Click **OK** to start the analysis
- Monitor progress in the Geneious progress dialog
- The workflow will:
  1. Validate Docker availability
  2. Create a temporary workspace
  3. Generate Nextflow configuration files
  4. Execute the TaxTriage workflow via Docker
  5. Import compatible results back to Geneious

## Workflow Details

### Generated Files
The plugin creates these files in a temporary workspace:
- `config/nextflow.config`: Nextflow execution configuration
- `config/params.json`: Workflow parameters
- `config/samplesheet.csv`: Input sample specification
- `input/`: Copies of your input files
- `work/`: Nextflow working directory
- `output/`: TaxTriage analysis results

### Docker Execution
The plugin uses the `nextflow/nextflow:latest` Docker image to execute:
```bash
nextflow run https://github.com/jhuapl-bio/taxtriage -r main -profile docker \
  --input /data/work/config/samplesheet.csv \
  --outdir /data/output \
  --workdir /data/work \
  [additional parameters based on your configuration]
```

### Volume Mounts
Docker containers use these volume mappings:
- Input files: `workspace/input -> /data/input`
- Output directory: `workspace/output -> /data/output`
- Working directory: `workspace -> /data/work`

## Troubleshooting

### Docker Issues
- **"Docker not available"**: Ensure Docker is installed and running
- **Permission errors**: Check file permissions on input/output directories
- **Out of disk space**: Ensure sufficient disk space for databases and results

### Workflow Failures
- Check the Geneious log for detailed error messages
- Verify input files are valid FASTQ/FASTA format
- Ensure database selections are compatible with your data type

### Performance
- Long-read data (ONT) requires more memory and time
- First run will download databases (several GB)
- Subsequent runs use cached databases for faster execution

## Advanced Configuration

### Custom Parameters
For advanced users, you can modify the generated configuration files in the temporary workspace before execution.

### Database Paths
If using custom databases, select "Custom" in the database dropdowns and ensure your Docker containers have access to the database paths.

## Output Interpretation

TaxTriage generates several types of output:
- **Kraken2 reports**: Taxonomic classification results
- **Bracken reports**: Species-level abundance estimates
- **Krona plots**: Interactive visualization of results
- **MultiQC report**: Quality control summary

## Support

For issues specific to:
- **Plugin functionality**: Check Geneious logs and plugin documentation
- **TaxTriage workflow**: Refer to https://github.com/jhuapl-bio/taxtriage
- **Docker/Nextflow**: Check respective documentation

## Implementation Notes

This implementation provides:
- ✅ Complete workflow execution pipeline
- ✅ Docker integration with proper volume mounting
- ✅ Configuration file generation
- ✅ Progress monitoring and error handling
- ✅ File validation and workspace management
- ⚠️ Basic result import (full import requires database service integration)

The plugin creates a complete, working TaxTriage execution environment that integrates seamlessly with Geneious while leveraging the full power of the TaxTriage Nextflow workflow.