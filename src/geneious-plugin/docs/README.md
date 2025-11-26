# TaxTriage Geneious Plugin

This directory contains the Geneious plugin implementation for TaxTriage, which provides taxonomic classification and analysis capabilities directly within the Geneious bioinformatics platform.

## Project Structure

```text
src/geneious-plugin/
├── src/                  # Java source code
│   └── com/
│       └── jhuapl/
│           └── taxtriage/
│               └── geneious/
├── test/                 # Unit tests
│   └── com/
│       └── jhuapl/
│           └── taxtriage/
│               └── geneious/
├── lib/                  # Required JAR dependencies
│   ├── GeneiousPublicAPI.jar  # Geneious plugin API
│   ├── jdom.jar               # XML processing
│   └── jebl.jar               # Bioinformatics utilities
├── resources/            # Plugin resources
│   ├── presets/          # Predefined analysis presets
│   └── templates/        # Configuration templates
├── plugin.properties     # Plugin metadata
├── build.xml             # Ant build script
└── README.md             # This file
```

## Dependencies

The plugin requires the following Geneious SDK components:

- **GeneiousPublicAPI.jar**: Core plugin API for Geneious integration
- **jdom.jar**: XML document processing
- **jebl.jar**: Java Evolutionary Biology Library

These dependencies are automatically copied from the Geneious devkit during project setup.

## Building the Plugin

### Prerequisites

- Java 11 or higher
- Apache Ant
- Geneious devkit (for API dependencies)

### Build Commands

Make sure you `cd src/geneious-plugin/` before running the following commands:

```bash
# Clean previous builds
ant clean

# Compile source code
ant compile

# Build the plugin (.gplugin file)
ant distribute

# Install plugin to local Geneious installation
ant install
```

### Build Targets

- **compile**: Compiles Java source code to classes
- **compile-tests**: Compiles unit tests
- **build**: Creates the plugin folder structure with JAR and resources
- **distribute**: Packages the plugin into a .gplugin file for distribution
- **install**: Copies the .gplugin file to the local Geneious plugins directory
- **clean**: Removes build artifacts

## Plugin Architecture

The TaxTriage Geneious plugin follows the standard Geneious plugin architecture:

1. **Main Plugin Class**: `com.jhuapl.taxtriage.geneious.TaxTriagePlugin`

   - Extends `GeneiousPlugin`
   - Defines plugin metadata and services

2. **Service Classes**: Implement specific analysis operations

   - Taxonomic classification services
   - Report generation services
   - Data import/export services

3. **UI Components**: Provide user interface elements
   - Analysis parameter dialogs
   - Result viewers
   - Configuration panels

## Development Guidelines

### Package Structure

- Use the base package: `com.jhuapl.taxtriage.geneious`
- Organize classes into logical subpackages:
  - `.services` - Analysis services
  - `.ui` - User interface components
  - `.util` - Utility classes
  - `.model` - Data models

### API Version Compatibility

- Target Geneious API version 4.1+
- Use compatible Java features (Java 11+)
- Follow Geneious plugin development best practices

### Testing

- Place unit tests in the `test/` directory
- Mirror the source package structure in tests
- Use standard JUnit testing framework

## Integration with TaxTriage

This plugin integrates with the main TaxTriage pipeline by:

1. Providing a Geneious-native interface for taxonomic analysis
2. Executing TaxTriage workflows on sequence data within Geneious
3. Displaying results in Geneious viewers and reports
4. Enabling export of results to TaxTriage formats

## Installation

1. Build the plugin using `ant distribute`
2. Copy the generated `.gplugin` file to your Geneious plugins directory:
   - macOS: `~/Library/Application Support/Geneious/plugins/`
   - Windows: `%APPDATA%/Geneious/plugins/`
   - Linux: `~/.geneious/plugins/`
3. Restart Geneious to load the plugin

## Development Status

This is the initial project structure setup. The following components need to be implemented:

- [ ] Main plugin class
- [ ] Taxonomic classification service
- [ ] UI components for parameter configuration
- [ ] Result viewers and reports
- [ ] Integration with TaxTriage core libraries
- [ ] Unit tests
- [ ] Documentation and help system
