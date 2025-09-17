/**
 * Docker execution management for TaxTriage Geneious plugin.
 *
 * This package provides comprehensive Docker integration for executing TaxTriage
 * workflows within Geneious. It handles container lifecycle management, volume
 * mounting, progress monitoring, and cross-platform compatibility.
 *
 * <h2>Core Components</h2>
 * <ul>
 *   <li>{@link com.jhuapl.taxtriage.geneious.docker.DockerManager} - Main Docker operations manager</li>
 *   <li>{@link com.jhuapl.taxtriage.geneious.docker.VolumeMapper} - Cross-platform path mapping</li>
 *   <li>{@link com.jhuapl.taxtriage.geneious.docker.ExecutionMonitor} - Nextflow progress tracking</li>
 *   <li>{@link com.jhuapl.taxtriage.geneious.docker.ExecutionResult} - Command execution results</li>
 *   <li>{@link com.jhuapl.taxtriage.geneious.docker.DockerException} - Docker-specific exceptions</li>
 * </ul>
 *
 * <h2>Features</h2>
 * <ul>
 *   <li>Automatic Docker availability detection</li>
 *   <li>Container image management and pulling</li>
 *   <li>Cross-platform volume mounting (Windows, macOS, Linux)</li>
 *   <li>Real-time Nextflow progress monitoring</li>
 *   <li>Proper error handling and cancellation support</li>
 *   <li>Path escaping for special characters and spaces</li>
 *   <li>Comprehensive test coverage with mock implementations</li>
 * </ul>
 *
 * <h2>Usage Example</h2>
 * <pre>{@code
 * // Create Docker manager
 * DockerManager dockerManager = new DockerManager();
 *
 * // Execute TaxTriage workflow
 * ExecutionResult result = dockerManager.executeNextflowCommand(
 *     "nextflow run main.nf --input /input --output /output",
 *     inputDirectory,
 *     outputDirectory,
 *     workDirectory,
 *     progressListener
 * );
 *
 * if (result.isSuccess()) {
 *     // Process results
 *     String output = result.getOutput();
 * } else {
 *     // Handle errors
 *     String error = result.getErrorOutput();
 * }
 * }</pre>
 *
 * <h2>Testing</h2>
 * The package includes comprehensive test suites:
 * <ul>
 *   <li>Unit tests with JUnit 5 and Mockito</li>
 *   <li>Mock implementations for testing without Docker</li>
 *   <li>Simple assertion-based tests for basic functionality</li>
 *   <li>Cross-platform compatibility tests</li>
 * </ul>
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 * @since 1.0.0
 */
package com.jhuapl.taxtriage.geneious.docker;