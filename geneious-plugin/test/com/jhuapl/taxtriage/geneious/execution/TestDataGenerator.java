package com.jhuapl.taxtriage.geneious.execution;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.sequence.SequenceDocument;
import com.jhuapl.taxtriage.geneious.TaxTriageOptions;
import com.jhuapl.taxtriage.geneious.config.TaxTriageConfig;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

/**
 * Utility class for generating test data for workflow execution tests.
 *
 * This class provides methods to create mock documents, test sequences,
 * and sample files needed for comprehensive testing of the TaxTriage workflow.
 */
public class TestDataGenerator {

    private static final Random random = new Random(12345); // Fixed seed for reproducible tests

    private static final String[] SAMPLE_NAMES = {
            "sample_001", "sample_002", "sample_003", "E_coli_strain",
            "human_sample", "environmental_sample", "mock_community"
    };

    private static final String[] DNA_BASES = {"A", "T", "C", "G"};

    /**
     * Creates a mock AnnotatedPluginDocument with a sequence document.
     *
     * @param name the document name
     * @param sequence the DNA sequence
     * @return mock AnnotatedPluginDocument
     */
    public static AnnotatedPluginDocument createMockDocument(String name, String sequence) {
        AnnotatedPluginDocument mockDoc = mock(AnnotatedPluginDocument.class);
        SequenceDocument mockSeqDoc = mock(SequenceDocument.class);

        when(mockDoc.getName()).thenReturn(name);
        when(mockDoc.getDocument()).thenReturn(mockSeqDoc);
        when(mockSeqDoc.getName()).thenReturn(name);
        when(mockSeqDoc.getSequenceString()).thenReturn(sequence);
        when(mockSeqDoc.getSequenceLength()).thenReturn(sequence.length());

        return mockDoc;
    }

    /**
     * Creates multiple mock documents with random sequences.
     *
     * @param count the number of documents to create
     * @param sequenceLength the length of each sequence
     * @return array of mock documents
     */
    public static AnnotatedPluginDocument[] createMockDocuments(int count, int sequenceLength) {
        AnnotatedPluginDocument[] documents = new AnnotatedPluginDocument[count];
        for (int i = 0; i < count; i++) {
            String name = SAMPLE_NAMES[i % SAMPLE_NAMES.length] + "_" + (i + 1);
            String sequence = generateRandomDNASequence(sequenceLength);
            documents[i] = createMockDocument(name, sequence);
        }
        return documents;
    }

    /**
     * Creates a mock TaxTriageOptions with typical settings.
     *
     * @param preset the sequencing preset to use
     * @return mock TaxTriageOptions
     */
    public static TaxTriageOptions createMockOptions(TaxTriageOptions.SequencingPreset preset) {
        TaxTriageOptions mockOptions = mock(TaxTriageOptions.class);

        when(mockOptions.getSequencingPreset()).thenReturn(preset);
        when(mockOptions.getQualityThreshold()).thenReturn(20);
        when(mockOptions.getMinReadLength()).thenReturn(50);
        when(mockOptions.getThreadCount()).thenReturn(4);
        when(mockOptions.getMemoryLimit()).thenReturn(8);
        when(mockOptions.isEnableSubsampling()).thenReturn(false);
        when(mockOptions.getSubsampleSize()).thenReturn(1000000);

        return mockOptions;
    }

    /**
     * Creates a mock TaxTriageConfig with default settings.
     *
     * @param preset the sequencing preset
     * @param outputDir the output directory
     * @return mock TaxTriageConfig
     */
    public static TaxTriageConfig createMockConfig(String preset, File outputDir) {
        TaxTriageConfig mockConfig = mock(TaxTriageConfig.class);

        when(mockConfig.getPreset()).thenReturn(preset);
        when(mockConfig.getQualityThreshold()).thenReturn(20);
        when(mockConfig.getMinReadLength()).thenReturn(50);
        when(mockConfig.getThreadCount()).thenReturn(4);
        when(mockConfig.getMemoryLimitGb()).thenReturn(8);
        when(mockConfig.getOutputDirectory()).thenReturn(outputDir);
        when(mockConfig.getKrakenDatabase()).thenReturn("/databases/kraken2");
        when(mockConfig.getBrackenDatabase()).thenReturn("/databases/bracken");
        when(mockConfig.getDockerProfile()).thenReturn("docker");
        when(mockConfig.isEnableKrona()).thenReturn(true);
        when(mockConfig.isEnableMultiQC()).thenReturn(true);
        when(mockConfig.getTopTaxa()).thenReturn(10);
        when(mockConfig.getConfidenceThreshold()).thenReturn(0.1);
        when(mockConfig.validate()).thenReturn(null); // Valid configuration

        return mockConfig;
    }

    /**
     * Generates a random DNA sequence of the specified length.
     *
     * @param length the length of the sequence
     * @return random DNA sequence
     */
    public static String generateRandomDNASequence(int length) {
        StringBuilder sequence = new StringBuilder();
        for (int i = 0; i < length; i++) {
            sequence.append(DNA_BASES[random.nextInt(DNA_BASES.length)]);
        }
        return sequence.toString();
    }

    /**
     * Creates sample FASTQ files in the specified directory.
     *
     * @param directory the directory to create files in
     * @param fileNames the names of the files to create
     * @param sequenceLength the length of sequences in each file
     * @throws IOException if file creation fails
     */
    public static void createSampleFastqFiles(Path directory, String[] fileNames, int sequenceLength)
            throws IOException {
        Files.createDirectories(directory);

        for (String fileName : fileNames) {
            Path filePath = directory.resolve(fileName);
            StringBuilder content = new StringBuilder();

            // Create 10 sequences per file
            for (int i = 0; i < 10; i++) {
                String seqName = fileName.replace(".fastq", "") + "_seq_" + (i + 1);
                String sequence = generateRandomDNASequence(sequenceLength);
                String quality = generateQualityString(sequenceLength);

                content.append("@").append(seqName).append("\n");
                content.append(sequence).append("\n");
                content.append("+").append("\n");
                content.append(quality).append("\n");
            }

            Files.writeString(filePath, content.toString());
        }
    }

    /**
     * Creates a sample samplesheet.csv file.
     *
     * @param filePath the path to create the samplesheet
     * @param sampleNames the sample names to include
     * @param preset the sequencing preset
     * @throws IOException if file creation fails
     */
    public static void createSampleSheet(Path filePath, String[] sampleNames, String preset) throws IOException {
        StringBuilder content = new StringBuilder();
        content.append("sample,fastq_1,fastq_2,long_fastq,fasta,single_end\n");

        for (String sampleName : sampleNames) {
            switch (preset.toUpperCase()) {
                case "ILLUMINA_PE":
                    content.append(sampleName).append(",")
                           .append(sampleName).append("_R1.fastq,")
                           .append(sampleName).append("_R2.fastq,")
                           .append(",,false\n");
                    break;
                case "ILLUMINA_SE":
                    content.append(sampleName).append(",")
                           .append(sampleName).append(".fastq,")
                           .append(",,,true\n");
                    break;
                case "ONT":
                    content.append(sampleName).append(",")
                           .append(",,")
                           .append(sampleName).append(".fastq,")
                           .append(",false\n");
                    break;
                default:
                    content.append(sampleName).append(",")
                           .append(sampleName).append(".fastq,")
                           .append(",,,true\n");
                    break;
            }
        }

        Files.writeString(filePath, content.toString());
    }

    /**
     * Creates a mock workflow report file.
     *
     * @param filePath the path to create the report
     * @param workflowId the workflow ID
     * @param successful whether the workflow was successful
     * @throws IOException if file creation fails
     */
    public static void createMockWorkflowReport(Path filePath, String workflowId, boolean successful)
            throws IOException {
        StringBuilder report = new StringBuilder();
        report.append("<!DOCTYPE html>\n<html>\n<head>\n");
        report.append("<title>TaxTriage Workflow Report</title>\n");
        report.append("</head>\n<body>\n");
        report.append("<h1>TaxTriage Workflow Report</h1>\n");
        report.append("<p>Workflow ID: ").append(workflowId).append("</p>\n");
        report.append("<p>Status: ").append(successful ? "Completed Successfully" : "Failed").append("</p>\n");

        if (successful) {
            report.append("<h2>Results Summary</h2>\n");
            report.append("<ul>\n");
            report.append("<li>Total reads processed: 50,000</li>\n");
            report.append("<li>Classified reads: 45,000 (90%)</li>\n");
            report.append("<li>Top genus: Escherichia (25%)</li>\n");
            report.append("<li>Top species: E. coli (20%)</li>\n");
            report.append("</ul>\n");

            report.append("<h2>Quality Metrics</h2>\n");
            report.append("<ul>\n");
            report.append("<li>Mean read quality: 32.5</li>\n");
            report.append("<li>Mean read length: 150 bp</li>\n");
            report.append("<li>Adapter contamination: 0.5%</li>\n");
            report.append("</ul>\n");
        } else {
            report.append("<h2>Error Information</h2>\n");
            report.append("<p>The workflow failed during execution. Please check the log files for details.</p>\n");
        }

        report.append("</body>\n</html>");
        Files.writeString(filePath, report.toString());
    }

    /**
     * Creates mock taxonomy results in CSV format.
     *
     * @param filePath the path to create the results file
     * @throws IOException if file creation fails
     */
    public static void createMockTaxonomyResults(Path filePath) throws IOException {
        StringBuilder csv = new StringBuilder();
        csv.append("taxonomy_id,scientific_name,common_name,rank,reads,percentage\n");
        csv.append("562,Escherichia coli,,species,10000,20.0\n");
        csv.append("511145,Escherichia coli str. K-12 substr. MG1655,,strain,8000,16.0\n");
        csv.append("1352,Enterococcus faecalis,,species,5000,10.0\n");
        csv.append("1280,Staphylococcus aureus,,species,3000,6.0\n");
        csv.append("287,Pseudomonas aeruginosa,,species,2500,5.0\n");
        csv.append("210,Helicobacter pylori,,species,2000,4.0\n");
        csv.append("1423,Bacillus subtilis,,species,1500,3.0\n");
        csv.append("9606,Homo sapiens,human,species,18000,36.0\n");

        Files.writeString(filePath, csv.toString());
    }

    /**
     * Creates a simple Nextflow configuration file.
     *
     * @param filePath the path to create the config file
     * @param preset the sequencing preset
     * @throws IOException if file creation fails
     */
    public static void createMockNextflowConfig(Path filePath, String preset) throws IOException {
        StringBuilder config = new StringBuilder();
        config.append("// TaxTriage Nextflow Configuration\n\n");
        config.append("params {\n");
        config.append("    preset = '").append(preset.toLowerCase()).append("'\n");
        config.append("    quality_threshold = 20\n");
        config.append("    min_read_length = 50\n");
        config.append("    enable_krona = true\n");
        config.append("    enable_multiqc = true\n");
        config.append("}\n\n");
        config.append("process {\n");
        config.append("    cpus = 4\n");
        config.append("    memory = '8.GB'\n");
        config.append("}\n\n");
        config.append("docker {\n");
        config.append("    enabled = true\n");
        config.append("    runOptions = '--rm'\n");
        config.append("}\n");

        Files.writeString(filePath, config.toString());
    }

    /**
     * Generates a quality string for FASTQ format.
     *
     * @param length the length of the quality string
     * @return quality string with random but realistic quality scores
     */
    private static String generateQualityString(int length) {
        StringBuilder quality = new StringBuilder();
        // Generate quality scores between 20-40 (ASCII 53-73)
        for (int i = 0; i < length; i++) {
            int qualityScore = 20 + random.nextInt(21); // 20-40
            char qualityChar = (char) (33 + qualityScore); // Convert to ASCII
            quality.append(qualityChar);
        }
        return quality.toString();
    }

    /**
     * Creates a realistic test workflow context.
     *
     * @param tempDir the temporary directory for the workflow
     * @param documentCount the number of input documents
     * @param preset the sequencing preset
     * @return configured WorkflowContext
     * @throws IOException if setup fails
     */
    public static WorkflowContext createRealisticWorkflowContext(Path tempDir, int documentCount,
                                                                TaxTriageOptions.SequencingPreset preset)
            throws IOException {
        // Create mock documents
        AnnotatedPluginDocument[] documents = createMockDocuments(documentCount, 150);
        List<AnnotatedPluginDocument> docList = List.of(documents);

        // Create mock options and config
        TaxTriageOptions options = createMockOptions(preset);
        TaxTriageConfig config = createMockConfig(preset.name(), tempDir.resolve("output").toFile());

        // Create workflow context
        WorkflowContext context = new WorkflowContext(docList, options, config, tempDir);

        // Create directory structure
        Files.createDirectories(context.getInputDirectory());
        Files.createDirectories(context.getOutputDirectory());
        Files.createDirectories(context.getConfigDirectory());

        // Create sample files
        String[] fileNames = new String[documentCount];
        for (int i = 0; i < documentCount; i++) {
            fileNames[i] = SAMPLE_NAMES[i % SAMPLE_NAMES.length] + "_" + (i + 1) + ".fastq";
        }

        createSampleFastqFiles(context.getInputDirectory(), fileNames, 150);

        // Create configuration files
        createSampleSheet(context.getConfigDirectory().resolve("samplesheet.csv"),
                         fileNames, preset.name());
        createMockNextflowConfig(context.getConfigDirectory().resolve("nextflow.config"),
                               preset.name());

        // Set files in context
        context.setSampleSheetFile(context.getConfigDirectory().resolve("samplesheet.csv").toFile());
        context.setConfigFile(context.getConfigDirectory().resolve("nextflow.config").toFile());
        context.setParamsFile(context.getConfigDirectory().resolve("params.json").toFile());

        List<File> exportedFiles = new ArrayList<>();
        for (String fileName : fileNames) {
            exportedFiles.add(context.getInputDirectory().resolve(fileName).toFile());
        }
        context.setExportedFiles(exportedFiles);

        return context;
    }
}