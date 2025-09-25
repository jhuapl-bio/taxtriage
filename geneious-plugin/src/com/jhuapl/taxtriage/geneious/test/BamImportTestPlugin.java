package com.jhuapl.taxtriage.geneious.test;

import com.biomatters.geneious.publicapi.databaseservice.DatabaseServiceException;
import com.biomatters.geneious.publicapi.databaseservice.Query;
import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.documents.sequence.SequenceDocument;
import com.biomatters.geneious.publicapi.plugin.*;
import com.jhuapl.taxtriage.geneious.importers.BamReferenceExtractor;
import com.jhuapl.taxtriage.geneious.importers.GenBankFileFixer;
import jebl.util.ProgressListener;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

/**
 * Test plugin that tries multiple approaches to import BAM files with automatic reference linking.
 * This will cascade through different strategies until one succeeds without showing a dialog.
 */
public class BamImportTestPlugin extends DocumentOperation {

    private static final Logger logger = Logger.getLogger(BamImportTestPlugin.class.getName());

    // Static cache for references across import attempts
    private static final Map<String, AnnotatedPluginDocument> REFERENCE_CACHE = new ConcurrentHashMap<>();

    @Override
    public String getHelp() {
        return "Tests multiple BAM import strategies to find one that works without dialogs";
    }

    @Override
    public DocumentSelectionSignature[] getSelectionSignatures() {
        return new DocumentSelectionSignature[0];
    }

    @Override
    public GeneiousActionOptions getActionOptions() {
        return new GeneiousActionOptions("Test BAM Import Strategies");
    }

    /**
     * Main test method that tries all import strategies
     */
    public static List<AnnotatedPluginDocument> testAllImportStrategies(
            Path minimap2Dir,
            WritableDatabaseService targetFolder,
            ProgressListener progressListener) {

        System.out.println("\n========================================");
        System.out.println("BAM IMPORT STRATEGY TESTER");
        System.out.println("========================================\n");

        // Find BAM and GenBank files
        List<File> bamFiles = new ArrayList<>();
        List<File> gbFiles = new ArrayList<>();

        try {
            Files.list(minimap2Dir)
                .forEach(path -> {
                    String name = path.getFileName().toString().toLowerCase();
                    if (name.endsWith(".bam")) {
                        bamFiles.add(path.toFile());
                    } else if (name.endsWith(".gb") || name.endsWith(".gbk")) {
                        gbFiles.add(path.toFile());
                    }
                });
        } catch (IOException e) {
            System.err.println("Error scanning directory: " + e.getMessage());
            return new ArrayList<>();
        }

        System.out.println("Found " + bamFiles.size() + " BAM file(s)");
        System.out.println("Found " + gbFiles.size() + " GenBank reference file(s)");

        if (bamFiles.isEmpty()) {
            System.err.println("No BAM files found in " + minimap2Dir);
            return new ArrayList<>();
        }

        File bamFile = bamFiles.get(0); // Test with first BAM file
        System.out.println("\nTesting with BAM file: " + bamFile.getName());

        // List of strategies to try
        List<ImportStrategy> strategies = Arrays.asList(
            new Strategy1_DirectToSameFolder(),
            new Strategy2_BatchImportTogether(),
            new Strategy3_PrewarmDatabase(),
            new Strategy4_CoLocationTemporary(),
            new Strategy5_ImportWithOptions(),
            new Strategy6_TransactionWrapped(),
            new Strategy7_BackgroundThreadTiming(),
            new Strategy8_SystemProperties(),
            new Strategy9_ForceIndexing(),
            new Strategy10_ReflectionContext(),
            new Strategy11_CustomProgressListener(),
            new Strategy12_SymlinkApproach(),
            new Strategy13_MultiStepWithCache(),
            new Strategy14_DirectMemoryImport(),
            new Strategy15_CombinedFolderImport()
        );

        // Try each strategy
        for (ImportStrategy strategy : strategies) {
            System.out.println("\n----------------------------------------");
            System.out.println("STRATEGY: " + strategy.getName());
            System.out.println("----------------------------------------");

            try {
                List<AnnotatedPluginDocument> result = strategy.importBam(
                    bamFile, gbFiles, minimap2Dir, targetFolder, progressListener);

                if (result != null && !result.isEmpty()) {
                    System.out.println("✓ SUCCESS: Imported " + result.size() + " document(s)");
                    System.out.println("Using strategy: " + strategy.getName());
                    return result;
                } else {
                    System.out.println("✗ Failed: No documents imported");
                }
            } catch (Exception e) {
                System.out.println("✗ Error: " + e.getMessage());
                e.printStackTrace();
            }

            // Small delay between attempts
            try {
                Thread.sleep(1000);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }
        }

        System.out.println("\n✗✗✗ All strategies failed ✗✗✗");
        return new ArrayList<>();
    }

    // ========== IMPORT STRATEGIES ==========

    interface ImportStrategy {
        String getName();
        List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                 Path minimap2Dir, WritableDatabaseService targetFolder,
                                                 ProgressListener progressListener) throws Exception;
    }

    /**
     * Strategy 1: Import references and BAM to same folder
     */
    static class Strategy1_DirectToSameFolder implements ImportStrategy {
        public String getName() { return "Direct Import to Same Folder"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            List<AnnotatedPluginDocument> allDocs = new ArrayList<>();

            // Import references first
            System.out.println("Importing " + gbFiles.size() + " reference files...");
            for (File gbFile : gbFiles) {
                GenBankFileFixer.fixGenBankFile(gbFile);
                List<AnnotatedPluginDocument> refs = PluginUtilities.importDocumentsToDatabase(
                    gbFile, targetFolder, progressListener);
                allDocs.addAll(refs);
            }

            // Wait for indexing
            Thread.sleep(2000);

            // Import BAM to same folder
            System.out.println("Importing BAM to same folder as references...");
            List<AnnotatedPluginDocument> bamDocs = PluginUtilities.importDocumentsToDatabase(
                bamFile, targetFolder, progressListener);
            allDocs.addAll(bamDocs);

            return bamDocs;
        }
    }

    /**
     * Strategy 2: Batch import all files together
     */
    static class Strategy2_BatchImportTogether implements ImportStrategy {
        public String getName() { return "Batch Import All Files Together"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            // Fix GenBank files first
            for (File gbFile : gbFiles) {
                GenBankFileFixer.fixGenBankFile(gbFile);
            }

            // Create array of all files
            List<File> allFiles = new ArrayList<>();
            allFiles.addAll(gbFiles);
            allFiles.add(bamFile);

            System.out.println("Batch importing " + allFiles.size() + " files together...");

            // Import all files
            List<AnnotatedPluginDocument> docs = new ArrayList<>();
            for (File file : allFiles) {
                docs.addAll(PluginUtilities.importDocuments(file, progressListener));
            }

            // Add to target folder
            List<AnnotatedPluginDocument> addedDocs = new ArrayList<>();
            for (AnnotatedPluginDocument doc : docs) {
                AnnotatedPluginDocument added = targetFolder.addDocumentCopy(doc, ProgressListener.EMPTY);
                if (added != null) {
                    addedDocs.add(added);
                }
            }

            return addedDocs;
        }
    }

    /**
     * Strategy 3: Pre-warm database with reference queries
     */
    static class Strategy3_PrewarmDatabase implements ImportStrategy {
        public String getName() { return "Database Pre-warming with Queries"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            // Import references
            System.out.println("Importing references...");
            List<AnnotatedPluginDocument> refs = new ArrayList<>();
            for (File gbFile : gbFiles) {
                GenBankFileFixer.fixGenBankFile(gbFile);
                refs.addAll(PluginUtilities.importDocumentsToDatabase(gbFile, targetFolder, progressListener));
            }

            // Pre-warm database cache
            System.out.println("Pre-warming database cache...");
            for (AnnotatedPluginDocument ref : refs) {
                try {
                    Query query = Query.Factory.createQuery("Name=" + ref.getName());
                    targetFolder.retrieve(query, ProgressListener.EMPTY);
                } catch (Exception e) {
                    // Ignore query errors
                }
            }

            Thread.sleep(1000);

            // Import BAM
            System.out.println("Importing BAM with pre-warmed cache...");
            return PluginUtilities.importDocumentsToDatabase(bamFile, targetFolder, progressListener);
        }
    }

    /**
     * Strategy 4: Co-location with temporary copy
     */
    static class Strategy4_CoLocationTemporary implements ImportStrategy {
        public String getName() { return "Temporary Co-location"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            // Create temp directory
            Path tempDir = Files.createTempDirectory("bam_import_");

            try {
                // Copy all files to temp directory
                System.out.println("Creating co-located temporary files...");
                for (File gbFile : gbFiles) {
                    GenBankFileFixer.fixGenBankFile(gbFile);
                    Files.copy(gbFile.toPath(), tempDir.resolve(gbFile.getName()), StandardCopyOption.REPLACE_EXISTING);
                }
                Path tempBam = tempDir.resolve(bamFile.getName());
                Files.copy(bamFile.toPath(), tempBam, StandardCopyOption.REPLACE_EXISTING);

                // Import from temp directory
                System.out.println("Importing from co-located temp directory...");
                List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(
                    tempDir.toFile(), progressListener);

                // Add to database
                List<AnnotatedPluginDocument> addedDocs = new ArrayList<>();
                for (AnnotatedPluginDocument doc : docs) {
                    addedDocs.add(targetFolder.addDocumentCopy(doc, ProgressListener.EMPTY));
                }

                return addedDocs;

            } finally {
                // Clean up temp files
                Files.list(tempDir).forEach(p -> {
                    try { Files.delete(p); } catch (Exception e) { }
                });
                Files.delete(tempDir);
            }
        }
    }

    /**
     * Strategy 5: Import with extensive options
     */
    static class Strategy5_ImportWithOptions implements ImportStrategy {
        public String getName() { return "Import with Reference Options"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            // Import references first
            System.out.println("Importing references...");
            List<AnnotatedPluginDocument> refs = new ArrayList<>();
            StringBuilder refUrns = new StringBuilder();
            StringBuilder refNames = new StringBuilder();

            for (File gbFile : gbFiles) {
                GenBankFileFixer.fixGenBankFile(gbFile);
                List<AnnotatedPluginDocument> imported = PluginUtilities.importDocumentsToDatabase(
                    gbFile, targetFolder, progressListener);
                refs.addAll(imported);

                for (AnnotatedPluginDocument ref : imported) {
                    if (refUrns.length() > 0) refUrns.append(",");
                    refUrns.append(ref.getURN());
                    if (refNames.length() > 0) refNames.append(",");
                    refNames.append(ref.getName());
                }
            }

            Thread.sleep(2000);

            // Build extensive options
            Map<String, String> options = new HashMap<>();
            options.put("referenceSequences", refUrns.toString());
            options.put("referenceURNs", refUrns.toString());
            options.put("referenceNames", refNames.toString());
            options.put("findReferences", "true");
            options.put("autoFindReferences", "true");
            options.put("suppressDialogs", "true");
            options.put("referenceDatabase", targetFolder.getName());
            options.put("targetDatabase", targetFolder.getUniqueID());

            System.out.println("Importing BAM with " + options.size() + " option parameters...");
            List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(bamFile, options, progressListener);

            // Add to database
            List<AnnotatedPluginDocument> addedDocs = new ArrayList<>();
            for (AnnotatedPluginDocument doc : docs) {
                addedDocs.add(targetFolder.addDocumentCopy(doc, ProgressListener.EMPTY));
            }

            return addedDocs;
        }
    }

    /**
     * Strategy 6: Transaction-wrapped import
     */
    static class Strategy6_TransactionWrapped implements ImportStrategy {
        public String getName() { return "Transaction-wrapped Import"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            System.out.println("Starting database transaction...");

            // Try to use transaction methods via reflection
            try {
                Method beginTx = targetFolder.getClass().getMethod("beginTransaction");
                Method commitTx = targetFolder.getClass().getMethod("commitTransaction");
                Method rollbackTx = targetFolder.getClass().getMethod("rollbackTransaction");

                beginTx.invoke(targetFolder);

                try {
                    // Import in transaction
                    List<AnnotatedPluginDocument> allDocs = new ArrayList<>();

                    for (File gbFile : gbFiles) {
                        GenBankFileFixer.fixGenBankFile(gbFile);
                        allDocs.addAll(PluginUtilities.importDocumentsToDatabase(
                            gbFile, targetFolder, progressListener));
                    }

                    commitTx.invoke(targetFolder);
                    Thread.sleep(3000);

                    beginTx.invoke(targetFolder);
                    List<AnnotatedPluginDocument> bamDocs = PluginUtilities.importDocumentsToDatabase(
                        bamFile, targetFolder, progressListener);
                    commitTx.invoke(targetFolder);

                    return bamDocs;

                } catch (Exception e) {
                    rollbackTx.invoke(targetFolder);
                    throw e;
                }

            } catch (NoSuchMethodException e) {
                System.out.println("Transaction methods not available, using standard import");
                return new Strategy1_DirectToSameFolder().importBam(bamFile, gbFiles, minimap2Dir, targetFolder, progressListener);
            }
        }
    }

    /**
     * Strategy 7: Background thread with timing
     */
    static class Strategy7_BackgroundThreadTiming implements ImportStrategy {
        public String getName() { return "Background Thread with Timing"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            System.out.println("Starting background reference import...");

            CompletableFuture<List<AnnotatedPluginDocument>> refFuture = CompletableFuture.supplyAsync(() -> {
                List<AnnotatedPluginDocument> refs = new ArrayList<>();
                try {
                    for (File gbFile : gbFiles) {
                        GenBankFileFixer.fixGenBankFile(gbFile);
                        refs.addAll(PluginUtilities.importDocumentsToDatabase(
                            gbFile, targetFolder, ProgressListener.EMPTY));
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
                return refs;
            });

            // Wait for references
            List<AnnotatedPluginDocument> refs = refFuture.get(30, TimeUnit.SECONDS);
            System.out.println("References imported: " + refs.size());

            // Extended settling time
            System.out.println("Waiting 5 seconds for database indexing...");
            Thread.sleep(5000);

            // Import BAM
            System.out.println("Importing BAM after background reference import...");
            return PluginUtilities.importDocumentsToDatabase(bamFile, targetFolder, progressListener);
        }
    }

    /**
     * Strategy 8: System properties configuration
     */
    static class Strategy8_SystemProperties implements ImportStrategy {
        public String getName() { return "System Properties Configuration"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            // Import references first
            System.out.println("Importing references...");
            for (File gbFile : gbFiles) {
                GenBankFileFixer.fixGenBankFile(gbFile);
                PluginUtilities.importDocumentsToDatabase(gbFile, targetFolder, progressListener);
            }

            Thread.sleep(2000);

            // Set system properties
            System.out.println("Setting system properties for import...");
            System.setProperty("geneious.import.references.autoFind", "true");
            System.setProperty("geneious.import.bam.referencePath", targetFolder.getName());
            System.setProperty("geneious.import.suppressDialogs", "true");
            System.setProperty("geneious.import.references.folder", targetFolder.getUniqueID());

            try {
                return PluginUtilities.importDocumentsToDatabase(bamFile, targetFolder, progressListener);
            } finally {
                // Clean up properties
                System.clearProperty("geneious.import.references.autoFind");
                System.clearProperty("geneious.import.bam.referencePath");
                System.clearProperty("geneious.import.suppressDialogs");
                System.clearProperty("geneious.import.references.folder");
            }
        }
    }

    /**
     * Strategy 9: Force database indexing
     */
    static class Strategy9_ForceIndexing implements ImportStrategy {
        public String getName() { return "Force Database Indexing"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            // Import references
            System.out.println("Importing references...");
            List<AnnotatedPluginDocument> refs = new ArrayList<>();
            for (File gbFile : gbFiles) {
                GenBankFileFixer.fixGenBankFile(gbFile);
                refs.addAll(PluginUtilities.importDocumentsToDatabase(gbFile, targetFolder, progressListener));
            }

            // Try to force indexing via reflection
            System.out.println("Attempting to force database indexing...");
            try {
                Method[] methods = targetFolder.getClass().getMethods();
                for (Method m : methods) {
                    if (m.getName().contains("index") || m.getName().contains("flush") ||
                        m.getName().contains("refresh") || m.getName().contains("sync")) {
                        try {
                            System.out.println("  Calling: " + m.getName());
                            m.invoke(targetFolder);
                        } catch (Exception e) {
                            // Ignore individual method failures
                        }
                    }
                }
            } catch (Exception e) {
                System.out.println("Could not force indexing: " + e.getMessage());
            }

            Thread.sleep(3000);

            // Import BAM
            System.out.println("Importing BAM after forced indexing...");
            return PluginUtilities.importDocumentsToDatabase(bamFile, targetFolder, progressListener);
        }
    }

    /**
     * Strategy 10: Reflection-based context injection
     */
    static class Strategy10_ReflectionContext implements ImportStrategy {
        public String getName() { return "Reflection-based Context Injection"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            // Import references
            System.out.println("Importing references...");
            List<AnnotatedPluginDocument> refs = new ArrayList<>();
            for (File gbFile : gbFiles) {
                GenBankFileFixer.fixGenBankFile(gbFile);
                refs.addAll(PluginUtilities.importDocumentsToDatabase(gbFile, targetFolder, progressListener));
            }

            Thread.sleep(2000);

            // Try to create import context via reflection
            System.out.println("Attempting reflection-based import context...");
            try {
                // Look for import context classes
                Class<?>[] classes = {
                    Class.forName("com.biomatters.geneious.publicapi.plugin.ImportContext"),
                    Class.forName("com.biomatters.geneious.publicapi.plugin.DocumentImportContext"),
                    Class.forName("com.biomatters.plugins.fileimportexport.ImportContext")
                };

                for (Class<?> contextClass : classes) {
                    try {
                        Object context = contextClass.newInstance();

                        // Try to set references
                        for (Method m : contextClass.getMethods()) {
                            if (m.getName().contains("Reference") || m.getName().contains("Database")) {
                                try {
                                    if (m.getParameterTypes().length == 1) {
                                        if (m.getParameterTypes()[0].isAssignableFrom(targetFolder.getClass())) {
                                            m.invoke(context, targetFolder);
                                        } else if (m.getParameterTypes()[0].isAssignableFrom(List.class)) {
                                            m.invoke(context, refs);
                                        }
                                    }
                                } catch (Exception e) {
                                    // Ignore
                                }
                            }
                        }

                        // Try to use context in import
                        Method importMethod = PluginUtilities.class.getMethod("importDocuments",
                            File.class, contextClass, ProgressListener.class);
                        List<AnnotatedPluginDocument> docs = (List<AnnotatedPluginDocument>)
                            importMethod.invoke(null, bamFile, context, progressListener);

                        if (docs != null && !docs.isEmpty()) {
                            return docs;
                        }

                    } catch (Exception e) {
                        // Try next context class
                    }
                }
            } catch (ClassNotFoundException e) {
                System.out.println("Import context classes not found");
            }

            // Fallback
            return PluginUtilities.importDocumentsToDatabase(bamFile, targetFolder, progressListener);
        }
    }

    /**
     * Strategy 11: Custom progress listener with reference context
     */
    static class Strategy11_CustomProgressListener implements ImportStrategy {
        public String getName() { return "Custom Progress Listener"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            // Import references
            System.out.println("Importing references...");
            List<AnnotatedPluginDocument> refs = new ArrayList<>();
            for (File gbFile : gbFiles) {
                GenBankFileFixer.fixGenBankFile(gbFile);
                refs.addAll(PluginUtilities.importDocumentsToDatabase(gbFile, targetFolder, progressListener));
            }

            // Cache references
            REFERENCE_CACHE.clear();
            for (AnnotatedPluginDocument ref : refs) {
                REFERENCE_CACHE.put(ref.getName(), ref);
            }

            Thread.sleep(2000);

            // Create custom progress listener wrapper
            System.out.println("Using custom progress listener wrapper...");

            // Since ProgressListener methods are final, we can't override them
            // Just use the standard listener but log reference cache status
            System.out.println("  Available references in cache: " + REFERENCE_CACHE.keySet());

            ProgressListener customListener = progressListener != null ? progressListener : ProgressListener.EMPTY;

            return PluginUtilities.importDocumentsToDatabase(bamFile, targetFolder, customListener);
        }
    }

    /**
     * Strategy 12: Symlink approach
     */
    static class Strategy12_SymlinkApproach implements ImportStrategy {
        public String getName() { return "Symbolic Link Co-location"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            Path tempDir = Files.createTempDirectory("bam_symlink_");

            try {
                System.out.println("Creating symbolic links for co-location...");

                // Create symlinks to GenBank files
                for (File gbFile : gbFiles) {
                    GenBankFileFixer.fixGenBankFile(gbFile);
                    Path linkPath = tempDir.resolve(gbFile.getName());
                    try {
                        Files.createSymbolicLink(linkPath, gbFile.toPath());
                    } catch (Exception e) {
                        // Fall back to copying if symlinks not supported
                        Files.copy(gbFile.toPath(), linkPath, StandardCopyOption.REPLACE_EXISTING);
                    }
                }

                // Create symlink to BAM file
                Path bamLink = tempDir.resolve(bamFile.getName());
                try {
                    Files.createSymbolicLink(bamLink, bamFile.toPath());
                } catch (Exception e) {
                    Files.copy(bamFile.toPath(), bamLink, StandardCopyOption.REPLACE_EXISTING);
                }

                // Import from symlinked directory
                System.out.println("Importing from symlinked directory...");
                List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(tempDir.toFile(), progressListener);

                // Add to database
                List<AnnotatedPluginDocument> addedDocs = new ArrayList<>();
                for (AnnotatedPluginDocument doc : docs) {
                    addedDocs.add(targetFolder.addDocumentCopy(doc, ProgressListener.EMPTY));
                }

                return addedDocs;

            } finally {
                // Clean up
                Files.list(tempDir).forEach(p -> {
                    try { Files.delete(p); } catch (Exception e) { }
                });
                Files.delete(tempDir);
            }
        }
    }

    /**
     * Strategy 13: Multi-step with reference cache
     */
    static class Strategy13_MultiStepWithCache implements ImportStrategy {
        public String getName() { return "Multi-step with Reference Cache"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            // Step 1: Import references to memory
            System.out.println("Step 1: Importing references to memory...");
            List<AnnotatedPluginDocument> memRefs = new ArrayList<>();
            for (File gbFile : gbFiles) {
                GenBankFileFixer.fixGenBankFile(gbFile);
                List<AnnotatedPluginDocument> docs = PluginUtilities.importDocuments(gbFile, progressListener);
                memRefs.addAll(docs);
            }

            // Step 2: Add references to database
            System.out.println("Step 2: Adding references to database...");
            List<AnnotatedPluginDocument> dbRefs = new ArrayList<>();
            for (AnnotatedPluginDocument doc : memRefs) {
                AnnotatedPluginDocument added = targetFolder.addDocumentCopy(doc, ProgressListener.EMPTY);
                if (added != null) {
                    dbRefs.add(added);
                    REFERENCE_CACHE.put(added.getName(), added);
                }
            }

            // Step 3: Extract BAM requirements
            System.out.println("Step 3: Extracting BAM requirements...");
            List<BamReferenceExtractor.ReferenceInfo> bamRefs = BamReferenceExtractor.extractReferences(bamFile);
            System.out.println("  BAM requires " + bamRefs.size() + " reference(s)");

            // Step 4: Match references
            Map<String, AnnotatedPluginDocument> refMap = new HashMap<>();
            for (BamReferenceExtractor.ReferenceInfo bamRef : bamRefs) {
                for (AnnotatedPluginDocument dbRef : dbRefs) {
                    if (dbRef.getName().equals(bamRef.name) ||
                        dbRef.getName().equals(bamRef.accession) ||
                        (bamRef.accession != null && dbRef.getName().contains(bamRef.accession))) {
                        refMap.put(bamRef.name, dbRef);
                        System.out.println("  Matched: " + bamRef.name + " -> " + dbRef.getName());
                        break;
                    }
                }
            }

            Thread.sleep(3000);

            // Step 5: Import BAM with matched references available
            System.out.println("Step 5: Importing BAM with " + refMap.size() + " matched reference(s)...");
            return PluginUtilities.importDocumentsToDatabase(bamFile, targetFolder, progressListener);
        }
    }

    /**
     * Strategy 14: Direct memory import without database
     */
    static class Strategy14_DirectMemoryImport implements ImportStrategy {
        public String getName() { return "Direct Memory Import"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            System.out.println("Importing everything to memory first...");

            // Import all to memory
            List<AnnotatedPluginDocument> allMemDocs = new ArrayList<>();

            for (File gbFile : gbFiles) {
                GenBankFileFixer.fixGenBankFile(gbFile);
                allMemDocs.addAll(PluginUtilities.importDocuments(gbFile, progressListener));
            }

            allMemDocs.addAll(PluginUtilities.importDocuments(bamFile, progressListener));

            // Now add all to database at once
            System.out.println("Adding " + allMemDocs.size() + " documents to database...");
            List<AnnotatedPluginDocument> dbDocs = new ArrayList<>();

            for (AnnotatedPluginDocument doc : allMemDocs) {
                AnnotatedPluginDocument added = targetFolder.addDocumentCopy(doc, ProgressListener.EMPTY);
                if (added != null) {
                    dbDocs.add(added);
                }
            }

            // Return only BAM documents
            List<AnnotatedPluginDocument> bamDocs = new ArrayList<>();
            for (AnnotatedPluginDocument doc : dbDocs) {
                if (doc.getName().toLowerCase().contains(".bam") ||
                    doc.getDocumentClass().getName().contains("Alignment")) {
                    bamDocs.add(doc);
                }
            }

            return bamDocs;
        }
    }

    /**
     * Strategy 15: Import to combined folder
     */
    static class Strategy15_CombinedFolderImport implements ImportStrategy {
        public String getName() { return "Combined Folder Import"; }

        public List<AnnotatedPluginDocument> importBam(File bamFile, List<File> gbFiles,
                                                       Path minimap2Dir, WritableDatabaseService targetFolder,
                                                       ProgressListener progressListener) throws Exception {
            // Create a subfolder specifically for this import
            System.out.println("Creating combined import folder...");
            WritableDatabaseService combinedFolder = null;

            try {
                // Try to create subfolder
                Method createFolder = targetFolder.getClass().getMethod("createSubFolder", String.class);
                combinedFolder = (WritableDatabaseService) createFolder.invoke(targetFolder, "BAM_Import_" + System.currentTimeMillis());
            } catch (Exception e) {
                System.out.println("Could not create subfolder, using main folder");
                combinedFolder = targetFolder;
            }

            // Import everything to the combined folder
            System.out.println("Importing all files to combined folder...");
            List<AnnotatedPluginDocument> allDocs = new ArrayList<>();

            // Import references
            for (File gbFile : gbFiles) {
                GenBankFileFixer.fixGenBankFile(gbFile);
                allDocs.addAll(PluginUtilities.importDocumentsToDatabase(gbFile, combinedFolder, progressListener));
            }

            Thread.sleep(2000);

            // Import BAM
            List<AnnotatedPluginDocument> bamDocs = PluginUtilities.importDocumentsToDatabase(
                bamFile, combinedFolder, progressListener);

            return bamDocs;
        }
    }

    @Override
    public List<AnnotatedPluginDocument> performOperation(AnnotatedPluginDocument[] documents,
                                                           ProgressListener progressListener,
                                                           Options options) throws DocumentOperationException {
        // This would be called when user runs the operation
        // For testing, use testAllImportStrategies() directly
        return new ArrayList<>();
    }
}