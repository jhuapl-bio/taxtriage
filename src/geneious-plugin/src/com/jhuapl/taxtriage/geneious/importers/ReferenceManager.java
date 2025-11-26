package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.databaseservice.DatabaseServiceException;
import com.biomatters.geneious.publicapi.databaseservice.Query;
import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.URN;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import jebl.util.ProgressListener;

import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Manages reference sequences for BAM import, ensuring they are properly
 * committed to the Geneious database before BAM files are imported.
 */
public class ReferenceManager {

    private static final Logger logger = Logger.getLogger(ReferenceManager.class.getName());

    /**
     * Prepares reference documents for BAM import by ensuring they have proper names
     * and are ready to be used. Since we can't directly commit to database from here,
     * we ensure the documents are properly configured.
     *
     * @param referenceDocuments The reference documents to prepare
     * @param bamReferences The BAM reference information for name matching
     * @param progressListener Progress listener
     * @return The prepared reference documents
     */
    public static List<AnnotatedPluginDocument> prepareReferencesForBAM(
            List<AnnotatedPluginDocument> referenceDocuments,
            List<BamReferenceExtractor.ReferenceInfo> bamReferences,
            ProgressListener progressListener) {

        logger.info("===== PREPARING REFERENCES FOR BAM =====");

        if (referenceDocuments == null || referenceDocuments.isEmpty()) {
            logger.warning("No reference documents to prepare");
            return referenceDocuments;
        }

        logger.info("Input: " + referenceDocuments.size() + " reference document(s)");
        logger.info("Input: " + bamReferences.size() + " BAM reference(s)");

        logger.info("BAM references expecting:");
        for (BamReferenceExtractor.ReferenceInfo bamRef : bamReferences) {
            logger.info("  - BAM expects: '" + bamRef.name + "' (accession: " + bamRef.accession + ")");
        }

        logger.info("Current document names:");
        for (AnnotatedPluginDocument doc : referenceDocuments) {
            logger.info("  - Doc has: '" + doc.getName() + "'");
        }

        if (progressListener != null) {
            progressListener.setMessage("Preparing reference sequences for BAM import...");
        }

        // Ensure each reference document has the correct name for BAM matching
        Map<String, String> nameMapping = createNameMappingForBAM(bamReferences, referenceDocuments);

        logger.info("Name mapping created with " + nameMapping.size() + " entries:");
        for (Map.Entry<String, String> entry : nameMapping.entrySet()) {
            logger.info("  - Map: '" + entry.getKey() + "' -> '" + entry.getValue() + "'");
        }

        int renameCount = 0;
        for (AnnotatedPluginDocument doc : referenceDocuments) {
            String currentName = doc.getName();
            String targetName = nameMapping.get(currentName);

            if (targetName != null && !targetName.equals(currentName)) {
                try {
                    doc.setName(targetName);
                    logger.info("RENAMED: '" + currentName + "' -> '" + targetName + "'");
                    renameCount++;
                } catch (Exception e) {
                    logger.log(Level.WARNING, "FAILED to rename: '" + currentName + "' -> '" + targetName + "'", e);
                }
            } else if (targetName == null) {
                logger.info("NO MAPPING for: '" + currentName + "'");
            } else {
                logger.info("ALREADY CORRECT: '" + currentName + "'");
            }
        }

        logger.info("Renamed " + renameCount + " document(s)");

        // Log the final prepared references
        logger.info("Final prepared references:");
        for (AnnotatedPluginDocument doc : referenceDocuments) {
            logger.info("  - '" + doc.getName() + "' (URN: " +
                       (doc.getURN() != null ? doc.getURN() : "pending") + ")");
        }

        logger.info("===== REFERENCE PREPARATION COMPLETE =====");

        return referenceDocuments;
    }

    /**
     * Creates a name mapping to ensure references match BAM expectations.
     */
    private static Map<String, String> createNameMappingForBAM(
            List<BamReferenceExtractor.ReferenceInfo> bamReferences,
            List<AnnotatedPluginDocument> referenceDocuments) {

        Map<String, String> mapping = new HashMap<>();

        // For each reference document, find the best matching BAM reference name
        for (AnnotatedPluginDocument doc : referenceDocuments) {
            String docName = doc.getName();
            if (docName == null) continue;

            // Extract accession from document name
            String docAccession = extractAccession(docName);

            // Find matching BAM reference
            for (BamReferenceExtractor.ReferenceInfo bamRef : bamReferences) {
                // Check if this BAM reference matches the document
                if (bamRef.accession != null && bamRef.accession.equals(docAccession)) {
                    // Map document name to BAM reference name
                    mapping.put(docName, bamRef.name);
                    logger.fine("Mapped " + docName + " to BAM name: " + bamRef.name);
                    break;
                } else if (bamRef.name.equals(docName)) {
                    // Exact match
                    mapping.put(docName, bamRef.name);
                    break;
                }
            }
        }

        return mapping;
    }

    /**
     * Creates a reference mapping between BAM reference names and Geneious document names.
     * This helps ensure BAM import can find the correct reference sequences.
     */
    public static Map<String, String> createReferenceMapping(List<BamReferenceExtractor.ReferenceInfo> bamReferences,
                                                            List<AnnotatedPluginDocument> referenceDocuments) {
        Map<String, String> mapping = new HashMap<>();

        // Create a map of document names for quick lookup
        Map<String, AnnotatedPluginDocument> docsByName = new HashMap<>();
        Map<String, AnnotatedPluginDocument> docsByAccession = new HashMap<>();

        for (AnnotatedPluginDocument doc : referenceDocuments) {
            String name = doc.getName();
            if (name != null) {
                docsByName.put(name, doc);

                // Also try to extract accession from name
                String accession = extractAccession(name);
                if (accession != null) {
                    docsByAccession.put(accession, doc);
                }
            }
        }

        // Try to match BAM references to documents
        for (BamReferenceExtractor.ReferenceInfo bamRef : bamReferences) {
            // First try exact name match
            if (docsByName.containsKey(bamRef.name)) {
                mapping.put(bamRef.name, bamRef.name);
                logger.fine("Exact match for reference: " + bamRef.name);
            }
            // Then try accession match
            else if (bamRef.accession != null && docsByAccession.containsKey(bamRef.accession)) {
                AnnotatedPluginDocument doc = docsByAccession.get(bamRef.accession);
                mapping.put(bamRef.name, doc.getName());
                logger.fine("Accession match for reference: " + bamRef.name + " -> " + doc.getName());
            }
            // Try partial matching
            else {
                for (String docName : docsByName.keySet()) {
                    if (docName.contains(bamRef.accession) || bamRef.name.contains(docName)) {
                        mapping.put(bamRef.name, docName);
                        logger.fine("Partial match for reference: " + bamRef.name + " -> " + docName);
                        break;
                    }
                }
            }

            if (!mapping.containsKey(bamRef.name)) {
                logger.warning("No match found for BAM reference: " + bamRef.name);
            }
        }

        return mapping;
    }

    private static String extractAccession(String name) {
        if (name == null) return null;

        // Common patterns for NCBI accessions
        java.util.regex.Pattern pattern = java.util.regex.Pattern.compile("([A-Z]{1,4}_?\\d+\\.?\\d*)");
        java.util.regex.Matcher matcher = pattern.matcher(name);

        if (matcher.find()) {
            return matcher.group(1);
        }

        return null;
    }

}