package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.databaseservice.DatabaseService;
import com.biomatters.geneious.publicapi.databaseservice.Query;
import com.biomatters.geneious.publicapi.databaseservice.QueryField;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import jebl.util.ProgressListener;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Helper class to check if reference sequences are available in Geneious database
 * before attempting BAM import.
 */
public class ReferenceAvailabilityChecker {

    private static final Logger logger = Logger.getLogger(ReferenceAvailabilityChecker.class.getName());

    /**
     * Checks if all required references for a BAM file are available in the database.
     *
     * @param bamFile The BAM file that requires references
     * @return Map of reference names to their availability status
     */
    public static Map<String, Boolean> checkReferencesAvailable(File bamFile) {
        return checkReferencesAvailable(bamFile, null);
    }

    /**
     * Checks if all required references for a BAM file are available in the database.
     *
     * @param bamFile The BAM file that requires references
     * @param knownReferences Optional list of known reference documents to check against
     * @return Map of reference names to their availability status
     */
    public static Map<String, Boolean> checkReferencesAvailable(File bamFile, List<AnnotatedPluginDocument> knownReferences) {
        Map<String, Boolean> availability = new HashMap<>();

        System.out.println("===== CHECKING REFERENCE AVAILABILITY =====");

        try {
            // Extract reference names from BAM
            List<BamReferenceExtractor.ReferenceInfo> bamRefs = BamReferenceExtractor.extractReferences(bamFile);
            System.out.println("BAM file requires " + bamRefs.size() + " reference(s)");

            // Get database service
            DatabaseService database = null;
            try {
                Object service = PluginUtilities.getGeneiousService("DatabaseService");
                if (service instanceof DatabaseService) {
                    database = (DatabaseService) service;
                    System.out.println("Got DatabaseService for searching");
                }
            } catch (Exception e) {
                System.out.println("Could not get DatabaseService: " + e.getMessage());
            }

            // Try alternative ways to search for references
            for (BamReferenceExtractor.ReferenceInfo ref : bamRefs) {
                boolean found = false;

                System.out.println("\nSearching for reference: '" + ref.name + "'");

                // Method 1: Check in provided known references first
                if (knownReferences != null && !knownReferences.isEmpty()) {
                    for (AnnotatedPluginDocument doc : knownReferences) {
                        if (doc.getName() != null && doc.getName().equals(ref.name)) {
                            System.out.println("  FOUND in provided references (URN: " + doc.getURN() + ")");
                            found = true;
                            break;
                        }
                    }
                }

                // Method 2: Try database query if we have a database service
                if (!found && database != null) {
                    try {
                        // Create a query for the reference name
                        Query query = Query.Factory.createQuery("Name=\"" + ref.name + "\"");
                        List<AnnotatedPluginDocument> results = database.retrieve(
                            query,
                            ProgressListener.EMPTY
                        );

                        if (!results.isEmpty()) {
                            System.out.println("  FOUND via database query: " + results.size() + " result(s)");
                            found = true;
                        }
                    } catch (Exception e) {
                        System.out.println("  Database query error: " + e.getMessage());
                    }
                }

                // Method 3: Check if found via other means
                if (!found) {
                    System.out.println("  Could not find reference via available search methods");
                }

                availability.put(ref.name, found);
                System.out.println("  Result: " + (found ? "AVAILABLE" : "NOT AVAILABLE"));
            }

        } catch (Exception e) {
            System.out.println("Error checking reference availability: " + e.getMessage());
            e.printStackTrace();
            logger.log(Level.WARNING, "Failed to check reference availability", e);
        }

        // Summary
        System.out.println("\n===== REFERENCE AVAILABILITY SUMMARY =====");
        int available = 0;
        int missing = 0;
        for (Map.Entry<String, Boolean> entry : availability.entrySet()) {
            if (entry.getValue()) {
                available++;
                System.out.println("  ✓ " + entry.getKey());
            } else {
                missing++;
                System.out.println("  ✗ " + entry.getKey());
            }
        }
        System.out.println("Total: " + available + " available, " + missing + " missing");
        System.out.println("========================================\n");

        return availability;
    }

    /**
     * Waits for references to become available, checking periodically.
     *
     * @param bamFile The BAM file requiring references
     * @param maxWaitMs Maximum time to wait in milliseconds
     * @param checkIntervalMs How often to check in milliseconds
     * @return true if all references became available, false if timeout
     */
    public static boolean waitForReferencesAvailable(File bamFile, int maxWaitMs, int checkIntervalMs) {
        System.out.println("Waiting for references to become available (max " + (maxWaitMs/1000) + " seconds)...");

        long startTime = System.currentTimeMillis();
        long endTime = startTime + maxWaitMs;

        while (System.currentTimeMillis() < endTime) {
            Map<String, Boolean> availability = checkReferencesAvailable(bamFile);

            // Check if all are available
            boolean allAvailable = true;
            for (Boolean available : availability.values()) {
                if (!available) {
                    allAvailable = false;
                    break;
                }
            }

            if (allAvailable) {
                System.out.println("All references are now available!");
                return true;
            }

            // Wait before checking again
            try {
                long elapsed = System.currentTimeMillis() - startTime;
                System.out.println("References not yet available after " + (elapsed/1000) + " seconds, waiting...");
                Thread.sleep(checkIntervalMs);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                break;
            }
        }

        System.out.println("Timeout: References did not become available within " + (maxWaitMs/1000) + " seconds");
        return false;
    }
}