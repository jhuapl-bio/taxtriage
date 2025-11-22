package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.databaseservice.DatabaseService;
import com.biomatters.geneious.publicapi.databaseservice.DatabaseServiceException;
import com.biomatters.geneious.publicapi.databaseservice.Query;
import com.biomatters.geneious.publicapi.databaseservice.WritableDatabaseService;
import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.documents.URN;
import com.biomatters.geneious.publicapi.plugin.PluginUtilities;
import jebl.util.ProgressListener;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Helper to force indexing of reference documents in Geneious database.
 */
public class ForceIndexingHelper {

    private static final Logger logger = Logger.getLogger(ForceIndexingHelper.class.getName());

    /**
     * Forces Geneious to index reference documents by triggering a database refresh.
     */
    public static void forceIndexing(List<AnnotatedPluginDocument> documents) {
        System.out.println("===== FORCING DATABASE INDEXING =====");

        try {
            // Get the database service
            DatabaseService databaseService = null;
            WritableDatabaseService writableService = null;

            try {
                Object service = PluginUtilities.getGeneiousService("DatabaseService");
                if (service instanceof DatabaseService) {
                    databaseService = (DatabaseService) service;
                    System.out.println("Got DatabaseService");
                }

                service = PluginUtilities.getGeneiousService("WritableDatabaseService");
                if (service instanceof WritableDatabaseService) {
                    writableService = (WritableDatabaseService) service;
                    System.out.println("Got WritableDatabaseService");
                }
            } catch (Exception e) {
                System.out.println("Could not get database service: " + e.getMessage());
            }

            // Try to trigger a refresh/reindex
            if (writableService != null) {
                try {
                    System.out.println("Attempting to trigger database refresh...");

                    // Try to get each document by URN to force caching
                    for (AnnotatedPluginDocument doc : documents) {
                        if (doc.getURN() != null) {
                            try {
                                System.out.println("  Forcing cache for: " + doc.getName() + " (URN: " + doc.getURN() + ")");

                                // Try to retrieve by URN to force indexing
                                if (databaseService != null) {
                                    try {
                                        // Use DocumentUtilities to get by URN
                                        AnnotatedPluginDocument retrieved = DocumentUtilities.getDocumentByURN(doc.getURN());
                                        if (retrieved != null) {
                                            System.out.println("    Successfully retrieved from database");
                                        }
                                    } catch (Exception ex) {
                                        // Try with database query
                                        try {
                                            Query query = Query.Factory.createQuery("Name=\"" + doc.getName() + "\"");
                                            List<AnnotatedPluginDocument> results = databaseService.retrieve(
                                                query,
                                                ProgressListener.EMPTY
                                            );
                                            if (!results.isEmpty()) {
                                                System.out.println("    Found via query: " + results.size() + " result(s)");
                                            }
                                        } catch (Exception qe) {
                                            System.out.println("    Query failed: " + qe.getMessage());
                                        }
                                    }
                                }
                            } catch (Exception e) {
                                System.out.println("    Could not retrieve: " + e.getMessage());
                            }
                        }
                    }

                    // Add a longer wait for indexing to complete
                    System.out.println("Waiting 5 seconds for database indexing to complete...");
                    Thread.sleep(5000);
                    System.out.println("Indexing wait complete");

                } catch (Exception e) {
                    System.out.println("Error during indexing: " + e.getMessage());
                    logger.log(Level.WARNING, "Failed to force indexing", e);
                }
            }

        } catch (Exception e) {
            System.out.println("Unexpected error in forceIndexing: " + e.getMessage());
            logger.log(Level.WARNING, "Error in forceIndexing", e);
        }

        System.out.println("=====================================\n");
    }

    /**
     * Waits for documents to be fully indexed and searchable.
     */
    public static void waitForIndexing(List<AnnotatedPluginDocument> documents, int maxWaitMs) {
        System.out.println("===== WAITING FOR DOCUMENT INDEXING =====");
        System.out.println("Waiting up to " + (maxWaitMs/1000) + " seconds for " + documents.size() + " document(s) to be indexed");

        long startTime = System.currentTimeMillis();
        long endTime = startTime + maxWaitMs;

        // Initial wait
        try {
            Thread.sleep(2000);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }

        // Check periodically
        while (System.currentTimeMillis() < endTime) {
            boolean allIndexed = true;

            // Check if all documents have URNs (indicating they're in the database)
            for (AnnotatedPluginDocument doc : documents) {
                if (doc.getURN() == null) {
                    allIndexed = false;
                    System.out.println("  Document not yet indexed: " + doc.getName());
                    break;
                }
            }

            if (allIndexed) {
                System.out.println("All documents appear to be indexed");

                // Extra wait to ensure search index is updated
                try {
                    System.out.println("Waiting additional 3 seconds for search index update...");
                    Thread.sleep(3000);
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                }

                break;
            }

            // Wait before checking again
            try {
                Thread.sleep(1000);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                break;
            }
        }

        long elapsed = System.currentTimeMillis() - startTime;
        System.out.println("Indexing wait completed after " + (elapsed/1000) + " seconds");
        System.out.println("=========================================\n");
    }
}