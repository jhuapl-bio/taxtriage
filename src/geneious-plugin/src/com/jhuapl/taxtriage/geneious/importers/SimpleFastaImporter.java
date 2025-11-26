package com.jhuapl.taxtriage.geneious.importers;

import com.biomatters.geneious.publicapi.documents.AnnotatedPluginDocument;
import com.biomatters.geneious.publicapi.documents.DocumentUtilities;
import com.biomatters.geneious.publicapi.documents.sequence.SequenceDocument;
import com.biomatters.geneious.publicapi.implementations.sequence.DefaultNucleotideSequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Simple FASTA importer for TaxTriage reference sequences.
 */
public class SimpleFastaImporter {

    private static final Logger logger = Logger.getLogger(SimpleFastaImporter.class.getName());

    /**
     * Imports a FASTA file and returns a list of AnnotatedPluginDocuments.
     *
     * @param fastaFile The FASTA file to import
     * @return List of imported sequence documents
     */
    public static List<AnnotatedPluginDocument> importFasta(File fastaFile) {
        List<AnnotatedPluginDocument> documents = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(fastaFile))) {
            String line;
            StringBuilder currentSequence = new StringBuilder();
            String currentName = "";

            while ((line = reader.readLine()) != null) {
                if (line.length() == 0) {
                    continue;
                }

                if (line.startsWith(">")) {
                    // Save previous sequence if exists
                    if (currentSequence.length() > 0 && !currentName.isEmpty()) {
                        try {
                            SequenceDocument seq = new DefaultNucleotideSequence(
                                currentName,
                                currentSequence.toString()
                            );
                            AnnotatedPluginDocument doc = DocumentUtilities.createAnnotatedPluginDocument(seq);
                            documents.add(doc);
                        } catch (Exception e) {
                            logger.log(Level.WARNING, "Failed to create sequence: " + currentName, e);
                        }
                        currentSequence = new StringBuilder();
                    }
                    // Start new sequence
                    currentName = line.substring(1).trim();
                } else {
                    // Add to current sequence
                    currentSequence.append(line.trim());
                }
            }

            // Don't forget the last sequence
            if (currentSequence.length() > 0 && !currentName.isEmpty()) {
                try {
                    SequenceDocument seq = new DefaultNucleotideSequence(
                        currentName,
                        currentSequence.toString()
                    );
                    AnnotatedPluginDocument doc = DocumentUtilities.createAnnotatedPluginDocument(seq);
                    documents.add(doc);
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Failed to create sequence: " + currentName, e);
                }
            }

            logger.info("Imported " + documents.size() + " sequences from " + fastaFile.getName());

        } catch (IOException e) {
            logger.log(Level.SEVERE, "Error reading FASTA file: " + fastaFile, e);
        }

        return documents;
    }
}