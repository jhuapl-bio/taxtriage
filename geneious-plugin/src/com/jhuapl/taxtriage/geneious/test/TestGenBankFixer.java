package com.jhuapl.taxtriage.geneious.test;

import com.jhuapl.taxtriage.geneious.importers.GenBankFileFixer;
import java.io.File;

/**
 * Test program for GenBankFileFixer
 */
public class TestGenBankFixer {
    public static void main(String[] args) {
        File testFile = new File("/tmp/test.gb");

        System.out.println("Testing GenBankFileFixer on: " + testFile.getAbsolutePath());
        System.out.println("Valid before fix: " + GenBankFileFixer.validateGenBankFile(testFile));

        boolean fixed = GenBankFileFixer.fixGenBankFile(testFile);
        System.out.println("Fix successful: " + fixed);

        System.out.println("Valid after fix: " + GenBankFileFixer.validateGenBankFile(testFile));
    }
}