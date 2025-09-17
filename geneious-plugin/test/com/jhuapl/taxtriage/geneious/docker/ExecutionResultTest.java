package com.jhuapl.taxtriage.geneious.docker;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Test suite for ExecutionResult class.
 *
 * Tests all methods of the ExecutionResult class including construction,
 * accessors, and utility methods.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class ExecutionResultTest {

    @Test
    void testSuccessfulExecution() {
        ExecutionResult result = new ExecutionResult(0, "Success output", "");

        assertEquals(0, result.getExitCode());
        assertEquals("Success output", result.getOutput());
        assertEquals("", result.getErrorOutput());
        assertTrue(result.isSuccess());
    }

    @Test
    void testFailedExecution() {
        ExecutionResult result = new ExecutionResult(1, "Partial output", "Error occurred");

        assertEquals(1, result.getExitCode());
        assertEquals("Partial output", result.getOutput());
        assertEquals("Error occurred", result.getErrorOutput());
        assertFalse(result.isSuccess());
    }

    @Test
    void testConstructorWithNullValues() {
        ExecutionResult result = new ExecutionResult(0, null, null);

        assertEquals(0, result.getExitCode());
        assertEquals("", result.getOutput());
        assertEquals("", result.getErrorOutput());
        assertTrue(result.isSuccess());
    }

    @Test
    void testGetAllOutputWithBothOutputs() {
        ExecutionResult result = new ExecutionResult(0, "Standard output", "Error output");

        String allOutput = result.getAllOutput();
        assertTrue(allOutput.contains("STDOUT:"));
        assertTrue(allOutput.contains("Standard output"));
        assertTrue(allOutput.contains("STDERR:"));
        assertTrue(allOutput.contains("Error output"));
    }

    @Test
    void testGetAllOutputWithOnlyStdout() {
        ExecutionResult result = new ExecutionResult(0, "Only stdout", "");

        String allOutput = result.getAllOutput();
        assertTrue(allOutput.contains("STDOUT:"));
        assertTrue(allOutput.contains("Only stdout"));
        assertFalse(allOutput.contains("STDERR:"));
    }

    @Test
    void testGetAllOutputWithOnlyStderr() {
        ExecutionResult result = new ExecutionResult(1, "", "Only stderr");

        String allOutput = result.getAllOutput();
        assertFalse(allOutput.contains("STDOUT:"));
        assertTrue(allOutput.contains("STDERR:"));
        assertTrue(allOutput.contains("Only stderr"));
    }

    @Test
    void testGetAllOutputWithNoOutput() {
        ExecutionResult result = new ExecutionResult(0, "", "");

        String allOutput = result.getAllOutput();
        assertEquals("", allOutput);
    }

    @Test
    void testToString() {
        ExecutionResult result = new ExecutionResult(1, "Some output", "Some error");

        String toString = result.toString();
        assertTrue(toString.contains("exitCode=1"));
        assertTrue(toString.contains("outputLength=11"));
        assertTrue(toString.contains("errorOutputLength=10"));
        assertTrue(toString.contains("success=false"));
    }

    @Test
    void testToStringWithLongOutput() {
        String longOutput = "A".repeat(1000);
        String longError = "B".repeat(500);
        ExecutionResult result = new ExecutionResult(0, longOutput, longError);

        String toString = result.toString();
        assertTrue(toString.contains("outputLength=1000"));
        assertTrue(toString.contains("errorOutputLength=500"));
        assertTrue(toString.contains("success=true"));
    }

    @Test
    void testIsSuccessWithVariousExitCodes() {
        assertTrue(new ExecutionResult(0, "", "").isSuccess());
        assertFalse(new ExecutionResult(1, "", "").isSuccess());
        assertFalse(new ExecutionResult(-1, "", "").isSuccess());
        assertFalse(new ExecutionResult(127, "", "").isSuccess());
        assertFalse(new ExecutionResult(255, "", "").isSuccess());
    }

    @Test
    void testMultilineOutputs() {
        String multilineOutput = "Line 1\nLine 2\nLine 3";
        String multilineError = "Error line 1\nError line 2";

        ExecutionResult result = new ExecutionResult(0, multilineOutput, multilineError);

        assertEquals(multilineOutput, result.getOutput());
        assertEquals(multilineError, result.getErrorOutput());

        String allOutput = result.getAllOutput();
        assertTrue(allOutput.contains("Line 1"));
        assertTrue(allOutput.contains("Line 2"));
        assertTrue(allOutput.contains("Line 3"));
        assertTrue(allOutput.contains("Error line 1"));
        assertTrue(allOutput.contains("Error line 2"));
    }

    @Test
    void testSpecialCharactersInOutput() {
        String specialOutput = "Output with \"quotes\" and 'apostrophes' and \t tabs \n newlines";
        String specialError = "Error with unicode: \u2603 and symbols: $@#%";

        ExecutionResult result = new ExecutionResult(0, specialOutput, specialError);

        assertEquals(specialOutput, result.getOutput());
        assertEquals(specialError, result.getErrorOutput());

        String allOutput = result.getAllOutput();
        assertTrue(allOutput.contains("\"quotes\""));
        assertTrue(allOutput.contains("\u2603"));
        assertTrue(allOutput.contains("$@#%"));
    }
}