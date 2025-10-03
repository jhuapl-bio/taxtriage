package com.jhuapl.taxtriage.geneious.docker;

import org.junit.jupiter.api.Test;

import java.time.LocalDateTime;

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
        LocalDateTime now = LocalDateTime.now();
        ExecutionResult result = new ExecutionResult("test command", 0, "Success output", "", now, now);

        assertEquals(0, result.getExitCode());
        assertEquals("Success output", result.getStandardOutput());
        assertEquals("", result.getErrorOutput());
        assertTrue(result.isSuccess());
        assertTrue(result.isSuccessful());
    }

    @Test
    void testFailedExecution() {
        LocalDateTime now = LocalDateTime.now();
        ExecutionResult result = new ExecutionResult("test command", 1, "Partial output", "Error occurred", now, now);

        assertEquals(1, result.getExitCode());
        assertEquals("Partial output", result.getStandardOutput());
        assertEquals("Error occurred", result.getErrorOutput());
        assertFalse(result.isSuccess());
        assertFalse(result.isSuccessful());
        assertTrue(result.isFailed());
    }

    @Test
    void testConstructorWithNullValues() {
        LocalDateTime now = LocalDateTime.now();
        ExecutionResult result = new ExecutionResult("test command", 0, null, null, now, now);

        assertEquals(0, result.getExitCode());
        assertEquals("", result.getStandardOutput());
        assertEquals("", result.getErrorOutput());
        assertTrue(result.isSuccess());
    }

    @Test
    void testStaticFactoryMethods() {
        ExecutionResult success = ExecutionResult.success("test command", "Success output");
        assertTrue(success.isSuccessful());
        assertEquals(0, success.getExitCode());
        assertEquals("Success output", success.getStandardOutput());

        ExecutionResult failure = ExecutionResult.failure("test command", 1, "Error output");
        assertTrue(failure.isFailed());
        assertEquals(1, failure.getExitCode());
        assertEquals("Error output", failure.getErrorOutput());

        ExecutionResult fromException = ExecutionResult.fromException("test command", new RuntimeException("Test error"));
        assertTrue(fromException.isFailed());
        assertEquals(-1, fromException.getExitCode());
        assertTrue(fromException.getErrorOutput().contains("Test error"));
    }

    @Test
    void testOutputLines() {
        LocalDateTime now = LocalDateTime.now();
        String multilineOutput = "Line 1\nLine 2\nLine 3";
        String multilineError = "Error line 1\nError line 2";

        ExecutionResult result = new ExecutionResult("test command", 0, multilineOutput, multilineError, now, now);

        assertEquals(3, result.getOutputLines().size());
        assertEquals(2, result.getErrorLines().size());
        assertEquals("Line 1", result.getOutputLines().get(0));
        assertEquals("Error line 1", result.getErrorLines().get(0));
    }

    @Test
    void testToString() {
        LocalDateTime now = LocalDateTime.now();
        ExecutionResult result = new ExecutionResult("test command", 1, "Some output", "Some error", now, now);

        String toString = result.toString();
        assertTrue(toString.contains("exitCode=1"));
        assertTrue(toString.contains("successful=false"));
        assertTrue(toString.contains("test command"));
    }

    @Test
    void testIsSuccessWithVariousExitCodes() {
        LocalDateTime now = LocalDateTime.now();
        assertTrue(new ExecutionResult("cmd", 0, "", "", now, now).isSuccess());
        assertFalse(new ExecutionResult("cmd", 1, "", "", now, now).isSuccess());
        assertFalse(new ExecutionResult("cmd", -1, "", "", now, now).isSuccess());
        assertFalse(new ExecutionResult("cmd", 127, "", "", now, now).isSuccess());
        assertFalse(new ExecutionResult("cmd", 255, "", "", now, now).isSuccess());
    }

    @Test
    void testContainsInOutput() {
        LocalDateTime now = LocalDateTime.now();
        String output = "Line 1\nLine 2 with pattern\nLine 3";
        ExecutionResult result = new ExecutionResult("test command", 0, output, "", now, now);

        assertTrue(result.containsInOutput("pattern"));
        assertFalse(result.containsInOutput("notfound"));
    }

    @Test
    void testContainsInError() {
        LocalDateTime now = LocalDateTime.now();
        String error = "Error line 1\nError line 2 with error pattern\nError line 3";
        ExecutionResult result = new ExecutionResult("test command", 1, "", error, now, now);

        assertTrue(result.containsInError("error pattern"));
        assertFalse(result.containsInError("notfound"));
    }

    @Test
    void testFindInOutput() {
        LocalDateTime now = LocalDateTime.now();
        String output = "Line 1\nLine 2 with pattern\nLine 3";
        ExecutionResult result = new ExecutionResult("test command", 0, output, "", now, now);

        String found = result.findInOutput("pattern");
        assertNotNull(found);
        assertTrue(found.contains("pattern"));

        String notFound = result.findInOutput("notfound");
        assertNull(notFound);
    }

    @Test
    void testFindInError() {
        LocalDateTime now = LocalDateTime.now();
        String error = "Error line 1\nError line 2 with error pattern\nError line 3";
        ExecutionResult result = new ExecutionResult("test command", 1, "", error, now, now);

        String found = result.findInError("error pattern");
        assertNotNull(found);
        assertTrue(found.contains("error pattern"));

        String notFound = result.findInError("notfound");
        assertNull(notFound);
    }

    @Test
    void testHasErrorOutput() {
        LocalDateTime now = LocalDateTime.now();
        ExecutionResult withError = new ExecutionResult("test command", 1, "", "Error", now, now);
        ExecutionResult withoutError = new ExecutionResult("test command", 0, "Output", "", now, now);

        assertTrue(withError.hasErrorOutput());
        assertFalse(withoutError.hasErrorOutput());
    }

    @Test
    void testExecutionTime() {
        LocalDateTime start = LocalDateTime.now();
        LocalDateTime end = start.plusSeconds(10);
        ExecutionResult result = new ExecutionResult("test command", 0, "", "", start, end);

        assertEquals(10, result.getExecutionSeconds());
        assertTrue(result.getExecutionMillis() >= 10000);
    }

    @Test
    void testGetSummary() {
        LocalDateTime now = LocalDateTime.now();
        ExecutionResult result = new ExecutionResult("test command", 0, "Output", "Error", now, now);

        String summary = result.getSummary();
        assertTrue(summary.contains("test command"));
        assertTrue(summary.contains("Exit Code: 0"));
        assertTrue(summary.contains("SUCCESS"));
    }
}