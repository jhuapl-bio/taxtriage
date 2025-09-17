package com.jhuapl.taxtriage.geneious.docker;

import org.junit.jupiter.api.Test;

import java.io.IOException;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Test suite for DockerException class.
 *
 * Tests all constructors and inheritance behavior of the DockerException class.
 *
 * @author Johns Hopkins University Applied Physics Laboratory
 * @version 1.0.0
 */
class DockerExceptionTest {

    @Test
    void testConstructorWithMessage() {
        String message = "Docker operation failed";
        DockerException exception = new DockerException(message);

        assertEquals(message, exception.getMessage());
        assertNull(exception.getCause());
    }

    @Test
    void testConstructorWithMessageAndCause() {
        String message = "Docker operation failed";
        IOException cause = new IOException("IO error");
        DockerException exception = new DockerException(message, cause);

        assertEquals(message, exception.getMessage());
        assertEquals(cause, exception.getCause());
    }

    @Test
    void testConstructorWithCause() {
        IOException cause = new IOException("IO error");
        DockerException exception = new DockerException(cause);

        assertEquals(IOException.class.getName() + ": IO error", exception.getMessage());
        assertEquals(cause, exception.getCause());
    }

    @Test
    void testInheritanceFromException() {
        DockerException exception = new DockerException("Test message");

        assertTrue(exception instanceof Exception);
        assertTrue(exception instanceof Throwable);
    }

    @Test
    void testExceptionChaining() {
        IOException rootCause = new IOException("Root cause");
        RuntimeException intermediateCause = new RuntimeException("Intermediate", rootCause);
        DockerException exception = new DockerException("Final message", intermediateCause);

        assertEquals("Final message", exception.getMessage());
        assertEquals(intermediateCause, exception.getCause());
        assertEquals(rootCause, exception.getCause().getCause());
    }

    @Test
    void testNullMessage() {
        DockerException exception = new DockerException((String) null);

        assertNull(exception.getMessage());
        assertNull(exception.getCause());
    }

    @Test
    void testNullCause() {
        DockerException exception = new DockerException((Throwable) null);

        assertEquals("null", exception.getMessage());
        assertNull(exception.getCause());
    }

    @Test
    void testEmptyMessage() {
        String emptyMessage = "";
        DockerException exception = new DockerException(emptyMessage);

        assertEquals(emptyMessage, exception.getMessage());
        assertNull(exception.getCause());
    }

    @Test
    void testLongMessage() {
        String longMessage = "A".repeat(1000);
        DockerException exception = new DockerException(longMessage);

        assertEquals(longMessage, exception.getMessage());
        assertEquals(1000, exception.getMessage().length());
    }

    @Test
    void testMessageWithSpecialCharacters() {
        String specialMessage = "Error with unicode: \u2603 and symbols: $@#% and newlines:\nLine 2";
        DockerException exception = new DockerException(specialMessage);

        assertEquals(specialMessage, exception.getMessage());
        assertTrue(exception.getMessage().contains("\u2603"));
        assertTrue(exception.getMessage().contains("$@#%"));
        assertTrue(exception.getMessage().contains("\n"));
    }

    @Test
    void testStackTrace() {
        DockerException exception = new DockerException("Test exception");

        assertNotNull(exception.getStackTrace());
        assertTrue(exception.getStackTrace().length > 0);

        // The top stack trace element should be from this test method
        StackTraceElement topElement = exception.getStackTrace()[0];
        assertEquals("testStackTrace", topElement.getMethodName());
        assertEquals(getClass().getName(), topElement.getClassName());
    }

    @Test
    void testSuppressedExceptions() {
        DockerException mainException = new DockerException("Main exception");
        IOException suppressedException = new IOException("Suppressed exception");

        mainException.addSuppressed(suppressedException);

        Throwable[] suppressed = mainException.getSuppressed();
        assertEquals(1, suppressed.length);
        assertEquals(suppressedException, suppressed[0]);
    }

    @Test
    void testExceptionAsString() {
        DockerException exception = new DockerException("Test message");
        String exceptionString = exception.toString();

        assertTrue(exceptionString.contains("DockerException"));
        assertTrue(exceptionString.contains("Test message"));
    }

    @Test
    void testCausePreservation() {
        Exception originalCause = new RuntimeException("Original");
        DockerException wrappedException = new DockerException("Wrapped", originalCause);

        // Verify that the original cause is preserved exactly
        assertSame(originalCause, wrappedException.getCause());
    }
}