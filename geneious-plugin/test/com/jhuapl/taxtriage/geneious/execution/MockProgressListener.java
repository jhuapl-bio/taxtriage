package com.jhuapl.taxtriage.geneious.execution;

import jebl.util.ProgressListener;

import java.util.ArrayList;
import java.util.List;

/**
 * Mock implementation of ProgressListener for testing.
 *
 * This class captures all progress updates and messages for verification
 * in unit tests, and allows simulation of user cancellation.
 */
public class MockProgressListener implements ProgressListener {

    private final List<Double> progressValues = new ArrayList<>();
    private final List<String> messages = new ArrayList<>();
    private double currentProgress = 0.0;
    private String currentMessage = "";
    private boolean cancelled = false;

    @Override
    public void setProgress(double progress) {
        this.currentProgress = progress;
        this.progressValues.add(progress);
    }

    @Override
    public void setMessage(String message) {
        this.currentMessage = message;
        this.messages.add(message);
    }

    @Override
    public boolean isCancelled() {
        return cancelled;
    }

    /**
     * Simulates user cancellation.
     */
    public void cancel() {
        this.cancelled = true;
    }

    /**
     * Resets the progress listener to initial state.
     */
    public void reset() {
        progressValues.clear();
        messages.clear();
        currentProgress = 0.0;
        currentMessage = "";
        cancelled = false;
    }

    /**
     * Gets the current progress value.
     *
     * @return current progress (0.0 to 1.0)
     */
    public double getCurrentProgress() {
        return currentProgress;
    }

    /**
     * Gets the current message.
     *
     * @return current message
     */
    public String getCurrentMessage() {
        return currentMessage;
    }

    /**
     * Gets all progress values that were set.
     *
     * @return list of progress values
     */
    public List<Double> getProgressValues() {
        return new ArrayList<>(progressValues);
    }

    /**
     * Gets all messages that were set.
     *
     * @return list of messages
     */
    public List<String> getMessages() {
        return new ArrayList<>(messages);
    }

    /**
     * Gets the number of progress updates.
     *
     * @return number of progress updates
     */
    public int getProgressUpdateCount() {
        return progressValues.size();
    }

    /**
     * Gets the number of message updates.
     *
     * @return number of message updates
     */
    public int getMessageUpdateCount() {
        return messages.size();
    }

    /**
     * Checks if progress reached the specified value.
     *
     * @param expectedProgress the expected progress value
     * @return true if progress reached the value
     */
    public boolean progressReached(double expectedProgress) {
        return progressValues.stream().anyMatch(p -> Math.abs(p - expectedProgress) < 0.001);
    }

    /**
     * Checks if a message containing the specified text was set.
     *
     * @param expectedText the expected text in message
     * @return true if message containing text was found
     */
    public boolean messageContains(String expectedText) {
        return messages.stream().anyMatch(m -> m.contains(expectedText));
    }

    /**
     * Gets the final (last) progress value.
     *
     * @return final progress value, or 0.0 if no progress was set
     */
    public double getFinalProgress() {
        return progressValues.isEmpty() ? 0.0 : progressValues.get(progressValues.size() - 1);
    }

    /**
     * Gets the final (last) message.
     *
     * @return final message, or empty string if no message was set
     */
    public String getFinalMessage() {
        return messages.isEmpty() ? "" : messages.get(messages.size() - 1);
    }

    /**
     * Checks if progress is monotonically increasing.
     *
     * @return true if progress values are non-decreasing
     */
    public boolean isProgressMonotonic() {
        for (int i = 1; i < progressValues.size(); i++) {
            if (progressValues.get(i) < progressValues.get(i - 1)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Gets a summary of all captured progress and messages.
     *
     * @return formatted summary
     */
    public String getSummary() {
        StringBuilder summary = new StringBuilder();
        summary.append("Progress Updates: ").append(progressValues.size()).append("\n");
        summary.append("Message Updates: ").append(messages.size()).append("\n");
        summary.append("Final Progress: ").append(String.format("%.1f%%", currentProgress * 100)).append("\n");
        summary.append("Final Message: ").append(currentMessage).append("\n");
        summary.append("Cancelled: ").append(cancelled).append("\n");

        if (!progressValues.isEmpty()) {
            summary.append("Progress Range: ")
                   .append(String.format("%.3f", progressValues.stream().mapToDouble(Double::doubleValue).min().orElse(0.0)))
                   .append(" - ")
                   .append(String.format("%.3f", progressValues.stream().mapToDouble(Double::doubleValue).max().orElse(0.0)))
                   .append("\n");
        }

        return summary.toString();
    }

    @Override
    public String toString() {
        return String.format("MockProgressListener{progress=%.3f, message='%s', cancelled=%s}",
                currentProgress, currentMessage, cancelled);
    }
}