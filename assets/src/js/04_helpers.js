      /* ═══════════════════════════════════════════════════════════════════════════
       -  §  HELPERS  (cross-tab small functions)
       -     rankLabel() — pull the right taxonomy field from a row given the
       -     currently selected rank. Used by Heatmap, TASS, Sunburst, Coverage.
═══════════════════════════════════════════════════════════════════════════ */
      /** Return the rank label string for a row given the current rank setting. */
      function rankLabel(row, rank) {
        if (!row) return "";
        const field = rank === "Species" ? "Detected Organism" : rank;
        return (row[field] || "") + (rank !== "Species" ? ` (${rank})` : "");
      }

