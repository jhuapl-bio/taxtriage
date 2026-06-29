      /* ═══════════════════════════════════════════════════════════════════════════
       -  §  PROTEIN COLUMN RENAME MAP
       -     Static lookup mapping the noisy upstream column names from VFDB /
       -     CARD / RecPort exports to short, friendly labels used in the
       -     VF/AMR tab UI. Lookup is exact first, then case-insensitive.
       -     Columns in PROT_COL_IGNORE are dropped entirely. Extend
       -     window.PROT_COL_REMAP at runtime to add more aliases.
═══════════════════════════════════════════════════════════════════════════ */
      window.PROT_COL_REMAP = {
        evalue: "E-value",
        bitscore: "Bitscore",
        source_id: "Source ID",
        source: "Source",
        gene_name: "Gene",
        pident: "%id",
        property: "Property",
        genus: "Genus",
        species: "Species",
        level: "Level",
        class: "Property",
        product: "Product",
        antibiotics_class: "Antibiotics Class",
        antibiotics: "Antibiotics",
        organism: "Reference Organism",
        sample: "Specimen ID",
        sample_name: "Specimen ID",
        "sample name": "Specimen ID",
        "sample id": "Specimen ID",
        "specimen id": "Specimen ID",
      };

      window.PROT_COL_IGNORE = new Set([
        "reference organism",
        "classification",
        "antibiotics class",
        "antibiotics_class",
        "Antibiotics Class",
        "sseqid",
        "host_name",
      ]);

      /** Apply PROT_COL_REMAP and drop PROT_COL_IGNORE columns from an array of rows. */
      function _applyProtColRemap(rows) {
        const map = window.PROT_COL_REMAP || {};
        const ignore = window.PROT_COL_IGNORE || new Set();
        return rows.map((row) => {
          const out = {};
          for (const [k, v] of Object.entries(row)) {
            if (ignore.has(k.toLowerCase())) continue;
            const mapped = map[k] ?? map[k.toLowerCase()];
            out[mapped !== undefined ? mapped : k] = v;
          }
          return out;
        });
      }

