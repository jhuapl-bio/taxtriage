      const BOOT = window.HEATMAP_BOOT || {};
      let DATA = BOOT.records || [];
      let ALL_COLS = (BOOT.all_cols || []).filter((c) => c !== "High ANI Matches" && c !== "ANI Annotated");
      // Analysis-only fields carried per-record (ANI list + capability flag) that
      // must never be exposed as detections-table columns.
      const _NON_DISPLAY_COLS = new Set(["High ANI Matches", "ANI Annotated"]);
      let NUMERIC = new Set(BOOT.numeric_cols || []);
      let SAMPLE_META = BOOT.sample_meta || {};
      const PROT = BOOT.prot_data || {};
      let HAS_PROT = BOOT.has_prot || false;
      const PROT_HIDDEN_PROPS = new Set();
      // Reference-free novelty payload from NOVELTY_COLLECT (per-sample summary + candidate taxa).
      const NOVELTY = BOOT.novelty || { samples: {} };
      const NOVELTY_DL = BOOT.novelty_downloads || [];
      let HAS_NOVELTY = BOOT.has_novelty || false;
      // Pathogen reference lookups (assets/pathogen_sheet.csv) used to flag listed
      // pathogens that have NO reference alignment but DO appear in novelty candidates
      // (pre-stamped server-side as candidate.pathogen) or in VF/AMR hits (matched here
      // by species/genus name). Shape: {by_taxid:{}, by_name:{}, by_genus:{}}.
      const PATHO = BOOT.pathogens || { by_taxid: {}, by_name: {}, by_genus: {} };
      const HAS_PATHO = !!BOOT.has_pathogens;

      // Which backend produced the novelty results (mmseqs2 | kaiju | bracken). Prefer the
      // top-level field; fall back to any sample's summary.classifier for older feeds.
      function _noveltyClassifier() {
        let c = (NOVELTY && NOVELTY.classifier) || "";
        if (!c && NOVELTY && NOVELTY.samples) {
          for (const s of Object.values(NOVELTY.samples)) {
            const cc = s && s.summary && s.summary.classifier;
            if (cc) {
              c = cc;
              break;
            }
          }
        }
        return String(c || "").toLowerCase();
      }
      // Gene mode (default on; off under --disable_gene): the query fed to the backend was
      // Pyrodigal-predicted genes, not whole de novo contigs. Read the top-level flag, falling
      // back to any sample summary.
      function _novGeneMode() {
        if (NOVELTY && NOVELTY.gene_mode != null) return !!+NOVELTY.gene_mode;
        if (NOVELTY && NOVELTY.samples) {
          for (const s of Object.values(NOVELTY.samples)) {
            const g = s && s.summary && s.summary.gene_mode;
            if (g != null && g !== "") return !!+g;
          }
        }
        return false;
      }
      // Per-backend display label + one-line method description for the Novelty tab.
      const _NOVELTY_METHODS = {
        mmseqs2: {
          label: "MMseqs2",
          short: "MMseqs",
          desc: "translated-search LCA (mmseqs taxonomy) on the de novo contigs — open-set, reports a best-hit %identity tail.",
        },
        kaiju: {
          label: "Kaiju",
          short: "Kaiju",
          desc: "greedy translated search (Kaiju) on the de novo contigs — open-set; no per-hit %identity, so the low-identity tail component is 0.",
        },
        bracken: {
          label: "Kraken2 + Bracken",
          short: "Bracken",
          desc: "k-mer classification (Kraken2) with Bracken abundance re-estimation on the de novo contigs — closed-set; counts are weighted per taxon.",
        },
      };
      function _noveltyMethodInfo() {
        const c = _noveltyClassifier();
        const info = _NOVELTY_METHODS[c] || { label: c ? c : "—", short: c ? c : "novelty", desc: "" };
        // In gene mode the query is Pyrodigal-predicted genes, not whole contigs — reflect that in
        // the method blurb so the description matches the count unit shown elsewhere on the tab.
        if (_novGeneMode() && info.desc) {
          return Object.assign({}, info, {
            desc: info.desc.replace(/the de novo contigs/g, "Pyrodigal-predicted genes (de novo contigs → ORFs)"),
          });
        }
        return info;
      }
      // Short, capitalized backend label used inline across the Novelty tab (column headers,
      // bucket names, tooltips) so everything tracks the active --novelty mode instead of
      // always reading "mmseqs". Falls back to "novelty" when no classifier is known.
      function _novClsShort() {
        return _noveltyMethodInfo().short || "novelty";
      }
      // Lowercase form for mid-sentence prose ("placed by kaiju", "X kaiju-only").
      function _novClsLc() {
        return _novClsShort().toLowerCase();
      }
      // The unit label for what the novelty classifier counted: "reads" for bracken (which
      // classifies reads directly via Kraken2 k-mers); otherwise "genes" by default (Pyrodigal-
      // predicted genes), or "contigs" under --disable_gene (kaiju/mmseqs2 on de novo contigs).
      function _novUnit() {
        const c = _noveltyClassifier();
        if (c === "bracken") return "reads";
        return _novGeneMode() ? "genes" : "contigs";
      }
      function _novUnitCap() {
        const u = _novUnit();
        return u.charAt(0).toUpperCase() + u.slice(1);
      }
      // Short abbreviation for tight UI spaces ("ctgs" / "genes" → "seqs", "reads" → "rds").
      function _novUnitShort() {
        const u = _novUnit();
        if (u === "contigs") return "ctgs";
        if (u === "genes") return "seqs";
        return "rds";
      }
      let CONTIG_DATA = BOOT.contig_data || []; // [{sample, organism, taxon_id, contigs, depth_histogram}]
      // Decode any base64-encoded breadth histograms baked into the BOOT data.
      // (Uploaded JSON entries are decoded later via _decodeBreadthHist in _extractContigs.)
      CONTIG_DATA.forEach(function (cd) {
        var bh = cd.breadth_histogram;
        if (bh && bh.b64 && !bh.bins) {
          try {
            var raw = atob(bh.b64);
            bh.bins = Array.from({ length: raw.length }, function (_, i) {
              return raw.charCodeAt(i);
            });
          } catch (e) {}
        }
      });
      // best_cutoffs.subkey.best_threshold — stored as 0–100
      const BEST_TASS_THRESH = ((BOOT.best_cutoffs || {}).subkey || {}).best_threshold;
      // ── Run metadata ───────────────────────────────────────────────────────
      const RUN_META = BOOT.run_metadata_records || []; // [{sample_name, run_id, latitude, longitude, ...}]
      const HAS_GEO = BOOT.has_geo || false;

      // ── Multi-value metadata normalization ──────────────────────────────────
      // Some metadata fields (e.g. host_disease) pack a comma-separated list of
      // values into one cell — "runny stool, cramps, cold-like symptoms". Split
      // those into arrays so the aggregation/comparison code (which already
      // handles arrays) counts each value separately instead of treating the
      // whole string as a single group. Idempotent: arrays and empties pass
      // through unchanged, so it is safe to run on boot data, uploaded records,
      // and re-loaded saved states alike.
      const _MULTI_VALUE_META_FIELDS = ["host_disease"];
      function _splitMultiMetaValue(v) {
        if (v == null) return v;
        if (Array.isArray(v)) return v.map((s) => String(s).trim()).filter(Boolean);
        const s = String(v).trim();
        if (s === "") return v;
        if (s.indexOf(",") === -1) return s; // single value — keep as a plain string
        return s
          .split(",")
          .map((p) => p.trim())
          .filter(Boolean);
      }
      function _normalizeMetaRecord(rec) {
        if (!rec || typeof rec !== "object") return rec;
        _MULTI_VALUE_META_FIELDS.forEach((f) => {
          if (f in rec) rec[f] = _splitMultiMetaValue(rec[f]);
        });
        return rec;
      }
      RUN_META.forEach(_normalizeMetaRecord);

      // ── User annotations ────────────────────────────────────────────────────
      // Captures metadata the user types directly into the report (rather than
      // via the init samplesheet or an uploaded CSV). Two surfaces:
      //   • rowCols / rowData  — custom columns + per-row values added to the
      //     main detection table (keyed by `${Specimen ID}||${Detected Organism}||${Taxonomic ID #}`).
      //   • metaCols           — custom columns added to the Run Metadata grid
      //     (the values themselves live on the RUN_META records so they export
      //     and persist with the rest of the metadata).
      // The whole object is serialized into the exported state (see ttExportState)
      // and restored on load, so edits survive a save → reload round-trip.
      const TT_ANNOT =
        BOOT.annotations && typeof BOOT.annotations === "object"
          ? {
              rowCols: BOOT.annotations.rowCols || [],
              rowData: BOOT.annotations.rowData || {},
              metaCols: BOOT.annotations.metaCols || [],
            }
          : { rowCols: [], rowData: {}, metaCols: [] };

      // Mirror stored row annotations back onto the live DATA rows so that the
      // existing sort / filter / DOM-based table export paths see them as plain
      // column values. Safe to call repeatedly (after every merge / redraw).
      function _ttApplyRowAnnotations() {
        if (!Array.isArray(DATA) || !TT_ANNOT.rowCols.length) return;
        DATA.forEach((r) => {
          // Same stable key as _tblRowKey() — inlined so this stays self-contained.
          const k = `${r["Specimen ID"]}||${r["Detected Organism"]}||${r["Taxonomic ID #"]}`;
          const a = TT_ANNOT.rowData[k];
          TT_ANNOT.rowCols.forEach((col) => {
            if (a && a[col] != null) r[col] = a[col];
            else if (r[col] == null) r[col] = "";
          });
        });
      }
