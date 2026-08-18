/* ═══════════════════════════════════════════════════════════════════════════
       -  §  FILE UPLOAD
       -     Client-side ingestion of .paths.json / .xlsx / .tsv files via
       -     drag-and-drop or the file picker. Mirrors the Python
       -     _flatten_organism logic so fresh uploads go through the same
       -     normalization that make_report.py applies on the server.
       -     End result: appends rows to DATA, then calls redraw().
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  /* ── flatten one organism entry (mirrors Python _flatten_organism) ─────── */
  function _flattenOrg(org, sampleName, sampleType, totalReads) {
    const strain_reads = parseFloat(org.numreads || 0) || 0;
    const pct = (strain_reads / Math.max(1, totalReads)) * 100;
    const covered = parseInt(org.covered_bases || 0) || 0;
    const genome_len = parseInt(org.length || 0) || 0;
    const breadth_pct = genome_len > 0 ? Math.round((covered / genome_len) * 10000) / 100 : 0;
    const tass = parseFloat(org.tass_score || 0) || 0;
    const tax = org.taxonomy || {};
    return {
      "Specimen ID": sampleName,
      "Sample Type": sampleType || "unknown",
      "Detected Organism": org.name || "Unknown",
      "TASS Score": Math.round(tass * 10) / 10, // 4th col
      "Taxonomic ID #": String(org.key || ""),
      Subkey: String(org.subkey || org.key || ""),
      "Microbial Category": org.microbial_category || "Unknown",
      "Ann Class": org.annClass || "",
      IsAnnotated: org.is_annotated === "Yes" ? "Yes" : "No",
      "High Consequence": !!org.high_cons,
      "Mol Type": (org.mol_type || "").trim().toLowerCase(),
      Status: org.status || "",
      "# Reads Aligned": Math.round(strain_reads),
      "% Reads": Math.round(pct * 10000) / 10000,
      Coverage: Math.round((parseFloat(org.coverage || 0) || 0) * 1000) / 10,
      "Covered Bases": covered,
      "Genome Length (bp)": genome_len,
      "Breadth %": breadth_pct,
      "Mean Depth": Math.round((parseFloat(org.meandepth || 0) || 0) * 100) / 100,
      "Gini Coefficient": Math.round((parseFloat(org.gini_coefficient || 0) || 0) * 1000) / 1000,
      "Mean MapQ": Math.round((parseFloat(org.meanmapq || 0) || 0) * 10) / 10,
      "Mean BaseQ": Math.round((parseFloat(org.meanbaseq || 0) || 0) * 10) / 10,
      "Minhash Score": Math.round((parseFloat(org.minhash_reduction || 0) || 0) * 1000) / 1000,
      "Breadth Score": Math.round((parseFloat(org.breadth_log_score || 0) || 0) * 1000) / 1000,
      "MapQ Score": Math.round((parseFloat(org.mapq_score || 0) || 0) * 1000) / 1000,
      "Disparity Score": Math.round((parseFloat(org.disparity || 0) || 0) * 1000) / 1000,
      "Diamond Identity": Math.round((parseFloat(org.diamond_identity || 0) || 0) * 10) / 10,
      "K2 Reads": parseInt(org.k2_reads || 0) || 0,
      RPM: Math.round((parseFloat(org.rpm || 0) || 0) * 100) / 100,
      RPKM: Math.round((parseFloat(org.rpkm || 0) || 0) * 10000) / 10000,
      "Passes Threshold": !!org.passes_threshold,
      // ANI annotation (mirrors Python _flatten_organism). Presence of the
      // 'high_ani_matches' key — even empty — means ANI was computed for this
      // sample; absence ⇒ the upload predates ANI support (unsupported).
      "ANI Annotated": Object.prototype.hasOwnProperty.call(org, "high_ani_matches"),
      "High ANI Matches": (Array.isArray(org.high_ani_matches) ? org.high_ani_matches : []).map((m) => ({
        key: String((m && m.key) || ""),
        ani_pct: Math.round((parseFloat((m && m.ani_pct) || 0) || 0) * 100) / 100,
      })),
      Domain: tax.domain || "",
      Superkingdom: tax.superkingdom || "",
      Kingdom: tax.kingdom || "",
      Phylum: tax.phylum || "",
      Class: tax["class"] || "",
      Order: tax.order || "",
      Family: tax.family || "",
      Genus: tax.genus || "",
    };
  }

  /* ── decode a base64-encoded Uint8 breadth histogram ───────────────────── */
  function _decodeBreadthHist(bh) {
    if (!bh) return null;
    if (Array.isArray(bh.bins)) return bh; // already decoded (legacy plain array)
    if (!bh.b64) return null;
    try {
      const raw = atob(bh.b64);
      const bins = Array.from({ length: raw.length }, (_, i) => raw.charCodeAt(i));
      return { bin_size: bh.bin_size, total_len: bh.total_len, bins, breaks: bh.breaks || [0] };
    } catch (e) {
      return null;
    }
  }

  /* ── extract contig data from a strain entry ────────────────────────────── */
  function _extractContigs(strain, sampleName) {
    if (!strain.contigs && !strain.depth_histogram && !strain.breadth_histogram) return null;
    return {
      sample: sampleName,
      organism: strain.name || "Unknown",
      taxon_id: String(strain.key || ""),
      contigs: strain.contigs || [],
      depth_histogram: strain.depth_histogram || {},
      breadth_histogram: _decodeBreadthHist(strain.breadth_histogram),
    };
  }

  /* ── parse a .paths.json file ────────────────────────────────────────────── */
  function _parseJson(data, filename) {
    const rows = [];
    const contigs = [];
    const protHits = [];
    const meta = data.metadata || {};
    const sampleName = meta.sample_name || filename.replace(/\.paths\.json$|\.json$/, "");
    const sampleType = meta.sample_type || "unknown";
    const totalReads = parseInt(meta.total_reads || 1) || 1;

    function _collectAnnotations(annotations, organism, topLevel, taxid, level) {
      for (const ann of annotations || []) {
        protHits.push({
          "Specimen ID": sampleName,
          _organism_name: organism || "Unknown",
          _toplevel_name: topLevel || organism || "Unknown",
          _taxid: String(taxid || ""),
          _level: level,
          ...ann,
        });
      }
    }

    function _collectStrain(strain, topLevel) {
      rows.push(_flattenOrg(strain, sampleName, sampleType, totalReads));
      const cd = _extractContigs(strain, sampleName);
      if (cd) contigs.push(cd);
      _collectAnnotations(strain.protein_annotations, strain.name, topLevel, strain.key, "strain");
    }

    for (const grp of data.organisms || []) {
      const topLevel = grp.toplevelname || grp.name || "Unknown";
      _collectAnnotations(grp.protein_annotations_genus, topLevel, topLevel, grp.toplevelkey || grp.key, "genus");
      for (const sk_m of grp.members || []) {
        const speciesName = sk_m.subkeyname || sk_m.name || topLevel;
        _collectAnnotations(sk_m.protein_annotations, speciesName, topLevel, sk_m.subkey || sk_m.key, "species");
        for (const strain of sk_m.members || []) {
          _collectStrain(strain, topLevel);
        }
        if (!(sk_m.members && sk_m.members.length)) {
          rows.push(_flattenOrg(sk_m, sampleName, sampleType, totalReads));
          const cd = _extractContigs(sk_m, sampleName);
          if (cd) contigs.push(cd);
        }
      }
      if (!(grp.members && grp.members.length)) {
        rows.push(_flattenOrg(grp, sampleName, sampleType, totalReads));
        const cd = _extractContigs(grp, sampleName);
        if (cd) contigs.push(cd);
      }
    }

    const prot_data = protHits.length ? { per_gene_hits: _applyProtColRemap(protHits) } : {};

    return { rows, contigs, prot_data, sampleName, meta };
  }

  /* ── novelty JSON (NOVELTY_COLLECT output) ───────────────────────────────
           Two shapes ship from the pipeline and neither carries an `organisms`
           array, so _parseJson() returned zero rows for them and the upload was
           silently discarded — no sample in the sidebar, no banner change, no
           Novelty tab:
             all.novelty.json  → {taxtriage_novelty, classifier, gene_mode,
                                  samples: {<sample>: {summary, candidates}}}
             <sample>.novelty.json → {summary: {sample, …}, candidates: [...]}
           Returns null when `data` is not a novelty payload. */
  function _looksLikeNoveltySample(v) {
    return !!v && typeof v === "object" && (v.summary != null || Array.isArray(v.candidates));
  }

  function _parseNoveltyJson(data, filename) {
    if (!data || typeof data !== "object") return null;

    // ── combined all.novelty.json ────────────────────────────────────────
    const bag = data.samples;
    const isCombined =
      data.taxtriage_novelty === true ||
      (bag && typeof bag === "object" && !Array.isArray(bag) && Object.values(bag).some(_looksLikeNoveltySample));
    if (isCombined && bag && typeof bag === "object" && !Array.isArray(bag)) {
      const samples = {};
      Object.keys(bag).forEach((s) => {
        if (s && _looksLikeNoveltySample(bag[s])) samples[s] = bag[s];
      });
      if (!Object.keys(samples).length) return null;
      return {
        novelty: {
          samples,
          classifier: data.classifier || "",
          gene_mode: data.gene_mode != null ? data.gene_mode : null,
        },
      };
    }

    // ── per-sample <sample>.novelty.json ─────────────────────────────────
    if (_looksLikeNoveltySample(data) && !data.organisms) {
      const summary = data.summary || {};
      const sampleName =
        summary.sample ||
        String(filename || "")
          .replace(/\.novelty\.json$/i, "")
          .replace(/\.json$/i, "");
      if (!sampleName) return null;
      return {
        novelty: {
          samples: { [sampleName]: { summary, candidates: data.candidates || [] } },
          classifier: summary.classifier || "",
          gene_mode: summary.gene_mode != null ? summary.gene_mode : null,
        },
      };
    }
    return null;
  }

  /* ── standalone annotate_report.tsv / .xlsx (annotate_report.py output) ──
           These carry VF *and* AMR hits in one flat table. Without routing they
           fell through to the detections parser, which produced junk rows and
           never lit up the VF/AMR tab. Mirrors _is_amr_annotation() +
           load_standalone_annotations() in bin/make_report.py. */
  const _ANNOT_GENE_KEYS = ["gene_name", "Gene", "Gene Name", "gene"];
  const _ANNOT_SIGNAL_KEYS = [
    "property",
    "Property",
    "source",
    "Source",
    "pident",
    "%id",
    "antibiotics",
    "Antibiotics",
    "antibiotics_class",
    "Antibiotics Class",
    "source_id",
    "Source ID",
  ];

  function _looksLikeAnnotateReport(rows, filename) {
    if (/\.annotate_report\b/i.test(String(filename || ""))) return true;
    if (!rows || !rows.length) return false;
    const keys = Object.keys(rows[0]);
    // A detections export always has this; an annotate report never does.
    if (keys.includes("Detected Organism") || keys.includes("TASS Score")) return false;
    const hasGene = _ANNOT_GENE_KEYS.some((k) => keys.includes(k));
    const hasSignal = _ANNOT_SIGNAL_KEYS.some((k) => keys.includes(k));
    return hasGene && hasSignal;
  }

  function _isAmrAnnotation(prop, abxClass, abx, source) {
    const s = (v) => (v == null ? "" : String(v)).trim().toLowerCase();
    const p = s(prop);
    if (["resist", "amr", "antibiotic", "antimicrobial"].some((w) => p.includes(w))) return true;
    if (s(abxClass) || s(abx)) return true;
    if (["card", "ncbiamr", "resfinder", "argannot", "amrfinder"].some((w) => s(source).includes(w))) return true;
    return false;
  }

  /* Turn annotate-report rows into the {per_gene_hits, amr_genes} shape PROT
           expects. `fallbackSample` recovers "<sample>.annotate_report.tsv" names
           for files whose rows carry no sample column. */
  function _parseAnnotateReport(rows, filename) {
    const fallbackSample = String(filename || "")
      .replace(/\.[^.]*$/, "")
      .split(".annotate_report")[0];
    const pick = (r, ...keys) => {
      for (const k of keys) {
        const v = r[k];
        if (v !== undefined && v !== null && v !== "") return v;
      }
      return null;
    };
    const per_gene_hits = [];
    const amr_genes = [];
    (rows || []).forEach((r) => {
      const prop = pick(r, "property", "Property");
      const abxClass = pick(r, "antibiotics_class", "Antibiotics Class");
      const abx = pick(r, "antibiotics", "Antibiotics");
      const source = pick(r, "source", "Source");
      const out = {
        "Specimen ID": pick(r, "Specimen ID", "sample", "sample_name", "Sample", "Sample ID") || fallbackSample,
        Genus: pick(r, "genus", "Genus") || "Unknown",
        Species: pick(r, "species", "Species"),
        Gene: pick(r, ..._ANNOT_GENE_KEYS),
        Product: pick(r, "product", "Product"),
        Property: prop,
        Description: pick(r, "function", "product", "Description"),
        "Antibiotics Class": abxClass,
        Antibiotics: abx,
        Source: source,
        "Source ID": pick(r, "source_id", "Source ID"),
        "%id": pick(r, "pident", "%id"),
        "E-value": pick(r, "evalue", "E-value"),
        Bitscore: pick(r, "bitscore", "Bitscore"),
        "Reference Organism": pick(r, "organism", "Reference Organism"),
        Level: pick(r, "level", "Level"),
      };
      if (_isAmrAnnotation(prop, abxClass, abx, source)) amr_genes.push(out);
      else per_gene_hits.push(out);
    });
    return { per_gene_hits, amr_genes };
  }

  /* ── remap alternative-format column names to canonical names ──────────── */
  // Some TXT/TSV exports (e.g. TASS Challenge reports) use slightly different
  // column headers.  This map converts those names to the canonical ones the
  // visualisations expect so every tab renders correctly.
  const _ALT_COL_MAP = {
    AnnClass: "Ann Class",
    "Protein Identity Score": "Diamond Identity",
    "Sample ID": "Specimen ID",
    SampleID: "Specimen ID",
    "Sample Name": "Specimen ID",
    Sample: "Specimen ID",
  };

  function _normalizeAltFormatRows(rows) {
    if (!rows.length) return rows;
    const firstKeys = Object.keys(rows[0]);
    const needsRemap = firstKeys.some((k) => Object.prototype.hasOwnProperty.call(_ALT_COL_MAP, k));
    // Also detect when "Breadth %" is absent but "Coverage" is present —
    // in that case treat "Coverage" as the breadth percentage.
    const hasBreadth = firstKeys.includes("Breadth %");
    const hasCoverage = firstKeys.includes("Coverage");
    const fillBreadth = !hasBreadth && hasCoverage;
    if (!needsRemap && !fillBreadth) return rows;
    return rows.map((row) => {
      const out = {};
      for (const [k, v] of Object.entries(row)) {
        out[_ALT_COL_MAP[k] !== undefined ? _ALT_COL_MAP[k] : k] = v;
      }
      ["Coverage", "Breadth %", "% Reads"].forEach((k) => {
        if (typeof out[k] === "string" && out[k].includes("%")) {
          const n = parseFloat(out[k].replace(/%/g, ""));
          out[k] = isNaN(n) ? out[k] : n;
        }
      });
      if (out["Coverage"] === undefined && out["Mean Coverage"] !== undefined) {
        out["Coverage"] = out["Mean Coverage"];
      }
      if (out["% Reads"] === undefined && out["% Reads Aligned"] !== undefined) {
        out["% Reads"] = out["% Reads Aligned"];
      }
      if (fillBreadth && out["Breadth %"] === undefined) {
        out["Breadth %"] = out["Coverage"];
      }
      if ((out["Genus"] === undefined || out["Genus"] === "") && out["Detected Organism"]) {
        const org = String(out["Detected Organism"]).trim();
        const first = org.split(/\s+/)[0] || "";
        out["Genus"] = first;
      }
      return out;
    });
  }

  /* ── parse an XLSX / TSV file using SheetJS ─────────────────────────────── */
  function _parseXlsx(ab, filename) {
    const wb = XLSX.read(ab, { type: "array" });

    // Known protein-annotation sheet names and their prot_data keys
    const PROT_SHEET_MAP = {
      "Genus Summary": "genus_summary",
      "Per-Gene Hits": "per_gene_hits",
      "Sample Overview": "sample_overview",
      "AMR Genes": "amr_genes",
    };
    const PROT_SHEET_NAMES = new Set(Object.keys(PROT_SHEET_MAP));

    const HEADER_KEYS = ["Detected Organism", "Specimen ID", "TASS Score", "# Reads Aligned"];
    function _scoreRows(rows) {
      if (!rows.length) return 0;
      let score = 0;
      HEADER_KEYS.forEach((k) => {
        if (rows[0][k] !== undefined) score += 1;
      });
      return score;
    }
    function _getRowsWithFallback(ws) {
      const rows0 = XLSX.utils.sheet_to_json(ws, { defval: null });
      const rows1 = XLSX.utils.sheet_to_json(ws, { defval: null, range: 1 });
      return _scoreRows(rows1) > _scoreRows(rows0) ? rows1 : rows0;
    }

    // ── Find organism sheet (prefer header match, else "Organisms", else first non-prot sheet) ──
    let bestSheetName = null;
    let bestRows = [];
    let bestScore = -1;
    wb.SheetNames.forEach((name) => {
      const ws = wb.Sheets[name];
      const rows = _getRowsWithFallback(ws);
      const score = _scoreRows(rows);
      if (score > bestScore) {
        bestScore = score;
        bestRows = rows;
        bestSheetName = name;
      }
    });
    const orgSheetName =
      bestScore > 0
        ? bestSheetName
        : wb.SheetNames.find((n) => n === "Organisms") ||
          wb.SheetNames.find((n) => !PROT_SHEET_NAMES.has(n)) ||
          wb.SheetNames[0];

    const ws = wb.Sheets[orgSheetName];
    const rows = bestScore > 0 && bestSheetName === orgSheetName ? bestRows : _getRowsWithFallback(ws);
    if (rows.length && rows[0]["Detected Organism"] !== undefined) {
      // Fix degree symbol in organism names
      rows.forEach((r) => {
        if (r["Detected Organism"]) r["Detected Organism"] = String(r["Detected Organism"]).replace(/°/g, "").trim();
        if (r["High Consequence"] !== undefined) r["High Consequence"] = isTruthy(r["High Consequence"]);
        if (r["Passes Threshold"] !== undefined) r["Passes Threshold"] = isTruthy(r["Passes Threshold"]);
      });
    }

    // ── Extract protein annotation sheets ────────────────────────────────
    const prot_data = {};
    for (const [sheetName, key] of Object.entries(PROT_SHEET_MAP)) {
      if (wb.Sheets[sheetName]) {
        const sheetRows = XLSX.utils.sheet_to_json(wb.Sheets[sheetName], { defval: null });
        prot_data[key] = _applyProtColRemap(sheetRows);
      }
    }

    // No organism sheet was recognised, but the workbook looks like a
    // standalone annotate_report.xlsx — route it to VF/AMR instead of letting
    // its rows masquerade as detections.
    if (bestScore <= 0 && _looksLikeAnnotateReport(rows, filename)) {
      const ann = _parseAnnotateReport(rows, filename);
      _PROT_KEYS.forEach((k) => {
        if (ann[k] && ann[k].length) prot_data[k] = [...(prot_data[k] || []), ...ann[k]];
      });
      return { rows: [], contigs: [], prot_data };
    }

    return { rows: _normalizeAltFormatRows(rows), contigs: [], prot_data };
  }

  /* ── derive ALL_COLS and NUMERIC from rows ──────────────────────────────── */
  function _deriveSchema(rows) {
    if (!rows.length) return { cols: [], numeric: new Set() };
    const cols = Object.keys(rows[0]);
    const numeric = new Set();
    for (const col of cols) {
      const sample = rows
        .slice(0, 50)
        .map((r) => r[col])
        .filter((v) => v !== null && v !== undefined && v !== "");
      if (
        sample.length &&
        sample.every((v) => typeof v === "number" || (typeof v === "string" && !isNaN(parseFloat(v))))
      ) {
        numeric.add(col);
      }
    }
    return { cols, numeric };
  }

  /* ── merge uploaded data into globals ───────────────────────────────────── */
  let _uploadedRows = [];
  let _uploadedContigs = [];
  let _uploadedProtData = { genus_summary: [], per_gene_hits: [], sample_overview: [], amr_genes: [] };
  let _uploadedNames = [];
  let _uploadedMetaRecords = []; // run metadata extracted from uploaded .paths.json files
  let _uploadedNovelty = { samples: {}, classifier: "", gene_mode: null };

  const _PROT_KEYS = ["genus_summary", "per_gene_hits", "sample_overview", "amr_genes"];

  function _uploadedNoveltyCount() {
    return Object.keys(_uploadedNovelty.samples || {}).length;
  }

  /* Fold uploaded novelty into the live payload. NOVELTY / NOVELTY_DL are const
           bindings from early.js, so they are mutated in place; BOOT.novelty is
           updated alongside them so clearUploadedData() can restore correctly. */
  function _mergeNovelty() {
    if (!_uploadedNoveltyCount()) return;
    if (NOVELTY && (!NOVELTY.samples || typeof NOVELTY.samples !== "object")) NOVELTY.samples = {};
    if (BOOT) {
      if (!BOOT.novelty || typeof BOOT.novelty !== "object") BOOT.novelty = { samples: {} };
      if (!BOOT.novelty.samples) BOOT.novelty.samples = {};
    }
    Object.keys(_uploadedNovelty.samples).forEach((s) => {
      const payload = _uploadedNovelty.samples[s];
      if (NOVELTY && NOVELTY.samples) NOVELTY.samples[s] = payload;
      if (BOOT && BOOT.novelty && BOOT.novelty.samples && BOOT.novelty.samples !== (NOVELTY || {}).samples) {
        BOOT.novelty.samples[s] = payload;
      }
    });
    if (_uploadedNovelty.classifier) {
      if (NOVELTY && !NOVELTY.classifier) NOVELTY.classifier = _uploadedNovelty.classifier;
      if (BOOT && BOOT.novelty && !BOOT.novelty.classifier) BOOT.novelty.classifier = _uploadedNovelty.classifier;
    }
    if (_uploadedNovelty.gene_mode != null) {
      if (NOVELTY && NOVELTY.gene_mode == null) NOVELTY.gene_mode = _uploadedNovelty.gene_mode;
      if (BOOT && BOOT.novelty && BOOT.novelty.gene_mode == null) BOOT.novelty.gene_mode = _uploadedNovelty.gene_mode;
    }
    HAS_NOVELTY = !!(NOVELTY && NOVELTY.samples && Object.keys(NOVELTY.samples).length);
    if (BOOT) BOOT.has_novelty = HAS_NOVELTY;
    const novBtn = document.getElementById("novelty-tab-btn");
    if (novBtn) novBtn.classList.toggle("hidden", !HAS_NOVELTY);
  }

  function _refreshUploadedUi() {
    const safe = (label, fn) => {
      try {
        fn();
      } catch (e) {
        if (window.console && console.error) console.error(`upload refresh: ${label} failed`, e);
      }
    };

    safe("filter cache", () => {
      if (typeof _invalidateFilterCache === "function") _invalidateFilterCache();
    });
    safe("banner", () => {
      const banner = document.getElementById("banner-sub");
      if (banner) banner.textContent = _buildBannerSub();
    });
    safe("sample sidebar", () => buildSampleList());
    safe("table", () => {
      buildTable();
      renderTableHeaders();
    });
    safe("heatmap selector", () => buildHmValueSel());
    safe("active tab", () => redraw());

    // Optional tabs are lazy-rendered, but uploaded auxiliary data must also
    // invalidate any content they rendered before the upload. Draw both now so
    // their sample selectors/charts are current on the first click.
    if (HAS_PROT && typeof drawProteins === "function") safe("VF/AMR tab", () => drawProteins());
    if (HAS_NOVELTY && typeof drawNovelty === "function") safe("novelty tab", () => drawNovelty());
  }

  /* Rebuild the UI exactly once after an upload.
           The double-rAF exists so the upload-status paint lands first, but rAF
           is PAUSED while the tab is hidden or backgrounded — and dropping a file
           very often means the user just came from another window. When that
           happened the callback never fired, so the sidebar, the banner and the
           Summary tab kept showing the pre-upload state until some unrelated
           interaction (changing the TASS filter, switching tabs) forced a
           redraw. A timeout fallback races the rAF and whichever wins runs the
           refresh; _uiRefreshPending makes sure it only happens once. */
  let _uiRefreshPending = false;
  function _scheduleUploadedUiRefresh() {
    if (_uiRefreshPending) return;
    _uiRefreshPending = true;
    const run = () => {
      if (!_uiRefreshPending) return;
      _uiRefreshPending = false;
      _refreshUploadedUi();
    };
    if (typeof requestAnimationFrame === "function") {
      requestAnimationFrame(() => requestAnimationFrame(run));
    }
    // Fallback: fires if rAF is throttled (hidden tab) or unavailable. Timers
    // are throttled in background tabs too but, unlike rAF, still run.
    setTimeout(run, 150);
  }

  function _mergeAndRedraw() {
    // Novelty-only and VF/AMR-only uploads carry no detection rows. Counting
    // just rows here meant those files were parsed and then thrown away: the
    // sidebar, the banner totals and the Novelty / VF-AMR tabs never updated.
    const hasAnyData =
      _uploadedRows.length > 0 ||
      _PROT_KEYS.some((k) => (_uploadedProtData[k] || []).length > 0) ||
      _uploadedNoveltyCount() > 0;
    if (!hasAnyData) return;

    /* Every merge step below is individually guarded, the same way
             _emptyAllData() guards its teardown steps. Before this, the whole
             pipeline was one straight run of statements: if ANY step threw —
             most likely right after "Remove All", which leaves SAMPLE_META,
             SPECIMEN_OVERRIDE, TT_ANNOT and friends in an atypical, fully-empty
             shape that a helper like _applyMetaRecords()/_syncViewLevelOptions()
             may not expect — the function aborted before ever reaching
             _scheduleUploadedUiRefresh(). The upload-status line (built from
             _uploadedRows/_uploadedProtData/_uploadedNovelty, all populated
             BEFORE this function runs) still reported success, but
             buildSampleList()/buildTable()/drawProteins()/drawNovelty() never
             ran — the exact "loaded successfully but the UI never refreshed"
             symptom. Nothing here depends on a later step succeeding, so a
             failure in one must not skip the rest. */
    const step = (label, fn) => {
      try {
        fn();
      } catch (e) {
        if (window.console && console.error) console.error("_mergeAndRedraw: " + label + " failed", e);
      }
    };

    step("novelty merge", () => _mergeNovelty());

    const combined = [...(BOOT.records || []), ..._uploadedRows];
    DATA = combined;

    step("contig data", () => {
      // Merge contig data; _uploadedContigs are already decoded via _extractContigs
      CONTIG_DATA = [...(BOOT.contig_data || []), ..._uploadedContigs];
      CONTIG_DATA.forEach(function (cd) {
        var bh = cd.breadth_histogram;
        if (bh && bh.b64 && !bh.bins) {
          try {
            var r = atob(bh.b64);
            bh.bins = Array.from({ length: r.length }, function (_, i) {
              return r.charCodeAt(i);
            });
          } catch (e) {}
        }
      });
      if (typeof _invalidateSummaryHistMap === "function") _invalidateSummaryHistMap();
    });

    step("run metadata", () => {
      // Merge run metadata from uploaded JSON files (only fields from the JSON)
      if (_uploadedMetaRecords.length && window._applyMetaRecords) {
        window._applyMetaRecords(_uploadedMetaRecords, false);
      }
    });

    step("protein annotations", () => {
      // Merge protein annotation data into PROT
      _PROT_KEYS.forEach((k) => {
        PROT[k] = [...((BOOT.prot_data || {})[k] || []), ...(_uploadedProtData[k] || [])];
      });
      // New dataset — re-arm the "VF/AMR categories only" default so it is applied
      // against the merged properties rather than staying pinned to the old set.
      if (typeof _resetProtHiddenDefaults === "function") _resetProtHiddenDefaults();

      // Update proteins tab visibility
      HAS_PROT = hasProtAnnotations(PROT);
      const protBtn = document.getElementById("prot-tab-btn");
      if (protBtn) protBtn.classList.toggle("hidden", !HAS_PROT);
    });

    step("schema", () => {
      // Recompute schema from merged data
      const { cols, numeric } = _deriveSchema(combined);
      // Prefer boot-time cols list (keeps column order), union any new cols.
      // _NON_DISPLAY_COLS are analysis-only fields (ANI list/flag) that must
      // not surface as detections-table columns.
      const knownCols = new Set(BOOT.all_cols || []);
      const extraCols = cols.filter((c) => !knownCols.has(c) && !_NON_DISPLAY_COLS.has(c));
      ALL_COLS = [...(BOOT.all_cols || cols).filter((c) => !_NON_DISPLAY_COLS.has(c)), ...extraCols];
      NUMERIC = numeric;
    });

    step("passes-threshold checkbox", () => {
      // Auto-disable passes-threshold if nothing passes
      const passCount = DATA.filter((r) => isTruthy(r["Passes Threshold"])).length;
      const fpEl = document.getElementById("filter-pass");
      if (fpEl && passCount === 0) fpEl.checked = false;
    });

    step("view-level sync", () => {
      // Re-sync view-level options. Reset the synthesis guard first so
      // _synthesizeHierarchy re-runs on the combined data — uploaded rows
      // from _parseJson arrive without a Level field and need proxy rows.
      _HIERARCHY_SYNTHESIZED = false;
      if (typeof _syncViewLevelOptions === "function") _syncViewLevelOptions();
    });

    step("BSL levels", () => _computeBslLevels());

    // Invalidate immediately, then rebuild after the browser has committed the
    // upload-status paint. Each UI surface is guarded independently so a chart
    // error cannot prevent the banner or sidebar from showing the new samples.
    // This step MUST run even if everything above failed — it's what actually
    // reflects the merged data (however partial) in the sidebar and tables.
    step("schedule UI refresh", () => {
      if (typeof _invalidateFilterCache === "function") _invalidateFilterCache();
      _scheduleUploadedUiRefresh();
    });

    step("histogram tab", () => {
      // Show histogram tab (always)
      const histBtn = document.getElementById("hist-tab-btn");
      if (histBtn) histBtn.classList.remove("hidden");

      // Clear histogram selectors so they rebuild.
      if (window._resetHistSelectors) window._resetHistSelectors();
    });
  }

  /* ── status display helpers ─────────────────────────────────────────────── */
  function _setStatus(msg, isError) {
    const el = document.getElementById("upload-status");
    if (el) {
      el.textContent = msg;
      el.style.color = isError ? "#b33" : "#2a7a2a";
    }
  }

  /* ── process a single File object ──────────────────────────────────────── */
  function _processFile(file) {
    return new Promise((resolve, reject) => {
      const name = file.name.toLowerCase();
      const reader = new FileReader();

      if (name.endsWith(".json")) {
        reader.onload = (e) => {
          try {
            const data = JSON.parse(e.target.result);
            // ── Novelty payload (all.novelty.json / <sample>.novelty.json) ──
            // Checked before the detections parsers: it has no `organisms`
            // array, so those would return an empty result and drop the file.
            const nov = _parseNoveltyJson(data, file.name);
            if (nov) {
              resolve(nov);
              return;
            }
            // ── Combined all.samples.json: iterate over every sample ─────
            if (data.taxtriage_combined === true && Array.isArray(data.samples)) {
              const merged = { rows: [], contigs: [], prot_data: {}, metas: [] };
              const _PROT_KEYS_LOCAL = ["genus_summary", "per_gene_hits", "sample_overview", "amr_genes"];
              _PROT_KEYS_LOCAL.forEach((k) => {
                merged.prot_data[k] = [];
              });
              for (const sample of data.samples) {
                const r = _parseJson(sample, file.name);
                _ttPushAll(merged.rows, r.rows || []);
                _ttPushAll(merged.contigs, r.contigs || []);
                if (r.prot_data) {
                  _PROT_KEYS_LOCAL.forEach((k) => {
                    if (r.prot_data[k] && r.prot_data[k].length) _ttPushAll(merged.prot_data[k], r.prot_data[k]);
                  });
                }
                if (r.sampleName && r.meta) {
                  merged.metas.push({ sampleName: r.sampleName, meta: r.meta });
                }
              }
              // ── Blocks embedded by combine_samples_json.py (v1.1+) ────
              // all.odr.json can now carry the novelty payload and the
              // standalone VF/AMR annotations, so one drag lights up the
              // Novelty and VF/AMR tabs too. Older combined files simply
              // lack these keys.
              const embeddedNov = _parseNoveltyJson({ samples: (data.novelty || {}).samples }, file.name);
              if (embeddedNov) {
                merged.novelty = {
                  samples: embeddedNov.novelty.samples,
                  classifier: data.novelty.classifier || "",
                  gene_mode: data.novelty.gene_mode != null ? data.novelty.gene_mode : null,
                };
              }
              if (data.prot_data && typeof data.prot_data === "object") {
                _PROT_KEYS_LOCAL.forEach((k) => {
                  const rows = data.prot_data[k];
                  if (Array.isArray(rows) && rows.length) _ttPushAll(merged.prot_data[k], rows);
                });
              }
              resolve(merged);
            } else {
              resolve(_parseJson(data, file.name));
            }
          } catch (err) {
            reject(`${file.name}: JSON parse error – ${err.message}`);
          }
        };
        reader.readAsText(file);
      } else if (name.endsWith(".xlsx") || name.endsWith(".xls")) {
        reader.onload = (e) => {
          try {
            const ab = new Uint8Array(e.target.result);
            const result = _parseXlsx(ab, file.name);
            resolve(result);
          } catch (err) {
            reject(`${file.name}: XLSX parse error – ${err.message}`);
          }
        };
        reader.readAsArrayBuffer(file);
      } else if (name.endsWith(".tsv") || name.endsWith(".txt") || name.endsWith(".csv")) {
        reader.onload = (e) => {
          try {
            const sep = name.endsWith(".csv") ? "," : "\t";
            const lines = e.target.result.split(/\r\n|\n|\r/);
            const headers = lines[0].split(sep).map((h) => h.trim().replace(/^"|"$/g, ""));
            if (headers.length && headers[0].charAt(0) === "\ufeff") headers[0] = headers[0].slice(1);
            const rows = [];
            for (let i = 1; i < lines.length; i++) {
              if (!lines[i].trim()) continue;
              const vals = lines[i].split(sep).map((v) => v.trim().replace(/^"|"$/g, ""));
              const row = {};
              headers.forEach((h, j) => {
                row[h] = vals[j] !== undefined ? vals[j] : null;
              });
              rows.push(row);
            }
            // A standalone annotate_report is VF/AMR annotation, not detections.
            if (_looksLikeAnnotateReport(rows, file.name)) {
              resolve({ rows: [], contigs: [], prot_data: _parseAnnotateReport(rows, file.name) });
              return;
            }
            resolve({ rows: _normalizeAltFormatRows(rows), contigs: [] });
          } catch (err) {
            reject(`${file.name}: parse error – ${err.message}`);
          }
        };
        reader.readAsText(file);
      } else {
        reject(`${file.name}: unsupported file type`);
      }
    });
  }

  /* ── public: handle file input change ──────────────────────────────────── */
  window.handleUploadFiles = async function (files) {
    const fileArr = Array.from(files);
    if (!fileArr.length) return;

    // ── Per-file progress list ──────────────────────────────────────────
    // Show a spinning bacteria logo next to each file while it is being
    // read/parsed; swap to a check (or error) icon when it finishes.
    const _statEl = document.getElementById("upload-status");
    const _uplRows = {};
    if (_statEl) {
      _statEl.style.color = "#555";
      _statEl.textContent = "";
      fileArr.forEach((file) => {
        const row = document.createElement("div");
        row.style.cssText = "display:flex;align-items:center;gap:.45em;margin:.18em 0";
        const ic = document.createElement("i");
        ic.className = "fas fa-bacteria fa-spin";
        ic.style.cssText = "color:#1565c0;flex:0 0 auto";
        const nm = document.createElement("span");
        nm.textContent = file.name;
        nm.style.cssText = "flex:1 1 auto;overflow:hidden;text-overflow:ellipsis;white-space:nowrap";
        const st = document.createElement("span");
        st.textContent = "uploading…";
        st.style.cssText = "flex:0 0 auto;color:#888;font-size:.92em";
        row.appendChild(ic);
        row.appendChild(nm);
        row.appendChild(st);
        _statEl.appendChild(row);
        _uplRows[file.name] = { icon: ic, state: st };
      });
    }
    // Yield a frame so the spinners actually paint before the synchronous
    // FileReader / JSON.parse work starts.
    await new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r)));

    const _markRow = (name, ok, label) => {
      const r = _uplRows[name];
      if (!r) return;
      r.icon.className = ok ? "fas fa-check-circle" : "fas fa-triangle-exclamation";
      r.icon.style.color = ok ? "#2a7a2a" : "#b33";
      r.state.textContent = label;
      r.state.style.color = ok ? "#2a7a2a" : "#b33";
    };

    const errors = [];
    for (const file of fileArr) {
      try {
        const result = await _processFile(file);
        _ttPushAll(_uploadedRows, result.rows || []);
        _ttPushAll(_uploadedContigs, result.contigs || []);
        if (result.prot_data) {
          _PROT_KEYS.forEach((k) => {
            if (result.prot_data[k] && result.prot_data[k].length)
              _ttPushAll(_uploadedProtData[k], result.prot_data[k]);
          });
        }
        // Protein-only annotation files have no detection rows or metadata
        // object. Register their sample ids explicitly so the sidebar and
        // banner use the same authoritative sample set as the VF/AMR tab.
        const auxiliarySamples = new Set();
        _PROT_KEYS.forEach((k) => {
          ((result.prot_data || {})[k] || []).forEach((row) => {
            const sample = row && (row["Specimen ID"] || row.Sample || row.sample);
            if (sample) auxiliarySamples.add(String(sample));
          });
        });
        auxiliarySamples.forEach((sample) => {
          if (!SAMPLE_META[sample]) SAMPLE_META[sample] = { sample_name: sample };
          if (BOOT && BOOT.sample_meta) BOOT.sample_meta[sample] = SAMPLE_META[sample];
        });
        if (auxiliarySamples.size && typeof _noteSampleMetaChanged === "function") _noteSampleMetaChanged();
        // ── Novelty payload ────────────────────────────────────────────
        // Registering each novelty sample in SAMPLE_META is what puts it in
        // the sidebar (_allSampleIds treats SAMPLE_META as authoritative) and
        // what feeds the banner's read counters, which read total_reads /
        // total_organism_reads from there.
        if (result.novelty && result.novelty.samples) {
          const nv = result.novelty;
          Object.keys(nv.samples).forEach((sn) => {
            _uploadedNovelty.samples[sn] = nv.samples[sn];
            const summary = (nv.samples[sn] || {}).summary || {};
            const existing = SAMPLE_META[sn] || {};
            const totalReads = parseFloat(summary.total_reads);
            // ref_aligned is the novelty summary's count of reads that hit a
            // reference. Only used when the sample has no detections-derived
            // figure, so a real .paths.json upload always wins.
            const refAligned = parseFloat(summary.ref_aligned);
            SAMPLE_META[sn] = Object.assign({}, existing, {
              sample_name: sn,
              total_reads: existing.total_reads || (isNaN(totalReads) ? 0 : totalReads),
              total_organism_reads: existing.total_organism_reads || (isNaN(refAligned) ? 0 : refAligned),
            });
            if (BOOT && BOOT.sample_meta) BOOT.sample_meta[sn] = SAMPLE_META[sn];
          });
          if (nv.classifier && !_uploadedNovelty.classifier) _uploadedNovelty.classifier = nv.classifier;
          if (nv.gene_mode != null && _uploadedNovelty.gene_mode == null) _uploadedNovelty.gene_mode = nv.gene_mode;
          if (typeof _noteSampleMetaChanged === "function") _noteSampleMetaChanged();
        }
        // Extract run metadata fields — handles both single (.paths.json)
        // and combined (all.samples.json) results.
        const _metaEntries = result.metas
          ? result.metas // combined: array of {sampleName, meta}
          : result.meta && result.sampleName
          ? [{ sampleName: result.sampleName, meta: result.meta }]
          : [];
        const _META_NUM = new Set(["latitude", "longitude", "depth", "salinity"]);
        for (const { sampleName: sn, meta: m } of _metaEntries) {
          const metaRec = { sample_name: sn };
          ["run_id", "latitude", "longitude", "depth", "salinity", "collection_time", "location"].forEach((k) => {
            const v = m[k];
            if (v === undefined || v === null || v === "") {
              metaRec[k] = null;
            } else if (_META_NUM.has(k)) {
              const n = parseFloat(v);
              metaRec[k] = isNaN(n) ? null : n;
            } else {
              metaRec[k] = String(v);
            }
          });
          if (Object.values(metaRec).some((v, i) => i > 0 && v != null)) {
            _uploadedMetaRecords.push(metaRec);
          }
          // Always update SAMPLE_META with read-count fields so the
          // banner % reads classified stays accurate after a JSON upload.
          const existing = SAMPLE_META[sn] || {};
          SAMPLE_META[sn] = Object.assign({}, existing, {
            sample_name: sn,
            total_reads: parseFloat(m.total_reads) || existing.total_reads || 0,
            aligned_reads: parseFloat(m.aligned_reads) || existing.aligned_reads || 0,
            total_organism_reads: parseFloat(m.total_organism_reads) || existing.total_organism_reads || 0,
          });
          if (BOOT && BOOT.sample_meta) {
            BOOT.sample_meta[sn] = SAMPLE_META[sn];
          }
          if (typeof _noteSampleMetaChanged === "function") _noteSampleMetaChanged();
        }
        _uploadedNames.push(file.name);
        _markRow(file.name, true, "loaded");
      } catch (err) {
        errors.push(err);
        _markRow(file.name, false, "failed");
      }
    }

    // Append a summary line beneath the per-file rows (don't overwrite them).
    if (_statEl) {
      const summary = document.createElement("div");
      summary.style.cssText = "margin-top:.4em;font-weight:600";
      if (errors.length) {
        summary.style.color = "#b33";
        summary.textContent = "Errors: " + errors.join("; ");
      } else {
        summary.style.color = "#2a7a2a";
        // Break the count out per feed. "Loaded N row(s)" alone could not
        // distinguish "this file has no VF/AMR block" from "the VF/AMR block
        // was ingested but the tab failed to render" — the two look identical
        // from the outside and need opposite fixes.
        const _protCount = _PROT_KEYS.reduce((n, k) => n + (_uploadedProtData[k] || []).length, 0);
        const _novCount = _uploadedNoveltyCount();
        const parts = [`${_uploadedRows.length} detection row(s)`];
        parts.push(_novCount ? `${_novCount} novelty sample(s)` : "no novelty data");
        parts.push(_protCount ? `${_protCount} VF/AMR row(s)` : "no VF/AMR data");
        summary.textContent = `Loaded from ${_uploadedNames.length} file(s): ` + parts.join(" · ");
      }
      _statEl.appendChild(summary);
    }

    const hasUpload =
      _uploadedRows.length > 0 ||
      _PROT_KEYS.some((k) => (_uploadedProtData[k] || []).length > 0) ||
      _uploadedNoveltyCount() > 0;
    if (hasUpload) {
      const clearBtn = document.getElementById("upload-clear-btn");
      if (clearBtn) clearBtn.classList.remove("hidden");
      _mergeAndRedraw();
    }

    // Reset file input so same file can be re-loaded
    const inp = document.getElementById("file-upload-input");
    if (inp) inp.value = "";
  };

  /* ── public: drag-drop handler ──────────────────────────────────────────── */
  window.handleUploadDrop = function (event) {
    event.preventDefault();
    document.getElementById("upload-zone").classList.remove("drag-over");
    const files = event.dataTransfer.files;
    if (files.length) window.handleUploadFiles(files);
  };

  /* ── public: attach a file to a specific sample row ──────────────────────
           Triggered by the per-sample paperclip button in the sidebar. Marks the
           row as loading (spinning bacteria icon) BEFORE the synchronous parse so
           the spinner actually paints, then reuses the normal upload pipeline so
           the new JSON / XLSX / TSV merges into the dataset. */
  window.addFileToSample = async function (sampleId, fileList) {
    const files = Array.from(fileList || []);
    if (!files.length) return;
    _loadingSampleIds.add(sampleId);
    if (typeof buildSampleList === "function") buildSampleList(); // paint spinner on this row
    // Yield two frames so the spinner is visible before the blocking parse.
    await new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r)));
    try {
      await window.handleUploadFiles(files);
    } finally {
      _loadingSampleIds.delete(sampleId);
      if (typeof buildSampleList === "function") buildSampleList();
    }
  };

  /* ── public: clear ONLY the upload buffers (no UI rebuild).
           Called from removeAllSamples() so that the next file drop won't
           re-merge previously-uploaded rows on top of a freshly-emptied
           dataset. */
  window._clearUploadBuffers = function () {
    _uploadedRows = [];
    _uploadedContigs = [];
    _uploadedProtData = { genus_summary: [], per_gene_hits: [], sample_overview: [], amr_genes: [] };
    _uploadedNames = [];
    _uploadedMetaRecords = [];
    _uploadedNovelty = { samples: {}, classifier: "", gene_mode: null };
    const clearBtn = document.getElementById("upload-clear-btn");
    if (clearBtn) clearBtn.classList.add("hidden");
    const stat = document.getElementById("upload-status");
    if (stat) stat.textContent = "";
  };

  /* ── public: prune upload buffers for a single sample id.
           Called from _deleteSampleData(id) so the next file drop won't
           resurrect a sample the user explicitly deleted. */
  window._pruneUploadBuffersForSample = function (id) {
    _uploadedRows = _uploadedRows.filter((r) => r["Specimen ID"] !== id);
    _uploadedContigs = _uploadedContigs.filter((c) => c.sample !== id);
    _PROT_KEYS.forEach((k) => {
      _uploadedProtData[k] = (_uploadedProtData[k] || []).filter((r) => r["Specimen ID"] !== id);
    });
    _uploadedMetaRecords = _uploadedMetaRecords.filter((m) => m.sample_name !== id);
    // Novelty-only samples are resurrected the same way — the buffer must be
    // pruned too, or a later drop merges the deleted sample's novelty payload
    // back in (see _mergeNovelty()).
    if (_uploadedNovelty && _uploadedNovelty.samples) delete _uploadedNovelty.samples[id];
  };

  /* ── public: rename a sample id inside the upload buffers.
           Called from _renameSample so a subsequent file drop (which
           re-runs _mergeAndRedraw) doesn't resurrect the OLD name from
           the accumulated buffers. */
  window._renameInUploadBuffers = function (oldName, newName) {
    if (!oldName || !newName || oldName === newName) return;
    _uploadedRows.forEach((r) => {
      if (r["Specimen ID"] === oldName) r["Specimen ID"] = newName;
    });
    _uploadedContigs.forEach((c) => {
      if (c.sample === oldName) c.sample = newName;
    });
    _PROT_KEYS.forEach((k) => {
      (_uploadedProtData[k] || []).forEach((r) => {
        if (r["Specimen ID"] === oldName) r["Specimen ID"] = newName;
      });
    });
    _uploadedMetaRecords.forEach((m) => {
      if (m.sample_name === oldName) m.sample_name = newName;
    });
  };

  /* ── public: clear uploaded data ────────────────────────────────────────── */
  window.clearUploadedData = function () {
    _uploadedRows = [];
    _uploadedContigs = [];
    _uploadedProtData = { genus_summary: [], per_gene_hits: [], sample_overview: [], amr_genes: [] };
    _uploadedNames = [];
    _uploadedMetaRecords = [];
    _uploadedNovelty = { samples: {}, classifier: "", gene_mode: null };
    DATA = BOOT.records || [];
    ALL_COLS = (BOOT.all_cols || []).filter((c) => !_NON_DISPLAY_COLS.has(c));
    NUMERIC = new Set(BOOT.numeric_cols || []);
    CONTIG_DATA = BOOT.contig_data || [];
    CONTIG_DATA.forEach(function (cd) {
      var bh = cd.breadth_histogram;
      if (bh && bh.b64 && !bh.bins) {
        try {
          var r = atob(bh.b64);
          bh.bins = Array.from({ length: r.length }, function (_, i) {
            return r.charCodeAt(i);
          });
        } catch (e) {}
      }
    });
    if (typeof _invalidateSummaryHistMap === "function") _invalidateSummaryHistMap();
    // Restore PROT to boot-time data
    _PROT_KEYS.forEach((k) => {
      PROT[k] = (BOOT.prot_data || {})[k] || [];
    });
    HAS_PROT = BOOT.has_prot || false;
    const protBtnClear = document.getElementById("prot-tab-btn");
    if (protBtnClear) protBtnClear.classList.toggle("hidden", !HAS_PROT);

    const passCount = DATA.filter((r) => isTruthy(r["Passes Threshold"])).length;
    const fpEl = document.getElementById("filter-pass");
    if (fpEl) fpEl.checked = passCount > 0;

    const samples = uniq(DATA.map((r) => r["Specimen ID"] || "")).filter(Boolean);
    document.getElementById("banner-sub").textContent = _buildBannerSub();

    const clearBtn = document.getElementById("upload-clear-btn");
    if (clearBtn) clearBtn.classList.add("hidden");
    document.getElementById("upload-status").textContent = "";
    _setStatus("", false);

    // Recompute BSL levels from restored BOOT data
    _computeBslLevels();

    buildSampleList();
    buildTable();
    renderTableHeaders();
    buildHmValueSel();
    if (window._resetHistSelectors) window._resetHistSelectors();
    redraw();
  };
})();
