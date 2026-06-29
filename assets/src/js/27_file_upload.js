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

          function _collectStrain(strain) {
            rows.push(_flattenOrg(strain, sampleName, sampleType, totalReads));
            const cd = _extractContigs(strain, sampleName);
            if (cd) contigs.push(cd);
            for (const ann of strain.protein_annotations || []) {
              protHits.push({ "Specimen ID": sampleName, ...ann });
            }
          }

          for (const grp of data.organisms || []) {
            for (const sk_m of grp.members || []) {
              for (const strain of sk_m.members || []) {
                _collectStrain(strain);
              }
              if (!(sk_m.members && sk_m.members.length)) {
                _collectStrain(sk_m);
              }
            }
            if (!(grp.members && grp.members.length)) {
              _collectStrain(grp);
            }
          }

          const prot_data = protHits.length ? { per_gene_hits: _applyProtColRemap(protHits) } : {};

          return { rows, contigs, prot_data, sampleName, meta };
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
              if (r["Detected Organism"])
                r["Detected Organism"] = String(r["Detected Organism"]).replace(/°/g, "").trim();
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

        const _PROT_KEYS = ["genus_summary", "per_gene_hits", "sample_overview", "amr_genes"];

        function _mergeAndRedraw() {
          const hasAnyData =
            _uploadedRows.length > 0 || _PROT_KEYS.some((k) => (_uploadedProtData[k] || []).length > 0);
          if (!hasAnyData) return;

          const combined = [...(BOOT.records || []), ..._uploadedRows];
          DATA = combined;

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

          // Merge run metadata from uploaded JSON files (only fields from the JSON)
          if (_uploadedMetaRecords.length && window._applyMetaRecords) {
            window._applyMetaRecords(_uploadedMetaRecords, false);
          }

          // Merge protein annotation data into PROT
          _PROT_KEYS.forEach((k) => {
            PROT[k] = [...((BOOT.prot_data || {})[k] || []), ...(_uploadedProtData[k] || [])];
          });

          // Update proteins tab visibility
          HAS_PROT = _PROT_KEYS.some((k) => Array.isArray(PROT[k]) && PROT[k].length > 0);
          const protBtn = document.getElementById("prot-tab-btn");
          if (protBtn) protBtn.classList.toggle("hidden", !HAS_PROT);

          // Recompute schema from merged data
          const { cols, numeric } = _deriveSchema(combined);
          // Prefer boot-time cols list (keeps column order), union any new cols.
          // _NON_DISPLAY_COLS are analysis-only fields (ANI list/flag) that must
          // not surface as detections-table columns.
          const knownCols = new Set(BOOT.all_cols || []);
          const extraCols = cols.filter((c) => !knownCols.has(c) && !_NON_DISPLAY_COLS.has(c));
          ALL_COLS = [...(BOOT.all_cols || cols).filter((c) => !_NON_DISPLAY_COLS.has(c)), ...extraCols];
          NUMERIC = numeric;

          // Auto-disable passes-threshold if nothing passes
          const passCount = DATA.filter((r) => isTruthy(r["Passes Threshold"])).length;
          const fpEl = document.getElementById("filter-pass");
          if (fpEl && passCount === 0) fpEl.checked = false;

          // Re-sync view-level options. Reset the synthesis guard first so
          // _synthesizeHierarchy re-runs on the combined data — uploaded rows
          // from _parseJson arrive without a Level field and need proxy rows.
          _HIERARCHY_SYNTHESIZED = false;
          if (typeof _syncViewLevelOptions === "function") _syncViewLevelOptions();

          // Recompute BSL levels with merged data
          _computeBslLevels();

          // Rebuild UI
          const samples = uniq(DATA.map((r) => r["Specimen ID"] || "")).filter(Boolean);
          document.getElementById("banner-sub").textContent = _buildBannerSub();

          buildSampleList();
          buildTable();
          renderTableHeaders();
          buildHmValueSel();
          redraw();

          // Show histogram tab (always)
          const histBtn = document.getElementById("hist-tab-btn");
          if (histBtn) histBtn.classList.remove("hidden");

          // Clear histogram selectors so they rebuild
          if (window._resetHistSelectors) window._resetHistSelectors();
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
                  // ── Combined all.samples.json: iterate over every sample ─────
                  if (data.taxtriage_combined === true && Array.isArray(data.samples)) {
                    const merged = { rows: [], contigs: [], prot_data: {}, metas: [] };
                    const _PROT_KEYS_LOCAL = ["genus_summary", "per_gene_hits", "sample_overview", "amr_genes"];
                    _PROT_KEYS_LOCAL.forEach((k) => {
                      merged.prot_data[k] = [];
                    });
                    for (const sample of data.samples) {
                      const r = _parseJson(sample, file.name);
                      merged.rows.push(...(r.rows || []));
                      merged.contigs.push(...(r.contigs || []));
                      if (r.prot_data) {
                        _PROT_KEYS_LOCAL.forEach((k) => {
                          if (r.prot_data[k] && r.prot_data[k].length) merged.prot_data[k].push(...r.prot_data[k]);
                        });
                      }
                      if (r.sampleName && r.meta) {
                        merged.metas.push({ sampleName: r.sampleName, meta: r.meta });
                      }
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
              _uploadedRows.push(...(result.rows || []));
              _uploadedContigs.push(...(result.contigs || []));
              if (result.prot_data) {
                _PROT_KEYS.forEach((k) => {
                  if (result.prot_data[k] && result.prot_data[k].length)
                    _uploadedProtData[k].push(...result.prot_data[k]);
                });
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
              summary.textContent = `Loaded ${_uploadedRows.length} row(s) from ${_uploadedNames.length} file(s)`;
            }
            _statEl.appendChild(summary);
          }

          const hasUpload = _uploadedRows.length > 0 || _PROT_KEYS.some((k) => (_uploadedProtData[k] || []).length > 0);
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

