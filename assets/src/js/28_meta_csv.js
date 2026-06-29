      /* ═══════════════════════════════════════════════════════════════════════════
       -  §  META CSV
       -     Drag-and-drop metadata injection. Parses CSV/TSV with the
       -     columns: sample (or sample_name), run_id, latitude, longitude,
       -     depth, salinity, collection_time, location. Numeric fields are
       -     coerced. Matched rows merge into existing DATA records by sample.
═══════════════════════════════════════════════════════════════════════════ */
      (function () {
        const _META_NUM_FIELDS = new Set(["latitude", "longitude", "depth", "salinity"]);

        // ── Parse a CSV/TSV string → array of objects ────────────────────────
        function _parseCsvText(text, filename) {
          const sep = (filename || "").toLowerCase().endsWith(".tsv") ? "\t" : ",";
          const lines = text.split(/\r?\n/);
          if (!lines.length) return [];
          const headers = lines[0].split(sep).map((h) => h.trim().replace(/^"|"$/g, "").toLowerCase());
          const rows = [];
          for (let i = 1; i < lines.length; i++) {
            if (!lines[i].trim()) continue;
            const vals = lines[i].split(sep).map((v) => v.trim().replace(/^"|"$/g, ""));
            const row = {};
            headers.forEach((h, j) => {
              row[h] = vals[j] !== undefined ? vals[j] : "";
            });
            rows.push(row);
          }
          return rows;
        }

        // ── Convert a parsed CSV row → a RUN_META record (all columns kept) ──
        function _rowToMetaRecord(row) {
          // Accept "sample" or "sample_name" as the key column; normalize spaces → underscores
          const sampleName = (row["sample"] || row["sample_name"] || "").trim().replace(/\s+/g, "_");
          if (!sampleName) return null;
          const rec = { sample_name: sampleName };
          // Pass through EVERY column (except the sample key itself)
          Object.entries(row).forEach(([k, v]) => {
            if (k === "sample" || k === "sample_name") return;
            const vs = typeof v === "string" ? v.trim() : String(v ?? "").trim();
            if (vs === "" || vs.toLowerCase() === "null" || vs.toLowerCase() === "na") {
              rec[k] = null;
            } else if (_META_NUM_FIELDS.has(k)) {
              const n = parseFloat(vs);
              rec[k] = isNaN(n) ? null : n;
            } else {
              rec[k] = vs;
            }
          });
          return rec;
        }

        // ── Apply an array of RUN_META records into the live globals ─────────
        // fromBoot=true means we're restoring from BOOT (replaces everything)
        // fromBoot=false means user-uploaded CSV (merges / overrides by sample_name)
        window._applyMetaRecords = function (newRecords, fromBoot) {
          // Split comma-packed multi-value fields (e.g. host_disease) into arrays
          // so they are counted per value. Runs on every entry path (boot restore,
          // uploaded JSON, metadata CSV, re-loaded saved state).
          if (typeof _normalizeMetaRecord === "function") (newRecords || []).forEach(_normalizeMetaRecord);
          if (fromBoot) {
            // Full replace: restore to boot state
            RUN_META.length = 0;
            (newRecords || []).forEach((r) => RUN_META.push(r));
            // Restore SAMPLE_META run fields from boot
            Object.assign(SAMPLE_META, BOOT.sample_meta || {});
          } else {
            // Merge: upsert by sample_name
            newRecords.forEach((rec) => {
              if (!rec.sample_name) return;
              // Update RUN_META
              const existing = RUN_META.findIndex((r) => r.sample_name === rec.sample_name);
              if (existing >= 0) {
                Object.assign(RUN_META[existing], rec);
              } else {
                RUN_META.push(rec);
              }
              // Update SAMPLE_META — copy all keys from the meta record
              if (!SAMPLE_META[rec.sample_name]) {
                SAMPLE_META[rec.sample_name] = { sample_name: rec.sample_name };
              }
              Object.entries(rec).forEach(([k, v]) => {
                if (k !== "sample_name") SAMPLE_META[rec.sample_name][k] = v;
              });
            });
          }

          // Recompute HAS_GEO and update tab visibility
          const hasGeo = RUN_META.some((r) => r.latitude != null && r.longitude != null);
          const mapBtn = document.getElementById("map-tab-btn");
          const rmBtn = document.getElementById("runmeta-tab-btn");
          if (mapBtn) mapBtn.classList.toggle("hidden", !hasGeo);
          const _rmHasSamples = (DATA || []).some((r) => r["Specimen ID"]);
          if (rmBtn) rmBtn.classList.toggle("hidden", RUN_META.length === 0 && !_rmHasSamples);

          // Rebuild metadata table + update sub-tab enabled states
          _buildRunMetaTable();
          if (typeof _updateMetaSubTabStates === "function") _updateMetaSubTabStates();

          // Rebuild / refresh map if it was already initialized
          if (_leafletMap) {
            _markerLayer && _markerLayer.clearLayers();
            _markerObjects = {};
            _selectedSample = null;
            closeMapPanel();
            // Re-add markers
            const geoRows = RUN_META.filter((r) => r.latitude != null && r.longitude != null);
            const sampleNames = [...new Set(geoRows.map((r) => r.sample_name))];
            sampleNames.forEach((n, i) => {
              if (!sampleColors[n]) sampleColors[n] = PALETTE[i % PALETTE.length];
            });
            const bounds = [];
            const _lg = {};
            geoRows.forEach((rec) => {
              const lat = parseFloat(rec.latitude),
                lon = parseFloat(rec.longitude);
              if (isNaN(lat) || isNaN(lon)) return;
              const k = `${lat.toFixed(3)}_${lon.toFixed(3)}`;
              if (!_lg[k]) _lg[k] = { lat, lon, recs: [] };
              _lg[k].recs.push(rec);
            });
            Object.entries(_lg).forEach(([locKey, grp]) => {
              const { lat, lon, recs } = grp;
              const isSingle = recs.length === 1;
              const colors = recs.map((r) => sampleColors[r.sample_name] || "#1565C0");
              const mk = L.marker([lat, lon], { icon: isSingle ? _svgDot(colors[0], false) : _pieSvg(colors, false) });
              mk.on("click", () => {
                Object.values(_markerObjects).forEach((obj) => {
                  if (!obj.marker) return;
                  if (obj.recs && obj.recs.length > 1) {
                    obj.marker.setIcon(
                      _pieSvg(
                        obj.recs.map((r) => sampleColors[r.sample_name] || "#1565C0"),
                        false,
                      ),
                    );
                  } else {
                    obj.marker.setIcon(_svgDot(obj.color || "#1565C0", false));
                  }
                });
                _selectedSample = recs[0].sample_name;
                mk.setIcon(isSingle ? _svgDot(colors[0], true) : _pieSvg(colors, true));
                _renderMapGroupPanel(recs);
              });
              mk.addTo(_markerLayer);
              _markerObjects[locKey] = { marker: mk, recs, color: colors[0] };
              if (isSingle) _markerObjects[recs[0].sample_name] = _markerObjects[locKey];
              bounds.push([lat, lon]);
            });
            if (bounds.length > 0) {
              _leafletMap.fitBounds(bounds, { padding: [40, 40], maxZoom: 9 });
            }
            setTimeout(() => _leafletMap.invalidateSize(), 50);
          } else if (hasGeo) {
            // Map tab is now visible but not yet initialized — will init on tab click
          }
        };

        // ── Status helpers ────────────────────────────────────────────────────
        function _setMetaStatus(msg, isError) {
          const el = document.getElementById("meta-upload-status");
          if (el) {
            el.textContent = msg;
            el.style.color = isError ? "#b33" : "#2a7a2a";
          }
        }

        // ── Process a single meta file (CSV, TSV, or XLSX) ───────────────────
        function _processMetaFile(file) {
          return new Promise((resolve, reject) => {
            const name = file.name.toLowerCase();
            const reader = new FileReader();

            if (name.endsWith(".xlsx") || name.endsWith(".xls")) {
              // Parse XLSX — SheetJS is already loaded
              reader.onload = (e) => {
                try {
                  if (typeof XLSX === "undefined") throw new Error("SheetJS not loaded");
                  const wb = XLSX.read(new Uint8Array(e.target.result), { type: "array" });
                  // Use first sheet, or one named "Metadata" / "metadata" if present
                  const sheetName = wb.SheetNames.find((n) => /^metadata$/i.test(n)) || wb.SheetNames[0];
                  const ws = wb.Sheets[sheetName];
                  const rawRows = XLSX.utils.sheet_to_json(ws, { defval: "" });
                  // Normalise column names to lowercase
                  const rows = rawRows.map((r) => {
                    const out = {};
                    Object.entries(r).forEach(([k, v]) => {
                      out[k.toLowerCase().trim()] = String(v ?? "");
                    });
                    return out;
                  });
                  const records = rows.map(_rowToMetaRecord).filter(Boolean);
                  resolve(records);
                } catch (err) {
                  reject(`${file.name}: ${err.message}`);
                }
              };
              reader.onerror = () => reject(`${file.name}: could not read file`);
              reader.readAsArrayBuffer(file);
            } else {
              // CSV / TSV / TXT
              reader.onload = (e) => {
                try {
                  const rows = _parseCsvText(e.target.result, file.name);
                  const records = rows.map(_rowToMetaRecord).filter(Boolean);
                  resolve(records);
                } catch (err) {
                  reject(`${file.name}: ${err.message}`);
                }
              };
              reader.onerror = () => reject(`${file.name}: could not read file`);
              reader.readAsText(file);
            }
          });
        }

        // ── Public handlers ────────────────────────────────────────────────────
        window.handleMetaFiles = async function (files) {
          const fileArr = Array.from(files);
          if (!fileArr.length) return;
          _setMetaStatus(`Loading metadata…`, false);
          try {
            let allRecords = [];
            for (const f of fileArr) {
              const recs = await _processMetaFile(f);
              allRecords = allRecords.concat(recs);
            }
            if (!allRecords.length) {
              _setMetaStatus("No valid rows found in CSV.", true);
              return;
            }
            window._applyMetaRecords(allRecords, false);
            _setMetaStatus(`✓ Metadata loaded: ${allRecords.length} sample(s)`, false);
            const clearBtn = document.getElementById("meta-clear-btn");
            if (clearBtn) clearBtn.classList.remove("hidden");
          } catch (err) {
            _setMetaStatus(String(err), true);
          }
          const inp = document.getElementById("meta-upload-input");
          if (inp) inp.value = "";
        };

        window.handleMetaDrop = function (event) {
          event.preventDefault();
          document.getElementById("meta-upload-zone").classList.remove("drag-over");
          const files = event.dataTransfer.files;
          if (files.length) window.handleMetaFiles(files);
        };

        window.clearMetaData = function () {
          // Restore RUN_META to boot-time values
          window._applyMetaRecords(BOOT.run_metadata_records || [], true);
          _setMetaStatus("", false);
          const clearBtn = document.getElementById("meta-clear-btn");
          if (clearBtn) clearBtn.classList.add("hidden");
        };
      })();

      // ── Run Metadata analysis sub-tab state ──────────────────────────────
      let _activeMetaSub = null; // e.g. "longi" | "geo" | "host" | "cmp"

      // Heavy report build. Defined as a named function (was an immediately-invoked
      // IIFE) so the deferred scheduler below can run it AFTER the loading overlay
      // has painted, instead of blocking the first paint synchronously.
      function __ttRunInit() {
        // Banner subtitle
        document.getElementById("banner-sub").textContent = _buildBannerSub();

        // Auto-disable passes-threshold filter if it would show nothing
        const passCount = DATA.filter((r) => isTruthy(r["Passes Threshold"])).length;
        if (passCount === 0) {
          const fpEl = document.getElementById("filter-pass");
          if (fpEl) fpEl.checked = false;
        }

        // Sync the View-level dropdown to what the data actually contains.
        _syncViewLevelOptions();

        // Pre-populate min TASS filter from the recommended best_cutoffs threshold.
        // CRITICAL: also push the value into the paired range slider — otherwise the
        // slider knob stays at its HTML default (70) while the number reads e.g. 20,
        // and the very next drag snaps the number to wherever the knob is, which
        // looks like the slider is producing random values.
        {
          const minEl = document.getElementById("filter-min");
          const rangeEl = document.getElementById("filter-min-range");
          if (minEl) {
            const fromBestCutoff =
              BEST_TASS_THRESH != null && !isNaN(BEST_TASS_THRESH) ? Math.round(BEST_TASS_THRESH) : null;
            if (fromBestCutoff != null) {
              const clamped = Math.max(0, Math.min(100, fromBestCutoff));
              minEl.value = clamped;
              if (rangeEl) rangeEl.value = clamped;
            } else if (rangeEl) {
              // No best_cutoff override — still keep the two controls in sync with whatever
              // the number input was authored with in the HTML.
              rangeEl.value = minEl.value;
            }
          }
        }

        // Show proteins tab only if annotation data present
        const hasProtData =
          HAS_PROT &&
          [PROT.genus_summary, PROT.per_gene_hits, PROT.sample_overview, PROT.amr_genes, PROT.genus_by_property].some(
            (arr) => Array.isArray(arr) && arr.length > 0,
          );
        HAS_PROT = hasProtData;
        const protBtn = document.getElementById("prot-tab-btn");
        if (protBtn) {
          protBtn.classList.toggle("hidden", !hasProtData);
        }

        // Always show histogram tab (shows "no data" message gracefully when empty)
        const histBtn = document.getElementById("hist-tab-btn");
        if (histBtn) histBtn.classList.remove("hidden");

        HAS_NOVELTY = !!(NOVELTY && NOVELTY.samples && Object.keys(NOVELTY.samples).length);
        const noveltyBtn = document.getElementById("novelty-tab-btn");
        if (noveltyBtn) noveltyBtn.classList.toggle("hidden", !HAS_NOVELTY);

        // Show Map / Run Metadata tabs — always derived from live RUN_META, never from
        // the stale HAS_GEO constant, so real reports without metadata stay hidden.
        const mapBtn = document.getElementById("map-tab-btn");
        const runmetaBtn = document.getElementById("runmeta-tab-btn");
        const _initHasGeo = RUN_META.some((r) => r.latitude != null && r.longitude != null);
        if (mapBtn) mapBtn.classList.toggle("hidden", !_initHasGeo);
        // Run Metadata tab is available whenever there are samples (or metadata),
        // so users can add metadata directly in the report even when none was
        // supplied via the samplesheet / CSV.
        const _hasSamples = (DATA || []).some((r) => r["Specimen ID"]);
        if (runmetaBtn) runmetaBtn.classList.toggle("hidden", RUN_META.length === 0 && !_hasSamples);

        buildHmValueSel();
        buildSampleList();
        _computeBslLevels(); // inject BSL Level into DATA + ALL_COLS before table is built
        buildTable();
        renderTableHeaders();
        redraw();

        // ── Build Run Metadata table + initialize sub-tab states on load
        _buildRunMetaTable();
        if (typeof _updateMetaSubTabStates === "function") _updateMetaSubTabStates();
        // Auto-fill country / state / location from coordinates for any samples
        // that have lat/long but no geographic origin yet (deferred + no-op when
        // nothing needs filling, so reports without gaps make no network calls).
        if (typeof _ttAutofillGeoFromCoords === "function") {
          setTimeout(() => _ttAutofillGeoFromCoords(), 600);
        }
        _initExportEnhancer();
        const pdfBtn = document.getElementById("report-pdf-btn");
        if (pdfBtn) {
          pdfBtn.addEventListener("click", _openReportPdfModal);
          // Hover tooltip — same pattern as info-tip-btn icons
          pdfBtn.addEventListener("mouseenter", (ev) => showTip(pdfBtn.dataset.tip, ev));
          pdfBtn.addEventListener("mousemove", moveTip);
          pdfBtn.addEventListener("mouseleave", hideTip);
        }

        // ── Tooltip toggle button ──────────────────────────────────────────
        const _ttBtn = document.getElementById("tooltip-toggle-btn");
        const _ttLabel = document.getElementById("tooltip-toggle-label");
        if (_ttBtn) {
          const _updateTtBtn = () => {
            _ttBtn.title = tooltipsEnabled ? "Tooltips on — click to disable" : "Tooltips off — click to enable";
            _ttBtn.style.color = tooltipsEnabled ? "#1565c0" : "#999";
            _ttBtn.style.borderColor = tooltipsEnabled ? "#ccd6e8" : "#ddd";
            _ttBtn.querySelector("i").className = tooltipsEnabled ? "fas fa-comment-dots" : "fas fa-comment-slash";
            if (_ttLabel) {
              _ttLabel.textContent = tooltipsEnabled ? "Tooltips on" : "Tooltips off";
              _ttLabel.style.color = tooltipsEnabled ? "inherit" : "#aaa";
            }
          };
          _updateTtBtn();
          _ttBtn.addEventListener("click", () => {
            tooltipsEnabled = !tooltipsEnabled;
            if (!tooltipsEnabled) hideTip();
            _updateTtBtn();
          });
        }

        // ── Wire up sortable map-panel column headers ──────────────────────
        _initPanelSortHeaders();

        // ── Wire sidebar legend collapse toggle ────────────────────────────
        const _legendToggle = document.getElementById("sidebar-legend-toggle");
        const _legendBody = document.getElementById("sidebar-legend-body");
        if (_legendToggle && _legendBody) {
          _legendToggle.addEventListener("click", () => {
            const open = _legendBody.style.display !== "none";
            _legendBody.style.display = open ? "none" : "";
            _legendToggle.innerHTML = open ? "&#9660;" : "&#9650;";
            _legendToggle.title = open ? "Expand legend" : "Collapse legend";
          });
        }

        // Also update legend when the global TASS slider changes
        const _fmin = document.getElementById("filter-min");
        if (_fmin)
          _fmin.addEventListener("input", () => {
            if (typeof _updateSidebarLegend === "function") _updateSidebarLegend();
          });
        // Also update when a sample color picker changes
        document.getElementById("sample-list") &&
          document.getElementById("sample-list").addEventListener("input", (e) => {
            if (e.target && e.target.type === "color" && typeof _updateSidebarLegend === "function")
              _updateSidebarLegend();
          });
      }

      // ── Deferred init scheduler ───────────────────────────────────────────────
      // Paint the loading overlay first, then run the heavy report build on the
      // next frame, then fade the overlay out. Deferring this work (instead of
      // running it synchronously during parse) lets the browser show the loading
      // screen immediately rather than freezing on a blank page.
      (function __ttScheduleInit() {
        function __ttHideLoader() {
          var ov = document.getElementById("tt-loading-overlay");
          if (!ov) return;
          ov.classList.add("tt-hide");
          setTimeout(function () {
            if (ov && ov.parentNode) ov.parentNode.removeChild(ov);
          }, 450);
        }
        function run() {
          // Double rAF guarantees the overlay has been painted before the
          // (potentially long) synchronous build begins.
          requestAnimationFrame(function () {
            requestAnimationFrame(function () {
              try {
                __ttRunInit();
              } catch (e) {
                console.error("[taxtriage] report init failed:", e);
                var msg = document.getElementById("tt-loading-msg");
                if (msg) msg.textContent = "Something went wrong while building the report. See console for details.";
              } finally {
                __ttHideLoader();
              }
            });
          });
        }
        if (document.readyState === "loading") {
          document.addEventListener("DOMContentLoaded", run, { once: true });
        } else {
          run();
        }
        // Safety net: never leave the overlay stuck if something stalls.
        window.addEventListener("load", function () {
          setTimeout(__ttHideLoader, 8000);
        });
      })();

