      /* ═══════════════════════════════════════════════════════════════════════════
       -  §  SESSION STATE  (export / import) + GITHUB-PAGES STARTUP DIALOG
       -     ttExportState()      — serialize current data + all UI filters → JSON
       -     ttLoadStateFiles()   — read a state JSON file from <input type=file>
       -     ttLoadState(state)   — apply a parsed state envelope (data + filters)
       -     _ttApplyBoot(boot)   — replace all live data globals from a BOOT-shape
       -     _ttCaptureUi/_ttRestoreUi — snapshot & restore every id'd form control
       -     The GitHub-Pages build (which inlines its dataset as <script
       -     id="BOOTSTRAP">) shows a startup dialog: Demo (default) / Empty /
       -     Upload. Pipeline reports (window.HEATMAP_BOOT) skip it entirely.
═══════════════════════════════════════════════════════════════════════════ */
      (function () {
        function _ttSetStateStatus(msg, isError) {
          var el = document.getElementById("state-status");
          if (el) {
            el.textContent = msg;
            el.style.color = isError ? "#b33" : "#2a7a2a";
          }
        }

        // Decode any base64-packed breadth histograms in-place (same as boot).
        function _ttDecodeContigB64(list) {
          (list || []).forEach(function (cd) {
            var bh = cd && cd.breadth_histogram;
            if (bh && bh.b64 && !bh.bins) {
              try {
                var raw = atob(bh.b64);
                bh.bins = Array.from({ length: raw.length }, function (_, i) {
                  return raw.charCodeAt(i);
                });
              } catch (e) {}
            }
          });
        }

        function _ttEmptyBoot() {
          return {
            records: [],
            all_cols: [],
            numeric_cols: [],
            sample_meta: {},
            prot_data: {},
            has_prot: false,
            contig_data: [],
            best_cutoffs: {},
            run_metadata_records: [],
            has_geo: false,
          };
        }

        // Replace every live data global from a BOOT-shaped object, then rebuild
        // the whole UI by reusing the existing clear/restore code paths.
        function _ttApplyBoot(b) {
          b = b || {};
          // Rewrite the BOOT baseline so the well-tested clearUploadedData() /
          // _applyMetaRecords() restore paths rebuild from the *new* data rather
          // than resurrecting the original report.
          BOOT.records = Array.isArray(b.records) ? b.records : [];
          BOOT.all_cols = Array.isArray(b.all_cols) ? b.all_cols : [];
          BOOT.numeric_cols = Array.isArray(b.numeric_cols) ? b.numeric_cols : [];
          BOOT.sample_meta = b.sample_meta || {};
          BOOT.prot_data = b.prot_data || {};
          BOOT.has_prot = !!b.has_prot;
          BOOT.contig_data = Array.isArray(b.contig_data) ? b.contig_data : [];
          BOOT.best_cutoffs = b.best_cutoffs || {};
          BOOT.run_metadata_records = Array.isArray(b.run_metadata_records) ? b.run_metadata_records : [];
          BOOT.has_geo = !!b.has_geo;
          _ttDecodeContigB64(BOOT.contig_data);

          // Restore user annotations (detection-row notes + custom metadata cols).
          var _ann = b.annotations && typeof b.annotations === "object" ? b.annotations : {};
          TT_ANNOT.rowCols = Array.isArray(_ann.rowCols) ? _ann.rowCols : [];
          TT_ANNOT.rowData = _ann.rowData && typeof _ann.rowData === "object" ? _ann.rowData : {};
          TT_ANNOT.metaCols = Array.isArray(_ann.metaCols) ? _ann.metaCols : [];

          // Full-replace SAMPLE_META *up front* — before clearUploadedData()'s
          // redraw — so the summary KPIs (% Classified, Total Reads, which read
          // total_reads / total_organism_reads from SAMPLE_META) compute against
          // the new metadata instead of an empty object (which showed "N/A").
          // _applyMetaRecords(...,true) re-applies BOOT.sample_meta afterward
          // (idempotent).
          SAMPLE_META = Object.assign({}, b.sample_meta || {});

          // Re-synthesize the taxonomy hierarchy against the new records.
          try {
            _HIERARCHY_SYNTHESIZED = false;
          } catch (e) {}

          // Detections side: rebuild table/heatmap/etc. from the new BOOT.
          if (window.clearUploadedData) window.clearUploadedData();
          // Metadata side: replace RUN_META (also splits host_disease), refresh
          // map + run-meta tab + sub-tab states.
          if (window._applyMetaRecords) window._applyMetaRecords(BOOT.run_metadata_records || [], true);

          if (typeof _syncViewLevelOptions === "function") _syncViewLevelOptions();
          var histBtn = document.getElementById("hist-tab-btn");
          if (histBtn) histBtn.classList.remove("hidden");
        }
        window._ttApplyBoot = _ttApplyBoot;

        // ── Snapshot every id'd form control + sample colors + active tab ──────
        function _ttCaptureUi() {
          var ui = { controls: {}, sampleColors: {}, activeTab: null };
          document.querySelectorAll("input[id], select[id], textarea[id]").forEach(function (el) {
            if (el.type === "file") return;
            if (el.type === "checkbox" || el.type === "radio") ui.controls[el.id] = { checked: el.checked };
            else ui.controls[el.id] = { value: el.value };
          });
          try {
            if (typeof sampleColors !== "undefined" && sampleColors) ui.sampleColors = Object.assign({}, sampleColors);
          } catch (e) {}
          try {
            if (typeof activeTab !== "undefined") ui.activeTab = activeTab;
          } catch (e) {}
          return ui;
        }

        function _ttRestoreUi(ui) {
          if (!ui) return;
          // Sample colors first so the redraw below uses them.
          try {
            if (ui.sampleColors && typeof sampleColors !== "undefined" && sampleColors) {
              Object.keys(ui.sampleColors).forEach(function (k) {
                sampleColors[k] = ui.sampleColors[k];
              });
            }
          } catch (e) {}
          Object.keys(ui.controls || {}).forEach(function (id) {
            var el = document.getElementById(id);
            if (!el || el.type === "file") return;
            var v = ui.controls[id];
            if (v && Object.prototype.hasOwnProperty.call(v, "checked")) el.checked = v.checked;
            else if (v && v.value !== undefined) el.value = v.value;
          });
          // Keep the TASS min number<->range slider pair in sync.
          var minEl = document.getElementById("filter-min");
          var rangeEl = document.getElementById("filter-min-range");
          if (minEl && rangeEl) rangeEl.value = minEl.value;
          // Re-apply filters / rerender.
          try {
            if (typeof buildTable === "function") buildTable();
          } catch (e) {}
          try {
            if (typeof renderTableHeaders === "function") renderTableHeaders();
          } catch (e) {}
          try {
            if (typeof redraw === "function") redraw();
          } catch (e) {}
          try {
            if (typeof _updateMetaSubTabStates === "function") _updateMetaSubTabStates();
          } catch (e) {}
          // Restore the active tab (only if visible + enabled).
          if (ui.activeTab) {
            var btn = document.querySelector('.tab-btn[data-tab="' + ui.activeTab + '"]');
            if (btn && !btn.classList.contains("hidden") && !btn.classList.contains("tab-disabled")) btn.click();
          }
        }

        // ── Export: full data snapshot + UI filter state → downloadable JSON ──
        window.ttExportState = function () {
          try {
            var boot = {
              records: Array.isArray(DATA) ? DATA : [],
              all_cols: Array.isArray(ALL_COLS) ? ALL_COLS : [],
              numeric_cols: typeof NUMERIC !== "undefined" && NUMERIC ? Array.from(NUMERIC) : [],
              sample_meta: typeof SAMPLE_META !== "undefined" && SAMPLE_META ? SAMPLE_META : {},
              prot_data: typeof PROT !== "undefined" && PROT ? PROT : {},
              has_prot: !!(typeof HAS_PROT !== "undefined" && HAS_PROT),
              contig_data: Array.isArray(CONTIG_DATA) ? CONTIG_DATA : [],
              best_cutoffs: (BOOT && BOOT.best_cutoffs) || {},
              run_metadata_records: Array.isArray(RUN_META) ? RUN_META : [],
              has_geo:
                Array.isArray(RUN_META) &&
                RUN_META.some(function (r) {
                  return r.latitude != null && r.longitude != null;
                }),
              // User annotations typed into the report (detection-row notes +
              // custom metadata columns). Metadata cell *values* already ride
              // along on run_metadata_records above.
              annotations: {
                rowCols: TT_ANNOT.rowCols.slice(),
                rowData: TT_ANNOT.rowData,
                metaCols: TT_ANNOT.metaCols.slice(),
              },
            };
            var state = {
              __taxtriage_state__: true,
              version: 1,
              exported_at: new Date().toISOString(),
              boot: boot,
              ui: _ttCaptureUi(),
            };
            var json = JSON.stringify(state);
            var blob = new Blob([json], { type: "application/json" });
            var url = URL.createObjectURL(blob);
            var a = document.createElement("a");
            var stamp = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 19);
            a.href = url;
            a.download = "taxtriage_state_" + stamp + ".json";
            document.body.appendChild(a);
            a.click();
            setTimeout(function () {
              document.body.removeChild(a);
              URL.revokeObjectURL(url);
            }, 0);
            _ttSetStateStatus("✓ State exported (" + boot.records.length + " rows).", false);
          } catch (e) {
            _ttSetStateStatus("Export failed: " + (e && e.message ? e.message : e), true);
          }
        };

        // ── Apply a parsed state envelope (or a bare BOOT-shaped object) ──────
        window.ttLoadState = function (state) {
          if (!state || typeof state !== "object") throw new Error("Not a valid state file");
          var boot = state.boot || (Array.isArray(state.records) ? state : null);
          if (!boot) throw new Error("State file has no data payload");
          _ttApplyBoot(boot);
          if (state.ui) _ttRestoreUi(state.ui);
        };

        // ── Read a state JSON file from the sidebar <input type=file> ─────────
        window.ttLoadStateFiles = function (files) {
          var f = files && files[0];
          if (!f) return;
          _ttSetStateStatus("Loading state…", false);
          var reader = new FileReader();
          reader.onload = function (e) {
            try {
              var state = JSON.parse(e.target.result);
              window.ttLoadState(state);
              _ttSetStateStatus("✓ State loaded from " + f.name + ".", false);
            } catch (err) {
              _ttSetStateStatus("Load failed: " + (err && err.message ? err.message : err), true);
            }
          };
          reader.onerror = function () {
            _ttSetStateStatus("Could not read " + f.name, true);
          };
          reader.readAsText(f);
          var inp = document.getElementById("state-load-input");
          if (inp) inp.value = "";
        };

        // ── GitHub-Pages startup dialog: Demo (default) / Empty / Upload ──────
        function _ttShowPagesDialog() {
          if (document.getElementById("tt-pages-overlay")) return;
          var ov = document.createElement("div");
          ov.id = "tt-pages-overlay";
          ov.style.cssText =
            "position:fixed;inset:0;z-index:99999;display:flex;align-items:center;" +
            "justify-content:center;background:rgba(13,40,66,0.55);backdrop-filter:blur(2px)";
          ov.innerHTML =
            '<div role="dialog" aria-modal="true" aria-labelledby="tt-pages-title" ' +
            'style="background:#fff;border-radius:12px;max-width:460px;width:90%;' +
            'box-shadow:0 18px 50px rgba(0,0,0,.35);overflow:hidden;font-family:inherit">' +
            '<div id="tt-pages-title" style="background:#1565c0;color:#fff;padding:14px 18px;' +
            'font-weight:700;font-size:1.05em"><i class="fas fa-vials" style="margin-right:.5em"></i>' +
            "TaxTriage Interactive Report</div>" +
            '<div style="padding:18px;color:#334;line-height:1.5;font-size:.92em">' +
            '<p style="margin:0 0 14px">Choose how to start. You can change the data source at any ' +
            "time from the <b>Upload Data</b> panel in the sidebar.</p>" +
            '<div style="display:flex;flex-direction:column;gap:10px">' +
            '<button id="tt-pages-demo" style="width:100%;padding:.7em;font-size:.95em;border:none;' +
            'border-radius:8px;background:#1565c0;color:#fff;cursor:pointer;font-weight:700">' +
            '<i class="fas fa-flask" style="margin-right:.45em"></i>Load demo data ' +
            '<span style="opacity:.85;font-weight:500">(recommended)</span></button>' +
            '<button id="tt-pages-empty" style="width:100%;padding:.7em;font-size:.92em;' +
            "border:1px solid #b0cce8;border-radius:8px;background:#f5f9ff;color:#1565c0;" +
            'cursor:pointer;font-weight:600"><i class="fas fa-file" style="margin-right:.45em"></i>' +
            "Start empty</button>" +
            '<button id="tt-pages-upload" style="width:100%;padding:.7em;font-size:.92em;' +
            "border:1px solid #b0cce8;border-radius:8px;background:#fff;color:#1565c0;" +
            'cursor:pointer;font-weight:600"><i class="fas fa-upload" style="margin-right:.45em"></i>' +
            "Upload your own data…</button>" +
            "</div></div></div>";
          document.body.appendChild(ov);
          function close() {
            ov.remove();
          }
          var demoBtn = document.getElementById("tt-pages-demo");
          var emptyBtn = document.getElementById("tt-pages-empty");
          var upBtn = document.getElementById("tt-pages-upload");
          if (demoBtn)
            demoBtn.addEventListener("click", function () {
              close();
            }); // demo already loaded on boot
          if (emptyBtn)
            emptyBtn.addEventListener("click", function () {
              _ttApplyBoot(_ttEmptyBoot());
              close();
            });
          if (upBtn)
            upBtn.addEventListener("click", function () {
              _ttApplyBoot(_ttEmptyBoot());
              close();
              var zone = document.getElementById("upload-zone");
              if (zone && zone.scrollIntoView) zone.scrollIntoView({ behavior: "smooth", block: "center" });
              var inp = document.getElementById("file-upload-input");
              if (inp)
                setTimeout(function () {
                  inp.click();
                }, 250);
            });
          ov.addEventListener("keydown", function (e) {
            if (e.key === "Escape") close();
          });
          setTimeout(function () {
            if (demoBtn) demoBtn.focus();
          }, 30);
        }

        // Only the GitHub-Pages build inlines its dataset as <script id="BOOTSTRAP">.
        // Pipeline reports set window.HEATMAP_BOOT and never carry that element.
        if (document.getElementById("BOOTSTRAP")) {
          if (document.readyState === "complete" || document.readyState === "interactive") {
            setTimeout(_ttShowPagesDialog, 0);
          } else {
            window.addEventListener("DOMContentLoaded", function () {
              setTimeout(_ttShowPagesDialog, 0);
            });
          }
        }
      })();

