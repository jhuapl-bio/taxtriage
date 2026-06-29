      /* ═══════════════════════════════════════════════════════════════════════════
       -  §  SAMPLE SIDEBAR  (right-side panel, visible on every tab)
       -     Owns: per-sample row UI, hide/show/remove buttons, drag-to-reorder,
       -     "Hide all / Show all" + "Remove all" bulk actions, sample-name
       -     hover tooltip, and the inline rename flow.
═══════════════════════════════════════════════════════════════════════════ */
      function removeSample(id) {
        const n = DATA.filter((r) => r["Specimen ID"] === id).length;
        if (
          !confirm(`Delete "${id}" and its ${n} row${n !== 1 ? "s" : ""} from the dataset?\n\nThis cannot be undone.`)
        )
          return;
        _deleteSampleData(id);
      }

      function _deleteSampleData(id) {
        for (let i = DATA.length - 1; i >= 0; i--) {
          if (DATA[i]["Specimen ID"] === id) DATA.splice(i, 1);
        }
        // Prune RUN_META so map markers reflect the deletion
        for (let i = RUN_META.length - 1; i >= 0; i--) {
          if (RUN_META[i].sample_name === id) RUN_META.splice(i, 1);
        }
        // Prune CONTIG_DATA so the histogram tab loses this sample's bars
        for (let i = CONTIG_DATA.length - 1; i >= 0; i--) {
          if (CONTIG_DATA[i].sample === id) CONTIG_DATA.splice(i, 1);
        }
        if (typeof _invalidateSummaryHistMap === "function") _invalidateSummaryHistMap();
        // Prune PROT so VF/AMR loses this sample's hits
        Object.keys(PROT || {}).forEach((k) => {
          if (!Array.isArray(PROT[k])) return;
          for (let i = PROT[k].length - 1; i >= 0; i--) {
            if (PROT[k][i] && PROT[k][i]["Specimen ID"] === id) PROT[k].splice(i, 1);
          }
        });
        // Prune BOOT in place — the source of resurrection on subsequent
        // file drops if we don't filter it here.
        if (BOOT) {
          if (Array.isArray(BOOT.records)) {
            for (let i = BOOT.records.length - 1; i >= 0; i--) {
              if (BOOT.records[i]["Specimen ID"] === id) BOOT.records.splice(i, 1);
            }
          }
          if (Array.isArray(BOOT.contig_data)) {
            for (let i = BOOT.contig_data.length - 1; i >= 0; i--) {
              if (BOOT.contig_data[i].sample === id) BOOT.contig_data.splice(i, 1);
            }
          }
          if (BOOT.prot_data && typeof BOOT.prot_data === "object") {
            Object.keys(BOOT.prot_data).forEach((k) => {
              if (!Array.isArray(BOOT.prot_data[k])) return;
              for (let i = BOOT.prot_data[k].length - 1; i >= 0; i--) {
                if (BOOT.prot_data[k][i] && BOOT.prot_data[k][i]["Specimen ID"] === id) {
                  BOOT.prot_data[k].splice(i, 1);
                }
              }
            });
          }
        }
        // Drop this sample from accumulated upload state too
        if (window._pruneUploadBuffersForSample) window._pruneUploadBuffersForSample(id);
        delete sampleColors[id];
        delete sampleHidden[id];
        delete sampleRescale[id];
        _sampleOrder = _sampleOrder.filter((s) => s !== id);
        // Re-evaluate VF/AMR tab visibility after the prune
        HAS_PROT = Object.keys(PROT || {}).some((k) => Array.isArray(PROT[k]) && PROT[k].length > 0);
        const protBtn = document.getElementById("prot-tab-btn");
        if (protBtn) protBtn.classList.toggle("hidden", !HAS_PROT);
        // Re-evaluate histogram tab visibility too
        const histBtn = document.getElementById("hist-tab-btn");
        if (histBtn) histBtn.classList.toggle("hidden", CONTIG_DATA.length === 0);
        // Prune SAMPLE_META so % classified drops the deleted sample's reads
        delete SAMPLE_META[id];
        if (BOOT && BOOT.sample_meta) delete BOOT.sample_meta[id];
        buildSampleList();
        buildTable();
        _rebuildMapMarkers();
        if (window._resetHistSelectors) window._resetHistSelectors();
        document.getElementById("banner-sub").textContent = _buildBannerSub();
        redraw();
      }

      function removeAllSamples() {
        const samples = uniq(DATA.map((r) => r["Specimen ID"] || "")).filter(Boolean);
        if (!samples.length) return;
        const ok = confirm(
          `Permanently delete all ${samples.length} sample${
            samples.length !== 1 ? "s" : ""
          } and their data?\n\nThis cannot be undone.`,
        );
        if (!ok) return;
        // Batch delete — clear EVERY data source, including BOOT and the
        // upload buffers. Without this, dropping a new JSON afterwards would
        // re-merge BOOT.records (and any previously-uploaded rows) on top of
        // the freshly-emptied dataset, resurrecting the samples the user
        // just told us to delete.
        DATA.length = 0;
        RUN_META.length = 0;
        CONTIG_DATA.length = 0;
        if (typeof _invalidateSummaryHistMap === "function") _invalidateSummaryHistMap();
        // Reset BOOT in place so all references see the change.
        if (BOOT) {
          if (Array.isArray(BOOT.records)) BOOT.records.length = 0;
          if (Array.isArray(BOOT.contig_data)) BOOT.contig_data.length = 0;
          if (BOOT.prot_data && typeof BOOT.prot_data === "object") {
            Object.keys(BOOT.prot_data).forEach((k) => {
              if (Array.isArray(BOOT.prot_data[k])) BOOT.prot_data[k].length = 0;
            });
          }
          BOOT.has_prot = false;
        }
        // Empty PROT in place, keep the same key set so the rest of the code
        // that does `PROT.per_gene_hits` etc. still finds an array.
        Object.keys(PROT || {}).forEach((k) => {
          if (Array.isArray(PROT[k])) PROT[k].length = 0;
        });
        HAS_PROT = false;
        const protBtn = document.getElementById("prot-tab-btn");
        if (protBtn) protBtn.classList.add("hidden");
        // Drop accumulated upload state so the next drop is a clean start.
        if (window._clearUploadBuffers) window._clearUploadBuffers();
        samples.forEach((id) => {
          delete sampleColors[id];
          delete sampleHidden[id];
          delete sampleRescale[id];
        });
        _sampleOrder = [];
        // Clear SAMPLE_META so % classified resets to N/A with no data
        Object.keys(SAMPLE_META).forEach((k) => delete SAMPLE_META[k]);
        if (BOOT && BOOT.sample_meta) Object.keys(BOOT.sample_meta).forEach((k) => delete BOOT.sample_meta[k]);
        buildSampleList();
        buildTable();
        _rebuildMapMarkers();
        // Hide the histogram tab again (it was shown when contigs existed)
        const histBtn = document.getElementById("hist-tab-btn");
        if (histBtn) histBtn.classList.add("hidden");
        if (window._resetHistSelectors) window._resetHistSelectors();
        document.getElementById("banner-sub").textContent = _buildBannerSub();
        redraw();
      }

      /* ── Hide / show every sample at once.
         If any sample is currently visible, hide them all. Otherwise (all
         already hidden, or no samples loaded) show them all. We mutate
         sampleHidden in bulk and call redraw() ONCE — the per-sample eye
         icons are refreshed by buildSampleList(). */
      function toggleAllSamples() {
        const samples = Array.from(_allSampleIds()).filter(Boolean);
        if (!samples.length) return;
        const anyVisible = samples.some((id) => !sampleHidden[id]);
        const next = anyVisible; // true = hide all; false = show all
        samples.forEach((id) => {
          sampleHidden[id] = next;
        });
        buildSampleList();
        _syncToggleAllSamplesBtn();
        redraw();
      }

      /* Keep the Hide All / Show All button label + icon in sync with the
         current visibility state. Called from buildSampleList() so it
         reflects single-sample eye toggles too. */
      function _syncToggleAllSamplesBtn() {
        const btn = document.getElementById("toggle-all-samples-btn");
        if (!btn) return;
        const samples = Array.from(_allSampleIds()).filter(Boolean);
        if (!samples.length) {
          btn.disabled = true;
          btn.style.opacity = "0.5";
          btn.style.cursor = "not-allowed";
          btn.innerHTML = '<i class="fas fa-eye-slash"></i> Hide All';
          return;
        }
        btn.disabled = false;
        btn.style.opacity = "1";
        btn.style.cursor = "pointer";
        const anyVisible = samples.some((id) => !sampleHidden[id]);
        if (anyVisible) {
          btn.innerHTML = '<i class="fas fa-eye-slash"></i> Hide All';
          btn.title = "Hide every sample from all charts and tables";
        } else {
          btn.innerHTML = '<i class="fas fa-eye"></i> Show All';
          btn.title = "Show every sample again";
        }
      }

      // Authoritative list of every sample processed in the run — NOT just those
      // with detections. Samples with no Kraken2/alignment hits never appear in
      // DATA (records), but they are still present in SAMPLE_META (per-sample
      // metadata is emitted for every sample) and may carry novelty / VF-AMR
      // evidence. Union all sources so the right-panel lists every sample.
      function _allSampleIds() {
        const ids = new Set();
        (DATA || []).forEach((r) => {
          const s = r["Specimen ID"] || "";
          if (s) ids.add(s);
        });
        // Every sample that produced per-sample metadata (includes zero-hit samples).
        Object.keys(SAMPLE_META || {}).forEach((s) => {
          if (s) ids.add(s);
        });
        // Samples that only surfaced via novelty detection.
        try {
          const nov = (NOVELTY && NOVELTY.samples) || {};
          Object.keys(nov).forEach((s) => {
            if (s) ids.add(s);
          });
        } catch (e) {}
        // Samples that only carry VF/AMR annotation rows.
        try {
          const scanProt = (rows) =>
            (rows || []).forEach((r) => {
              const s = r["Specimen ID"] || r.Sample || r.sample;
              if (s) ids.add(s);
            });
          scanProt(PROT.per_gene_hits);
          scanProt(PROT.amr_genes);
        } catch (e) {}
        return ids;
      }

      // True when a sample has no usable detections — few reads, no alignment, or
      // no rows in DATA. These samples still appear in the list but their plots /
      // scores can be unreliable, so the row gets a warning icon. Returns a
      // {warn, msg} object; warn=false means the sample looks fine.
      const _LOW_READ_THRESHOLD = 1000; // total reads below this == "few reads"
      function _sampleDataWarning(id) {
        const meta = SAMPLE_META[id] || {};
        const totalReads = parseInt(meta.total_reads) || 0;
        const alignedReads = parseInt(meta.aligned_reads) || 0;
        const hasRows = (DATA || []).some((r) => (r["Specimen ID"] || "") === id);
        const noAlign = alignedReads === 0;
        const fewReads = totalReads > 0 && totalReads < _LOW_READ_THRESHOLD;
        if (!noAlign && hasRows && !fewReads) return { warn: false, msg: "" };
        const bits = [];
        if (!hasRows || noAlign) bits.push("no alignment hits");
        if (fewReads) bits.push(`only ${totalReads.toLocaleString()} reads`);
        if (!bits.length && !hasRows) bits.push("no detections");
        const detail = bits.length ? ` (${bits.join(" · ")})` : "";
        return {
          warn: true,
          msg:
            `<b>${id}</b>${detail}<br>` +
            `<span style="color:#ccc">Few reads or no alignments — heatmap, TASS and coverage plots for this sample may be empty or unreliable.</span>`,
        };
      }

      function buildSampleList() {
        const rawSamples = Array.from(_allSampleIds()).filter(Boolean).sort();

        // Initialise per-sample state for any new samples
        rawSamples.forEach((id, i) => {
          if (!sampleColors[id]) sampleColors[id] = PALETTE[i % PALETTE.length];
          if (sampleHidden[id] === undefined) sampleHidden[id] = false;
          if (sampleRescale[id] === undefined) sampleRescale[id] = false;
        });

        // Build / sync _sampleOrder: keep existing order for known samples, append new ones
        // New samples (or first load) are natural-sorted (numeric-aware, case-insensitive).
        const knownSet = new Set(_sampleOrder);
        const newOnes = rawSamples
          .filter((id) => !knownSet.has(id))
          .sort((a, b) => a.localeCompare(b, undefined, { numeric: true, sensitivity: "base" }));
        _sampleOrder = _sampleOrder.filter((id) => rawSamples.includes(id));
        _sampleOrder.push(...newOnes);

        const cont = document.getElementById("sample-list");
        if (!cont) return;
        cont.innerHTML = "";

        // ── drag state ──────────────────────────────────────────────────────
        let _dragSrc = null;

        function _rebuildFromDOM() {
          const rows = cont.querySelectorAll(".sample-entry[data-sid]");
          _sampleOrder = Array.from(rows).map((r) => r.dataset.sid);
          redraw();
        }

        _sampleOrder.forEach((id) => {
          const div = document.createElement("div");
          div.className = "sample-entry";
          div.dataset.sid = id;
          div.draggable = true;

          // ── drag-and-drop ──────────────────────────────────────────────
          div.addEventListener("dragstart", (e) => {
            _dragSrc = div;
            e.dataTransfer.effectAllowed = "move";
            e.dataTransfer.setData("text/plain", id);
            div.style.opacity = "0.5";
          });
          div.addEventListener("dragend", () => {
            div.style.opacity = "";
          });
          div.addEventListener("dragover", (e) => {
            e.preventDefault();
            e.dataTransfer.dropEffect = "move";
          });
          div.addEventListener("dragenter", (e) => {
            e.preventDefault();
            div.style.background = "#e8f0fe";
          });
          div.addEventListener("dragleave", () => {
            div.style.background = "";
          });
          div.addEventListener("drop", (e) => {
            e.preventDefault();
            div.style.background = "";
            if (_dragSrc && _dragSrc !== div) {
              const allRows = Array.from(cont.querySelectorAll(".sample-entry[data-sid]"));
              const srcIdx = allRows.indexOf(_dragSrc);
              const dstIdx = allRows.indexOf(div);
              if (srcIdx < dstIdx) cont.insertBefore(_dragSrc, div.nextSibling);
              else cont.insertBefore(_dragSrc, div);
              _rebuildFromDOM();
            }
          });

          // ── drag handle ────────────────────────────────────────────────
          const grip = document.createElement("div");
          grip.innerHTML = "&#8942;&#8942;"; // ⋮⋮ two narrow dot columns — visual hint only
          grip.title = "Drag row to reorder";
          grip.className = "drag-grip";
          grip.style.pointerEvents = "none"; // row handles drag, grip is decoration

          // ── Color swatch ───────────────────────────────────────────────
          const colorIn = document.createElement("input");
          colorIn.type = "color";
          colorIn.value = sampleColors[id] || "#888888";
          colorIn.title = id;
          colorIn.addEventListener("input", (e) => {
            sampleColors[id] = e.target.value;
            _refreshMapMarkerColors();
            redraw();
          });

          // ── Label ──────────────────────────────────────────────────────
          const span = document.createElement("span");
          span.textContent = id;
          span.style.cursor = "grab";
          span.style.flex = "1";
          span.style.overflow = "hidden";
          span.style.textOverflow = "ellipsis";
          span.style.whiteSpace = "nowrap";
          span.style.userSelect = "none";
          // Custom tooltip — shows the full sample name on hover. We use the
          // app's own showTip()/moveTip()/hideTip() so it appears instantly
          // (no native ~500ms delay) and matches the tooltip used elsewhere.
          span.addEventListener("mouseover", (ev) => showTip(`<b>${id}</b>`, ev));
          span.addEventListener("mousemove", moveTip);
          span.addEventListener("mouseout", hideTip);

          // ── Low-reads / no-alignment warning icon ──────────────────────
          // Surfaced for samples with few reads or no alignment hits so the user
          // knows their plots/scores may be unreliable. Uses the app's instant
          // tooltip (showTip/moveTip/hideTip) to explain the caveat on hover.
          const _warn = _sampleDataWarning(id);
          let warnIcon = null;
          if (_warn.warn) {
            warnIcon = document.createElement("i");
            warnIcon.className = "fas fa-triangle-exclamation sample-warn-icon";
            warnIcon.style.color = "#e0a800";
            warnIcon.style.fontSize = "0.82em";
            warnIcon.style.cursor = "help";
            warnIcon.style.marginRight = "2px";
            warnIcon.addEventListener("mouseover", (ev) => showTip(_warn.msg, ev));
            warnIcon.addEventListener("mousemove", moveTip);
            warnIcon.addEventListener("mouseout", hideTip);
          }

          // ── Eye toggle button ──────────────────────────────────────────
          const btn = document.createElement("button");
          btn.className = "eye-btn";
          btn.title = "Toggle visibility";
          btn.type = "button";
          const icon = document.createElement("i");
          icon.className = `fas ${sampleHidden[id] ? "fa-eye-slash" : "fa-eye"}`;
          btn.appendChild(icon);
          btn.addEventListener("click", () => {
            sampleHidden[id] = !sampleHidden[id];
            icon.className = `fas ${sampleHidden[id] ? "fa-eye-slash" : "fa-eye"}`;
            btn.style.opacity = sampleHidden[id] ? "0.35" : "1";
            _syncToggleAllSamplesBtn();
            redraw();
          });
          btn.style.opacity = sampleHidden[id] ? "0.35" : "1";

          // ── ×100 rescale toggle ────────────────────────────────────────
          const rescaleBtn = document.createElement("button");
          rescaleBtn.className = "eye-btn";
          rescaleBtn.title = "Scale TASS Score & Coverage ×100 (legacy 0–1 files)";
          rescaleBtn.type = "button";
          rescaleBtn.textContent = "×100";
          rescaleBtn.style.fontSize = "0.68em";
          rescaleBtn.style.padding = "0 3px";
          rescaleBtn.style.opacity = sampleRescale[id] ? "1" : "0.3";
          rescaleBtn.style.color = sampleRescale[id] ? "#c06000" : "inherit";
          rescaleBtn.addEventListener("click", () => {
            sampleRescale[id] = !sampleRescale[id];
            rescaleBtn.style.opacity = sampleRescale[id] ? "1" : "0.3";
            rescaleBtn.style.color = sampleRescale[id] ? "#c06000" : "inherit";
            redraw();
          });

          // ── Rename button ──────────────────────────────────────────────
          const editBtn = document.createElement("button");
          editBtn.className = "eye-btn";
          editBtn.title = "Rename sample";
          editBtn.type = "button";
          editBtn.innerHTML = '<i class="fas fa-pen"></i>';
          editBtn.addEventListener("click", () => editSampleNameFromSidebar(id));

          // ── Attach-file button ─────────────────────────────────────────
          // Lets the user drop a new .json / .xlsx / .tsv / .txt / .csv next to
          // this specific sample. While the file is ingesting, the icon swaps to
          // a spinning bacteria so it's obvious the row is loading.
          const addBtn = document.createElement("button");
          addBtn.className = "eye-btn";
          addBtn.title = "Add a file (JSON / XLSX / TSV …) for this sample";
          addBtn.type = "button";
          const addIcon = document.createElement("i");
          const _isLoading = _loadingSampleIds.has(id);
          addIcon.className = _isLoading ? "fas fa-bacteria fa-spin" : "fas fa-paperclip";
          addIcon.style.color = _isLoading ? "#1565c0" : "inherit";
          addBtn.appendChild(addIcon);
          // Hidden file input dedicated to this row.
          const addInput = document.createElement("input");
          addInput.type = "file";
          addInput.accept = ".json,.xlsx,.xls,.tsv,.txt,.csv";
          addInput.multiple = true;
          addInput.style.display = "none";
          addInput.addEventListener("change", (e) => {
            if (e.target.files && e.target.files.length) window.addFileToSample(id, e.target.files);
            e.target.value = "";
          });
          addBtn.addEventListener("click", (e) => {
            e.stopPropagation();
            if (_loadingSampleIds.has(id)) return; // already busy
            addInput.click();
          });
          // Allow dragging a file straight onto the row's attach button.
          addBtn.addEventListener("dragover", (e) => {
            e.preventDefault();
            e.stopPropagation();
            e.dataTransfer.dropEffect = "copy";
          });
          addBtn.addEventListener("drop", (e) => {
            e.preventDefault();
            e.stopPropagation();
            if (e.dataTransfer.files && e.dataTransfer.files.length) window.addFileToSample(id, e.dataTransfer.files);
          });
          if (_isLoading) addBtn.style.opacity = "1";

          // ── Delete button ──────────────────────────────────────────────
          const delBtn = document.createElement("button");
          delBtn.className = "eye-btn";
          delBtn.title = "Delete sample";
          delBtn.type = "button";
          delBtn.innerHTML = '<i class="fas fa-trash"></i>';
          delBtn.style.color = "#c0392b";
          delBtn.style.opacity = "0.5";
          delBtn.addEventListener("mouseenter", () => {
            delBtn.style.opacity = "1";
          });
          delBtn.addEventListener("mouseleave", () => {
            delBtn.style.opacity = "0.5";
          });
          delBtn.addEventListener("click", () => removeSample(id));

          div.appendChild(colorIn);
          div.appendChild(span);
          if (warnIcon) div.appendChild(warnIcon);
          div.appendChild(editBtn);
          div.appendChild(rescaleBtn);
          div.appendChild(btn);
          div.appendChild(addBtn);
          div.appendChild(addInput);
          div.appendChild(delBtn);
          div.appendChild(grip);
          cont.appendChild(div);
        });
        // Keep the Hide All / Show All button label aligned with the
        // current visibility state on every list rebuild.
        _syncToggleAllSamplesBtn();
        // Rebuild per-sample-type TASS cutoff UI
        buildPerTypeTassUI();
        // Update the sidebar legend (sample swatches + TASS cutoffs)
        _updateSidebarLegend();
      }

      /** Rebuild the sidebar legend: sample color swatches + active TASS cutoffs.
       *  Called whenever buildSampleList() runs or TASS thresholds change. */
      function _updateSidebarLegend() {
        const panel = document.getElementById("sidebar-legend");
        const samplesEl = document.getElementById("sidebar-legend-samples");
        const tassEl = document.getElementById("sidebar-legend-tass");
        if (!panel || !samplesEl || !tassEl) return;

        // ── sample swatches ──────────────────────────────────────────────
        if (_sampleOrder.length > 0) {
          let swatchHtml =
            '<div style="font-size:0.72em;font-weight:600;color:#546e7a;margin-bottom:0.25em">Sample Colors</div>';
          swatchHtml += '<div style="display:flex;flex-direction:column;gap:0.18em">';
          _sampleOrder.forEach((id) => {
            const col = sampleColors[id] || "#888";
            const hidden = sampleHidden[id];
            swatchHtml +=
              `<div style="display:flex;align-items:center;gap:0.4em;opacity:${hidden ? 0.4 : 1}">` +
              `<span style="display:inline-block;width:11px;height:11px;border-radius:2px;background:${col};border:1px solid rgba(0,0,0,0.18);flex-shrink:0"></span>` +
              `<span style="font-size:0.78em;color:#37474f;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;max-width:130px" title="${id}">${id}</span>` +
              `</div>`;
          });
          swatchHtml += "</div>";
          samplesEl.innerHTML = swatchHtml;
          panel.style.display = "";
        } else {
          samplesEl.innerHTML = "";
          panel.style.display = "none";
        }

        // ── TASS cutoffs ─────────────────────────────────────────────────
        const cs = _tassCutoffSummary();
        let tassHtml =
          '<div style="font-size:0.72em;font-weight:600;color:#546e7a;margin:0.45em 0 0.25em">TASS Cutoffs</div>';
        if (cs.mode === "byType" && cs.items.length) {
          tassHtml += '<div style="display:flex;flex-direction:column;gap:0.15em">';
          cs.items.forEach(({ type, applied }) => {
            tassHtml +=
              `<div style="display:flex;align-items:center;gap:0.4em">` +
              `<span style="display:inline-block;width:11px;height:2px;background:#f59f00;flex-shrink:0;border-radius:1px"></span>` +
              `<span style="font-size:0.78em;color:#37474f">${type}: <b>${applied}</b></span>` +
              `</div>`;
          });
          tassHtml += "</div>";
        } else {
          const v = cs.global || 0;
          tassHtml +=
            `<div style="display:flex;align-items:center;gap:0.4em">` +
            `<span style="display:inline-block;width:11px;height:2px;background:#f59f00;flex-shrink:0;border-radius:1px"></span>` +
            `<span style="font-size:0.78em;color:#37474f">Global: <b>${v}</b></span>` +
            `</div>`;
        }
        tassEl.innerHTML = tassHtml;
      }

      /* ── Per-sample-type TASS helpers ────────────────────────────────────── */

      /** Return the recommended TASS threshold for a given sample type by
       *  inspecting SAMPLE_META entries that share the same type, then taking
       *  the median of their per-sample best_cutoffs thresholds. */
      function _defaultTassForType(sampleType) {
        const t = (sampleType || "").toLowerCase().trim();
        const vals = Object.values(SAMPLE_META)
          .filter((m) => (m.sample_type || "").toLowerCase().trim() === t)
          .map((m) => {
            const gran = m.preferred_granularity || "subkey";
            return ((m.best_cutoffs || {})[gran] || {}).best_threshold;
          })
          .filter((v) => v != null && !isNaN(v));
        if (!vals.length) {
          return BEST_TASS_THRESH != null && !isNaN(BEST_TASS_THRESH)
            ? Math.round(BEST_TASS_THRESH)
            : parseInt(document.getElementById("filter-min").value, 10) || 70;
        }
        vals.sort((a, b) => a - b);
        const mid = Math.floor(vals.length / 2);
        return Math.round(vals.length % 2 ? vals[mid] : (vals[mid - 1] + vals[mid]) / 2);
      }

      /* ── Per-sample-type TASS cutoff UI ──────────────────────────────────── */
      function buildPerTypeTassUI() {
        const wrap = document.getElementById("per-type-tass-wrap");
        const rowsEl = document.getElementById("per-type-tass-rows");
        const resetAllBtn = document.getElementById("per-type-reset-all-btn");
        const setAllBtn = document.getElementById("per-type-set-all-btn");
        if (!wrap || !rowsEl) return;

        const types = Array.from(
          new Set(DATA.map((r) => (r["Sample Type"] || "").trim().toLowerCase()).filter((t) => t && t !== "unknown")),
        ).sort();

        if (!types.length) {
          wrap.style.display = "none";
          // Show global slider as normal fallback
          const lbl = document.getElementById("filter-min-label");
          const wrapEl = document.getElementById("filter-min-wrap");
          if (lbl) lbl.style.display = "";
          if (wrapEl) wrapEl.style.display = "";
          return;
        }

        // Hide global slider — per-type controls replace it
        const globalLbl = document.getElementById("filter-min-label");
        const globalWrap = document.getElementById("filter-min-wrap");
        if (globalLbl) globalLbl.style.display = "none";
        if (globalWrap) globalWrap.style.display = "none";

        wrap.style.display = "";
        rowsEl.innerHTML = "";

        // Initialise defaults on first build; preserve user overrides on rebuild
        types.forEach((t) => {
          if (perTypeTass[t] == null) perTypeTass[t] = _defaultTassForType(t);
        });

        const singleType = types.length === 1;
        // which row is currently expanded (only for 2+ types)
        let activeType = singleType ? types[0] : null;

        // Hint for multi-type mode
        if (!singleType) {
          const hint = document.createElement("div");
          hint.style.cssText =
            "font-size:0.73em;color:#888;margin-bottom:0.4em;display:flex;align-items:center;gap:4px";
          hint.innerHTML =
            '<i class="fas fa-hand-pointer" style="font-size:0.9em"></i> Click a type to adjust its cutoff';
          rowsEl.appendChild(hint);
        }

        // One row per type
        types.forEach((t) => {
          const defaultVal = _defaultTassForType(t);

          /* ─ wrapper card ─ */
          const card = document.createElement("div");
          card.className = "pt-tass-card";
          card.dataset.type = t;
          card.style.cssText =
            "border:1px solid transparent;border-radius:5px;margin-bottom:3px;overflow:hidden;transition:border-color 0.15s,background 0.1s";

          /* ─ clickable header row ─ */
          const hdr = document.createElement("div");
          hdr.className = "pt-tass-hdr";
          hdr.style.cssText =
            "display:flex;align-items:center;gap:0.35em;padding:4px 7px;cursor:" +
            (singleType ? "default" : "pointer") +
            ";border-radius:5px;font-size:0.82em;transition:background 0.12s";

          const iconEl = document.createElement("i");
          iconEl.className = "fas fa-sliders";
          iconEl.style.cssText =
            "font-size:0.85em;color:#1565c0;opacity:" +
            (singleType ? "0.7" : "0.35") +
            ";transition:opacity 0.15s,transform 0.2s;flex-shrink:0;width:12px;text-align:center";

          const lblEl = document.createElement("span");
          lblEl.textContent = t;
          lblEl.style.cssText =
            "flex:1;text-transform:capitalize;font-weight:500;color:#333;overflow:hidden;text-overflow:ellipsis;white-space:nowrap";

          const valDisp = document.createElement("span");
          valDisp.className = "pt-val-disp";
          valDisp.textContent = (perTypeTass[t] ?? defaultVal).toFixed(0);
          valDisp.style.cssText =
            "font-size:0.9em;color:#1565c0;font-weight:600;min-width:2em;text-align:right;font-variant-numeric:tabular-nums";

          hdr.appendChild(iconEl);
          hdr.appendChild(lblEl);
          hdr.appendChild(valDisp);
          card.appendChild(hdr);

          /* ─ expandable slider panel ─ */
          const panel = document.createElement("div");
          panel.className = "pt-tass-panel";
          panel.style.cssText =
            "padding:5px 8px 7px 8px;background:#f0f4fc;border-top:1px solid #d0daf0;display:" +
            (singleType ? "block" : "none");

          const sliderRow = document.createElement("div");
          sliderRow.style.cssText = "display:flex;align-items:center;gap:0.4em";

          const rangeEl = document.createElement("input");
          rangeEl.type = "range";
          rangeEl.min = "0";
          rangeEl.max = "100";
          rangeEl.step = "1";
          rangeEl.value = perTypeTass[t] ?? defaultVal;
          rangeEl.style.cssText = "flex:1;cursor:pointer;accent-color:#1565c0";

          const numEl = document.createElement("input");
          numEl.type = "number";
          numEl.min = "0";
          numEl.max = "100";
          numEl.step = "1";
          numEl.value = perTypeTass[t] ?? defaultVal;
          numEl.style.cssText =
            "width:3.6em;padding:0.1em 0.2em;border:1px solid #b0bec5;border-radius:4px;font-size:0.88em;text-align:right";

          const resetBtn = document.createElement("button");
          resetBtn.innerHTML = '<i class="fas fa-rotate-left"></i>';
          resetBtn.title = `Reset to default (${defaultVal})`;
          resetBtn.style.cssText =
            "background:transparent;border:none;color:#90a4ae;cursor:pointer;font-size:0.95em;padding:0 2px;line-height:1;flex-shrink:0;transition:color 0.12s";
          resetBtn.addEventListener("mouseenter", () => (resetBtn.style.color = "#1565c0"));
          resetBtn.addEventListener("mouseleave", () => (resetBtn.style.color = "#90a4ae"));

          // Default label below slider
          const defLabel = document.createElement("div");
          defLabel.style.cssText = "font-size:0.7em;color:#90a4ae;margin-top:2px;text-align:right";
          defLabel.textContent = `default: ${defaultVal}`;

          const applyVal = (v) => {
            const clamped = Math.max(0, Math.min(100, Math.round(parseFloat(v))));
            perTypeTass[t] = clamped;
            rangeEl.value = clamped;
            numEl.value = clamped;
            valDisp.textContent = clamped.toFixed(0);
            if (setAllBtn) {
              setAllBtn.textContent = `Set All to ${clamped}`;
              setAllBtn.style.display = "";
            }
            _invalidateFilterCache();
            redraw();
            if (window._resetHistSelectors) window._resetHistSelectors();
          };

          rangeEl.addEventListener("input", () => applyVal(rangeEl.value));
          numEl.addEventListener("input", () => {
            if (!isNaN(parseFloat(numEl.value))) applyVal(numEl.value);
          });
          resetBtn.addEventListener("click", (e) => {
            e.stopPropagation();
            applyVal(defaultVal);
          });

          sliderRow.appendChild(rangeEl);
          sliderRow.appendChild(numEl);
          sliderRow.appendChild(resetBtn);
          panel.appendChild(sliderRow);
          panel.appendChild(defLabel);
          card.appendChild(panel);
          rowsEl.appendChild(card);

          /* ─ hover + click behaviour (multi-type mode only) ─ */
          if (!singleType) {
            hdr.addEventListener("mouseenter", () => {
              hdr.style.background = "#e8f0fe";
              iconEl.style.opacity = "1";
              iconEl.style.transform = "rotate(-25deg) scale(1.1)";
            });
            hdr.addEventListener("mouseleave", () => {
              if (activeType !== t) {
                hdr.style.background = "";
                iconEl.style.opacity = "0.35";
                iconEl.style.transform = "";
              }
            });
            hdr.addEventListener("click", () => {
              const isOpen = panel.style.display !== "none";
              // Close all panels first
              rowsEl.querySelectorAll(".pt-tass-card").forEach((c) => {
                c.querySelector(".pt-tass-panel").style.display = "none";
                c.style.borderColor = "transparent";
                c.style.background = "";
                const ic = c.querySelector(".pt-tass-hdr i");
                if (ic) {
                  ic.style.opacity = "0.35";
                  ic.style.transform = "";
                }
                c.querySelector(".pt-tass-hdr").style.background = "";
              });
              activeType = null;
              if (setAllBtn) setAllBtn.style.display = "none";

              if (!isOpen) {
                // Open this panel
                panel.style.display = "block";
                card.style.borderColor = "#90caf9";
                card.style.background = "#f7fbff";
                hdr.style.background = "#e8f0fe";
                iconEl.style.opacity = "1";
                iconEl.style.transform = "rotate(-25deg) scale(1.1)";
                activeType = t;
                const curVal = perTypeTass[t] ?? defaultVal;
                if (setAllBtn) {
                  setAllBtn.textContent = `Set All to ${Math.round(curVal)}`;
                  setAllBtn.style.display = "";
                }
              }
            });
          }
        });

        /* ─ Reset All button ─ */
        if (resetAllBtn) {
          // Re-attach to avoid duplicate listeners from rebuilds
          const newReset = resetAllBtn.cloneNode(true);
          resetAllBtn.parentNode.replaceChild(newReset, resetAllBtn);
          newReset.addEventListener("click", () => {
            types.forEach((t) => {
              const def = _defaultTassForType(t);
              perTypeTass[t] = def;
              const card = rowsEl.querySelector(`[data-type="${t}"]`);
              if (!card) return;
              const vd = card.querySelector(".pt-val-disp");
              if (vd) vd.textContent = def.toFixed(0);
              const panel = card.querySelector(".pt-tass-panel");
              if (panel) {
                const r = panel.querySelector("input[type=range]");
                const n = panel.querySelector("input[type=number]");
                if (r) r.value = def;
                if (n) n.value = def;
              }
            });
            _invalidateFilterCache();
            redraw();
            if (window._resetHistSelectors) window._resetHistSelectors();
          });
        }

        /* ─ Set All button ─ */
        if (setAllBtn) {
          const newSetAll = setAllBtn.cloneNode(true);
          setAllBtn.parentNode.replaceChild(newSetAll, setAllBtn);
          newSetAll.style.display = "none";
          newSetAll.addEventListener("click", () => {
            if (!activeType || perTypeTass[activeType] == null) return;
            const val = perTypeTass[activeType];
            types.forEach((t) => {
              perTypeTass[t] = val;
              const card = rowsEl.querySelector(`[data-type="${t}"]`);
              if (!card) return;
              const vd = card.querySelector(".pt-val-disp");
              if (vd) vd.textContent = val.toFixed(0);
              const panel = card.querySelector(".pt-tass-panel");
              if (panel) {
                const r = panel.querySelector("input[type=range]");
                const n = panel.querySelector("input[type=number]");
                if (r) r.value = val;
                if (n) n.value = val;
              }
            });
            _invalidateFilterCache();
            redraw();
            if (window._resetHistSelectors) window._resetHistSelectors();
            _updateSidebarLegend();
          });
        }
      }

