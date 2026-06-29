      /* ═══════════════════════════════════════════════════════════════════════════
       -  §  UTILITIES
       -     Pure helpers (isTruthy / num / uniq), the memoized filteredData()
       -     pipeline + cache invalidator, applyRescale(), and the shared
       -     tooltip (showTip / moveTip / hideTip).
       -     filteredData() is the hot path — every tab calls it to get the
       -     current filtered, sample-hidden, rescaled view of DATA.
═══════════════════════════════════════════════════════════════════════════ */
      function isTruthy(v) {
        if (v === true || v === 1) return true;
        if (typeof v === "string") {
          const s = v.trim().toLowerCase();
          return s === "true" || s === "yes" || s === "1" || s === "t";
        }
        return false;
      }

      /* ── filteredData() — memoized.
         A single redraw triggers ~10 calls; with 14k rows this was the
         dominant cost. We compute a cheap fingerprint of every input that
         affects the result and reuse the previously computed array if the
         fingerprint hasn't changed. The cache is invalidated automatically
         (no manual bookkeeping needed) because every UI change alters at
         least one of the inputs we hash. */
      let _FD_CACHE = { key: null, value: null };
      function _invalidateFilterCache() {
        _FD_CACHE.key = null;
        _FD_CACHE.value = null;
      }
      function _hasAnyRescale() {
        for (const k in sampleRescale) if (sampleRescale[k]) return true;
        return false;
      }
      function _hashHidden() {
        // sampleHidden / sampleRescale typically have <50 entries — cheap.
        const h = [];
        for (const k in sampleHidden) if (sampleHidden[k]) h.push(k);
        return h.sort().join("|");
      }
      function _hashRescale() {
        const r = [];
        for (const k in sampleRescale) if (sampleRescale[k]) r.push(k);
        return r.sort().join("|");
      }
      function _hashPerTypeTass() {
        return Object.entries(perTypeTass)
          .filter(([, v]) => v != null)
          .sort(([a], [b]) => a.localeCompare(b))
          .map(([k, v]) => `${k}:${v}`)
          .join("|");
      }
      // Classify a data row into a kingdom/domain filter bucket.
      function _rowKingdomBucket(r) {
        const domain = (r["Domain"] || "").trim();
        const kingdom = (r["Kingdom"] || "").trim();
        if (domain === "Bacteria") return "Bacteria";
        if (domain === "Archaea") return "Archaea";
        if (domain === "Viruses" || domain === "Virus") return "Viruses";
        if (domain === "Eukaryota") {
          const kl = kingdom.toLowerCase();
          if (kl.includes("fung") || kl.includes("mycota")) return "Fungi";
          if (kl.includes("plant") || kl.includes("viridiplant") || kl.includes("streptophyt")) return "Plantae";
          if (kl.includes("animal") || kl.includes("metazoa")) return "Animalia";
          return "OtherEuk";
        }
        // Fallback: infer from kingdom string when Domain is absent
        const kl = kingdom.toLowerCase();
        if (kl.includes("bacteri")) return "Bacteria";
        if (kl.includes("virus") || kl.includes("virae")) return "Viruses";
        if (kl.includes("archae")) return "Archaea";
        if (kl.includes("fung")) return "Fungi";
        if (kl.includes("plant") || kl.includes("viridiplant")) return "Plantae";
        if (kl.includes("animal") || kl.includes("metazoa")) return "Animalia";
        return null;
      }
      // Search scope + last-known regex validity (surfaced in the sidebar status line).
      let _searchInvalid = false;
      function filteredData() {
        const txt = document.getElementById("filter-text").value.trim();
        const ic = document.getElementById("filter-ic").checked;
        const scope = (document.getElementById("filter-scope") || {}).value || "both";
        const viewLevel = (document.getElementById("view-level") || {}).value || "Strain";
        const minV = parseFloat(document.getElementById("filter-min").value) || 0;
        const onlyHC = document.getElementById("filter-hc").checked;
        const onlyP = document.getElementById("filter-pass").checked;
        const mcSel = Array.from(document.getElementById("filter-mc").selectedOptions).map((o) => o.value);
        const mtDNA = document.getElementById("filter-mt-dna")
          ? document.getElementById("filter-mt-dna").checked
          : true;
        const mtRNA = document.getElementById("filter-mt-rna")
          ? document.getElementById("filter-mt-rna").checked
          : true;
        const mtBoth = document.getElementById("filter-mt-both")
          ? document.getElementById("filter-mt-both").checked
          : true;
        const kingdomSel = Array.from(document.querySelectorAll(".fk-cb:checked")).map((cb) => cb.value);

        // Fingerprint everything that can change the result. Note: DATA.length
        // covers most upload-driven changes; if you ever mutate DATA in place
        // without changing length, call _invalidateFilterCache() afterwards.
        const cacheKey = [
          txt,
          ic ? 1 : 0,
          scope,
          minV,
          onlyHC ? 1 : 0,
          onlyP ? 1 : 0,
          mcSel.join(","),
          mtDNA ? 1 : 0,
          mtRNA ? 1 : 0,
          mtBoth ? 1 : 0,
          kingdomSel.join(","),
          _hashHidden(),
          _hashRescale(),
          DATA.length,
          _hashPerTypeTass(),
          ROLLUP_PASS ? 1 : 0,
          viewLevel,
          watchFilterMode,
          watchFilterMode === "all" ? "" : Array.from(watchlist).sort().join(","),
        ].join("\u0001");
        if (_FD_CACHE.key === cacheKey && _FD_CACHE.value) return _FD_CACHE.value;

        // Compile the search pattern defensively. While the user is still
        // typing, an incomplete pattern (e.g. "(" or "[a-") would otherwise
        // throw and break every view. On error we fall back to "no text
        // filter" (results stay visible) and flag the box as invalid.
        let rx = null;
        _searchInvalid = false;
        if (txt) {
          try {
            rx = new RegExp(txt, ic ? "i" : "");
          } catch (e) {
            rx = null;
            _searchInvalid = true;
          }
        }
        const mcFilter = mcSel.length > 0 ? new Set(mcSel) : null; // null = no filter
        const kingdomFilter = kingdomSel.length > 0 ? new Set(kingdomSel) : null;
        const anyRescale = _hasAnyRescale();

        // Base (non-level, non-threshold) predicates shared by every view.
        function _basePass(r) {
          if (sampleHidden[r["Specimen ID"]]) return false;
          if (rx) {
            const inSample = rx.test(r["Specimen ID"] || "");
            const inOrg = rx.test(r["Detected Organism"] || "");
            if (scope === "sample") {
              if (!inSample) return false;
            } else if (scope === "organism") {
              if (!inOrg) return false;
            } else if (!inSample && !inOrg) {
              return false;
            }
          }
          if (onlyHC && !isTruthy(r["High Consequence"])) return false;
          if (mcFilter && !mcFilter.has(r["Microbial Category"] || "Unknown")) return false;
          if (kingdomFilter) {
            const bucket = _rowKingdomBucket(r);
            if (!bucket || !kingdomFilter.has(bucket)) return false;
          }
          const mt = (r["Mol Type"] || "").toLowerCase();
          if (mt === "dna" && !mtDNA) return false;
          if (mt === "rna" && !mtRNA) return false;
          if (mt !== "dna" && mt !== "rna" && !mtBoth) return false;
          // Follow-up (watchlist) view filter — applies across every tab.
          if (watchFilterMode !== "all" && watchlist.size) {
            const watched = watchlist.has(_watchKey(r));
            if (watchFilterMode === "only" && !watched) return false;
            if (watchFilterMode === "hide" && watched) return false;
          }
          return true;
        }
        const _lvlOf = (r) => r["Level"] || "Strain";

        const filtered = DATA.filter((r) => {
          if (_lvlOf(r) !== viewLevel) return false;
          if (!_basePass(r)) return false;
          const pi = rowPassInfo(r);
          // Threshold gate: when rollup is on, keep a row whose species or
          // genus aggregation passes even if the row itself is below cutoff.
          if (ROLLUP_PASS) {
            if (!isNaN(pi.strain) && !pi.effectivePass) return false;
          } else {
            if (!isNaN(pi.strain) && pi.strain < pi.thr) return false;
          }
          // "Passes threshold" filter stays STRICT (own-level) so users can
          // still isolate rows that pass on their own merits.
          if (onlyP && !pi.strainPass) return false;
          return true;
        });

        // ── Species rollup in the Strain view ──────────────────────────────
        // When every visible strain of a species is below the cutoff but the
        // species aggregation passes, collapse those strains into a single
        // Species summary row (flagged so the UI can badge it). Species where
        // at least one strain passes keep their detailed strain rows.
        let finalRows = filtered;
        if (viewLevel === "Strain" && ROLLUP_PASS) {
          const speciesRowByKey = new Map();
          for (const r of DATA) {
            if ((r["Level"] || "") === "Species") {
              speciesRowByKey.set((r["Specimen ID"] || "") + "" + (r["Subkey"] || ""), r);
            }
          }
          const groups = new Map();
          for (const r of filtered) {
            const k = (r["Specimen ID"] || "") + "" + (r["Subkey"] || "");
            if (!groups.has(k)) groups.set(k, []);
            groups.get(k).push(r);
          }
          finalRows = [];
          const handled = new Set();
          for (const r of filtered) {
            const k = (r["Specimen ID"] || "") + "" + (r["Subkey"] || "");
            if (handled.has(k)) continue;
            const grpRows = groups.get(k);
            const anyStrainPass = grpRows.some((g) => rowPassInfo(g).strainPass);
            const specRow = speciesRowByKey.get(k);
            if (!anyStrainPass && specRow && rowPassInfo(specRow).effectivePass) {
              // promote the species summary row in place of the failing strains
              const inj = Object.assign({}, specRow);
              inj.__rescuedSpecies = true;
              inj.__rescuedStrainCount = grpRows.length;
              finalRows.push(inj);
              handled.add(k);
            } else {
              finalRows.push(r);
            }
          }
        }
        // Skip the .map(applyRescale) allocation entirely when no sample is
        // currently rescaled — saves a 14k-object copy on every filter pass.
        const out = anyRescale ? finalRows.map(applyRescale) : finalRows;
        _FD_CACHE = { key: cacheKey, value: out };
        return out;
      }

      /* Return a shallow copy of row r with TASS Score + Coverage scaled ×100
        when sampleRescale is enabled for legacy 0–1 data. */
      const _RESCALE_COLS = ["TASS Score", "Coverage", "Breadth %"];
      function applyRescale(r) {
        if (!sampleRescale[r["Specimen ID"]]) return r;
        const copy = Object.assign({}, r);
        for (const col of _RESCALE_COLS) {
          const v = parseFloat(copy[col]);
          if (!isNaN(v)) copy[col] = Math.round(v * 1000) / 10;
        }
        return copy;
      }

      /* ── Threshold rollup (effective pass) ──────────────────────────────
         A strain may fall below the TASS cutoff while its parent species or
         genus aggregation passes. We compute an "effective pass" that rolls
         up: a row counts as effectively passing if the strain OR its species
         OR its genus clears the threshold. The strict strain-level pass is
         still kept for the "Passes Threshold" filter and badges.            */
      function _rescaleVal(r, v) {
        return _hasAnyRescale() && sampleRescale[r["Specimen ID"]] ? v * 100 : v;
      }
      function thresholdForRow(r) {
        const minV = parseFloat(document.getElementById("filter-min").value) || 0;
        const _st = (r["Sample Type"] || "").trim().toLowerCase();
        return _st && _st !== "unknown" && perTypeTass[_st] != null ? perTypeTass[_st] : minV;
      }
      function rowPassInfo(r) {
        const thr = thresholdForRow(r);
        const strain = _rescaleVal(r, parseFloat(r["TASS Score"]));
        const species = _rescaleVal(r, parseFloat(r["Species TASS"]));
        const genus = _rescaleVal(r, parseFloat(r["Genus TASS"]));
        const sP = !isNaN(strain) && strain >= thr;
        const spP = !isNaN(species) && species >= thr;
        const gP = !isNaN(genus) && genus >= thr;
        return {
          thr,
          strain,
          species,
          genus,
          strainPass: sP,
          speciesPass: spP,
          genusPass: gP,
          effectivePass: sP || spP || gP,
          rescued: !sP && (spP || gP),
          rescueLevel: sP ? null : spP ? "species" : gP ? "genus" : null,
        };
      }
      // Global toggle: when true, the threshold filter uses effective (rolled-up)
      // pass so species/genus-rescued strains stay visible. Default on.
      let ROLLUP_PASS = true;

      function uniq(arr) {
        return [...new Set(arr)];
      }
      function num(v) {
        const n = parseFloat(v);
        return isNaN(n) ? 0 : n;
      }

      let tooltipsEnabled = true;

      const tooltip = document.getElementById("tooltip");
      function showTip(html, event) {
        if (!tooltipsEnabled) return;
        tooltip.innerHTML = html;
        // Park at 0,0 before measuring so offsetWidth reflects the true
        // unconstrained width (not a clipped value from the previous position).
        tooltip.style.left = "0px";
        tooltip.style.top = "0px";
        tooltip.style.opacity = 1;
        moveTip(event);
      }
      function moveTip(event) {
        const x = event.clientX + 12,
          y = event.clientY - 28;
        // offsetWidth after parking at 0,0 (or from previous frame) is reliable.
        const tw = tooltip.offsetWidth || 520;
        const th = tooltip.offsetHeight || 0;
        const left = Math.min(x, window.innerWidth - tw - 16);
        const top = Math.min(Math.max(y, 4), window.innerHeight - th - 8);
        tooltip.style.left = left + "px";
        tooltip.style.top = top + "px";
      }
      function hideTip() {
        tooltip.style.opacity = 0;
      }

      /* ── Strain-level hover tooltip helpers ─────────────────────────────── */
      /**
       * Return up to `n` Strain-level rows from DATA that belong to the same
       * species (or genus) and sample as `row`, sorted by TASS desc.
       */
      function _getTopStrainsForRow(row, n) {
        n = n || 5;
        var lvl = row["Level"] || "Strain";
        if (lvl === "Strain") return [];
        var sampleId = row["Specimen ID"];
        var strains = DATA.filter(function (r) {
          if ((r["Level"] || "Strain") !== "Strain") return false;
          if (r["Specimen ID"] !== sampleId) return false;
          if (lvl === "Species") return r["Species Name"] === row["Species Name"];
          if (lvl === "Genus") return r["Genus Name"] === row["Genus Name"];
          return false;
        });
        strains.sort(function (a, b) {
          return (parseFloat(b["TASS Score"]) || 0) - (parseFloat(a["TASS Score"]) || 0);
        });
        return strains.slice(0, n);
      }

      /**
       * Build the inner HTML for a strain-breakdown tooltip table.
       */
      function _strainTooltipHtml(strains, rowLevel, parentName) {
        if (!strains || !strains.length) return null;
        var title = rowLevel === "Genus" ? "Top strains in genus (by TASS)" : "Top strains in species (by TASS)";
        // Build regex to strip shared parent-name prefix from strain labels
        var _prefixRe = null;
        if (parentName) {
          var _esc = parentName.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
          _prefixRe = new RegExp("^" + _esc + "\\s*", "i");
        }
        var html =
          '<div style="font-weight:700;margin-bottom:5px;font-size:0.9em;border-bottom:1px solid rgba(255,255,255,0.25);padding-bottom:4px">' +
          title +
          "</div>";
        html +=
          '<table style="border-collapse:collapse;font-size:0.82em;width:100%">' +
          "<thead><tr>" +
          '<th style="text-align:left;padding:2px 8px 3px 0;color:#ffd580;font-weight:600;white-space:nowrap">Strain</th>' +
          '<th style="text-align:right;padding:2px 6px 3px;color:#ffd580;font-weight:600;white-space:nowrap">TASS</th>' +
          '<th style="text-align:right;padding:2px 6px 3px;color:#ffd580;font-weight:600;white-space:nowrap">% Reads</th>' +
          '<th style="text-align:right;padding:2px 0 3px 6px;color:#ffd580;font-weight:600;white-space:nowrap"># Reads</th>' +
          "</tr></thead><tbody>";
        strains.forEach(function (s) {
          var tass = parseFloat(s["TASS Score"]);
          var pct = parseFloat(s["% Reads"]);
          var reads = s["# Reads Aligned"];
          var tassStr = isNaN(tass) ? "—" : tass.toFixed(1);
          var pctStr = isNaN(pct) ? "—" : pct.toFixed(2) + "%";
          var readsStr = _fmtBig(reads).short;
          var rawName = s["Detected Organism"] || "—";
          var displayName = rawName;
          if (_prefixRe && _prefixRe.test(rawName)) {
            var suffix = rawName.replace(_prefixRe, "").trim();
            displayName = suffix
              ? '<span style="color:rgba(255,255,255,0.4)" title="' + rawName + '">…</span> ' + suffix
              : rawName;
          }
          html +=
            "<tr>" +
            '<td style="padding:2px 8px 2px 0;max-width:220px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap"><i>' +
            displayName +
            "</i></td>" +
            '<td style="text-align:right;padding:2px 6px;font-weight:700">' +
            tassStr +
            "</td>" +
            '<td style="text-align:right;padding:2px 6px">' +
            pctStr +
            "</td>" +
            '<td style="text-align:right;padding:2px 0 2px 6px">' +
            readsStr +
            "</td>" +
            "</tr>";
        });
        html += "</tbody></table>";
        return html;
      }

      /* ── Plot + table exports ────────────────────────────────────────────── */
      const _EXPORT_STATE = { target: null, mode: null, observer: null, timer: null };
      const _TABLE_SETTINGS = {};

      function _slug(s) {
        return (
          String(s || "taxtriage-export")
            .trim()
            .toLowerCase()
            .replace(/[^a-z0-9]+/g, "-")
            .replace(/^-+|-+$/g, "")
            .slice(0, 80) || "taxtriage-export"
        );
      }

      function _downloadBlob(blob, filename) {
        const a = document.createElement("a");
        a.href = URL.createObjectURL(blob);
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        setTimeout(() => {
          URL.revokeObjectURL(a.href);
          a.remove();
        }, 250);
      }

      function _downloadText(text, filename, type) {
        _downloadBlob(new Blob([text], { type: type || "text/plain;charset=utf-8" }), filename);
      }

      function _nearestTitle(el) {
        const wrap = el.closest(".chart-wrap, .sun-panel, .hist-bar-track, .pane") || el.parentElement;
        const title =
          wrap?.querySelector(".chart-title, .sun-charthead, .hist-bar-title, .compare-title")?.textContent ||
          document.querySelector(`.tab-btn[data-tab="${activeTab}"]`)?.textContent ||
          document.title ||
          "TaxTriage export";
        return title.replace(/\s+/g, " ").trim();
      }

      function _ensureExportModal() {
        let overlay = document.getElementById("export-overlay");
        if (overlay) return overlay;
        overlay = document.createElement("div");
        overlay.id = "export-overlay";
        overlay.className = "export-modal-overlay";
        overlay.innerHTML = `
          <div class="export-modal" role="dialog" aria-modal="true" aria-labelledby="export-title">
            <header>
              <i class="fas fa-file-export"></i>
              <span id="export-title">Export</span>
              <button type="button" id="export-close" title="Close">x</button>
            </header>
            <div class="export-modal-body">
              <label>Format
                <select id="export-format"></select>
              </label>
              <label id="export-delim-wrap">Delimiter
                <select id="export-delim">
                  <option value=",">Comma (,)</option>
                  <option value="\\t">Tab</option>
                  <option value=";">Semicolon (;)</option>
                  <option value="|">Pipe (|)</option>
                  <option value="custom">Custom</option>
                </select>
              </label>
              <label id="export-custom-delim-wrap" style="display:none">Custom delimiter
                <input id="export-custom-delim" value="," maxlength="8" />
              </label>
              <label id="export-width-wrap">Width (px)
                <input id="export-width" type="number" min="50" step="50" />
              </label>
              <label id="export-height-wrap">Height (px)
                <input id="export-height" type="number" min="50" step="50" />
              </label>
            </div>
            <div class="export-modal-actions">
              <button type="button" id="export-cancel">Cancel</button>
              <button type="button" class="primary" id="export-save">Save</button>
            </div>
          </div>`;
        document.body.appendChild(overlay);
        const close = () => (overlay.style.display = "none");
        overlay.querySelector("#export-close").addEventListener("click", close);
        overlay.querySelector("#export-cancel").addEventListener("click", close);
        overlay.addEventListener("click", (e) => {
          if (e.target === overlay) close();
        });
        overlay.querySelector("#export-delim").addEventListener("change", (e) => {
          overlay.querySelector("#export-custom-delim-wrap").style.display = e.target.value === "custom" ? "" : "none";
        });
        overlay.querySelector("#export-format").addEventListener("change", () => _syncExportModalFields());
        overlay.querySelector("#export-save").addEventListener("click", () => _runExportFromModal());
        return overlay;
      }

      function _syncExportModalFields() {
        const overlay = document.getElementById("export-overlay");
        if (!overlay) return;
        const fmt = overlay.querySelector("#export-format").value;
        const isTable = _EXPORT_STATE.mode === "table";
        const isDelimited = fmt === "csv" || fmt === "tsv";
        overlay.querySelector("#export-delim-wrap").style.display = isTable && isDelimited ? "" : "none";
        overlay.querySelector("#export-custom-delim-wrap").style.display =
          isTable && isDelimited && overlay.querySelector("#export-delim").value === "custom" ? "" : "none";
        overlay.querySelector("#export-width-wrap").style.display = isTable ? "none" : "";
        overlay.querySelector("#export-height-wrap").style.display = isTable ? "none" : "";
      }

      function _openPlotExport(svg) {
        const overlay = _ensureExportModal();
        _EXPORT_STATE.mode = "plot";
        _EXPORT_STATE.target = svg;
        overlay.querySelector("#export-title").textContent = "Export Plot";
        overlay.querySelector("#export-format").innerHTML =
          '<option value="png">PNG</option><option value="jpeg">JPEG</option><option value="svg">SVG</option><option value="pdf">PDF</option><option value="html">HTML</option>';
        const size = _plotNaturalSize(svg);
        const w = Math.max(50, Math.round(size.width || 900));
        const h = Math.max(50, Math.round(size.height || 500));
        overlay.querySelector("#export-width").value = w;
        overlay.querySelector("#export-height").value = h;
        _syncExportModalFields();
        overlay.style.display = "flex";
      }

      function _openTableExport(table) {
        const overlay = _ensureExportModal();
        _EXPORT_STATE.mode = "table";
        _EXPORT_STATE.target = table;
        overlay.querySelector("#export-title").textContent = "Export Table";
        overlay.querySelector("#export-format").innerHTML =
          '<option value="xlsx">XLSX</option><option value="csv">CSV</option><option value="tsv">TSV</option>';
        overlay.querySelector("#export-delim").value = ",";
        _syncExportModalFields();
        overlay.style.display = "flex";
      }

      function _cloneSvgWithStyles(svg, width, height) {
        const clone = svg.cloneNode(true);
        clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
        clone.setAttribute("width", width);
        clone.setAttribute("height", height);
        if (!clone.getAttribute("viewBox")) {
          const ow = parseFloat(svg.getAttribute("width")) || width;
          const oh = parseFloat(svg.getAttribute("height")) || height;
          clone.setAttribute("viewBox", `0 0 ${ow} ${oh}`);
        }
        const src = [svg, ...svg.querySelectorAll("*")];
        const dst = [clone, ...clone.querySelectorAll("*")];
        const props = [
          "font-family",
          "font-size",
          "font-weight",
          "font-style",
          "fill",
          "stroke",
          "stroke-width",
          "opacity",
          "text-anchor",
          "dominant-baseline",
          "shape-rendering",
        ];
        for (let i = 0; i < src.length; i++) {
          const cs = getComputedStyle(src[i]);
          props.forEach((p) => {
            const v = cs.getPropertyValue(p);
            if (v) dst[i].style.setProperty(p, v);
          });
        }
        return new XMLSerializer().serializeToString(clone);
      }

      function _svgTextToImage(svgText) {
        return new Promise((resolve, reject) => {
          const img = new Image();
          const url = URL.createObjectURL(new Blob([svgText], { type: "image/svg+xml;charset=utf-8" }));
          img.onload = () => {
            URL.revokeObjectURL(url);
            resolve(img);
          };
          img.onerror = () => {
            URL.revokeObjectURL(url);
            reject(new Error("Could not render SVG for export."));
          };
          img.src = url;
        });
      }

      function _plotSvgs(target) {
        const svgs =
          target && target.tagName && target.tagName.toLowerCase() === "svg"
            ? [target]
            : Array.from(target.querySelectorAll("svg"));
        return svgs.filter((svg) => {
          const w = parseFloat(svg.getAttribute("width")) || svg.getBoundingClientRect().width;
          const h = parseFloat(svg.getAttribute("height")) || svg.getBoundingClientRect().height;
          return w >= 90 && h >= 60 && getComputedStyle(svg).display !== "none";
        });
      }

      function _plotNaturalSize(target) {
        if (target && target.tagName && target.tagName.toLowerCase() === "svg") {
          return {
            width:
              parseFloat(target.getAttribute("width")) ||
              target.viewBox?.baseVal?.width ||
              target.getBoundingClientRect().width ||
              900,
            height:
              parseFloat(target.getAttribute("height")) ||
              target.viewBox?.baseVal?.height ||
              target.getBoundingClientRect().height ||
              500,
          };
        }
        const rect = target.getBoundingClientRect();
        return {
          width: Math.max(target.scrollWidth || 0, rect.width || 0, 900),
          height: Math.max(target.scrollHeight || 0, rect.height || 0, 500),
        };
      }

      function _plotToSvgText(target, width, height) {
        if (target && target.tagName && target.tagName.toLowerCase() === "svg") {
          return _cloneSvgWithStyles(target, width, height);
        }
        const natural = _plotNaturalSize(target);
        const base = target.getBoundingClientRect();
        const pieces = _plotSvgs(target).map((svg) => {
          const r = svg.getBoundingClientRect();
          const x = Math.max(0, r.left - base.left + (target.scrollLeft || 0));
          const y = Math.max(0, r.top - base.top + (target.scrollTop || 0));
          const w = parseFloat(svg.getAttribute("width")) || r.width;
          const h = parseFloat(svg.getAttribute("height")) || r.height;
          const inner = _cloneSvgWithStyles(svg, w, h);
          return `<g transform="translate(${x.toFixed(2)},${y.toFixed(2)})">${inner}</g>`;
        });
        return `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${
          natural.width
        } ${natural.height}"><rect width="${natural.width}" height="${natural.height}" fill="#fff"/>${pieces.join(
          "",
        )}</svg>`;
      }

      async function _exportPlot(target, opts) {
        const title = _nearestTitle(target);
        const filename = _slug(title);
        const natural = _plotNaturalSize(target);
        const width = Math.max(50, parseInt(opts.width, 10) || natural.width || 900);
        const height = Math.max(50, parseInt(opts.height, 10) || natural.height || 500);
        const svgText = _plotToSvgText(target, width, height);

        if (opts.format === "svg") {
          _downloadText(svgText, `${filename}.svg`, "image/svg+xml;charset=utf-8");
          return;
        }
        if (opts.format === "html") {
          const html =
            '<!doctype html><html><head><meta charset="utf-8"><title>' +
            title.replace(/</g, "&lt;") +
            "</title><style>body{margin:24px;font-family:system-ui,sans-serif;background:#fff;color:#222}svg{max-width:100%;height:auto}</style></head><body>" +
            '<h1 style="font-size:18px;margin:0 0 14px">' +
            title.replace(/</g, "&lt;") +
            "</h1>" +
            svgText +
            "</body></html>";
          _downloadText(html, `${filename}.html`, "text/html;charset=utf-8");
          return;
        }

        const img = await _svgTextToImage(svgText);
        const canvas = document.createElement("canvas");
        canvas.width = width;
        canvas.height = height;
        const ctx = canvas.getContext("2d");
        ctx.fillStyle = "#ffffff";
        ctx.fillRect(0, 0, width, height);
        ctx.drawImage(img, 0, 0, width, height);

        if (opts.format === "pdf") {
          const jsPDF = window.jspdf && window.jspdf.jsPDF;
          if (!jsPDF) {
            alert("PDF export requires the jsPDF library. Try PNG/JPEG/SVG, or open the report with internet access.");
            return;
          }
          const orient = width >= height ? "landscape" : "portrait";
          const pdf = new jsPDF({ orientation: orient, unit: "pt", format: [width, height] });
          pdf.addImage(canvas.toDataURL("image/jpeg", 0.95), "JPEG", 0, 0, width, height);
          pdf.save(`${filename}.pdf`);
          return;
        }

        const mime = opts.format === "jpeg" ? "image/jpeg" : "image/png";
        canvas.toBlob(
          (blob) => {
            if (blob) _downloadBlob(blob, `${filename}.${opts.format === "jpeg" ? "jpg" : "png"}`);
          },
          mime,
          0.95,
        );
      }

      function _tableToRows(table) {
        const rows = [];
        table.querySelectorAll("tr").forEach((tr) => {
          if (tr.offsetParent === null && getComputedStyle(tr).display === "none") return;
          const cells = Array.from(tr.children).filter((c) => c.matches("th,td"));
          if (!cells.length) return;
          rows.push(
            cells
              .filter((c) => getComputedStyle(c).display !== "none")
              .map((c) => c.innerText.replace(/\s+/g, " ").trim()),
          );
        });
        return rows;
      }

      function _delimitedEscape(value, delimiter) {
        const s = value == null ? "" : String(value);
        return /["\r\n]/.test(s) || s.includes(delimiter) ? `"${s.replace(/"/g, '""')}"` : s;
      }

      function _exportTable(table, opts) {
        const rows = _tableToRows(table);
        if (!rows.length) {
          alert("No visible table rows to export.");
          return;
        }
        const title = _nearestTitle(table);
        const filename = _slug(title || table.id || "taxtriage-table");
        if (opts.format === "xlsx") {
          if (typeof XLSX === "undefined") {
            alert("XLSX export requires the SheetJS library. Try CSV or TSV.");
            return;
          }
          const ws = XLSX.utils.aoa_to_sheet(rows);
          const wb = XLSX.utils.book_new();
          XLSX.utils.book_append_sheet(wb, ws, "Table");
          XLSX.writeFile(wb, `${filename}.xlsx`);
          return;
        }
        let delimiter = opts.delimiter === "custom" ? opts.customDelimiter || "," : opts.delimiter || ",";
        if (delimiter === "\\t") delimiter = "\t";
        const ext = opts.format === "tsv" ? "tsv" : "csv";
        const text = rows.map((r) => r.map((v) => _delimitedEscape(v, delimiter)).join(delimiter)).join("\r\n");
        _downloadText(text, `${filename}.${ext}`, "text/plain;charset=utf-8");
      }

      function _runExportFromModal() {
        const overlay = document.getElementById("export-overlay");
        if (!overlay || !_EXPORT_STATE.target) return;
        const opts = {
          format: overlay.querySelector("#export-format").value,
          width: overlay.querySelector("#export-width").value,
          height: overlay.querySelector("#export-height").value,
          delimiter: overlay.querySelector("#export-delim").value,
          customDelimiter: overlay.querySelector("#export-custom-delim").value,
        };
        overlay.style.display = "none";
        if (_EXPORT_STATE.mode === "plot") _exportPlot(_EXPORT_STATE.target, opts);
        else _exportTable(_EXPORT_STATE.target, opts);
      }

      function _plotExportHost(svg) {
        return (
          svg.closest(
            "#heatmap-svg-wrap, #tass-svg-wrap, #coverage-svg-wrap, #summary-genus-wrap, #prot-genus-svg, #prot-prop-svg, #prot-compare-svg, #longi-chart-wrap, #hist-main, #explore-bubble-wrap, #explore-radar-wrap, #explore-chord-wrap, #explore-lolli-wrap, #explore-corr-wrap, .sun-panel-svg-wrap, #xs-dist-svg, #xs-pca-svg, #xs-cooc-svg, #cmp-par-svg, .chart-wrap, .sun-panel",
          ) || svg.parentElement
        );
      }

      function _tableExportHost(table) {
        return table.parentElement;
      }

      function _keepButtonInScrollport(host, btn) {
        if (!host || !btn || btn._scrollPinned) return;
        btn._scrollPinned = true;
        const update = () => {
          btn.style.transform = `translate(${host.scrollLeft || 0}px, ${host.scrollTop || 0}px)`;
        };
        host.addEventListener("scroll", update, { passive: true });
        update();
      }

      function _tableSettingsKey(table) {
        if (table.id) return `table:${table.id}`;
        const host = _tableExportHost(table);
        const hostId = host && host.id ? host.id : _slug(_nearestTitle(table));
        const headers = _tableColumns(table)
          .map((c) => c.label)
          .join("|");
        return `table:${hostId}:${headers}`;
      }

      function _tableColumns(table) {
        const headerRow =
          table.querySelector("thead tr") ||
          Array.from(table.querySelectorAll("tr")).find((tr) => tr.querySelector("th,td"));
        if (!headerRow) return [];
        return Array.from(headerRow.children)
          .filter((c) => c.matches("th,td"))
          .map((cell, index) => ({
            index,
            label: cell.innerText.replace(/\s+/g, " ").trim() || `Column ${index + 1}`,
          }));
      }

      function _tableSettingsFor(table) {
        const key = _tableSettingsKey(table);
        if (!_TABLE_SETTINGS[key]) _TABLE_SETTINGS[key] = { hidden: {}, fontSize: "", fontFamily: "" };
        return _TABLE_SETTINGS[key];
      }

      function _autoPrintTableFontSize(visibleCount) {
        if (visibleCount <= 4) return "13px";
        if (visibleCount <= 6) return "12px";
        if (visibleCount <= 8) return "10px";
        if (visibleCount <= 12) return "9px";
        if (visibleCount <= 16) return "8px";
        return "7px";
      }

      function _printFontFromSetting(fontSize, visibleCount) {
        if (!fontSize) return _autoPrintTableFontSize(visibleCount);
        const mapped = {
          "0.68em": "8px",
          "0.8em": "10px",
          "0.95em": "12px",
        };
        return mapped[fontSize] || fontSize;
      }

      function _applyTableSettings(table) {
        if (!table) return;
        const settings = _tableSettingsFor(table);
        table.style.fontSize = settings.fontSize || "";
        table.style.fontFamily = settings.fontFamily || "";
        const cols = _tableColumns(table);
        const hiddenIdx = new Set(cols.filter((c) => settings.hidden[c.label]).map((c) => c.index));
        const visibleCount = Math.max(1, cols.length - hiddenIdx.size);
        table.dataset.visibleCols = String(visibleCount);
        table.style.setProperty("--tt-table-print-font-size", _printFontFromSetting(settings.fontSize, visibleCount));
        table.querySelectorAll("colgroup col").forEach((col, idx) => {
          col.style.display = hiddenIdx.has(idx) ? "none" : "";
        });
        table.querySelectorAll("tr").forEach((tr) => {
          const cells = Array.from(tr.children).filter((c) => c.matches("th,td"));
          if (cells.length !== cols.length) return;
          cells.forEach((cell, idx) => {
            cell.style.display = hiddenIdx.has(idx) ? "none" : "";
          });
        });
      }

      function _ensureTableSettingsModal() {
        let overlay = document.getElementById("table-settings-overlay");
        if (overlay) return overlay;
        overlay = document.createElement("div");
        overlay.id = "table-settings-overlay";
        overlay.className = "export-modal-overlay";
        overlay.innerHTML = `
          <div class="export-modal" role="dialog" aria-modal="true" aria-labelledby="table-settings-title">
            <header>
              <i class="fas fa-pen-to-square"></i>
              <span id="table-settings-title">Table Display</span>
              <button type="button" id="table-settings-close" title="Close">x</button>
            </header>
            <div class="table-settings-body">
              <div style="color:#667;font-size:.86em;line-height:1.4">
                Select the columns to show. Deselected columns are hidden in the report, exports, and PDF print view.
              </div>
              <label>Shown columns
                <select id="table-settings-columns" multiple></select>
              </label>
              <div class="table-settings-tools">
                <button type="button" id="table-settings-all">Show all</button>
                <button type="button" id="table-settings-none">Hide all</button>
              </div>
              <div class="table-settings-inline">
                <label>Font size
                  <select id="table-settings-font-size">
                    <option value="">Default</option>
                    <option value="0.68em">Compact</option>
                    <option value="0.8em">Normal</option>
                    <option value="0.95em">Large</option>
                    <option value="11px">11 px</option>
                    <option value="12px">12 px</option>
                    <option value="14px">14 px</option>
                  </select>
                </label>
                <label>Font type
                  <select id="table-settings-font-family">
                    <option value="">Default</option>
                    <option value="system-ui, -apple-system, Segoe UI, Roboto, sans-serif">System Sans</option>
                    <option value="Arial, Helvetica, sans-serif">Arial</option>
                    <option value="Georgia, serif">Georgia</option>
                    <option value="'Courier New', monospace">Monospace</option>
                  </select>
                </label>
              </div>
            </div>
            <div class="export-modal-actions">
              <button type="button" id="table-settings-cancel">Cancel</button>
              <button type="button" class="primary" id="table-settings-apply">Apply</button>
            </div>
          </div>`;
        document.body.appendChild(overlay);
        const close = () => (overlay.style.display = "none");
        overlay.querySelector("#table-settings-close").addEventListener("click", close);
        overlay.querySelector("#table-settings-cancel").addEventListener("click", close);
        overlay.addEventListener("click", (e) => {
          if (e.target === overlay) close();
        });
        overlay.querySelector("#table-settings-all").addEventListener("click", () => {
          overlay.querySelectorAll("#table-settings-columns option").forEach((opt) => (opt.selected = true));
        });
        overlay.querySelector("#table-settings-none").addEventListener("click", () => {
          overlay.querySelectorAll("#table-settings-columns option").forEach((opt) => (opt.selected = false));
        });
        overlay.querySelector("#table-settings-apply").addEventListener("click", () => {
          const table = overlay._tableTarget;
          if (!table) return close();
          const settings = _tableSettingsFor(table);
          const shown = new Set(
            Array.from(overlay.querySelector("#table-settings-columns").selectedOptions).map((opt) => opt.value),
          );
          settings.hidden = {};
          _tableColumns(table).forEach((col) => {
            if (!shown.has(col.label)) settings.hidden[col.label] = true;
          });
          settings.fontSize = overlay.querySelector("#table-settings-font-size").value;
          settings.fontFamily = overlay.querySelector("#table-settings-font-family").value;
          _applyTableSettings(table);
          close();
        });
        return overlay;
      }

      function _openTableSettings(table) {
        const overlay = _ensureTableSettingsModal();
        overlay._tableTarget = table;
        const cols = _tableColumns(table);
        const settings = _tableSettingsFor(table);
        const colSel = overlay.querySelector("#table-settings-columns");
        colSel.innerHTML = "";
        cols.forEach((col) => {
          const opt = document.createElement("option");
          opt.value = col.label;
          opt.textContent = col.label;
          opt.selected = !settings.hidden[col.label];
          colSel.appendChild(opt);
        });
        overlay.querySelector("#table-settings-font-size").value = settings.fontSize || "";
        overlay.querySelector("#table-settings-font-family").value = settings.fontFamily || "";
        overlay.querySelector("#table-settings-title").textContent = `Table Display: ${_nearestTitle(table)}`;
        overlay.style.display = "flex";
      }

      function _enhanceExports() {
        document.querySelectorAll("svg").forEach((svg) => {
          if (svg.closest("#tooltip, #embed-detail, .export-modal-overlay, .export-control, #help-overlay")) return;
          const w = parseFloat(svg.getAttribute("width")) || svg.getBoundingClientRect().width;
          const h = parseFloat(svg.getAttribute("height")) || svg.getBoundingClientRect().height;
          if (w < 90 || h < 60) return;
          const host = _plotExportHost(svg);
          if (!host || host.querySelector(":scope > .plot-export-btn")) return;
          host.classList.add("plot-export-host");
          const target = host.id === "hist-main" || host.querySelectorAll("svg").length > 1 ? host : svg;
          const btn = document.createElement("button");
          btn.type = "button";
          btn.className = "plot-export-btn export-control";
          btn.title = "Export plot";
          btn.innerHTML = '<i class="fas fa-download"></i>';
          btn.addEventListener("click", (e) => {
            e.stopPropagation();
            _openPlotExport(target);
          });
          host.appendChild(btn);
          _keepButtonInScrollport(host, btn);
        });

        document.querySelectorAll("table").forEach((table) => {
          if (table.closest("#tooltip, #embed-detail, .export-modal-overlay, .export-control")) return;
          const host = _tableExportHost(table);
          if (!host) return;
          host.classList.add("table-export-host");
          _applyTableSettings(table);
          if (!host.querySelector(":scope > .table-export-btn")) {
            const btn = document.createElement("button");
            btn.type = "button";
            btn.className = "table-export-btn export-control";
            btn.title = "Export table";
            btn.innerHTML = '<i class="fas fa-file-export"></i>';
            btn.addEventListener("click", (e) => {
              e.stopPropagation();
              _openTableExport(table);
            });
            host.appendChild(btn);
            _keepButtonInScrollport(host, btn);
          }
          if (!host.querySelector(":scope > .table-settings-btn")) {
            const settingsBtn = document.createElement("button");
            settingsBtn.type = "button";
            settingsBtn.className = "table-settings-btn export-control";
            settingsBtn.title = "Edit table columns and font";
            settingsBtn.innerHTML = '<i class="fas fa-pen-to-square"></i>';
            settingsBtn.addEventListener("click", (e) => {
              e.stopPropagation();
              _openTableSettings(table);
            });
            host.appendChild(settingsBtn);
            _keepButtonInScrollport(host, settingsBtn);
          }
        });
      }

      function _scheduleExportEnhance() {
        clearTimeout(_EXPORT_STATE.timer);
        _EXPORT_STATE.timer = setTimeout(_enhanceExports, 60);
      }

      function _initExportEnhancer() {
        _ensureExportModal();
        _scheduleExportEnhance();
        if (_EXPORT_STATE.observer) return;
        _EXPORT_STATE.observer = new MutationObserver((records) => {
          for (const rec of records) {
            const target = rec.target && rec.target.nodeType === 1 ? rec.target : null;
            if (target && target.closest("#tooltip, #embed-detail, .export-modal-overlay")) continue;
            for (const node of rec.addedNodes || []) {
              if (node.nodeType !== 1) continue;
              if (node.closest && node.closest("#tooltip, #embed-detail, .export-modal-overlay")) continue;
              if (
                node.matches?.("svg,table") ||
                node.querySelector?.("svg,table") ||
                target?.matches?.("svg,table") ||
                target?.querySelector?.("svg,table")
              ) {
                _scheduleExportEnhance();
                return;
              }
            }
          }
        });
        _EXPORT_STATE.observer.observe(document.body, { childList: true, subtree: true });
      }

      function _ensureReportPdfModal() {
        let overlay = document.getElementById("report-pdf-overlay");
        if (overlay) return overlay;
        overlay = document.createElement("div");
        overlay.id = "report-pdf-overlay";
        overlay.className = "export-modal-overlay";
        overlay.innerHTML = `
          <div class="export-modal" role="dialog" aria-modal="true" aria-labelledby="report-pdf-title">
            <header>
              <i class="fas fa-file-pdf"></i>
              <span id="report-pdf-title">Export Report PDF</span>
              <button type="button" id="report-pdf-close" title="Close">x</button>
            </header>
            <div style="padding:1em 1.05em;color:#334;line-height:1.45;font-size:.88em">
              <p style="margin:0 0 .65em">
                This will prepare every report tab using the current filters, sample visibility, chart controls, and table state.
              </p>
              <p style="margin:0 0 .8em;color:#667;font-size:.86em">
                Your browser print dialog will open next. Select <b>Save as PDF</b> as the destination to write the report.
              </p>
              <label
                style="display:flex;align-items:center;gap:.5em;font-size:.85em;color:#334;font-weight:600;margin-bottom:.4em"
              >
                Sample-color &amp; TASS-cutoff legend
                <select id="report-pdf-legend-place" style="margin-left:auto;font-size:.95em;padding:2px 6px">
                  <option value="cover" selected>On cover page</option>
                  <option value="footer">Repeat on every page</option>
                  <option value="off">Don't include</option>
                </select>
              </label>
              <p style="margin:0;color:#8a93a3;font-size:.78em">
                The legend maps each plot color to its sample and records the applied TASS cutoffs.
              </p>
            </div>
            <div class="export-modal-actions">
              <button type="button" id="report-pdf-cancel">Cancel</button>
              <button type="button" class="primary" id="report-pdf-confirm">
                <i class="fas fa-print"></i> Prepare PDF
              </button>
            </div>
          </div>`;
        document.body.appendChild(overlay);
        const close = () => (overlay.style.display = "none");
        overlay.querySelector("#report-pdf-close").addEventListener("click", close);
        overlay.querySelector("#report-pdf-cancel").addEventListener("click", close);
        overlay.addEventListener("click", (e) => {
          if (e.target === overlay) close();
        });
        overlay.querySelector("#report-pdf-confirm").addEventListener("click", () => {
          overlay.style.display = "none";
          _exportReportPdf().catch((err) => {
            _hidePdfProgress();
            alert("PDF export failed: " + (err && err.message ? err.message : err));
          });
        });
        return overlay;
      }

      function _openReportPdfModal() {
        const overlay = _ensureReportPdfModal();
        overlay.style.display = "flex";
      }

      function _visibleReportTabs() {
        return Array.from(document.querySelectorAll(".tab-btn[data-tab]"))
          .filter((btn) => !btn.classList.contains("hidden"))
          .map((btn) => btn.dataset.tab)
          .filter((tab) => document.getElementById(`pane-${tab}`));
      }

      function _restoreActiveTab(tab) {
        document.querySelectorAll(".tab-btn").forEach((b) => b.classList.toggle("active", b.dataset.tab === tab));
        document.querySelectorAll(".pane").forEach((p) => p.classList.remove("active"));
        const pane = document.getElementById(`pane-${tab}`);
        if (pane) pane.classList.add("active");
        activeTab = tab;
      }

      function _pdfDelay(ms) {
        return new Promise((resolve) => setTimeout(resolve, ms));
      }

      function _setPdfProgress(message, done, total, subtext) {
        const overlay = document.getElementById("pdf-progress-overlay");
        const text = document.getElementById("pdf-progress-text");
        const fill = document.getElementById("pdf-progress-fill");
        const sub = document.getElementById("pdf-progress-sub");
        if (overlay) overlay.style.display = "flex";
        if (text) text.textContent = message || "Preparing PDF…";
        const pct = total ? Math.max(0, Math.min(100, Math.round((done / total) * 100))) : 0;
        if (fill) fill.style.width = pct + "%";
        if (sub) sub.textContent = subtext || "The print dialog will open when this is ready.";
      }

      function _hidePdfProgress() {
        const overlay = document.getElementById("pdf-progress-overlay");
        if (overlay) overlay.style.display = "none";
      }

      function _tabLabel(tab) {
        const btn = document.querySelector(`.tab-btn[data-tab="${tab}"]`);
        return btn ? btn.textContent.replace(/\s+/g, " ").trim() : tab;
      }

      /** Renders all enabled metadata sub-tabs so their charts appear in the PDF. */
      async function _renderMetaSubTabsForPdf() {
        const subIds = ["longi", "geo", "host", "cmp"];
        const originalSub = _activeMetaSub;
        for (const id of subIds) {
          const btn = document.querySelector(`.meta-subtab[data-metasub="${id}"]`);
          if (!btn || btn.disabled) continue;
          _switchMetaSub(id);
          // For the compare sub-tab, render all four analysis types
          if (id === "cmp" && typeof _buildComparison === "function") {
            const origAnalysis = _cmpAnalysis;
            const analyses = ["cosine", "jaccard", "enrichment", "profile"];
            const chart = document.getElementById("cmp-chart");
            // Create a print-only container to hold all analysis snapshots
            let printSnap = document.getElementById("cmp-print-snapshots");
            if (!printSnap) {
              printSnap = document.createElement("div");
              printSnap.id = "cmp-print-snapshots";
              printSnap.style.cssText = "display:none";
              const cmpPane = document.getElementById("meta-subpane-cmp");
              if (cmpPane) cmpPane.appendChild(printSnap);
            }
            printSnap.innerHTML = "";
            for (const analysis of analyses) {
              _cmpAnalysis = analysis;
              // Reflect active state on buttons
              document
                .querySelectorAll("#cmp-analysis-tabs .cmp-atab")
                .forEach((b) => b.classList.toggle("active", b.getAttribute("data-analysis") === analysis));
              _buildComparison();
              await _pdfDelay(300);
              if (chart && chart.querySelector("svg")) {
                const section = document.createElement("div");
                section.className = "meta-subtabs-pdf-section";
                section.style.cssText = "margin-bottom:1.2em";
                const heading = document.createElement("div");
                heading.style.cssText =
                  "font-size:0.82em;font-weight:700;color:#1565c0;margin-bottom:0.4em;text-transform:uppercase;letter-spacing:.04em";
                const labels = {
                  cosine: "Cosine Similarity",
                  jaccard: "Jaccard Overlap",
                  enrichment: "Enrichment",
                  profile: "Profile Heatmap",
                };
                heading.textContent = labels[analysis] || analysis;
                section.appendChild(heading);
                Array.from(chart.querySelectorAll("svg")).forEach((svg) => section.appendChild(svg.cloneNode(true)));
                printSnap.appendChild(section);
              }
            }
            // Restore active analysis
            _cmpAnalysis = origAnalysis;
            _buildComparison();
          } else {
            await _pdfDelay(id === "longi" ? 400 : 250);
          }
        }
        // Restore original sub-tab
        if (originalSub) _switchMetaSub(originalSub);
      }

      async function _renderTabsForPdf() {
        const originalTab = activeTab;
        const tabs = _visibleReportTabs();
        for (let i = 0; i < tabs.length; i++) {
          const tab = tabs[i];
          _setPdfProgress(`Rendering ${_tabLabel(tab)}…`, i, tabs.length);
          _restoreActiveTab(tab);
          if (tab === "runmeta") {
            if (typeof _buildRunMetaTable === "function") _buildRunMetaTable();
            if (typeof _updateMetaSubTabStates === "function") _updateMetaSubTabStates();
            await _pdfDelay(80);
            _setPdfProgress("Rendering metadata sub-tabs…", i, tabs.length);
            await _renderMetaSubTabsForPdf();
            // Make all rendered sub-tab snapshots visible for print
            const snap = document.getElementById("cmp-print-snapshots");
            if (snap) snap.style.display = "";
            await _pdfDelay(80);
          } else if (tab === "map") {
            if (typeof L !== "undefined" && typeof _initMap === "function") {
              _initMap();
              if (typeof _rebuildMapMarkers === "function") _rebuildMapMarkers();
              if (typeof _refreshMapMarkerColors === "function") _refreshMapMarkerColors();
              if (_leafletMap) {
                _leafletMap.invalidateSize();
                setTimeout(() => _leafletMap && _leafletMap.invalidateSize(), 180);
              }
              await _pdfDelay(900);
            }
          } else {
            const pane = document.getElementById(`pane-${tab}`);
            const hasRenderedContent = !!(pane && pane.querySelector("svg,table"));
            if (_TAB_DIRTY[tab] || !hasRenderedContent) _drawTab(tab);
            await _pdfDelay(tab === "sunburst" ? 950 : 120);
          }

          // If this is the summary tab, also render all cross-sample sub-views
          // and snapshot them into a print-only container so they appear in the PDF.
          if (tab === "summary") {
            _setPdfProgress("Rendering Cross-Sample plots…", i, tabs.length + 1);
            await _renderXsForPdf();
          }

          _setPdfProgress(`Finished ${_tabLabel(tab)}.`, i + 1, tabs.length);
        }
        _restoreActiveTab(originalTab);
        _setPdfProgress("Finalizing print layout…", tabs.length, tabs.length);
        await _pdfDelay(150);
      }

      /** Renders each xs sub-view and clones the resulting SVG into a print-only div. */
      async function _renderXsForPdf() {
        // Ensure the xs aggregate is populated. _drawCrossSample is the
        // canonical initializer; fall back to any legacy aliases.
        if (typeof _XS !== "undefined" && !_XS.lastAgg) {
          const fd = typeof filteredData === "function" ? filteredData() : [];
          const samps = typeof uniq === "function" ? uniq(fd.map((r) => r["Specimen ID"] || "")).filter(Boolean) : [];
          if (typeof _drawCrossSample === "function") _drawCrossSample(fd, samps);
          else if (typeof _drawXS === "function") _drawXS();
          else if (typeof _xsBuild === "function") _xsBuild();
          await _pdfDelay(200);
        }

        const printWrap = document.getElementById("xs-print-wrap");
        if (!printWrap) return;
        printWrap.innerHTML = "";

        // Guard: nothing to render if aggregate is still missing
        if (typeof _XS === "undefined" || !_XS.lastAgg) return;

        const xsViews = [
          { key: "dist", label: "Prevalence & TASS" },
          { key: "pca", label: "Sample Clustering" },
          { key: "cooc", label: "Co-occurrence" },
        ];

        const origView = _XS.view;
        for (const { key, label } of xsViews) {
          _XS.view = key;
          document.querySelectorAll(".xs-subtab").forEach((b) => b.classList.toggle("active", b.dataset.xs === key));
          if (typeof _xsRenderBody === "function") _xsRenderBody();
          // Longer delay for co-occurrence (O(n²)) and PCA
          await _pdfDelay(key === "cooc" || key === "pca" ? 500 : 350);
          const body = document.getElementById("xs-body");
          // Clone all SVGs in the view (some views render more than one)
          const svgs = body ? Array.from(body.querySelectorAll("svg")) : [];
          if (svgs.length) {
            const section = document.createElement("div");
            section.style.cssText = "margin-top:1em;break-inside:avoid;page-break-inside:avoid";
            const heading = document.createElement("div");
            heading.style.cssText = "font-size:0.88em;font-weight:600;color:#1565c0;margin-bottom:0.3em";
            heading.textContent = label;
            section.appendChild(heading);
            svgs.forEach((svg) => section.appendChild(svg.cloneNode(true)));
            printWrap.appendChild(section);
          }
        }

        // Restore the original view
        _XS.view = origView || "table";
        if (typeof _xsRenderBody === "function") _xsRenderBody();
      }

      /* ── TASS cutoff summary (shared by the KPI card logic + the PDF legend) ──
         Returns { mode:"byType"|"global", items:[{type,applied,def}], global, recommended }. */
      function _tassCutoffSummary() {
        const types = Array.from(
          new Set(DATA.map((r) => (r["Sample Type"] || "").trim().toLowerCase()).filter((t) => t && t !== "unknown")),
        ).sort();
        const globalFallback = parseFloat((document.getElementById("filter-min") || {}).value) || 0;
        const recommended = BEST_TASS_THRESH != null && !isNaN(BEST_TASS_THRESH) ? Number(BEST_TASS_THRESH) : null;
        if (types.length) {
          return {
            mode: "byType",
            items: types.map((t) => ({
              type: t,
              applied: perTypeTass[t] != null ? perTypeTass[t] : globalFallback,
              def: _defaultTassForType(t),
            })),
            global: globalFallback,
            recommended,
          };
        }
        return { mode: "global", items: [], global: globalFallback, recommended };
      }

      function _legendEsc(s) {
        return String(s == null ? "" : s)
          .replace(/&/g, "&amp;")
          .replace(/</g, "&lt;")
          .replace(/>/g, "&gt;");
      }

      // Ordered, currently-visible sample list with resolved swatch colors.
      function _legendSamples() {
        const visible = uniq(DATA.map((r) => r["Specimen ID"] || ""))
          .filter(Boolean)
          .filter((id) => !sampleHidden[id]);
        return _orderedSamples(visible).map((id) => ({ id, color: sampleColors[id] || "#90a4ae" }));
      }

      /* Build the printable legend (sample colors + TASS cutoffs + marker key).
         placement: "cover" (full block on page 1) | "footer" (compact strip
         repeated on every page) | "off" (clear both). */
      function _buildPdfLegend(placement) {
        const coverEl = document.getElementById("pdf-legend");
        const footEl = document.getElementById("pdf-legend-footer");
        document.body.classList.remove("pdf-legend-footer-on");
        if (coverEl) coverEl.innerHTML = "";
        if (footEl) footEl.innerHTML = "";
        if (placement === "off") return;

        const samples = _legendSamples();
        const cut = _tassCutoffSummary();
        const hasWatch = typeof watchlist !== "undefined" && watchlist.size > 0;

        // Human-readable cutoff strings.
        const cutStrings =
          cut.mode === "byType"
            ? cut.items.map((it) => `${it.type}: ${Number(it.applied).toFixed(0)}`)
            : [
                `${
                  cut.global > 0 ? cut.global.toFixed(0) : cut.recommended != null ? cut.recommended.toFixed(0) : "—"
                }`,
              ];

        if (placement === "footer") {
          if (!footEl) return;
          const swatches = samples
            .map(
              (s) =>
                `<span class="pdf-foot-item"><span class="pdf-swatch" style="background:${s.color}"></span>${_legendEsc(
                  s.id,
                )}</span>`,
            )
            .join("");
          const cutTxt =
            cut.mode === "byType"
              ? `TASS cutoff &mdash; ${cutStrings.map(_legendEsc).join(" · ")}`
              : `TASS cutoff: ${_legendEsc(cutStrings[0])}`;
          footEl.innerHTML =
            `<div class="pdf-foot-row"><span class="pdf-foot-label">Samples:</span>${swatches}` +
            `<span class="pdf-foot-item"><span class="pdf-foot-label">${cutTxt}</span></span></div>`;
          document.body.classList.add("pdf-legend-footer-on");
          return;
        }

        // placement === "cover": full legend block.
        if (!coverEl) return;
        const sampleItems = samples
          .map(
            (s) =>
              `<span class="pdf-legend-item"><span class="pdf-swatch" style="background:${s.color}"></span>${_legendEsc(
                s.id,
              )}</span>`,
          )
          .join("");

        let cutItems;
        if (cut.mode === "byType") {
          cutItems = cut.items
            .map(
              (it) =>
                `<span class="pdf-legend-item"><b>${_legendEsc(it.type)}</b>: ${Number(it.applied).toFixed(1)}${
                  it.def != null && !isNaN(it.def)
                    ? ` <span style="color:#888">(default ${Number(it.def).toFixed(0)})</span>`
                    : ""
                }</span>`,
            )
            .join("");
        } else {
          const g = cut.global;
          cutItems = `<span class="pdf-legend-item">${
            g > 0 ? `Applied cutoff: <b>${g.toFixed(1)}</b>` : "No minimum cutoff applied"
          }${
            cut.recommended != null
              ? ` <span style="color:#888">(recommended ${cut.recommended.toFixed(1)})</span>`
              : ""
          }</span>`;
        }

        // Marker key — explains the badges/markers that appear across the report.
        const keyParts = [];
        if (hasWatch)
          keyParts.push(`<span class="pdf-legend-item"><span style="color:#f5a623">&#9733;</span> Follow-up</span>`);
        keyParts.push(
          `<span class="pdf-legend-item"><span style="color:#c62828">&#9888;</span> High consequence</span>`,
        );
        keyParts.push(
          `<span class="pdf-legend-item"><span style="display:inline-flex;width:12px;height:12px;border-radius:50%;background:#1565c0;color:#fff;font-size:8px;font-weight:700;align-items:center;justify-content:center">D</span> DNA</span>`,
        );
        keyParts.push(
          `<span class="pdf-legend-item"><span style="display:inline-flex;width:12px;height:12px;border-radius:50%;background:#6a1b9a;color:#fff;font-size:8px;font-weight:700;align-items:center;justify-content:center">R</span> RNA</span>`,
        );
        keyParts.push(
          `<span class="pdf-legend-item"><span style="font-size:8px;font-weight:700;padding:0 3px;border-radius:3px;background:#e3f2fd;color:#1565c0">Strain</span>/<span style="font-size:8px;font-weight:700;padding:0 3px;border-radius:3px;background:#f3e5f5;color:#6a1b9a">Species</span>/<span style="font-size:8px;font-weight:700;padding:0 3px;border-radius:3px;background:#e8f5e9;color:#2e7d32">Genus</span> level</span>`,
        );

        coverEl.innerHTML =
          `<div class="pdf-legend-section">` +
          `<div class="pdf-legend-h">Sample colors</div>` +
          `<div class="pdf-legend-items">${
            sampleItems || '<span class="pdf-legend-item">No visible samples</span>'
          }</div>` +
          `</div>` +
          `<div class="pdf-legend-section">` +
          `<div class="pdf-legend-h">TASS cutoff${cut.mode === "byType" ? "s (by sample type)" : ""}</div>` +
          `<div class="pdf-legend-items">${cutItems}</div>` +
          `</div>` +
          `<div class="pdf-legend-section">` +
          `<div class="pdf-legend-h">Marker key</div>` +
          `<div class="pdf-legend-items">${keyParts.join("")}</div>` +
          `</div>`;
      }

      // Tag fixed-size SVGs with a viewBox so print CSS can scale them to the
      // page width without distortion. Idempotent + harmless on screen.
      function _optimizeSvgForPrint() {
        // Inline a minimal font-face so SVG text prints correctly even in
        // offline / sandboxed PDF renderers that strip external CSS.
        const _fontFaceStyle = `
          @font-face {
            font-family: system-ui;
            src: local("system-ui"), local("-apple-system"), local("Segoe UI"),
                 local("Helvetica Neue"), local("Arial");
          }`;
        document.querySelectorAll("#content svg, #xs-print-wrap svg").forEach((svg) => {
          // Skip Leaflet's map overlay — it must keep its own coordinate system.
          if (svg.closest(".leaflet-pane") || svg.classList.contains("leaflet-zoom-animated")) return;
          const wAttr = parseFloat(svg.getAttribute("width"));
          const hAttr = parseFloat(svg.getAttribute("height"));
          // Add viewBox from width/height attrs if missing.
          if (!svg.getAttribute("viewBox") && wAttr > 0 && hAttr > 0) {
            svg.setAttribute("viewBox", `0 0 ${wAttr} ${hAttr}`);
            svg.setAttribute("preserveAspectRatio", "xMinYMin meet");
          }
          if (svg.getAttribute("viewBox")) svg.classList.add("pdf-fit");
          // Inject a <style> with font-face into each SVG so text renders
          // correctly when the SVG is rasterized by the browser print engine.
          if (!svg.querySelector("style.pdf-font")) {
            const st = document.createElementNS("http://www.w3.org/2000/svg", "style");
            st.className = "pdf-font";
            st.textContent = _fontFaceStyle + " text { font-family: system-ui, Arial, sans-serif; }";
            svg.insertBefore(st, svg.firstChild);
          }
          // Ensure text elements have an explicit font-family attribute so
          // browsers that ignore injected SVG <style> still get the right font.
          svg.querySelectorAll("text").forEach((t) => {
            if (!t.getAttribute("font-family")) {
              t.setAttribute("font-family", "system-ui, Arial, sans-serif");
            }
          });
        });
      }

      /* ── PDF section metadata ─────────────────────────────────────────────
         Per-tab grouping + plain-language explainers ("what this shows" /
         "how to read it") injected at the top of each printed pane so the
         exported report reads as an organized, annotated document rather than
         a raw dump of plots. */
      const PDF_SECTION_INFO = {
        summary: {
          group: "Overview",
          what: "Headline run metrics plus the core detection tables — sample count, total reads, percent of reads classified, unique organisms/genera, and the high-confidence calls.",
          how: "Read the KPI cards left-to-right for run scale and quality. In the detection table, each row is one organism call with its TASS confidence score; red rows are high-consequence pathogens.",
        },
        heatmap: {
          group: "Taxonomic Composition",
          what: "A samples × organisms grid in which each cell's colour encodes the selected metric (TASS score, reads, or coverage).",
          how: "Darker cells are higher values. Scan across a row to see which samples share an organism; scan down a column to profile one sample. Empty cells mean the organism was not detected there.",
        },
        tass: {
          group: "Detection Confidence",
          what: "TASS confidence scores for each detection drawn side-by-side as bars, coloured by sample.",
          how: "Taller bars are more confident detections. Compare each bar against the TASS cutoff — bars below it are weak calls. Bar colours match the sample legend.",
        },
        sunburst: {
          group: "Taxonomic Composition",
          what: "A radial view of the taxonomy: inner rings are higher ranks (kingdom→genus), outer rings are species/strains, sized by the chosen metric.",
          how: "Wider wedges represent a larger share of the metric. Trace a wedge outward to follow a lineage from genus down to strain; the centre is the root of the tree.",
        },
        coverage: {
          group: "Confidence & Coverage",
          what: "A scatter of genome coverage against TASS score (or read count); each point is one detection, coloured by sample.",
          how: "Points high and to the right are strong, well-covered detections. Outliers — e.g. many reads but low breadth — sit away from the main cluster and deserve a closer look.",
        },
        proteins: {
          group: "Virulence & Resistance",
          what: "Virulence-factor (VF) and antimicrobial-resistance (AMR) gene hits detected across the run, summarised by genus with a detailed gene table.",
          how: "Bars rank genera by gene-hit count; AMR hits flag possible drug resistance. Use the table to see individual gene names and which samples carried them.",
        },
        histogram: {
          group: "Per-Genome Coverage",
          what: "Per-contig / per-genome read-depth distributions showing how sequencing depth is spread across each reference genome.",
          how: "Even, full-width coverage indicates a genuine, well-covered genome; spiky or sparse coverage suggests partial or ambiguous mapping. Bin size varies with genome length.",
        },
        explore: {
          group: "Cross-Sample Patterns",
          what: "Free-form cross-sample charts for spotting multi-metric trends, batch effects, and sample-type patterns across the whole run.",
          how: "Pick the axes/metrics that answer your question; clusters and trends reveal samples that behave alike or stand apart.",
        },
        table: {
          group: "Detection Records",
          what: "The complete detections table behind every chart — every organism call with its metrics, in the current filtered view.",
          how: "Each row is one detection. Key columns: TASS Score (confidence), Breadth % (genome coverage), and # Reads Aligned. Red rows are high-consequence; amber-barred rows were rescued by species/genus roll-up.",
        },
        map: {
          group: "Provenance & Geography",
          what: "Geographic placement of samples that carry latitude/longitude metadata.",
          how: "Each pin is a sampling site, coloured by sample. Use it to track spread across locations for surveillance.",
        },
        runmeta: {
          group: "Provenance & Geography",
          what: "Per-sample run provenance — sample type, sequencing platform, and any uploaded metadata fields — collected in one table.",
          how: "Use it to confirm each sample's origin and sequencing context when interpreting the detections above.",
        },
      };

      function _esc(s) {
        return String(s == null ? "" : s)
          .replace(/&/g, "&amp;")
          .replace(/</g, "&lt;")
          .replace(/>/g, "&gt;");
      }

      /* Compute the concise KPI set for the PDF cover — mirrors the Summary
         tab's KPI logic but reads from the current filtered data so the cover
         reflects exactly what was exported. */
      function _pdfKpiCards() {
        const fd = typeof filteredData === "function" ? filteredData() : [];
        const samples = uniq(fd.map((r) => r["Specimen ID"] || "")).filter(Boolean);
        const totalInput = samples.reduce((s, sn) => s + (parseFloat((SAMPLE_META[sn] || {}).total_reads) || 0), 0);
        const totalOrg =
          samples.reduce((s, sn) => s + (parseFloat((SAMPLE_META[sn] || {}).total_organism_reads) || 0), 0) ||
          fd.reduce((s, r) => s + (parseFloat(r["# Reads Aligned"]) || 0), 0);
        let pctClass = totalInput > 0 ? ((totalOrg / totalInput) * 100).toFixed(1) + "%" : "N/A";
        if (pctClass === "0.0%" && totalOrg > 0) pctClass = "<0.1%";
        const uniqOrgs = new Set(fd.map((r) => r["Taxonomic ID #"] || "").filter(Boolean)).size;
        const uniqGenera = new Set(fd.map((r) => r["Genus"]).filter(Boolean)).size;
        const hcCount = new Set(
          fd.filter((r) => isTruthy(r["High Consequence"])).map((r) => r["Taxonomic ID #"] || r["Detected Organism"]),
        ).size;
        const platforms = uniq(
          samples.map((sn) => (SAMPLE_META[sn] || {}).platform).filter((p) => p && p !== "unknown"),
        );
        const cut = _tassCutoffSummary();
        let cutVal, cutSub;
        if (cut.mode === "byType" && cut.items.length) {
          const vals = cut.items.map((it) => Number(it.applied));
          const same = vals.every((v) => v === vals[0]);
          if (same) {
            cutVal = vals[0].toFixed(1);
            cutSub = cut.items.length === 1 ? cut.items[0].type : "all types";
          } else {
            cutVal = `${Math.min(...vals).toFixed(0)}–${Math.max(...vals).toFixed(0)}`;
            cutSub = "by sample type";
          }
        } else {
          cutVal = cut.global > 0 ? cut.global.toFixed(1) : cut.recommended != null ? cut.recommended.toFixed(1) : "—";
          cutSub = cut.recommended != null && cut.global <= 0 ? "recommended" : "applied filter";
        }
        const tReads = _fmtBig(totalInput);
        const aReads = _fmtBig(totalOrg);
        return [
          { label: "Samples", value: samples.length, sub: platforms.length ? platforms.join(", ") : "in report" },
          { label: "Total Reads", value: totalInput > 0 ? tReads.short : "N/A", sub: "input reads" },
          { label: "% Classified", value: pctClass, sub: aReads.short + " aligned" },
          { label: "Organisms", value: _fmtInt(uniqOrgs), sub: "unique taxa" },
          { label: "Genera", value: _fmtInt(uniqGenera), sub: "unique genera" },
          { label: "High Consequence", value: _fmtInt(hcCount), sub: "flagged", hc: true },
          { label: "TASS Cutoff", value: cutVal, sub: cutSub },
        ];
      }

      /* Build the PDF cover: KPI banner + an ordered "contents" list of the
         sections that will follow. */
      function _buildPdfCover() {
        const bannerEl = document.getElementById("pdf-kpi-banner");
        if (bannerEl) {
          bannerEl.innerHTML = _pdfKpiCards()
            .map(
              (c) =>
                `<div class="pdf-kpi${c.hc && c.value !== "0" ? " hc" : ""}">` +
                `<div class="pdf-kpi-label">${_esc(c.label)}</div>` +
                `<div class="pdf-kpi-value">${_esc(c.value)}</div>` +
                (c.sub ? `<div class="pdf-kpi-sub">${_esc(c.sub)}</div>` : "") +
                `</div>`,
            )
            .join("");
        }
        const contentsEl = document.getElementById("pdf-contents");
        if (contentsEl) {
          const tabs = _visibleReportTabs();
          const items = tabs
            .map((tab) => {
              const info = PDF_SECTION_INFO[tab] || {};
              const title = _tabLabel(tab);
              return `<li><b>${_esc(title)}</b>${
                info.group ? ` <span style="color:#8a93a3">— ${_esc(info.group)}</span>` : ""
              }</li>`;
            })
            .join("");
          contentsEl.innerHTML = `<div class="pdf-contents-h">In this report</div><ol>${items}</ol>`;
        }
      }

      // Track injected nodes so cleanup can remove them.
      let _pdfInjectedNodes = [];
      function _clearPdfSections() {
        _pdfInjectedNodes.forEach((n) => n.parentNode && n.parentNode.removeChild(n));
        _pdfInjectedNodes = [];
      }

      /* Prepend a numbered section header + "what / how" explainer to each
         visible pane so the printed output is grouped and annotated. */
      function _injectPdfSections() {
        _clearPdfSections();
        const tabs = _visibleReportTabs();
        tabs.forEach((tab, i) => {
          const pane = document.getElementById(`pane-${tab}`);
          if (!pane) return;
          const info = PDF_SECTION_INFO[tab] || {};
          const title = _tabLabel(tab);

          const head = document.createElement("div");
          head.className = "pdf-section-head";
          head.innerHTML =
            (info.group ? `<div class="pdf-section-group">${_esc(info.group)}</div>` : "") +
            `<div class="pdf-section-title"><span class="pdf-section-num">${i + 1}</span>${_esc(title)}</div>`;

          pane.insertBefore(head, pane.firstChild);
          _pdfInjectedNodes.push(head);

          if (info.what || info.how) {
            const exp = document.createElement("div");
            exp.className = "pdf-explainer";
            exp.innerHTML =
              (info.what ? `<div><div class="pdf-explainer-h">What this shows</div>${_esc(info.what)}</div>` : "") +
              (info.how ? `<div><div class="pdf-explainer-h">How to read it</div>${_esc(info.how)}</div>` : "");
            head.insertAdjacentElement("afterend", exp);
            _pdfInjectedNodes.push(exp);
          }
        });
      }

      async function _exportReportPdf() {
        _setPdfProgress("Collecting current report state…", 0, 1);
        const coverSub = document.getElementById("pdf-report-subtext");
        if (coverSub) {
          const stamp = new Date().toLocaleString();
          coverSub.textContent = `${_buildBannerSub()} · Exported ${stamp}`;
        }
        await _renderTabsForPdf();
        // Build the cover KPI banner + contents, then group/annotate each pane.
        _buildPdfCover();
        _injectPdfSections();
        _enhanceExports();
        // Scale every plot to the page (vector, no quality loss).
        _optimizeSvgForPrint();
        // Build the sample-color + TASS-cutoff legend per the chosen placement.
        const _legendPlace = (document.getElementById("report-pdf-legend-place") || {}).value || "cover";
        _buildPdfLegend(_legendPlace);
        document.body.classList.add("report-pdf-printing");
        // The print stylesheet shrinks #map-split from 560px → 320px. Leaflet
        // must recompute its container size and re-fit the data bounds at the
        // new height, otherwise it keeps the on-screen zoom and the CARTO
        // attribution control dominates the smaller printed frame.
        let _mapSettleDelay = 0;
        if (_leafletMap) {
          _leafletMap.invalidateSize();
          _fitMapToData();
          // Give the re-fit tiles a moment to load before the print dialog.
          _mapSettleDelay = 350;
        }
        const cleanup = () => {
          document.body.classList.remove("report-pdf-printing");
          document.body.classList.remove("pdf-legend-footer-on");
          _clearPdfSections();
          // Hide the cmp analysis snapshots again after printing
          const snap = document.getElementById("cmp-print-snapshots");
          if (snap) snap.style.display = "none";
          _hidePdfProgress();
          window.removeEventListener("afterprint", cleanup);
          _scheduleExportEnhance();
          // Re-render the longitudinal chart after the print layout is removed so
          // it measures the correct container width (sidebar back, normal layout).
          requestAnimationFrame(() =>
            requestAnimationFrame(() => {
              if (typeof _drawLongitudinalPlot === "function") {
                const longiWrap = document.getElementById("longi-chart-wrap");
                if (longiWrap && longiWrap.querySelector("svg")) _drawLongitudinalPlot();
              }
            }),
          );
        };
        window.addEventListener("afterprint", cleanup);
        setTimeout(() => {
          _setPdfProgress("Opening print dialog…", 1, 1, "Choose Save as PDF in your browser dialog.");
          window.print();
          setTimeout(() => {
            if (document.body.classList.contains("report-pdf-printing")) cleanup();
          }, 1500);
        }, 120 + _mapSettleDelay);
      }

