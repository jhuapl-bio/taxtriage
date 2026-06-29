      /* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB: HISTOGRAM          (data-tab="histogram"  —  hidden by default)
       -     Per-contig read / depth distribution, one panel per organism,
       -     wrapped in an IIFE that exposes its draw/clear callbacks. Driven
       -     by genome-coverage histogram data attached to organism rows.
═══════════════════════════════════════════════════════════════════════════ */
      (function () {
        function _getOrgKey(cd) {
          return `${cd.organism}||${cd.taxon_id}`;
        }

        // Set of `${Specimen ID}||${Taxonomic ID #}` for visible detections that
        // pass the TASS cutoff at their OWN level. Used to gate the histogram so
        // organisms below the cutoff (including roll-up-rescued ones) can't be
        // selected or charted.
        function _histPassingTaxa() {
          const set = new Set();
          filteredData().forEach((r) => {
            const pi = typeof rowPassInfo === "function" ? rowPassInfo(r) : null;
            // Keep rows with no numeric strain TASS (nothing to gate on), else
            // require a strict own-level pass.
            if (pi && !isNaN(pi.strain) && !pi.strainPass) return;
            set.add(`${r["Specimen ID"]}||${r["Taxonomic ID #"]}`);
          });
          return set;
        }

        function _buildSelectors() {
          const orgSel = document.getElementById("hist-org-sel");
          const sampleSel = document.getElementById("hist-sample-sel");
          const searchEl = document.getElementById("hist-org-search");
          const searchClear = document.getElementById("hist-org-search-clear");
          if (!orgSel || !sampleSel) return;

          function _getTextFilter() {
            const txt = (document.getElementById("filter-text") || { value: "" }).value.trim();
            const ic = (document.getElementById("filter-ic") || { checked: true }).checked;
            if (!txt) return null;
            try {
              return new RegExp(txt, ic ? "i" : "");
            } catch (err) {
              return null;
            }
          }

          function _getSearchTerm() {
            return (searchEl ? searchEl.value : "").trim().toLowerCase();
          }

          function _buildOrgKeys() {
            const orgKeys = {};
            const rx = _getTextFilter();
            const q = _getSearchTerm();
            // Only organisms whose detection clears the TASS cutoff at its own
            // level are selectable — a below-cutoff organism (even one kept by
            // species/genus roll-up) must not be pickable here.
            const visibleTaxa = _histPassingTaxa();
            CONTIG_DATA.forEach((cd) => {
              const orgName = cd.organism || "";
              const sample = cd.sample || "";
              if (sampleHidden[sample]) return;
              if (!visibleTaxa.has(`${sample}||${cd.taxon_id}`)) return;
              if (!cd.contigs || !cd.contigs.length) return;
              if (rx && !rx.test(orgName) && !rx.test(sample)) return;
              if (q && !orgName.toLowerCase().includes(q)) return;
              const k = _getOrgKey(cd);
              if (!orgKeys[k]) orgKeys[k] = cd.organism;
            });
            return orgKeys;
          }

          function _refreshOrgs() {
            const prev = orgSel.value;
            const orgKeys = _buildOrgKeys();
            orgSel.innerHTML = "";
            Object.entries(orgKeys)
              .sort((a, b) => a[1].localeCompare(b[1]))
              .forEach(([k, name]) => {
                const opt = document.createElement("option");
                opt.value = k;
                opt.textContent = name;
                orgSel.appendChild(opt);
              });
            if (prev && Array.from(orgSel.options).some((o) => o.value === prev)) {
              orgSel.value = prev;
            }
          }

          function _refreshSamples() {
            const orgKey = orgSel.value;
            sampleSel.innerHTML = '<option value="">All samples</option>';
            const rx = _getTextFilter();
            const visibleTaxa = _histPassingTaxa();
            CONTIG_DATA.filter(
              (cd) =>
                _getOrgKey(cd) === orgKey &&
                !sampleHidden[cd.sample] &&
                visibleTaxa.has(`${cd.sample}||${cd.taxon_id}`) &&
                cd.contigs &&
                cd.contigs.length,
            ).forEach((cd) => {
              if (rx && !rx.test(cd.organism || "") && !rx.test(cd.sample || "")) return;
              const opt = document.createElement("option");
              opt.value = cd.sample;
              opt.textContent = cd.sample;
              sampleSel.appendChild(opt);
            });
          }

          _refreshOrgs();
          _refreshSamples();

          if (!orgSel._wired) {
            orgSel._wired = true;
            orgSel.addEventListener("change", () => {
              _refreshSamples();
              window.drawHistogram();
            });
          }
          if (searchEl && !searchEl._wired) {
            searchEl._wired = true;
            searchEl.addEventListener("input", () => {
              _refreshOrgs();
              _refreshSamples();
              window.drawHistogram();
            });
          }
          if (searchClear && !searchClear._wired) {
            searchClear._wired = true;
            searchClear.addEventListener("click", () => {
              if (searchEl) searchEl.value = "";
              _refreshOrgs();
              _refreshSamples();
              window.drawHistogram();
            });
          }
          ["filter-text", "filter-ic", "filter-mc", "filter-min", "filter-min-range", "filter-rollup"].forEach((id) => {
            const el = document.getElementById(id);
            if (el && !el._histWired) {
              el._histWired = true;
              const _onHistFilterChange = () => {
                _invalidateFilterCache();
                _refreshOrgs();
                _refreshSamples();
                window.drawHistogram();
              };
              el.addEventListener("input", _onHistFilterChange);
              el.addEventListener("change", _onHistFilterChange);
            }
          });
        }

        // Persistent filter state across drawHistogram calls
        let histSelectedContig = null;
        let _histLastOrgKey = null;

        window.drawHistogram = function () {
          const chartWrap = document.getElementById("hist-chart-wrap");
          const depthWrap = document.getElementById("hist-depth-chart");
          const noData = document.getElementById("hist-no-data");
          if (!chartWrap) return;

          chartWrap.innerHTML = "";
          if (depthWrap) depthWrap.innerHTML = "";
          const breadthChartEl = document.getElementById("hist-breadth-chart");
          if (breadthChartEl) breadthChartEl.innerHTML = "";

          const orgSel = document.getElementById("hist-org-sel");
          const sampleSel = document.getElementById("hist-sample-sel");
          const metricSel = document.getElementById("hist-metric-sel");
          const top20El = document.getElementById("hist-top20");
          if (!orgSel) return;

          const orgKey = orgSel.value;
          // Reset per-contig filter when the organism changes
          if (orgKey !== _histLastOrgKey) {
            histSelectedContig = null;
            _histLastOrgKey = orgKey;
          }
          const sample = sampleSel ? sampleSel.value : "";
          const metric = metricSel ? metricSel.value : "reads";
          const limitTo = top20El ? parseInt(top20El.value, 10) || 15 : 9999;

          if (!orgKey) {
            noData.style.display = "block";
            chartWrap.style.display = "none";
            return;
          }

          const _fdTaxaKeys = _histPassingTaxa();
          let entries = CONTIG_DATA.filter(
            (cd) =>
              _getOrgKey(cd) === orgKey && !sampleHidden[cd.sample] && _fdTaxaKeys.has(`${cd.sample}||${cd.taxon_id}`),
          );
          if (sample) entries = entries.filter((cd) => cd.sample === sample);
          // Sort entries to match the user-defined sample order from the main panel
          if (_sampleOrder.length) {
            const _idx = Object.fromEntries(_sampleOrder.map((id, i) => [id, i]));
            entries = [...entries].sort((a, b) => {
              const ia = _idx[a.sample] !== undefined ? _idx[a.sample] : 9999;
              const ib = _idx[b.sample] !== undefined ? _idx[b.sample] : 9999;
              return ia - ib;
            });
          }
          // Always show the sample name in chart titles when the dropdown is
          // on "All samples" — even if just one entry matches. Previously this
          // required `entries.length > 1`, which meant uploaded samples that
          // were the sole match for an organism rendered with no name at all.
          const showAll = !sample;

          if (!entries.length || entries.every((cd) => !cd.contigs || !cd.contigs.length)) {
            noData.style.display = "block";
            chartWrap.style.display = "none";
            return;
          }
          noData.style.display = "none";
          chartWrap.style.display = "block";

          const METRIC_LABELS = {
            reads: "Reads Aligned",
            mean_depth: "Mean Depth",
            coverage: "Coverage %",
            covered_bases: "Covered Bases",
          };
          const mLabel = METRIC_LABELS[metric] || metric;

          // ── Shared constants ───────────────────────────────────────────────────
          const barH = 18,
            mT_reads = 24,
            mB_reads = 30,
            mL_reads = 160,
            mR_reads = 36;
          const mT_brd = 28,
            mB_brd = 36,
            iH_brd_min = 110;
          const mT_dep = 22,
            mB_dep = 28,
            mExtra_dep = 14;
          const BRD_SVG_MIN = mT_brd + iH_brd_min + mB_brd; // 174
          const DEP_SVG_MIN = 130 + mExtra_dep; // 144

          // ── Per-row heights: all three columns sync to the tallest ─────────────
          const rowHeights = entries.map((cd) => {
            const nC = (cd.contigs || []).length > 0 ? Math.min(cd.contigs.length, limitTo) : 1;
            return Math.max(mT_reads + nC * barH + mB_reads, BRD_SVG_MIN, DEP_SVG_MIN);
          });

          function _rowSep() {
            const d = document.createElement("div");
            d.style.cssText = "height:1px;background:#d0d5e0;margin:5px 0;";
            return d;
          }
          function fmtLen(bp) {
            if (bp >= 1e6) return (bp / 1e6).toFixed(2) + " Mbp";
            if (bp >= 1e3) return (bp / 1e3).toFixed(1) + " Kbp";
            return bp + " bp";
          }
          // Returns the contig index (into bh.contig_names) that owns a given bin
          function _binContigIdx(binIdx, breaks) {
            let ci = 0;
            for (let i = breaks.length - 1; i >= 0; i--) {
              if (binIdx >= breaks[i]) {
                ci = i;
                break;
              }
            }
            return ci;
          }

          const hasBreadth = entries.some(
            (cd) =>
              cd.breadth_histogram && Array.isArray(cd.breadth_histogram.bins) && cd.breadth_histogram.bins.length > 0,
          );
          const hasDepth = entries.some(
            (cd) => cd.depth_histogram && Object.values(cd.depth_histogram).some((v) => v > 0),
          );

          // ── Per-entry breadth metadata (populated after breadth renders, used by bar hover) ──
          const _brdMeta = {};
          function _showBrdLabel(ei, cName) {
            const lbl = document.getElementById(`brd-hlbl-${ei}`);
            if (!lbl) return;
            const m = _brdMeta[ei];
            if (!m || !cName) {
              lbl.textContent = "";
              return;
            }
            const ci = m.contigNames.indexOf(cName);
            if (ci < 0) {
              lbl.textContent = "";
              return;
            }
            const brkStart = m.breaks[ci] || 0;
            const brkEnd = ci + 1 < m.breaks.length ? m.breaks[ci + 1] : m.nBins;
            const midBin = (brkStart + brkEnd) / 2;
            const xPos = m.xSc(midBin);
            if (xPos < -2 || xPos > m.iW + 2) {
              lbl.textContent = "";
              return;
            }
            lbl.setAttribute("x", Math.max(2, Math.min(xPos, m.iW - 70)));
            lbl.textContent = cName.length > 16 ? cName.slice(0, 15) + "\u2026" : cName;
          }
          function _clearBrdLabel(ei) {
            const lbl = document.getElementById(`brd-hlbl-${ei}`);
            if (lbl) lbl.textContent = "";
          }
          function _showBrdHighlight(ei, cName) {
            const rect = document.getElementById(`brd-hover-region-${ei}`);
            if (!rect) return;
            const m = _brdMeta[ei];
            if (!m || !cName) {
              rect.setAttribute("display", "none");
              return;
            }
            const ci = m.contigNames.indexOf(cName);
            if (ci < 0) {
              rect.setAttribute("display", "none");
              return;
            }
            const brkStart = m.breaks[ci] || 0;
            const brkEnd = ci + 1 < m.breaks.length ? m.breaks[ci + 1] : m.nBins;
            const x1 = m.xSc(brkStart);
            const x2 = m.xSc(brkEnd);
            if (x2 <= 0 || x1 >= m.iW) {
              rect.setAttribute("display", "none");
              return;
            }
            rect.setAttribute("x", Math.max(0, x1));
            rect.setAttribute("width", Math.min(m.iW, x2) - Math.max(0, x1));
            rect.setAttribute("display", "");
          }
          function _clearBrdHighlight(ei) {
            const rect = document.getElementById(`brd-hover-region-${ei}`);
            if (rect) rect.setAttribute("display", "none");
          }

          // ── Reads bar chart — clickable, filters breadth & depth ──────────────
          entries.forEach((cd, ei) => {
            if (!cd.contigs || !cd.contigs.length) return;
            const contigs = [...cd.contigs]
              .sort((a, b) => {
                const va = metric === "coverage" ? b[metric] * 100 : b[metric];
                const vb = metric === "coverage" ? a[metric] * 100 : a[metric];
                return va - vb;
              })
              .slice(0, limitTo);

            const rowH = rowHeights[ei];
            const W = Math.max(260, chartWrap.clientWidth || 320) - 4;
            const iW = W - mL_reads - mR_reads;
            const nTicks = Math.max(2, Math.min(8, Math.floor(iW / 55)));

            const naturalBarRegionH = contigs.length * barH;
            const availForBars = rowH - mT_reads - mB_reads;
            const barOffset = Math.max(0, Math.floor((availForBars - naturalBarRegionH) / 2));
            const barsTop = mT_reads + barOffset;

            const valFn = (c) => (metric === "coverage" ? c[metric] * 100 : c[metric]);
            const xMax = d3.max(contigs, valFn) || 1;
            const xScale = d3.scaleLinear().domain([0, xMax]).range([0, iW]).nice();
            const yScale = d3
              .scaleBand()
              .domain(contigs.map((c) => c.name))
              .range([barsTop, barsTop + naturalBarRegionH])
              .padding(0.1);
            const colorFn = d3.scaleSequential(d3.interpolateBlues).domain([0, xMax]);

            const svgWrap = document.createElement("div");
            chartWrap.appendChild(svgWrap);
            // NO title — organism is shown once in breadth panel above
            const svg = d3.select(svgWrap).append("svg").attr("width", W).attr("height", rowH);
            const g = svg.append("g").attr("transform", `translate(${mL_reads},0)`);

            // Bars — clickable to filter; dim unselected
            g.selectAll("rect.bar-rect")
              .data(contigs)
              .enter()
              .append("rect")
              .attr("class", "bar-rect")
              .attr("x", 0)
              .attr("y", (c) => yScale(c.name))
              .attr("height", yScale.bandwidth())
              .attr("width", (c) => Math.max(2, xScale(valFn(c))))
              .attr("fill", (c) => colorFn(valFn(c)))
              .attr("rx", 2)
              .attr("opacity", (c) => (!histSelectedContig || histSelectedContig === c.name ? 1 : 0.22))
              .attr("stroke", (c) => (histSelectedContig === c.name ? "#e65100" : null))
              .attr("stroke-width", (c) => (histSelectedContig === c.name ? 1.8 : null))
              .style("cursor", "pointer")
              .on("click", function (ev, c) {
                histSelectedContig = histSelectedContig === c.name ? null : c.name;
                window.drawHistogram();
              })
              .on("mouseover", function (ev, c) {
                showTip(
                  `<b>${c.name}</b><br>` +
                    `<span style="color:#aaa;font-size:0.9em">Click to filter breadth &amp; depth</span><br>` +
                    `Reads: ${c.reads.toLocaleString()}<br>Mean Depth: ${c.mean_depth}x<br>` +
                    `Coverage: ${(c.coverage * 100).toFixed(1)}%<br>` +
                    `Covered: ${c.covered_bases.toLocaleString()} / ${c.length.toLocaleString()} bp`,
                  ev,
                );
                _showBrdLabel(ei, c.name);
                _showBrdHighlight(ei, c.name);
              })
              .on("mousemove", moveTip)
              .on("mouseout", function () {
                hideTip();
                _clearBrdLabel(ei);
                _clearBrdHighlight(ei);
              });

            g.selectAll("text.val")
              .data(contigs)
              .enter()
              .append("text")
              .attr("class", "val")
              .attr("x", (c) => Math.max(2, xScale(valFn(c))) + 4)
              .attr("y", (c) => yScale(c.name) + yScale.bandwidth() / 2 + 4)
              .attr("font-size", 9)
              .attr("fill", "#555")
              .attr("opacity", (c) => (!histSelectedContig || histSelectedContig === c.name ? 1 : 0.3))
              .text((c) => {
                const v = valFn(c);
                if (metric === "coverage") return v.toFixed(1) + "%";
                if (metric === "mean_depth") return v.toFixed(1) + "x";
                return v >= 1000 ? `${(v / 1000).toFixed(1)}k` : v;
              });

            g.append("g")
              .call(d3.axisLeft(yScale).tickSize(0))
              .selectAll("text")
              .attr("font-size", 9)
              .text((d) => (d.length > 22 ? d.slice(0, 21) + "…" : d));

            g.append("g")
              .attr("transform", `translate(0,${barsTop + naturalBarRegionH})`)
              .call(d3.axisBottom(xScale).ticks(nTicks))
              .append("text")
              .attr("x", iW / 2)
              .attr("y", 26)
              .attr("text-anchor", "middle")
              .attr("fill", "#666")
              .attr("font-size", 10)
              .text(mLabel);

            if (ei < entries.length - 1) chartWrap.appendChild(_rowSep());
          });

          // ── Positional breadth histogram (samtools-style) ─────────────────────
          // X = genomic bin position, Y = % of that bin covered (0–100 %)
          // When a contig is selected, x-axis zooms to that contig's bin range.
          (function () {
            const breadthWrap = document.getElementById("hist-breadth-chart");
            const breadthSection = document.getElementById("hist-breadth-wrap");
            if (!breadthWrap) return;
            breadthWrap.innerHTML = "";

            if (!hasBreadth) {
              if (breadthSection) breadthSection.style.display = "none";
              return;
            }
            if (breadthSection) breadthSection.style.display = "";

            entries.forEach((cd, ei) => {
              const rowH = rowHeights[ei];
              const bh = cd.breadth_histogram;
              const hasBh = bh && Array.isArray(bh.bins) && bh.bins.length > 0;

              if (!hasBh) {
                const spacer = document.createElement("div");
                spacer.style.cssText = `height:${rowH}px;display:flex;align-items:center;justify-content:center;color:#ccc;font-size:9px;`;
                spacer.textContent = "—";
                breadthWrap.appendChild(spacer);
                if (ei < entries.length - 1) breadthWrap.appendChild(_rowSep());
                return;
              }

              const bins = bh.bins;
              const binSz = bh.bin_size || 2000;
              const totLen = bh.total_len || bins.length * binSz;
              const breaks = bh.breaks || [0];
              const contigNames = Array.isArray(bh.contig_names) ? bh.contig_names : [];
              const sLbl = showAll ? ` — ${cd.sample}` : "";
              const nBins = bins.length;
              const maxPct = 100;

              // Top-N set — build filtered bins/breaks/names restricted to visible contigs
              const topContigSet = new Set((cd.contigs || []).slice(0, limitTo).map((c) => c.name));

              // Remap breadth histogram to only the top-N contigs so grey sections disappear
              let filteredBins = bins,
                filteredBreaks = breaks,
                filteredContigNames = contigNames;
              if (topContigSet.size > 0 && contigNames.length > 0) {
                const remappedBins = [];
                const remappedBreaks = [];
                const remappedNames = [];
                contigNames.forEach((cName, ci) => {
                  if (!topContigSet.has(cName)) return;
                  const bStart = breaks[ci] || 0;
                  const bEnd = ci + 1 < breaks.length ? breaks[ci + 1] : nBins;
                  remappedBreaks.push(remappedBins.length);
                  remappedNames.push(cName);
                  for (let b = bStart; b < bEnd; b++) remappedBins.push(bins[b]);
                });
                if (remappedBins.length > 0) {
                  filteredBins = remappedBins;
                  filteredBreaks = remappedBreaks;
                  filteredContigNames = remappedNames;
                }
              }
              const fNBins = filteredBins.length;
              const avgPct = fNBins ? (filteredBins.reduce((a, b) => a + b, 0) / fNBins).toFixed(1) : 0;
              const maxBin = fNBins ? Math.max(...filteredBins).toFixed(1) : 0;

              // If histSelectedContig is no longer in top-N, clear it
              if (
                histSelectedContig &&
                filteredContigNames.length &&
                !filteredContigNames.includes(histSelectedContig)
              ) {
                histSelectedContig = null;
              }

              // Selected contig bin range (within filtered bins)
              const selIdx =
                histSelectedContig && filteredContigNames.length ? filteredContigNames.indexOf(histSelectedContig) : -1;
              const selBinStart = selIdx >= 0 ? filteredBreaks[selIdx] || 0 : 0;
              const selBinEnd =
                selIdx >= 0 ? (selIdx + 1 < filteredBreaks.length ? filteredBreaks[selIdx + 1] : fNBins) : fNBins;

              // X domain: zoom to selected contig when filtered
              const xDomStart = selIdx >= 0 ? selBinStart : 0;
              const xDomEnd = selIdx >= 0 ? selBinEnd : fNBins;
              const visRange = xDomEnd - xDomStart;

              const panelW = Math.max(200, (breadthWrap.clientWidth || 400) - 4);
              const mL = 42,
                mR = 8;
              const iW = panelW - mL - mR;
              const iH = Math.max(iH_brd_min, rowH - mT_brd - mB_brd);

              const svg = d3
                .select(breadthWrap)
                .append("svg")
                .attr("width", panelW)
                .attr("height", rowH)
                .style("display", "block")
                .style("margin-bottom", "0");

              svg
                .append("text")
                .attr("x", mL + iW / 2)
                .attr("y", 14)
                .attr("text-anchor", "middle")
                .attr("font-size", 10)
                .attr("font-weight", "600")
                .attr("fill", "#333")
                .text(`${cd.organism}${sLbl}`);

              const g = svg.append("g").attr("transform", `translate(${mL},${mT_brd})`);
              const xSc = d3.scaleLinear().domain([xDomStart, xDomEnd]).range([0, iW]);
              const ySc = d3.scaleLinear().domain([0, maxPct]).range([iH, 0]);

              // Store metadata so bar-chart hover can drive this label
              _brdMeta[ei] = { xSc, contigNames: filteredContigNames, breaks: filteredBreaks, nBins: fNBins, iW };

              // ClipPath so zoomed bars don't bleed outside the panel
              const clipId = `brd-clip-${ei}`;
              svg
                .append("defs")
                .append("clipPath")
                .attr("id", clipId)
                .append("rect")
                .attr("width", iW)
                .attr("height", iH);

              g.append("rect")
                .attr("x", 0)
                .attr("y", 0)
                .attr("width", iW)
                .attr("height", iH)
                .attr("fill", "#f8f9fa")
                .attr("stroke", "#ddd");

              [25, 50, 75, 100].forEach((pct) => {
                const y = ySc(pct);
                g.append("line")
                  .attr("x1", 0)
                  .attr("x2", iW)
                  .attr("y1", y)
                  .attr("y2", y)
                  .attr("stroke", pct === 100 ? "#bbb" : "#e0e0e0")
                  .attr("stroke-dasharray", pct === 100 ? "none" : "3,3");
              });

              const binPx = Math.max(1, iW / visRange);
              const colorFn = d3.scaleSequential(d3.interpolateBlues).domain([0, 100]);

              // Bind as {v, i} objects so index is always reliable
              const binsData = filteredBins.map((v, i) => ({ v, i }));

              // ── Hover highlight band (sits behind bars, driven by hover in either panel) ──
              g.append("rect")
                .attr("id", `brd-hover-region-${ei}`)
                .attr("x", 0)
                .attr("y", 0)
                .attr("width", 0)
                .attr("height", iH)
                .attr("fill", "rgba(230,81,0,0.13)")
                .attr("stroke", "#e65100")
                .attr("stroke-width", 0.8)
                .attr("display", "none")
                .style("pointer-events", "none");

              const barsG = g.append("g").attr("clip-path", `url(#${clipId})`);
              barsG
                .selectAll("rect.bh-bar")
                .data(binsData)
                .enter()
                .append("rect")
                .attr("class", "bh-bar")
                .attr("x", ({ i }) => xSc(i))
                .attr("y", ({ v }) => ySc(Math.min(maxPct, v)))
                .attr("width", binPx + 0.5)
                .attr("height", ({ v }) => iH - ySc(Math.min(maxPct, v)))
                .attr("fill", ({ v, i }) => {
                  return v > 0 ? colorFn(v) : "#ececec";
                })
                .attr("opacity", ({ i }) => (selIdx < 0 ? 1 : i >= selBinStart && i < selBinEnd ? 1 : 0.15))
                .attr("shape-rendering", "crispEdges")
                .style("cursor", "pointer")
                .on("click", function (ev, { v, i }) {
                  const ci = _binContigIdx(i, filteredBreaks);
                  const cName = filteredContigNames[ci] || null;
                  if (cName) {
                    histSelectedContig = histSelectedContig === cName ? null : cName;
                    window.drawHistogram();
                  }
                })
                .on("mouseover", function (ev, { v, i }) {
                  const ci = _binContigIdx(i, filteredBreaks);
                  const cName = filteredContigNames[ci] || "";
                  const posStart = i * binSz;
                  const posEnd = Math.min(totLen, posStart + binSz);
                  showTip(
                    `<b>${cName || `Bin ${i + 1}`}</b><br>` +
                      `<span style="color:#aaa;font-size:0.9em">Click to filter</span><br>` +
                      `Pos: ${fmtLen(posStart)} \u2013 ${fmtLen(posEnd)}<br>` +
                      `Covered: <b>${v.toFixed(1)}%</b>`,
                    ev,
                  );
                  d3.select(this).attr("stroke", "#333").attr("stroke-width", 0.8);
                  _showBrdLabel(ei, cName);
                  _showBrdHighlight(ei, cName);
                })
                .on("mousemove", moveTip)
                .on("mouseout", function () {
                  hideTip();
                  d3.select(this).attr("stroke", null);
                  _clearBrdLabel(ei);
                  _clearBrdHighlight(ei);
                });

              // ── Hover label: one floating text updated dynamically ─────────────
              g.append("text")
                .attr("id", `brd-hlbl-${ei}`)
                .attr("y", 9)
                .attr("font-size", 7)
                .attr("font-weight", "600")
                .attr("fill", "#333")
                .style("pointer-events", "none")
                .text("");

              // ── Contig break lines (only when NOT filtered — zoomed view has 1 contig) ──
              if (selIdx < 0) {
                filteredBreaks.slice(1).forEach((brkIdx) => {
                  const x = xSc(brkIdx);
                  if (x <= 0 || x >= iW) return;
                  g.append("line")
                    .attr("x1", x)
                    .attr("x2", x)
                    .attr("y1", 0)
                    .attr("y2", iH)
                    .attr("stroke", "#aaa")
                    .attr("stroke-width", 0.7)
                    .attr("stroke-dasharray", "3,2")
                    .attr("opacity", 0.7);
                });
              }

              // ── Orange border when a contig is selected (zoomed, full width) ───
              if (selIdx >= 0) {
                g.append("rect")
                  .attr("x", 0)
                  .attr("y", 0)
                  .attr("width", iW)
                  .attr("height", iH)
                  .attr("fill", "none")
                  .attr("stroke", "#e65100")
                  .attr("stroke-width", 2)
                  .attr("rx", 1);
              }

              // ── X axis ────────────────────────────────────────────────────────
              const _xTickOffset = selIdx >= 0 ? xDomStart : 0;
              const _xRefLen = selIdx >= 0 ? (selBinEnd - selBinStart) * binSz : totLen;
              const xAxis = d3
                .axisBottom(xSc)
                .ticks(Math.min(10, visRange))
                .tickFormat((d) => {
                  const bp = (d - _xTickOffset) * binSz;
                  if (_xRefLen >= 1e6) return (bp / 1e6).toFixed(1) + "M";
                  if (_xRefLen >= 1e3) return (bp / 1e3).toFixed(0) + "K";
                  return bp;
                });
              g.append("g").attr("transform", `translate(0,${iH})`).call(xAxis).selectAll("text").attr("font-size", 8);

              const xLbl =
                selIdx >= 0
                  ? `${filteredContigNames[selIdx] || "Contig"}  (bin\u202f=\u202f${fmtLen(
                      binSz,
                    )}, len\u202f\u2248\u202f${fmtLen(visRange * binSz)})`
                  : `Genomic position  (bin\u202f=\u202f${fmtLen(binSz)}, genome\u202f=\u202f${fmtLen(totLen)})`;
              g.append("text")
                .attr("x", iW / 2)
                .attr("y", iH + mB_brd - 6)
                .attr("text-anchor", "middle")
                .attr("font-size", 9)
                .attr("fill", selIdx >= 0 ? "#e65100" : "#666")
                .text(xLbl);

              g.append("g")
                .call(
                  d3
                    .axisLeft(ySc)
                    .ticks(5)
                    .tickFormat((d) => d + "%"),
                )
                .selectAll("text")
                .attr("font-size", 8);
              svg
                .append("text")
                .attr("transform", "rotate(-90)")
                .attr("x", -(mT_brd + iH / 2))
                .attr("y", 10)
                .attr("text-anchor", "middle")
                .attr("font-size", 9)
                .attr("fill", "#555")
                .text("% covered");

              const subLine =
                selIdx >= 0
                  ? `Filtered: ${filteredContigNames[selIdx] || ""} \u00b7 max=${maxBin}% \u00b7 mean=${avgPct}%`
                  : `bin=${fmtLen(binSz)} \u00b7 genome=${fmtLen(
                      totLen,
                    )} \u00b7 max=${maxBin}% \u00b7 mean=${avgPct}%` +
                    (filteredBreaks.length > 1 ? ` \u00b7 ${filteredBreaks.length} contigs` : "");
              svg
                .append("text")
                .attr("x", mL + iW / 2)
                .attr("y", mT_brd - 4)
                .attr("text-anchor", "middle")
                .attr("font-size", 8)
                .attr("fill", "#888")
                .text(subLine);

              if (ei < entries.length - 1) breadthWrap.appendChild(_rowSep());
            });
          })();

          // ── Depth histogram — per-contig when selected, else organism-level ─────
          if (!depthWrap) return;
          if (!hasDepth) return;

          const buckets = ["0x", "1-5x", "5-10x", "10-50x", ">50x"];
          const bColors = ["#e0e0e0", "#c6dbef", "#6baed6", "#2171b5", "#08306b"];

          entries.forEach((cd, ei) => {
            const rowH = rowHeights[ei];

            // Prefer per-contig depth if the selected contig has its own histogram
            const selContig = histSelectedContig ? (cd.contigs || []).find((c) => c.name === histSelectedContig) : null;
            const hist = selContig && selContig.depth_histogram ? selContig.depth_histogram : cd.depth_histogram;
            const isFiltered = !!(selContig && selContig.depth_histogram);
            const isAggregate = !!(histSelectedContig && selContig && !selContig.depth_histogram);
            const hasDh = hist && Object.values(hist).some((v) => v > 0);

            if (!hasDh) {
              const spacer = document.createElement("div");
              spacer.style.cssText = `height:${rowH}px;display:flex;align-items:center;justify-content:center;color:#ccc;font-size:9px;`;
              spacer.textContent = "—";
              depthWrap.appendChild(spacer);
              if (ei < entries.length - 1) depthWrap.appendChild(_rowSep());
              return;
            }

            const total = buckets.reduce((s, b) => s + (hist[b] || 0), 0) || 1;
            const sLbl = showAll ? ` — ${cd.sample}` : "";
            let depthTitle;
            if (isFiltered) {
              const n = histSelectedContig;
              depthTitle = n.length > 13 ? n.slice(0, 12) + "…" : n;
            } else if (isAggregate) {
              depthTitle = `Depth (org.)${sLbl}`;
            } else {
              depthTitle = `Depth dist.${sLbl}`;
            }

            const dW = Math.max(160, (depthWrap.clientWidth || 200) - 4);
            const mL = 50,
              mT = mT_dep,
              mB = mB_dep;
            const iW2 = dW - mL - 8;
            const iH2 = Math.max(60, rowH - mExtra_dep - mT - mB);

            const svg2 = d3.select(depthWrap).append("svg").attr("width", dW).attr("height", rowH);
            svg2
              .append("text")
              .attr("x", dW / 2)
              .attr("y", 12)
              .attr("text-anchor", "middle")
              .attr("font-size", 9)
              .attr("fill", isFiltered ? "#e65100" : "#555")
              .text(depthTitle);

            const g2 = svg2.append("g").attr("transform", `translate(${mL},${mT})`);
            const xS = d3.scaleBand().domain(buckets).range([0, iW2]).padding(0.1);
            const yS = d3
              .scaleLinear()
              .domain([0, d3.max(buckets, (b) => hist[b] || 0)])
              .range([iH2, 0])
              .nice();

            g2.selectAll("rect")
              .data(buckets)
              .enter()
              .append("rect")
              .attr("x", (b) => xS(b))
              .attr("y", (b) => yS(hist[b] || 0))
              .attr("width", xS.bandwidth())
              .attr("height", (b) => iH2 - yS(hist[b] || 0))
              .attr("fill", (b, i) => bColors[i])
              .attr("rx", 2)
              .on("mouseover", (ev, b) => {
                const cnt = hist[b] || 0;
                showTip(`<b>${b}</b><br>${cnt.toLocaleString()} bases (${((cnt / total) * 100).toFixed(1)}%)`, ev);
              })
              .on("mousemove", moveTip)
              .on("mouseout", hideTip);

            g2.append("g")
              .attr("transform", `translate(0,${iH2})`)
              .call(d3.axisBottom(xS).tickSize(0))
              .selectAll("text")
              .attr("font-size", 8);
            g2.append("g")
              .call(d3.axisLeft(yS).ticks(3).tickFormat(d3.format(".2s")))
              .selectAll("text")
              .attr("font-size", 8);

            if (ei < entries.length - 1) depthWrap.appendChild(_rowSep());
          });

          // ── Filter badge ────────────────────────────────────────────────────────
          (function () {
            const badge = document.getElementById("hist-filter-badge");
            const nameEl = document.getElementById("hist-filter-name");
            const clearEl = document.getElementById("hist-filter-clear");
            if (!badge) return;
            if (histSelectedContig) {
              badge.style.display = "flex";
              if (nameEl) nameEl.textContent = histSelectedContig;
              if (clearEl)
                clearEl.onclick = function (e) {
                  e.preventDefault();
                  histSelectedContig = null;
                  window.drawHistogram();
                };
            } else {
              badge.style.display = "none";
            }
          })();
        };

        // Wire histogram controls after DOM is ready (scripts are at end of body)
        _buildSelectors();
        ["hist-sample-sel", "hist-metric-sel"].forEach((id) => {
          const el = document.getElementById(id);
          if (el) el.addEventListener("change", window.drawHistogram);
        });
        const top20El = document.getElementById("hist-top20");
        if (top20El) top20El.addEventListener("input", window.drawHistogram);

        // Expose selector reset for file-upload integration
        window._resetHistSelectors = function () {
          _buildSelectors();
          window.drawHistogram && window.drawHistogram();
        };
      })();

      /* ── Table pagination state ─────────────────────────────────────────────── */
      let _tblPage = 0; // 0-based current page index

      function _tblPageSize() {
        const sel = document.getElementById("tbl-page-size");
        return sel ? parseInt(sel.value) || 0 : 100;
      }

      /* Reset to page 0 — call whenever filters/sort change so the user
         doesn't land on an empty page after narrowing results. */
      function _tblResetPage() {
        _tblPage = 0;
      }

      const _tblPinned = new Set();
      function _tblRowKey(row) {
        return `${row["Specimen ID"]}||${row["Detected Organism"]}||${row["Taxonomic ID #"]}`;
      }
      function _renderOneRow(row) {
        // Novelty-only indicator row: sample has no passing detections but carries
        // novelty genus evidence from translated-search / closed-set compare.
        if (row.__noveltyOnly) {
          const itr = document.createElement("tr");
          itr.className = "novelty-only-row";
          const itd = document.createElement("td");
          itd.colSpan = visibleCols.length || 1;
          itd.innerHTML = _noveltyOnlyMessageHTML(row["Specimen ID"], row.__novelty);
          const _s = row.__novelty;
          const _msg = itd.querySelector(".novelty-only-msg");
          if (_msg) {
            _msg.addEventListener("mouseover", (ev) => showTip(_noveltyOnlyTip(row["Specimen ID"], _s), ev));
            _msg.addEventListener("mousemove", moveTip);
            _msg.addEventListener("mouseout", hideTip);
          }
          const _lnk = itd.querySelector(".novelty-only-link");
          if (_lnk)
            _lnk.addEventListener("click", (e) => {
              e.stopPropagation();
              hideTip();
              _noveltyOnlyOpen(row["Specimen ID"]);
            });
          itr.appendChild(itd);
          return itr;
        }
        // Empty-sample indicator row: no detections, no VF/AMR, no novelty evidence.
        if (row.__emptyOnly) {
          const itr = document.createElement("tr");
          itr.className = "empty-only-row";
          const itd = document.createElement("td");
          itd.colSpan = visibleCols.length || 1;
          itd.innerHTML = _emptyOnlyMessageHTML(row["Specimen ID"], row.__emptyMeta);
          itr.appendChild(itd);
          return itr;
        }
        // VF/AMR-only indicator row: sample has no passing detections but carries
        // virulence / resistance gene hits. Render a single spanning cell.
        if (row.__vfamrOnly) {
          const itr = document.createElement("tr");
          itr.className = "vfamr-only-row";
          const itd = document.createElement("td");
          itd.colSpan = visibleCols.length || 1;
          itd.innerHTML = _vfamrOnlyMessageHTML(row["Specimen ID"], row.__vfamr);
          const _s = row.__vfamr;
          const _msg = itd.querySelector(".vfamr-only-msg");
          if (_msg) {
            _msg.addEventListener("mouseover", (ev) => showTip(_vfamrOnlyTip(row["Specimen ID"], _s), ev));
            _msg.addEventListener("mousemove", moveTip);
            _msg.addEventListener("mouseout", hideTip);
          }
          const _lnk = itd.querySelector(".vfamr-only-link");
          if (_lnk)
            _lnk.addEventListener("click", (e) => {
              e.stopPropagation();
              hideTip();
              _vfamrOnlyOpen(_s);
            });
          itr.appendChild(itd);
          return itr;
        }
        const tr = document.createElement("tr");
        if (isTruthy(row["High Consequence"])) tr.classList.add("hc-row");
        if (row.__belowCutoffVFAMR) tr.classList.add("below-cutoff-row");
        if (row.__belowCutoffAligned) tr.classList.add("below-cutoff-aligned-row");
        if (row.__noveltyNoAlign) tr.classList.add("novelty-noalign-row");
        const _rk = _tblRowKey(row);
        if (_tblPinned.has(_rk)) tr.classList.add("row-pinned");
        tr.addEventListener("click", (e) => {
          if (e.target.closest("a") || e.target.closest(".watch-star")) return;
          if (_tblPinned.has(_rk)) {
            _tblPinned.delete(_rk);
            tr.classList.remove("row-pinned");
          } else {
            _tblPinned.add(_rk);
            tr.classList.add("row-pinned");
          }
          _updateComparePinnedBtn();
        });
        let _orgTd = null;
        visibleCols.forEach((c) => {
          const td = document.createElement("td");
          const val = row[c] != null ? row[c] : "";
          if (c === "Detected Organism") {
            _orgTd = td;
            const isHC = row["High Consequence"];
            const mt = (row["Mol Type"] || "").toLowerCase();
            const rowLevel = row["Level"] || "Strain";
            td.style.cssText =
              "position:relative;white-space:nowrap;overflow:hidden;padding-right:70px;padding-left:4px;" +
              (isHC ? "color:#c62828;font-weight:600;" : "");
            // Star icon for watchlist toggle — sits inline right after the organism
            // text (appended below). Click star to watch; click row to pin.
            const _wk = _watchKey(row);
            const _starred = watchlist.has(_wk);
            const starEl = document.createElement("i");
            starEl.className = "watch-star" + (_starred ? " on" : "");
            starEl.title = _starred ? "On the follow-up list — click to remove" : "Add to follow-up list";
            starEl.style.cssText = "margin-left:5px;font-size:11px;cursor:pointer;vertical-align:middle;";
            starEl.innerHTML = `<i class="${_starred ? "fas" : "far"} fa-star"></i>`;
            starEl.addEventListener("click", (e) => {
              e.stopPropagation();
              _toggleWatchKey(_wk);
            });
            const contentSpan = document.createElement("span");
            contentSpan.style.cssText =
              "display:inline-block;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;vertical-align:middle;max-width:100%";
            if (isHC) {
              const badge = document.createElement("span");
              badge.textContent = "⚠ ";
              badge.title = "High Consequence Pathogen";
              badge.style.cssText = "color:#c62828;margin-right:2px";
              contentSpan.appendChild(badge);
            }
            contentSpan.appendChild(document.createTextNode(val));
            if (mt === "dna" || mt === "rna") {
              const mtBadge = document.createElement("span");
              mtBadge.textContent = mt === "dna" ? "D" : "R";
              mtBadge.title = mt === "dna" ? "DNA pathogen" : "RNA pathogen";
              mtBadge.style.cssText = `display:inline-flex;align-items:center;justify-content:center;width:14px;height:14px;border-radius:50%;background:${
                mt === "dna" ? "#1565c0" : "#6a1b9a"
              };color:#fff;font-size:8px;font-weight:700;vertical-align:middle;margin-left:3px;line-height:1;flex-shrink:0`;
              contentSpan.appendChild(mtBadge);
            }
            if (row.__belowCutoffVFAMR) {
              const _bc = document.createElement("span");
              _bc.className = "below-cutoff-badge";
              _bc.textContent = "↓ below cutoff";
              const _pi = rowPassInfo(row);
              _bc.title =
                `TASS ${isNaN(_pi.strain) ? "—" : _pi.strain.toFixed(1)} is below the cutoff (${_pi.thr}) — ` +
                `hidden by score, but shown because VF/AMR genes were detected for this organism's genus in this sample.`;
              contentSpan.appendChild(_bc);
            }
            if (row.__belowCutoffAligned || row.__noveltyNoAlign) {
              const _stb = document.createElement("span");
              _stb.innerHTML = _subThresholdBadgeHTML(row);
              while (_stb.firstChild) contentSpan.appendChild(_stb.firstChild);
            }
            td.appendChild(contentSpan);
            // Star sits inline, right after the organism text.
            td.appendChild(starEl);
            // Right-anchored group: level badge (flush right when not hovered), then
            // the pin-hint thumbtack which only appears on hover/pin and pushes the
            // badge left. Absolute so overflow:hidden on the cell doesn't clip it.
            const rightGroup = document.createElement("span");
            rightGroup.style.cssText =
              "position:absolute;right:6px;top:50%;transform:translateY(-50%);display:inline-flex;align-items:center;white-space:nowrap;";
            // Taxonomic level badge: always shown so users can tell whether
            // a row is a direct Strain call, a Species aggregate, or a Genus aggregate.
            const _lvlStyles = {
              Strain: { bg: "#e3f2fd", fg: "#1565c0", label: "Strain" },
              Species: { bg: "#f3e5f5", fg: "#6a1b9a", label: "Species" },
              Genus: { bg: "#e8f5e9", fg: "#2e7d32", label: "Genus" },
            };
            const _ls = _lvlStyles[rowLevel] || { bg: "#f5f5f5", fg: "#555", label: rowLevel };
            const lvlBadge = document.createElement("span");
            lvlBadge.textContent = _ls.label;
            lvlBadge.title =
              rowLevel === "Strain"
                ? "Strain-level detection — TASS reflects this specific strain's alignments."
                : rowLevel === "Species"
                ? "Species-level aggregate — TASS reflects all reads assigned to this species."
                : "Genus-level aggregate — TASS reflects all reads assigned to this genus.";
            lvlBadge.style.cssText =
              `display:inline-block;font-size:9px;font-weight:700;padding:0 4px;border-radius:3px;line-height:1.5;` +
              `background:${_ls.bg};color:${_ls.fg};border:1px solid ${_ls.fg}44`;
            rightGroup.appendChild(lvlBadge);
            // Pin-hint thumbtack — last, so it sits at the very end of the cell.
            const pinHint = document.createElement("span");
            pinHint.className = "row-pin-hint";
            rightGroup.appendChild(pinHint);
            td.appendChild(rightGroup);
          } else {
            td.textContent = val;
          }
          // User annotation column: make the cell editable. Commits write back to
          // TT_ANNOT.rowData (keyed by the stable row key) AND onto the row object
          // so sort / export stay consistent.
          if (TT_ANNOT.rowCols.includes(c)) {
            td.contentEditable = "true";
            td.spellcheck = false;
            td.dataset.annotCol = c;
            td.classList.add("tt-annot-cell");
            td.style.cssText = "min-width:90px;white-space:pre-wrap;cursor:text;background:#fffdf5;outline:none";
            td.title = "Click to edit — saved with exported state";
            const _commit = () => {
              const v = td.innerText.trim();
              const store = (TT_ANNOT.rowData[_rk] = TT_ANNOT.rowData[_rk] || {});
              store[c] = v;
              row[c] = v;
            };
            td.addEventListener("click", (e) => e.stopPropagation()); // don't toggle pin while editing
            td.addEventListener("blur", _commit);
            td.addEventListener("keydown", (e) => {
              if (e.key === "Enter" && !e.shiftKey) {
                e.preventDefault();
                td.blur();
              }
            });
          }
          // Taxonomic ID #: link out to the NCBI Taxonomy browser.
          if (c === "Taxonomic ID #" && /^\d+$/.test(String(val).trim())) {
            const tid = String(val).trim();
            td.textContent = "";
            const a = document.createElement("a");
            a.href = `https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=${tid}`;
            a.target = "_blank";
            a.rel = "noopener";
            a.textContent = tid;
            a.title = `View taxon ${tid} on NCBI`;
            a.style.color = "#1565c0";
            a.addEventListener("click", (e) => e.stopPropagation());
            td.appendChild(a);
          }
          if (c === "BSL Level" && val) {
            const n = parseInt(val.replace("BSL-", "")) || 0;
            const bg = { 1: "#e8f5e9", 2: "#fff9c4", 3: "#ffe0b2", 4: "#ffcdd2" }[n] || "";
            const fg = { 1: "#2e7d32", 2: "#b45309", 3: "#c2410c", 4: "#b91c1c" }[n] || "";
            td.style.cssText = `background:${bg};color:${fg};font-weight:700;text-align:center;white-space:nowrap`;
          }
          // TASS Score: max 2 decimal places + linear background 0–100, low-alpha blue → green gradient
          if (c === "TASS Score") {
            const score = parseFloat(val);
            if (!isNaN(score)) td.textContent = +score.toFixed(2);
            if (!isNaN(score) && score > 0) {
              const t = Math.min(1, Math.max(0, score / 100));
              // interpolate hue 210 (blue) → 145 (green) as score rises
              const hue = Math.round(210 - t * 65);
              const alpha = (0.08 + t * 0.27).toFixed(3); // 0.08 @ 0 → 0.35 @ 100
              td.style.background = `hsla(${hue},70%,45%,${alpha})`;
            }
            td.style.textAlign = "center";
            td.style.fontWeight = "600";
            // Species-rollup hit: mark the row with the rescued-row style.
            // The level badge in the Detected Organism column already identifies
            // this as a Species row — no duplicate text badge needed here.
            if (row.__rescuedSpecies) {
              tr.classList.add("rescued-row");
            }
            // Pass / rescue indicator. Passing rows get a green cell background;
            // no text badge. Rescued rows get the ↑ rollup badge. ✗ below is
            // shown in Species/Genus view when the aggregation itself doesn't pass.
            const _vlBadge = (document.getElementById("view-level") || {}).value || "Strain";
            const _pi = rowPassInfo(row);
            if (!row.__rescuedSpecies && !isNaN(_pi.strain)) {
              if (_vlBadge !== "Strain") {
                // Species / Genus view: green for pass, rescue badge if rescued, ✗ below otherwise.
                if (_pi.strainPass) {
                  td.style.background = "#c8e6c9";
                  td.title = _vlBadge + " passes the TASS threshold (≥ " + _pi.thr + ")";
                } else if (_pi.rescued) {
                  const badge = document.createElement("span");
                  badge.style.cssText =
                    "display:inline-block;margin-left:5px;font-size:9px;font-weight:700;" +
                    "padding:0 4px;border-radius:3px;vertical-align:middle;line-height:1.5;" +
                    "background:#ffe0b2;color:#c2410c";
                  badge.textContent = _pi.rescueLevel === "species" ? "↑ species" : "↑ genus";
                  badge.title =
                    "Strain is BELOW the cutoff (" +
                    _pi.strain.toFixed(1) +
                    " < " +
                    _pi.thr +
                    ") but kept visible because its " +
                    _pi.rescueLevel +
                    " aggregation passes (" +
                    (_pi.rescueLevel === "species"
                      ? "Species TASS " + (isNaN(_pi.species) ? "?" : _pi.species.toFixed(1))
                      : "Genus TASS " + (isNaN(_pi.genus) ? "?" : _pi.genus.toFixed(1))) +
                    " ≥ " +
                    _pi.thr +
                    '). Turn off "Roll up threshold" to hide it.';
                  tr.classList.add("rescued-row");
                  td.appendChild(badge);
                } else {
                  const badge = document.createElement("span");
                  badge.style.cssText =
                    "display:inline-block;margin-left:5px;font-size:9px;font-weight:700;" +
                    "padding:0 4px;border-radius:3px;vertical-align:middle;line-height:1.5;" +
                    "background:#ffcdd2;color:#c62828";
                  badge.textContent = "✗ below";
                  badge.title = _vlBadge + " is below the TASS cutoff (" + _pi.thr + ")";
                  td.appendChild(badge);
                }
              } else if (_pi.strainPass) {
                // Strain view, passes on its own: green cell, no badge.
                td.style.background = "#c8e6c9";
                td.title = "Strain passes the TASS threshold (≥ " + _pi.thr + ")";
              } else if (_pi.rescued) {
                const badge = document.createElement("span");
                badge.style.cssText =
                  "display:inline-block;margin-left:5px;font-size:9px;font-weight:700;" +
                  "padding:0 4px;border-radius:3px;vertical-align:middle;line-height:1.5;" +
                  "background:#ffe0b2;color:#c2410c";
                badge.textContent = _pi.rescueLevel === "species" ? "↑ species" : "↑ genus";
                badge.title =
                  "Strain is BELOW the cutoff (" +
                  _pi.strain.toFixed(1) +
                  " < " +
                  _pi.thr +
                  ") but kept visible because its " +
                  _pi.rescueLevel +
                  " aggregation passes (" +
                  (_pi.rescueLevel === "species"
                    ? "Species TASS " + (isNaN(_pi.species) ? "?" : _pi.species.toFixed(1))
                    : "Genus TASS " + (isNaN(_pi.genus) ? "?" : _pi.genus.toFixed(1))) +
                  " ≥ " +
                  _pi.thr +
                  '). Turn off "Roll up threshold" to hide it.';
                tr.classList.add("rescued-row");
                td.appendChild(badge);
              }
            }
          }
          // Species TASS / Genus TASS: same gradient styling as TASS Score,
          // 2-dp, so parent confidence is legible alongside the strain score.
          if (c === "Species TASS" || c === "Genus TASS") {
            const pscore = parseFloat(val);
            if (!isNaN(pscore)) {
              td.textContent = +pscore.toFixed(2);
              const t = Math.min(1, Math.max(0, pscore / 100));
              const hue = Math.round(210 - t * 65);
              const alpha = (0.05 + t * 0.18).toFixed(3);
              td.style.background = `hsla(${hue},55%,45%,${alpha})`;
              td.style.textAlign = "center";
            }
          }
          tr.appendChild(td);
        });
        // Strain-breakdown hover tooltip — only on the Detected Organism cell
        const _trLevel = row["Level"] || "Strain";
        if ((_trLevel === "Species" || _trLevel === "Genus") && _orgTd) {
          _orgTd.style.cursor = "help";
          _orgTd.addEventListener("mouseenter", function (e) {
            const _tipStrains = _getTopStrainsForRow(row);
            if (_tipStrains.length) {
              const _tipHtml = _strainTooltipHtml(_tipStrains, _trLevel, row["Detected Organism"]);
              if (_tipHtml) showTip(_tipHtml, e);
            }
          });
          _orgTd.addEventListener("mousemove", moveTip);
          _orgTd.addEventListener("mouseleave", hideTip);
        }
        return tr;
      }

      function populateTable() {
        let fd = filteredData();
        // Append faded below-cutoff rows that carry VF/AMR hits (toggle-gated).
        const _belowExtra = _belowCutoffExtraRows();
        if (_belowExtra.length) fd = fd.concat(_belowExtra);
        // Append faded sub-threshold (aligned) + no-alignment novelty rows (toggle-gated).
        const _novSubExtra = _noveltySubThresholdExtraRows();
        if (_novSubExtra.length) fd = fd.concat(_novSubExtra);
        // Append indicator rows for samples with VF/AMR hits but no detection rows.
        const _vfOnly = _vfamrOnlySampleRows(fd);
        if (_vfOnly.length) fd = fd.concat(_vfOnly);
        // Append novelty-only indicator rows for samples lacking visible detections.
        const _novOnly = _noveltyOnlySampleRows(fd);
        if (_novOnly.length) fd = fd.concat(_novOnly);
        // Append empty-sample indicator rows for samples with no detections at all.
        const _emptyOnly = _emptyOnlySampleRows(fd);
        if (_emptyOnly.length) fd = fd.concat(_emptyOnly);
        const _tblGrouped = document.getElementById("tbl-group-sample")?.checked;
        if (_tblGrouped) {
          // Group by Sample: primary sort follows _sampleOrder (right-panel order),
          // secondary sort on the active column (or TASS desc by default) within each group.
          const sc = sortCol || "TASS Score";
          const _grpIdx = _sampleOrder.length ? Object.fromEntries(_sampleOrder.map((id, i) => [id, i])) : {};
          fd = [...fd].sort((a, b) => {
            const sia = String(a["Specimen ID"] || "");
            const sib = String(b["Specimen ID"] || "");
            const ia = _grpIdx[sia] !== undefined ? _grpIdx[sia] : 9999;
            const ib = _grpIdx[sib] !== undefined ? _grpIdx[sib] : 9999;
            if (ia !== ib) return ia - ib;
            // Fallback within same group: active sort col
            const na = parseFloat(a[sc]),
              nb = parseFloat(b[sc]);
            if (!isNaN(na) && !isNaN(nb)) return sortCol ? (sortAsc ? na - nb : nb - na) : nb - na;
            return sortAsc
              ? String(a[sc] || "").localeCompare(String(b[sc] || ""))
              : String(b[sc] || "").localeCompare(String(a[sc] || ""));
          });
        } else if (sortCol) {
          fd = [...fd].sort((a, b) => {
            const va = a[sortCol],
              vb = b[sortCol];
            const na = parseFloat(va),
              nb = parseFloat(vb);
            if (!isNaN(na) && !isNaN(nb)) return sortAsc ? na - nb : nb - na;
            return sortAsc
              ? String(va || "").localeCompare(String(vb || ""))
              : String(vb || "").localeCompare(String(va || ""));
          });
        }

        const total = fd.length;
        const pageSize = _tblPageSize();
        const pages = pageSize > 0 ? Math.max(1, Math.ceil(total / pageSize)) : 1;

        // Clamp page index in case total shrank
        _tblPage = Math.max(0, Math.min(_tblPage, pages - 1));

        const start = pageSize > 0 ? _tblPage * pageSize : 0;
        const end = pageSize > 0 ? Math.min(start + pageSize, total) : total;
        const slice = pageSize > 0 ? fd.slice(start, end) : fd;

        // ── render rows ──────────────────────────────────────────────────────
        const tbody = document.getElementById("table-body");
        tbody.innerHTML = "";
        const frag = document.createDocumentFragment();
        let _lastGrp = null;
        slice.forEach((row) => {
          if (_tblGrouped && row["Specimen ID"] !== _lastGrp) {
            _lastGrp = row["Specimen ID"];
            const gr = document.createElement("tr");
            gr.className = "grp-row";
            const gtd = document.createElement("td");
            gtd.colSpan = visibleCols.length;
            gtd.innerHTML = `<i class="fas fa-vial"></i> ${_lastGrp}`;
            gr.appendChild(gtd);
            frag.appendChild(gr);
          }
          frag.appendChild(_renderOneRow(row));
        });
        tbody.appendChild(frag);

        // ── update toolbar ───────────────────────────────────────────────────
        const countEl = document.getElementById("tbl-row-count");
        if (countEl) {
          if (total === 0) {
            countEl.textContent = "No rows";
          } else if (pageSize > 0 && pages > 1) {
            countEl.textContent = `Rows ${start + 1}–${end} of ${total.toLocaleString()}`;
          } else {
            countEl.textContent = `${total.toLocaleString()} row${total !== 1 ? "s" : ""}`;
          }
        }

        const labelEl = document.getElementById("tbl-page-label");
        const firstBtn = document.getElementById("tbl-first");
        const prevBtn = document.getElementById("tbl-prev");
        const nextBtn = document.getElementById("tbl-next");
        const lastBtn = document.getElementById("tbl-last");
        const pager = document.getElementById("tbl-pager");

        const showPager = pageSize > 0 && pages > 1;
        if (pager) pager.style.display = showPager ? "flex" : "none";

        if (showPager) {
          if (labelEl) labelEl.textContent = `Page ${_tblPage + 1} of ${pages}`;
          if (firstBtn) firstBtn.disabled = _tblPage === 0;
          if (prevBtn) prevBtn.disabled = _tblPage === 0;
          if (nextBtn) nextBtn.disabled = _tblPage >= pages - 1;
          if (lastBtn) lastBtn.disabled = _tblPage >= pages - 1;
        }
      }

      /* ── Wire pagination controls (runs once after DOM is ready) ────────────── */
      (function _initTablePagination() {
        function _goPage(n) {
          _tblPage = n;
          populateTable();
        }

        document.getElementById("tbl-first").addEventListener("click", () => _goPage(0));
        document.getElementById("tbl-prev").addEventListener("click", () => _goPage(_tblPage - 1));
        document.getElementById("tbl-next").addEventListener("click", () => _goPage(_tblPage + 1));
        document.getElementById("tbl-last").addEventListener("click", () => {
          const total = filteredData().length;
          const pageSize = _tblPageSize();
          const pages = pageSize > 0 ? Math.max(1, Math.ceil(total / pageSize)) : 1;
          _goPage(pages - 1);
        });
        document.getElementById("tbl-page-size").addEventListener("change", () => {
          _tblPage = 0;
          populateTable();
        });
        const _tblGrpEl = document.getElementById("tbl-group-sample");
        if (_tblGrpEl)
          _tblGrpEl.addEventListener("change", () => {
            _tblPage = 0;
            buildTable();
            renderTableHeaders();
            populateTable();
          });

        // ── Add an editable annotation column to the detection table ──────────
        const _addNoteBtn = document.getElementById("tbl-add-note-col");
        if (_addNoteBtn)
          _addNoteBtn.addEventListener("click", () => {
            const name = (window.prompt("New column name (e.g. Notes, Reviewer, Action):", "Notes") || "").trim();
            if (!name) return;
            if (ALL_COLS.includes(name) || TT_ANNOT.rowCols.includes(name)) {
              alert("A column named “" + name + "” already exists.");
              return;
            }
            TT_ANNOT.rowCols.push(name);
            buildTable();
            renderTableHeaders();
            populateTable();
          });
      })();

