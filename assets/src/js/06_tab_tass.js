      /* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB: TASS COMPARISON    (data-tab="tass")
       -     drawTassChart()  →  D3 grouped bar chart. One group per taxonomy
       -     rank value, one bar per sample, height = max TASS Score.
       -     UI controls: rank selector, Y-scale (linear / log / sqrt).
       -     Sample legend is intentionally NOT drawn — colors are shown in
       -     the right sidebar.
═══════════════════════════════════════════════════════════════════════════ */
      function drawTassChart() {
        const wrap = document.getElementById("tass-svg-wrap");
        wrap.innerHTML = "";
        const fd = filteredData();
        const rank = document.getElementById("tass-rank").value || "Genus";
        const rankField = rank === "Species" ? "Detected Organism" : rank;

        const samples = _orderedSamples(uniq(fd.map((r) => r["Specimen ID"])).filter(Boolean));
        const groupKeys = uniq(fd.map((r) => r[rankField] || "Unknown"))
          .filter(Boolean)
          .sort();

        if (!groupKeys.length) {
          wrap.innerHTML =
            '<p style="color:#999;padding:1em">No detection (taxonomy ranks may not be populated; ensure taxdump is provided).</p>';
          return;
        }

        // TASS level: plot strain score, or the rolled-up species/genus TASS.
        const _tassLevel = (document.getElementById("tass-level") || {}).value || "strain";
        const _tassCol =
          _tassLevel === "species" ? "Species TASS" : _tassLevel === "genus" ? "Genus TASS" : "TASS Score";

        // Build grouped data: {group: {sample: maxTASS}}
        const grouped = {};
        fd.forEach((r) => {
          const g = r[rankField] || "Unknown";
          const sp = r["Specimen ID"] || "";
          // Fall back to strain TASS when the chosen parent column is absent.
          let v = num(r[_tassCol]);
          if (_tassCol !== "TASS Score" && (isNaN(v) || v === 0) && r[_tassCol] == null) v = num(r["TASS Score"]);
          if (!grouped[g]) grouped[g] = {};
          grouped[g][sp] = Math.max(grouped[g][sp] || 0, v);
        });

        const marginL = 130,
          marginT = 16,
          marginR = 20;
        // Dynamic bottom margin: enough room for rotated (-35°) x-axis labels
        const _maxKeyLen = groupKeys.reduce((m, k) => Math.max(m, k.length), 0);
        const marginB = Math.max(80, Math.min(260, Math.ceil(_maxKeyLen * 6.2 * Math.sin((35 * Math.PI) / 180)) + 30));
        const W = Math.max(600, wrap.clientWidth || 900);
        const H = 420 + Math.max(0, marginB - 100); // grow SVG height with labels
        const innerW = W - marginL - marginR;
        const innerH = H - marginT - marginB;

        const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H);

        // Clip bars to plot area
        svg
          .append("defs")
          .append("clipPath")
          .attr("id", "tass-clip")
          .append("rect")
          .attr("width", innerW)
          .attr("height", innerH + 4);

        const g = svg.append("g").attr("transform", `translate(${marginL},${marginT})`);

        const x0 = d3.scaleBand().domain(groupKeys).range([0, innerW]).paddingInner(0.15).paddingOuter(0.1);
        const x1 = d3.scaleBand().domain(samples).range([0, x0.bandwidth()]).padding(0.05);
        // Y-axis scale type — linear default; log clamps to a positive minimum
        // since TASS scores can be 0 and log(0) is undefined.
        const tassScaleType = (document.getElementById("tass-scale") || {}).value || "linear";
        let y;
        if (tassScaleType === "log") {
          y = d3.scaleLog().domain([0.1, 100]).range([innerH, 0]).clamp(true);
        } else if (tassScaleType === "sqrt") {
          y = d3.scaleSqrt().domain([0, 100]).range([innerH, 0]).nice();
        } else {
          y = d3.scaleLinear().domain([0, 100]).range([innerH, 0]).nice();
        }

        // Gridlines (static)
        g.append("g")
          .attr("class", "grid")
          .selectAll("line")
          .data(y.ticks(5))
          .enter()
          .append("line")
          .attr("x1", 0)
          .attr("x2", innerW)
          .attr("y1", (d) => y(d))
          .attr("y2", (d) => y(d))
          .attr("stroke", "#e0e0e0")
          .attr("stroke-dasharray", "3,3");

        // Y axis (static) — use a SI-friendly tick format on log so we get
        // 0.1, 1, 10, 100 rather than scientific notation.
        const yAxisCall =
          tassScaleType === "log"
            ? d3
                .axisLeft(y)
                .ticks(5, "~s")
                .tickFormat((d) => d + "%")
            : d3
                .axisLeft(y)
                .ticks(5)
                .tickFormat((d) => d + "%");
        g.append("g").attr("class", "axis").call(yAxisCall);

        // X axis group — updated on zoom
        const xAxisG = g.append("g").attr("class", "axis").attr("transform", `translate(0,${innerH})`);

        // Zoom-capture rect (behind bars)
        const zoomRect = g
          .append("rect")
          .attr("width", innerW)
          .attr("height", innerH)
          .attr("fill", "none")
          .attr("pointer-events", "all")
          .attr("cursor", "grab");

        // Bars group (clipped)
        const barsG = g.append("g").attr("clip-path", "url(#tass-clip)");

        function _drawBars(x0c, x1c) {
          barsG.selectAll("g.tass-grp").remove();
          groupKeys.forEach((gk) => {
            const gData = grouped[gk] || {};
            const gg = barsG
              .append("g")
              .attr("class", "tass-grp")
              .attr("transform", `translate(${x0c(gk)},0)`);
            samples.forEach((sp) => {
              const v = gData[sp] || 0;
              const fill = sampleColors[sp] || "#90a4ae";
              gg.append("rect")
                .attr("x", x1c(sp))
                .attr("width", Math.max(0, x1c.bandwidth()))
                .attr("y", y(v))
                .attr("height", innerH - y(v))
                .attr("fill", fill)
                .attr("rx", 2)
                .on("mouseover", (ev) => showTip(`<b>${gk}</b> — ${sp}<br>TASS: <b>${v.toFixed(1)}</b>`, ev))
                .on("mousemove", moveTip)
                .on("mouseout", hideTip);
            });
          });
        }

        function _updateXAxis(x0c) {
          xAxisG
            .call(d3.axisBottom(x0c).tickSizeOuter(0))
            .selectAll("text")
            .attr("transform", "rotate(-35)")
            .style("text-anchor", "end");
        }

        _updateXAxis(x0);
        _drawBars(x0, x1);

        // ── Threshold overlay ────────────────────────────────────────────
        // Horizontal dashed line at the active TASS cutoff so it is obvious
        // which rank groups clear the threshold (a species/genus bar above the
        // line "rescues" its strains under the rollup). Uses the global min;
        // per-sample-type cutoffs vary, so this is the baseline reference.
        // Toggled by #tass-show-cutoff checkbox.
        (function () {
          const _showCutoff = (document.getElementById("tass-show-cutoff") || { checked: true }).checked;
          if (!_showCutoff) return;
          const _thr = parseFloat(document.getElementById("filter-min").value) || 0;
          if (_thr <= 0 || _thr > 100) return;
          const yT = y(_thr);
          const tg = g.append("g").attr("class", "tass-threshold");
          // White halo behind the dashed line so it stays legible over any bar color
          tg.append("line")
            .attr("x1", 0)
            .attr("x2", innerW)
            .attr("y1", yT)
            .attr("y2", yT)
            .attr("stroke", "white")
            .attr("stroke-width", 4)
            .attr("pointer-events", "none");
          // Dashed cutoff line on top
          tg.append("line")
            .attr("x1", 0)
            .attr("x2", innerW)
            .attr("y1", yT)
            .attr("y2", yT)
            .attr("stroke", "#e53935")
            .attr("stroke-width", 2)
            .attr("stroke-dasharray", "8,4")
            .attr("pointer-events", "none");
          // Pill background behind label so it's readable over any bar
          const _labelText = `TASS cutoff ${_thr}`;
          const _lblPadX = 5,
            _lblPadY = 3,
            _lblFontSize = 10;
          const _lblW = _labelText.length * 6 + _lblPadX * 2;
          const _lblH = _lblFontSize + _lblPadY * 2;
          tg.append("rect")
            .attr("x", innerW - 4 - _lblW)
            .attr("y", yT - _lblH - 3)
            .attr("width", _lblW)
            .attr("height", _lblH)
            .attr("fill", "#e53935")
            .attr("rx", 3)
            .attr("opacity", 0.9)
            .attr("pointer-events", "none");
          tg.append("text")
            .attr("x", innerW - 4 - _lblW / 2)
            .attr("y", yT - 3 - _lblPadY)
            .attr("text-anchor", "middle")
            .attr("dominant-baseline", "middle")
            .attr("font-size", _lblFontSize)
            .attr("font-weight", "700")
            .attr("fill", "white")
            .attr("pointer-events", "none")
            .text(_labelText);
        })();

        // D3 zoom — rescale band x0 range
        const zoom = d3
          .zoom()
          .scaleExtent([1, 10])
          .translateExtent([
            [0, 0],
            [innerW, innerH],
          ])
          .on("zoom", (event) => {
            const t = event.transform;
            x0.range([t.applyX(0), t.applyX(innerW)]);
            x1.range([0, x0.bandwidth()]);
            _updateXAxis(x0);
            _drawBars(x0, x1);
          });

        zoomRect.call(zoom);

        // Reset zoom button
        svg
          .append("text")
          .attr("x", marginL + innerW - 2)
          .attr("y", marginT - 8)
          .attr("text-anchor", "end")
          .attr("font-size", 10)
          .attr("fill", "#1565c0")
          .style("cursor", "pointer")
          .text("⟳ Reset zoom")
          .on("click", () => zoomRect.call(zoom.transform, d3.zoomIdentity));

        // (Legend removed — sample colors are already shown in the right panel.)

        // Y label
        svg
          .append("text")
          .attr("transform", `translate(12,${H / 2}) rotate(-90)`)
          .attr("text-anchor", "middle")
          .attr("font-size", 12)
          .text("Max TASS Score");
      }

