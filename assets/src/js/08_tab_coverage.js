      /* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB: COVERAGE           (data-tab="coverage")
       -     drawCoverage()  →  zoomable bubble scatter; X = breadth, Y = depth
       -     by default, bubble size = TASS Score. Each axis (X, Y, R) and the
       -     scale type (linear / log / sqrt) is user-selectable. Colored per
       -     sample, with legend rendered above the SVG.
═══════════════════════════════════════════════════════════════════════════ */
      function drawCoverage() {
        const wrap = document.getElementById("coverage-svg-wrap");
        wrap.innerHTML = "";
        const fd = filteredData();

        // Wire axis selectors once
        ["cov-x-sel", "cov-y-sel", "cov-r-sel"].forEach((sid) => {
          const el = document.getElementById(sid);
          if (el && !el._wired) {
            el._wired = true;
            el.addEventListener("change", drawCoverage);
          }
        });
        const xField = (document.getElementById("cov-x-sel") || {}).value || "Breadth %";
        const yField = (document.getElementById("cov-y-sel") || {}).value || "Mean Depth";
        const rField = (document.getElementById("cov-r-sel") || {}).value || "TASS Score";
        const covScaleType = (document.getElementById("cov-scale") || {}).value || "linear";

        // Legend — above SVG
        const legWrap = document.getElementById("coverage-legend-wrap");
        if (legWrap) {
          const legSamples = _orderedSamples(uniq(fd.map((r) => r["Specimen ID"])).filter(Boolean));
          legWrap.innerHTML = legSamples
            .map((sp) => {
              const col = sampleColors[sp] || "#90a4ae";
              return `<span style="display:flex;align-items:center;gap:.35em"><svg width="10" height="10" style="flex-shrink:0"><circle cx="5" cy="5" r="5" fill="${col}"/></svg><span>${sp}</span></span>`;
            })
            .join("");
        }

        if (!fd.length) {
          wrap.innerHTML = '<p style="color:#999;padding:1em">No detection.</p>';
          return;
        }

        const marginL = 65,
          marginT = 30,
          marginR = 40,
          marginB = 60;
        const W = Math.max(600, wrap.clientWidth || 900);
        const H = 420;
        const iW = W - marginL - marginR;
        const iH = H - marginT - marginB;

        const xMax = d3.max(fd, (r) => num(r[xField])) || 100;
        const yMax = d3.max(fd, (r) => num(r[yField])) || 10;
        const rMax = d3.max(fd, (r) => num(r[rField])) || 100;

        // Add 6% padding so dots at the max value are never clipped by the axis edge.
        const xPad = xMax * 0.06;
        const yPad = yMax * 0.06;

        // Build x/y scales according to the user-selected scale type. For
        // log we clamp to a positive minimum so that zero-valued points
        // collapse to the axis edge instead of erroring out.
        function _buildAxisScale(type, max, pad, range) {
          if (type === "log")
            return d3
              .scaleLog()
              .domain([0.01, Math.max(max + pad, 0.02)])
              .range(range)
              .clamp(true);
          if (type === "sqrt")
            return d3
              .scaleSqrt()
              .domain([0, max + pad])
              .range(range)
              .nice();
          return d3
            .scaleLinear()
            .domain([0, max + pad])
            .range(range)
            .nice();
        }
        const x = _buildAxisScale(covScaleType, xMax, xPad, [0, iW]);
        const y = _buildAxisScale(covScaleType, yMax, yPad, [iH, 0]);
        const rScale = d3.scaleSqrt().domain([0, rMax]).range([3, 14]);

        const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H);

        // Clip path to prevent dots drawing outside the plot area
        svg
          .append("defs")
          .append("clipPath")
          .attr("id", "cov-clip")
          .append("rect")
          .attr("width", iW)
          .attr("height", iH);

        const g = svg.append("g").attr("transform", `translate(${marginL},${marginT})`);

        // Axis groups (re-rendered on zoom)
        const xAxisG = g.append("g").attr("class", "axis").attr("transform", `translate(0,${iH})`);
        const yAxisG = g.append("g").attr("class", "axis");
        const xGrid = g.append("g").attr("class", "x-grid");
        const yGrid = g.append("g").attr("class", "y-grid");

        // Zoom-capture rect — kept BEHIND circles so hover events reach dots
        const zoomRect = g
          .append("rect")
          .attr("width", iW)
          .attr("height", iH)
          .attr("fill", "none")
          .attr("pointer-events", "all")
          .attr("cursor", "grab");

        // Dots group (clipped) — rendered on top of zoom rect so circles receive hover events
        const dotsG = g.append("g").attr("clip-path", "url(#cov-clip)");

        function renderAxes(xt, yt) {
          // Use compact SI tick format on log so values like 1k, 10k render
          // cleanly instead of scientific notation.
          const xAxis = covScaleType === "log" ? d3.axisBottom(xt).ticks(6, "~s") : d3.axisBottom(xt).ticks(8);
          const yAxis = covScaleType === "log" ? d3.axisLeft(yt).ticks(6, "~s") : d3.axisLeft(yt).ticks(8);
          xAxisG.call(xAxis);
          yAxisG.call(yAxis);
          xGrid
            .selectAll("line")
            .data(xt.ticks(covScaleType === "log" ? 6 : 8))
            .join("line")
            .attr("x1", (d) => xt(d))
            .attr("x2", (d) => xt(d))
            .attr("y1", 0)
            .attr("y2", iH)
            .attr("stroke", "#eee")
            .attr("stroke-dasharray", "3,3");
          yGrid
            .selectAll("line")
            .data(yt.ticks(covScaleType === "log" ? 6 : 8))
            .join("line")
            .attr("x1", 0)
            .attr("x2", iW)
            .attr("y1", (d) => yt(d))
            .attr("y2", (d) => yt(d))
            .attr("stroke", "#eee")
            .attr("stroke-dasharray", "3,3");
        }

        renderAxes(x, y);

        const circles = dotsG
          .selectAll("circle")
          .data(fd)
          .enter()
          .append("circle")
          .attr("cx", (r) => x(num(r[xField])))
          .attr("cy", (r) => y(num(r[yField])))
          .attr("r", (r) => rScale(num(r[rField])))
          .attr("fill", (r) => sampleColors[r["Specimen ID"]] || "#90a4ae")
          .attr("opacity", 0.75)
          .attr("stroke", "#fff")
          .attr("stroke-width", 0.8)
          .style("cursor", "pointer")
          .on("mouseover", (ev, r) =>
            showTip(
              `<b>${r["Detected Organism"]}</b><br>Sample: ${r["Specimen ID"]}<br>` +
                `${xField}: ${num(r[xField]).toFixed(2)} &nbsp; ${yField}: ${num(r[yField]).toFixed(2)}<br>` +
                `${rField}: ${num(r[rField]).toFixed(2)}<br>` +
                `Covered: ${r["Covered Bases"]} bp &nbsp; Length: ${r["Genome Length (bp)"]} bp<br>` +
                `HC: ${isTruthy(r["High Consequence"]) ? "Yes" : "No"}`,
              ev,
            ),
          )
          .on("mousemove", moveTip)
          .on("mouseout", hideTip);

        // ── D3 zoom ─────────────────────────────────────────────────────────────
        const zoom = d3
          .zoom()
          .scaleExtent([0.5, 50])
          .extent([
            [0, 0],
            [iW, iH],
          ])
          .translateExtent([
            [-iW, -iH],
            [2 * iW, 2 * iH],
          ])
          .on("zoom", (event) => {
            const xt = event.transform.rescaleX(x);
            const yt = event.transform.rescaleY(y);
            renderAxes(xt, yt);
            circles.attr("cx", (r) => xt(num(r[xField]))).attr("cy", (r) => yt(num(r[yField])));
          });

        // Attach zoom behaviour to the behind-circles rect
        zoomRect.call(zoom);

        // Reset zoom button
        svg
          .append("text")
          .attr("x", marginL + iW - 2)
          .attr("y", marginT - 8)
          .attr("text-anchor", "end")
          .attr("font-size", 10)
          .attr("fill", "#1565c0")
          .style("cursor", "pointer")
          .text("⟳ Reset zoom")
          .on("click", () => zoomRect.call(zoom.transform, d3.zoomIdentity));

        // Axes labels
        svg
          .append("text")
          .attr("transform", `translate(${marginL + iW / 2},${H - 14})`)
          .attr("text-anchor", "middle")
          .attr("font-size", 12)
          .text(xField);
        svg
          .append("text")
          .attr("transform", `translate(14,${marginT + iH / 2}) rotate(-90)`)
          .attr("text-anchor", "middle")
          .attr("font-size", 12)
          .text(yField);
        svg
          .append("text")
          .attr("x", marginL + iW / 2)
          .attr("y", marginT - 10)
          .attr("text-anchor", "middle")
          .attr("font-size", 9)
          .attr("fill", "#999")
          .text("Scroll to zoom · drag to pan");
      }

