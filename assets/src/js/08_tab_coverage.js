/* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB: COVERAGE           (data-tab="coverage")
       -     drawCoverage()  →  zoomable bubble scatter; X = breadth, Y = depth
       -     by default, bubble size = TASS Score. Each axis (X, Y, R) and the
       -     scale type (linear / log / sqrt) is user-selectable. Colored per
       -     sample; the colour legend sits BELOW the plot as a collapsible,
       -     one-row-per-sample list that isolates a sample's dots on hover.
═══════════════════════════════════════════════════════════════════════════ */
// Whether the coverage colour legend is expanded. null = not chosen yet, so the
// first render picks a default from the sample count; after that the user's
// choice survives redraws.
let _covLegendOpen = null;
// Legend-local search text. Scoped to the legend on purpose: it hides legend
// entries and their dots only, and never touches filteredData(), so the rest of
// the report (tables, other tabs, exports) is unaffected by typing here.
let _covLegendQuery = "";
// Sample whose dots are locked in isolation by a click on its legend row.
// Hover is transient; this survives redraws so a filter tweak doesn't drop the
// sample you were tracking. Cleared when that sample leaves the filtered set.
let _covLegendPin = null;

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

  if (!fd.length) {
    wrap.innerHTML = '<p style="color:#999;padding:1em">No detection.</p>';
    _renderCovLegend(null, fd);
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
  svg.append("defs").append("clipPath").attr("id", "cov-clip").append("rect").attr("width", iW).attr("height", iH);

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

  _renderCovLegend(circles, fd);
}

/* ── Coverage colour legend ────────────────────────────────────────────────
   Lives BELOW the plot (#coverage-legend-wrap follows #coverage-svg-wrap in
   the markup): a colour key is a lookup you consult after reading the chart,
   so it should never be what you scroll past to reach it — on a 164-specimen
   run the old wrapped chip row pushed the plot itself off screen.

   Entries run left-to-right and wrap, numbered so the reading order is still
   explicit, in the same sample order used everywhere else in the report. A
   search box filters the legend itself — matching entries and their dots stay,
   the rest drop out of the plot — which is deliberately NOT a report filter:
   nothing outside this legend and this scatter changes.

   Hovering (or keyboard-focusing) an entry isolates that sample: its dots
   come forward at full strength and every other dot fades to a faint ghost, so
   you can see where one sample sits without losing the shape of the whole
   cloud. Click locks that isolation so you can move the mouse to the plot. */
function _renderCovLegend(circlesSel, fd) {
  const legWrap = document.getElementById("coverage-legend-wrap");
  if (!legWrap) return;
  const legSamples = _orderedSamples(uniq(fd.map((r) => r["Specimen ID"])).filter(Boolean));
  const unit = typeof specimenMergeEnabled !== "undefined" && specimenMergeEnabled ? "specimen" : "sample";
  const n = legSamples.length;
  const AUTO_OPEN_MAX = 8;
  // Remember the user's choice across redraws (filters redraw constantly);
  // otherwise the legend would snap shut on every keystroke.
  if (_covLegendOpen === null) _covLegendOpen = n > 0 && n <= AUTO_OPEN_MAX;
  if (_covLegendPin && legSamples.indexOf(_covLegendPin) === -1) _covLegendPin = null;

  if (!n) {
    legWrap.innerHTML = "";
    return;
  }

  const nPts = {};
  fd.forEach((r) => {
    const sp = r["Specimen ID"];
    if (sp) nPts[sp] = (nPts[sp] || 0) + 1;
  });

  const rows = legSamples
    .map((sp, i) => {
      const col = sampleColors[sp] || "#90a4ae";
      const cnt = nPts[sp] || 0;
      return (
        `<li><button type="button" class="cov-leg-item${_covLegendPin === sp ? " pinned" : ""}"` +
        ` data-sp="${_esc(sp)}" title="Hover to isolate this ${unit}'s points · click to lock">` +
        `<span class="cov-leg-idx">${i + 1}</span>` +
        `<span class="cov-leg-dot" style="background:${_esc(col)}"></span>` +
        `<span class="cov-leg-name">${_esc(sp)}</span>` +
        `<span class="cov-leg-n">${cnt.toLocaleString()}</span>` +
        `</button></li>`
      );
    })
    .join("");

  legWrap.innerHTML =
    `<details id="cov-legend-details"${_covLegendOpen ? " open" : ""}>` +
    `<summary><b>${n}</b> ${unit}${n === 1 ? "" : "s"} plotted` +
    `<span style="color:#8a97a4;font-weight:400"> — colour legend</span>` +
    `<span class="cov-leg-hint">hover an entry to isolate its points</span></summary>` +
    `<div class="cov-leg-tools">` +
    `<input type="search" id="cov-leg-search" placeholder="Search ${unit}s in this legend…"` +
    ` value="${_esc(_covLegendQuery)}" autocomplete="off"` +
    ` title="Filters this legend and its dots only — the rest of the report is untouched" />` +
    `<span id="cov-leg-count"></span>` +
    `</div>` +
    `<ol id="cov-legend-list">${rows}</ol>` +
    `</details>`;

  const det = document.getElementById("cov-legend-details");
  if (det) det.addEventListener("toggle", () => (_covLegendOpen = det.open));

  // ── Isolate-on-hover ────────────────────────────────────────────────
  // Attribute writes on the existing selection, not a redraw: the zoom
  // transform and every dot position stay exactly as they are.
  const list = document.getElementById("cov-legend-list");
  if (!list || !circlesSel) return;

  const paint = (sp) => {
    if (!sp) {
      circlesSel.attr("opacity", 0.75).attr("stroke", "#fff").attr("stroke-width", 0.8);
    } else {
      circlesSel
        .attr("opacity", (r) => (r["Specimen ID"] === sp ? 0.95 : 0.07))
        .attr("stroke", (r) => (r["Specimen ID"] === sp ? "#263238" : "#fff"))
        .attr("stroke-width", (r) => (r["Specimen ID"] === sp ? 1.4 : 0.8));
      // Bring the isolated dots above the ghosts they were drawn under.
      circlesSel.filter((r) => r["Specimen ID"] === sp).raise();
    }
    list.querySelectorAll(".cov-leg-item").forEach((b) => {
      b.classList.toggle("dim", !!sp && b.dataset.sp !== sp);
      b.classList.toggle("active", !!sp && b.dataset.sp === sp);
    });
  };

  const enter = (ev) => {
    const btn = ev.target.closest(".cov-leg-item");
    if (btn) paint(btn.dataset.sp);
  };
  const leave = () => paint(_covLegendPin);

  list.addEventListener("mouseover", enter);
  list.addEventListener("focusin", enter);
  list.addEventListener("mouseleave", leave);
  list.addEventListener("focusout", (ev) => {
    if (!list.contains(ev.relatedTarget)) leave();
  });
  list.addEventListener("click", (ev) => {
    const btn = ev.target.closest(".cov-leg-item");
    if (!btn) return;
    _covLegendPin = _covLegendPin === btn.dataset.sp ? null : btn.dataset.sp;
    list.querySelectorAll(".cov-leg-item").forEach((b) => b.classList.toggle("pinned", b.dataset.sp === _covLegendPin));
    paint(_covLegendPin || btn.dataset.sp);
  });

  // ── Legend-scoped search ────────────────────────────────────────────
  // Hides non-matching entries and their dots. Nothing here calls a redraw or
  // touches the global filter state — clearing the box restores everything.
  const searchEl = document.getElementById("cov-leg-search");
  const countEl = document.getElementById("cov-leg-count");
  const _matches = (sp) => {
    const q = (_covLegendQuery || "").trim().toLowerCase();
    return !q || String(sp).toLowerCase().indexOf(q) !== -1;
  };
  const applySearch = () => {
    let shown = 0;
    list.querySelectorAll(".cov-leg-item").forEach((b) => {
      const ok = _matches(b.dataset.sp);
      b.parentElement.style.display = ok ? "" : "none";
      if (ok) shown++;
    });
    circlesSel.style("display", (r) => (_matches(r["Specimen ID"]) ? null : "none"));
    if (countEl) {
      const q = (_covLegendQuery || "").trim();
      countEl.textContent = q ? `${shown} of ${n} shown` : "";
      countEl.classList.toggle("none", !!q && shown === 0);
    }
    // An isolation locked on a sample the search just hid would leave every
    // remaining dot ghosted with nothing lit — drop it instead.
    if (_covLegendPin && !_matches(_covLegendPin)) {
      _covLegendPin = null;
      list.querySelectorAll(".cov-leg-item.pinned").forEach((b) => b.classList.remove("pinned"));
      paint(null);
    }
  };
  if (searchEl) {
    searchEl.addEventListener("input", () => {
      _covLegendQuery = searchEl.value;
      applySearch();
    });
    // Escape clears the box without bubbling out to any global key handler.
    searchEl.addEventListener("keydown", (ev) => {
      if (ev.key === "Escape") {
        ev.stopPropagation();
        searchEl.value = "";
        _covLegendQuery = "";
        applySearch();
      }
    });
  }
  applySearch();

  // A pin set before this redraw is re-applied to the fresh circles.
  if (_covLegendPin) paint(_covLegendPin);
}
