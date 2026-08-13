/* ═══════════════════════════════════════════════════════════════════════════
       -  §  LONGITUDINAL ANALYSIS  —  SECOND DEFINITION  (this is the active one)
       -     Lives inside the Run Metadata tab. Renders a line plot:
       -        – X axis: collection_time
       -        – Y axis: user-selected metric
       -        – One line per selected organism (chip selector)
       -        – Per-plot sample visibility independent of the main sidebar
       -          sampleHidden map.
       -     Functions: _parseLongiDate, _buildLongitudinalSection,
       -     _populateLongiOrgs, _buildLongiSamplePanel, _drawLongitudinalPlot.
       -     State variables are declared above the init IIFE to avoid TDZ
       -     errors when build is called during page bootstrap.
═══════════════════════════════════════════════════════════════════════════ */

// ── Run-metadata table pagination state ──────────────────────────────
// The metadata grid holds one row per sample; for large runs we page the
// rows so the tab stays responsive. All rows stay in the DOM (so edit
// indices, highlight-scroll and full-table export keep working) —
// pagination only toggles per-row visibility.
let _runmetaPage = 1;
const _RUNMETA_PAGE_SIZE = 25;

// Standard metadata columns always shown in the Run Metadata table, even
// when the samplesheet supplied no value for any of them. This gives the
// user empty, editable cells to fill in and guarantees the geo /
// longitudinal / host sub-tabs have the fields they analyse.
const _STANDARD_META_COLS = [
  "sample_origin_country",
  "sample_origin_state_province_territory",
  "latitude",
  "longitude",
  "collection_time",
  "host_scientific_name",
  "host_disease",
  "environmental_site",
];

function _parseLongiDate(s) {
  if (!s) return null;
  const str = String(s).trim();
  // Try native Date parse first (handles ISO 8601 and many standard formats)
  let d = new Date(str);
  if (!isNaN(d)) return d;
  // M/D/YY h:mm  or  M/D/YYYY h:mm
  let m = str.match(/^(\d{1,2})\/(\d{1,2})\/(\d{2,4})\s+(\d{1,2}):(\d{2})/);
  if (m) {
    let yr = parseInt(m[3]);
    if (yr < 100) yr += 2000;
    d = new Date(yr, parseInt(m[1]) - 1, parseInt(m[2]), parseInt(m[4]), parseInt(m[5]));
    if (!isNaN(d)) return d;
  }
  // M/D/YY  or  M/D/YYYY  (date only)
  m = str.match(/^(\d{1,2})\/(\d{1,2})\/(\d{2,4})$/);
  if (m) {
    let yr = parseInt(m[3]);
    if (yr < 100) yr += 2000;
    d = new Date(yr, parseInt(m[1]) - 1, parseInt(m[2]));
    if (!isNaN(d)) return d;
  }
  return null;
}

function _buildLongitudinalSection() {
  const section = document.getElementById("longi-section");
  if (!section) return;

  // Build time map from live RUN_META
  const timeMap = {};
  (RUN_META || []).forEach((r) => {
    const d = _parseLongiDate(r.collection_time);
    if (d) timeMap[r.sample_name] = d;
  });

  const timedSamples = Object.keys(timeMap);
  if (typeof _updateMetaSubTabStates === "function") _updateMetaSubTabStates();
  if (!timedSamples.length) return;

  // Organisms that appear in at least one timed sample
  const eligible = new Set(timedSamples);
  const orgSet = new Set(
    filteredData()
      .filter((r) => eligible.has(r["Specimen ID"]))
      .map((r) => r["Detected Organism"] || "")
      .filter(Boolean),
  );
  _longiOrgList = [...orgSet].sort((a, b) => a.localeCompare(b));

  // Assign stable colours to every organism
  _longiOrgList.forEach((o, i) => {
    if (!_longiOrgColors[o]) _longiOrgColors[o] = _LONGI_ORG_PAL[i % _LONGI_ORG_PAL.length];
  });

  // Default: select the first organism if nothing is selected yet
  if (_longiSelectedOrgs.size === 0 && _longiOrgList.length > 0) {
    _longiSelectedOrgs.add(_longiOrgList[0]);
  }

  _populateLongiOrgs(_longiOrgList);
  _buildLongiSamplePanel(timedSamples);

  // Wire events once
  if (!_longiBuilt) {
    _longiBuilt = true;
    const orgSearch = document.getElementById("longi-org-search");
    const ySel = document.getElementById("longi-y-sel");
    const scaleSel = document.getElementById("longi-scale-sel");
    if (ySel) ySel.addEventListener("change", _drawLongitudinalPlot);
    if (scaleSel) scaleSel.addEventListener("change", _drawLongitudinalPlot);
    if (orgSearch)
      orgSearch.addEventListener("input", () => {
        const q = (orgSearch.value || "").trim().toLowerCase();
        _populateLongiOrgs(q ? _longiOrgList.filter((o) => o.toLowerCase().includes(q)) : _longiOrgList);
      });

    // Show All / Hide All (None) buttons in the sample panel header
    const showAllBtn = document.getElementById("longi-show-all");
    const hideAllBtn = document.getElementById("longi-hide-all");
    if (showAllBtn) {
      showAllBtn.addEventListener("click", () => {
        Object.keys(_longiHidden).forEach((id) => {
          _longiHidden[id] = false;
        });
        // Sync checkboxes
        document.querySelectorAll("#longi-sample-list input[type=checkbox]").forEach((cb) => {
          cb.checked = true;
        });
        _drawLongitudinalPlot();
      });
    }
    if (hideAllBtn) {
      hideAllBtn.addEventListener("click", () => {
        Object.keys(_longiHidden).forEach((id) => {
          _longiHidden[id] = true;
        });
        // Sync checkboxes
        document.querySelectorAll("#longi-sample-list input[type=checkbox]").forEach((cb) => {
          cb.checked = false;
        });
        _drawLongitudinalPlot();
      });
    }
  }

  _drawLongitudinalPlot();
}

function _populateLongiOrgs(visibleOrgs) {
  const container = document.getElementById("longi-org-list");
  const counter = document.getElementById("longi-org-count");
  if (!container) return;

  container.innerHTML = "";
  visibleOrgs.forEach((org) => {
    const selected = _longiSelectedOrgs.has(org);
    const color = _longiOrgColors[org] || "#1565c0";

    const chip = document.createElement("button");
    chip.type = "button";
    chip.title = org;
    chip.style.cssText = [
      "display:inline-flex;align-items:center;gap:5px",
      "padding:3px 8px 3px 6px",
      "border-radius:14px",
      "font-size:0.76em",
      "cursor:pointer",
      "transition:all .12s",
      "white-space:nowrap",
      "max-width:190px",
      selected
        ? `background:${color};color:#fff;border:1.5px solid ${color};font-weight:600`
        : "background:#f0f4f9;color:#455a64;border:1.5px solid #ccd6e8;font-weight:400",
    ].join(";");

    // Colour dot
    const dot = document.createElement("span");
    dot.style.cssText = `display:inline-block;width:8px;height:8px;border-radius:50%;flex-shrink:0;
            background:${selected ? "#fff" : color};border:1px solid rgba(0,0,0,.18)`;

    // Label (truncated)
    const lbl = document.createElement("span");
    lbl.style.cssText = "overflow:hidden;text-overflow:ellipsis;max-width:160px;display:block";
    lbl.textContent = org;

    chip.appendChild(dot);
    chip.appendChild(lbl);

    chip.addEventListener("click", () => {
      if (_longiSelectedOrgs.has(org)) {
        _longiSelectedOrgs.delete(org);
      } else {
        _longiSelectedOrgs.add(org);
      }
      // Re-render chips (keeping search filter)
      const q = (document.getElementById("longi-org-search") || {}).value || "";
      _populateLongiOrgs(
        q.trim() ? _longiOrgList.filter((o) => o.toLowerCase().includes(q.trim().toLowerCase())) : _longiOrgList,
      );
      _drawLongitudinalPlot();
    });

    container.appendChild(chip);
  });

  // Update counter label
  if (counter) {
    const n = _longiSelectedOrgs.size;
    counter.textContent = n ? `— ${n} selected` : "";
  }
}

function _buildLongiSamplePanel(sampleNames) {
  const list = document.getElementById("longi-sample-list");
  if (!list) return;
  list.innerHTML = "";
  sampleNames
    .slice()
    .sort()
    .forEach((id) => {
      if (_longiHidden[id] === undefined) _longiHidden[id] = false;
      const color = sampleColors[id] || "#1565c0";
      const div = document.createElement("div");
      div.style.cssText = "display:flex;align-items:center;gap:5px;margin-bottom:5px";

      const swatch = document.createElement("span");
      swatch.style.cssText = `display:inline-block;width:10px;height:10px;border-radius:50%;
            background:${color};flex-shrink:0;border:1.5px solid rgba(0,0,0,.15)`;

      const lbl = document.createElement("label");
      lbl.style.cssText =
        "font-size:0.78em;cursor:pointer;flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;display:flex;align-items:center;gap:4px";
      lbl.title = id;

      const cb = document.createElement("input");
      cb.type = "checkbox";
      cb.checked = !_longiHidden[id];
      cb.style.cssText = "flex-shrink:0;cursor:pointer";
      cb.addEventListener("change", (e) => {
        _longiHidden[id] = !e.target.checked;
        _drawLongitudinalPlot();
      });

      const txt = document.createTextNode(id.length > 18 ? id.slice(0, 17) + "…" : id);
      lbl.appendChild(cb);
      lbl.appendChild(txt);
      div.appendChild(swatch);
      div.appendChild(lbl);
      list.appendChild(div);
    });
}

function _drawLongitudinalPlot() {
  const wrap = document.getElementById("longi-chart-wrap");
  const noData = document.getElementById("longi-no-data");
  if (!wrap) return;

  wrap.innerHTML = "";

  const ySel = document.getElementById("longi-y-sel");
  const yField = (ySel && ySel.value) || "TASS Score";

  const selectedOrgs = [..._longiSelectedOrgs];
  if (!selectedOrgs.length) {
    if (noData) noData.style.display = "block";
    return;
  }

  // Rebuild time map and run map from live RUN_META
  const timeMap = {};
  const runMap = {}; // sample_name → run label (if present)
  (RUN_META || []).forEach((r) => {
    const d = _parseLongiDate(r.collection_time);
    if (d) timeMap[r.sample_name] = d;
    // Accept "run" (from CSV/xlsx) or "run_id" (from JSON metadata)
    const rv = r.run || r.run_id || null;
    if (rv != null) runMap[r.sample_name] = String(rv);
  });
  const hasRunInfo = Object.keys(runMap).length > 0;

  // All visible samples that have a collection_time (used for zero-fill)
  const allTimedSamples = Object.keys(timeMap).filter((s) => !_longiHidden[s]);

  // Build mol_type lookup for organisms
  const _longiMolType = {};
  filteredData().forEach((r) => {
    if (r["Mol Type"]) _longiMolType[r["Detected Organism"]] = r["Mol Type"];
  });

  // Build one series per selected organism: [{org, pts:[{sample,date,y,run,isZero}]}]
  const series = selectedOrgs
    .map((org) => {
      // Real hits
      const hitPts = filteredData()
        .filter(
          (r) => (r["Detected Organism"] || "") === org && timeMap[r["Specimen ID"]] && !_longiHidden[r["Specimen ID"]],
        )
        .map((r) => ({
          org,
          sample: r["Specimen ID"] || "",
          date: timeMap[r["Specimen ID"]],
          y: parseFloat(r[yField]),
          run: runMap[r["Specimen ID"]] || null,
          isZero: false,
        }))
        .filter((p) => !isNaN(p.y));

      // For runs that have ≥1 real hit, add zero points for same-run samples with no hit
      const hitSamples = new Set(hitPts.map((p) => p.sample));
      const runsWithHits = new Set(hitPts.map((p) => p.run).filter((r) => r != null));

      const zeroPts = [];
      if (hasRunInfo && runsWithHits.size > 0) {
        allTimedSamples.forEach((s) => {
          if (hitSamples.has(s)) return; // already has a real hit
          const r = runMap[s] || null;
          if (!runsWithHits.has(r)) return; // not in a run that detected this org
          zeroPts.push({
            org,
            sample: s,
            date: timeMap[s],
            y: 0,
            run: r,
            isZero: true,
          });
        });
      }

      const pts = [...hitPts, ...zeroPts].sort((a, b) => a.date - b.date || (a.run || "").localeCompare(b.run || ""));
      const molType = (_longiMolType[org] || "").toLowerCase();
      return { org, color: _longiOrgColors[org] || "#1565c0", pts, molType };
    })
    .filter((s) => s.pts.some((p) => !p.isZero)); // keep series only if it has ≥1 real hit

  if (!series.length) {
    if (noData) noData.style.display = "block";
    return;
  }
  if (noData) noData.style.display = "none";

  const allPts = series.flatMap((s) => s.pts.filter((p) => !p.isZero));
  const allPtsInc = series.flatMap((s) => s.pts); // includes zero pts for x-domain
  const mo = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];

  // ── Legend (above chart) ──
  const legendDiv = document.createElement("div");
  legendDiv.style.cssText = "display:flex;flex-wrap:wrap;gap:6px 12px;margin-bottom:6px;font-size:0.76em";
  series.forEach((s) => {
    const mtBadge =
      s.molType === "dna"
        ? ' <span style="display:inline-flex;align-items:center;justify-content:center;width:14px;height:14px;border-radius:50%;background:#1565c0;color:#fff;font-size:8px;font-weight:700;vertical-align:middle;line-height:1" title="DNA pathogen">D</span>'
        : s.molType === "rna"
        ? ' <span style="display:inline-flex;align-items:center;justify-content:center;width:14px;height:14px;border-radius:50%;background:#6a1b9a;color:#fff;font-size:8px;font-weight:700;vertical-align:middle;line-height:1" title="RNA pathogen">R</span>'
        : "";
    const item = document.createElement("span");
    item.style.cssText = "display:inline-flex;align-items:center;gap:4px;color:#333";
    item.innerHTML = `<span style="display:inline-block;width:18px;height:3px;background:${s.color};border-radius:2px;vertical-align:middle"></span>
            <span style="max-width:220px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="${s.org}">${s.org}${mtBadge}</span>`;
    legendDiv.appendChild(item);
  });
  wrap.appendChild(legendDiv);

  // ── SVG chart ──
  // If the container has no layout width yet (tab hidden on first render),
  // defer one animation frame so the browser can compute the real clientWidth.
  const _rawW = wrap.clientWidth || wrap.offsetWidth || 0;
  if (!_rawW) {
    wrap.innerHTML = "";
    requestAnimationFrame(_drawLongitudinalPlot);
    return;
  }
  const W = Math.max(320, _rawW);
  const H = 290;
  const ML = 62,
    MR = 20,
    MT = 14,
    MB = 68;
  const iW = W - ML - MR;
  const iH = H - MT - MB;

  const svg = d3
    .select(wrap)
    .append("svg")
    .attr("viewBox", `0 0 ${W} ${H}`)
    .attr("height", H)
    .style("width", "100%")
    .style("display", "block")
    .style("overflow", "visible");

  const g = svg.append("g").attr("transform", `translate(${ML},${MT})`);

  // Scales (shared across all series)
  // Use allPtsInc (includes zero pts) for x domain so early-run no-detects stay in view
  const xExt = d3.extent(allPtsInc, (p) => p.date);
  const xPad = Math.max((xExt[1] - xExt[0]) * 0.08, 6 * 3600 * 1000);
  const xSc = d3
    .scaleTime()
    .domain([new Date(+xExt[0] - xPad), new Date(+xExt[1] + xPad)])
    .range([0, iW]);

  const yMax = d3.max(allPts, (p) => p.y) || 1;
  const scaleType = (document.getElementById("longi-scale-sel") || {}).value || "linear";
  let ySc;
  if (scaleType === "log") {
    const yMinPos = Math.max(d3.min(allPts, (p) => p.y) || 0.001, 0.001);
    ySc = d3
      .scaleLog()
      .domain([yMinPos * 0.9, yMax * 1.5])
      .range([iH, 0])
      .clamp(true);
  } else if (scaleType === "sqrt") {
    ySc = d3
      .scaleSqrt()
      .domain([0, yMax * 1.15])
      .range([iH, 0])
      .nice();
  } else {
    ySc = d3
      .scaleLinear()
      .domain([0, yMax * 1.15])
      .range([iH, 0])
      .nice();
  }

  // Horizontal grid
  g.append("g")
    .call(d3.axisLeft(ySc).tickSize(-iW).tickFormat("").ticks(5))
    .call((ax) => ax.select(".domain").remove())
    .call((ax) => ax.selectAll("line").attr("stroke", "#e5eaf2").attr("stroke-dasharray", "4,3"));

  // X axis
  g.append("g")
    .attr("transform", `translate(0,${iH})`)
    .call(
      d3
        .axisBottom(xSc)
        .ticks(Math.min(allPts.length + 1, 7))
        .tickFormat((d) => `${mo[d.getMonth()]} ${d.getDate()}, ${d.getFullYear()}`),
    )
    .call((ax) =>
      ax.selectAll("text").attr("transform", "rotate(-35)").style("text-anchor", "end").style("font-size", "0.71em"),
    )
    .call((ax) => ax.select(".domain").attr("stroke", "#ccc"))
    .call((ax) => ax.selectAll(".tick line").attr("stroke", "#ccc"));

  // Y axis
  const yAxis = scaleType === "log" ? d3.axisLeft(ySc).ticks(5, "~s") : d3.axisLeft(ySc).ticks(5);
  g.append("g")
    .call(yAxis)
    .call((ax) => ax.selectAll("text").style("font-size", "0.71em"))
    .call((ax) => ax.select(".domain").attr("stroke", "#ccc"))
    .call((ax) => ax.selectAll(".tick line").attr("stroke", "#ccc"));

  g.append("text")
    .attr("transform", "rotate(-90)")
    .attr("x", -iH / 2)
    .attr("y", -ML + 14)
    .attr("text-anchor", "middle")
    .style("font-size", "0.74em")
    .style("fill", "#555")
    .text(yField);

  // Tooltip
  const tip = d3
    .select(wrap)
    .append("div")
    .style("position", "absolute")
    .style("background", "rgba(22,35,58,0.95)")
    .style("color", "#fff")
    .style("border-radius", "8px")
    .style("padding", "8px 13px")
    .style("font-size", "0.78em")
    .style("pointer-events", "none")
    .style("display", "none")
    .style("white-space", "nowrap")
    .style("z-index", "9999")
    .style("box-shadow", "0 2px 10px rgba(0,0,0,.4)");

  const lineGen = d3
    .line()
    .defined((p) => scaleType !== "log" || p.y > 0)
    .x((p) => xSc(p.date))
    .y((p) => ySc(scaleType === "log" ? Math.max(ySc.domain()[0], p.y) : p.y));

  // Draw one line + dots per organism
  series.forEach((s) => {
    // Lines: if run info is present, segment so points from different runs don't connect
    if (s.pts.length > 1) {
      if (hasRunInfo && s.pts.some((p) => p.run != null)) {
        // Split into contiguous same-run segments and draw each separately
        let i = 0;
        while (i < s.pts.length) {
          const segRun = s.pts[i].run;
          let j = i + 1;
          while (j < s.pts.length && s.pts[j].run === segRun) j++;
          const seg = s.pts.slice(i, j);
          if (seg.length > 1) {
            g.append("path")
              .datum(seg)
              .attr("fill", "none")
              .attr("stroke", s.color)
              .attr("stroke-width", 2.2)
              .attr("opacity", 0.75)
              .attr("d", lineGen);
          }
          i = j;
        }
      } else {
        g.append("path")
          .datum(s.pts)
          .attr("fill", "none")
          .attr("stroke", s.color)
          .attr("stroke-width", 2.2)
          .attr("opacity", 0.75)
          .attr("d", lineGen);
      }
    }

    // Dots — bullseye: outer=organism color, white ring, inner=sample color
    // Zero/absent points stay as hollow dashed circles
    s.pts.forEach((p) => {
      const isZ = p.isZero;
      const cx = xSc(p.date);
      const cy = isZ ? ySc(ySc.domain()[0]) : ySc(p.y);
      const sampleCol = sampleColors[p.sample] || "#90a4ae";

      const dotG = g.append("g").attr("transform", `translate(${cx},${cy})`).style("cursor", "pointer");

      if (isZ) {
        // Hollow dashed circle for "not detected"
        dotG
          .append("circle")
          .attr("r", 5)
          .attr("fill", "#fff")
          .attr("stroke", s.color)
          .attr("stroke-width", 1.5)
          .attr("stroke-dasharray", "3,2")
          .attr("opacity", 0.65);
      } else {
        // Outer: organism color
        dotG
          .append("circle")
          .attr("r", 7)
          .attr("fill", s.color)
          .attr("stroke", "none")
          .style("filter", "drop-shadow(0 1px 3px rgba(0,0,0,.28))");
        // White ring gap
        dotG
          .append("circle")
          .attr("r", 4.5)
          .attr("fill", "#fff")
          .attr("stroke", "none")
          .style("pointer-events", "none");
        // Inner: sample color
        dotG
          .append("circle")
          .attr("r", 2.5)
          .attr("fill", sampleCol)
          .attr("stroke", "none")
          .style("pointer-events", "none");
      }

      dotG
        .on("mouseover", function (event) {
          d3.select(this)
            .transition()
            .duration(80)
            .attr("transform", `translate(${cx},${cy}) scale(${isZ ? 1.4 : 1.35})`);
          const dateStr = `${mo[p.date.getMonth()]} ${p.date.getDate()}, ${p.date.getFullYear()}`;
          const fmt = (v) => (v % 1 === 0 ? v.toLocaleString() : v.toFixed(3));
          tip
            .html(
              `<span style="color:#90caf9;font-weight:700">${p.sample}</span>` +
                (p.run ? `<span style="color:#78909c;font-size:0.88em"> · ${p.run}</span>` : "") +
                `<br><span style="color:#b0bec5;font-size:0.9em">${dateStr}</span><br>` +
                `<span style="color:#cfd8dc;font-size:0.88em">${
                  p.org.length > 35 ? p.org.slice(0, 34) + "…" : p.org
                }</span><br>` +
                (isZ
                  ? `<span style="color:#ef9a9a">Not detected</span>`
                  : `<span style="color:#fff">${yField}:</span> <b>${fmt(p.y)}</b>`),
            )
            .style("display", "block");
        })
        .on("mousemove", function (event) {
          const box = wrap.getBoundingClientRect();
          const tipW = 200;
          let lft = event.clientX - box.left + 14;
          if (lft + tipW > box.width) lft = event.clientX - box.left - tipW - 10;
          tip.style("left", lft + "px").style("top", event.clientY - box.top - 42 + "px");
        })
        .on("mouseout", function () {
          d3.select(this).transition().duration(80).attr("transform", `translate(${cx},${cy}) scale(1)`);
          tip.style("display", "none");
        });
    });
  });

  // Sample name labels — one per unique sample (real hits only), de-overlapped
  // when multiple samples share the same (or very close) x position
  const labelMap = new Map();
  series.forEach((s) =>
    s.pts
      .filter((p) => !p.isZero)
      .forEach((p) => {
        const py = ySc(p.y);
        if (!labelMap.has(p.sample) || py < labelMap.get(p.sample).topY) {
          labelMap.set(p.sample, { date: p.date, topY: py });
        }
      }),
  );

  // Resolve overlapping labels: sort by x, then stagger y for nearby neighbours
  const lblData = [...labelMap.entries()].map(([sample, v]) => ({
    sample,
    x: xSc(v.date),
    labelY: v.topY - 13, // default label y (above topmost dot)
  }));
  lblData.sort((a, b) => a.x - b.x);
  const MIN_X_GAP = 68; // px — below this horizontal distance, stagger vertically
  const Y_STEP = 12; // px — vertical shift per collision level
  for (let i = 1; i < lblData.length; i++) {
    if (lblData[i].x - lblData[i - 1].x < MIN_X_GAP) {
      lblData[i].labelY = lblData[i - 1].labelY - Y_STEP;
    }
  }

  g.selectAll(".longi-lbl-sample")
    .data(lblData)
    .enter()
    .append("text")
    .attr("x", (d) => d.x)
    .attr("y", (d) => d.labelY)
    .attr("text-anchor", "middle")
    .style("font-size", "0.66em")
    .style("fill", "#455a64")
    .style("pointer-events", "none")
    .text((d) => (d.sample.length > 16 ? d.sample.slice(0, 15) + "…" : d.sample));
}

// ── Metadata editing helpers ──────────────────────────────────────────
const _META_NUM_KEYS = new Set(["latitude", "longitude", "depth", "salinity"]);

// Seed a RUN_META record for every sample in the dataset that doesn't have
// one yet, so the user can type metadata in even when none was supplied.
function _ttSeedMetaRowsForAllSamples() {
  const have = new Set(RUN_META.map((r) => r.sample_name));
  const ids = (typeof uniq === "function" ? uniq((DATA || []).map((r) => r["Specimen ID"] || "")) : [])
    .filter(Boolean)
    .sort((a, b) => a.localeCompare(b, undefined, { numeric: true, sensitivity: "base" }));
  let added = 0;
  ids.forEach((id) => {
    if (!have.has(id)) {
      RUN_META.push({ sample_name: id });
      have.add(id);
      added++;
    }
  });
  return added;
}

// Re-render every Run-Metadata-derived view so edited cells flow through
// to the choropleth, geographic/host bar charts, longitudinal & cross-entry
// plots, and the Leaflet marker map. Deliberately does NOT rebuild the
// metadata table itself — that would steal focus and break tab-to-next-cell
// editing. The edited value is already shown in the cell.
function _ttRefreshMetaDerived(changedKey) {
  const _geoCoord = changedKey === "latitude" || changedKey === "longitude";
  const _geoOrigin = changedKey === "sample_origin_country" || changedKey === "sample_origin_state_province_territory";
  if (_geoCoord && typeof _rebuildMapMarkers === "function") {
    try {
      _rebuildMapMarkers();
    } catch (e) {}
    // Refresh the precise-map view (Mapping & Geography sub-tab) so newly
    // added / edited coordinates appear without re-selecting the level.
    if (typeof _geoRedraw === "function") {
      try {
        _geoRedraw();
      } catch (e) {}
    }
  }
  // Backfill country / state / location from coordinates when either the
  // coordinates change OR a geo-origin field is cleared (so deleting "MA"
  // with lat/long present re-detects it). Self-guards: only fills blanks
  // and makes no network call when nothing needs filling.
  if ((_geoCoord || _geoOrigin) && typeof _ttAutofillGeoFromCoords === "function") {
    _ttAutofillGeoFromCoords();
  }
  // Enable/disable the analysis sub-tabs based on the now-current metadata…
  if (typeof _updateMetaSubTabStates === "function") _updateMetaSubTabStates();
  // …and re-draw whichever analysis sub-tab is open (choropleth, geographic
  // comparison, host & disease, longitudinal, cross-entry).
  if (typeof _activeMetaSub !== "undefined" && _activeMetaSub && typeof _switchMetaSub === "function") {
    try {
      _switchMetaSub(_activeMetaSub);
    } catch (e) {}
  }
}

// Commit an edited Run Metadata cell back into its RUN_META record.
function _ttCommitMetaCell(td) {
  const idx = parseInt(td.dataset.metaIdx, 10);
  const key = td.dataset.metaKey;
  const rec = RUN_META[idx];
  if (!rec || !key) return;
  let v = td.innerText.replace(/ /g, " ").trim();
  if (v === "" || v === "—") {
    rec[key] = null;
  } else if (_META_NUM_KEYS.has(key)) {
    const n = parseFloat(v);
    rec[key] = isNaN(n) ? null : n;
  } else {
    rec[key] = v;
  }
  if (typeof _normalizeMetaRecord === "function") _normalizeMetaRecord(rec);
  // An in-place cell edit changes a value without changing the record count,
  // which the grouping model's cache token cannot see. Drop the cache and
  // re-profile so a newly-typed category shows up in the "Group by" bar and
  // the group views immediately.
  try {
    if (window.metaGrouping) {
      window.metaGrouping.refresh();
      if (typeof _mgSyncGroupBar === "function") _mgSyncGroupBar();
      if (typeof _mgBroadcastRedraw === "function") _mgBroadcastRedraw();
    }
  } catch (e) {}
  // Keep SAMPLE_META in sync so any KPI / lookup reading from it sees the edit.
  try {
    if (typeof SAMPLE_META !== "undefined" && SAMPLE_META) {
      SAMPLE_META[rec.sample_name] = SAMPLE_META[rec.sample_name] || { sample_name: rec.sample_name };
      SAMPLE_META[rec.sample_name][key] = rec[key];
      // Editing a specimen cell in the metadata table regroups the run, so the
      // specimen caches must be dropped (the entry count is unchanged, which a
      // size-only signature would not notice).
      if (typeof _noteSampleMetaChanged === "function") _noteSampleMetaChanged();
    }
  } catch (e) {}
  // When lat/lon changes and both coordinates are now valid, clear the
  // auto-derived geographic fields so _ttAutofillGeoFromCoords will
  // refill them from the new position instead of skipping because they
  // are already set.  Without this, moving a sample to a different
  // country never updates the choropleth.
  if (key === "latitude" || key === "longitude") {
    const _lat = parseFloat(rec.latitude),
      _lon = parseFloat(rec.longitude);
    if (!isNaN(_lat) && !isNaN(_lon)) {
      rec.sample_origin_country = null;
      rec.sample_origin_state_province_territory = null;
      rec.location = null;
      try {
        if (typeof SAMPLE_META !== "undefined" && SAMPLE_META && SAMPLE_META[rec.sample_name]) {
          SAMPLE_META[rec.sample_name].sample_origin_country = null;
          SAMPLE_META[rec.sample_name].sample_origin_state_province_territory = null;
          SAMPLE_META[rec.sample_name].location = null;
        }
      } catch (e) {}
    }
  }
  // Propagate the edit to all metadata plots + the marker map WITHOUT
  // rebuilding the table (so the user can keep editing adjacent cells).
  _ttRefreshMetaDerived(key);
}

// Wire the metadata toolbar buttons exactly once.
let _runmetaToolbarWired = false;
function _wireRunMetaToolbar() {
  if (_runmetaToolbarWired) return;
  _runmetaToolbarWired = true;
  const addColBtn = document.getElementById("runmeta-add-col");
  const addRowsBtn = document.getElementById("runmeta-add-rows");
  const exportBtn = document.getElementById("runmeta-export-xlsx");
  const fillGeoBtn = document.getElementById("runmeta-fill-geo");
  // Opt-in "only filtered samples" scoping for this table.
  const scopeCb = document.getElementById("runmeta-filter-scope");
  if (scopeCb)
    scopeCb.addEventListener("change", () => {
      if (typeof ttBusyRun === "function") ttBusyRun("Filtering metadata…", _buildRunMetaTable);
      else _buildRunMetaTable();
    });
  if (fillGeoBtn)
    fillGeoBtn.addEventListener("click", () => {
      if (typeof _ttAutofillGeoFromCoords === "function") _ttAutofillGeoFromCoords({ notify: true });
    });
  if (addColBtn)
    addColBtn.addEventListener("click", () => {
      const raw = (window.prompt("New metadata column name (e.g. notes, collected_by):", "") || "").trim();
      if (!raw) return;
      const key = raw.toLowerCase().replace(/\s+/g, "_");
      const _isHidden = (TT_ANNOT.hiddenCols || []).includes(key);
      const exists = !_isHidden && (RUN_META.some((r) => key in r) || TT_ANNOT.metaCols.includes(key));
      if (exists) {
        alert("A metadata column “" + raw + "” already exists.");
        return;
      }
      // Re-adding a previously removed column simply un-hides it.
      if (_isHidden) TT_ANNOT.hiddenCols = TT_ANNOT.hiddenCols.filter((k) => k !== key);
      if (!TT_ANNOT.metaCols.includes(key)) TT_ANNOT.metaCols.push(key);
      if (RUN_META.length === 0) _ttSeedMetaRowsForAllSamples();
      // Ensure the key is present on records so the column renders.
      RUN_META.forEach((r) => {
        if (!(key in r)) r[key] = null;
      });
      const rmBtn = document.getElementById("runmeta-tab-btn");
      if (rmBtn) rmBtn.classList.toggle("hidden", RUN_META.length === 0);
      _buildRunMetaTable();
    });
  if (addRowsBtn)
    addRowsBtn.addEventListener("click", () => {
      const added = _ttSeedMetaRowsForAllSamples();
      const rmBtn = document.getElementById("runmeta-tab-btn");
      if (rmBtn) rmBtn.classList.toggle("hidden", RUN_META.length === 0);
      _buildRunMetaTable();
      if (typeof _updateMetaSubTabStates === "function") _updateMetaSubTabStates();
      if (!added) alert("Every sample already has a metadata row.");
    });
  if (exportBtn)
    exportBtn.addEventListener("click", () => {
      const tbl = document.getElementById("runmeta-table");
      if (!tbl || !document.getElementById("runmeta-body").children.length) {
        alert("No metadata to export yet. Add rows / columns first.");
        return;
      }
      // Pagination hides off-page rows with display:none, which the table
      // exporter skips. Clone the table with every row made visible so the
      // export always contains all samples, not just the current page.
      const exportTbl = tbl.cloneNode(true);
      exportTbl.querySelectorAll("tr").forEach((tr) => (tr.style.display = ""));
      // Drop the header ✕ remove buttons so they don't leak into column names.
      exportTbl.querySelectorAll(".runmeta-col-remove").forEach((el) => el.remove());
      if (typeof _openTableExport === "function") _openTableExport(exportTbl);
    });
}

function _buildRunMetaTable() {
  _wireRunMetaToolbar();
  const hdrRow = document.getElementById("runmeta-header-row");
  const tbody = document.getElementById("runmeta-body");
  // No metadata (e.g. after "Start empty" or loading an empty state):
  // clear any previously-rendered rows instead of leaving them stale.
  if (!RUN_META || RUN_META.length === 0) {
    if (hdrRow) hdrRow.innerHTML = "";
    if (tbody) tbody.innerHTML = "";
    const pager = document.getElementById("runmeta-pager");
    if (pager) {
      pager.innerHTML = "";
      pager.style.display = "none";
    }
    return;
  }
  if (!hdrRow || !tbody) return;

  // Discover all keys present across every record, keeping sample_name first
  const keySet = new Set();
  RUN_META.forEach((r) =>
    Object.keys(r).forEach((k) => {
      if (k !== "sample_name") keySet.add(k);
    }),
  );
  // Always surface user-added metadata columns, even if still empty.
  TT_ANNOT.metaCols.forEach((k) => keySet.add(k));
  // Always surface the standard metadata columns (country, lat/long,
  // collection_time, host / disease, …) even when the samplesheet supplied
  // none of them, so users have empty cells to fill in.
  _STANDARD_META_COLS.forEach((k) => keySet.add(k));
  // Preferred order for known fields, then any extras alphabetically
  const KNOWN_ORDER = [
    "run_id",
    "sample_id",
    "organism",
    "host_scientific_name",
    "host_disease",
    "environmental_site",
    "sample_origin_country",
    "sample_origin_state_province_territory",
    "latitude",
    "longitude",
    "depth",
    "salinity",
    "collection_time",
    "location",
    "submitter_organization_name",
    "library_preparation_kit",
    "sequencing_instrument",
    "sequencing_platform",
    "sequencing_protocol_primer_set",
  ];
  const ordered = [
    "sample_name",
    ...KNOWN_ORDER.filter((k) => keySet.has(k)),
    ...[...keySet].filter((k) => !KNOWN_ORDER.includes(k)).sort(),
  ];
  // Only show columns that have at least one non-null value — EXCEPT the
  // sample_name key, the standard metadata columns, and any user-added
  // columns (TT_ANNOT.metaCols), which are always shown so the user has
  // empty cells to type into.
  const _alwaysShow = new Set(["sample_name", ..._STANDARD_META_COLS, ...TT_ANNOT.metaCols]);
  const _hidden = new Set(TT_ANNOT.hiddenCols || []);
  const activeCols = ordered.filter(
    (k) =>
      // sample_name is never removable; every other column can be hidden
      // via the header ✕ (stored in TT_ANNOT.hiddenCols).
      (k === "sample_name" || !_hidden.has(k)) &&
      (_alwaysShow.has(k) ||
        RUN_META.some((r) => {
          const _v = r[k];
          return _v != null && _v !== "" && (Array.isArray(_v) || typeof _v !== "object");
        })),
  );
  hdrRow.innerHTML = activeCols
    .map((k) => {
      // Every column except the sample_name key gets a ✕ to remove it.
      const removeBtn =
        k !== "sample_name"
          ? `<span class="runmeta-col-remove" data-col="${k}" title="Remove this column" role="button" tabindex="0">&times;</span>`
          : "";
      return `<th style="background:#1565C0;color:#fff;padding:6px 10px;text-align:left;white-space:nowrap;font-weight:600"><span class="runmeta-col-label">${_metaKeyLabel(
        k,
      )}</span>${removeBtn}</th>`;
    })
    .join("");

  // Wire the per-column remove (✕) buttons.
  hdrRow.querySelectorAll(".runmeta-col-remove").forEach((el) => {
    const _do = () => _ttRemoveMetaColumn(el.dataset.col);
    el.addEventListener("click", (e) => {
      e.stopPropagation();
      _do();
    });
    el.addEventListener("keydown", (e) => {
      if (e.key === "Enter" || e.key === " ") {
        e.preventDefault();
        _do();
      }
    });
  });

  // ── Optional filter scoping ────────────────────────────────────────
  // Off by default (this table is the run inventory). When on, keep only the
  // samples surviving the sidebar filters. CRITICAL: the original RUN_META
  // index travels with each row — the editable cells address records by index
  // via data-meta-idx, so renumbering them would write edits into the wrong
  // record.
  const _scopeEl = document.getElementById("runmeta-filter-scope");
  const _scoped = !!(_scopeEl && _scopeEl.checked);
  const _matched = _scoped && typeof _ttFilterMatchedKeys === "function" ? _ttFilterMatchedKeys() : null;
  const _visible = RUN_META.map((rec, i) => ({ rec, i })).filter(
    ({ rec }) =>
      !_matched || (typeof _ttSampleMatchesFilter === "function" && _ttSampleMatchesFilter(rec.sample_name, _matched)),
  );
  const _cntEl = document.getElementById("runmeta-filter-scope-count");
  if (_cntEl) _cntEl.textContent = _scoped ? `(${_visible.length} of ${RUN_META.length})` : "";

  let _rowN = 0; // display position, for zebra striping (RUN_META index may skip)
  tbody.innerHTML = _visible
    .map(({ rec, i }) => {
      const isHighlighted = rec.sample_name === _runmetaHighlightSample;
      const bg = isHighlighted ? "#fff9c4" : _rowN++ % 2 === 0 ? "#f0f6ff" : "#fff";
      const border = isHighlighted ? "2px solid #f9a825" : "none";
      const cells = activeCols
        .map((k) => {
          const v = rec[k];
          const disp = Array.isArray(v) ? v.join(", ") : v != null && v !== "" ? v : "";
          const title = Array.isArray(v)
            ? v.join(", ").replace(/"/g, "&quot;")
            : v != null
            ? String(v).replace(/"/g, "&quot;")
            : "";
          // sample_name is the record key — not editable. Every other column
          // is contenteditable; edits commit on blur via _ttCommitMetaCell.
          const editable = k !== "sample_name";
          const cls = editable ? ' class="runmeta-editable"' : "";
          const editAttrs = editable
            ? ` contenteditable="true" spellcheck="false" data-meta-idx="${i}" data-meta-key="${k}"`
            : "";
          const dispEsc = String(disp).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
          const mergedBadge =
            k === "sample_name" && typeof _mergedSampleBadgeHTML === "function"
              ? _mergedSampleBadgeHTML(rec.sample_name)
              : "";
          return `<td${cls}${editAttrs} style="padding:5px 10px;border-bottom:1px solid #ddd;color:#222;max-width:220px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="${title}">${dispEsc}${mergedBadge}</td>`;
        })
        .join("");
      return `<tr id="runmeta-row-${CSS.escape(
        rec.sample_name,
      )}" style="background:${bg};outline:${border}">${cells}</tr>`;
    })
    .join("");

  // Attach edit-commit handlers to the editable cells.
  tbody.querySelectorAll("td.runmeta-editable").forEach((td) => {
    td.addEventListener("blur", () => _ttCommitMetaCell(td));
    td.addEventListener("keydown", (e) => {
      if (e.key === "Enter") {
        e.preventDefault();
        td.blur();
      } else if (e.key === "Escape") {
        td.blur();
      }
    });
  });

  // ── Pagination ──────────────────────────────────────────────────────
  // If a sample is highlighted, jump to the page that contains it so the
  // scroll-into-view below lands on a visible row.
  if (_runmetaHighlightSample) {
    const _hi = RUN_META.findIndex((r) => r.sample_name === _runmetaHighlightSample);
    if (_hi >= 0) _runmetaPage = Math.floor(_hi / _RUNMETA_PAGE_SIZE) + 1;
  }
  _applyRunMetaPage();

  // Scroll highlighted row into view
  if (_runmetaHighlightSample) {
    const row = document.getElementById(`runmeta-row-${CSS.escape(_runmetaHighlightSample)}`);
    if (row) setTimeout(() => row.scrollIntoView({ behavior: "smooth", block: "center" }), 50);
  }

  // Point the user at the "Rows for all samples" button whenever some
  // samples still have no metadata row.
  _updateRunMetaAddRowsHint();
}

// Count samples in DATA that don't yet have a metadata row (non-mutating).
function _runMetaMissingSampleCount() {
  const have = new Set((RUN_META || []).map((r) => r.sample_name));
  const ids = (typeof uniq === "function" ? uniq((DATA || []).map((r) => r["Specimen ID"] || "")) : []).filter(Boolean);
  return ids.reduce((n, id) => (have.has(id) ? n : n + 1), 0);
}

// Pulse-highlight the "Rows for all samples" button (and show a callout)
// when one or more samples are missing a metadata row. Clears once every
// sample has a row.
function _updateRunMetaAddRowsHint() {
  const btn = document.getElementById("runmeta-add-rows");
  const hint = document.getElementById("runmeta-add-rows-hint");
  const missing = _runMetaMissingSampleCount();
  if (btn) btn.classList.toggle("runmeta-pulse", missing > 0);
  if (hint) {
    hint.style.display = missing > 0 ? "inline-flex" : "none";
    if (missing > 0) {
      hint.querySelector("span").textContent = `${missing} sample${
        missing === 1 ? "" : "s"
      } missing a row — click to add`;
    }
  }
}

// Remove (hide) a metadata column from the table. sample_name can't be
// removed. User-added columns are also dropped from TT_ANNOT.metaCols.
function _ttRemoveMetaColumn(key) {
  if (!key || key === "sample_name") return;
  const label = typeof _metaKeyLabel === "function" ? _metaKeyLabel(key) : key;
  if (!window.confirm(`Remove the “${label}” column from the metadata table?`)) return;
  TT_ANNOT.hiddenCols = TT_ANNOT.hiddenCols || [];
  if (!TT_ANNOT.hiddenCols.includes(key)) TT_ANNOT.hiddenCols.push(key);
  // Drop from user-added columns so it isn't re-surfaced as always-shown.
  TT_ANNOT.metaCols = (TT_ANNOT.metaCols || []).filter((k) => k !== key);
  _buildRunMetaTable();
}

// Show only the current page of metadata rows and (re)render the pager.
// All rows remain in the DOM; off-page rows are simply hidden.
function _applyRunMetaPage() {
  const tbody = document.getElementById("runmeta-body");
  if (!tbody) return;
  const rows = Array.from(tbody.children);
  const total = rows.length;
  const totalPages = Math.max(1, Math.ceil(total / _RUNMETA_PAGE_SIZE));
  if (_runmetaPage > totalPages) _runmetaPage = totalPages;
  if (_runmetaPage < 1) _runmetaPage = 1;
  const start = (_runmetaPage - 1) * _RUNMETA_PAGE_SIZE;
  const end = start + _RUNMETA_PAGE_SIZE;
  rows.forEach((tr, i) => {
    tr.style.display = i >= start && i < end ? "" : "none";
  });
  _renderRunMetaPager(total, totalPages, start, end);
}

// Navigate to a specific page of the metadata table.
function _gotoRunMetaPage(n) {
  _runmetaPage = n;
  _applyRunMetaPage();
}

// Build the pager control below the metadata table. Hidden when a single
// page holds every row.
function _renderRunMetaPager(total, totalPages, start, end) {
  const pager = document.getElementById("runmeta-pager");
  if (!pager) return;
  if (total <= _RUNMETA_PAGE_SIZE) {
    pager.innerHTML = "";
    pager.style.display = "none";
    return;
  }
  pager.style.display = "flex";
  const shownFrom = start + 1;
  const shownTo = Math.min(end, total);
  const btn = (label, page, disabled, title) =>
    `<button type="button" class="runmeta-page-btn" data-page="${page}"${disabled ? " disabled" : ""}${
      title ? ` title="${title}"` : ""
    }>${label}</button>`;
  pager.innerHTML =
    `<span style="color:#555;margin-right:auto">Showing ${shownFrom}–${shownTo} of ${total} samples</span>` +
    btn('<i class="fas fa-angle-double-left"></i>', 1, _runmetaPage === 1, "First page") +
    btn('<i class="fas fa-angle-left"></i>', _runmetaPage - 1, _runmetaPage === 1, "Previous page") +
    `<span style="padding:0 0.6em;color:#333">Page ${_runmetaPage} / ${totalPages}</span>` +
    btn('<i class="fas fa-angle-right"></i>', _runmetaPage + 1, _runmetaPage === totalPages, "Next page") +
    btn('<i class="fas fa-angle-double-right"></i>', totalPages, _runmetaPage === totalPages, "Last page");
  pager.querySelectorAll(".runmeta-page-btn").forEach((b) => {
    b.addEventListener("click", () => {
      const p = parseInt(b.dataset.page, 10);
      if (!isNaN(p)) _gotoRunMetaPage(p);
    });
  });
}
