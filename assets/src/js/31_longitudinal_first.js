/* ═══════════════════════════════════════════════════════════════════════════
       -  §  LONGITUDINAL ANALYSIS  —  FIRST DEFINITION  (shadowed below)
       -     This block is overridden by the second copy further down (~line
       -     9290+); since both define the same function names, the later
       -     definitions win. Keeping it as a visual landmark until somebody
       -     dedupes.
       -     Functions: _parseLongiDate, _buildLongitudinalSection,
       -     _populateLongiOrgs, _buildLongiSamplePanel, _drawLongitudinalPlot.
═══════════════════════════════════════════════════════════════════════════ */

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
