/* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB: EXPLORE            (data-tab="explore")
       -     drawExplore() — top-level dispatcher; reads the explore-mode
       -     selector (bubble / chord / radar / lollipop) and calls one of:
       -        _drawExploreBubble()   — coverage % vs TASS scatter, sized by
       -                                 reads aligned / k2 reads. Scale type
       -                                 (linear / log / sqrt) shared by both
       -                                 axes.
       -        _drawExploreChord()    — sample × organism flow chart.
       -        _drawExploreRadar()    — multi-axis polygon per sample.
       -        _drawExploreLollipop() — top-N sticks comparing samples.
       -        _drawExploreCorrelogram() — Pearson correlation grid across
       -                                    user-picked numeric metrics; uses
       -                                    the same `fd` so it auto-syncs to
       -                                    sidebar filters and sample renames.
       -     _exploreColorScale() centralizes the d3 color domain shared by
       -     the sub-plots.
═══════════════════════════════════════════════════════════════════════════ */
function drawExplore() {
  const fd = filteredData();
  if (!fd.length) return;
  _drawExploreBubble(fd);
  _drawExploreRadar(fd);
  _drawExploreChord(fd);
  _drawExploreLollipop(fd);
  _drawExploreCorrelogram(fd);
}

function _exploreColorScale(fd) {
  const by = (document.getElementById("explore-color-by") || { value: "sample" }).value || "sample";
  const field = by === "sample" ? "Specimen ID" : by === "genus" ? "Genus" : "Microbial Category";
  const keys = [...new Set(fd.map((r) => r[field] || "Unknown"))];
  const pal = [
    "#4e79a7",
    "#f28e2b",
    "#e15759",
    "#76b7b2",
    "#59a14f",
    "#edc948",
    "#b07aa1",
    "#ff9da7",
    "#9c755f",
    "#bab0ac",
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
  ];
  return { field, scale: d3.scaleOrdinal().domain(keys).range(pal), keys };
}

function _drawExploreBubble(fd) {
  const wrap = document.getElementById("explore-bubble-wrap");
  if (!wrap) return;
  wrap.innerHTML = "";
  wrap.style.position = "relative"; // needed for fixed tip positioning

  const W = wrap.clientWidth || 700,
    H = 400;
  const mL = 60,
    mT = 20,
    mR = 165,
    mB = 55;
  const iW = W - mL - mR,
    iH = H - mT - mB;

  const { field, scale: col, keys } = _exploreColorScale(fd);
  const sizeBy = (document.getElementById("explore-size-by") || { value: "reads" }).value;
  const sizeField = sizeBy === "k2" ? "K2 Reads" : "# Reads Aligned";

  const points = fd
    .map((r) => ({
      x: parseFloat(r["Breadth %"] || r["Coverage"] || 0),
      y: parseFloat(r["TASS Score"] || 0),
      r: Math.max(4, Math.sqrt(Math.max(0, parseFloat(r[sizeField] || 0) + 1)) * 0.4),
      label: r["Detected Organism"] || "",
      color: col(r[field] || "Unknown"),
      raw: r,
    }))
    .filter((p) => !isNaN(p.x) && !isNaN(p.y));

  // Axis scale type for both x (Coverage %) and y (TASS Score).
  // Default linear; log/sqrt are available for spreading out the
  // dense low-end common in metagenomics data.
  const expScaleType = (document.getElementById("explore-bubble-scale") || {}).value || "linear";
  const xDataMax = Math.max(d3.max(points, (p) => p.x) || 1, 1);
  function _buildExpScale(type, max, range) {
    if (type === "log")
      return d3
        .scaleLog()
        .domain([0.01, Math.max(max, 0.02)])
        .range(range)
        .clamp(true);
    if (type === "sqrt") return d3.scaleSqrt().domain([0, max]).range(range).nice();
    return d3.scaleLinear().domain([0, max]).range(range).nice();
  }
  const xBase = _buildExpScale(expScaleType, xDataMax, [0, iW]);
  const yBase = _buildExpScale(expScaleType, 100, [iH, 0]);

  const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H);
  const g = svg.append("g").attr("transform", `translate(${mL},${mT})`);

  // Clip path so circles don't overflow on zoom
  svg.append("defs").append("clipPath").attr("id", "bubble-clip").append("rect").attr("width", iW).attr("height", iH);

  // Axes (will be rescaled on zoom)
  const xAxisG = g
    .append("g")
    .attr("transform", `translate(0,${iH})`)
    .call(
      d3
        .axisBottom(xBase)
        .ticks(6)
        .tickFormat((d) => d + "%"),
    );
  const yAxisG = g.append("g").call(d3.axisLeft(yBase).ticks(6));

  svg
    .append("text")
    .attr("x", mL + iW / 2)
    .attr("y", H - 6)
    .attr("text-anchor", "middle")
    .style("font-size", "11px")
    .text("Coverage (Breadth %)");
  svg
    .append("text")
    .attr("transform", `translate(14,${mT + iH / 2})rotate(-90)`)
    .attr("text-anchor", "middle")
    .style("font-size", "11px")
    .text("TASS Score");

  // Gridlines group (rescaled on zoom)
  const gridX = g
    .append("g")
    .attr("class", "grid")
    .attr("transform", `translate(0,${iH})`)
    .call(d3.axisBottom(xBase).ticks(6).tickSize(-iH).tickFormat(""));
  gridX.selectAll("line").style("stroke", "#eee").style("stroke-dasharray", "3,3");
  gridX.select(".domain").remove();
  const gridY = g.append("g").attr("class", "grid").call(d3.axisLeft(yBase).ticks(6).tickSize(-iW).tickFormat(""));
  gridY.selectAll("line").style("stroke", "#eee").style("stroke-dasharray", "3,3");
  gridY.select(".domain").remove();

  // Clipped circle group
  const circleG = g.append("g").attr("clip-path", "url(#bubble-clip)");

  // On-chart hover label (follows current zoom scales)
  let curX = xBase;
  let curY = yBase;
  const hoverG = g.append("g").style("pointer-events", "none").style("display", "none");
  const hoverText = hoverG
    .append("text")
    .style("font-size", "10px")
    .attr("fill", "#222")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "central");
  const hoverBg = hoverG
    .insert("rect", "text")
    .attr("fill", "rgba(255,255,255,0.9)")
    .attr("stroke", "#cbd5e1")
    .attr("stroke-width", 0.8)
    .attr("rx", 3)
    .attr("ry", 3);

  // Fixed tooltip anchored to document body (fixes offset on nested scroll)
  const tip = d3
    .select("body")
    .selectAll(".explore-bubble-tip")
    .data([1])
    .join("div")
    .attr("class", "explore-bubble-tip")
    .style("position", "fixed")
    .style("background", "rgba(20,20,30,.88)")
    .style("color", "#fff")
    .style("padding", ".4em .7em")
    .style("border-radius", "5px")
    .style("font-size", "11.5px")
    .style("pointer-events", "none")
    .style("display", "none")
    .style("max-width", "240px")
    .style("line-height", "1.6")
    .style("z-index", "9999")
    .style("box-shadow", "0 2px 8px rgba(0,0,0,.35)");

  circleG
    .selectAll("circle")
    .data(points)
    .join("circle")
    .attr("cx", (p) => xBase(p.x))
    .attr("cy", (p) => yBase(p.y))
    .attr("r", (p) => p.r)
    .attr("fill", (p) => p.color)
    .attr("opacity", 0.75)
    .attr("stroke", "#fff")
    .attr("stroke-width", 0.8)
    .style("cursor", "crosshair")
    .on("mouseover", (ev, p) => {
      d3.select(ev.currentTarget)
        .attr("stroke", "#111")
        .attr("stroke-width", 1.4)
        .attr("opacity", 1)
        .attr("r", p.r + 2);
      hoverText.text(p.label || "");
      const bb = hoverText.node().getBBox();
      hoverBg
        .attr("width", bb.width + 10)
        .attr("height", bb.height + 6)
        .attr("x", bb.x - 5)
        .attr("y", bb.y - 3);
      const hx = curX(p.x);
      const hy = curY(p.y) - 14;
      hoverG.attr("transform", `translate(${hx},${hy})`);
      hoverG.style("display", "block");
      const r = p.raw;
      tip
        .style("display", "block")
        .html(
          `<b style="font-size:12px">${p.label}</b><br>` +
            `<span style="color:#adf">Sample:</span> ${r["Specimen ID"] || "—"}<br>` +
            `<span style="color:#adf">TASS Score:</span> ${p.y.toFixed(2)}<br>` +
            `<span style="color:#adf">Coverage:</span> ${p.x.toFixed(2)}%<br>` +
            `<span style="color:#adf">${sizeField}:</span> ${(+r[sizeField] || 0).toLocaleString()}<br>` +
            `<span style="color:#adf">Category:</span> ${r["Microbial Category"] || "—"}<br>` +
            `<span style="color:#adf">Mean Depth:</span> ${parseFloat(r["Mean Depth"] || 0).toFixed(3)}<br>` +
            `<span style="color:#adf">Passes Threshold:</span> ${r["Passes Threshold"] || "—"}`,
        );
    })
    .on("mousemove", (ev) => {
      tip.style("left", ev.clientX + 14 + "px").style("top", ev.clientY - 10 + "px");
      const target = d3.select(ev.currentTarget).datum();
      if (target) {
        const hx = curX(target.x);
        const hy = curY(target.y) - 14;
        hoverG.attr("transform", `translate(${hx},${hy})`);
      }
    })
    .on("mouseout", (ev, p) => {
      d3.select(ev.currentTarget).attr("stroke", "#fff").attr("stroke-width", 0.8).attr("opacity", 0.75).attr("r", p.r);
      hoverG.style("display", "none");
      tip.style("display", "none");
    });

  // D3 zoom — scroll to zoom, drag to pan, dbl-click to reset
  const zoom = d3
    .zoom()
    .scaleExtent([0.5, 20])
    .on("zoom", (event) => {
      const t = event.transform;
      const newX = t.rescaleX(xBase);
      const newY = t.rescaleY(yBase);
      curX = newX;
      curY = newY;
      xAxisG.call(
        d3
          .axisBottom(newX)
          .ticks(6)
          .tickFormat((d) => d + "%"),
      );
      yAxisG.call(d3.axisLeft(newY).ticks(6));
      gridX.call(d3.axisBottom(newX).ticks(6).tickSize(-iH).tickFormat(""));
      gridX.selectAll("line").style("stroke", "#eee").style("stroke-dasharray", "3,3");
      gridX.select(".domain").remove();
      gridY.call(d3.axisLeft(newY).ticks(6).tickSize(-iW).tickFormat(""));
      gridY.selectAll("line").style("stroke", "#eee").style("stroke-dasharray", "3,3");
      gridY.select(".domain").remove();
      circleG
        .selectAll("circle")
        .attr("cx", (p) => newX(p.x))
        .attr("cy", (p) => newY(p.y));
    });
  const zoomRect = g
    .insert("rect", ":first-child")
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", iW)
    .attr("height", iH)
    .attr("fill", "none")
    .attr("pointer-events", "all");
  zoomRect.call(zoom);
  svg.on("dblclick.zoom", () => svg.transition().duration(350).call(zoom.transform, d3.zoomIdentity));

  // Zoom hint
  svg
    .append("text")
    .attr("x", mL + iW)
    .attr("y", mT - 5)
    .attr("text-anchor", "end")
    .style("font-size", "9px")
    .attr("fill", "#bbb")
    .text("scroll to zoom · drag to pan · dbl-click reset");

  // Legend (right side)
  const lg = svg.append("g").attr("transform", `translate(${W - mR + 12},${mT})`);
  keys.slice(0, 14).forEach((k, i) => {
    lg.append("circle")
      .attr("cx", 6)
      .attr("cy", i * 18 + 6)
      .attr("r", 5)
      .attr("fill", col(k));
    lg.append("text")
      .attr("x", 14)
      .attr("y", i * 18 + 10)
      .style("font-size", "10px")
      .text(k.length > 18 ? k.slice(0, 17) + "…" : k);
  });
}

function _drawExploreChord(fd) {
  const wrap = document.getElementById("explore-chord-wrap");
  const legWrap = document.getElementById("explore-chord-legend");
  if (!wrap) return;
  wrap.innerHTML = "";
  if (legWrap) legWrap.innerHTML = "";

  const groupBy = (document.getElementById("explore-chord-group") || { value: "category" }).value;
  const valueBy = (document.getElementById("explore-chord-value") || { value: "reads" }).value;
  const valField = valueBy === "tass" ? "TASS Score" : "# Reads Aligned";
  const catField = groupBy === "genus" ? "Genus" : "Microbial Category";

  const samples = _orderedSamples([...new Set(fd.map((r) => r["Specimen ID"] || ""))].filter(Boolean));
  const cats = [...new Set(fd.map((r) => r[catField] || "Other"))].filter(Boolean).sort();
  if (!samples.length || !cats.length) return;

  // Combined node list: samples first, then a spacer gap, then categories
  const GAP_KEY = "__gap__";
  const nodes = [...samples, GAP_KEY, ...cats];
  const gapIndex = samples.length;
  const n = nodes.length;
  const idxOf = (k) => nodes.indexOf(k);

  // Build n×n flow matrix
  const mat = Array.from({ length: n }, () => new Array(n).fill(0));
  fd.forEach((r) => {
    const s = r["Specimen ID"] || "";
    const c = r[catField] || "Other";
    if (!s || !c) return;
    const si = idxOf(s),
      ci = idxOf(c);
    if (si < 0 || ci < 0) return;
    const v = Math.max(0, parseFloat(r[valField]) || 0);
    mat[si][ci] += v;
    mat[ci][si] += v;
  });

  const wrapW = wrap.clientWidth || 360;
  const SIZE = Math.max(280, Math.min(520, Math.floor(wrapW * 0.92)));
  const padX = 80;
  const padY = 30;
  const SVGW = Math.max(wrapW, SIZE + padX * 2);
  const SVGH = SIZE + padY * 2;
  const cx = SVGW / 2;
  const cy = SVGH / 2;
  const R = SIZE * 0.34;
  const innerR = SIZE * 0.28;
  const pal20 = [
    "#4e79a7",
    "#f28e2b",
    "#e15759",
    "#76b7b2",
    "#59a14f",
    "#edc948",
    "#b07aa1",
    "#ff9da7",
    "#9c755f",
    "#bab0ac",
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
  ];
  const colorOf = (i) => pal20[i % pal20.length];

  // Add a small self-value to the gap node so it renders a subtle spacer
  const totalFlow = mat.reduce((sum, row) => sum + row.reduce((s, v) => s + v, 0), 0);
  const gapVal = Math.max(1, totalFlow * 0.03);
  mat[gapIndex][gapIndex] = gapVal;

  const chord = d3.chord().padAngle(0.04).sortSubgroups(d3.descending)(mat);
  const arc = d3.arc().innerRadius(innerR).outerRadius(R);
  const ribbon = d3.ribbon().radius(innerR);

  const svg = d3
    .select(wrap)
    .append("svg")
    .attr("width", SVGW)
    .attr("height", SVGH)
    .style("display", "block")
    .style("margin", "0 auto");
  const g = svg.append("g").attr("transform", `translate(${cx},${cy})`);

  // Tooltip (body-fixed)
  const chTip = d3
    .select("body")
    .selectAll(".explore-chord-tip")
    .data([1])
    .join("div")
    .attr("class", "explore-chord-tip")
    .style("position", "fixed")
    .style("background", "rgba(20,20,30,.88)")
    .style("color", "#fff")
    .style("padding", ".35em .65em")
    .style("border-radius", "5px")
    .style("font-size", "11.5px")
    .style("pointer-events", "none")
    .style("display", "none")
    .style("z-index", "9999")
    .style("max-width", "220px")
    .style("line-height", "1.6")
    .style("box-shadow", "0 2px 8px rgba(0,0,0,.35)");
  const chMove = (ev) => chTip.style("left", ev.clientX + 14 + "px").style("top", ev.clientY - 10 + "px");
  const chOut = () => chTip.style("display", "none");

  // Arcs (outer ring)
  let labelSel = null;
  const arcSel = g
    .append("g")
    .selectAll("path")
    .data(chord.groups.filter((d) => d.index !== gapIndex))
    .join("path")
    .attr("d", arc)
    .attr("fill", (d) => colorOf(d.index))
    .attr("stroke", "#fff")
    .attr("stroke-width", 0.5)
    .attr("opacity", 0.85)
    .on("mouseover", (ev, d) => {
      if (d.index === gapIndex) return;
      chTip
        .style("display", "block")
        .html(`<b>${nodes[d.index]}</b><br>Sum (of Filtered) ${valField}: ${d3.format(",.0f")(d.value)}`);
      // Dim unrelated ribbons
      g.selectAll(".chord-ribbon").attr("opacity", (ch) =>
        ch.source.index === d.index || ch.target.index === d.index ? 0.75 : 0.05,
      );
      const connected = new Set([d.index]);
      chord.forEach((ch) => {
        if (ch.source.index === d.index || ch.target.index === d.index) {
          connected.add(ch.source.index);
          connected.add(ch.target.index);
        }
      });
      arcSel.attr("opacity", (a) => (connected.has(a.index) ? 1 : 0.2));
      if (labelSel) labelSel.attr("opacity", (a) => (connected.has(a.d.index) ? 1 : 0.15));
    })
    .on("mousemove", chMove)
    .on("mouseout", () => {
      g.selectAll(".chord-ribbon").attr("opacity", 0.6);
      arcSel.attr("opacity", 0.85);
      if (labelSel) labelSel.attr("opacity", 1);
      chOut();
    });

  // Arc labels
  const labelMinAngle = 0.18;
  const labelGap = 0.16;
  const labelR = R + 14;
  const maxX = SVGW / 2 - 6;
  const minX = -SVGW / 2 + 6;
  const maxY = SVGH / 2 - 6;
  const minY = -SVGH / 2 + 6;
  const labelData = chord.groups
    .map((d) => {
      const mid = (d.startAngle + d.endAngle) / 2 - Math.PI / 2;
      const span = d.endAngle - d.startAngle;
      return { d, mid, span };
    })
    .filter((l) => l.span >= labelMinAngle && l.d.index !== gapIndex)
    .sort((a, b) => a.mid - b.mid);

  const placed = [];
  labelData.forEach((l) => {
    const tooClose = placed.some((p) => Math.abs(p.mid - l.mid) < labelGap);
    if (tooClose) return;
    placed.push(l);
  });

  labelSel = g
    .append("g")
    .attr("class", "chord-labels")
    .selectAll("text")
    .data(placed)
    .join("text")
    .attr("x", (l) => {
      const x = labelR * Math.cos(l.mid);
      return Math.max(minX, Math.min(maxX, x));
    })
    .attr("y", (l) => {
      const y = labelR * Math.sin(l.mid);
      return Math.max(minY, Math.min(maxY, y));
    })
    .attr("text-anchor", (l) => (Math.cos(l.mid) > 0 ? "start" : "end"))
    .attr("dominant-baseline", "central")
    .style("font-size", "9px")
    .attr("fill", "#444")
    .text((l) => {
      const label = nodes[l.d.index];
      return label.length > 12 ? label.slice(0, 11) + "…" : label;
    });

  // Ribbons (chords)
  g.append("g")
    .attr("fill-opacity", 0.6)
    .selectAll("path")
    .data(chord)
    .join("path")
    .attr("class", "chord-ribbon")
    .attr("d", ribbon)
    .attr("fill", (d) => colorOf(d.target.index))
    .attr("stroke", (d) => colorOf(d.target.index))
    .attr("stroke-width", 0.3)
    .attr("opacity", 0.6)
    .on("mouseover", (ev, d) => {
      const src = nodes[d.source.index],
        tgt = nodes[d.target.index];
      chTip
        .style("display", "block")
        .html(`<b>${src}</b> ↔ <b>${tgt}</b><br>` + `${valField}: ${d3.format(",.0f")(d.source.value)}`);
    })
    .on("mousemove", chMove)
    .on("mouseout", chOut);

  // Legend
  if (legWrap) {
    const legBox = document.createElement("div");
    legBox.style.cssText = "display:flex;gap:1em;align-items:flex-start";

    const mkCol = (title) => {
      const col = document.createElement("div");
      col.style.cssText = "display:flex;flex-direction:column;gap:.2em;min-width:140px";
      const hdr = document.createElement("div");
      hdr.style.cssText = "font-size:.72em;color:#666;text-transform:uppercase;letter-spacing:.05em";
      hdr.textContent = title;
      col.appendChild(hdr);
      return col;
    };

    const sampleCol = mkCol("Samples");
    const groupCol = mkCol(groupBy === "genus" ? "Genera" : "Categories");
    groupCol.style.borderLeft = "1px solid #e5e7eb";
    groupCol.style.paddingLeft = "0.8em";

    function addLegItem(col, name, idx) {
      const div = document.createElement("div");
      div.style.cssText = "display:flex;align-items:center;gap:.4em;line-height:1.7;cursor:pointer";
      div.innerHTML =
        `<span style="display:inline-block;width:12px;height:12px;border-radius:2px;background:${colorOf(
          idx,
        )};flex-shrink:0"></span>` + `${name.length > 28 ? name.slice(0, 27) + "…" : name}`;
      div.addEventListener("mouseover", () => {
        const arcData = chord.groups.find((g) => g.index === idx);
        if (!arcData) return;
        arcSel.filter((a) => a.index === idx).dispatch("mouseover", { detail: arcData });
      });
      div.addEventListener("mouseout", () => {
        arcSel.filter((a) => a.index === idx).dispatch("mouseout", {});
      });
      col.appendChild(div);
    }

    samples.forEach((name, i) => addLegItem(sampleCol, name, i));
    cats.forEach((name, i) => addLegItem(groupCol, name, i + gapIndex + 1));

    legBox.appendChild(sampleCol);
    legBox.appendChild(groupCol);
    legWrap.appendChild(legBox);
  }
}

function _drawExploreRadar(fd) {
  const wrap = document.getElementById("explore-radar-wrap");
  const legWrap = document.getElementById("explore-radar-legend");
  if (!wrap) return;
  wrap.innerHTML = "";
  if (legWrap) legWrap.innerHTML = "";

  const axes = ["TASS Score", "Breadth %", "Mean Depth", "Minhash Score", "# Reads Aligned"];
  const labels = ["TASS Score", "Coverage", "Mean Depth", "Minhash", "Reads"];
  const wrapW = wrap.clientWidth || 280;
  const SIZE = Math.max(240, Math.min(420, Math.floor(wrapW * 0.9)));
  const SVGW = Math.max(wrapW, SIZE);
  const cx = SVGW / 2,
    cy = SIZE / 2,
    R = SIZE * 0.36;

  // Get top organisms by TASS Score
  const orgs = [...new Set(fd.map((r) => r["Detected Organism"] || ""))].filter(Boolean);
  const orgSel = document.getElementById("explore-radar-orgs");
  if (orgSel && !orgSel.options.length) {
    const topOrgs = orgs.slice().sort((a, b) => {
      const ta =
        d3.max(
          fd.filter((r) => r["Detected Organism"] === a),
          (r) => parseFloat(r["TASS Score"]) || 0,
        ) || 0;
      const tb =
        d3.max(
          fd.filter((r) => r["Detected Organism"] === b),
          (r) => parseFloat(r["TASS Score"]) || 0,
        ) || 0;
      return tb - ta;
    });
    topOrgs.forEach((o, i) => {
      const opt = document.createElement("option");
      opt.value = opt.textContent = o;
      opt.selected = i < 5;
      orgSel.appendChild(opt);
    });
    orgSel.addEventListener("change", () => _drawExploreRadar(filteredData()));
  }
  const selected = orgSel ? [...orgSel.selectedOptions].map((o) => o.value) : orgs.slice(0, 5);

  // Compute per-axis max for normalization
  const maxVal = {};
  axes.forEach((a) => {
    maxVal[a] = d3.max(fd, (r) => parseFloat(r[a]) || 0) || 1;
  });

  const pal = ["#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc948", "#b07aa1", "#9c755f"];
  const svg = d3
    .select(wrap)
    .append("svg")
    .attr("width", SVGW)
    .attr("height", SIZE)
    .style("display", "block")
    .style("margin", "0 auto");
  const n = axes.length;
  const angleStep = (2 * Math.PI) / n;

  // Grid rings
  [0.25, 0.5, 0.75, 1].forEach((frac) => {
    const pts = axes.map((_, i) => {
      const a = i * angleStep - Math.PI / 2;
      return [cx + frac * R * Math.cos(a), cy + frac * R * Math.sin(a)];
    });
    svg
      .append("polygon")
      .attr("points", pts.map((p) => p.join(",")).join(" "))
      .attr("fill", "none")
      .attr("stroke", "#ddd")
      .attr("stroke-width", 1);
  });

  // Axis lines + labels
  axes.forEach((ax, i) => {
    const a = i * angleStep - Math.PI / 2;
    const x2 = cx + R * Math.cos(a),
      y2 = cy + R * Math.sin(a);
    svg
      .append("line")
      .attr("x1", cx)
      .attr("y1", cy)
      .attr("x2", x2)
      .attr("y2", y2)
      .attr("stroke", "#ccc")
      .attr("stroke-width", 1);
    const lx = cx + (R + 16) * Math.cos(a),
      ly = cy + (R + 16) * Math.sin(a);
    svg
      .append("text")
      .attr("x", lx)
      .attr("y", ly)
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "central")
      .style("font-size", "9px")
      .attr("fill", "#555")
      .text(labels[i]);
  });

  // Data polygons
  selected.forEach((org, oi) => {
    const rows = fd.filter((r) => r["Detected Organism"] === org);
    if (!rows.length) return;
    const vals = axes.map((ax) => d3.mean(rows, (r) => parseFloat(r[ax]) || 0) || 0);
    const pts = vals.map((v, i) => {
      const norm = Math.min(1, v / maxVal[axes[i]]);
      const a = i * angleStep - Math.PI / 2;
      return [cx + norm * R * Math.cos(a), cy + norm * R * Math.sin(a)];
    });
    const color = pal[oi % pal.length];
    svg
      .append("polygon")
      .attr("points", pts.map((p) => p.join(",")).join(" "))
      .attr("fill", color)
      .attr("fill-opacity", 0.18)
      .attr("stroke", color)
      .attr("stroke-width", 1.5);
    if (legWrap) {
      const div = document.createElement("div");
      div.style.cssText = `display:flex;align-items:center;gap:.4em;line-height:1.5`;
      div.innerHTML = `<span style="display:inline-block;width:12px;height:12px;border-radius:2px;background:${color}"></span>${
        org.length > 30 ? org.slice(0, 29) + "…" : org
      }`;
      legWrap.appendChild(div);
    }
  });
}

function _drawExploreLollipop(fd) {
  const wrap = document.getElementById("explore-lolli-wrap");
  if (!wrap) return;
  wrap.innerHTML = "";
  const metricEl = document.getElementById("explore-lolli-metric");
  const multi = (document.getElementById("explore-lolli-multi") || { checked: true }).checked;
  const sampleSel = document.getElementById("explore-lolli-sample");
  const sampleLbl = document.getElementById("explore-lolli-sample-label");
  const metric = metricEl ? metricEl.value : "TASS Score";

  // Aggregate: max per organism (or per organism+sample if multi)
  const samples = _orderedSamples([...new Set(fd.map((r) => r["Specimen ID"] || ""))].filter(Boolean));
  const orgs = [...new Set(fd.map((r) => r["Detected Organism"] || ""))].filter(Boolean);
  const pal = d3
    .scaleOrdinal(["#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc948", "#b07aa1", "#9c755f"])
    .domain(samples);

  if (sampleSel && sampleLbl) {
    sampleLbl.classList.toggle("hidden", multi || samples.length <= 1);
    sampleSel.classList.toggle("hidden", multi || samples.length <= 1);
    if (!sampleSel._wired) {
      sampleSel._wired = true;
      sampleSel.addEventListener("change", () => drawExplore());
    }
    if (!sampleSel.options.length || sampleSel.options.length !== samples.length) {
      const prev = sampleSel.value;
      sampleSel.innerHTML = "";
      samples.forEach((s) => {
        const opt = document.createElement("option");
        opt.value = s;
        opt.textContent = s;
        sampleSel.appendChild(opt);
      });
      if (prev && samples.includes(prev)) sampleSel.value = prev;
    }
  }

  let fdUse = fd;
  let singleSample = "";
  if (!multi && samples.length > 0) {
    singleSample = sampleSel && sampleSel.value ? sampleSel.value : samples[0];
    if (sampleSel && sampleSel.value !== singleSample) sampleSel.value = singleSample;
    fdUse = fd.filter((r) => (r["Specimen ID"] || "") === singleSample);
  }

  // Build data: for each organism, max value across samples
  const orgVals = orgs
    .map((o) => {
      const rows = fdUse.filter((r) => r["Detected Organism"] === o);
      const val = d3.max(rows, (r) => parseFloat(r[metric]) || 0) || 0;
      return { org: o, val, rows };
    })
    .sort((a, b) => b.val - a.val)
    .slice(0, 25);

  const mL = Math.max(120, orgs.reduce((m, o) => Math.max(m, o.length), 0) * 5.5 + 10);
  const mT = 10,
    mR = 40,
    mB = 30,
    rowH = 22;
  const H = orgVals.length * rowH + mT + mB;
  const W = wrap.clientWidth || 640;
  const iW = W - mL - mR;
  const xMax = d3.max(orgVals, (d) => d.val) || 1;
  const xS = d3.scaleLinear().domain([0, xMax]).range([0, iW]).nice();
  const yS = d3
    .scaleBand()
    .domain(orgVals.map((d) => d.org))
    .range([0, H - mT - mB])
    .padding(0.3);

  const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H);
  const g = svg.append("g").attr("transform", `translate(${mL},${mT})`);

  g.append("g").call(d3.axisLeft(yS).tickSize(0)).select(".domain").remove();
  g.selectAll("text").style("font-size", "10px");
  g.append("g")
    .attr("transform", `translate(0,${H - mT - mB})`)
    .call(d3.axisBottom(xS).ticks(5));
  svg
    .append("text")
    .attr("x", mL + iW / 2)
    .attr("y", H - 4)
    .attr("text-anchor", "middle")
    .style("font-size", "10px")
    .attr("fill", "#666")
    .text(metric);

  // Shared rich tooltip anchored to body
  const lTip = d3
    .select("body")
    .selectAll(".explore-lolli-tip")
    .data([1])
    .join("div")
    .attr("class", "explore-lolli-tip")
    .style("position", "fixed")
    .style("background", "rgba(20,20,30,.88)")
    .style("color", "#fff")
    .style("padding", ".4em .7em")
    .style("border-radius", "5px")
    .style("font-size", "11.5px")
    .style("pointer-events", "none")
    .style("display", "none")
    .style("max-width", "260px")
    .style("line-height", "1.65")
    .style("z-index", "9999")
    .style("box-shadow", "0 2px 8px rgba(0,0,0,.35)");

  const _tipHtml = (r, val) =>
    `<b style="font-size:12px">${r["Detected Organism"] || "—"}</b><br>` +
    `<span style="color:#adf">Sample:</span> ${r["Specimen ID"] || "—"}<br>` +
    `<span style="color:#adf">${metric}:</span> ${(+val || 0).toFixed(2)}<br>` +
    `<span style="color:#adf"># Reads Aligned:</span> ${(+(r["# Reads Aligned"] || 0)).toLocaleString()}<br>` +
    `<span style="color:#adf">K2 Reads:</span> ${(+(r["K2 Reads"] || 0)).toLocaleString()}<br>` +
    `<span style="color:#adf">Category:</span> ${r["Microbial Category"] || "—"}<br>` +
    `<span style="color:#adf">Coverage:</span> ${parseFloat(r["Breadth %"] || r["Coverage"] || 0).toFixed(2)}%<br>` +
    `<span style="color:#adf">Mean Depth:</span> ${parseFloat(r["Mean Depth"] || 0).toFixed(3)}<br>` +
    `<span style="color:#adf">Passes Threshold:</span> ${r["Passes Threshold"] || "—"}`;

  const _onMove = (ev) => lTip.style("left", ev.clientX + 14 + "px").style("top", ev.clientY - 10 + "px");
  const _onOut = () => lTip.style("display", "none");

  if (multi && samples.length > 1) {
    const bw = yS.bandwidth() / samples.length;
    orgVals.forEach((d) => {
      d.rows.forEach((r) => {
        const sp = r["Specimen ID"] || "";
        const si = samples.indexOf(sp);
        const val = parseFloat(r[metric]) || 0;
        const yPos = yS(d.org) + si * bw + bw / 2;
        g.append("line")
          .attr("x1", 0)
          .attr("x2", xS(val))
          .attr("y1", yPos)
          .attr("y2", yPos)
          .attr("stroke", pal(sp))
          .attr("stroke-width", 1.5)
          .attr("opacity", 0.7);
        g.append("circle")
          .attr("cx", xS(val))
          .attr("cy", yPos)
          .attr("r", 4)
          .attr("fill", pal(sp))
          .attr("stroke", "#fff")
          .attr("stroke-width", 0.5)
          .style("cursor", "crosshair")
          .on("mouseover", (ev) => {
            lTip.style("display", "block").html(_tipHtml(r, val));
            _onMove(ev);
          })
          .on("mousemove", _onMove)
          .on("mouseout", _onOut);
      });
    });
  } else {
    orgVals.forEach((d) => {
      const yPos = yS(d.org) + yS.bandwidth() / 2;
      const rep = d.rows[0] || {};
      g.append("line")
        .attr("x1", 0)
        .attr("x2", xS(d.val))
        .attr("y1", yPos)
        .attr("y2", yPos)
        .attr("stroke", "#1565c0")
        .attr("stroke-width", 1.5);
      g.append("circle")
        .attr("cx", xS(d.val))
        .attr("cy", yPos)
        .attr("r", 5)
        .attr("fill", "#1565c0")
        .attr("stroke", "#fff")
        .attr("stroke-width", 1)
        .style("cursor", "crosshair")
        .on("mouseover", (ev) => {
          lTip.style("display", "block").html(_tipHtml(rep, d.val));
          _onMove(ev);
        })
        .on("mousemove", _onMove)
        .on("mouseout", _onOut);
      g.append("text")
        .attr("x", xS(d.val) + 7)
        .attr("y", yPos + 4)
        .style("font-size", "9px")
        .attr("fill", "#333")
        .text(typeof d.val === "number" ? d.val.toFixed(1) : d.val);
    });
  }
}

/* ───────────────────────────────────────────────────────────────────────
         _drawExploreCorrelogram(fd)
           Pearson correlation matrix across user-selected numeric metrics, drawn
           over the same `fd` rows that every other Explore plot uses — so it
           is automatically synced to the sidebar filters, sample renames, and
           hide/show toggles.

           Each cell shows r ∈ [-1, +1] color-encoded with a diverging palette
           (red = negative, white = zero, blue = positive). Diagonal cells are
           marked self-correlation. Only pairs with ≥3 finite observations are
           computed — sparser pairs render as a hatched "n/a" cell.

           Optional log10 transform helps when a metric (e.g. read counts) spans
           many orders of magnitude — we add 1 before log to swallow zeros.
         ─────────────────────────────────────────────────────────────────── */
function _drawExploreCorrelogram(fd) {
  const wrap = document.getElementById("explore-corr-wrap");
  const empty = document.getElementById("explore-corr-empty");
  if (!wrap) return;
  wrap.innerHTML = "";
  if (empty) empty.style.display = "none";

  // Which metric columns to include (multi-select).
  const sel = document.getElementById("explore-corr-metrics");
  let metrics = [];
  if (sel) metrics = Array.from(sel.selectedOptions).map((o) => o.value);
  if (!metrics.length) {
    metrics = ["TASS Score", "# Reads Aligned", "Breadth %", "Mean Depth", "% Reads"];
  }
  const showVals = (document.getElementById("explore-corr-show-vals") || { checked: true }).checked;
  const useLog = (document.getElementById("explore-corr-log") || { checked: false }).checked;

  if (!fd || !fd.length || metrics.length < 2) {
    if (empty) empty.style.display = "block";
    return;
  }

  // Build vectors for each metric — coerce to finite numbers; non-numeric
  // entries are kept as NaN and ignored pairwise. log10 is applied with a
  // +1 offset for "read-like" metrics where 0 is meaningful.
  const N = fd.length;
  const cols = {};
  const looksLogish = (m) =>
    /reads|depth|count/i.test(m) || m === "# Reads Aligned" || m === "K2 Reads" || m === "Mean Depth";
  metrics.forEach((m) => {
    const v = new Array(N);
    for (let i = 0; i < N; i++) {
      const raw = fd[i] ? fd[i][m] : null;
      const n = parseFloat(raw);
      if (!isFinite(n)) {
        v[i] = NaN;
      } else if (useLog && looksLogish(m)) {
        v[i] = Math.log10(Math.max(0, n) + 1);
      } else {
        v[i] = n;
      }
    }
    cols[m] = v;
  });

  // Pearson correlation between two vectors, ignoring index pairs where
  // either side is NaN. Returns null if fewer than 3 valid pairs or if
  // one side has zero variance.
  function _pearson(a, b) {
    let n = 0,
      sa = 0,
      sb = 0,
      saa = 0,
      sbb = 0,
      sab = 0;
    const len = Math.min(a.length, b.length);
    for (let i = 0; i < len; i++) {
      const x = a[i],
        y = b[i];
      if (!isFinite(x) || !isFinite(y)) continue;
      n++;
      sa += x;
      sb += y;
      saa += x * x;
      sbb += y * y;
      sab += x * y;
    }
    if (n < 3) return { r: null, n };
    const meanA = sa / n,
      meanB = sb / n;
    const varA = saa / n - meanA * meanA;
    const varB = sbb / n - meanB * meanB;
    const cov = sab / n - meanA * meanB;
    if (varA <= 0 || varB <= 0) return { r: null, n };
    const r = cov / Math.sqrt(varA * varB);
    // Clamp tiny FP overshoot.
    return { r: Math.max(-1, Math.min(1, r)), n };
  }

  // Compute the full matrix (symmetric; diagonal = 1).
  const M = metrics.length;
  const mat = []; // [{i,j,r,n}]
  for (let i = 0; i < M; i++) {
    for (let j = 0; j < M; j++) {
      if (i === j) {
        mat.push({ i, j, r: 1, n: cols[metrics[i]].filter(isFinite).length });
      } else if (j < i) {
        // mirror the upper-triangle entry already computed
        const prior = mat.find((e) => e.i === j && e.j === i);
        mat.push({ i, j, r: prior ? prior.r : null, n: prior ? prior.n : 0 });
      } else {
        const { r, n } = _pearson(cols[metrics[i]], cols[metrics[j]]);
        mat.push({ i, j, r, n });
      }
    }
  }

  // Layout: square grid; size each cell so the whole thing fits the wrap.
  const wrapW = wrap.clientWidth || 700;
  const labelW = 130; // left + top label gutter
  const cellMax = 56;
  const cell = Math.max(28, Math.min(cellMax, Math.floor((wrapW - labelW - 40) / M)));
  const W = labelW + cell * M + 60; // +60 for the colour-key strip
  const H = labelW + cell * M + 30;

  const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H);
  const gMat = svg.append("g").attr("transform", `translate(${labelW},${labelW})`);

  // Diverging colour scale: r=-1 red → r=0 white → r=+1 blue.
  const color = d3.scaleLinear().domain([-1, 0, 1]).range(["#c0392b", "#ffffff", "#1f4e79"]).clamp(true);

  // Tooltip — re-uses the same styling as the bubble tip.
  const tip = d3
    .select("body")
    .selectAll(".explore-corr-tip")
    .data([1])
    .join("div")
    .attr("class", "explore-corr-tip")
    .style("position", "fixed")
    .style("background", "rgba(20,20,30,.88)")
    .style("color", "#fff")
    .style("padding", ".4em .7em")
    .style("border-radius", "5px")
    .style("font-size", "11.5px")
    .style("pointer-events", "none")
    .style("display", "none")
    .style("max-width", "240px")
    .style("line-height", "1.6")
    .style("z-index", "9999")
    .style("box-shadow", "0 2px 8px rgba(0,0,0,.35)");

  // Cells.
  gMat
    .selectAll("rect.cell")
    .data(mat)
    .join("rect")
    .attr("class", "cell")
    .attr("x", (d) => d.j * cell)
    .attr("y", (d) => d.i * cell)
    .attr("width", cell - 1)
    .attr("height", cell - 1)
    .attr("fill", (d) => (d.r == null ? "#f1f1f1" : color(d.r)))
    .attr("stroke", "#fff")
    .attr("stroke-width", 1)
    .style("cursor", "crosshair")
    .on("mouseover", (ev, d) => {
      const a = metrics[d.i],
        b = metrics[d.j];
      const rTxt = d.r == null ? "n/a (insufficient data)" : d.r.toFixed(3);
      tip
        .style("display", "block")
        .html(
          `<b style="font-size:12px">${a}</b><br>` +
            `<span style="color:#adf">vs:</span> ${b}<br>` +
            `<span style="color:#adf">Pearson r:</span> ${rTxt}<br>` +
            `<span style="color:#adf">n:</span> ${d.n}` +
            (useLog && (looksLogish(a) || looksLogish(b))
              ? '<br><span style="color:#9bd">log10 transform applied</span>'
              : ""),
        );
    })
    .on("mousemove", (ev) => tip.style("left", ev.clientX + 14 + "px").style("top", ev.clientY - 10 + "px"))
    .on("mouseout", () => tip.style("display", "none"));

  // Cell text — the r value (or "—" for n/a). Auto colour for legibility.
  if (showVals) {
    gMat
      .selectAll("text.val")
      .data(mat)
      .join("text")
      .attr("class", "val")
      .attr("x", (d) => d.j * cell + cell / 2)
      .attr("y", (d) => d.i * cell + cell / 2)
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "central")
      .style("font-size", Math.max(8, Math.min(12, cell * 0.28)) + "px")
      .style("pointer-events", "none")
      .attr("fill", (d) => (d.r == null ? "#888" : Math.abs(d.r) > 0.55 ? "#fff" : "#222"))
      .text((d) => (d.r == null ? "—" : d.r.toFixed(2)));
  }

  // Hatched overlay for "no data" cells.
  const defs = svg.append("defs");
  defs
    .append("pattern")
    .attr("id", "corr-hatch")
    .attr("patternUnits", "userSpaceOnUse")
    .attr("width", 6)
    .attr("height", 6)
    .append("path")
    .attr("d", "M0,6 l6,-6")
    .attr("stroke", "#bbb")
    .attr("stroke-width", 1);
  gMat
    .selectAll("rect.nan")
    .data(mat.filter((d) => d.r == null))
    .join("rect")
    .attr("class", "nan")
    .attr("x", (d) => d.j * cell)
    .attr("y", (d) => d.i * cell)
    .attr("width", cell - 1)
    .attr("height", cell - 1)
    .attr("fill", "url(#corr-hatch)")
    .attr("pointer-events", "none");

  // Row + column labels.
  svg
    .append("g")
    .attr("transform", `translate(${labelW - 6},${labelW})`)
    .selectAll("text")
    .data(metrics)
    .join("text")
    .attr("y", (_, i) => i * cell + cell / 2)
    .attr("text-anchor", "end")
    .attr("dominant-baseline", "central")
    .style("font-size", "11px")
    .style("fill", "#222")
    .text((m) => m);
  svg
    .append("g")
    .attr("transform", `translate(${labelW},${labelW - 6})`)
    .selectAll("text")
    .data(metrics)
    .join("text")
    .attr("transform", (_, j) => `translate(${j * cell + cell / 2},0) rotate(-40)`)
    .attr("text-anchor", "start")
    .attr("dominant-baseline", "alphabetic")
    .style("font-size", "11px")
    .style("fill", "#222")
    .text((m) => m);

  // Color-key strip on the right (vertical, -1 → +1).
  const keyX = labelW + cell * M + 18;
  const keyY = labelW;
  const keyH = cell * M;
  const keyW = 12;
  const grad = defs
    .append("linearGradient")
    .attr("id", "corr-grad")
    .attr("x1", "0")
    .attr("x2", "0")
    .attr("y1", "0")
    .attr("y2", "1");
  [
    { o: 0, c: color(1) },
    { o: 0.5, c: color(0) },
    { o: 1, c: color(-1) },
  ].forEach((s) => grad.append("stop").attr("offset", s.o).attr("stop-color", s.c));
  svg
    .append("rect")
    .attr("x", keyX)
    .attr("y", keyY)
    .attr("width", keyW)
    .attr("height", keyH)
    .attr("fill", "url(#corr-grad)")
    .attr("stroke", "#bbb");
  const keyScale = d3.scaleLinear().domain([1, -1]).range([0, keyH]);
  const keyAxis = d3.axisRight(keyScale).ticks(5).tickFormat(d3.format(".1f"));
  svg
    .append("g")
    .attr("transform", `translate(${keyX + keyW},${keyY})`)
    .call(keyAxis)
    .selectAll("text")
    .style("font-size", "10px");
  svg
    .append("text")
    .attr("x", keyX + keyW / 2)
    .attr("y", keyY - 6)
    .attr("text-anchor", "middle")
    .style("font-size", "10px")
    .style("fill", "#444")
    .text("Pearson r");
}

/* ── Explore event listeners ─────────────────────────────────────────── */
[
  "explore-color-by",
  "explore-size-by",
  "explore-lolli-metric",
  "explore-chord-group",
  "explore-chord-value",
  "explore-corr-metrics",
].forEach((id) => {
  const el = document.getElementById(id);
  if (el) el.addEventListener("change", () => drawExplore());
});
document.getElementById("explore-lolli-multi")?.addEventListener("change", () => drawExplore());
document.getElementById("explore-lolli-sample")?.addEventListener("change", () => drawExplore());
document.getElementById("explore-corr-show-vals")?.addEventListener("change", () => drawExplore());
document.getElementById("explore-corr-log")?.addEventListener("change", () => drawExplore());

/* ── Map state — declared here so redraw() and init() can reference these
         without hitting the Temporal Dead Zone (let/const vars aren't hoisted) ── */
let _leafletMap = null; // Leaflet map instance
let _markerLayer = null; // L.layerGroup holding all markers
let _selectedSample = null; // currently-selected sample name (single-sample panel)
let _selectedGroup = null; // currently-selected group recs (multi-sample panel)
let _markerObjects = {}; // sampleName → {marker, rec, color}
let _panelSortCol = "TASS Score";
let _panelSortAsc = false;
let _mapDrag = null; // drag-resize state
// Mirrors PALETTE (colorblind-safe) so map markers match sample colors.
const _MAP_PALETTE = [
  "#0072B2",
  "#E69F00",
  "#009E73",
  "#CC79A7",
  "#56B4E9",
  "#D55E00",
  "#117733",
  "#882255",
  "#DDCC77",
  "#AA4499",
];
// Highlight state for "View" button — must live here so init() can call _buildRunMetaTable
let _runmetaHighlightSample = null;
// Category colour map — must live here so init() finishing doesn't leave it in TDZ
const _CAT_COLORS = {
  Primary: "#c62828",
  Commensal: "#1565c0",
  Opportunistic: "#f57c00",
  Potential: "#6a1b9a",
  Unknown: "#555",
};
// Longitudinal analysis state — declared here (before init IIFE) to avoid TDZ
let _longiOrgList = [];
let _longiSelectedOrgs = new Set();
let _longiOrgColors = {};
let _longiHidden = {};
let _longiBuilt = false;
const _LONGI_ORG_PAL = [
  "#1565c0",
  "#c62828",
  "#2e7d32",
  "#f57c00",
  "#6a1b9a",
  "#00838f",
  "#ad1457",
  "#4527a0",
  "#558b2f",
  "#e65100",
  "#0277bd",
  "#d84315",
  "#00695c",
  "#4a148c",
  "#827717",
];

// ── Map panel drag-resize listeners ───────────────────────────────────
document.addEventListener("mousedown", (e) => {
  if (e.target && e.target.id === "map-resize-handle") {
    const panel = document.getElementById("map-panel");
    if (!panel) return;
    _mapDrag = { startX: e.clientX, startPanelW: panel.offsetWidth };
    document.body.style.cursor = "ew-resize";
    document.body.style.userSelect = "none";
    e.preventDefault();
  }
});
document.addEventListener("mousemove", (e) => {
  if (!_mapDrag) return;
  const dx = _mapDrag.startX - e.clientX;
  const minW = 200,
    maxW = window.innerWidth * 0.72;
  const newW = Math.max(minW, Math.min(maxW, _mapDrag.startPanelW + dx));
  const panel = document.getElementById("map-panel");
  if (panel) panel.style.width = newW + "px";
});
document.addEventListener("mouseup", () => {
  if (_mapDrag) {
    _mapDrag = null;
    document.body.style.cursor = "";
    document.body.style.userSelect = "";
    if (_leafletMap) setTimeout(() => _leafletMap.invalidateSize(), 50);
  }
});

/* ── Tab-aware redraw.
         Previously every filter change re-rendered every chart on every
         tab. With 14k organisms × ~20 samples that meant the heatmap, TASS
         chart, sunburst, coverage, table, etc. all redrew even though only
         one is visible — this dominated the wall-clock cost.
         Now: filter changes invalidate the cache, redraw the currently
         visible tab immediately, and mark every other tab "dirty". When the
         user switches tabs, the new tab is rendered only if dirty. */
const _TAB_DIRTY = {
  summary: false,
  heatmap: false,
  tass: false,
  sunburst: false,
  coverage: false,
  proteins: false,
  histogram: false,
  explore: false,
  table: false,
  map: false,
  runmeta: false,
  novelty: false,
};
// Tracks whether each pane has EVER been drawn. The dirty flags above only get
// flipped true by redraw() (inside __ttRunInit); if init throws before reaching
// redraw() — e.g. an environment quirk in buildTable/_computeBslLevels — a lazily
// drawn pane (novelty, embed, …) would stay clean and a click would render nothing.
// This guard guarantees a first click always draws the pane regardless.
const _TAB_RENDERED = {};
