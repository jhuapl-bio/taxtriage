/* ═══════════════════════════════════════════════════════════════════════════
       -  §  GROUP VIEWS  (group × organism heatmap, group network / chord)
       -     Both views are driven entirely by the shared metaGrouping engine
       -     (38_meta_grouping.js) so whatever columns the user unions in the
       -     "Group by" bar become the rows / nodes here. With no grouping
       -     selected both fall back to individual samples, which keeps them
       -     useful on runs with no categorical metadata at all.
       -
       -       _mgGroupProfiles()   — group → {taxon: value} matrices, the shared
       -                              aggregation both views consume.
       -       _mgBuildGroupHeatmap() — group × organism grid (sub-tab "ghm").
       -       _mgBuildNetwork()      — force-directed graph or chord diagram of
       -                              inter-group organism-profile similarity
       -                              (sub-tab "net").
       -     Similarity reuses the cosine / Jaccard definitions from the
       -     Cross-Entry Comparison tab so the two tabs never disagree.
═══════════════════════════════════════════════════════════════════════════ */

/* ── Append a group glyph as a real SVG node ────────────────────────────
   The map builds its markers from HTML strings (Leaflet divIcon), but D3
   selections need actual elements. This mirrors _mgShapeEl's geometry so a
   group's shape is identical in the map, the heatmap rows, and the network. */
function _mgAppendShape(sel, shape, cx, cy, r, fill, opts) {
  opts = opts || {};
  const pts = (n, rot) => {
    const a = [];
    for (let i = 0; i < n; i++) {
      const t = rot + (i * 2 * Math.PI) / n;
      a.push(`${(cx + r * Math.cos(t)).toFixed(2)},${(cy + r * Math.sin(t)).toFixed(2)}`);
    }
    return a.join(" ");
  };
  let node;
  switch (shape) {
    case "square":
      node = sel
        .append("rect")
        .attr("x", cx - r * 0.88)
        .attr("y", cy - r * 0.88)
        .attr("width", r * 1.76)
        .attr("height", r * 1.76);
      break;
    case "triangle":
      node = sel.append("polygon").attr("points", pts(3, -Math.PI / 2));
      break;
    case "diamond":
      node = sel.append("polygon").attr("points", pts(4, -Math.PI / 2));
      break;
    case "hexagon":
      node = sel.append("polygon").attr("points", pts(6, -Math.PI / 2));
      break;
    case "star": {
      const a = [];
      for (let i = 0; i < 10; i++) {
        const rr = i % 2 ? r * 0.45 : r;
        const t = -Math.PI / 2 + (i * Math.PI) / 5;
        a.push(`${(cx + rr * Math.cos(t)).toFixed(2)},${(cy + rr * Math.sin(t)).toFixed(2)}`);
      }
      node = sel.append("polygon").attr("points", a.join(" "));
      break;
    }
    case "cross": {
      const w = r * 0.4;
      node = sel
        .append("polygon")
        .attr(
          "points",
          `${cx - w},${cy - r} ${cx + w},${cy - r} ${cx + w},${cy - w} ${cx + r},${cy - w} ` +
            `${cx + r},${cy + w} ${cx + w},${cy + w} ${cx + w},${cy + r} ${cx - w},${cy + r} ` +
            `${cx - w},${cy + w} ${cx - r},${cy + w} ${cx - r},${cy - w} ${cx - w},${cy - w}`,
        );
      break;
    }
    case "wedge":
      node = sel
        .append("polygon")
        .attr("points", `${cx},${cy + r} ${cx + r},${cy - r * 0.6} ${cx - r},${cy - r * 0.6}`);
      break;
    default:
      node = sel.append("circle").attr("cx", cx).attr("cy", cy).attr("r", r);
  }
  return node
    .attr("fill", fill)
    .attr("stroke", opts.stroke || "#37474f")
    .attr("stroke-width", opts.strokeWidth == null ? 1 : opts.strokeWidth)
    .attr("fill-opacity", opts.fillOpacity == null ? 1 : opts.fillOpacity);
}

/* ── Shared aggregation ─────────────────────────────────────────────────
   Returns:
     { units:[{id,label,color,shape,samples:[…]}],
       taxa:[…], matrix:{id:{taxon:value}},
       totals:{id:{detections,reads,samples,tassSum,tassN}} }
   `unit` is "group" when a grouping is active, else "sample".                */
function _mgGroupProfiles(opts) {
  opts = opts || {};
  const level = opts.level || "genus";
  const metric = opts.metric || "detections";
  const mg = window.metaGrouping;
  const grouped = !!(mg && mg.active());

  const taxonOf = (r) =>
    level === "genus"
      ? String(r["Genus Name"] || r["Genus"] || "").trim()
      : String(r["Species Name"] || r["Detected Organism"] || "").trim();

  const valOf = (r) => {
    if (metric === "tass") {
      const t = level === "genus" ? parseFloat(r["Genus TASS"]) : parseFloat(r["Species TASS"]);
      const fb = parseFloat(r["TASS Score"]);
      return Math.max(isNaN(t) ? 0 : t, isNaN(fb) ? 0 : fb);
    }
    if (metric === "coverage") {
      const b = parseFloat(r["Breadth %"]);
      const c = parseFloat(r["Coverage"]);
      return Math.max(isNaN(b) ? 0 : b, isNaN(c) ? 0 : c);
    }
    if (metric === "reads") {
      const n = parseFloat(r["# Reads Aligned"]);
      return isNaN(n) ? 0 : n;
    }
    return 1; // detections / prevalence
  };

  // Max for score-like metrics (a group's value is its best hit, not a sum);
  // sum for count-like metrics.
  const isMax = metric === "tass" || metric === "coverage";

  const matrix = {};
  const totals = {};
  const members = {};
  const taxaScore = new Map();
  // prevalence needs the distinct sample count per (unit, taxon)
  const prevalence = {};

  const rows = typeof filteredData === "function" ? filteredData() : [];
  rows.forEach((r) => {
    const sample = r["Specimen ID"] || "";
    if (!sample) return;
    if (typeof sampleHidden !== "undefined" && sampleHidden[sample]) return;
    const taxon = taxonOf(r);
    if (!taxon) return;

    // Groups toggled off in the legend drop out here, which is what keeps the
    // heatmap's Top-N organism ranking and the network's similarity edges
    // consistent with what is actually on screen — filtering only at render
    // time would leave the ranking computed over hidden data.
    const ids = grouped ? mg.groupsOf(sample).filter((k) => !mg.isHidden(k)) : [sample];
    if (!ids.length) return;
    const v = valOf(r);

    ids.forEach((id) => {
      if (!matrix[id]) {
        matrix[id] = {};
        totals[id] = { detections: 0, reads: 0, tassSum: 0, tassN: 0, samples: new Set() };
        members[id] = new Set();
        prevalence[id] = {};
      }
      matrix[id][taxon] = matrix[id][taxon] == null ? v : isMax ? Math.max(matrix[id][taxon], v) : matrix[id][taxon] + v;
      (prevalence[id][taxon] = prevalence[id][taxon] || new Set()).add(sample);
      const t = totals[id];
      t.detections += 1;
      const rd = parseFloat(r["# Reads Aligned"]);
      if (!isNaN(rd)) t.reads += rd;
      const ts = parseFloat(r["TASS Score"]);
      if (!isNaN(ts)) {
        t.tassSum += ts;
        t.tassN += 1;
      }
      t.samples.add(sample);
      members[id].add(sample);
      taxaScore.set(taxon, (taxaScore.get(taxon) || 0) + 1);
    });
  });

  // prevalence = fraction of the unit's samples carrying the taxon
  if (metric === "prevalence") {
    Object.keys(matrix).forEach((id) => {
      const n = totals[id].samples.size || 1;
      Object.keys(matrix[id]).forEach((tx) => {
        matrix[id][tx] = (prevalence[id][tx].size / n) * 100;
      });
    });
  }

  const model = window.metaGrouping ? window.metaGrouping.model() : { groups: {}, order: [] };
  const ids = Object.keys(matrix);
  // Group order follows the engine's (largest first); sample order follows the
  // report's existing sample ordering so the view matches the sidebar.
  const order = grouped
    ? model.visible.filter((k) => matrix[k])
    : typeof _orderedSamples === "function"
      ? _orderedSamples(ids)
      : ids.sort();
  // Anything the engine didn't know about (shouldn't happen, but be safe)
  ids.forEach((id) => {
    if (order.indexOf(id) === -1) order.push(id);
  });

  const units = order.map((id) => ({
    id,
    label: id,
    color: grouped
      ? mg.color(id)
      : (typeof sampleColors !== "undefined" && sampleColors[id]) || "#0072B2",
    shape: grouped ? mg.shape(id) : "circle",
    samples: Array.from(members[id] || []),
    totals: totals[id],
    // Highlight state, so these views fade the same groups the map fades.
    highlighted: grouped ? mg.isHighlighted(id) : false,
    dimmed: grouped ? mg.isDimmed(id) : false,
  }));

  const taxa = Array.from(taxaScore.keys()).sort((a, b) => (taxaScore.get(b) || 0) - (taxaScore.get(a) || 0));

  return { grouped, unit: grouped ? "group" : "sample", units, taxa, matrix, totals };
}

/* ═══════════════ GROUP × ORGANISM HEATMAP ═══════════════════════════════ */

let _mgHeatSort = { col: null, asc: false };

function _mgBuildGroupHeatmap() {
  const wrap = document.getElementById("mg-heat-chart");
  const nodata = document.getElementById("mg-heat-no-data");
  if (!wrap) return;

  const level = (document.getElementById("mg-heat-level") || {}).value || "genus";
  const metric = (document.getElementById("mg-heat-metric") || {}).value || "detections";
  const topN = parseInt((document.getElementById("mg-heat-topn") || {}).value || "25", 10);

  const grouped = !!(window.metaGrouping && window.metaGrouping.active());
  const hint = document.getElementById("mg-heat-hint");
  if (hint) {
    hint.textContent = grouped
      ? `Rows are ${window.metaGrouping.describe()} groups. Click a group below to hide it, shift-click to show only that one.`
      : 'No grouping selected — rows are individual samples. Pick one or more columns in the "Group by" bar to aggregate.';
  }
  // Per-tab toggle row. Same chips as the shared bar, so hiding a group here
  // also updates the map and the network. No onPick: the engine's change
  // broadcast already rebuilds whichever group view is active, and adding one
  // here would render this chart twice per click.
  if (typeof window._mgRenderLegend === "function") window._mgRenderLegend("mg-heat-legend", {});

  const P = _mgGroupProfiles({ level, metric });
  wrap.innerHTML = "";
  if (!P.units.length || !P.taxa.length) {
    // Distinguish "your filters excluded everything" from "you switched every
    // group off" — the fix is completely different and the empty grid alone
    // does not say which happened.
    if (nodata) {
      nodata.innerHTML = _mgEmptyMessage();
      nodata.style.display = "block";
    }
    return;
  }
  if (nodata) nodata.style.display = "none";

  const taxa = P.taxa.slice(0, Math.max(1, topN));

  // Optional column sort: clicking an organism header sorts rows by that
  // organism's value, so a user can rank sites by a single pathogen.
  let units = P.units.slice();
  if (_mgHeatSort.col && taxa.indexOf(_mgHeatSort.col) !== -1) {
    const c = _mgHeatSort.col;
    units.sort((a, b) => {
      const va = P.matrix[a.id][c] || 0;
      const vb = P.matrix[b.id][c] || 0;
      return _mgHeatSort.asc ? va - vb : vb - va;
    });
  }

  let maxV = 0;
  units.forEach((u) => taxa.forEach((t) => (maxV = Math.max(maxV, P.matrix[u.id][t] || 0))));
  if (!maxV) maxV = 1;

  const cellW = 26;
  const cellH = 22;
  const labelW = Math.min(
    300,
    Math.max(120, ...units.map((u) => Math.min(300, u.label.length * 6.6 + 34))),
  );
  const headH = Math.min(220, Math.max(80, ...taxa.map((t) => t.length * 6.2)) + 12);
  const legendH = 46;
  const W = labelW + taxa.length * cellW + 24;
  const H = headH + units.length * cellH + legendH;

  const color = d3.scaleSequential(d3.interpolateYlGnBu).domain([0, maxV]);

  const svg = d3
    .select(wrap)
    .append("svg")
    .attr("width", W)
    .attr("height", H)
    .attr("font-family", "inherit");

  // ── column headers (rotated organism names, clickable to sort) ──────────
  const headSel = svg
    .append("g")
    .selectAll("text")
    .data(taxa)
    .join("text")
    .attr("transform", (d, i) => `translate(${labelW + i * cellW + cellW / 2},${headH - 6}) rotate(-60)`)
    .attr("font-size", 10)
    .attr("fill", (d) => (_mgHeatSort.col === d ? "#0072B2" : "#455a64"))
    .attr("font-weight", (d) => (_mgHeatSort.col === d ? 700 : 400))
    .style("cursor", "pointer")
    .text((d) => (d.length > 34 ? d.slice(0, 33) + "…" : d));
  headSel.append("title").text((d) => d + " — click to sort rows by this organism");
  // Toggle direction when re-clicking the same column, else sort descending.
  headSel.on("click", (ev, d) => {
    _mgHeatSort = { col: d, asc: _mgHeatSort.col === d ? !_mgHeatSort.asc : false };
    _mgBuildGroupHeatmap();
  });

  // ── rows ────────────────────────────────────────────────────────────────
  const rowG = svg
    .append("g")
    .selectAll("g")
    .data(units)
    .join("g")
    .attr("transform", (d, i) => `translate(0,${headH + i * cellH})`)
    // A highlighted group keeps full contrast; the rest fade, matching the map.
    .attr("opacity", (d) => (d.dimmed ? 0.32 : 1));

  // Group glyph — the same shape the map marker and network node use.
  rowG.each(function (u) {
    _mgAppendShape(d3.select(this), u.shape, 9, cellH / 2, 6.5, u.color);
  });

  rowG
    .append("text")
    .attr("x", 22)
    .attr("y", cellH / 2 + 4)
    .attr("font-size", 11)
    .attr("fill", "#263238")
    .style("cursor", "pointer")
    .text((d) => (d.label.length > 40 ? d.label.slice(0, 39) + "…" : d.label))
    .on("click", (ev, d) => _mgFocusUnit(d))
    .append("title")
    .text((d) => `${d.label}\n${d.samples.length} sample(s): ${d.samples.slice(0, 20).join(", ")}`);

  rowG
    .selectAll("rect.cell")
    .data((u) => taxa.map((t) => ({ u, t, v: P.matrix[u.id][t] || 0 })))
    .join("rect")
    .attr("class", "cell")
    .attr("x", (d, i) => labelW + i * cellW)
    .attr("y", 1)
    .attr("width", cellW - 1.5)
    .attr("height", cellH - 2)
    .attr("rx", 2)
    .attr("fill", (d) => (d.v ? color(d.v) : "#f2f5f8"))
    .attr("stroke", "#fff")
    .attr("stroke-width", 0.5)
    .style("cursor", "pointer")
    .append("title")
    .text((d) => `${d.u.label}\n${d.t}\n${_mgMetricLabel(metric)}: ${_mgFmt(d.v, metric)}`);

  // ── colour legend ───────────────────────────────────────────────────────
  const lg = svg.append("g").attr("transform", `translate(${labelW},${headH + units.length * cellH + 16})`);
  const gradId = "mg-heat-grad";
  const defs = svg.append("defs");
  const grad = defs.append("linearGradient").attr("id", gradId);
  d3.range(0, 1.01, 0.1).forEach((t) =>
    grad.append("stop").attr("offset", `${t * 100}%`).attr("stop-color", color(t * maxV)),
  );
  lg.append("rect").attr("width", 180).attr("height", 10).attr("rx", 2).attr("fill", `url(#${gradId})`);
  lg.append("text").attr("x", 0).attr("y", 24).attr("font-size", 10).attr("fill", "#607d8b").text("0");
  lg.append("text")
    .attr("x", 180)
    .attr("y", 24)
    .attr("text-anchor", "end")
    .attr("font-size", 10)
    .attr("fill", "#607d8b")
    .text(_mgFmt(maxV, metric));
  lg.append("text")
    .attr("x", 196)
    .attr("y", 9)
    .attr("font-size", 10.5)
    .attr("fill", "#455a64")
    .text(_mgMetricLabel(metric));
}

/* Empty-state copy. When the view is blank because every group is toggled off,
   say so and offer the way back — otherwise the user is left widening filters
   that were never the problem. */
function _mgEmptyMessage(needTwo) {
  const mg = window.metaGrouping;
  const hidden = mg && mg.active() ? mg.hiddenCount() : 0;
  if (hidden) {
    const total = mg.model().order.length;
    return (
      '<i class="fas fa-eye-slash" style="margin-right:.35em"></i>' +
      `${hidden} of ${total} groups are switched off. ` +
      'Use <b>Show all</b> in the toggle row above to bring them back.'
    );
  }
  return (
    '<i class="fas fa-info-circle" style="margin-right:.3em"></i>' +
    (needTwo
      ? "Need at least two groups (or two samples) with organism hits. Widen the filters or pick a grouping with more than one value."
      : "No organism hits in the current filter. Widen the sidebar filters or pick a different grouping.")
  );
}

function _mgMetricLabel(m) {
  return (
    {
      detections: "Detections",
      reads: "# Reads Aligned",
      tass: "TASS Score (max)",
      coverage: "Breadth % (max)",
      prevalence: "Prevalence (% of samples)",
    }[m] || m
  );
}

function _mgFmt(v, metric) {
  if (v == null || isNaN(v)) return "—";
  if (metric === "prevalence" || metric === "coverage") return v.toFixed(1) + "%";
  if (metric === "tass") return v.toFixed(1);
  if (v >= 1e6) return (v / 1e6).toFixed(1) + "M";
  if (v >= 1e3) return (v / 1e3).toFixed(1) + "k";
  return String(Math.round(v));
}

/* Clicking a group's row label / node highlights it across every grouped view
   — the same state the legend's first click produces, so the map picks it up
   and draws that group's similarity network. Clicking it again clears the
   highlight rather than advancing to "hidden": accidentally hiding a group by
   double-clicking its label would be a nasty surprise, and the legend chip is
   the deliberate place to reach the hidden state. */
function _mgFocusUnit(u) {
  if (!u) return;
  const mg = window.metaGrouping;
  if (!mg || !mg.active()) return;
  try {
    mg.setState(u.id, mg.isHighlighted(u.id) ? mg.STATE_NORMAL : mg.STATE_HIGHLIGHT);
    if (typeof _mgMapFitVisible === "function") _mgMapFitVisible();
  } catch (e) {}
}

/* ═══════════════ GROUP NETWORK / CHORD ══════════════════════════════════ */

let _mgNetSim = null;

function _mgNetSimilarity(a, b, mode, taxa) {
  // a, b: {taxon: value}
  if (mode === "jaccard") {
    let inter = 0;
    let uni = 0;
    taxa.forEach((t) => {
      const x = (a[t] || 0) > 0;
      const y = (b[t] || 0) > 0;
      if (x || y) uni += 1;
      if (x && y) inter += 1;
    });
    return uni ? inter / uni : 0;
  }
  if (mode === "shared") {
    let inter = 0;
    taxa.forEach((t) => {
      if ((a[t] || 0) > 0 && (b[t] || 0) > 0) inter += 1;
    });
    return inter;
  }
  // cosine on log1p-scaled values so a single huge read count can't dominate
  let dot = 0;
  let na = 0;
  let nb = 0;
  taxa.forEach((t) => {
    const x = Math.log1p(a[t] || 0);
    const y = Math.log1p(b[t] || 0);
    dot += x * y;
    na += x * x;
    nb += y * y;
  });
  return na && nb ? dot / Math.sqrt(na * nb) : 0;
}

function _mgBuildNetwork() {
  const wrap = document.getElementById("mg-net-chart");
  const nodata = document.getElementById("mg-net-no-data");
  if (!wrap) return;

  const layout = (document.getElementById("mg-net-layout") || {}).value || "force";
  const level = (document.getElementById("mg-net-level") || {}).value || "genus";
  const basis = (document.getElementById("mg-net-basis") || {}).value || "presence";
  const simMode = (document.getElementById("mg-net-sim") || {}).value || "cosine";
  const thresh = parseFloat((document.getElementById("mg-net-threshold") || {}).value || "0.3");
  const threshLbl = document.getElementById("mg-net-threshold-val");
  if (threshLbl) threshLbl.textContent = simMode === "shared" ? `≥ ${Math.round(thresh * 20)} taxa` : thresh.toFixed(2);

  if (typeof window._mgRenderLegend === "function") window._mgRenderLegend("mg-net-legend", {});

  const metric = basis === "tass" ? "tass" : basis === "coverage" ? "coverage" : "reads";
  const P = _mgGroupProfiles({ level, metric });

  if (_mgNetSim) {
    _mgNetSim.stop();
    _mgNetSim = null;
  }
  wrap.innerHTML = "";

  if (P.units.length < 2) {
    if (nodata) {
      nodata.innerHTML = _mgEmptyMessage(true);
      nodata.style.display = "block";
    }
    return;
  }
  if (nodata) nodata.style.display = "none";

  const hint = document.getElementById("mg-net-hint");
  if (hint) {
    hint.textContent =
      (window.metaGrouping && window.metaGrouping.active()
        ? `Nodes are ${window.metaGrouping.describe()} groups — click a group below to hide it, shift-click to isolate it`
        : "Nodes are individual samples") +
      `. An edge means the two share an organism profile above the link threshold.`;
  }

  // Cap the node count — a force layout of 400 samples is unreadable and slow.
  const MAXN = 80;
  const units = P.units.slice(0, MAXN);
  const truncated = P.units.length > MAXN;
  const taxa = P.taxa;

  // Pairwise similarity
  const links = [];
  let maxSim = 0;
  for (let i = 0; i < units.length; i++) {
    for (let j = i + 1; j < units.length; j++) {
      const s = _mgNetSimilarity(P.matrix[units[i].id], P.matrix[units[j].id], simMode, taxa);
      if (s > maxSim) maxSim = s;
      links.push({ source: units[i].id, target: units[j].id, value: s });
    }
  }
  const cut = simMode === "shared" ? Math.round(thresh * 20) : thresh;
  const kept = links.filter((l) => l.value >= cut && l.value > 0);

  const W = Math.max(560, wrap.clientWidth || 720);
  const H = 520;
  const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H).attr("font-family", "inherit");

  const stats = document.getElementById("mg-net-stats");
  if (stats) {
    stats.textContent =
      `${units.length} node${units.length === 1 ? "" : "s"}, ${kept.length} link${kept.length === 1 ? "" : "s"} ` +
      `of ${links.length} possible` +
      (truncated ? ` — showing the ${MAXN} largest (group your samples to see all of them)` : "");
  }

  if (layout === "chord") {
    _mgDrawChord(svg, units, links, W, H, simMode);
  } else {
    _mgDrawForce(svg, units, kept, W, H, simMode, P);
  }
}

function _mgDrawForce(svg, units, links, W, H, simMode, P) {
  const byId = new Map(units.map((u) => [u.id, u]));
  const nodes = units.map((u) => ({
    id: u.id,
    label: u.label,
    color: u.color,
    shape: u.shape,
    dimmed: !!u.dimmed,
    n: (u.totals && u.totals.samples ? u.totals.samples.size : u.samples.length) || 1,
    det: (u.totals && u.totals.detections) || 0,
  }));
  const dimOf = (id) => {
    const n = nodes.find((x) => x.id === id);
    return n ? n.dimmed : false;
  };
  const maxN = d3.max(nodes, (d) => d.n) || 1;
  const rOf = (d) => 7 + 13 * Math.sqrt(d.n / maxN);
  const maxV = d3.max(links, (l) => l.value) || 1;

  const linkSel = svg
    .append("g")
    .attr("stroke", "#90a4ae")
    .selectAll("line")
    .data(links.map((l) => Object.assign({}, l)))
    .join("line")
    .attr("stroke-width", (d) => 0.6 + 3.4 * (d.value / maxV))
    // An edge fades if either endpoint is faded, so a highlighted group's
    // links stay readable against a quiet background.
    .attr("stroke-opacity", (d) => {
      const base = 0.2 + 0.5 * (d.value / maxV);
      const a = typeof d.source === "object" ? d.source.id : d.source;
      const b = typeof d.target === "object" ? d.target.id : d.target;
      return dimOf(a) || dimOf(b) ? base * 0.25 : base;
    });
  linkSel.append("title").text((d) => `${d.source} ↔ ${d.target}\n${simMode}: ${d.value.toFixed(3)}`);

  const nodeSel = svg
    .append("g")
    .selectAll("g")
    .data(nodes)
    .join("g")
    .style("cursor", "pointer")
    .attr("opacity", (d) => (d.dimmed ? 0.3 : 1));

  nodeSel.each(function (d) {
    _mgAppendShape(d3.select(this), d.shape, 0, 0, rOf(d), d.color, {
      strokeWidth: 1.4,
      fillOpacity: 0.88,
    });
  });

  nodeSel
    .append("title")
    .text((d) => `${d.label}\n${d.n} sample(s), ${d.det} detection(s)`);

  nodeSel
    .append("text")
    .attr("y", (d) => rOf(d) + 11)
    .attr("text-anchor", "middle")
    .attr("font-size", 10)
    .attr("fill", "#37474f")
    .attr("pointer-events", "none")
    .text((d) => (d.label.length > 22 ? d.label.slice(0, 21) + "…" : d.label));

  nodeSel.on("click", (ev, d) => _mgFocusUnit(byId.get(d.id)));

  _mgNetSim = d3
    .forceSimulation(nodes)
    .force(
      "link",
      d3
        .forceLink(linkSel.data())
        .id((d) => d.id)
        .distance((d) => 190 - 130 * (d.value / maxV))
        .strength((d) => 0.08 + 0.5 * (d.value / maxV)),
    )
    .force("charge", d3.forceManyBody().strength(-260))
    .force("center", d3.forceCenter(W / 2, H / 2))
    .force(
      "collide",
      d3.forceCollide().radius((d) => rOf(d) + 16),
    )
    .on("tick", () => {
      nodes.forEach((d) => {
        d.x = Math.max(28, Math.min(W - 28, d.x));
        d.y = Math.max(24, Math.min(H - 24, d.y));
      });
      linkSel
        .attr("x1", (d) => d.source.x)
        .attr("y1", (d) => d.source.y)
        .attr("x2", (d) => d.target.x)
        .attr("y2", (d) => d.target.y);
      nodeSel.attr("transform", (d) => `translate(${d.x},${d.y})`);
    });

  nodeSel.call(
    d3
      .drag()
      .on("start", (ev, d) => {
        if (!ev.active && _mgNetSim) _mgNetSim.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
      })
      .on("drag", (ev, d) => {
        d.fx = ev.x;
        d.fy = ev.y;
      })
      .on("end", (ev, d) => {
        if (!ev.active && _mgNetSim) _mgNetSim.alphaTarget(0);
        d.fx = null;
        d.fy = null;
      }),
  );
}

function _mgDrawChord(svg, units, links, W, H, simMode) {
  // Chord needs a dense square matrix; cap harder than the force layout since
  // every pair gets a ribbon.
  const MAXC = 24;
  const us = units.slice(0, MAXC);
  const idx = new Map(us.map((u, i) => [u.id, i]));
  const n = us.length;
  const mat = Array.from({ length: n }, () => new Array(n).fill(0));
  links.forEach((l) => {
    const a = idx.get(l.source);
    const b = idx.get(l.target);
    if (a == null || b == null) return;
    mat[a][b] = l.value;
    mat[b][a] = l.value;
  });

  const outer = Math.min(W, H) / 2 - 110;
  const inner = outer - 14;
  const g = svg.append("g").attr("transform", `translate(${W / 2},${H / 2})`);

  const chords = d3.chord().padAngle(0.045).sortSubgroups(d3.descending)(mat);

  const arc = d3.arc().innerRadius(inner).outerRadius(outer);
  const ribbon = d3.ribbon().radius(inner);

  g.append("g")
    .selectAll("path")
    .data(chords.groups)
    .join("path")
    .attr("d", arc)
    .attr("fill", (d) => us[d.index].color)
    .attr("stroke", "#fff")
    .append("title")
    .text((d) => us[d.index].label);

  g.append("g")
    .attr("fill-opacity", 0.62)
    .selectAll("path")
    .data(chords)
    .join("path")
    .attr("d", ribbon)
    .attr("fill", (d) => us[d.source.index].color)
    .attr("stroke", "#ffffff")
    .attr("stroke-opacity", 0.5)
    .append("title")
    .text(
      (d) =>
        `${us[d.source.index].label} ↔ ${us[d.target.index].label}\n${simMode}: ${d.source.value.toFixed(3)}`,
    );

  g.append("g")
    .selectAll("text")
    .data(chords.groups)
    .join("text")
    .attr("transform", (d) => {
      const a = (d.startAngle + d.endAngle) / 2;
      return `rotate(${(a * 180) / Math.PI - 90}) translate(${outer + 8}) ${a > Math.PI ? "rotate(180)" : ""}`;
    })
    .attr("text-anchor", (d) => ((d.startAngle + d.endAngle) / 2 > Math.PI ? "end" : "start"))
    .attr("dy", "0.35em")
    .attr("font-size", 10.5)
    .attr("fill", "#37474f")
    .text((d) => {
      const l = us[d.index].label;
      return l.length > 24 ? l.slice(0, 23) + "…" : l;
    });

  if (units.length > MAXC) {
    svg
      .append("text")
      .attr("x", 10)
      .attr("y", H - 10)
      .attr("font-size", 11)
      .attr("fill", "#78909c")
      .text(`Showing the ${MAXC} largest of ${units.length} — chord ribbons get unreadable past that.`);
  }
}

/* ── Wiring ──────────────────────────────────────────────────────────────── */
(function () {
  function wire(ids, fn) {
    ids.forEach((id) => {
      const el = document.getElementById(id);
      if (el) el.addEventListener(el.tagName === "INPUT" && el.type === "range" ? "input" : "change", fn);
    });
  }
  function ready() {
    wire(["mg-heat-level", "mg-heat-metric", "mg-heat-topn"], () => {
      _mgHeatSort = { col: null, asc: false };
      _mgBuildGroupHeatmap();
    });
    wire(
      ["mg-net-layout", "mg-net-level", "mg-net-basis", "mg-net-sim", "mg-net-threshold"],
      () => _mgBuildNetwork(),
    );
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", ready);
  else ready();
})();
