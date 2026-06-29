/* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB: SUNBURST           (data-tab="sunburst")
       -     D3 v7 half-dome partition chart, MBON Dashboard style. Renders
       -     ONE panel per visible sample inside an IIFE that registers panels
       -     in a private `_panels` map keyed by id. Zoom path + search are
       -     synced across all panels via shared globals.
       -     Per-panel state: hroot, focusedNode, colorScale, legendEl, etc.
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  // ── Per-panel state registry ──────────────────────────────────────────
  const _panels = {};
  let _nextId = 1;

  const _sun20 = [
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

  // ── Global sync state (zoom path + search shared across all panels) ───
  let _globalZoomPath = []; // node name path from root, e.g. ['Bacteria','Proteobacteria']
  let _globalSearch = "";

  function _getZoomPath(node) {
    if (!node || node.data.name === "root") return [];
    return node
      .ancestors()
      .reverse()
      .slice(1)
      .map((n) => n.data.name);
  }

  // Walk hroot following path; stop at deepest ancestor found in this hierarchy
  function _findNodeByPath(hroot, path) {
    let cur = hroot;
    for (const name of path) {
      if (!cur.children) return cur;
      const child = cur.children.find((c) => c.data.name === name);
      if (!child) return cur;
      cur = child;
    }
    return cur;
  }

  // ── Color: zoom-aware — distinct color per direct child of focused node ──
  function _getNodeColor(st, d) {
    if (!st.colorScale) return "#999";
    // When zoomed into a non-root node, color by which direct-child "bucket" d belongs to
    const focused = st.focusedNode;
    if (focused && focused !== st.hroot && st._zoomColorScale) {
      let n = d;
      while (n.parent && n.parent !== focused) n = n.parent;
      if (n.parent === focused) return st._zoomColorScale(n.data.name);
    }
    // Default: color by Superkingdom (depth-1 ancestor)
    let n = d;
    while (n.depth > 1 && n.parent) n = n.parent;
    return st.colorScale(n.data.name);
  }

  // Fill-opacity given current focus
  function _opac(d, focused, hroot) {
    if (!focused || focused === hroot) return Math.max(0.3, 1 - d.depth * 0.12);
    let n = d;
    while (n) {
      if (n === focused) return d === focused ? 0 : Math.max(0.35, 1 - (d.depth - focused.depth) * 0.15);
      n = n.parent;
    }
    return 0.06;
  }

  // ── Arc label path (half-dome): flip right-side arcs so text reads L→R ──
  function _arcLabelPath(st, d) {
    const r = (st.ySc(d.y0) + st.ySc(d.y1)) / 2;
    const a0 = Math.max(-Math.PI / 2, Math.min(Math.PI / 2, st.xSc(d.x0)));
    const a1 = Math.max(-Math.PI / 2, Math.min(Math.PI / 2, st.xSc(d.x1)));
    const aMid = (a0 + a1) / 2;
    // Right side of dome (aMid > 0): reverse so text isn't written downward
    const flip = aMid > 0;
    const sa = flip ? a1 : a0;
    const ea = flip ? a0 : a1;
    const sx = r * Math.sin(sa),
      sy = -r * Math.cos(sa);
    const ex = r * Math.sin(ea),
      ey = -r * Math.cos(ea);
    const largeArc = Math.abs(ea - sa) > Math.PI ? 1 : 0;
    return `M${sx.toFixed(2)},${sy.toFixed(2)} A${r.toFixed(2)},${r.toFixed(2)} 0 ${largeArc} ${
      flip ? 0 : 1
    } ${ex.toFixed(2)},${ey.toFixed(2)}`;
  }

  // Arc label visible when arc-length > ~45px
  function _arcLabelVisible(st, d) {
    const r = (st.ySc(d.y0) + st.ySc(d.y1)) / 2;
    const a0 = Math.max(-Math.PI / 2, Math.min(Math.PI / 2, st.xSc(d.x0)));
    const a1 = Math.max(-Math.PI / 2, Math.min(Math.PI / 2, st.xSc(d.x1)));
    return r * Math.abs(a1 - a0) > 45;
  }

  // Refresh label paths + visibility after a zoom transition completes
  function _updateArcLabels(st) {
    if (!st.svg) return;
    st.svg
      .select("defs")
      .selectAll("path.arc-lp")
      .attr("d", (d) => _arcLabelPath(st, d));
    if (st.g) st.g.selectAll("text.arc-lbl").style("opacity", (d) => (_arcLabelVisible(st, d) ? 1 : 0));
  }

  // ── Populate per-panel sample dropdown ───────────────────────────────
  function _populateSampleSel(st) {
    if (!st.sampleSel) return;
    const samples = _orderedSamples(uniq(filteredData().map((r) => r["Specimen ID"] || "")).filter(Boolean));
    const prev = st.sampleSel.value || st._defaultSample || "";
    st._defaultSample = null; // consume one-time default
    st.sampleSel.innerHTML = '<option value="">All Samples</option>';
    samples.forEach((s) => {
      const opt = document.createElement("option");
      opt.value = s;
      opt.textContent = s;
      st.sampleSel.appendChild(opt);
    });
    if (samples.includes(prev)) st.sampleSel.value = prev;
  }

  // ── Build D3 hierarchy ────────────────────────────────────────────────
  function _buildHier(metric, searchTerm, sampleFilter) {
    let fd = filteredData();
    if (sampleFilter) fd = fd.filter((r) => r["Specimen ID"] === sampleFilter);
    if (!fd.length) return null;
    const q = searchTerm ? searchTerm.trim().toLowerCase() : "";
    const TAX_PATH = ["Domain", "Phylum", "Class", "Genus"];
    const rootNode = { name: "root", children: {}, value: 0 };
    fd.forEach((r) => {
      const org = r["Detected Organism"] || "Unknown";
      if (q) {
        const hay = [...TAX_PATH.map((k) => (r[k] || "").toLowerCase()), org.toLowerCase()].join(" ");
        if (!hay.includes(q)) return;
      }
      const v = Math.max(num(r[metric]), 0);
      let node = rootNode.children;
      TAX_PATH.forEach((rank) => {
        const key = r[rank] || `Unknown ${rank}`;
        if (!node[key]) node[key] = { name: key, children: {}, value: 0 };
        node[key].value += v;
        node = node[key].children;
      });
      if (!node[org]) node[org] = { name: org, children: {}, value: 0 };
      node[org].value += v;
    });
    function toH(n) {
      const ch = Object.values(n.children).map(toH);
      return ch.length ? { name: n.name, children: ch } : { name: n.name, value: Math.max(n.value, 0.01) };
    }
    const hier = { name: "root", children: Object.values(rootNode.children).map(toH) };
    return hier.children.length ? hier : null;
  }

  // ── Scrollable children list ──────────────────────────────────────────
  function _updateList(st, node) {
    const el = st.listEl;
    if (!el) return;
    el.innerHTML = "";
    const children = (node ? node.children : st.hroot && st.hroot.children) || [];
    const nm = !node || node === st.hroot ? "Root" : node.data.name;
    const hdr = document.createElement("div");
    hdr.className = "sun-list-hdr";
    hdr.textContent = children.length
      ? `${nm} · ${children.length} item${children.length !== 1 ? "s" : ""}`
      : nm + " (leaf)";
    el.appendChild(hdr);
    if (!children.length) {
      const leafRow = document.createElement("div");
      leafRow.className = "sun-list-row";
      leafRow.style.cssText = "font-style:italic;color:#999";
      const isInt = st.metric === "# Reads Aligned" || st.metric === "K2 Reads";
      leafRow.textContent = isInt ? Math.round(node.value).toLocaleString() + " reads" : (node.value || 0).toFixed(2);
      el.appendChild(leafRow);
      return;
    }
    const total = node ? node.value || 1 : 1;
    const isInt = st.metric === "# Reads Aligned" || st.metric === "K2 Reads";
    const sorted = [...children].sort((a, b) => b.value - a.value);
    sorted.forEach((c) => {
      const color = _getNodeColor(st, c); // zoom-aware: distinct per child of focused
      const pct = ((c.value / total) * 100).toFixed(1);
      const valTxt = isInt ? Math.round(c.value).toLocaleString() : c.value.toFixed(2);
      const row = document.createElement("div");
      row.className = "sun-list-row";
      row.onclick = () => _zoomTo(st, c);
      const dot = document.createElement("span");
      dot.className = "sun-list-dot";
      dot.style.background = color;
      const nm2 = document.createElement("span");
      nm2.className = "sun-list-name";
      nm2.title = c.data.name;
      nm2.textContent = c.data.name + (c.children && c.children.length ? " ›" : "");
      const val = document.createElement("span");
      val.className = "sun-list-val";
      val.textContent = `${valTxt} (${pct}%)`;
      row.append(dot, nm2, val);
      el.appendChild(row);
    });
  }

  // ── Legend — zoom-aware (shows children of focused, or Superkingdoms at root) ──
  function _updateLegend(st) {
    const el = st.legendEl;
    if (!el || !st.hroot) return;
    el.innerHTML = "";
    const focused = st.focusedNode;
    const items = focused && focused !== st.hroot && focused.children ? focused.children : st.hroot.children || [];
    items.forEach((n) => {
      const row = document.createElement("div");
      row.className = "sun-legend-row";
      row.style.cursor = "pointer";
      const sw = document.createElement("span");
      sw.className = "sun-legend-sw";
      sw.style.background = _getNodeColor(st, n);
      const lbl = document.createElement("span");
      lbl.textContent = n.data.name;
      row.append(sw, lbl);
      row.onclick = () => _zoomTo(st, n);
      el.appendChild(row);
    });
  }

  // ── Zoom to a node ────────────────────────────────────────────────────
  function _zoomTo(st, node, skipBroadcast) {
    if (!st.paths || !st.xSc || !st.ySc || !st.arc) return;
    st.focusedNode = node;
    const { xSc, ySc, maxR, arc, hroot, paths } = st;
    const p = node || hroot;
    // Build per-zoom color scale for direct children of focused node
    if (p && p !== hroot && p.children && p.children.length) {
      st._zoomColorScale = d3.scaleOrdinal(_sun20).domain(p.children.map((c) => c.data.name));
    } else {
      st._zoomColorScale = null;
    }
    if (st.bcEl) {
      st.bcEl.textContent =
        !node || node === hroot
          ? "root"
          : node
              .ancestors()
              .reverse()
              .slice(1)
              .map((n) => n.data.name)
              .join(" › ");
    }
    _updateList(st, node === hroot ? hroot : node);
    _updateLegend(st);
    // Hide arc labels immediately — reappear correctly when transition ends
    if (st.g) st.g.selectAll("text.arc-lbl").style("opacity", 0);
    // Apply new colors immediately before transition
    paths.attr("fill", (d) => _getNodeColor(st, d));
    paths
      .transition()
      .duration(750)
      .ease(d3.easeCubicInOut)
      .tween("scale", () => {
        const xd = d3.interpolate(xSc.domain(), [p.x0, p.x1]);
        const yd = d3.interpolate(ySc.domain(), [p.y0, 1]);
        const yr = d3.interpolate(ySc.range(), [p.depth ? maxR * 0.12 : maxR * 0.04, maxR]);
        return (t) => {
          xSc.domain(xd(t));
          ySc.domain(yd(t)).range(yr(t));
        };
      })
      .attrTween("d", function (d) {
        return () => arc(d);
      })
      .attr("fill-opacity", (d) => _opac(d, p, hroot))
      .on("end", () => _updateArcLabels(st));
    // Broadcast zoom to all other panels (skip when already responding to a broadcast)
    if (!skipBroadcast) {
      _globalZoomPath = _getZoomPath(!node || node === hroot ? null : node);
      Object.keys(_panels).forEach((oid) => {
        const ost = _panels[Number(oid)];
        if (!ost || ost === st || !ost.hroot || !ost.paths) return;
        const target = _findNodeByPath(ost.hroot, _globalZoomPath);
        _zoomTo(ost, target, true);
      });
    }
  }

  // ── Render one panel ──────────────────────────────────────────────────
  function _renderPanel(id) {
    const st = _panels[id];
    if (!st) return;
    st.svgWrapEl.innerHTML = "";
    const metric = document.getElementById("sun-metric")?.value || "# Reads Aligned";
    const searchTerm = st.searchEl ? st.searchEl.value : "";
    st.metric = metric;
    _populateSampleSel(st);
    const sampleFilter = st.sampleSel ? st.sampleSel.value : "";
    const hier = _buildHier(metric, searchTerm, sampleFilter);
    if (!hier) {
      st.svgWrapEl.innerHTML = '<p style="color:#999;padding:1em;font-size:.82em">No data for this selection.</p>';
      if (st.listEl) st.listEl.innerHTML = "";
      if (st.badgeEl) st.badgeEl.textContent = "0";
      if (st.legendEl) st.legendEl.innerHTML = "";
      return;
    }
    // Normalized partition [0,1] × [0,1]
    st.hroot = d3
      .hierarchy(hier)
      .sum((d) => d.value || 0)
      .sort((a, b) => b.value - a.value);
    d3.partition().size([1, 1])(st.hroot);
    const visNodes = st.hroot.descendants().filter((d) => d.depth > 0 && d.x1 > d.x0);
    if (st.badgeEl) st.badgeEl.textContent = visNodes.length;
    // Reset zoom color scale on full re-render
    st._zoomColorScale = null;
    // Color by Superkingdom at root level
    st.colorScale = d3.scaleOrdinal(_sun20).domain((st.hroot.children || []).map((n) => n.data.name));
    // Sizing
    const panelW = st.svgWrapEl.closest(".sun-panel")?.clientWidth || 360;
    const W = Math.max(panelW - 16, 260);
    const maxR = W / 2 - 8;
    st.maxR = maxR;
    // Scales — half-dome: flat base, arcs fan upward from center
    const innerR = Math.max(10, maxR * 0.025);
    st.xSc = d3
      .scaleLinear()
      .domain([0, 1])
      .range([-Math.PI / 2, Math.PI / 2]);
    st.ySc = d3.scaleSqrt().domain([0, 1]).range([innerR, maxR]);
    st.arc = d3
      .arc()
      .startAngle((d) => Math.max(-Math.PI / 2, Math.min(Math.PI / 2, st.xSc(d.x0))))
      .endAngle((d) => Math.max(-Math.PI / 2, Math.min(Math.PI / 2, st.xSc(d.x1))))
      .innerRadius((d) => Math.max(0, st.ySc(d.y0)))
      .outerRadius((d) => Math.max(0, st.ySc(d.y1)))
      .padAngle(0.003)
      .padRadius(maxR / 2)
      .cornerRadius(3);
    // SVG height = half-dome radius + label padding below baseline
    const SVGh = maxR + 40;
    const svgSel = d3
      .select(st.svgWrapEl)
      .append("svg")
      .attr("width", W)
      .attr("height", SVGh)
      .style("overflow", "visible");
    st.svg = svgSel;
    // ── Defs: arc label paths ──────────────────────────────────────────
    visNodes.forEach((d, i) => {
      d._lid = `arc-lp-${id}-${i}`;
    });
    const defs = svgSel.append("defs");
    defs
      .selectAll("path.arc-lp")
      .data(visNodes)
      .join("path")
      .attr("class", "arc-lp")
      .attr("id", (d) => d._lid)
      .attr("d", (d) => _arcLabelPath(st, d));
    const g = svgSel.append("g").attr("transform", `translate(${W / 2},${maxR + 8})`);
    st.g = g;
    // ── Arcs ──────────────────────────────────────────────────────────
    st.paths = g
      .selectAll("path")
      .data(visNodes)
      .join("path")
      .attr("fill", (d) => _getNodeColor(st, d))
      .attr("fill-opacity", (d) => Math.max(0.3, 1 - d.depth * 0.12))
      .attr("stroke", "#fff")
      .attr("stroke-width", 0.5)
      .style("cursor", "pointer")
      .attr("d", (d) => st.arc(d))
      .on("mouseover", (ev, d) => {
        const chain = d
          .ancestors()
          .reverse()
          .slice(1)
          .map((n) => n.data.name)
          .join(" › ");
        const isInt = metric === "# Reads Aligned" || metric === "K2 Reads";
        const valStr = isInt ? Math.round(d.value).toLocaleString() : d.value.toFixed(2);
        showTip(`<b>${d.data.name}</b><br><small>${chain}</small><br>${metric}: <b>${valStr}</b>`, ev);
      })
      .on("mousemove", moveTip)
      .on("mouseout", hideTip)
      .on("click", (ev, d) => {
        ev.stopPropagation();
        _zoomTo(st, d);
      });
    // ── Arc labels (curved textPath, visible when arc-length allows) ───
    g.selectAll("text.arc-lbl")
      .data(visNodes)
      .join("text")
      .attr("class", "arc-lbl")
      .style("pointer-events", "none")
      .style("font-size", "10px")
      .style("font-family", "system-ui,sans-serif")
      .style("fill", "#fff")
      .style("opacity", (d) => (_arcLabelVisible(st, d) ? 1 : 0))
      .append("textPath")
      .attr("href", (d) => `#${d._lid}`)
      .attr("startOffset", "50%")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "central")
      .text((d) => {
        const name = d.data.name;
        return name.length > 14 ? name.slice(0, 13) + "\u2026" : name;
      });
    // ── Center circle — click goes UP one level ────────────────────────
    const cg = g.append("g").style("cursor", "pointer");
    cg.append("circle")
      .attr("r", innerR)
      .attr("fill", "#fff")
      .attr("fill-opacity", 0.95)
      .attr("stroke", "#dee2e6")
      .attr("stroke-width", 1.5);
    if (st.centerLabelEl) st.centerLabelEl.remove();
    const lbl = sampleFilter || "All";
    const centerLabel = document.createElement("div");
    centerLabel.className = "sun-center-label";
    centerLabel.style.cssText = "text-align:center;margin-top:4px;line-height:1.1";
    centerLabel.innerHTML =
      `<div style="font-size:10px;font-weight:700;color:#333">${lbl.length > 18 ? lbl.slice(0, 17) + "…" : lbl}</div>` +
      `<div style="font-size:9px;color:#888">${metric.length > 22 ? metric.slice(0, 21) + "…" : metric}</div>`;
    st.svgWrapEl.appendChild(centerLabel);
    st.centerLabelEl = centerLabel;
    cg.on("click", () => {
      const focused = st.focusedNode;
      const target = focused && focused !== st.hroot && focused.parent ? focused.parent : st.hroot;
      _zoomTo(st, target);
    })
      .on("mouseover", (ev) => showTip("Click to go up one level", ev))
      .on("mousemove", moveTip)
      .on("mouseout", hideTip);
    // Entrance animation — arcs fan in from center
    st.paths.each(function (d) {
      this._x0 = d.x0;
      this._x1 = d.x0;
      this._y0 = d.y0;
      this._y1 = d.y1;
    });
    st.paths
      .transition()
      .duration(700)
      .attrTween("d", function (d) {
        const i = d3.interpolate({ x0: this._x0, x1: this._x1, y0: this._y0, y1: this._y1 }, d);
        return (t) => st.arc(i(t));
      });
    // Initial state
    st.focusedNode = st.hroot;
    if (st.titleEl) st.titleEl.textContent = sampleFilter || "All Samples";
    if (st.bcEl) st.bcEl.textContent = "root";
    _updateList(st, st.hroot);
    _updateLegend(st);
  }

  // ── Create a new panel card ───────────────────────────────────────────
  function _createPanel(sampleValue) {
    const container = document.getElementById("sun-panels-container");
    if (!container) return;
    const id = _nextId++;
    const card = document.createElement("div");
    card.className = "sun-panel";
    card.id = `sun-panel-${id}`;
    card.innerHTML = `
      <div class="sun-panel-filterbar">
        <select class="sun-sample-sel"><option value="">All Samples</option></select>
        <button class="sun-btn-ghost sun-reset-btn" title="Reset to root">↺</button>
        <input class="sun-panel-search" type="text" placeholder="🔍 Search…">
        <span class="sun-badge sun-badge-el">0</span>
        <button class="sun-remove-btn sun-remove-btn-el">✕</button>
      </div>
      <div class="sun-charthead sun-title-el">All Samples</div>
      <div class="sun-breadcrumb sun-bc-el">root</div>
      <div class="sun-panel-svg-wrap"></div>
      <div class="sun-legend sun-legend-el"></div>
      <div class="sun-panel-list sun-list-el"></div>`;
    container.appendChild(card);
    const st = {
      id,
      panelEl: card,
      svgWrapEl: card.querySelector(".sun-panel-svg-wrap"),
      sampleSel: card.querySelector(".sun-sample-sel"),
      resetBtn: card.querySelector(".sun-reset-btn"),
      searchEl: card.querySelector(".sun-panel-search"),
      badgeEl: card.querySelector(".sun-badge-el"),
      removeBtn: card.querySelector(".sun-remove-btn-el"),
      titleEl: card.querySelector(".sun-title-el"),
      bcEl: card.querySelector(".sun-bc-el"),
      legendEl: card.querySelector(".sun-legend-el"),
      listEl: card.querySelector(".sun-list-el"),
      hroot: null,
      focusedNode: null,
      colorScale: null,
      _zoomColorScale: null,
      _defaultSample: sampleValue || null,
      metric: null,
      maxR: null,
      paths: null,
      xSc: null,
      ySc: null,
      arc: null,
      svg: null,
      g: null,
    };
    _panels[id] = st;
    st.resetBtn.addEventListener("click", () => {
      if (st.hroot) _zoomTo(st, st.hroot, false);
    });
    st.removeBtn.addEventListener("click", () => _removePanel(id));
    st.sampleSel.addEventListener("change", () => _renderPanel(id));
    let _deb;
    st.searchEl.addEventListener("input", () => {
      clearTimeout(_deb);
      _deb = setTimeout(() => {
        _globalSearch = st.searchEl.value;
        // Sync search text to all other panels
        Object.keys(_panels).forEach((oid) => {
          const ost = _panels[Number(oid)];
          if (ost && ost !== st && ost.searchEl) ost.searchEl.value = _globalSearch;
        });
        // Re-render all panels (search reshapes the hierarchy)
        Object.keys(_panels).forEach((oid) => _renderPanel(Number(oid)));
        // Restore global zoom position in each panel after re-render
        if (_globalZoomPath.length > 0) {
          Object.keys(_panels).forEach((oid) => {
            const ost = _panels[Number(oid)];
            if (ost && ost.hroot) {
              const target = _findNodeByPath(ost.hroot, _globalZoomPath);
              if (target && target !== ost.hroot) _zoomTo(ost, target, true);
            }
          });
        }
      }, 250);
    });
    _refreshRemoveBtns();
    _renderPanel(id);
    return id;
  }

  function _removePanel(id) {
    if (Object.keys(_panels).length <= 1) return;
    _panels[id]?.panelEl.remove();
    delete _panels[id];
    _refreshRemoveBtns();
  }

  function _refreshRemoveBtns() {
    const count = Object.keys(_panels).length;
    Object.values(_panels).forEach((st) => {
      if (st.removeBtn) st.removeBtn.style.display = count <= 1 ? "none" : "";
    });
  }

  // ── Public: redraw all panels (called by global redraw()) ─────────────
  window.drawSunburst = function () {
    if (Object.keys(_panels).length === 0) {
      // Auto-create one panel per sample; fall back to a single panel if no samples found
      const samples = _orderedSamples(uniq(filteredData().map((r) => r["Specimen ID"] || "")).filter(Boolean));
      if (samples.length > 0) {
        samples.forEach((s) => _createPanel(s));
      } else {
        _createPanel();
      }
    } else {
      Object.keys(_panels).forEach((id) => _renderPanel(Number(id)));
    }
  };

  // ── Wire toolbar ──────────────────────────────────────────────────────
  const _addBtn = document.getElementById("sun-add-panel-btn");
  const _metricSel = document.getElementById("sun-metric");
  if (_addBtn) _addBtn.addEventListener("click", _createPanel);
  if (_metricSel) _metricSel.addEventListener("change", () => window.drawSunburst && window.drawSunburst());
})();
