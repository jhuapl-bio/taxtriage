/* ═══════════════════════════════════════════════════════════════════════════
       -  §  MAP GROUP NETWORK  (highlight → on-map similarity links + table)
       -     When a group is HIGHLIGHTED in the legend (first click of the
       -     normal → highlight → hidden cycle), the map stops being a plain
       -     scatter of markers and shows the structure inside that category:
       -
       -       • every member of the highlighted group(s) stays vivid while the
       -         rest of the map fades back (handled by _mapIsDimmed),
       -       • an edge is drawn between each pair of members whose organism
       -         profiles are similar enough, weighted by that similarity,
       -       • a table below the map lists those pairs — sortable, and
       -         collapsible so it never covers the map.
       -
       -     Edge definition reuses _mgNetSimilarity from the Group Network tab
       -     (cosine / Jaccard / shared-taxa) so the two views can never
       -     disagree about which sites are linked.
       -
       -     At the Country / State choropleth levels there are no per-sample
       -     points, so the same idea is applied to regions: non-matching
       -     regions fade, matching ones are outlined, and edges connect region
       -     centroids. See _mgGeoHighlightOverlay.
       -
       -       _mgMapNetworkModel()  — nodes + edges for the highlighted groups
       -       _mgRefreshMapNetwork() — redraw / tear down everything below
       -       _mgDrawMapEdges()      — the Leaflet polyline layer
       -       _mgRenderNetworkTable() — the panel under the map
═══════════════════════════════════════════════════════════════════════════ */

let _mgEdgeLayer = null; // L.LayerGroup holding the polylines
let _mgNetTableSort = { col: "value", asc: false };
let _mgNetTableOpen = true;

/* Per-sample organism profile at the requested level, restricted to samples
   that are actually on the map. Mirrors _mgGroupProfiles' value rules but
   keyed by sample rather than by group, since the on-map network links
   individual sites. */
function _mgSampleProfiles(sampleIds, level, basis) {
  const want = new Set(sampleIds);
  const taxonOf = (r) =>
    level === "genus"
      ? String(r["Genus Name"] || r["Genus"] || "").trim()
      : String(r["Species Name"] || r["Detected Organism"] || "").trim();
  const valOf = (r) => {
    if (basis === "tass") {
      const t = level === "genus" ? parseFloat(r["Genus TASS"]) : parseFloat(r["Species TASS"]);
      const fb = parseFloat(r["TASS Score"]);
      return Math.max(isNaN(t) ? 0 : t, isNaN(fb) ? 0 : fb);
    }
    if (basis === "coverage") {
      const b = parseFloat(r["Breadth %"]);
      const c = parseFloat(r["Coverage"]);
      return Math.max(isNaN(b) ? 0 : b, isNaN(c) ? 0 : c);
    }
    const n = parseFloat(r["# Reads Aligned"]);
    return isNaN(n) ? 1 : Math.max(n, 1);
  };
  const isMax = basis === "tass" || basis === "coverage";

  const prof = {};
  const taxa = new Set();
  const rows = typeof filteredData === "function" ? filteredData() : [];
  rows.forEach((r) => {
    const s = r["Specimen ID"] || "";
    if (!want.has(s)) return;
    const t = taxonOf(r);
    if (!t) return;
    const v = valOf(r);
    if (!prof[s]) prof[s] = {};
    prof[s][t] = prof[s][t] == null ? v : isMax ? Math.max(prof[s][t], v) : prof[s][t] + v;
    taxa.add(t);
  });
  return { profiles: prof, taxa: Array.from(taxa) };
}

/* Great-circle distance in km — shown in the table so "these two are similar"
   can be read against "and they are 900 km apart", which is the whole point of
   putting the network on a map rather than in an abstract layout. */
function _mgHaversineKm(a, b) {
  const R = 6371;
  const dLat = ((b[0] - a[0]) * Math.PI) / 180;
  const dLon = ((b[1] - a[1]) * Math.PI) / 180;
  const la1 = (a[0] * Math.PI) / 180;
  const la2 = (b[0] * Math.PI) / 180;
  const h = Math.sin(dLat / 2) ** 2 + Math.cos(la1) * Math.cos(la2) * Math.sin(dLon / 2) ** 2;
  return 2 * R * Math.asin(Math.min(1, Math.sqrt(h)));
}

/* The current on-map network. Returns null when nothing is highlighted, which
   is the signal to tear the overlay down entirely. */
function _mgMapNetworkModel() {
  const mg = window.metaGrouping;
  if (!mg || !mg.active() || !mg.anyHighlighted()) return null;
  if (!_markerObjects) return null;

  const level = (document.getElementById("mg-map-net-level") || {}).value || "genus";
  const basis = (document.getElementById("mg-map-net-basis") || {}).value || "presence";
  const simMode = (document.getElementById("mg-map-net-sim") || {}).value || "cosine";
  const thresh = parseFloat((document.getElementById("mg-map-net-threshold") || {}).value || "0.5");

  // Members of the highlighted groups that are actually plotted and visible.
  const nodes = [];
  const seen = new Set();
  Object.keys(_markerObjects).forEach((sn) => {
    const o = _markerObjects[sn];
    if (!o || !_mapMarkerVisible(sn, o.members)) return;
    if (!mg.sampleHighlighted(sn)) return;
    if (seen.has(sn)) return;
    seen.add(sn);
    const gk = mg.groupsOf(sn).filter((k) => mg.isHighlighted(k));
    nodes.push({
      id: sn,
      lat: o.lat,
      lon: o.lon,
      group: gk[0] || null,
      groups: gk,
      color: (gk[0] && mg.color(gk[0])) || o.color,
    });
  });

  const { profiles, taxa } = _mgSampleProfiles(
    nodes.map((n) => n.id),
    level,
    basis,
  );

  const cut = simMode === "shared" ? Math.max(1, Math.round(thresh * 20)) : thresh;
  const edges = [];
  for (let i = 0; i < nodes.length; i++) {
    for (let j = i + 1; j < nodes.length; j++) {
      const a = profiles[nodes[i].id] || {};
      const b = profiles[nodes[j].id] || {};
      const v = typeof _mgNetSimilarity === "function" ? _mgNetSimilarity(a, b, simMode, taxa) : 0;
      if (!(v >= cut) || v <= 0) continue;
      const shared = taxa.filter((t) => (a[t] || 0) > 0 && (b[t] || 0) > 0);
      edges.push({
        a: nodes[i],
        b: nodes[j],
        value: v,
        shared,
        km: _mgHaversineKm([nodes[i].lat, nodes[i].lon], [nodes[j].lat, nodes[j].lon]),
        crossGroup: nodes[i].group !== nodes[j].group,
      });
    }
  }
  edges.sort((x, y) => y.value - x.value);

  return { nodes, edges, simMode, level, basis, cut, groups: mg.highlighted(), taxaCount: taxa.length };
}

/* ── Leaflet overlay ────────────────────────────────────────────────────── */
function _mgClearMapEdges() {
  if (_mgEdgeLayer && _leafletMap && _leafletMap.hasLayer(_mgEdgeLayer)) _leafletMap.removeLayer(_mgEdgeLayer);
  _mgEdgeLayer = null;
}

function _mgDrawMapEdges(model) {
  _mgClearMapEdges();
  if (!_leafletMap || !model || !model.edges.length) return;
  const maxV = model.edges[0].value || 1;
  _mgEdgeLayer = L.layerGroup();

  model.edges.forEach((e) => {
    const w = 1 + 4 * (e.value / maxV);
    const op = 0.25 + 0.5 * (e.value / maxV);
    // A link between two DIFFERENT highlighted groups is the interesting case
    // (two categories sharing a signal), so it is dashed and darker rather
    // than taking one side's colour and implying it belongs to that group.
    const line = L.polyline(
      [
        [e.a.lat, e.a.lon],
        [e.b.lat, e.b.lon],
      ],
      {
        color: e.crossGroup ? "#37474f" : e.a.color,
        weight: w,
        opacity: op,
        dashArray: e.crossGroup ? "6 5" : null,
        interactive: true,
        className: "mg-map-edge",
      },
    );
    line.bindTooltip(
      `<b>${_mgEsc(e.a.id)} ↔ ${_mgEsc(e.b.id)}</b><br>` +
        `${model.simMode}: ${e.value.toFixed(3)}<br>` +
        `${e.shared.length} shared ${model.level === "genus" ? "genera" : "species"}<br>` +
        `${_mgKm(e.km)} apart`,
      { sticky: true },
    );
    line.on("click", () => _mgHighlightTableRow(e));
    _mgEdgeLayer.addLayer(line);
  });

  _mgEdgeLayer.addTo(_leafletMap);
  // Keep the links under the markers so a dot is always clickable.
  if (_mgEdgeLayer.eachLayer) _mgEdgeLayer.eachLayer((l) => l.bringToBack && l.bringToBack());
}

function _mgKm(km) {
  if (km == null || isNaN(km)) return "—";
  if (km < 1) return `${Math.round(km * 1000)} m`;
  if (km < 100) return `${km.toFixed(1)} km`;
  return `${Math.round(km).toLocaleString()} km`;
}

/* ── The table under the map ─────────────────────────────────────────────
   Deliberately NOT a popup over the map: it lives in its own panel below the
   precise view, is collapsible, and disappears entirely when no group is
   highlighted, so it never covers what the user is looking at. */
function _mgRenderNetworkTable(model) {
  const panel = document.getElementById("mg-map-net-panel");
  if (!panel) return;
  if (!model) {
    panel.style.display = "none";
    return;
  }
  panel.style.display = "";

  const title = document.getElementById("mg-map-net-title");
  if (title) {
    const gl = model.groups.map((g) => `<span class="mg-net-chip">${_mgEsc(g)}</span>`).join("");
    title.innerHTML =
      `<i class="fas fa-circle-nodes"></i> Network for ${gl} ` +
      `<span class="mg-net-sub">${model.nodes.length} site${model.nodes.length === 1 ? "" : "s"}, ` +
      `${model.edges.length} link${model.edges.length === 1 ? "" : "s"}</span>`;
  }

  const body = document.getElementById("mg-map-net-body");
  if (body) body.style.display = _mgNetTableOpen ? "" : "none";
  const caret = document.getElementById("mg-map-net-caret");
  if (caret) caret.className = _mgNetTableOpen ? "fas fa-chevron-down" : "fas fa-chevron-right";

  const wrap = document.getElementById("mg-map-net-table-wrap");
  if (!wrap) return;

  if (!model.nodes.length) {
    wrap.innerHTML =
      '<p class="mg-net-empty">No plotted samples in the highlighted group — it may have no ' +
      "latitude / longitude, or the sidebar filters may have excluded it.</p>";
    return;
  }
  if (!model.edges.length) {
    wrap.innerHTML =
      `<p class="mg-net-empty">${model.nodes.length} site${model.nodes.length === 1 ? "" : "s"} highlighted, ` +
      `but no pair reaches the link threshold (${model.simMode === "shared" ? model.cut + " shared taxa" : model.cut}). ` +
      "Lower the threshold above to see weaker links.</p>";
    return;
  }

  const rows = model.edges.slice();
  const s = _mgNetTableSort;
  rows.sort((x, y) => {
    let d = 0;
    if (s.col === "value") d = x.value - y.value;
    else if (s.col === "km") d = x.km - y.km;
    else if (s.col === "shared") d = x.shared.length - y.shared.length;
    else if (s.col === "a") d = String(x.a.id).localeCompare(String(y.a.id));
    else if (s.col === "b") d = String(x.b.id).localeCompare(String(y.b.id));
    return s.asc ? d : -d;
  });

  const th = (col, label, help) =>
    `<th data-col="${col}" class="mg-net-th${s.col === col ? " mg-net-sorted" : ""}" title="${help}">` +
    `${label}${s.col === col ? (s.asc ? " ▲" : " ▼") : ""}</th>`;

  const levelWord = model.level === "genus" ? "genera" : "species";
  // At the aggregated levels the rows are region pairs, not sample pairs, and
  // the distance column is meaningless (centroids are projected screen points).
  const isRegion = model.unit === "region";
  const A = isRegion ? "Region A" : "Site A";
  const B = isRegion ? "Region B" : "Site B";
  wrap.innerHTML =
    `<table class="mg-net-table"><thead><tr>` +
    th("a", A, "Sort by the first name") +
    th("b", B, "Sort by the second name") +
    th("value", model.simMode === "shared" ? "Shared taxa" : "Similarity", "Sort by link strength") +
    th("shared", `Shared ${levelWord}`, "Sort by how many taxa the pair has in common") +
    (isRegion ? "" : th("km", "Distance", "Sort by great-circle distance")) +
    `<th>Top shared ${levelWord}</th></tr></thead><tbody>` +
    rows
      .map((e) => {
        const top = e.shared.slice(0, 4).map((t) => `<span class="mg-net-taxon">${_mgEsc(t)}</span>`).join("");
        const more = e.shared.length > 4 ? `<span class="mg-net-more">+${e.shared.length - 4}</span>` : "";
        return (
          `<tr data-a="${_mgEsc(e.a.id)}" data-b="${_mgEsc(e.b.id)}"${e.crossGroup ? ' class="mg-net-cross"' : ""}>` +
          `<td>${_mgEsc(e.a.id)}</td><td>${_mgEsc(e.b.id)}</td>` +
          `<td class="mg-net-num">${model.simMode === "shared" ? e.value : e.value.toFixed(3)}</td>` +
          `<td class="mg-net-num">${e.shared.length}</td>` +
          (isRegion ? "" : `<td class="mg-net-num">${_mgKm(e.km)}</td>`) +
          `<td>${top}${more}</td></tr>`
        );
      })
      .join("") +
    `</tbody></table>` +
    (model.edges.some((e) => e.crossGroup)
      ? '<p class="mg-net-note"><span class="mg-net-dash"></span> Dashed rows / links join two <em>different</em> ' +
        "highlighted groups — a shared signal across categories.</p>"
      : "");

  wrap.querySelectorAll(".mg-net-th").forEach((h) => {
    h.addEventListener("click", () => {
      const c = h.getAttribute("data-col");
      _mgNetTableSort = { col: c, asc: _mgNetTableSort.col === c ? !_mgNetTableSort.asc : c === "a" || c === "b" };
      _mgRenderNetworkTable(model);
    });
  });
  // Hovering a row pans nothing but flashes the pair on the map, so the table
  // and the map stay legible as one view rather than two.
  wrap.querySelectorAll("tbody tr").forEach((tr) => {
    tr.addEventListener("mouseenter", () => _mgFlashPair(tr.getAttribute("data-a"), tr.getAttribute("data-b")));
    tr.addEventListener("mouseleave", () => _mgFlashPair(null, null));
  });
}

// Emphasise one edge's endpoints on hover by fading the other polylines.
function _mgFlashPair(a, b) {
  if (!_mgEdgeLayer || !_mgEdgeLayer.eachLayer) return;
  _mgEdgeLayer.eachLayer((l) => {
    if (!l.setStyle || !l.getLatLngs) return;
    if (!a) {
      l.setStyle({ opacity: l.options.opacity });
      return;
    }
    const ll = l.getLatLngs();
    const o = _markerObjects || {};
    const match =
      o[a] && o[b] && ll.length === 2
        ? (Math.abs(ll[0].lat - o[a].lat) < 1e-9 && Math.abs(ll[1].lat - o[b].lat) < 1e-9) ||
          (Math.abs(ll[0].lat - o[b].lat) < 1e-9 && Math.abs(ll[1].lat - o[a].lat) < 1e-9)
        : false;
    l.setStyle({ opacity: match ? 0.95 : 0.08 });
  });
}

function _mgHighlightTableRow(e) {
  const wrap = document.getElementById("mg-map-net-table-wrap");
  if (!wrap) return;
  wrap.querySelectorAll("tbody tr").forEach((tr) => tr.classList.remove("mg-net-rowsel"));
  const tr = wrap.querySelector(`tbody tr[data-a="${CSS.escape(e.a.id)}"][data-b="${CSS.escape(e.b.id)}"]`);
  if (tr) {
    tr.classList.add("mg-net-rowsel");
    tr.scrollIntoView({ block: "nearest" });
  }
}

/* ═══════════════ CHOROPLETH HIGHLIGHT (Country / State levels) ═══════════
   The aggregated levels have no per-sample points, so the same three ideas are
   expressed on regions instead:
     • regions containing no highlighted sample fade back,
     • regions that do are outlined in the group's colour,
     • edges join region CENTROIDS, so a link between two countries reads as
       "these two territories share a signal" rather than implying a
       sample-level path.
   The table below the map switches to region pairs to match.                 */
function _mgGeoRegionModel() {
  const ctx = window._geoLastDraw;
  const mg = window.metaGrouping;
  if (!ctx || !mg || !mg.active() || !mg.anyHighlighted()) return null;

  const level = (document.getElementById("mg-map-net-level") || {}).value || "genus";
  const basis = (document.getElementById("mg-map-net-basis") || {}).value || "presence";
  const simMode = (document.getElementById("mg-map-net-sim") || {}).value || "cosine";
  const thresh = parseFloat((document.getElementById("mg-map-net-threshold") || {}).value || "0.5");

  // Region → the highlighted samples inside it.
  const regions = [];
  (ctx.features || []).forEach((f) => {
    const samples =
      typeof window._geoFeatureSamples === "function"
        ? window._geoFeatureSamples(ctx.field, f, ctx.metaByNorm, ctx.stateCountries)
        : [];
    if (!samples.length) return;
    const hot = samples.filter((s) => mg.sampleHighlighted(s));
    if (!hot.length) return;
    const c = ctx.path && ctx.path.centroid ? ctx.path.centroid(f) : null;
    if (!c || isNaN(c[0]) || isNaN(c[1])) return;
    const p = f.properties || {};
    const gk = mg.groupsOf(hot[0]).filter((k) => mg.isHighlighted(k));
    regions.push({
      id: p.name || p.NAME || p.admin || "region",
      cx: c[0],
      cy: c[1],
      samples: hot,
      group: gk[0] || null,
      color: (gk[0] && mg.color(gk[0])) || "#1565c0",
      feature: f,
    });
  });

  if (!regions.length) return { ctx, regions: [], edges: [], simMode, level, basis, cut: thresh };

  // Aggregate each region's member samples into one profile, then compare.
  const allSamples = [];
  regions.forEach((r) => r.samples.forEach((s) => allSamples.push(s)));
  const { profiles, taxa } = _mgSampleProfiles(allSamples, level, basis);
  const regionProf = regions.map((r) => {
    const p = {};
    r.samples.forEach((s) => {
      const sp = profiles[s] || {};
      Object.keys(sp).forEach((t) => (p[t] = Math.max(p[t] || 0, sp[t])));
    });
    return p;
  });

  const cut = simMode === "shared" ? Math.max(1, Math.round(thresh * 20)) : thresh;
  const edges = [];
  for (let i = 0; i < regions.length; i++) {
    for (let j = i + 1; j < regions.length; j++) {
      const v =
        typeof _mgNetSimilarity === "function" ? _mgNetSimilarity(regionProf[i], regionProf[j], simMode, taxa) : 0;
      if (!(v >= cut) || v <= 0) continue;
      const shared = taxa.filter((t) => (regionProf[i][t] || 0) > 0 && (regionProf[j][t] || 0) > 0);
      edges.push({
        a: regions[i],
        b: regions[j],
        value: v,
        shared,
        km: NaN, // projected centroids: a screen-space line, not a measured distance
        crossGroup: regions[i].group !== regions[j].group,
      });
    }
  }
  edges.sort((x, y) => y.value - x.value);
  return { ctx, regions, edges, simMode, level, basis, cut, groups: mg.highlighted() };
}

function _mgGeoHighlightOverlay() {
  const ctx = window._geoLastDraw;
  if (!ctx || !ctx.svg) return;
  // Remove any previous overlay before re-deciding.
  ctx.svg.selectAll("g.mg-geo-links").remove();

  const mg = window.metaGrouping;
  const on = !!(mg && mg.active() && mg.anyHighlighted());
  if (!on) {
    // Restore the plain choropleth look.
    if (ctx.sel) ctx.sel.attr("opacity", 1).attr("stroke-dasharray", null);
    return;
  }

  const model = _mgGeoRegionModel();
  const hot = new Map((model && model.regions ? model.regions : []).map((r) => [r.feature, r]));

  ctx.sel
    .attr("opacity", (f) => (hot.has(f) ? 1 : 0.25))
    .attr("stroke", (f) => (hot.has(f) ? hot.get(f).color : "#ffffff"))
    .attr("stroke-width", (f) => (hot.has(f) ? 2.4 : 0.4));

  if (!model || !model.edges.length) return;
  const maxV = model.edges[0].value || 1;
  const g = ctx.svg.append("g").attr("class", "mg-geo-links").attr("pointer-events", "none");
  model.edges.forEach((e) => {
    g.append("line")
      .attr("x1", e.a.cx)
      .attr("y1", e.a.cy)
      .attr("x2", e.b.cx)
      .attr("y2", e.b.cy)
      .attr("stroke", e.crossGroup ? "#37474f" : e.a.color)
      .attr("stroke-width", 1 + 3.5 * (e.value / maxV))
      .attr("stroke-opacity", 0.3 + 0.5 * (e.value / maxV))
      .attr("stroke-dasharray", e.crossGroup ? "6 5" : null);
  });
  model.regions.forEach((r) => {
    g.append("circle")
      .attr("cx", r.cx)
      .attr("cy", r.cy)
      .attr("r", 4)
      .attr("fill", r.color)
      .attr("stroke", "#fff")
      .attr("stroke-width", 1.2);
  });
}
window._mgGeoHighlightOverlay = _mgGeoHighlightOverlay;

/* ── Entry point ─────────────────────────────────────────────────────────
   Called from the grouping engine's change broadcast and from the overlay's
   own controls. Safe to call when nothing is highlighted: it tears down. */
function _mgRefreshMapNetwork() {
  const level = (document.getElementById("geo-cmp-level") || {}).value || "precise";
  const aggregated = level !== "precise";

  // Point map: markers + polylines. Aggregated levels: regions + centroid
  // links, drawn into the choropleth's own SVG.
  const pointModel = aggregated ? null : _mgMapNetworkModel();
  _mgDrawMapEdges(pointModel);
  _mgGeoHighlightOverlay();

  const model = aggregated ? _mgGeoRegionModel() : pointModel;
  // The table is the same panel either way; only the row nouns change.
  _mgRenderNetworkTable(
    model && aggregated
      ? Object.assign({}, model, { nodes: model.regions, unit: "region" })
      : model && Object.assign({}, model, { unit: "site" }),
  );
  return model;
}
window._mgRefreshMapNetwork = _mgRefreshMapNetwork;

/* ── Controls ───────────────────────────────────────────────────────────── */
(function () {
  function ready() {
    ["mg-map-net-level", "mg-map-net-basis", "mg-map-net-sim", "mg-map-net-threshold"].forEach((id) => {
      const el = document.getElementById(id);
      if (!el) return;
      el.addEventListener(el.type === "range" ? "input" : "change", () => {
        const lbl = document.getElementById("mg-map-net-threshold-val");
        const sim = (document.getElementById("mg-map-net-sim") || {}).value || "cosine";
        const v = parseFloat((document.getElementById("mg-map-net-threshold") || {}).value || "0.5");
        if (lbl) lbl.textContent = sim === "shared" ? `≥ ${Math.max(1, Math.round(v * 20))} taxa` : v.toFixed(2);
        _mgRefreshMapNetwork();
      });
    });
    const head = document.getElementById("mg-map-net-head");
    if (head) {
      head.addEventListener("click", () => {
        _mgNetTableOpen = !_mgNetTableOpen;
        _mgRefreshMapNetwork();
      });
    }
    const clr = document.getElementById("mg-map-net-clear");
    if (clr) {
      clr.addEventListener("click", (ev) => {
        ev.stopPropagation(); // the header toggles collapse; this must not
        if (window.metaGrouping) window.metaGrouping.clearHighlight();
      });
    }
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", ready);
  else ready();
})();
