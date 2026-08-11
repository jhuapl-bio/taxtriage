/* ═══════════════════════════════════════════════════════════════════════════
       -  §  MAP  (Precise level of the Mapping & Geography sub-tab)
       -     Leaflet map with CARTO Voyager tiles (no Referer restrictions).
       -     Renders one marker per sample with lat/lon metadata; clicking a
       -     marker opens a side-panel with a filtered TASS table for that
       -     sample. Group mode (multi-marker selection) shows a combined
       -     panel with sortable headers.
       -     Key functions: _initMap / _doInitMap, _rebuildMapMarkers,
       -     _refreshMapMarkerColors, _renderMapPanel / _renderMapGroupPanel,
       -     _setMapPanelGroupMode, _refreshMapPanelTable, closeMapPanel,
       -     viewSampleInMetaTab, editSampleNameFromSidebar / _renameSample.
       -     State variables and drag listeners are declared earlier (before
       -     redraw()) to avoid the `let` Temporal Dead Zone.
═══════════════════════════════════════════════════════════════════════════ */
/* `dimmed` renders a sample that does not match the active filters: smaller,
   faint grey, no shadow. It stays on the map on purpose — seeing that the four
   RSV-A hits sit in one place is only meaningful next to the samples that had
   no hit, so hiding the non-matches would throw away the denominator. */
function _svgDot(color, selected, dimmed, shape) {
  const r = dimmed ? 5 : selected ? 11 : 8;
  const fill = dimmed ? "#b9c3cd" : color;
  const cx = r + 5;
  const ring =
    selected && !dimmed
      ? `<circle cx="${cx}" cy="${cx}" r="${
          r + 3
        }" fill="none" stroke="${color}" stroke-width="2" opacity=".45"/>`
      : "";
  // Shape coding: when a metadata grouping is active every group gets its own
  // colour AND its own glyph, so the categories stay separable past the
  // 15-colour palette and in greyscale print.
  const glyph =
    shape && shape !== "circle" && typeof window._mgShapeEl === "function"
      ? `${window._mgShapeEl(shape, cx, cx, r)} fill="${fill}" ` +
        `stroke="rgba(255,255,255,${dimmed ? ".7" : ".9"})" stroke-width="${dimmed ? 1 : selected ? 2.5 : 1.5}" ` +
        `${dimmed ? "" : 'filter="drop-shadow(0 1px 3px rgba(0,0,0,.38))"'} />`
      : `<circle cx="${cx}" cy="${cx}" r="${r}" fill="${fill}"
              stroke="rgba(255,255,255,${dimmed ? ".7" : ".9"})" stroke-width="${dimmed ? 1 : selected ? 2.5 : 1.5}"
              ${dimmed ? "" : 'filter="drop-shadow(0 1px 3px rgba(0,0,0,.38))"'}/>`;
  return L.divIcon({
    className: "",
    html: `<svg xmlns="http://www.w3.org/2000/svg" width="${cx * 2}" height="${cx * 2}"${
      dimmed ? ' opacity="0.45"' : ""
    }>
            ${glyph}
            ${ring}
          </svg>`,
    iconSize: [cx * 2, cx * 2],
    iconAnchor: [cx, cx],
  });
}

/* ── Group styling for markers ───────────────────────────────────────────
   With no grouping selected the map behaves exactly as before: one colour per
   sample, drawn from `sampleColors`. Once the user picks columns in the shared
   "Group by" bar, markers instead take the colour and shape of the group their
   sample belongs to, which is what turns the map into a categorical comparison
   (Agricultural vs Water vs Soil) rather than a per-sample scatter. */
function _mapGroupStyle(sn) {
  const mg = window.metaGrouping;
  if (mg && mg.active()) {
    const keys = mg.groupsOf(sn);
    if (keys.length) {
      return {
        grouped: true,
        group: keys[0],
        groups: keys,
        color: mg.color(keys[0]),
        shape: mg.shape(keys[0]),
      };
    }
    // Sample carries no value for any grouping column — draw it neutrally so
    // it is visibly "ungrouped" rather than silently miscoloured.
    return { grouped: true, group: null, groups: [], color: "#b0bec5", shape: "circle" };
  }
  return {
    grouped: false,
    group: null,
    groups: [],
    color: (typeof sampleColors !== "undefined" && sampleColors[sn]) || "#1565C0",
    shape: "circle",
  };
}

/* ── Filter scoping ──────────────────────────────────────────────────────
   The map used to read RUN_META directly, so sidebar filters never reached it.
   When _mapFilterScoped is on (default) a sample whose detections are all
   filtered out is drawn dimmed. Only applied while a filter is actually
   narrowing something — otherwise every marker would trivially "match". */
let _mapFilterScoped = true;
let _mapMatchedKeys = null;

function _mapRefreshMatchSet() {
  const cb = document.getElementById("map-filter-scope");
  if (cb) _mapFilterScoped = !!cb.checked;
  const active = typeof _ttAnyFilterActive === "function" && _ttAnyFilterActive();
  _mapMatchedKeys =
    _mapFilterScoped && active && typeof _ttFilterMatchedKeys === "function" ? _ttFilterMatchedKeys() : null;
  // Small "(4 of 164 samples)" readout beside the toggle: without it a fully
  // dimmed map is ambiguous between "nothing matched" and "scoping is off".
  const out = document.getElementById("map-filter-scope-count");
  if (out) {
    if (!_mapMatchedKeys) {
      out.textContent = _mapFilterScoped ? "(no active filter)" : "";
    } else {
      const geo = (typeof RUN_META !== "undefined" ? RUN_META : []).filter(
        (r) => r.latitude != null && r.longitude != null,
      );
      const hit = geo.filter((r) => _ttSampleMatchesFilter(r.sample_name, _mapMatchedKeys)).length;
      out.textContent = `(${hit} of ${geo.length} mapped sample(s))`;
    }
  }
  return _mapMatchedKeys;
}

/* Wire the map's filter-scope checkbox once. Re-icons the existing markers
   rather than rebuilding the layer, so the current zoom / selection survive. */
function _wireMapFilterScope() {
  const cb = document.getElementById("map-filter-scope");
  if (!cb || cb._wired) return;
  cb._wired = true;
  cb.addEventListener("change", () => {
    _mapFilterScoped = !!cb.checked;
    _refreshMapMarkerColors();
  });
}

/* ── Map-local sample visibility ─────────────────────────────────────────
   A show/hide list that affects THIS MAP ONLY. Deliberately separate from the
   sidebar's `sampleHidden`: that one removes a sample from every tab, which is
   the wrong tool when you just want to declutter the map while keeping the
   sample in the table, heatmap and cross-sample views.

   Two independent mechanisms end up on the same markers, and they mean
   different things:
     • filter scoping (above) DIMS a sample the sidebar filters exclude —
       context you still want to see.
     • this list HIDES a sample outright — decluttering you asked for.
   `sampleHidden` still wins over both, since that is a global "this sample is
   out of the analysis" statement.

   Stores the HIDDEN names, so the default (empty set) shows everything and new
   samples arriving from an upload are visible without touching this state. */
const _mapHiddenSamples = new Set();

/** Is this marker hidden by the map-local list? Keyed on the marker's own name,
 *  which is the specimen when merge is on and the raw sample otherwise. */
function _mapLocallyHidden(sn) {
  return _mapHiddenSamples.has(String(sn));
}

/** Everything that decides whether a marker is on the map at all. */
function _mapMarkerVisible(sn, members) {
  const list = Array.isArray(members) && members.length ? members : [sn];
  if (list.every((s) => sampleHidden[s])) return false; // global sidebar hide
  if (_mapLocallyHidden(sn)) return false; // map-only hide
  // Group toggled off in the shared legend. A sample in several groups stays
  // on the map while any one of them is still switched on.
  try {
    if (window.metaGrouping && window.metaGrouping.active() && !window.metaGrouping.sampleVisible(sn)) return false;
  } catch (e) {}
  return true;
}

/** True when this marker should be drawn faded. Two independent reasons:
 *   • the sample falls outside the active sidebar filters, or
 *   • a group is highlighted and this sample is in none of the highlighted
 *     groups — the "click once to highlight, fade all others" state.
 *  A merged specimen counts as matching if ANY of its member libraries does. */
function _mapIsDimmed(sn, members) {
  try {
    if (window.metaGrouping && window.metaGrouping.sampleDimmed(sn)) return true;
  } catch (e) {}
  if (!_mapMatchedKeys) return false;
  const list = Array.isArray(members) && members.length ? members : [sn];
  if (typeof _ttSampleMatchesFilter !== "function") return false;
  if (_ttSampleMatchesFilter(sn, _mapMatchedKeys)) return false;
  return !list.some((m) => _ttSampleMatchesFilter(m, _mapMatchedKeys));
}

/* ── Pie-chart marker for co-located samples ──────────────────────── */
function _pieSvg(colors, selected) {
  const r = selected ? 11 : 9;
  const cx = r + 5,
    cy = r + 5,
    n = colors.length;
  if (n === 1) return _svgDot(colors[0], selected);
  let paths = "";
  for (let i = 0; i < n; i++) {
    const a0 = (i / n) * 2 * Math.PI - Math.PI / 2;
    const a1 = ((i + 1) / n) * 2 * Math.PI - Math.PI / 2;
    const x0 = cx + r * Math.cos(a0),
      y0 = cy + r * Math.sin(a0);
    const x1 = cx + r * Math.cos(a1),
      y1 = cy + r * Math.sin(a1);
    const large = a1 - a0 > Math.PI ? 1 : 0;
    paths += `<path d="M${cx},${cy} L${x0.toFixed(2)},${y0.toFixed(2)} A${r},${r} 0 ${large} 1 ${x1.toFixed(
      2,
    )},${y1.toFixed(2)} Z" fill="${colors[i]}" stroke="white" stroke-width="0.8"/>`;
  }
  const ring = selected
    ? `<circle cx="${cx}" cy="${cy}" r="${r + 3}" fill="none" stroke="#333" stroke-width="1.5" opacity=".35"/>`
    : "";
  const size = (r + 5) * 2;
  return L.divIcon({
    className: "",
    html: `<svg xmlns="http://www.w3.org/2000/svg" width="${size}" height="${size}" style="filter:drop-shadow(0 1px 3px rgba(0,0,0,.38))">${paths}${ring}</svg>`,
    iconSize: [size, size],
    iconAnchor: [cx, cy],
  });
}

/* ── Cluster bubble icon ─────────────────────────────────────────────
   Overlapping sample dots merge into one bubble showing the sample count.
   Size grows with the count so dense clusters read clearly; the bubble
   splits back into individual dots as the user zooms in. */
function _clusterIcon(cluster) {
  const n = cluster.getChildCount();
  const size = n < 10 ? 34 : n < 100 ? 42 : 50;
  const fs = n < 100 ? 14 : 12;
  // With a filter active, a bubble reading just "40" hides the thing you are
  // looking for. Show "4/40" and grey out bubbles with no match, so a cluster
  // containing the hits is identifiable without zooming in first.
  let label = String(n);
  let bg = "rgba(21,101,192,.88)";
  let title = `${n} sample(s)`;
  if (_mapMatchedKeys) {
    let hit = 0;
    const kids = typeof cluster.getAllChildMarkers === "function" ? cluster.getAllChildMarkers() : [];
    kids.forEach((m) => {
      if (m && !m.ttDimmed) hit++;
    });
    label = `${hit}/${n}`;
    bg = hit ? "rgba(21,101,192,.92)" : "rgba(150,163,176,.55)";
    title = `${hit} of ${n} sample(s) match the active filters`;
  }
  // ── Group-composition bubble ──────────────────────────────────────────
  // With a grouping active, a flat count bubble throws away the one thing the
  // user asked the map to show. Draw the cluster as a donut whose wedges are
  // the group mix inside it, so "this cluster is 3 Water / 1 Soil" is legible
  // without zooming in.
  const mg = window.metaGrouping;
  if (mg && mg.active()) {
    const kids = typeof cluster.getAllChildMarkers === "function" ? cluster.getAllChildMarkers() : [];
    const tally = new Map();
    kids.forEach((m) => {
      const g = m && m.ttGroup != null ? m.ttGroup : "(ungrouped)";
      tally.set(g, (tally.get(g) || 0) + 1);
    });
    const parts = Array.from(tally.entries()).sort((a, b) => b[1] - a[1]);
    if (parts.length) {
      const R = size / 2;
      const cx = R;
      const cy = R;
      const rOuter = R - 2;
      const rInner = rOuter * 0.56;
      let acc = 0;
      let wedges = "";
      parts.forEach(([g, cnt]) => {
        const a0 = (acc / n) * 2 * Math.PI - Math.PI / 2;
        acc += cnt;
        const a1 = (acc / n) * 2 * Math.PI - Math.PI / 2;
        const col = g === "(ungrouped)" ? "#b0bec5" : mg.color(g);
        const large = a1 - a0 > Math.PI ? 1 : 0;
        // A single-group cluster has no arc to draw — use a full ring instead.
        if (parts.length === 1) {
          wedges = `<circle cx="${cx}" cy="${cy}" r="${(rOuter + rInner) / 2}" fill="none" stroke="${col}" stroke-width="${rOuter - rInner}"/>`;
          return;
        }
        const p = (r, a) => `${(cx + r * Math.cos(a)).toFixed(2)},${(cy + r * Math.sin(a)).toFixed(2)}`;
        wedges +=
          `<path d="M${p(rOuter, a0)} A${rOuter},${rOuter} 0 ${large} 1 ${p(rOuter, a1)} ` +
          `L${p(rInner, a1)} A${rInner},${rInner} 0 ${large} 0 ${p(rInner, a0)} Z" ` +
          `fill="${col}" stroke="#fff" stroke-width="0.9"/>`;
      });
      const tip = parts.map(([g, c]) => `${g}: ${c}`).join("\n");
      const gh =
        `<div title="${_mgEscAttr(title + "\n" + tip)}" style="position:relative;width:${size}px;height:${size}px">` +
        `<svg xmlns="http://www.w3.org/2000/svg" width="${size}" height="${size}" ` +
        `style="filter:drop-shadow(0 1px 4px rgba(0,0,0,.4))">` +
        `<circle cx="${R}" cy="${R}" r="${rInner}" fill="rgba(255,255,255,.96)"/>${wedges}</svg>` +
        `<div style="position:absolute;inset:0;display:flex;align-items:center;justify-content:center;` +
        `color:#263238;font-weight:700;font-size:${fs - 1}px;font-family:system-ui,sans-serif">${label}</div></div>`;
      return L.divIcon({ html: gh, className: "tt-cluster", iconSize: [size, size] });
    }
  }

  const html =
    `<div title="${title}" style="width:${size}px;height:${size}px;border-radius:50%;` +
    `background:${bg};border:2px solid rgba(255,255,255,.95);` +
    `box-shadow:0 1px 4px rgba(0,0,0,.4);display:flex;align-items:center;` +
    `justify-content:center;color:#fff;font-weight:700;font-size:${fs}px;` +
    `font-family:system-ui,sans-serif">${label}</div>`;
  return L.divIcon({ html, className: "tt-cluster", iconSize: [size, size] });
}

function _mgEscAttr(s) {
  return String(s == null ? "" : s).replace(/&/g, "&amp;").replace(/"/g, "&quot;").replace(/</g, "&lt;");
}

/* Re-fit the map to whatever is currently visible. Called after a group
   toggle so isolating a category actually zooms to it instead of leaving the
   user looking at a mostly-empty world map. */
function _mgMapFitVisible() {
  if (!_leafletMap || !_markerObjects) return;
  const pts = Object.keys(_markerObjects)
    .filter((n) => {
      const o = _markerObjects[n];
      return o && _mapMarkerVisible(n, o.members);
    })
    .map((n) => [_markerObjects[n].lat, _markerObjects[n].lon]);
  if (pts.length) _leafletMap.fitBounds(pts, { padding: [40, 40], maxZoom: 11 });
}
window._mgMapFitVisible = _mgMapFitVisible;

/* Isolate one group on the map. Delegates to the shared engine so the Group
   Heatmap and Group Network isolate with it — the grouping is one shared
   selection, so its visibility should be shared too. Re-soloing the group
   that is already alone restores the rest. */
function _mgMapFocusGroup(key) {
  const mg = window.metaGrouping;
  if (!mg || !mg.active() || !key) return;
  mg.solo(key); // broadcasts → markers, legends and the group views re-render
  _mgMapFitVisible();
  const picker = document.getElementById("map-sample-picker");
  if (picker && picker.open && typeof _mapRenderSampleList === "function") _mapRenderSampleList();
}
window._mgMapFocusGroup = _mgMapFocusGroup;

/* Create the marker container. Uses Leaflet.markercluster when available
   (zoom-based clustering with count bubbles); falls back to a plain layer
   group if the plugin failed to load (e.g. offline). */
function _makeMarkerLayer() {
  if (typeof L.markerClusterGroup === "function") {
    return L.markerClusterGroup({
      maxClusterRadius: 45, // px — dots closer than this merge into a bubble
      // We handle cluster clicks ourselves (see the "clusterclick" listener in
      // _doInitMap) so the behaviour can switch between list / expand modes.
      spiderfyOnMaxZoom: false,
      showCoverageOnHover: false,
      zoomToBoundsOnClick: false,
      iconCreateFunction: _clusterIcon,
    });
  }
  return L.layerGroup();
}

/* Set what happens when a cluster bubble is clicked (see the "clusterclick"
   handler in _doInitMap, which reads this):
     "list"   — open the side panel with every sample in the cluster (the
                pre-clustering behaviour, default).
     "expand" — zoom in to split the cluster, or fan out (spiderfy) dots that
                share exact coordinates. */
function _setMapClusterMode(mode) {
  _mapClusterMode = mode === "expand" ? "expand" : "list";
}

/* Small on-map control to switch between the two cluster-click modes. */
function _addMapClusterControl() {
  if (typeof L.Control === "undefined" || _mapClusterCtlAdded) return;
  const Ctl = L.Control.extend({
    options: { position: "topright" },
    onAdd: function () {
      const div = L.DomUtil.create("div", "leaflet-bar tt-cluster-mode");
      div.style.cssText =
        "background:#fff;padding:5px 8px;font:12px system-ui,sans-serif;color:#333;" +
        "line-height:1.3;box-shadow:0 1px 4px rgba(0,0,0,.25);border-radius:4px";
      div.innerHTML =
        '<div style="font-weight:600;margin-bottom:3px">On cluster click</div>' +
        '<label style="display:block;cursor:pointer"><input type="radio" name="tt-cluster-mode" value="list" checked> List samples</label>' +
        '<label style="display:block;cursor:pointer"><input type="radio" name="tt-cluster-mode" value="expand"> Zoom / expand</label>';
      L.DomEvent.disableClickPropagation(div);
      div
        .querySelectorAll('input[name="tt-cluster-mode"]')
        .forEach((el) => el.addEventListener("change", (e) => _setMapClusterMode(e.target.value)));
      return div;
    },
  });
  _leafletMap.addControl(new Ctl());
  _mapClusterCtlAdded = true;
}

/* Build one colored dot per geo-located sample and add it to the marker
   layer. Shared by _doInitMap and _rebuildMapMarkers. Returns lat/lon
   bounds for fitting the view. */
function _addSampleMarkers() {
  _wireMapFilterScope();
  _wireMapSamplePicker();
  _mapRefreshMatchSet();
  const geoRows = RUN_META.filter((r) => r.latitude != null && r.longitude != null);
  const mergeOn =
    typeof specimenMergeEnabled !== "undefined" && specimenMergeEnabled && typeof specimenOf === "function";
  const markerRows = [];
  if (mergeOn) {
    const bySpec = new Map();
    geoRows.forEach((r) => {
      const spec = specimenOf(r.sample_name);
      if (!bySpec.has(spec)) bySpec.set(spec, []);
      bySpec.get(spec).push(r);
    });
    bySpec.forEach((recs, spec) => {
      if (recs.length === 1 && spec === recs[0].sample_name) {
        markerRows.push(recs[0]);
        return;
      }
      const resolved =
        typeof SPECIMEN_META_RESOLVED !== "undefined" && SPECIMEN_META_RESOLVED[spec]
          ? SPECIMEN_META_RESOLVED[spec]
          : {};
      const rec = Object.assign({}, recs[0], resolved, {
        sample_name: spec,
        __mergedMembers:
          typeof specimenGroups === "function" && specimenGroups().has(spec)
            ? specimenGroups().get(spec).slice()
            : recs.map((r) => r.sample_name),
      });
      markerRows.push(rec);
    });
  } else {
    markerRows.push(...geoRows);
  }

  const sampleNames = [...new Set(markerRows.map((r) => r.sample_name))];
  sampleNames.forEach((n, i) => {
    if (!sampleColors[n]) sampleColors[n] = PALETTE[i % PALETTE.length];
  });

  const bounds = [];
  markerRows.forEach((rec) => {
    const lat = parseFloat(rec.latitude);
    const lon = parseFloat(rec.longitude);
    if (isNaN(lat) || isNaN(lon)) return;
    const sn = rec.sample_name;
    const style = _mapGroupStyle(sn);
    const color = style.color;
    const selected = _selectedSample === sn;
    const members0 = Array.isArray(rec.__mergedMembers) ? rec.__mergedMembers : [sn];
    const dimmed0 = _mapIsDimmed(sn, members0);
    const mk = L.marker([lat, lon], { icon: _svgDot(color, selected, dimmed0, style.shape) });
    mk.ttRec = rec; // used to list a cluster's samples in "list" mode
    mk.ttDimmed = dimmed0; // read by _clusterIcon to show the matching share
    mk.ttGroup = style.group; // read by _clusterIcon to build the group pie
    mk.on("click", () => {
      // Deselect any previously-selected dot, select this one, open its panel
      _selectedGroup = null;
      _selectedSample = sn;
      _refreshMapMarkerColors();
      _renderMapPanel(rec);
    });
    const members = Array.isArray(rec.__mergedMembers) ? rec.__mergedMembers : [sn];
    _markerObjects[sn] = { marker: mk, rec, color, lat, lon, members };
    if (_mapMarkerVisible(sn, members)) _markerLayer.addLayer(mk);
    bounds.push([lat, lon]);
  });
  // Drop picker entries for markers that no longer exist (e.g. after toggling
  // specimen merge, where the names change from libraries to specimens) so a
  // stale name can't keep something invisible with no way to switch it back on.
  const _live = new Set(Object.keys(_markerObjects));
  [..._mapHiddenSamples].forEach((n) => {
    if (!_live.has(n)) _mapHiddenSamples.delete(n);
  });
  const _picker = document.getElementById("map-sample-picker");
  if (_picker && _picker.open) _mapRenderSampleList();
  return bounds;
}

function _refreshMapMarkerColors() {
  if (!_markerObjects || !_leafletMap || !_markerLayer) return;
  _mapRefreshMatchSet();
  Object.entries(_markerObjects).forEach(([sn, obj]) => {
    if (!obj.marker) return;
    const members = Array.isArray(obj.members) && obj.members.length ? obj.members : [sn];
    if (!_mapMarkerVisible(sn, members)) {
      if (_markerLayer.hasLayer(obj.marker)) _markerLayer.removeLayer(obj.marker);
      return;
    }
    if (!_markerLayer.hasLayer(obj.marker)) _markerLayer.addLayer(obj.marker);
    const style = _mapGroupStyle(sn);
    const color = style.color;
    obj.color = color;
    obj.group = style.group;
    const dimmed = _mapIsDimmed(sn, members);
    obj.marker.ttDimmed = dimmed;
    obj.marker.ttGroup = style.group;
    obj.marker.setIcon(_svgDot(color, _selectedSample === sn, dimmed, style.shape));
  });
  // Keep the map's own group legend in step with the markers. onPick refits
  // the view after a toggle — the hide/show and shift-to-isolate behaviour is
  // handled by the legend itself, shared with the other grouped views.
  if (typeof window._mgRenderLegend === "function") {
    window._mgRenderLegend("map-group-legend", { onPick: _mgMapFitVisible });
  }
  // Recompute cluster bubbles (counts / membership may have changed)
  if (typeof _markerLayer.refreshClusters === "function") _markerLayer.refreshClusters();
  // The picker annotates rows with "not in filter"; keep that in step.
  const _pk = document.getElementById("map-sample-picker");
  if (_pk && _pk.open) _mapRenderSampleList();
}

/* Move the reusable map markup (#map-split, parked in the hidden #map-host)
   into the "Precise (lat/long)" view of the Mapping & Geography sub-tab. Runs
   once — after that the node already lives in the slot. */
function _ensureMapInPreciseWrap() {
  const slot = document.getElementById("geo-precise-slot");
  const split = document.getElementById("map-split");
  if (slot && split && split.parentElement !== slot) {
    slot.appendChild(split);
    if (_leafletMap) setTimeout(() => _leafletMap.invalidateSize(), 60);
  }
}

function _initMap() {
  _ensureMapInPreciseWrap();
  if (_leafletMap) {
    _leafletMap.invalidateSize();
    return;
  }
  _doInitMap();
}

/* Re-fit the map view to the current geo markers. Used after the
         container is resized (e.g. when the PDF print stylesheet shrinks the
         map) so the printed map matches the on-screen framing instead of
         keeping a stale zoom level. */
function _fitMapToData() {
  if (!_leafletMap) return;
  const bounds = [];
  RUN_META.forEach((r) => {
    if (r.latitude == null || r.longitude == null) return;
    const lat = parseFloat(r.latitude);
    const lon = parseFloat(r.longitude);
    if (!isNaN(lat) && !isNaN(lon)) bounds.push([lat, lon]);
  });
  if (bounds.length > 0) {
    _leafletMap.fitBounds(bounds, { padding: [40, 40], maxZoom: 9 });
  } else {
    _leafletMap.setView([20, 0], 2);
  }
}

function _doInitMap() {
  const container = document.getElementById("map-container");
  if (!container || _leafletMap) return;

  // Guard: Leaflet must be loaded (requires internet for CDN)
  if (typeof L === "undefined") {
    container.innerHTML =
      '<p style="padding:2em;text-align:center;color:#888">' +
      '<i class="fas fa-triangle-exclamation" style="font-size:2em;display:block;margin-bottom:.5em;opacity:.4"></i>' +
      "Map requires an internet connection to load the Leaflet library.<br>" +
      "<small>Open this report from a web server or ensure internet access.</small></p>";
    return;
  }

  const geoRows = RUN_META.filter((r) => r.latitude != null && r.longitude != null);
  if (geoRows.length === 0) return; // no geo data

  _leafletMap = L.map("map-container", { zoomControl: true });

  // CARTO Voyager — no Referer/policy restrictions
  L.tileLayer("https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png", {
    attribution:
      '© <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors © <a href="https://carto.com/attributions">CARTO</a>',
    subdomains: "abcd",
    maxZoom: 19,
  }).addTo(_leafletMap);

  _markerLayer = _makeMarkerLayer();
  _markerLayer.addTo(_leafletMap);

  // Cluster-click behaviour + its toggle control (only when clustering is on).
  if (typeof _markerLayer.on === "function" && typeof L.markerClusterGroup === "function") {
    _markerLayer.on("clusterclick", (a) => {
      const cluster = a.layer;
      if (_mapClusterMode === "list") {
        // Open the panel listing every sample in the cluster.
        const recs = cluster
          .getAllChildMarkers()
          .map((m) => m.ttRec)
          .filter(Boolean);
        if (recs.length) {
          _selectedSample = null;
          _renderMapGroupPanel(recs);
        }
      } else {
        // Expand: fan out if the dots share the same point, else zoom to fit.
        const lls = cluster.getAllChildMarkers().map((m) => m.getLatLng());
        const allSame = lls.length > 0 && lls.every((ll) => ll.equals(lls[0]));
        if (allSame) cluster.spiderfy();
        else cluster.zoomToBounds({ padding: [40, 40] });
      }
    });
    _addMapClusterControl();
  }

  // One colored dot per geo-located sample; overlapping dots auto-cluster
  // into a numbered bubble (see _clusterIcon / _makeMarkerLayer).
  const bounds = _addSampleMarkers();

  if (bounds.length > 0) {
    _leafletMap.fitBounds(bounds, { padding: [40, 40], maxZoom: 9 });
  } else {
    _leafletMap.setView([20, 0], 2);
  }
  // Force Leaflet to recalculate container size in case flex layout settled late
  setTimeout(() => {
    if (_leafletMap) _leafletMap.invalidateSize();
  }, 150);
}

/* ── Rebuild all map markers from current RUN_META ─────────────────── */
function _rebuildMapMarkers() {
  if (!_leafletMap || !_markerLayer) return;
  _markerLayer.clearLayers();
  _markerObjects = {};
  _selectedSample = null;
  _selectedGroup = null;
  closeMapPanel();

  const geoRows = RUN_META.filter((r) => r.latitude != null && r.longitude != null);
  if (geoRows.length === 0) return;

  _addSampleMarkers();

  setTimeout(() => _leafletMap.invalidateSize(), 50);
}

/* ── Map side panel ─────────────────────────────────────────────────── */

/* ── Swap the map-panel thead between single-sample and group modes ─── */
function _setMapPanelGroupMode(isGroup) {
  const thead = document.getElementById("map-panel-thead");
  const tbl = document.getElementById("map-panel-table");
  if (!thead) return;
  if (isGroup) {
    if (tbl) tbl.style.tableLayout = "auto";
    thead.innerHTML = `
            <th style="padding:5px 6px;text-align:left;white-space:nowrap;background:#f0f6ff" colspan="2">Sample / Organism</th>
            <th style="padding:5px 8px;text-align:right;white-space:nowrap;background:#f0f6ff">TASS&nbsp;%</th>
            <th style="padding:5px 8px;text-align:right;white-space:nowrap;background:#f0f6ff">#&nbsp;Reads</th>
            <th style="padding:5px 8px;text-align:right;white-space:nowrap;background:#f0f6ff">Cov&nbsp;%</th>
            <th style="padding:5px 6px;text-align:left;background:#f0f6ff">Category</th>`;
  } else {
    if (tbl) tbl.style.tableLayout = "fixed";
    thead.innerHTML = `
            <th id="mph-organism" data-col="Detected Organism" data-num="0" class="mph sortable-col"
                style="padding:5px 8px;text-align:left;cursor:pointer;user-select:none">
              Organism <span class="sort-arrow"></span>
            </th>
            <th id="mph-tass" data-col="TASS Score" data-num="1" class="mph sortable-col"
                style="padding:5px 8px;text-align:right;cursor:pointer;user-select:none">
              TASS&nbsp;% <span class="sort-arrow">▼</span>
            </th>
            <th id="mph-reads" data-col="% Reads" data-num="1" class="mph sortable-col"
                style="padding:5px 8px;text-align:right;cursor:pointer;user-select:none">
              %&nbsp;Reads <span class="sort-arrow"></span>
            </th>
            <th id="mph-cov" data-col="Coverage" data-num="1" class="mph sortable-col"
                style="padding:5px 8px;text-align:right;cursor:pointer;user-select:none">
              Cov&nbsp;% <span class="sort-arrow"></span>
            </th>
            <th id="mph-cat" data-col="Microbial Category" data-num="0" class="mph sortable-col"
                style="padding:5px 8px;text-align:left;cursor:pointer;user-select:none">
              Category <span class="sort-arrow"></span>
            </th>`;
    // Re-bind sort click handlers after rebuilding the thead
    document.querySelectorAll("#map-panel-thead .sortable-col").forEach((th) => {
      th.addEventListener("click", () => {
        const col = th.getAttribute("data-col");
        if (_panelSortCol === col) _panelSortAsc = !_panelSortAsc;
        else {
          _panelSortCol = col;
          _panelSortAsc = col !== "TASS Score";
        }
        _refreshMapPanelTable();
      });
    });
  }
}

function _renderMapPanel(rec) {
  const panel = document.getElementById("map-panel");
  if (!panel) return;

  // Switch to single-sample panel mode (only rebuild thead if coming from group mode)
  const _wasGroup = !!_selectedGroup;
  _selectedGroup = null;
  if (_wasGroup) _setMapPanelGroupMode(false);

  // Make panel and resize handle visible
  panel.style.display = "flex";
  const handle = document.getElementById("map-resize-handle");
  if (handle) handle.style.display = "block";

  // Header: sample name (truncate long names, full name on hover)
  const title = document.getElementById("map-panel-title");
  if (title) {
    const _full = rec.sample_name || "";
    title.textContent = _truncSampleName(_full, _sampleNameCap());
    title.title = _full;
  }

  // Show the "View in metadata table" button only if runmeta tab is available
  const viewBtn = document.getElementById("map-panel-view-btn");
  if (viewBtn) {
    const hasRunMeta = RUN_META.length > 0 && !Array.isArray(rec.__mergedMembers);
    viewBtn.style.display = hasRunMeta ? "inline-block" : "none";
  }

  // Meta info strip — show lat/lon first, then all other non-null fields dynamically
  const metaEl = document.getElementById("map-panel-meta");
  if (metaEl) {
    const skip = new Set(["sample_name", "latitude", "longitude", "__mergedMembers"]);
    const lines = [];

    // Lat/Lon always first (it's why the dot is on the map)
    if (rec.latitude != null && rec.longitude != null) {
      lines.push(`<b>Lat/Lon:</b> ${parseFloat(rec.latitude).toFixed(4)}, ${parseFloat(rec.longitude).toFixed(4)}`);
    }

    // All remaining non-null fields in order
    Object.entries(rec).forEach(([k, v]) => {
      if (skip.has(k) || v == null || v === "") return;
      lines.push(`<b>${_metaKeyLabel(k)}:</b> ${v}`);
    });

    metaEl.innerHTML = lines.length
      ? lines.join(" &nbsp;·&nbsp; ")
      : '<span style="color:#bbb">No additional metadata</span>';
  }

  _refreshMapPanelTable();
}

/* ── Group panel: single or multi-sample at same lat/lon ───────────── */
function _renderMapGroupPanel(recs) {
  if (!recs || recs.length === 0) return;
  if (recs.length === 1) {
    _renderMapPanel(recs[0]);
    return;
  }

  // Sort group by custom sample order so the map panel respects reordering
  _selectedGroup = _orderedSamples(recs.map((r) => r.sample_name))
    .map((sn) => recs.find((r) => r.sample_name === sn))
    .filter(Boolean);
  _selectedSample = null;
  _setMapPanelGroupMode(true);

  const panel = document.getElementById("map-panel");
  if (!panel) return;
  panel.style.display = "flex";
  const handle = document.getElementById("map-resize-handle");
  if (handle) handle.style.display = "block";

  // Title
  const title = document.getElementById("map-panel-title");
  if (title) title.textContent = `${recs.length} samples at this location`;

  // Hide "View in metadata table" button (ambiguous for group)
  const viewBtn = document.getElementById("map-panel-view-btn");
  if (viewBtn) viewBtn.style.display = "none";

  // Meta strip: lat/lon + color swatches for all samples
  const metaEl = document.getElementById("map-panel-meta");
  if (metaEl) {
    const r0 = recs[0];
    const coords =
      r0.latitude != null
        ? `<b>Lat/Lon:</b> ${parseFloat(r0.latitude).toFixed(4)}, ${parseFloat(r0.longitude).toFixed(4)} &nbsp;·&nbsp; `
        : "";
    // Keep the strip to a single compact line — just the location and the
    // sample count. Individual names appear as header rows in the table below.
    metaEl.innerHTML = `<div style="display:flex;align-items:center;flex-wrap:wrap">${coords}<b>${recs.length} samples</b></div>`;
  }

  // Populate the standard panel tbody with organisms for all samples
  _refreshMapGroupPanelTable();
}

/* ── Organism table for a group of co-located samples ────────────── */
function _refreshMapGroupPanelTable() {
  if (!_selectedGroup) return;
  // Re-sort to reflect any reordering in the main sample panel
  const _sgNames = _orderedSamples(_selectedGroup.map((r) => r.sample_name));
  _selectedGroup = _sgNames.map((sn) => _selectedGroup.find((r) => r.sample_name === sn)).filter(Boolean);
  const tbody = document.getElementById("map-panel-tbody");
  const empty = document.getElementById("map-panel-empty");
  const footer = document.getElementById("map-panel-footer");
  if (!tbody) return;

  const searchVal = (document.getElementById("map-panel-search")?.value || "").trim().toLowerCase();
  const sampleNames = new Set(_selectedGroup.map((r) => r.sample_name));

  // Respect the same filters the rest of the UI uses
  let allRows = filteredData().filter((r) => sampleNames.has(r["Specimen ID"]));
  if (searchVal) {
    allRows = allRows.filter(
      (r) =>
        (r["Detected Organism"] || "").toLowerCase().includes(searchVal) ||
        (r["Microbial Category"] || "").toLowerCase().includes(searchVal),
    );
  }

  if (allRows.length === 0) {
    tbody.innerHTML = "";
    if (empty) empty.style.display = "block";
    if (footer) footer.textContent = "No results";
    return;
  }
  if (empty) empty.style.display = "none";

  // Group rows by sample, preserving the order from _selectedGroup
  const grouped = {};
  _selectedGroup.forEach((r) => {
    grouped[r.sample_name] = [];
  });
  allRows.forEach((r) => {
    const sn = r["Specimen ID"];
    if (grouped[sn]) grouped[sn].push(r);
    else grouped[sn] = [r];
  });
  // Sort each sample's organisms by TASS desc
  Object.values(grouped).forEach((arr) =>
    arr.sort((a, b) => parseFloat(b["TASS Score"] || 0) - parseFloat(a["TASS Score"] || 0)),
  );

  let html = "";
  _selectedGroup.forEach((meta) => {
    const sn = meta.sample_name;
    const orgRows = grouped[sn] || [];
    if (orgRows.length === 0) return;
    const col = sampleColors[sn] || "#1565C0";
    const swatch = `<span title="${sn}" style="display:inline-block;width:12px;height:12px;
            border-radius:2px;background:${col};margin-right:5px;vertical-align:middle;
            box-shadow:0 0 0 1px rgba(0,0,0,.18);flex-shrink:0"></span>`;

    // ── Sample header row (truncate long names, full name on hover) ──
    const _snShort = _truncSampleName(sn, _sampleNameCap());
    const _snTitle = String(sn).replace(/"/g, "&quot;");
    html += `<tr style="background:#e8f0fe;border-top:2px solid ${col}">
            <td colspan="6" style="padding:5px 8px;font-weight:700;font-size:0.88em;color:#1a237e">
              ${swatch}<span title="${_snTitle}">${_snShort}</span>
              <span style="font-weight:400;color:#666;margin-left:6px;font-size:0.9em">
                ${orgRows.length} organism${orgRows.length !== 1 ? "s" : ""}
              </span>
            </td>
          </tr>`;

    // ── Organism rows ────────────────────────────────────────────────
    orgRows.forEach((r, i) => {
      const bg = i % 2 === 0 ? "#fafbff" : "#fff";
      const cat = r["Microbial Category"] || "Unknown";
      const catColor = _CAT_COLORS[cat] || "#555";
      const tass = parseFloat(r["TASS Score"] || 0);
      const reads = parseFloat(r["# Reads Aligned"] || r["% Reads"] || 0);
      const cov = parseFloat(r["Coverage"] || 0);
      const org = r["Detected Organism"] || "";
      const bar = `<div style="display:inline-block;width:${Math.round(tass * 0.38)}px;
              height:4px;background:${catColor};border-radius:2px;opacity:.5;
              margin-left:3px;vertical-align:middle"></div>`;
      html += `<tr style="background:${bg}">
              <td colspan="2" style="padding:3px 8px 3px 20px;border-bottom:1px solid #eef;
                max-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap"
                title="${org.replace(/"/g, "&quot;")}">${org}</td>
              <td style="padding:3px 8px;border-bottom:1px solid #eef;text-align:right;white-space:nowrap">
                ${tass.toFixed(1)}${bar}</td>
              <td style="padding:3px 8px;border-bottom:1px solid #eef;text-align:right;white-space:nowrap">
                ${reads >= 1000 ? reads.toLocaleString() : reads.toFixed(reads < 1 ? 2 : 0)}</td>
              <td style="padding:3px 8px;border-bottom:1px solid #eef;text-align:right;white-space:nowrap">
                ${cov.toFixed(1)}</td>
              <td style="padding:3px 8px;border-bottom:1px solid #eef">
                <span style="background:${catColor}22;color:${catColor};font-size:0.72em;font-weight:700;
                  padding:1px 5px;border-radius:4px">${cat}</span>
              </td>
            </tr>`;
    });
  });

  tbody.innerHTML = html;

  const total = allRows.length;
  if (footer) footer.textContent = `${total} organism${total !== 1 ? "s" : ""} across ${_selectedGroup.length} samples`;
}

function _renameSample(oldName, newName) {
  if (!oldName || !newName || oldName === newName) return;

  // Update DATA rows
  DATA.forEach((r) => {
    if (r["Specimen ID"] === oldName) r["Specimen ID"] = newName;
  });

  // Update protein rows
  [PROT.genus_summary, PROT.per_gene_hits, PROT.sample_overview, PROT.amr_genes].forEach((arr) => {
    if (!Array.isArray(arr)) return;
    arr.forEach((r) => {
      if (r["Specimen ID"] === oldName) r["Specimen ID"] = newName;
    });
  });

  // Update run metadata
  RUN_META.forEach((r) => {
    if (r.sample_name === oldName) r.sample_name = newName;
  });

  // Update CONTIG_DATA — drives the histogram tab. Without this the
  // histogram filter (visibleTaxa.has(`${cd.sample}||${cd.taxon_id}`))
  // misses the renamed sample because cd.sample still has the old name
  // while DATA's Specimen ID now has the new name.
  CONTIG_DATA.forEach((cd) => {
    if (cd && cd.sample === oldName) cd.sample = newName;
  });
  if (typeof _invalidateSummaryHistMap === "function") _invalidateSummaryHistMap();

  // Mirror the rename into BOOT so a later "Clear uploaded data" or
  // any code path that re-reads BOOT doesn't bring the old name back.
  if (BOOT) {
    if (Array.isArray(BOOT.records)) {
      BOOT.records.forEach((r) => {
        if (r["Specimen ID"] === oldName) r["Specimen ID"] = newName;
      });
    }
    if (Array.isArray(BOOT.contig_data)) {
      BOOT.contig_data.forEach((cd) => {
        if (cd && cd.sample === oldName) cd.sample = newName;
      });
    }
    if (BOOT.prot_data && typeof BOOT.prot_data === "object") {
      Object.keys(BOOT.prot_data).forEach((k) => {
        if (!Array.isArray(BOOT.prot_data[k])) return;
        BOOT.prot_data[k].forEach((r) => {
          if (r && r["Specimen ID"] === oldName) r["Specimen ID"] = newName;
        });
      });
    }
  }

  // Mirror the rename into the upload buffers so that the next file
  // drop (which calls _mergeAndRedraw → DATA = BOOT + _uploadedRows)
  // doesn't resurrect the OLD name from the buffers.
  if (window._renameInUploadBuffers) window._renameInUploadBuffers(oldName, newName);

  // Update per-sample UI state
  if (sampleColors[oldName] !== undefined) sampleColors[newName] = sampleColors[oldName];
  if (sampleHidden[oldName] !== undefined) sampleHidden[newName] = sampleHidden[oldName];
  if (sampleRescale[oldName] !== undefined) sampleRescale[newName] = sampleRescale[oldName];
  delete sampleColors[oldName];
  delete sampleHidden[oldName];
  delete sampleRescale[oldName];

  // Keep custom sample order in sync
  if (Array.isArray(_sampleOrder)) {
    _sampleOrder = _sampleOrder.map((s) => (s === oldName ? newName : s));
  }

  // Update map markers index
  if (_markerObjects[oldName]) {
    const obj = _markerObjects[oldName];
    obj.rec.sample_name = newName;
    _markerObjects[newName] = obj;
    delete _markerObjects[oldName];
  }

  if (_selectedSample === oldName) _selectedSample = newName;

  buildSampleList();
  _computeBslLevels();
  _buildRunMetaTable();
  if (typeof _updateMetaSubTabStates === "function") _updateMetaSubTabStates();
  // Rebuild histogram selectors so the rename is visible immediately
  if (window._resetHistSelectors) window._resetHistSelectors();
  redraw();

  const rec = RUN_META.find((r) => r.sample_name === newName);
  if (rec) _renderMapPanel(rec);
}

function editSampleName() {
  if (!_selectedSample) return;
  const current = _selectedSample;
  const next = prompt("Rename sample", current);
  if (next == null) return;
  const trimmed = String(next).trim();
  if (!trimmed) return;
  if (trimmed === current) return;

  const exists = DATA.some((r) => r["Specimen ID"] === trimmed) || RUN_META.some((r) => r.sample_name === trimmed);
  if (exists && !confirm(`Sample "${trimmed}" already exists. Merge under this name?`)) return;

  _renameSample(current, trimmed);
}

function editSampleNameFromSidebar(current) {
  if (!current) return;
  const next = prompt("Rename sample", current);
  if (next == null) return;
  const trimmed = String(next).trim();
  if (!trimmed) return;
  if (trimmed === current) return;

  const exists = DATA.some((r) => r["Specimen ID"] === trimmed) || RUN_META.some((r) => r.sample_name === trimmed);
  if (exists && !confirm(`Sample "${trimmed}" already exists. Merge under this name?`)) return;

  _renameSample(current, trimmed);
}

// Wire up sortable column headers (called once after init)
function _initPanelSortHeaders() {
  document.querySelectorAll("th.sortable-col").forEach((th) => {
    th.addEventListener("click", () => {
      const col = th.getAttribute("data-col");
      if (_panelSortCol === col) {
        _panelSortAsc = !_panelSortAsc;
      } else {
        _panelSortCol = col;
        _panelSortAsc = col === "Detected Organism" || col === "Microbial Category";
      }
      _refreshMapPanelTable();
    });
  });
}

function _refreshMapPanelTable() {
  if (_selectedGroup) {
    _refreshMapGroupPanelTable();
    return;
  }
  if (!_selectedSample) return;
  const tbody = document.getElementById("map-panel-tbody");
  const empty = document.getElementById("map-panel-empty");
  const footer = document.getElementById("map-panel-footer");
  if (!tbody) return;

  // Search filter text
  const searchVal = (document.getElementById("map-panel-search")?.value || "").trim().toLowerCase();

  // Apply same filters as filteredData() but restrict to the selected sample
  let rows = filteredData().filter((r) => r["Specimen ID"] === _selectedSample);

  // Additional organism name search
  if (searchVal) {
    rows = rows.filter(
      (r) =>
        (r["Detected Organism"] || "").toLowerCase().includes(searchVal) ||
        (r["Microbial Category"] || "").toLowerCase().includes(searchVal),
    );
  }

  // Sort by _panelSortCol / _panelSortAsc
  const isNum = !["Detected Organism", "Microbial Category"].includes(_panelSortCol);
  rows = rows.slice().sort((a, b) => {
    const av = isNum ? parseFloat(a[_panelSortCol] || 0) : a[_panelSortCol] || "";
    const bv = isNum ? parseFloat(b[_panelSortCol] || 0) : b[_panelSortCol] || "";
    const cmp = isNum ? av - bv : av.localeCompare(bv);
    return _panelSortAsc ? cmp : -cmp;
  });

  // Update sort arrow indicators
  document.querySelectorAll("th.sortable-col").forEach((th) => {
    const arrowEl = th.querySelector(".sort-arrow");
    if (!arrowEl) return;
    if (th.getAttribute("data-col") === _panelSortCol) {
      th.style.color = "#0d47a1";
      arrowEl.textContent = _panelSortAsc ? "▲" : "▼";
    } else {
      th.style.color = "";
      arrowEl.textContent = "";
    }
  });

  if (rows.length === 0) {
    tbody.innerHTML = "";
    if (empty) empty.style.display = "block";
    if (footer) footer.textContent = "No results";
    return;
  }
  if (empty) empty.style.display = "none";

  tbody.innerHTML = rows
    .map((r, i) => {
      const bg = i % 2 === 0 ? "#fafbff" : "#fff";
      const cat = r["Microbial Category"] || "Unknown";
      const catColor = _CAT_COLORS[cat] || "#555";
      const tass = parseFloat(r["TASS Score"] || 0);
      const bar = `<div style="display:inline-block;width:${Math.round(
        tass * 0.5,
      )}px;height:5px;background:${catColor};border-radius:3px;opacity:.6;margin-left:4px;vertical-align:middle"></div>`;
      const org = r["Detected Organism"] || "";
      return `<tr style="background:${bg}">
            <td style="padding:4px 8px;border-bottom:1px solid #eef" title="${org.replace(/"/g, "&quot;")}">${org}</td>
            <td style="padding:4px 8px;border-bottom:1px solid #eef;text-align:right">${tass.toFixed(1)}${bar}</td>
            <td style="padding:4px 8px;border-bottom:1px solid #eef;text-align:right">${parseFloat(
              r["% Reads"] || 0,
            ).toFixed(2)}</td>
            <td style="padding:4px 8px;border-bottom:1px solid #eef;text-align:right">${parseFloat(
              r["Coverage"] || 0,
            ).toFixed(1)}</td>
            <td style="padding:4px 8px;border-bottom:1px solid #eef"><span style="background:${catColor}22;color:${catColor};font-size:0.72em;font-weight:700;padding:1px 5px;border-radius:4px">${cat}</span></td>
          </tr>`;
    })
    .join("");

  const sortLabel =
    _panelSortCol === "TASS Score"
      ? "TASS"
      : _panelSortCol === "% Reads"
      ? "% Reads"
      : _panelSortCol === "Coverage"
      ? "Coverage"
      : _panelSortCol === "Microbial Category"
      ? "Category"
      : "Organism";
  const dirLabel = _panelSortAsc ? "↑" : "↓";
  if (footer)
    footer.textContent = `${rows.length} organism${rows.length !== 1 ? "s" : ""} · sorted by ${sortLabel} ${dirLabel}`;
}

function closeMapPanel() {
  const panel = document.getElementById("map-panel");
  const handle = document.getElementById("map-resize-handle");
  if (panel) panel.style.display = "none";
  if (handle) handle.style.display = "none";
  // Clear search bar
  const search = document.getElementById("map-panel-search");
  if (search) search.value = "";
  // Deselect single-sample marker
  if (_selectedSample && _markerObjects[_selectedSample]) {
    const obj = _markerObjects[_selectedSample];
    obj.marker.setIcon(_svgDot(obj.color, false));
  }
  _selectedSample = null;
  // Deselect group marker (reset all pie icons to non-selected)
  if (_selectedGroup) {
    _selectedGroup = null;
    _refreshMapMarkerColors();
    // Restore thead to single-sample layout for next open
    _setMapPanelGroupMode(false);
  }
  // Resize map to reclaim space
  if (_leafletMap) setTimeout(() => _leafletMap.invalidateSize(), 50);
}

// "View" button: jump to Run Metadata tab and highlight this sample's row
function viewSampleInMetaTab() {
  if (!_selectedSample) return;
  _runmetaHighlightSample = _selectedSample;

  // Switch to the runmeta tab programmatically
  const rmBtn = document.getElementById("runmeta-tab-btn");
  if (rmBtn) rmBtn.click();

  // Rebuild table with highlight (click triggers _buildRunMetaTable via tab handler)
  // But also call it directly in case the tab was already active
  _buildRunMetaTable();
}

/* ═══════════════════════════════════════════════════════════════════════════
   -  §  MAP SAMPLE PICKER
   -     The collapsible "Samples on map" panel: one checkbox per mapped
   -     sample, a search box to find one in a long list, and bulk actions.
   -     Everything here writes only to _mapHiddenSamples, so nothing it does
   -     leaks into the other tabs.
   ═══════════════════════════════════════════════════════════════════════════ */

/** Marker names currently on the map, in the order they are drawn. Derived
 *  from _markerObjects so it automatically follows the specimen-merge state
 *  (specimen names when merged, raw sample ids otherwise). */
function _mapPickerNames() {
  return Object.keys(_markerObjects || {}).sort((a, b) => a.localeCompare(b, undefined, { numeric: true }));
}

/** Re-render the checkbox list. `query` narrows which rows are LISTED; it does
 *  not hide anything on the map — that is what the checkboxes are for. */
function _mapRenderSampleList(query) {
  const list = document.getElementById("map-sample-list");
  if (!list) return;
  const q = String(query == null ? (document.getElementById("map-sample-search") || {}).value || "" : query)
    .trim()
    .toLowerCase();
  const names = _mapPickerNames();
  const shown = q ? names.filter((n) => n.toLowerCase().includes(q)) : names;

  if (!names.length) {
    list.innerHTML = '<div style="padding:6px 8px;color:#8a97a4;font-size:.8em">No mapped samples yet.</div>';
  } else if (!shown.length) {
    list.innerHTML = `<div style="padding:6px 8px;color:#8a97a4;font-size:.8em">No sample matches “${q}”.</div>`;
  } else {
    list.innerHTML = shown
      .map((n) => {
        const hidden = _mapLocallyHidden(n);
        const obj = _markerObjects[n] || {};
        const dot = obj.color || sampleColors[n] || "#1565C0";
        // Flag samples the sidebar filters exclude, so "hide everything that
        // isn't a hit" is an informed click rather than a guess.
        const off = obj.marker && obj.marker.ttDimmed ? " · not in filter" : "";
        const globallyHidden = (obj.members || [n]).every((s) => sampleHidden[s]);
        return (
          `<label class="map-sample-row" title="${n}${globallyHidden ? " — hidden globally in the sidebar" : ""}" ` +
          `style="display:flex;align-items:center;gap:7px;padding:3px 8px;cursor:pointer;` +
          `${globallyHidden ? "opacity:.45;" : ""}">` +
          `<input type="checkbox" data-map-sample="${String(n).replace(/"/g, "&quot;")}"${hidden ? "" : " checked"}${
            globallyHidden ? " disabled" : ""
          } />` +
          `<span style="width:9px;height:9px;border-radius:50%;background:${dot};flex:0 0 auto;` +
          `border:1px solid rgba(0,0,0,.2)"></span>` +
          `<span style="overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-size:.8em">${n}</span>` +
          `<span style="margin-left:auto;color:#9aa6b4;font-size:.72em;white-space:nowrap">${off}</span>` +
          `</label>`
        );
      })
      .join("");
  }

  const count = document.getElementById("map-sample-count");
  if (count) {
    const visible = names.filter((n) => !_mapLocallyHidden(n)).length;
    count.textContent = names.length ? `${visible} of ${names.length} shown` : "";
  }
}

/** Apply a bulk action to the CURRENTLY LISTED samples (so a search term also
 *  scopes the bulk action — "type RSV, click Only these" is the fast path). */
function _mapBulkSamples(action) {
  const q = ((document.getElementById("map-sample-search") || {}).value || "").trim().toLowerCase();
  const names = _mapPickerNames();
  const scoped = q ? names.filter((n) => n.toLowerCase().includes(q)) : names;
  const scopedSet = new Set(scoped);

  if (action === "all") {
    scoped.forEach((n) => _mapHiddenSamples.delete(n));
  } else if (action === "none") {
    scoped.forEach((n) => _mapHiddenSamples.add(n));
  } else if (action === "invert") {
    scoped.forEach((n) => (_mapHiddenSamples.has(n) ? _mapHiddenSamples.delete(n) : _mapHiddenSamples.add(n)));
  } else if (action === "only") {
    // Keep only what the search listed — everything else off the map.
    names.forEach((n) => (scopedSet.has(n) ? _mapHiddenSamples.delete(n) : _mapHiddenSamples.add(n)));
  } else if (action === "matches") {
    // Keep only samples the sidebar filters currently match.
    const matched = typeof _ttFilterMatchedKeys === "function" ? _ttFilterMatchedKeys() : null;
    names.forEach((n) => {
      const obj = _markerObjects[n] || {};
      const members = obj.members && obj.members.length ? obj.members : [n];
      const hit =
        matched && (_ttSampleMatchesFilter(n, matched) || members.some((m) => _ttSampleMatchesFilter(m, matched)));
      if (hit) _mapHiddenSamples.delete(n);
      else _mapHiddenSamples.add(n);
    });
  }
  _refreshMapMarkerColors();
  _mapRenderSampleList();
}

/** Wire the picker once. Uses event delegation for the checkboxes so the list
 *  can be re-rendered freely without rebinding. */
function _wireMapSamplePicker() {
  const panel = document.getElementById("map-sample-picker");
  if (!panel || panel._wired) return;
  panel._wired = true;

  const list = document.getElementById("map-sample-list");
  if (list)
    list.addEventListener("change", (e) => {
      const cb = e.target.closest("input[data-map-sample]");
      if (!cb) return;
      const n = cb.getAttribute("data-map-sample");
      if (cb.checked) _mapHiddenSamples.delete(n);
      else _mapHiddenSamples.add(n);
      _refreshMapMarkerColors();
      const count = document.getElementById("map-sample-count");
      if (count) {
        const names = _mapPickerNames();
        count.textContent = `${names.filter((x) => !_mapLocallyHidden(x)).length} of ${names.length} shown`;
      }
    });

  const search = document.getElementById("map-sample-search");
  if (search) search.addEventListener("input", () => _mapRenderSampleList());

  panel.querySelectorAll("[data-map-bulk]").forEach((btn) =>
    btn.addEventListener("click", (e) => {
      e.preventDefault();
      _mapBulkSamples(btn.getAttribute("data-map-bulk"));
    }),
  );

  // Populate lazily: building 160+ rows is wasted work until it's opened.
  panel.addEventListener("toggle", () => {
    if (panel.open) _mapRenderSampleList();
  });
}
