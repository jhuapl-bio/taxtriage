/* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB: MAP                (data-tab="map"  —  hidden if no lat/lon)
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
function _svgDot(color, selected) {
  const r = selected ? 11 : 8;
  const ring = selected
    ? `<circle cx="${r + 5}" cy="${r + 5}" r="${r + 3}" fill="none" stroke="${color}" stroke-width="2" opacity=".45"/>`
    : "";
  return L.divIcon({
    className: "",
    html: `<svg xmlns="http://www.w3.org/2000/svg" width="${(r + 5) * 2}" height="${(r + 5) * 2}">
            <circle cx="${r + 5}" cy="${r + 5}" r="${r}" fill="${color}"
              stroke="rgba(255,255,255,.9)" stroke-width="${selected ? 2.5 : 1.5}"
              filter="drop-shadow(0 1px 3px rgba(0,0,0,.38))"/>
            ${ring}
          </svg>`,
    iconSize: [(r + 5) * 2, (r + 5) * 2],
    iconAnchor: [r + 5, r + 5],
  });
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

function _refreshMapMarkerColors() {
  if (!_markerObjects || !_leafletMap || !_markerLayer) return;
  const _groupNames = _selectedGroup ? new Set(_selectedGroup.map((r) => r.sample_name)) : null;

  // Track which locKeys we've already processed (single-sample entries are indexed
  // both by locKey and by sample_name — skip the duplicate sample_name entry)
  const processed = new Set();

  Object.entries(_markerObjects).forEach(([key, obj]) => {
    if (!obj.marker || !obj.recs) return;
    // Skip duplicate sample_name aliases pointing to same obj as a locKey entry
    const locKey =
      obj.recs.length === 1
        ? `${parseFloat(obj.recs[0].latitude).toFixed(3)}_${parseFloat(obj.recs[0].longitude).toFixed(3)}`
        : null;
    if (locKey && locKey !== key && processed.has(locKey)) return;
    if (locKey) processed.add(locKey);

    const visibleRecs = obj.recs.filter((r) => !sampleHidden[r.sample_name]);
    const allHidden = visibleRecs.length === 0;

    if (allHidden) {
      // Remove marker from layer if present
      if (_markerLayer.hasLayer(obj.marker)) _markerLayer.removeLayer(obj.marker);
      return;
    }

    // Re-add to layer if it was previously removed
    if (!_markerLayer.hasLayer(obj.marker)) obj.marker.addTo(_markerLayer);

    const isSingle = visibleRecs.length === 1;
    let selected;
    if (obj.recs.length > 1) {
      selected = _groupNames
        ? visibleRecs.some((r) => _groupNames.has(r.sample_name))
        : visibleRecs.some((r) => r.sample_name === _selectedSample);
      const cols = visibleRecs.map((r) => sampleColors[r.sample_name] || "#1565C0");
      obj.marker.setIcon(isSingle ? _svgDot(cols[0], selected) : _pieSvg(cols, selected));
    } else {
      selected = _selectedSample === key || _selectedSample === obj.recs[0]?.sample_name;
      const color = sampleColors[key] || obj.color || "#1565C0";
      obj.color = color;
      obj.marker.setIcon(_svgDot(color, selected));
    }
  });
}

function _initMap() {
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

  _markerLayer = L.layerGroup().addTo(_leafletMap);

  // Assign colours (prefer sidebar sample colors)
  const sampleNames = [...new Set(geoRows.map((r) => r.sample_name))];
  sampleNames.forEach((n, i) => {
    if (!sampleColors[n]) sampleColors[n] = PALETTE[i % PALETTE.length];
  });

  // ── Group rows by lat/lon (rounded to 3 dp) for pie-chart markers ──
  const _locGroups = {};
  geoRows.forEach((rec) => {
    const lat = parseFloat(rec.latitude);
    const lon = parseFloat(rec.longitude);
    if (isNaN(lat) || isNaN(lon)) return;
    const key = `${lat.toFixed(3)}_${lon.toFixed(3)}`;
    if (!_locGroups[key]) _locGroups[key] = { lat, lon, recs: [] };
    _locGroups[key].recs.push(rec);
  });

  const bounds = [];
  Object.entries(_locGroups).forEach(([locKey, grp]) => {
    const { lat, lon, recs } = grp;
    const isSingle = recs.length === 1;
    const colors = recs.map((r) => sampleColors[r.sample_name] || "#1565C0");
    const icon = isSingle ? _svgDot(colors[0], false) : _pieSvg(colors, false);
    const mk = L.marker([lat, lon], { icon });

    mk.on("click", () => {
      // Deselect previous group marker
      Object.values(_markerObjects).forEach((obj) => {
        if (!obj.marker) return;
        const sel = false;
        if (obj.recs && obj.recs.length > 1) {
          const cols = obj.recs.map((r) => sampleColors[r.sample_name] || "#1565C0");
          obj.marker.setIcon(_pieSvg(cols, sel));
        } else {
          obj.marker.setIcon(_svgDot(obj.color || "#1565C0", sel));
        }
      });
      _selectedSample = recs[0].sample_name;
      mk.setIcon(isSingle ? _svgDot(colors[0], true) : _pieSvg(colors, true));
      _renderMapGroupPanel(recs);
    });

    mk.addTo(_markerLayer);
    _markerObjects[locKey] = { marker: mk, recs, color: colors[0] };
    // Also index by sample_name for single markers (backward compat)
    if (isSingle) _markerObjects[recs[0].sample_name] = _markerObjects[locKey];
    bounds.push([lat, lon]);
  });

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

  const _lg = {};
  geoRows.forEach((rec) => {
    const lat = parseFloat(rec.latitude),
      lon = parseFloat(rec.longitude);
    if (isNaN(lat) || isNaN(lon)) return;
    const key = `${lat.toFixed(3)}_${lon.toFixed(3)}`;
    if (!_lg[key]) _lg[key] = { lat, lon, recs: [] };
    _lg[key].recs.push(rec);
  });

  Object.entries(_lg).forEach(([locKey, grp]) => {
    const { lat, lon, recs } = grp;
    const isSingle = recs.length === 1;
    const colors = recs.map((r) => sampleColors[r.sample_name] || "#1565C0");
    const mk = L.marker([lat, lon], { icon: isSingle ? _svgDot(colors[0], false) : _pieSvg(colors, false) });
    mk.on("click", () => {
      Object.values(_markerObjects).forEach((obj) => {
        if (!obj.marker) return;
        if (obj.recs && obj.recs.length > 1)
          obj.marker.setIcon(
            _pieSvg(
              obj.recs.map((r) => sampleColors[r.sample_name] || "#1565C0"),
              false,
            ),
          );
        else obj.marker.setIcon(_svgDot(obj.color || "#1565C0", false));
      });
      _selectedSample = recs[0].sample_name;
      mk.setIcon(isSingle ? _svgDot(colors[0], true) : _pieSvg(colors, true));
      _renderMapGroupPanel(recs);
    });
    mk.addTo(_markerLayer);
    _markerObjects[locKey] = { marker: mk, recs, color: colors[0] };
    if (isSingle) _markerObjects[recs[0].sample_name] = _markerObjects[locKey];
  });

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

  // Header: sample name
  const title = document.getElementById("map-panel-title");
  if (title) title.textContent = rec.sample_name;

  // Show the "View in metadata table" button only if runmeta tab is available
  const viewBtn = document.getElementById("map-panel-view-btn");
  if (viewBtn) {
    const hasRunMeta = RUN_META.length > 0;
    viewBtn.style.display = hasRunMeta ? "inline-block" : "none";
  }

  // Meta info strip — show lat/lon first, then all other non-null fields dynamically
  const metaEl = document.getElementById("map-panel-meta");
  if (metaEl) {
    const skip = new Set(["sample_name", "latitude", "longitude"]);
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
    const swatches = recs
      .map((r) => {
        const col = sampleColors[r.sample_name] || "#1565C0";
        return `<span title="${r.sample_name}" style="display:inline-block;width:11px;height:11px;
              border-radius:2px;background:${col};margin-right:3px;vertical-align:middle;
              box-shadow:0 0 0 1px rgba(0,0,0,.2)"></span>`;
      })
      .join("");
    metaEl.innerHTML =
      coords +
      swatches +
      `<span style="color:#555;font-size:0.9em;margin-left:2px">${recs.map((r) => r.sample_name).join(", ")}</span>`;
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

    // ── Sample header row ────────────────────────────────────────────
    html += `<tr style="background:#e8f0fe;border-top:2px solid ${col}">
            <td colspan="6" style="padding:5px 8px;font-weight:700;font-size:0.88em;color:#1a237e">
              ${swatch}${sn}
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
