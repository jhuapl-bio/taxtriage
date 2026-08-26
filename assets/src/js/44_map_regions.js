/* ═══════════════════════════════════════════════════════════════════════════
       -  §  MAP REGIONS  —  draw an area, get a report for it
       -
       -     Optional, off until the user arms a tool. Two ways to cut out an
       -     area of the precise (lat/long) map:
       -
       -       lasso      hold the mouse down and trace a shape freehand
       -       rectangle  click and drag a box
       -
       -     Every region the user draws is kept in `_ttRegions` — the shape,
       -     its name and its colour — INDEPENDENTLY of whether it is currently
       -     drawn on the map. That separation is the whole point of the "Show
       -     regions" switch: flicking it off strips the overlays so the user
       -     can see the pins underneath, and flicking it back on restores
       -     every region exactly as it was, with its summaries intact.
       -
       -     What a region produces:
       -       • a live summary card — samples, detections, organisms, reads,
       -         mean TASS, and the organisms that dominate inside it
       -       • a cross-region comparison table once two or more exist,
       -         including shared vs. region-only organisms
       -       • optionally a `map_region` metadata column ("Group by region"),
       -         which makes the drawn regions a first-class grouping dimension
       -         for the shared Group by bar, the group heatmap, the group
       -         network and the cross-entry comparison
       -
       -     Membership is filter-aware: a region's stats describe the samples
       -     inside it that survive the current sidebar filters, so narrowing
       -     the run narrows every region report with it.
       -
       -     Public surface:
       -       ttMapRegions.list() / .capture() / .restore(state)
       -       _ttRegionRefresh()      — recompute + re-render (filters changed)
       -       _ttRegionOnMapShown()   — the Mapping tab just became visible
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  "use strict";

  /* ── state ──────────────────────────────────────────────────────────── */

  // {id, name, type:"lasso"|"rect", pts:[[lat,lng],…], color}
  let _regions = [];
  let _seq = 0;
  let _visible = true; // the "Show regions" switch
  let _armed = null; // null | "lasso" | "rect" — tool waiting for a drag
  let _asGrouping = false; // "Group by region" — publishes map_region
  const _layers = new Map(); // region id → Leaflet layer (only while visible)
  let _draft = null; // the shape being dragged right now
  let _abortDraft = null; // set once the map is wired — drops a half-drawn shape
  let _floatOpen = true; // the on-map card stack: expanded or collapsed to its bar
  let _chartRegion = null; // region id whose chart overlay is open, or null
  let _chartView = "bars"; // "bars" | "heat"
  let _wired = false;

  // Region colours are deliberately NOT the sample palette: a region outline
  // has to stay legible on top of coloured pins, so this is a short, high
  // contrast set that reads as "annotation", not "data".
  const REGION_COLORS = ["#e65100", "#6a1b9a", "#00695c", "#c2185b", "#1565c0", "#4e342e", "#33691e", "#ad1457"];

  const _esc = (s) =>
    String(s == null ? "" : s).replace(
      /[&<>"']/g,
      (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c],
    );
  const _num = (v) => {
    const n = parseFloat(v);
    return isNaN(n) ? 0 : n;
  };
  const _fmt = (n) =>
    n >= 1e9
      ? (n / 1e9).toFixed(1) + "B"
      : n >= 1e6
      ? (n / 1e6).toFixed(1) + "M"
      : n >= 1e3
      ? (n / 1e3).toFixed(1) + "k"
      : String(Math.round(n));

  /* ── geometry ───────────────────────────────────────────────────────── */

  /* Standard ray-casting point-in-polygon on [lat, lng] pairs. Longitudes are
     used as given: a run whose samples straddle the antimeridian would need a
     normalisation pass, but every lasso the user can physically draw on a
     single map view stays on one side of it. */
  function _inPoly(lat, lng, pts) {
    let inside = false;
    for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
      const yi = pts[i][0],
        xi = pts[i][1],
        yj = pts[j][0],
        xj = pts[j][1];
      const hit = yi > lat !== yj > lat && lng < ((xj - xi) * (lat - yi)) / (yj - yi || 1e-12) + xi;
      if (hit) inside = !inside;
    }
    return inside;
  }

  /* ── sample membership ──────────────────────────────────────────────── */

  /* Every geo-located sample currently ON the map, as {name, lat, lon}. Reads
     _markerObjects when the map has been built (so specimen merging and the
     map's own show/hide picker are already accounted for) and falls back to
     RUN_META when it hasn't. */
  function _geoSamples() {
    const out = [];
    if (typeof _markerObjects !== "undefined" && _markerObjects && Object.keys(_markerObjects).length) {
      Object.keys(_markerObjects).forEach((sn) => {
        const o = _markerObjects[sn];
        if (!o || !isFinite(o.lat) || !isFinite(o.lon)) return;
        out.push({ name: sn, lat: o.lat, lon: o.lon, members: o.members || [sn] });
      });
      return out;
    }
    (typeof RUN_META !== "undefined" ? RUN_META : []).forEach((r) => {
      const lat = parseFloat(r.latitude),
        lon = parseFloat(r.longitude);
      if (isNaN(lat) || isNaN(lon)) return;
      out.push({ name: r.sample_name, lat, lon, members: [r.sample_name] });
    });
    return out;
  }

  function _samplesIn(region) {
    return _geoSamples()
      .filter((s) => _inPoly(s.lat, s.lon, region.pts))
      .map((s) => s.name);
  }

  /* Detection rows for a set of samples, under the active sidebar filters.
     Merged specimens are matched through their member libraries so a region
     over a merged pin still finds its rows. */
  function _rowsFor(names) {
    const want = new Set(names);
    if (typeof specimenOf === "function" && typeof specimenMergeEnabled !== "undefined" && specimenMergeEnabled) {
      names.forEach((n) => want.add(n));
    }
    const fd = typeof filteredData === "function" ? filteredData() : [];
    return fd.filter((r) => {
      const sid = r["Specimen ID"];
      if (want.has(sid)) return true;
      if (typeof specimenOf === "function") {
        try {
          return want.has(specimenOf(sid));
        } catch (e) {}
      }
      return false;
    });
  }

  /* Everything a region card / comparison row needs, computed once. */
  function _stats(region) {
    const samples = _samplesIn(region);
    const rows = _rowsFor(samples);
    const byOrg = new Map();
    let reads = 0,
      tassSum = 0,
      tassN = 0,
      pass = 0;
    rows.forEach((r) => {
      const org = r["Detected Organism"] || "(unnamed)";
      const rd = _num(r["# Reads Aligned"]);
      const ts = _num(r["TASS Score"]);
      reads += rd;
      if (r["TASS Score"] != null && r["TASS Score"] !== "") {
        tassSum += ts;
        tassN++;
      }
      if (typeof isTruthy === "function" && isTruthy(r["Passes Threshold"])) pass++;
      const cur = byOrg.get(org) || { org, reads: 0, n: 0, tass: 0, samples: new Set() };
      cur.reads += rd;
      cur.n++;
      cur.tass = Math.max(cur.tass, ts);
      cur.samples.add(r["Specimen ID"]);
      byOrg.set(org, cur);
    });
    const orgs = Array.from(byOrg.values()).sort((a, b) => b.reads - a.reads || b.tass - a.tass);
    return {
      samples,
      nSamples: samples.length,
      nDetections: rows.length,
      nPassing: pass,
      nOrganisms: byOrg.size,
      reads,
      meanTass: tassN ? tassSum / tassN : 0,
      orgs,
      orgSet: new Set(byOrg.keys()),
    };
  }

  /* ── Leaflet overlay ────────────────────────────────────────────────── */

  function _map() {
    return typeof _leafletMap !== "undefined" ? _leafletMap : null;
  }

  function _drawRegion(reg) {
    const m = _map();
    if (!m || typeof L === "undefined") return;
    if (_layers.has(reg.id)) return;
    const layer = L.polygon(reg.pts, {
      color: reg.color,
      weight: 2,
      opacity: 0.95,
      fillColor: reg.color,
      fillOpacity: 0.1,
      // Click-through: a region must never swallow a click meant for a pin.
      interactive: false,
    });
    layer.addTo(m);
    // A small always-on label so a region is identifiable without the chips.
    const b = layer.getBounds();
    const tip = L.marker(b.getNorthWest ? b.getNorthWest() : reg.pts[0], {
      interactive: false,
      icon: L.divIcon({
        className: "mapdraw-tip",
        html:
          '<span style="background:' +
          reg.color +
          ";color:#fff;font:600 10px system-ui,sans-serif;padding:1px 6px;" +
          'border-radius:4px;white-space:nowrap;box-shadow:0 1px 3px rgba(0,0,0,.3)">' +
          _esc(reg.name) +
          "</span>",
        iconSize: null,
      }),
    });
    tip.addTo(m);
    _layers.set(reg.id, { shape: layer, tip });
  }

  function _undrawRegion(id) {
    const m = _map();
    const l = _layers.get(id);
    if (!l) return;
    if (m) {
      if (l.shape) m.removeLayer(l.shape);
      if (l.tip) m.removeLayer(l.tip);
    }
    _layers.delete(id);
  }

  /* Bring the overlay in line with `_visible` + the current region list. This
     is the ONLY place layers are added or removed, which is what keeps the
     show/hide switch cheap and lossless. */
  function _syncLayers() {
    if (!_map()) return;
    if (!_visible) {
      Array.from(_layers.keys()).forEach(_undrawRegion);
      return;
    }
    const live = new Set(_regions.map((r) => r.id));
    Array.from(_layers.keys()).forEach((id) => {
      if (!live.has(id)) _undrawRegion(id);
    });
    _regions.forEach(_drawRegion);
  }

  /* ── drawing interaction ────────────────────────────────────────────── */

  function _arm(tool) {
    const m = _map();
    if (!m) return;
    _armed = _armed === tool ? null : tool;
    // Switching tools (or switching off) mid-stroke drops the half-drawn shape.
    if (!_armed && typeof _abortDraft === "function") _abortDraft();
    const c = document.getElementById("map-container");
    if (c) c.classList.toggle("mapdraw-arming", !!_armed);
    // The floating cards and the chart overlay sit on top of the map, so while
    // a tool is armed they fade and stop taking pointer events — a drag that
    // starts over them still draws on the map underneath instead of silently
    // doing nothing.
    const split = document.getElementById("map-split");
    if (split) split.classList.toggle("mapdraw-arming", !!_armed);
    // While armed the map must not pan under the pointer, or a lasso would
    // drag the world instead of tracing a shape.
    if (_armed) {
      m.dragging.disable();
      m.boxZoom.disable();
    } else {
      m.dragging.enable();
      m.boxZoom.enable();
    }
    _renderBar();
  }

  function _disarm() {
    if (_armed) _arm(_armed); // toggles it off, restores dragging + cursor
  }

  /* The lowest slot number not currently in use — NOT a running counter.
     Delete "Region 1" and the next shape drawn is "Region 1" again, in the
     same colour, instead of climbing to "Region 3" and leaving a gap that
     reads like a region went missing. A region the user has renamed keeps
     holding its slot, so its colour is never handed to something else while
     it is still on the map. */
  function _freeSlot() {
    const taken = new Set(_regions.map((r) => r.slot).filter((n) => n));
    let n = 1;
    while (taken.has(n)) n++;
    return n;
  }

  function _finishDraft(pts, type) {
    const m = _map();
    if (_draft && m) m.removeLayer(_draft);
    _draft = null;
    // Guard against a stray click registering as a zero-area region.
    if (!pts || pts.length < 3) return;
    const id = "reg" + ++_seq;
    const slot = _freeSlot();
    const reg = {
      id,
      name: "Region " + slot,
      slot,
      type,
      pts,
      color: REGION_COLORS[(slot - 1) % REGION_COLORS.length],
    };
    _regions.push(reg);
    // Drawing something implies you want to see it.
    _visible = true;
    const t = document.getElementById("mapdraw-toggle");
    if (t) t.checked = true;
    _syncLayers();
    _publishGrouping();
    _render();
  }

  /* All pointer handling is done with DOM listeners on the map container, in
     the CAPTURE phase, rather than with Leaflet's own map events. Leaflet's
     map-level mousedown / mousemove never fire when the pointer is over a
     marker, a cluster bubble or a control — and a drag that starts on a pin,
     or crosses one, is precisely the drag a user makes when circling samples.
     Capturing on the container sees every one of them. */
  function _wireMapDrawing() {
    const m = _map();
    if (!m || _wired || typeof L === "undefined") return;
    const cont = m.getContainer();
    if (!cont) return;
    _wired = true;

    let start = null; // latlng where the drag began
    let last = null; // latlng under the pointer right now
    let pts = null; // lasso trace

    const ll = (e) => {
      try {
        return m.mouseEventToLatLng(e);
      } catch (err) {
        return null;
      }
    };

    function _begin(e) {
      if (!_armed || start) return;
      const p0 = ll(e);
      if (!p0) return;
      // Swallow the event so a pin under the cursor doesn't also select itself.
      e.preventDefault();
      e.stopPropagation();
      start = p0;
      last = p0;
      pts = [[p0.lat, p0.lng]];
      if (_draft) m.removeLayer(_draft);
      _draft =
        _armed === "rect"
          ? L.rectangle(L.latLngBounds(p0, p0), {
              color: "#1565c0",
              weight: 2,
              dashArray: "4 3",
              fillOpacity: 0.08,
              interactive: false,
            })
          : L.polyline([p0], { color: "#1565c0", weight: 2, dashArray: "4 3", interactive: false });
      _draft.addTo(m);
    }

    function _extend(e) {
      if (!_armed || !start || !_draft) return;
      const p1 = ll(e);
      if (!p1) return;
      last = p1;
      if (_armed === "rect") {
        _draft.setBounds(L.latLngBounds(start, p1));
        return;
      }
      // Thin the trace: a raw mousemove stream produces thousands of points for
      // one lasso, which makes the polygon (and the saved session) far heavier
      // than the shape needs.
      const prev = pts[pts.length - 1];
      const b = m.getBounds();
      const tol = Math.abs(b.getNorth() - b.getSouth()) / 400;
      if (Math.abs(prev[0] - p1.lat) < tol && Math.abs(prev[1] - p1.lng) < tol) return;
      pts.push([p1.lat, p1.lng]);
      _draft.setLatLngs(pts);
    }

    /* Close the shape. Bound to a document-level mouseup so a drag that ends
       outside the map — or on top of a control — still finishes cleanly. */
    function _finish() {
      if (!_armed || !start) return;
      const tool = _armed;
      const end = last || start;
      let out;
      if (tool === "rect") {
        const bb = L.latLngBounds(start, end);
        // Ignore an accidental click with no drag behind it.
        if (Math.abs(bb.getNorth() - bb.getSouth()) < 1e-9 && Math.abs(bb.getEast() - bb.getWest()) < 1e-9) {
          out = null;
        } else {
          out = [
            [bb.getNorth(), bb.getWest()],
            [bb.getNorth(), bb.getEast()],
            [bb.getSouth(), bb.getEast()],
            [bb.getSouth(), bb.getWest()],
          ];
        }
      } else {
        out = pts && pts.length >= 3 ? pts.slice() : null;
      }
      start = null;
      last = null;
      pts = null;
      _finishDraft(out, tool);
      // One shape per arming: the tool disarms itself so the very next drag
      // pans the map again, which is what people expect after drawing.
      _disarm();
    }

    function _abort() {
      if (!start) return;
      start = null;
      last = null;
      pts = null;
      if (_draft) {
        m.removeLayer(_draft);
        _draft = null;
      }
    }

    cont.addEventListener("mousedown", _begin, true);
    cont.addEventListener("mousemove", _extend, true);
    // A click that follows the drag would otherwise reach the marker under the
    // release point and open its side panel.
    cont.addEventListener(
      "click",
      (e) => {
        if (_armed) {
          e.preventDefault();
          e.stopPropagation();
        }
      },
      true,
    );
    document.addEventListener("mouseup", () => {
      if (_armed && start) _finish();
    });
    // Deliberately NO abort-on-mouseleave: a trace that strays past the edge of
    // the map — which happens constantly when circling samples near the border
    // — keeps the points it already has, and the document-level mouseup above
    // closes the shape wherever the button is finally released.
    _abortDraft = _abort;
  }

  /* ── "Group by region" — regions as a metadata column ───────────────── */

  /* Write (or clear) `map_region` on every RUN_META record, then let the
     grouping engine re-profile. Publishing it as an ordinary metadata column
     is what makes the drawn regions work everywhere a metadata column works:
     the Group by bar, the group heatmap, the network, cross-entry comparison,
     the map's own colouring, and the exported metadata table. */
  function _publishGrouping() {
    if (typeof RUN_META === "undefined") return;
    if (!_asGrouping || !_regions.length) {
      let touched = false;
      RUN_META.forEach((r) => {
        if ("map_region" in r) {
          delete r.map_region;
          touched = true;
        }
      });
      if (touched) _afterGroupingChange(true);
      return;
    }
    const byName = new Map();
    _regions.forEach((reg) => {
      _samplesIn(reg).forEach((sn) => {
        // A sample inside two overlapping regions belongs to both — the
        // grouping engine already handles multi-valued cells by fanning the
        // sample into every combination.
        const cur = byName.get(sn);
        if (cur) cur.push(reg.name);
        else byName.set(sn, [reg.name]);
      });
    });
    RUN_META.forEach((r) => {
      const hit = byName.get(r.sample_name);
      r.map_region = hit ? (hit.length === 1 ? hit[0] : hit.slice()) : "Outside regions";
    });
    _afterGroupingChange(false);
  }

  function _afterGroupingChange(cleared) {
    try {
      if (typeof _normalizeMetaRecord === "function") RUN_META.forEach(_normalizeMetaRecord);
    } catch (e) {}
    try {
      if (window.metaGrouping) window.metaGrouping.refresh();
    } catch (e) {}
    try {
      if (typeof _mgSyncGroupBar === "function") _mgSyncGroupBar();
    } catch (e) {}
    try {
      if (typeof _updateMetaSubTabStates === "function") _updateMetaSubTabStates();
    } catch (e) {}
    try {
      if (typeof _buildRunMetaTable === "function") _buildRunMetaTable();
    } catch (e) {}
    // Select the column the first time it is published, so "Group by region"
    // visibly does something instead of only making a column available.
    try {
      if (!cleared && window.metaGrouping) {
        const f = window.metaGrouping.fields();
        if (f.indexOf("map_region") === -1) window.metaGrouping.setFields(f.concat(["map_region"]));
      } else if (cleared && window.metaGrouping) {
        const f = window.metaGrouping.fields().filter((x) => x !== "map_region");
        window.metaGrouping.setFields(f);
      }
    } catch (e) {}
  }

  /* ── rendering: toolbar chips ───────────────────────────────────────── */

  function _renderBar() {
    const lasso = document.getElementById("mapdraw-lasso");
    const rect = document.getElementById("mapdraw-rect");
    if (lasso) lasso.classList.toggle("armed", _armed === "lasso");
    if (rect) rect.classList.toggle("armed", _armed === "rect");

    const hint = document.getElementById("mapdraw-hint");
    if (hint) {
      hint.textContent = _armed
        ? _armed === "lasso"
          ? "Hold the mouse down on the map and trace an area — release to finish."
          : "Click and drag a box on the map — release to finish."
        : _regions.length
        ? ""
        : "Draw an area to summarise the samples inside it.";
    }
    const count = document.getElementById("mapdraw-count");
    if (count) count.textContent = _regions.length ? `(${_regions.length})` : "";
    const clear = document.getElementById("mapdraw-clear");
    if (clear) clear.style.display = _regions.length ? "inline-flex" : "none";
    const grpWrap = document.getElementById("mapdraw-group-wrap");
    if (grpWrap) grpWrap.style.display = _regions.length ? "inline-flex" : "none";

    const chips = document.getElementById("mapdraw-chips");
    if (!chips) return;
    chips.innerHTML = _regions
      .map((r) => {
        const n = _samplesIn(r).length;
        return (
          `<span class="mapdraw-chip${_visible ? "" : " off"}" style="border-left-color:${r.color}" data-id="${
            r.id
          }">` +
          `<span class="mapdraw-chip-name" contenteditable="true" spellcheck="false" ` +
          `title="Click to rename">${_esc(r.name)}</span>` +
          `<span class="mapdraw-chip-n">${n} sample${n === 1 ? "" : "s"}</span>` +
          `<button type="button" class="mapdraw-chip-zoom" title="Zoom the map to this region"><i class="fas fa-crosshairs"></i></button>` +
          `<button type="button" class="mapdraw-chip-del" title="Delete this region">&times;</button>` +
          `</span>`
        );
      })
      .join("");

    chips.querySelectorAll(".mapdraw-chip").forEach((chip) => {
      const id = chip.getAttribute("data-id");
      const reg = _regions.find((r) => r.id === id);
      if (!reg) return;
      const nameEl = chip.querySelector(".mapdraw-chip-name");
      if (nameEl) {
        const commit = () => {
          const v = nameEl.textContent.trim();
          reg.name = v || reg.name;
          nameEl.textContent = reg.name;
          _undrawRegion(reg.id); // redraw so the on-map label follows the name
          _syncLayers();
          _publishGrouping();
          _render();
        };
        nameEl.addEventListener("blur", commit);
        nameEl.addEventListener("keydown", (e) => {
          if (e.key === "Enter") {
            e.preventDefault();
            nameEl.blur();
          }
        });
      }
      const zoom = chip.querySelector(".mapdraw-chip-zoom");
      if (zoom)
        zoom.addEventListener("click", () => {
          const m = _map();
          if (m) m.fitBounds(reg.pts, { padding: [40, 40] });
        });
      const del = chip.querySelector(".mapdraw-chip-del");
      if (del)
        del.addEventListener("click", () => {
          _undrawRegion(reg.id);
          _regions = _regions.filter((r) => r.id !== reg.id);
          _syncLayers();
          _publishGrouping();
          _render();
        });
    });
  }

  /* ── rendering: per-region cards + cross-region comparison ──────────── */

  function _card(reg, st, compact) {
    const top = st.orgs.slice(0, compact ? 3 : 5);
    return (
      `<div class="mapdraw-card${compact ? " compact" : ""}" data-id="${reg.id}" ` +
      `style="border-top-color:${reg.color}">` +
      `<h4><span class="mapdraw-swatch" style="background:${reg.color}"></span>` +
      `<span class="mapdraw-card-name">${_esc(reg.name)}</span>` +
      // Corner controls: the graphic for this region, and a jump to it.
      `<button type="button" class="mapdraw-card-zoom" title="Zoom the map to this region">` +
      `<i class="fas fa-crosshairs"></i></button>` +
      `<button type="button" class="mapdraw-chart-btn" title="Show a chart of what is in this region">` +
      `<i class="fas fa-chart-simple"></i></button></h4>` +
      `<div class="mapdraw-kv">` +
      `<span>Samples</span><b>${st.nSamples}</b>` +
      `<span>Detections</span><b>${st.nDetections}</b>` +
      `<span>Passing threshold</span><b>${st.nPassing}</b>` +
      `<span>Distinct organisms</span><b>${st.nOrganisms}</b>` +
      `<span>Reads aligned</span><b>${_fmt(st.reads)}</b>` +
      `<span>Mean TASS</span><b>${st.meanTass.toFixed(1)}</b>` +
      `</div>` +
      (top.length
        ? `<div class="mapdraw-top">Top organisms<ol>` +
          top
            .map(
              (o) =>
                `<li>${_esc(o.org)} <span style="color:#90a4ae">— ${_fmt(o.reads)} reads, ` +
                `${o.samples.size} sample${o.samples.size === 1 ? "" : "s"}</span></li>`,
            )
            .join("") +
          `</ol></div>`
        : `<div class="mapdraw-top mapdraw-empty">No detections inside this region under the current filters.</div>`) +
      `</div>`
    );
  }

  function _comparison(regs, stats) {
    // Metric bars, scaled per column so regions are comparable at a glance.
    const cols = [
      { key: "nSamples", label: "Samples", fmt: (v) => v },
      { key: "nDetections", label: "Detections", fmt: (v) => v },
      { key: "nOrganisms", label: "Organisms", fmt: (v) => v },
      { key: "reads", label: "Reads aligned", fmt: _fmt },
      { key: "meanTass", label: "Mean TASS", fmt: (v) => v.toFixed(1) },
    ];
    const max = {};
    cols.forEach((c) => (max[c.key] = Math.max(1, ...stats.map((s) => s[c.key] || 0))));

    let html =
      `<div class="mapdraw-sec-title"><i class="fas fa-scale-balanced"></i> Region comparison</div>` +
      `<div class="mapdraw-cmp"><table><thead><tr><th>Region</th>` +
      cols.map((c) => `<th>${c.label}</th>`).join("") +
      `</tr></thead><tbody>` +
      regs
        .map((r, i) => {
          const s = stats[i];
          return (
            `<tr><td><span class="mapdraw-swatch" style="background:${r.color}"></span>${_esc(r.name)}</td>` +
            cols
              .map((c) => {
                const v = s[c.key] || 0;
                const pct = Math.round((v / max[c.key]) * 100);
                return (
                  `<td class="mapdraw-bar-cell">` +
                  `<span class="mapdraw-bar-fill" style="width:${pct}%;background:${r.color}"></span>` +
                  `<span>${c.fmt(v)}</span></td>`
                );
              })
              .join("") +
            `</tr>`
          );
        })
        .join("") +
      `</tbody></table></div>`;

    // Shared vs region-only organisms — the question a two-region draw is
    // usually asking ("what is in this area that isn't in that one?").
    const shared = [];
    const all = new Set();
    stats.forEach((s) => s.orgSet.forEach((o) => all.add(o)));
    all.forEach((o) => {
      if (stats.every((s) => s.orgSet.has(o))) shared.push(o);
    });
    const only = regs.map((r, i) =>
      Array.from(stats[i].orgSet).filter((o) => stats.every((s, j) => j === i || !s.orgSet.has(o))),
    );

    html +=
      `<div class="mapdraw-sec-title"><i class="fas fa-circle-nodes"></i> Organism overlap</div>` +
      `<div class="mapdraw-cards">` +
      `<div class="mapdraw-card" style="border-top-color:#607d8b">` +
      `<h4><i class="fas fa-object-group" style="color:#607d8b"></i> In every region (${shared.length})</h4>` +
      (shared.length
        ? `<div class="mapdraw-top"><ol>${shared
            .slice(0, 12)
            .map((o) => `<li>${_esc(o)}</li>`)
            .join("")}</ol>${shared.length > 12 ? `<div>…and ${shared.length - 12} more</div>` : ""}</div>`
        : `<div class="mapdraw-empty">No organism is present in all of them.</div>`) +
      `</div>` +
      regs
        .map(
          (r, i) =>
            `<div class="mapdraw-card" style="border-top-color:${r.color}">` +
            `<h4><span class="mapdraw-swatch" style="background:${r.color}"></span>Only in ${_esc(r.name)} (${
              only[i].length
            })</h4>` +
            (only[i].length
              ? `<div class="mapdraw-top"><ol>${only[i]
                  .slice(0, 12)
                  .map((o) => `<li>${_esc(o)}</li>`)
                  .join("")}</ol>${only[i].length > 12 ? `<div>…and ${only[i].length - 12} more</div>` : ""}</div>`
              : `<div class="mapdraw-empty">Nothing unique to this region.</div>`) +
            `</div>`,
        )
        .join("") +
      `</div>`;

    return html;
  }

  /* ── rendering: the floating card stack ─────────────────────────────────
     The per-region cards sit ON the map rather than under it. Two reasons:
     the numbers describe a shape you are looking at, so they belong next to
     it, and the map is tall enough that a summary below the fold was a
     scroll away from the thing it summarises.

     The stack collapses to its header bar (click the header, or the caret),
     so it never permanently covers the corner of the map it sits in. That is
     separate from the "Show regions" switch: this hides the NUMBERS, that
     hides the SHAPES.                                                       */

  function _floatHost() {
    // #map-split is moved into the Mapping tab at runtime, so the panel is
    // (re)attached to wherever it currently lives rather than to a fixed spot.
    const split = document.getElementById("map-split");
    if (!split) return null;
    let el = document.getElementById("mapdraw-float");
    if (!el) {
      el = document.createElement("div");
      el.id = "mapdraw-float";
      split.appendChild(el);
    } else if (el.parentElement !== split) {
      split.appendChild(el);
    }
    return el;
  }

  function _renderFloat() {
    const host = _floatHost();
    if (!host) return;
    if (!_regions.length) {
      host.style.display = "none";
      host.innerHTML = "";
      _closeChart();
      return;
    }
    host.style.display = "";
    host.classList.toggle("collapsed", !_floatOpen);

    const stats = _regions.map(_stats);
    host.innerHTML =
      `<div id="mapdraw-float-head" title="Click to collapse / expand the region summaries">` +
      `<i class="fas fa-draw-polygon"></i>` +
      `<span>Regions <b>${_regions.length}</b></span>` +
      `<i class="fas ${_floatOpen ? "fa-chevron-up" : "fa-chevron-down"}" id="mapdraw-float-caret"></i>` +
      `</div>` +
      `<div id="mapdraw-float-body">` +
      _regions.map((r, i) => _card(r, stats[i], true)).join("") +
      `</div>`;

    const head = host.querySelector("#mapdraw-float-head");
    if (head)
      head.addEventListener("click", () => {
        _floatOpen = !_floatOpen;
        _renderFloat();
      });
    _wireCards(host);
    // Keep the map from panning when the pointer is working inside the panel.
    if (typeof L !== "undefined" && L.DomEvent) {
      L.DomEvent.disableClickPropagation(host);
      L.DomEvent.disableScrollPropagation(host);
    }
  }

  /* Shared by the floating stack and the comparison section below the map. */
  function _wireCards(root) {
    root.querySelectorAll(".mapdraw-card[data-id]").forEach((card) => {
      const reg = _regions.find((r) => r.id === card.getAttribute("data-id"));
      if (!reg) return;
      const chart = card.querySelector(".mapdraw-chart-btn");
      if (chart)
        chart.addEventListener("click", (e) => {
          e.stopPropagation();
          _openChart(reg.id);
        });
      const zoom = card.querySelector(".mapdraw-card-zoom");
      if (zoom)
        zoom.addEventListener("click", (e) => {
          e.stopPropagation();
          const m = _map();
          if (m) m.fitBounds(reg.pts, { padding: [40, 40] });
        });
    });
  }

  /* ── rendering: the chart overlay ───────────────────────────────────────
     A "what is actually in there" panel for ONE region, opened from the chart
     button in the corner of its card. Two views:

       bars     the region's top organisms by reads — the quick read
       heatmap  regions × organisms, only offered with two or more regions,
                because with one region every cell is the same row

     Drawn as plain inline SVG rather than through the charting stack: it has
     to render inside a floating panel over a Leaflet map, at whatever size
     the panel happens to be, without waiting on a layout pass.              */

  const HEAT_RAMP = ["#eef4fb", "#cfe1f6", "#9ec6ec", "#66a5de", "#3b82c9", "#1565c0"];

  function _chartHost() {
    const split = document.getElementById("map-split");
    if (!split) return null;
    let el = document.getElementById("mapdraw-chart");
    if (!el) {
      el = document.createElement("div");
      el.id = "mapdraw-chart";
      el.style.display = "none";
      split.appendChild(el);
    } else if (el.parentElement !== split) {
      split.appendChild(el);
    }
    return el;
  }

  function _openChart(id) {
    _chartRegion = id;
    if (_chartView === "heat" && _regions.length < 2) _chartView = "bars";
    _renderChart();
  }

  function _closeChart() {
    _chartRegion = null;
    const el = document.getElementById("mapdraw-chart");
    if (el) el.style.display = "none";
  }

  function _barsSvg(reg, st) {
    const top = st.orgs.slice(0, 8);
    if (!top.length) return '<div class="mapdraw-empty">Nothing detected inside this region right now.</div>';
    const W = 430,
      rowH = 26,
      padL = 8,
      padR = 8,
      labelH = 13;
    const H = top.length * rowH + 8;
    const max = Math.max(...top.map((o) => o.reads), 1);
    const barW = W - padL - padR;
    const rows = top
      .map((o, i) => {
        const y = i * rowH + 4;
        const w = Math.max(2, (o.reads / max) * barW);
        return (
          `<g><title>${_esc(o.org)} — ${_fmt(o.reads)} reads in ${o.samples.size} sample${
            o.samples.size === 1 ? "" : "s"
          }</title>` +
          `<rect x="${padL}" y="${y + labelH}" width="${w}" height="8" rx="3" fill="${reg.color}" opacity="0.85"/>` +
          `<text x="${padL}" y="${y + labelH - 3}" font-size="10" fill="#37474f">${_esc(
            o.org.length > 46 ? o.org.slice(0, 44) + "…" : o.org,
          )}</text>` +
          `<text x="${W - padR}" y="${y + labelH - 3}" font-size="10" text-anchor="end" fill="#90a4ae">${_fmt(
            o.reads,
          )}</text></g>`
        );
      })
      .join("");
    return (
      `<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}" role="img" ` +
      `aria-label="Top organisms in ${_esc(reg.name)} by reads aligned">${rows}</svg>` +
      `<div class="mapdraw-chart-foot">Bar length is reads aligned, scaled to the biggest organism in this region.</div>`
    );
  }

  function _heatSvg(highlightId) {
    const stats = _regions.map(_stats);
    // Column set: the organisms that matter across ALL regions, so the same
    // columns line up for every row and a gap is visible as an empty cell.
    const totals = new Map();
    stats.forEach((s) => s.orgs.forEach((o) => totals.set(o.org, (totals.get(o.org) || 0) + o.reads)));
    const cols = Array.from(totals.entries())
      .sort((a, b) => b[1] - a[1])
      .slice(0, 8)
      .map((e) => e[0]);
    if (!cols.length) return '<div class="mapdraw-empty">No detections in any region right now.</div>';

    const byRegion = stats.map((s) => {
      const m = new Map();
      s.orgs.forEach((o) => m.set(o.org, o.reads));
      return m;
    });
    // Shaded per COLUMN: the question is "which region has more of this
    // organism", and a single global scale would flatten every rare taxon.
    const colMax = cols.map((c) => Math.max(1, ...byRegion.map((m) => m.get(c) || 0)));

    const rowLabelW = 96,
      cellW = 42,
      cellH = 26,
      headH = 74;
    const W = rowLabelW + cols.length * cellW + 8;
    const H = headH + _regions.length * cellH + 6;

    const head = cols
      .map((c, j) => {
        const x = rowLabelW + j * cellW + cellW / 2;
        const short = c.length > 20 ? c.slice(0, 19) + "…" : c;
        return (
          `<g transform="translate(${x},${headH - 6}) rotate(-55)">` +
          `<title>${_esc(c)}</title>` +
          `<text font-size="9.5" fill="#546e7a">${_esc(short)}</text></g>`
        );
      })
      .join("");

    const body = _regions
      .map((r, i) => {
        const y = headH + i * cellH;
        const cells = cols
          .map((c, j) => {
            const v = byRegion[i].get(c) || 0;
            const t = v / colMax[j];
            const shade =
              v === 0 ? "#fbfdff" : HEAT_RAMP[Math.min(HEAT_RAMP.length - 1, Math.ceil(t * (HEAT_RAMP.length - 1)))];
            return (
              `<g><title>${_esc(r.name)} × ${_esc(c)} — ${v ? _fmt(v) + " reads" : "not detected"}</title>` +
              `<rect x="${rowLabelW + j * cellW + 1}" y="${y + 1}" width="${cellW - 2}" height="${
                cellH - 2
              }" rx="3" fill="${shade}" stroke="#e6edf5"/>` +
              (v
                ? `<text x="${rowLabelW + j * cellW + cellW / 2}" y="${
                    y + cellH / 2 + 3
                  }" font-size="9" text-anchor="middle" fill="${t > 0.6 ? "#fff" : "#546e7a"}">${_fmt(v)}</text>`
                : "") +
              `</g>`
            );
          })
          .join("");
        const name = r.name.length > 14 ? r.name.slice(0, 13) + "…" : r.name;
        return (
          `<g${r.id === highlightId ? ' class="mapdraw-heat-hi"' : ""}>` +
          `<rect x="0" y="${y + 1}" width="4" height="${cellH - 2}" rx="2" fill="${r.color}"/>` +
          `<text x="10" y="${y + cellH / 2 + 3}" font-size="10" fill="#37474f"${
            r.id === highlightId ? ' font-weight="700"' : ""
          }><title>${_esc(r.name)}</title>${_esc(name)}</text>` +
          cells +
          `</g>`
        );
      })
      .join("");

    return (
      `<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}" role="img" ` +
      `aria-label="Reads aligned per region and organism">${head}${body}</svg>` +
      `<div class="mapdraw-chart-foot">Reads aligned, shaded within each column — darkest cell is the region ` +
      `with most of that organism.</div>`
    );
  }

  function _renderChart() {
    const host = _chartHost();
    if (!host) return;
    const reg = _regions.find((r) => r.id === _chartRegion);
    if (!reg) {
      host.style.display = "none";
      return;
    }
    const st = _stats(reg);
    const canHeat = _regions.length > 1;
    host.style.display = "";
    host.innerHTML =
      `<div id="mapdraw-chart-head">` +
      `<span class="mapdraw-swatch" style="background:${reg.color}"></span>` +
      `<span class="mapdraw-chart-title">${_esc(reg.name)}</span>` +
      `<span class="mapdraw-chart-sub">${st.nSamples} sample${st.nSamples === 1 ? "" : "s"} · ${
        st.nOrganisms
      } organism${st.nOrganisms === 1 ? "" : "s"}</span>` +
      `<span class="mapdraw-chart-views">` +
      `<button type="button" class="mapdraw-view-btn${_chartView === "bars" ? " active" : ""}" data-view="bars">` +
      `<i class="fas fa-chart-simple"></i> Bars</button>` +
      (canHeat
        ? `<button type="button" class="mapdraw-view-btn${
            _chartView === "heat" ? " active" : ""
          }" data-view="heat"><i class="fas fa-border-all"></i> Heatmap</button>`
        : "") +
      `</span>` +
      `<button type="button" id="mapdraw-chart-close" title="Close">&times;</button>` +
      `</div>` +
      `<div id="mapdraw-chart-body">${_chartView === "heat" && canHeat ? _heatSvg(reg.id) : _barsSvg(reg, st)}</div>`;

    host.querySelectorAll(".mapdraw-view-btn").forEach((b) =>
      b.addEventListener("click", () => {
        _chartView = b.getAttribute("data-view");
        _renderChart();
      }),
    );
    const close = host.querySelector("#mapdraw-chart-close");
    if (close) close.addEventListener("click", _closeChart);
    if (typeof L !== "undefined" && L.DomEvent) {
      L.DomEvent.disableClickPropagation(host);
      L.DomEvent.disableScrollPropagation(host);
    }
  }

  /* ── rendering: the comparison section under the map ────────────────── */

  function _renderPanel() {
    const host = document.getElementById("mapdraw-report");
    if (!host) return;
    // The per-region numbers now live on the map; what stays below it is the
    // side-by-side comparison, which is too wide to float over anything.
    if (_regions.length < 2) {
      host.style.display = "none";
      host.innerHTML = "";
      return;
    }
    host.style.display = "";
    host.innerHTML = _comparison(_regions, _regions.map(_stats));
    _wireCards(host);
  }

  function _render() {
    _renderBar();
    _renderFloat();
    _renderPanel();
    if (_chartRegion) _renderChart();
  }

  /* ── wiring ─────────────────────────────────────────────────────────── */

  function _wireToolbar() {
    const bar = document.getElementById("mapdraw-bar");
    if (!bar || bar.getAttribute("data-wired")) return;
    bar.setAttribute("data-wired", "1");

    bar.querySelectorAll("[data-draw]").forEach((b) => {
      b.addEventListener("click", () => {
        if (!_map()) {
          const h = document.getElementById("mapdraw-hint");
          if (h) h.textContent = "The map isn't ready yet — it needs coordinates and an internet connection.";
          return;
        }
        _wireMapDrawing();
        _arm(b.getAttribute("data-draw"));
      });
    });

    const toggle = document.getElementById("mapdraw-toggle");
    if (toggle)
      toggle.addEventListener("change", () => {
        // Regions themselves are untouched — only the overlay comes and goes,
        // so flicking this back on restores exactly what was there.
        _visible = !!toggle.checked;
        _syncLayers();
        _renderBar();
      });

    const group = document.getElementById("mapdraw-group");
    if (group)
      group.addEventListener("change", () => {
        _asGrouping = !!group.checked;
        _publishGrouping();
        _render();
      });

    const clear = document.getElementById("mapdraw-clear");
    if (clear)
      clear.addEventListener("click", () => {
        if (_regions.length > 1 && !confirm(`Delete all ${_regions.length} drawn regions?`)) return;
        Array.from(_layers.keys()).forEach(_undrawRegion);
        _regions = [];
        _publishGrouping();
        _render();
      });
  }

  /* Called when the Mapping tab becomes visible: the map may only now exist. */
  window._ttRegionOnMapShown = function () {
    _wireToolbar();
    _wireMapDrawing();
    _syncLayers();
    _render();
  };

  /* Called after any filter change — membership is fixed, but the stats that
     describe it are filter-scoped, so the cards have to be recomputed. */
  window._ttRegionRefresh = function () {
    if (!_regions.length) return;
    _render();
  };

  window.ttMapRegions = {
    list: () => _regions.map((r) => Object.assign({}, r, { pts: r.pts.slice() })),
    statsFor: (id) => {
      const r = _regions.find((x) => x.id === id);
      return r ? _stats(r) : null;
    },
    capture: () => ({
      regions: _regions.map((r) => ({
        id: r.id,
        name: r.name,
        slot: r.slot,
        type: r.type,
        color: r.color,
        pts: r.pts,
      })),
      visible: _visible,
      asGrouping: _asGrouping,
      seq: _seq,
    }),
    restore: (state) => {
      state = state && typeof state === "object" ? state : {};
      Array.from(_layers.keys()).forEach(_undrawRegion);
      _regions = Array.isArray(state.regions)
        ? state.regions
            .filter((r) => r && Array.isArray(r.pts) && r.pts.length >= 3)
            .map((r, i) => ({
              id: r.id || "reg" + (i + 1),
              name: r.name || "Region " + (i + 1),
              // Slot drives both the auto-name and the colour, and a session
              // saved before slots existed simply falls back to position.
              slot: Number(r.slot) > 0 ? Number(r.slot) : i + 1,
              type: r.type === "rect" ? "rect" : "lasso",
              color: r.color || REGION_COLORS[i % REGION_COLORS.length],
              pts: r.pts,
            }))
        : [];
      _seq = Math.max(Number(state.seq) || 0, _regions.length);
      _visible = state.visible !== false;
      _asGrouping = !!state.asGrouping;
      const t = document.getElementById("mapdraw-toggle");
      if (t) t.checked = _visible;
      const g = document.getElementById("mapdraw-group");
      if (g) g.checked = _asGrouping;
      _syncLayers();
      if (_asGrouping) _publishGrouping();
      _render();
    },
  };

  // The toolbar is static markup, so it can be wired as soon as the DOM is up;
  // the map itself is wired lazily the first time the Mapping tab is opened.
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", _wireToolbar);
  } else {
    _wireToolbar();
  }
})();
