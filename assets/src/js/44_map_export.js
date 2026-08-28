/* ═══════════════════════════════════════════════════════════════════════════
       -  §  MAP EXPORT  (Precise lat/long map — Mapping & Geography sub-tab)
       -     Every other plot in the report is an <svg>, so the generic
       -     enhancer in 02_utilities.js can clone it and hand it to the
       -     export modal. The Leaflet map is not: it is a stack of <img>
       -     raster tiles from the CARTO CDN plus absolutely positioned HTML
       -     marker icons. This module composes an equivalent SVG — tiles
       -     embedded as base64 <image>, markers redrawn from their live DOM —
       -     and then reuses _renderSvgExport so the map offers the same
       -     PNG / JPEG / SVG / PDF / HTML formats and the same width/height
       -     controls as the rest of the report.
       -     Entry points: _ensureMapExportControl (adds the on-map button,
       -     called from _doInitMap) and _exportMapView (called by
       -     _runExportFromModal when the modal is in "map" mode).
═══════════════════════════════════════════════════════════════════════════ */

/* Tiles are re-fetched with crossOrigin="anonymous" so the canvas they are
   drawn on stays untainted and can be read back as a data URI. CARTO sends
   Access-Control-Allow-Origin, so this normally succeeds; when it does not
   (offline report, blocked CDN) the tile is dropped and the export falls back
   to markers on a plain background, with one warning to the user. */
const _MAP_TILE_CACHE = new Map();
let _mapExportTileWarned = false;

function _mapTileDataUrl(url) {
  if (!url) return Promise.resolve(null);
  if (_MAP_TILE_CACHE.has(url)) return _MAP_TILE_CACHE.get(url);
  const p = new Promise((resolve) => {
    const img = new Image();
    img.crossOrigin = "anonymous";
    img.onload = () => {
      try {
        const c = document.createElement("canvas");
        c.width = img.naturalWidth || 256;
        c.height = img.naturalHeight || 256;
        c.getContext("2d").drawImage(img, 0, 0);
        resolve(c.toDataURL("image/png"));
      } catch (e) {
        resolve(null); // tainted canvas — tile server did not allow CORS
      }
    };
    img.onerror = () => resolve(null);
    img.src = url;
  });
  _MAP_TILE_CACHE.set(url, p);
  return p;
}

function _mapExpEsc(s) {
  return String(s == null ? "" : s)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;");
}

function _mapExpAttr(s) {
  return _mapExpEsc(s).replace(/"/g, "&quot;");
}

/* A CSS colour that actually paints something (not `transparent` / alpha 0). */
function _mapExpColor(v) {
  if (!v) return null;
  const s = String(v).trim();
  if (s === "transparent" || s === "none") return null;
  const m = s.match(/rgba?\(([^)]+)\)/i);
  if (m) {
    const parts = m[1].split(",").map((x) => parseFloat(x));
    if (parts.length > 3 && !(parts[3] > 0)) return null;
  }
  return s;
}

/* Redraw one live DOM node (a Leaflet marker icon, or anything inside it) as
   SVG, positioned relative to the map container.

   Two cases cover every icon the map builds (see _svgDot / _pieIcon /
   _clusterIcon in 35_tab_map.js):
     • an <svg> icon  → cloned with its computed styles, exactly like the
       other plots are exported;
     • an HTML box    → its background (circle when fully rounded, else rect),
       its border, and its centred text are re-emitted as SVG primitives, then
       its children are walked the same way. That is what turns a cluster
       bubble ("4/40") or a group donut with a count label into real vector
       output instead of a blank gap. */
function _mapNodeToSvg(el, base) {
  if (!el || el.nodeType !== 1) return "";
  const tag = el.tagName ? el.tagName.toLowerCase() : "";
  if (tag === "script" || tag === "style") return "";
  const cs = getComputedStyle(el);
  if (cs.display === "none" || cs.visibility === "hidden") return "";
  const nodeOpacity = parseFloat(cs.opacity);
  if (nodeOpacity === 0) return "";
  const r = el.getBoundingClientRect();
  if (r.width <= 0 || r.height <= 0) return "";
  // Off-canvas markers (Leaflet keeps a margin of them around the viewport).
  if (r.right < base.left - 2 || r.left > base.right + 2 || r.bottom < base.top - 2 || r.top > base.bottom + 2) {
    return "";
  }
  const x = r.left - base.left;
  const y = r.top - base.top;
  const op = nodeOpacity >= 0 && nodeOpacity < 1 ? ` opacity="${nodeOpacity}"` : "";

  if (tag === "svg") {
    const inner = _cloneSvgWithStyles(el, r.width, r.height);
    return `<g transform="translate(${x.toFixed(2)},${y.toFixed(2)})"${op}>${inner}</g>`;
  }

  let out = "";
  const bg = _mapExpColor(cs.backgroundColor);
  const bw = parseFloat(cs.borderTopWidth) || 0;
  const bc = bw > 0 ? _mapExpColor(cs.borderTopColor) : null;
  if (bg || bc) {
    const radius = parseFloat(cs.borderTopLeftRadius) || 0;
    const fill = bg || "none";
    const stroke = bc ? ` stroke="${bc}" stroke-width="${bw}"` : "";
    if (radius >= Math.min(r.width, r.height) / 2 - 1) {
      out +=
        `<ellipse cx="${(x + r.width / 2).toFixed(2)}" cy="${(y + r.height / 2).toFixed(2)}" ` +
        `rx="${(r.width / 2).toFixed(2)}" ry="${(r.height / 2).toFixed(2)}" fill="${fill}"${stroke}${op}/>`;
    } else {
      out +=
        `<rect x="${x.toFixed(2)}" y="${y.toFixed(2)}" width="${r.width.toFixed(2)}" ` +
        `height="${r.height.toFixed(2)}" rx="${radius.toFixed(2)}" fill="${fill}"${stroke}${op}/>`;
    }
  }

  const text = Array.from(el.childNodes)
    .filter((n) => n.nodeType === 3)
    .map((n) => n.textContent)
    .join("")
    .replace(/\s+/g, " ")
    .trim();
  if (text) {
    const fontSize = parseFloat(cs.fontSize) || 12;
    out +=
      `<text x="${(x + r.width / 2).toFixed(2)}" y="${(y + r.height / 2).toFixed(2)}" ` +
      `text-anchor="middle" dominant-baseline="central" ` +
      `font-family="${_mapExpAttr(String(cs.fontFamily).replace(/"/g, "'"))}" font-size="${fontSize}" ` +
      `font-weight="${_mapExpAttr(cs.fontWeight || "400")}" fill="${cs.color || "#222"}"${op}>` +
      `${_mapExpEsc(text)}</text>`;
  }

  Array.from(el.children).forEach((child) => {
    out += _mapNodeToSvg(child, base);
  });
  return out;
}

/* Compose the whole map view (basemap + markers + attribution) into one SVG
   string sized `width` x `height`, drawn in the container's own pixel
   coordinates via a viewBox so the caller can rasterise at any size. */
async function _mapExportSvgText(width, height) {
  const container = document.getElementById("map-container");
  if (!container) throw new Error("The map is not open.");
  const base = container.getBoundingClientRect();
  const W = Math.max(1, Math.round(base.width));
  const H = Math.max(1, Math.round(base.height));
  const body = [`<rect width="${W}" height="${H}" fill="#f2efe9"/>`];

  // ── Basemap tiles ──────────────────────────────────────────────────────
  const tileEls = Array.from(container.querySelectorAll(".leaflet-tile-pane img"));
  const tileJobs = tileEls.map(async (img) => {
    const r = img.getBoundingClientRect();
    if (r.width <= 0 || r.height <= 0) return null;
    if (r.right < base.left - 1 || r.left > base.right + 1 || r.bottom < base.top - 1 || r.top > base.bottom + 1) {
      return null;
    }
    if (parseFloat(getComputedStyle(img).opacity) === 0) return null;
    const data = await _mapTileDataUrl(img.currentSrc || img.src);
    if (!data) return null;
    return (
      `<image x="${(r.left - base.left).toFixed(2)}" y="${(r.top - base.top).toFixed(2)}" ` +
      `width="${r.width.toFixed(2)}" height="${r.height.toFixed(2)}" preserveAspectRatio="none" ` +
      `href="${data}" xlink:href="${data}"/>`
    );
  });
  const tiles = (await Promise.all(tileJobs)).filter(Boolean);
  body.push(...tiles);
  const tilesMissing = tileEls.length > 0 && tiles.length === 0;

  // ── Vector overlays, then marker icons (Leaflet's own paint order) ──────
  ["leaflet-overlay-pane", "leaflet-marker-pane"].forEach((cls) => {
    const pane = container.querySelector("." + cls);
    if (!pane) return;
    Array.from(pane.children).forEach((el) => {
      body.push(_mapNodeToSvg(el, base));
    });
  });

  // ── Attribution — required by the OpenStreetMap / CARTO tile licences,
  //    so it travels with every exported image rather than living only in
  //    the on-screen control. ───────────────────────────────────────────
  const attribution = "© OpenStreetMap contributors © CARTO";
  const attrW = attribution.length * 5.6 + 10;
  body.push(
    `<g><rect x="${(W - attrW - 4).toFixed(2)}" y="${(H - 20).toFixed(2)}" width="${attrW.toFixed(2)}" ` +
      `height="16" fill="rgba(255,255,255,.78)" rx="3"/>` +
      `<text x="${(W - 9).toFixed(2)}" y="${(H - 12).toFixed(2)}" text-anchor="end" dominant-baseline="central" ` +
      `font-family="system-ui, sans-serif" font-size="9.5" fill="#4a4a4a">${_mapExpEsc(attribution)}</text></g>`,
  );

  if (tilesMissing && !_mapExportTileWarned) {
    _mapExportTileWarned = true;
    alert(
      "The basemap tiles could not be embedded (the tile server refused a cross-origin read). " +
        "The exported map shows the sample markers on a plain background.",
    );
  }

  return (
    `<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" ` +
    `width="${width}" height="${height}" viewBox="0 0 ${W} ${H}">${body.join("")}</svg>`
  );
}

/* Called by _runExportFromModal when the export modal is in "map" mode. */
async function _exportMapView(target, opts) {
  const container = document.getElementById("map-container");
  if (!container) return;
  const title = _nearestTitle(container) || "Sample map";
  const filename = _slug(title || "sample-map");
  const width = Math.max(50, parseInt(opts.width, 10) || container.clientWidth || 900);
  const height = Math.max(50, parseInt(opts.height, 10) || container.clientHeight || 500);
  let svgText;
  try {
    svgText = await _mapExportSvgText(width, height);
  } catch (e) {
    alert("Could not export the map: " + (e && e.message ? e.message : e));
    return;
  }
  return _renderSvgExport(svgText, { title, filename, width, height, format: opts.format });
}

/* Open the shared export modal against the map. Same dialog, same formats and
   size fields the SVG plots use. */
function _openMapExport() {
  const container = document.getElementById("map-container");
  if (!container) return;
  const overlay = _ensureExportModal();
  _EXPORT_STATE.mode = "map";
  _EXPORT_STATE.target = container;
  overlay.querySelector("#export-title").textContent = "Export Map";
  overlay.querySelector("#export-format").innerHTML =
    '<option value="png">PNG</option><option value="jpeg">JPEG</option><option value="svg">SVG</option><option value="pdf">PDF</option><option value="html">HTML</option>';
  overlay.querySelector("#export-width").value = Math.max(50, Math.round(container.clientWidth || 900));
  overlay.querySelector("#export-height").value = Math.max(50, Math.round(container.clientHeight || 500));
  _syncExportModalFields();
  overlay.style.display = "flex";
}

/* The download button itself. It is a Leaflet control rather than a plain
   absolutely-positioned button so it sits in the same top-right stack as the
   cluster-mode box and cannot be re-laid-out by Leaflet's panes. Styling and
   the hover-fade come from .tt-map-export in geo.css, and the shared
   .export-control class keeps it out of the printed report PDF. */
function _ensureMapExportControl() {
  if (typeof L === "undefined" || !L.Control || !_leafletMap) return;
  const container = document.getElementById("map-container");
  if (!container || container.querySelector(".tt-map-export")) return;
  const Ctl = L.Control.extend({
    options: { position: "topright" },
    onAdd: function () {
      const div = L.DomUtil.create("div", "tt-map-export");
      const btn = L.DomUtil.create("button", "export-control", div);
      btn.type = "button";
      btn.title = "Export map";
      btn.setAttribute("aria-label", "Export map");
      btn.innerHTML = '<i class="fas fa-download"></i>';
      L.DomEvent.disableClickPropagation(div);
      L.DomEvent.on(btn, "click", (e) => {
        L.DomEvent.stop(e);
        _openMapExport();
      });
      return div;
    },
  });
  _leafletMap.addControl(new Ctl());
}
