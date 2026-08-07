/* ═══════════════════════════════════════════════════════════════════════
       -  §  FEATURE COMPARE — per-cell sample tooltips + per-sample popup
       -     Hovering any matrix cell shows the relevant sample breakdown for
       -     that organism × metric. Clicking a cell opens a popup table listing
       -     every sample in view with that organism's per-sample values.
       ═══════════════════════════════════════════════════════════════════════ */
const _XS_SAMPLE_PALETTE = [
  "#1565c0",
  "#2e7d32",
  "#c62828",
  "#6a1b9a",
  "#ef6c00",
  "#00838f",
  "#ad1457",
  "#558b2f",
  "#4527a0",
  "#d84315",
  "#0277bd",
  "#9e9d24",
  "#5d4037",
  "#00695c",
  "#7b1fa2",
  "#827717",
];
function _cmpSampleTip(o, m, agg) {
  const all = (agg.sampleList || []).slice();
  const det = new Set(o.samples || []);
  // Lead with the exact value for this cell — it's no longer printed on the
  // matrix, so the tooltip is now the primary place to read it.
  const _cellVal = m.fmt(+o[m.key] || 0);
  const head =
    `<b><i>${o.name}</i></b><br>` +
    `<span style="color:#9bb">${m.label}:</span> <b>${_cellVal}</b>` +
    `<span style="color:#9bb;font-size:.85em"> — ${m.desc}</span><br>`;
  // Sample names are intentionally NOT listed on hover — show counts / summary
  // stats only. The full per-sample list is available via click (popup table).
  let body = "";
  if (m.key === "samplePct" || m.key === "sampleCount") {
    // Percent must use the run-wide specimen total (o.total) — the same
    // denominator the Pass / Below / Total column prints — not the in-view
    // sample list length, which is smaller and inflates the figure.
    const _tot = o.total != null ? o.total : all.length;
    body =
      `<span style="color:#7CFC9B"><b>Passing in ${
        o.passCount != null ? o.passCount : det.size
      } of ${_tot} specimen(s)</b></span>` +
      (m.key === "samplePct" ? ` (${(o.samplePct || 0).toFixed(1)}%)` : "") +
      (o.belowCount ? `<br><span style="color:#9bb">+${o.belowCount} detected below the cutoff</span>` : "");
  } else if (m.key === "reads") {
    body = `<b>Total aligned reads:</b> ${_fmtInt(o.reads || 0)} <span style="color:#9bb">across ${
      det.size
    } sample(s)</span>`;
  } else if (m.key === "aniGroupSize") {
    const members = (agg.aniGroups && agg.aniGroups.get(o.aniGroup)) || [o];
    body = `<b>ANI group — ${members.length} similar reference(s)</b>`;
  } else {
    const map = /tass/i.test(m.key) ? o.tassMap : /cov/i.test(m.key) ? o.covMap : null;
    if (map && map.size) {
      const vals = [...map.values()]
        .map(Number)
        .filter((v) => !isNaN(v))
        .sort((a, b) => a - b);
      const q = (p) => (typeof _xsQuart === "function" ? _xsQuart(vals, p) : vals[Math.floor((vals.length - 1) * p)]);
      body =
        `<span style="color:#9bb">Across ${vals.length} sample(s):</span> ` +
        `min ${vals[0].toFixed(1)} · med ${q(0.5).toFixed(1)} · max ${vals[vals.length - 1].toFixed(1)}`;
    } else {
      body = `<b>${m.label}:</b> ${m.fmt(+o[m.key] || 0)}`;
    }
  }
  return head + body + `<br><span style="color:#90caf9;font-size:.85em">click for full per-sample table</span>`;
}

function _xsOpenSamplePopup(o, m, agg) {
  const ov = document.getElementById("xspop-overlay");
  const body = document.getElementById("xspop-body");
  const titleEl = document.getElementById("xspop-title");
  if (!ov || !body) return;
  const all = (agg.sampleList || []).slice().sort();
  const det = new Set(o.samples || []);
  const tass = o.tassMap || new Map();
  const cov = o.covMap || new Map();
  const hlTass = /tass/i.test(m.key);
  const hlCov = /cov/i.test(m.key);
  const hlDet = m.key === "samplePct" || m.key === "sampleCount";
  if (titleEl) titleEl.innerHTML = `<i>${o.name}</i> — per-sample detail`;
  const hl = (on) => (on ? "background:#fff3bf" : "");
  const rows = all.map((s) => {
    const isDet = det.has(s);
    const tv = tass.has(s) ? (+tass.get(s)).toFixed(1) : "—";
    const cv = cov.has(s) ? (+cov.get(s)).toFixed(1) : "—";
    const html =
      `<tr style="${isDet ? "" : "color:#aaa"}">` +
      `<td style="padding:4px 10px;border-bottom:1px solid #eee;text-align:left">${s}</td>` +
      `<td style="padding:4px 10px;border-bottom:1px solid #eee;text-align:center;${hl(hlDet)}">${
        isDet ? '<span style="color:#2e7d32;font-weight:700">✓</span>' : '<span style="color:#c0392b">✗</span>'
      }</td>` +
      `<td style="padding:4px 10px;border-bottom:1px solid #eee;text-align:right;${hl(hlTass)}">${tv}</td>` +
      `<td style="padding:4px 10px;border-bottom:1px solid #eee;text-align:right;${hl(hlCov)}">${cv}</td>` +
      `</tr>`;
    return { text: s, html };
  });
  const nDet = det.size,
    nTot = all.length;
  const statsHtml =
    `<div style="font-size:.82em;color:#555;margin-bottom:.5em">` +
    `Detected in <b>${nDet}</b> / ${nTot} samples (${((nDet / Math.max(1, nTot)) * 100).toFixed(0)}%) · ` +
    `Mean TASS ${(o.meanTass || 0).toFixed(1)} · Max TASS ${(o.maxTass || 0).toFixed(1)} · ` +
    `Mean Cov ${(o.meanCov || 0).toFixed(1)} · Total reads ${_fmtInt(o.reads || 0)}</div>`;
  const headHtml =
    `<tr style="position:sticky;top:0;background:#f0f4fa">` +
    `<th style="padding:5px 10px;text-align:left;border-bottom:2px solid #1565c0">Sample</th>` +
    `<th style="padding:5px 10px;text-align:center;border-bottom:2px solid #1565c0;${hl(hlDet)}">Detected</th>` +
    `<th style="padding:5px 10px;text-align:right;border-bottom:2px solid #1565c0;${hl(hlTass)}">TASS</th>` +
    `<th style="padding:5px 10px;text-align:right;border-bottom:2px solid #1565c0;${hl(hlCov)}">Coverage %</th>` +
    `</tr>`;
  _xsRenderPagedTable(body, { statsHtml, headHtml, rows, colspan: 4, searchPlaceholder: "Search samples…" });
  ov.style.display = "flex";
}

// ── Reusable paginated + searchable table for detail popups ──────────────
// rows: array of { text, html } where text is the searchable string and html
// is the pre-rendered <tr>. Renders a search box, a scrollable table body, and
// prev/next pager into bodyEl and wires all interactions.
function _xsRenderPagedTable(bodyEl, opts) {
  opts = opts || {};
  const allRows = opts.rows || [];
  const pageSize = opts.pageSize || 20;
  const colspan = opts.colspan || 4;
  const uid = "xspg" + Math.random().toString(36).slice(2, 8);
  let page = 0,
    term = "";
  const btn =
    "padding:4px 12px;border:1px solid #1565c0;border-radius:5px;background:#fff;color:#1565c0;font-size:.82em;cursor:pointer";
  bodyEl.innerHTML =
    (opts.statsHtml || "") +
    `<div style="display:flex;align-items:center;gap:8px;margin-bottom:.5em">` +
    `<input id="${uid}-q" type="text" placeholder="${opts.searchPlaceholder || "Search…"}" ` +
    `style="flex:1;padding:5px 8px;border:1px solid #ccc;border-radius:5px;font-size:.85em"/>` +
    `<span id="${uid}-n" style="font-size:.8em;color:#789;white-space:nowrap"></span></div>` +
    `<div style="overflow:auto;max-height:50vh"><table style="border-collapse:collapse;width:100%;font-size:.85em">` +
    `<thead>${opts.headHtml || ""}</thead><tbody id="${uid}-tb"></tbody></table></div>` +
    `<div style="display:flex;align-items:center;justify-content:center;gap:12px;margin-top:.6em">` +
    `<button id="${uid}-prev" style="${btn}">‹ Prev</button>` +
    `<span id="${uid}-pg" style="font-size:.82em;color:#555"></span>` +
    `<button id="${uid}-next" style="${btn}">Next ›</button></div>`;
  const tb = document.getElementById(uid + "-tb"),
    pgEl = document.getElementById(uid + "-pg"),
    nEl = document.getElementById(uid + "-n"),
    prev = document.getElementById(uid + "-prev"),
    next = document.getElementById(uid + "-next"),
    q = document.getElementById(uid + "-q");
  const filtered = () => {
    const t = term.trim().toLowerCase();
    return t ? allRows.filter((r) => (r.text || "").toLowerCase().includes(t)) : allRows;
  };
  function render() {
    const f = filtered();
    const pages = Math.max(1, Math.ceil(f.length / pageSize));
    if (page >= pages) page = pages - 1;
    if (page < 0) page = 0;
    const slice = f.slice(page * pageSize, (page + 1) * pageSize);
    tb.innerHTML = slice.length
      ? slice.map((r) => r.html).join("")
      : `<tr><td colspan="${colspan}" style="padding:12px;text-align:center;color:#999">No matches</td></tr>`;
    pgEl.textContent = `Page ${page + 1} / ${pages}`;
    nEl.textContent = `${f.length} row(s)`;
    prev.disabled = page <= 0;
    next.disabled = page >= pages - 1;
    prev.style.opacity = prev.disabled ? 0.4 : 1;
    next.style.opacity = next.disabled ? 0.4 : 1;
  }
  prev.addEventListener("click", () => {
    if (page > 0) {
      page--;
      render();
    }
  });
  next.addEventListener("click", () => {
    page++;
    render();
  });
  q.addEventListener("input", () => {
    term = q.value;
    page = 0;
    render();
  });
  render();
}
