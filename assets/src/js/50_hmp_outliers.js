/* ═══════════════════════════════════════════════════════════════════════════
       -  §  HMP HEALTHY-ABUNDANCE OUTLIERS
       -     Port of the ODR PDF's "outlier" information: every detection whose
       -     sample type maps to an HMP body site is compared with how abundant
       -     that organism is in HEALTHY subjects at the same site.
       -
       -       z = (observed % reads − healthy mean %) / healthy std %
       -       z ≥ threshold          → Elevated            (shown normally)
       -       z <  threshold         → Within healthy range (faded + ◆ n (p%),
       -                                 exactly like the PDF's diamond rows)
       -       never seen at the site → Absent in healthy
       -
       -     make_report.py (bin/hmp_outliers.py) stamps the per-row figures
       -     ("HMP …" columns) and embeds the reference distributions it used in
       -     BOOT.hmp = { threshold, sites{site:N}, bins{lo,hi,step},
       -                  ref{"site:taxid": {m,s,n,N,h[],q[],name,rank}} }.
       -     This file only re-derives the status when the analyst moves the
       -     z threshold, re-computes z for merged specimens, and renders the
       -     markers / distribution plots / tooltips.
═══════════════════════════════════════════════════════════════════════════ */

const TT_HMP = (typeof BOOT !== "undefined" && BOOT && BOOT.hmp) || null;
let TT_HMP_Z = TT_HMP && isFinite(parseFloat(TT_HMP.threshold)) ? parseFloat(TT_HMP.threshold) : 2;
let TT_HMP_FADE = true;
const TT_HMP_ELEV = "Elevated";
const TT_HMP_WITHIN = "Within healthy range";
const TT_HMP_ABSENT = "Absent in healthy";
// Normally sterile site (blood / plasma / WB, CSF, ...): HMP has no healthy
// distribution because nothing is expected there, so every detection is outside
// the healthy expectation. Never faded.
const TT_HMP_STERILE = "Sterile site";

function ttHmpAvailable() {
  if (!TT_HMP) return false;
  for (let i = 0; i < DATA.length; i++) if (DATA[i] && DATA[i]["HMP Status"]) return true;
  return false;
}

function _hmpEsc(s) {
  return String(s == null ? "" : s)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}
function _hmpNum(v) {
  if (v === null || v === undefined || v === "") return NaN;
  const n = parseFloat(v);
  return isNaN(n) ? NaN : n;
}
function _hmpFmtPct(v) {
  if (!isFinite(v)) return "—";
  if (v === 0) return "0%";
  if (v >= 10) return v.toFixed(1) + "%";
  if (v >= 0.01) return v.toPrecision(3) + "%";
  return v.toExponential(1) + "%";
}
function _hmpFmtShare(n, N) {
  if (!(N > 0)) return "";
  const p = (100 * n) / N;
  return p > 0 && p < 0.001 ? "<0.001%" : p.toFixed(p < 1 ? 2 : 1) + "%";
}
function _hmpRefEntries(row) {
  if (!TT_HMP || !TT_HMP.ref) return [];
  return String(row["HMP Ref Key"] || "")
    .split("|")
    .filter(Boolean)
    .map((k) => ({ key: k, site: k.split(":")[0], taxid: k.split(":")[1], e: TT_HMP.ref[k] }))
    .filter((x) => x.e);
}

/* Approximate empirical percentile from the shipped log-histogram (used only
   when a merged specimen's abundance changes; the pipeline value is exact). */
function _hmpPercentileFromHist(entries, obs) {
  const b = (TT_HMP && TT_HMP.bins) || { lo: -6, hi: 2, step: 0.5 };
  let le = 0,
    tot = 0;
  entries.forEach(({ e }) => {
    const N = e.N || 0;
    tot += N;
    le += Math.max(0, N - (e.n || 0)); // absent healthy samples = 0 %
    if (!(obs > 0)) return;
    const lx = Math.log10(obs);
    (e.h || []).forEach((c, i) => {
      const lo = b.lo + i * b.step,
        hi = lo + b.step;
      if (lx >= hi) le += c;
      else if (lx > lo) le += (c * (lx - lo)) / b.step;
    });
  });
  return tot > 0 ? Math.min(100, (100 * le) / tot) : NaN;
}

/* (Re)derive a row's HMP status. `recompute` re-calculates z / percentile from
   the row's current "% Reads" (merged specimens); otherwise only the status is
   re-derived from the pipeline z against the active threshold. */
function ttHmpStamp(row, recompute) {
  if (!row || !row["HMP Status"]) return row;
  if (row["HMP Status"] === TT_HMP_ABSENT || row["HMP Status"] === TT_HMP_STERILE) {
    row.__hmpSortZ = 1e9;
    return row;
  }
  if (recompute) {
    const entries = _hmpRefEntries(row);
    const own = String(row["Taxonomic ID #"] || "");
    // Only rows compared at their own taxid; a strain compared at its species
    // keeps the species-level figures computed by the pipeline.
    if (entries.length && entries.every((x) => x.taxid === own)) {
      const obs = Math.max(0, _hmpNum(row["% Reads"]) || 0);
      const m = entries.reduce((s, x) => s + (x.e.m || 0), 0);
      const sd = entries.reduce((s, x) => s + (x.e.s || 0), 0);
      const z = sd > 0 ? (obs - m) / sd : obs > m ? 3 : 0;
      row["HMP Z-Score"] = Math.round(z * 1000) / 1000;
      const p = _hmpPercentileFromHist(entries, obs);
      if (isFinite(p)) row["HMP Healthy Percentile"] = Math.round(p * 100) / 100;
    }
  }
  const z = _hmpNum(row["HMP Z-Score"]);
  if (isFinite(z)) {
    row["HMP Status"] = z >= TT_HMP_Z ? TT_HMP_ELEV : TT_HMP_WITHIN;
    row.__hmpSortZ = z;
  }
  return row;
}
function ttHmpStampAll() {
  if (!TT_HMP) return;
  DATA.forEach((r) => ttHmpStamp(r, false));
}

function ttHmpIsWithin(row) {
  return !!row && row["HMP Status"] === TT_HMP_WITHIN;
}
/* Row class for the faded "within healthy range" look (PDF parity). */
function ttHmpRowClass(row) {
  return TT_HMP_FADE && ttHmpIsWithin(row) ? "hmp-within-row" : "";
}

/* The PDF's marker: ◆ n (p%) ◆ — n healthy reference samples carried the
   organism, p% of all reference samples for the body site. Only on rows that
   are within the healthy range. */
function ttHmpMarkerHTML(row) {
  if (!ttHmpIsWithin(row)) return "";
  const n = _hmpNum(row["HMP Healthy Samples"]);
  const N = _hmpNum(row["HMP Reference Samples"]);
  const lbl = isFinite(n) && N > 0 ? ` ${Math.round(n).toLocaleString()} (${_hmpFmtShare(n, N)})` : "";
  return `<span class="hmp-mark" data-hmp-mark="1">◆<small>${lbl}</small>◆</span>`;
}

/* Map a % abundance onto [0,1] across the log axis 1e-6 % … 100 %. */
function _hmpX(v) {
  const b = (TT_HMP && TT_HMP.bins) || { lo: -6, hi: 2 };
  if (!(v > 0)) return 0;
  return Math.min(1, Math.max(0, (Math.log10(v) - b.lo) / (b.hi - b.lo)));
}
function _hmpStatusColor(st) {
  return st === TT_HMP_ELEV
    ? "#e8590c"
    : st === TT_HMP_ABSENT || st === TT_HMP_STERILE
    ? "#c92a2a"
    : st === TT_HMP_WITHIN
    ? "#2b8a3e"
    : "#868e96";
}

/* Compact cell: healthy range bar (p5–p95 thin, IQR thick, median tick) on a
   log axis with the sample's abundance as a dot, plus z. Summary table column. */
function ttHmpCellHTML(row, key) {
  const st = row["HMP Status"];
  if (!st) return '<span style="color:#ccc" title="No healthy (HMP) reference for this sample type">—</span>';
  if (st === TT_HMP_STERILE)
    return (
      `<span class="hmp-cell" data-hmp-key="${_hmpEsc(
        key || "",
      )}" style="cursor:help;white-space:nowrap;font-size:0.82em;color:#c92a2a;font-weight:700">` +
      `<i class="fas fa-shield-virus" style="margin-right:3px"></i>sterile site</span>`
    );
  const W = 74,
    H = 14,
    pad = 3,
    iw = W - 2 * pad;
  const obs = Math.max(0, _hmpNum(row["% Reads"]) || 0);
  const col = _hmpStatusColor(st);
  let svg = `<svg width="${W}" height="${H}" style="vertical-align:middle" aria-hidden="true">`;
  svg += `<line x1="${pad}" y1="${H / 2}" x2="${W - pad}" y2="${H / 2}" stroke="#dee2e6" stroke-width="1"/>`;
  const ents = _hmpRefEntries(row);
  if (ents.length === 1 && ents[0].e.q && ents[0].e.q.length === 5) {
    const q = ents[0].e.q.map((v) => pad + _hmpX(v) * iw);
    svg += `<line x1="${q[0]}" y1="${H / 2}" x2="${q[4]}" y2="${H / 2}" stroke="#adb5bd" stroke-width="2"/>`;
    svg += `<rect x="${q[1]}" y="${H / 2 - 3}" width="${Math.max(
      1.5,
      q[3] - q[1],
    )}" height="6" fill="#ced4da" stroke="#868e96" stroke-width="0.6"/>`;
    svg += `<line x1="${q[2]}" y1="${H / 2 - 4}" x2="${q[2]}" y2="${H / 2 + 4}" stroke="#495057" stroke-width="1"/>`;
  } else if (ents.length) {
    const m = ents.reduce((s, x) => s + (x.e.m || 0), 0);
    const mx = pad + _hmpX(m) * iw;
    svg += `<line x1="${mx}" y1="${H / 2 - 4}" x2="${mx}" y2="${H / 2 + 4}" stroke="#495057" stroke-width="1"/>`;
  }
  const ox = pad + _hmpX(obs) * iw;
  svg += `<circle cx="${ox}" cy="${H / 2}" r="3.2" fill="${col}" stroke="#fff" stroke-width="0.8"/></svg>`;
  const z = _hmpNum(row["HMP Z-Score"]);
  const zTxt =
    st === TT_HMP_ABSENT
      ? '<span style="color:#c92a2a;font-weight:700" title="Never observed in healthy reference samples for this body site">new</span>'
      : isFinite(z)
      ? `<span style="color:${col};font-weight:${st === TT_HMP_ELEV ? 700 : 400}">z ${
          Math.abs(z) >= 1000 ? z.toExponential(1) : Math.abs(z) >= 100 ? z.toFixed(0) : z.toFixed(1)
        }</span>`
      : "";
  return (
    `<span class="hmp-cell" data-hmp-key="${_hmpEsc(
      key || "",
    )}" style="display:inline-flex;align-items:center;gap:4px;cursor:help;white-space:nowrap">` +
    `${svg}<span style="font-size:0.82em;min-width:38px;text-align:left">${zTxt}</span></span>`
  );
}

/* Full tooltip: figures + healthy distribution histogram with the sample marked. */
function ttHmpTipHTML(row) {
  const st = row["HMP Status"];
  const site = row["HMP Body Site"] || row["Sample Type"] || "";
  if (!st) {
    return (
      `<b>Healthy reference (HMP)</b><br><span style="color:#ccc">No healthy-subject reference for sample type ` +
      `“${_hmpEsc(
        row["Sample Type"] || "unknown",
      )}”. HMP distributions exist for stool, oral, nasal, skin, throat and vaginal sites.</span>`
    );
  }
  const obs = Math.max(0, _hmpNum(row["% Reads"]) || 0);
  const col = _hmpStatusColor(st);
  if (st === TT_HMP_STERILE) {
    return (
      `<div style="max-width:430px"><b>${_hmpEsc(row["Detected Organism"])}</b> · <span style="color:#ccc">${_hmpEsc(
        row["Specimen ID"],
      )}</span><br>` +
      `<span style="display:inline-block;margin:3px 0;padding:0 6px;border-radius:3px;background:${col};color:#fff;font-weight:700;font-size:0.85em">Normally sterile site</span>` +
      ` <span style="color:#ccc">sample type <b style="color:#fff">${_hmpEsc(row["Sample Type"])}</b> → ${_hmpEsc(
        site,
      )}</span>` +
      `<table style="border-collapse:collapse;font-size:0.85em;margin:5px 0"><tr><td style="color:#ffd580;padding-right:10px">This sample</td><td>${_hmpFmtPct(
        obs,
      )} of reads</td></tr><tr><td style="color:#ffd580;padding-right:10px">Healthy expectation</td><td>none — no organisms expected</td></tr></table>` +
      `<span style="color:#aaa;font-size:0.8em">${_hmpEsc(
        site,
      )} is normally sterile, so HMP has no healthy-abundance distribution for it: any detection is outside what a healthy subject carries (consider contamination / kit background via negative controls). Not faded.</span></div>`
    );
  }
  const ents = _hmpRefEntries(row);
  const n = _hmpNum(row["HMP Healthy Samples"]);
  const N = _hmpNum(row["HMP Reference Samples"]);
  const z = _hmpNum(row["HMP Z-Score"]);
  const pctl = _hmpNum(row["HMP Healthy Percentile"]);
  const mean = _hmpNum(row["HMP Healthy Mean %"]);
  const own = String(row["Taxonomic ID #"] || "");
  const via = ents.length && ents[0].taxid !== own ? ents[0] : null;
  let h = `<div style="max-width:430px"><b>${_hmpEsc(
    row["Detected Organism"],
  )}</b> · <span style="color:#ccc">${_hmpEsc(row["Specimen ID"])}</span><br>`;
  h += `<span style="display:inline-block;margin:3px 0;padding:0 6px;border-radius:3px;background:${col};color:#fff;font-weight:700;font-size:0.85em">${_hmpEsc(
    st,
  )}</span> <span style="color:#ccc">vs. healthy <b style="color:#fff">${_hmpEsc(site)}</b></span>`;
  if (via)
    h += `<br><span style="color:#ffd580;font-size:0.85em">No healthy reference at this taxid — compared at ${_hmpEsc(
      via.e.rank || "species",
    )} level: ${_hmpEsc(via.e.name || via.taxid)} (${_hmpEsc(via.taxid)}), using its abundance in this sample.</span>`;
  h += `<table style="border-collapse:collapse;font-size:0.85em;margin:5px 0">`;
  const tr = (k, v) =>
    `<tr><td style="color:#ffd580;padding-right:10px;white-space:nowrap">${k}</td><td>${v}</td></tr>`;
  h += tr("This sample", `${_hmpFmtPct(obs)} of reads`);
  if (st !== TT_HMP_ABSENT) {
    h += tr("Healthy mean ± sd", `${_hmpFmtPct(mean)} ± ${_hmpFmtPct(ents.reduce((s, x) => s + (x.e.s || 0), 0))}`);
    h += tr(
      "Z-score",
      `${isFinite(z) ? z.toFixed(2) : "—"} <span style="color:#aaa">(elevated at ≥ ${TT_HMP_Z})</span>`,
    );
  }
  if (isFinite(pctl))
    h += tr(
      "Healthy percentile",
      `${pctl.toFixed(
        pctl >= 99 ? 2 : 1,
      )} <span style="color:#aaa">— share of healthy samples at or below this abundance</span>`,
    );
  if (N > 0)
    h += tr(
      "Healthy prevalence",
      `${isFinite(n) ? Math.round(n).toLocaleString() : 0} of ${Math.round(
        N,
      ).toLocaleString()} reference samples (${_hmpFmtShare(n || 0, N)})`,
    );
  h += `</table>`;
  if (ents.length) h += _hmpHistSVG(ents, obs, col);
  h +=
    `<span style="color:#aaa;font-size:0.8em">Healthy distributions: HMP relative abundances per body site (mean/sd over samples where the organism was present). ` +
    (st === TT_HMP_WITHIN
      ? `Within the healthy range → less likely to be clinically significant; faded with ◆ as in the PDF report.`
      : st === TT_HMP_ABSENT
      ? `Never observed in healthy subjects at this site.`
      : `Above what healthy subjects at this site typically carry.`) +
    `</span></div>`;
  return h;
}

function _hmpHistSVG(ents, obs, col) {
  const b = (TT_HMP && TT_HMP.bins) || { lo: -6, hi: 2, step: 0.5 };
  const nb = Math.round((b.hi - b.lo) / b.step);
  const counts = new Array(nb).fill(0);
  ents.forEach(({ e }) => (e.h || []).forEach((c, i) => (counts[i] += c)));
  const max = Math.max(1, ...counts);
  const W = 400,
    H = 92,
    L = 28,
    R = 16,
    T = 6,
    B = 20,
    iw = W - L - R,
    ih = H - T - B,
    bw = iw / nb;
  let s = `<svg width="${W}" height="${H}" style="display:block;margin:2px 0 4px">`;
  counts.forEach((c, i) => {
    const bh = (c / max) * ih;
    if (c > 0)
      s += `<rect x="${L + i * bw + 0.5}" y="${T + ih - bh}" width="${
        bw - 1
      }" height="${bh}" fill="#74c0fc" opacity="0.85"><title>${c} healthy samples</title></rect>`;
  });
  s += `<line x1="${L}" y1="${T + ih}" x2="${W - R}" y2="${T + ih}" stroke="#868e96"/>`;
  for (let d = Math.ceil(b.lo); d <= b.hi; d += 2) {
    const x = L + ((d - b.lo) / (b.hi - b.lo)) * iw;
    const lab = d >= 0 ? Math.pow(10, d) + "%" : d >= -2 ? Math.pow(10, d).toFixed(-d) + "%" : "1e" + d + "%";
    s += `<line x1="${x}" y1="${T + ih}" x2="${x}" y2="${T + ih + 3}" stroke="#868e96"/>`;
    s += `<text x="${x}" y="${H - 5}" fill="#ced4da" font-size="9" text-anchor="middle">${lab}</text>`;
  }
  s += `<text x="${L - 4}" y="${T + 8}" fill="#ced4da" font-size="9" text-anchor="end">${max}</text>`;
  const ox = L + _hmpX(obs) * iw;
  s += `<line x1="${ox}" y1="${T - 2}" x2="${ox}" y2="${T + ih}" stroke="${col}" stroke-width="2"/>`;
  s += `<text x="${Math.min(W - R - 2, Math.max(L + 2, ox))}" y="${
    T + 8
  }" fill="${col}" font-size="9" font-weight="700" text-anchor="${ox > W - 60 ? "end" : "start"}"> this sample</text>`;
  s += `</svg><span style="color:#aaa;font-size:0.78em">Healthy samples by abundance (log scale). Samples where the organism was absent are not drawn.</span><br>`;
  return s;
}

/* Detections-table cell decoration for the "HMP Status" / "HMP Z-Score" columns. */
function ttHmpDecorateCell(td, row, col) {
  const st = row["HMP Status"];
  if (!st) return;
  if (col === "HMP Status") {
    const c = _hmpStatusColor(st);
    td.innerHTML = `<span style="display:inline-block;padding:0 5px;border-radius:3px;font-size:0.85em;font-weight:600;color:${c};border:1px solid ${c}55;background:${c}12;white-space:nowrap">${_hmpEsc(
      st,
    )}</span>`;
  } else if (col === "HMP Z-Score") {
    const z = _hmpNum(row[col]);
    if (isFinite(z)) td.textContent = Math.abs(z) >= 1000 ? z.toExponential(2) : z.toFixed(2);
  }
  td.style.cursor = "help";
  td.addEventListener("mouseover", (ev) => showTip(ttHmpTipHTML(row), ev));
  td.addEventListener("mousemove", moveTip);
  td.addEventListener("mouseout", hideTip);
}

/* Wire hover tooltips for markers / cells rendered as HTML strings. */
function ttHmpWireTips(root, resolveRow) {
  if (!root) return;
  root.querySelectorAll("[data-hmp-key]").forEach((el) => {
    const row = resolveRow(el.dataset.hmpKey);
    if (!row) return;
    el.addEventListener("mouseover", (ev) => {
      ev.stopPropagation();
      showTip(ttHmpTipHTML(row), ev);
    });
    el.addEventListener("mousemove", moveTip);
    el.addEventListener("mouseout", hideTip);
  });
}

/* ── Controls (injected only when the run carries HMP data) ──────────────── */
function _hmpControlHTML(idp) {
  return (
    `<label class="hmp-ctl" style="display:flex;align-items:center;gap:0.3em;white-space:nowrap" ` +
    `title="Fade organisms whose abundance is within the range seen in healthy subjects (HMP) at this body site — the PDF report's ◆ rows.">` +
    `<input type="checkbox" id="${idp}-hmp-fade" ${TT_HMP_FADE ? "checked" : ""}/> Fade healthy-range (HMP)` +
    `<span style="color:#868e96">z&lt;</span><input type="number" id="${idp}-hmp-z" value="${TT_HMP_Z}" step="0.5" ` +
    `style="width:3.6em;padding:0 0.2em;border:1px solid #ccc;border-radius:4px"/></label>`
  );
}
function _hmpSyncControls() {
  ["summary", "tbl"].forEach((p) => {
    const f = document.getElementById(p + "-hmp-fade");
    const z = document.getElementById(p + "-hmp-z");
    if (f) f.checked = TT_HMP_FADE;
    if (z && document.activeElement !== z) z.value = TT_HMP_Z;
  });
}
function _hmpRefresh() {
  ttHmpStampAll();
  _hmpSyncControls();
  if (typeof redraw === "function") redraw();
}
function ttHmpMountControls() {
  if (!ttHmpAvailable()) return;
  const spots = [
    ["summary", document.getElementById("summary-group-sample")],
    ["tbl", document.getElementById("tbl-group-sample")],
  ];
  spots.forEach(([p, anchor]) => {
    if (!anchor || document.getElementById(p + "-hmp-fade")) return;
    const lab = anchor.closest("label");
    if (!lab) return;
    lab.insertAdjacentHTML("beforebegin", _hmpControlHTML(p));
    document.getElementById(p + "-hmp-fade").addEventListener("change", (e) => {
      TT_HMP_FADE = !!e.target.checked;
      _hmpRefresh();
    });
    document.getElementById(p + "-hmp-z").addEventListener("change", (e) => {
      const v = parseFloat(e.target.value);
      if (!isFinite(v)) return;
      TT_HMP_Z = v;
      _hmpRefresh();
    });
  });
  // Legend entry next to the summary table's other badge explanations.
  const leg = document.getElementById("sum-level-legend-rollup");
  if (leg && !document.getElementById("sum-hmp-legend")) {
    const N = TT_HMP && TT_HMP.sites ? Object.entries(TT_HMP.sites) : [];
    const ref = N.length
      ? ` Reference: ${N.map(([s, n]) => `${Number(n).toLocaleString()} healthy ${_hmpEsc(s)}`).join(
          ", ",
        )} samples (HMP).`
      : "";
    leg.parentElement.insertAdjacentHTML(
      "beforeend",
      ` <span id="sum-hmp-legend"><span class="hmp-mark">◆<small> n (p%)</small>◆</span> = abundance within the ` +
        `healthy (HMP) range for the body site (z &lt; threshold): faded, n healthy reference samples carried it ` +
        `(p% of all samples for the site). The <b>Healthy (HMP)</b> column plots the healthy range ` +
        `(5–95% / IQR / median) against this sample's abundance (dot); <b style="color:#c92a2a">sterile site</b> = blood/CSF-type sample, ` +
        `where no organism is expected in health.${ref}</span>`,
    );
  }
}

(function _hmpInit() {
  if (!TT_HMP) return;
  ttHmpStampAll();
  const go = () => ttHmpMountControls();
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", go);
  else go();
})();
