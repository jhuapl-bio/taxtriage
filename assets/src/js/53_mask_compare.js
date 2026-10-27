/* ═══════════════════════════════════════════════════════════════════════════
       -  §  REGION MASKING — BEFORE / AFTER TASS   (match_paths --mask_compare)
       -     With --mask_compare (on by default when a mask is used) every masked
       -     sample is also scored WITHOUT the mask. Rows then carry
       -     "TASS (unmasked)" / "TASS Change (mask)" / "Reads (unmasked)" /
       -     "Breadth % (unmasked)" and SAMPLE_META[s].masking.compare holds the
       -     per-sample summary (+ organisms only detected without the mask).
       -
       -     Kept deliberately quiet: the masked TASS stays the headline number.
       -       • hover: Summary / Table TASS cells and heatmap cells
       -       • right panel "Masking: before / after": heatmap TASS mode
       -         (masked / unmasked / change), optional inline arrows in the
       -         Summary table, a one-line summary and a before/after table.
═══════════════════════════════════════════════════════════════════════════ */
let _ttMaskCmpModeVal = "masked";
let _ttMaskCmpShowTables = false;
let _ttMaskCmpAvail = null;
let _ttMaskCmpAvailN = -1;

function ttMaskCmpHas(r) {
  return !!r && r["TASS (unmasked)"] !== null && r["TASS (unmasked)"] !== undefined && r["TASS (unmasked)"] !== "";
}
function ttMaskCmpAvailable() {
  const n = typeof DATA !== "undefined" ? DATA.length : 0;
  if (_ttMaskCmpAvailN !== n) {
    _ttMaskCmpAvailN = n;
    _ttMaskCmpAvail = n > 0 && DATA.some(ttMaskCmpHas);
  }
  return _ttMaskCmpAvail;
}
function ttMaskCmpMode() {
  return ttMaskCmpAvailable() ? _ttMaskCmpModeVal : "masked";
}
function ttMaskCmpShowInTables() {
  return _ttMaskCmpShowTables && ttMaskCmpAvailable();
}
// Value of a row under a mode. Rows from samples scored without a mask have
// no unmasked field: their TASS is the same either way (change 0).
function ttMaskCmpValue(r, mode) {
  const after = num(r["TASS Score"]);
  if (!ttMaskCmpHas(r)) return mode === "delta" ? 0 : after;
  const before = num(r["TASS (unmasked)"]);
  return mode === "unmasked" ? before : mode === "delta" ? after - before : after;
}
function _mcSigned(v, d) {
  if (!isFinite(v)) return "—";
  const a = Math.abs(v).toFixed(d === undefined ? 1 : d);
  return (v > 0 ? "+" : v < 0 ? "−" : "±") + a;
}
function _mcCall(r, before, after) {
  const thr = typeof thresholdForRow === "function" ? thresholdForRow(r) : NaN;
  if (!isFinite(thr)) return null;
  const b = before >= thr,
    a = after >= thr;
  return { thr, b, a, changed: a !== b };
}
// Tooltip block (heatmap cells and other d3 tips).
function ttMaskCmpTipHTML(r) {
  if (!ttMaskCmpHas(r)) return "";
  const before = num(r["TASS (unmasked)"]),
    after = num(r["TASS Score"]);
  const c = _mcCall(r, before, after);
  const rb = num(r["Reads (unmasked)"]),
    ra = num(r["# Reads Aligned"]);
  const bb = num(r["Breadth % (unmasked)"]),
    ba = num(r["Breadth %"]);
  return (
    `<div style="margin-top:4px;padding-top:3px;border-top:1px solid rgba(255,255,255,.25);color:#cfd8dc;font-size:.92em">` +
    `Masking: TASS <b>${before.toFixed(1)}</b> → <b>${after.toFixed(1)}</b> (${_mcSigned(after - before)})` +
    `<br>reads ${_fmtInt(rb)} → ${_fmtInt(ra)} · breadth ${bb.toFixed(1)}% → ${ba.toFixed(1)}%` +
    (c && c.changed
      ? `<br><b style="color:${c.a ? "#a5d6a7" : "#ffab91"}">${
          c.a ? "passes only with" : "passes only without"
        } the mask</b> (cutoff ${c.thr})`
      : "") +
    `</div>`
  );
}
// Plain-text version for native title= tooltips (Summary / Table cells).
function ttMaskCmpTitle(r) {
  if (!ttMaskCmpHas(r)) return "";
  const before = num(r["TASS (unmasked)"]),
    after = num(r["TASS Score"]);
  const c = _mcCall(r, before, after);
  return (
    `Before masking: TASS ${before.toFixed(1)}  →  after: ${after.toFixed(1)} (${_mcSigned(after - before)})\n` +
    `Reads ${_fmtInt(num(r["Reads (unmasked)"]))} → ${_fmtInt(num(r["# Reads Aligned"]))} · ` +
    `breadth ${num(r["Breadth % (unmasked)"]).toFixed(1)}% → ${num(r["Breadth %"]).toFixed(1)}%` +
    (c && c.changed ? `\n${c.a ? "Passes only WITH" : "Passes only WITHOUT"} the mask (cutoff ${c.thr})` : "")
  );
}
// Small grey arrow after a TASS value (only with "Show in tables" on).
function ttMaskCmpInlineHTML(r) {
  if (!ttMaskCmpShowInTables() || !ttMaskCmpHas(r)) return "";
  const before = num(r["TASS (unmasked)"]),
    d = num(r["TASS Score"]) - before;
  const c = _mcCall(r, before, num(r["TASS Score"]));
  const col = c && c.changed ? (c.a ? "#2e7d32" : "#c62828") : "#8a8f98";
  const arrow = Math.abs(d) < 0.05 ? "=" : d > 0 ? "▴" : "▾";
  return (
    `<span class="mcmp-inline" style="color:${col};font-size:.78em;margin-left:3px;white-space:nowrap">` +
    `${arrow}${before.toFixed(1)}</span>`
  );
}

/* ── Rows for the before / after table ─────────────────────────────────── */
function _mcRows(level) {
  const lvl = level || "Strain";
  const out = [];
  const seen = new Set();
  (typeof DATA !== "undefined" ? DATA : []).forEach((r) => {
    if (!ttMaskCmpHas(r)) return;
    if ((r["Level"] || "Strain") !== lvl) return;
    const before = num(r["TASS (unmasked)"]),
      after = num(r["TASS Score"]);
    const c = _mcCall(r, before, after);
    seen.add(r["Specimen ID"]);
    out.push({
      sample: r["Specimen ID"],
      organism: r["Detected Organism"],
      level: lvl,
      before,
      after,
      delta: after - before,
      readsB: num(r["Reads (unmasked)"]),
      readsA: num(r["# Reads Aligned"]),
      brdB: num(r["Breadth % (unmasked)"]),
      brdA: num(r["Breadth %"]),
      notScored: r["Masked bp %"] === null || r["Masked bp %"] === undefined ? null : num(r["Masked bp %"]),
      call: c,
      removed: false,
    });
  });
  // organisms only detected without the mask (all their evidence masked away)
  const meta = (typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {};
  Object.keys(meta).forEach((s) => {
    const cmp = ((meta[s] || {}).masking || {}).compare;
    if (!cmp || !Array.isArray(cmp.dropped)) return;
    cmp.dropped.forEach((d) => {
      const dl = (d.level || "").toLowerCase();
      // a dropped genus/species hides everything below it, so list it at every level view
      const want = lvl.toLowerCase();
      const order = { genus: 0, species: 1, strain: 2 };
      if (order[dl] === undefined || order[dl] > order[want]) return;
      const ref = (typeof DATA !== "undefined" ? DATA : []).find((r) => r["Specimen ID"] === s) || {
        "Sample Type": (meta[s] || {}).sample_type,
      };
      const thr = typeof thresholdForRow === "function" ? thresholdForRow(ref) : NaN;
      out.push({
        sample: s,
        organism: d.name + (dl !== want ? ` (${dl})` : ""),
        level: lvl,
        before: num(d.tass_unmasked),
        after: null,
        delta: -num(d.tass_unmasked),
        readsB: num(d.reads_unmasked),
        readsA: 0,
        brdB: null,
        brdA: null,
        notScored: null,
        call: isFinite(thr)
          ? { thr, b: num(d.tass_unmasked) >= thr, a: false, changed: num(d.tass_unmasked) >= thr }
          : null,
        removed: true,
      });
    });
  });
  return out;
}

function _mcSummaryLine() {
  const rows = _mcRows("Strain");
  const real = rows.filter((x) => !x.removed);
  if (!real.length && !rows.length) return "";
  const eps = 0.5;
  const up = real.filter((x) => x.delta >= eps).length;
  const down = real.filter((x) => x.delta <= -eps).length;
  const mean = real.length ? real.reduce((s, x) => s + x.delta, 0) / real.length : 0;
  const gained = rows.filter((x) => x.call && x.call.changed && x.call.a).length;
  const lost = rows.filter((x) => x.call && x.call.changed && !x.call.a).length;
  const removed = rows.filter((x) => x.removed).length;
  return (
    `<b>${real.length}</b> strain${real.length === 1 ? "" : "s"} rescored without the mask: ` +
    `<span style="color:#2e7d32">${up} up</span>, <span style="color:#c62828">${down} down</span>, ` +
    `${real.length - up - down} within ±${eps} (mean ${_mcSigned(mean)}).` +
    `<br>Calls at the current cutoffs: <span style="color:#2e7d32">${gained} gained</span>, ` +
    `<span style="color:#c62828">${lost} lost</span>` +
    (removed ? ` · ${removed} only detected without the mask` : "") +
    "."
  );
}

function ttMaskCmpRefreshPanel() {
  const sec = document.getElementById("mcmp-section");
  if (!sec) return;
  const on = ttMaskCmpAvailable();
  sec.hidden = !on;
  if (!on) return;
  const line = document.getElementById("mcmp-summary");
  if (line) line.innerHTML = _mcSummaryLine();
  const sel = document.getElementById("mcmp-mode");
  if (sel && sel.value !== _ttMaskCmpModeVal) sel.value = _ttMaskCmpModeVal;
  const badge = document.getElementById("hm-mcmp-badge");
  if (badge) {
    const m = ttMaskCmpMode();
    badge.hidden = m === "masked";
    badge.textContent = m === "unmasked" ? "TASS before masking" : m === "delta" ? "TASS change from masking" : "";
  }
}

/* ── Before / after modal ───────────────────────────────────────────────── */
let _mcSort = { key: "absdelta", dir: -1 };
function _mcOpenModal() {
  let ov = document.getElementById("mcmp-modal");
  if (!ov) {
    ov = document.createElement("div");
    ov.id = "mcmp-modal";
    ov.className = "mcmp-modal";
    ov.innerHTML =
      `<div class="mcmp-dialog" role="dialog" aria-modal="true" aria-label="Masking before / after">` +
      `<div class="mcmp-head"><div><div class="mcmp-title">Region masking — TASS before / after</div>` +
      `<div class="mcmp-sub">Each masked sample was scored twice: without the mask (<b>before</b>) and with it (<b>after</b>, the TASS used everywhere else in this report).</div></div>` +
      `<button type="button" class="mcmp-x" aria-label="Close">×</button></div>` +
      `<div class="mcmp-samples"></div>` +
      `<div class="mcmp-ctrls">` +
      `<label>Level <select class="mcmp-level"><option>Strain</option><option>Species</option><option>Genus</option></select></label>` +
      `<label>Sample <select class="mcmp-sample"><option value="">All</option></select></label>` +
      `<label><input type="checkbox" class="mcmp-only"> Only changed (|change| ≥ 0.5, call changed or removed)</label>` +
      `<button type="button" class="mcmp-csv"><i class="fas fa-download"></i> CSV</button></div>` +
      `<div class="mcmp-tablewrap"><table class="mcmp-table"><thead></thead><tbody></tbody></table></div></div>`;
    document.body.appendChild(ov);
    ov.addEventListener("click", (e) => {
      if (e.target === ov || e.target.closest(".mcmp-x")) ov.classList.remove("show");
    });
    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape") ov.classList.remove("show");
    });
    ov.querySelectorAll(".mcmp-level,.mcmp-sample,.mcmp-only").forEach((el) =>
      el.addEventListener("change", _mcRenderModal),
    );
    ov.querySelector(".mcmp-csv").addEventListener("click", _mcCsv);
  }
  const ssel = ov.querySelector(".mcmp-sample");
  const samples = uniq(_mcRows("Genus").map((x) => x.sample)).sort();
  ssel.innerHTML =
    `<option value="">All</option>` + samples.map((s) => `<option value="${_esc(s)}">${_esc(s)}</option>`).join("");
  _mcRenderModal();
  ov.classList.add("show");
}

const _MC_COLS = [
  { key: "sample", label: "Sample" },
  { key: "organism", label: "Organism" },
  { key: "before", label: "TASS before", num: true, tip: "Scored without the mask" },
  { key: "after", label: "TASS after", num: true, tip: "Scored with the mask (the report's TASS)" },
  { key: "delta", label: "Change", num: true },
  {
    key: "call",
    label: "Call (before → after)",
    tip: "Pass / fail against the current TASS cutoff for the sample type",
  },
  { key: "readsB", label: "Reads before", num: true },
  { key: "readsA", label: "Reads after", num: true },
  { key: "brdB", label: "Breadth % before", num: true },
  { key: "brdA", label: "Breadth % after", num: true, tip: "Over the scored bases only" },
  { key: "notScored", label: "Not scored %", num: true },
];

function _mcFiltered() {
  const ov = document.getElementById("mcmp-modal");
  const lvl = ov.querySelector(".mcmp-level").value;
  const smp = ov.querySelector(".mcmp-sample").value;
  const only = ov.querySelector(".mcmp-only").checked;
  let rows = _mcRows(lvl).filter((x) => !smp || x.sample === smp);
  if (only) rows = rows.filter((x) => x.removed || Math.abs(x.delta) >= 0.5 || (x.call && x.call.changed));
  const k = _mcSort.key,
    dir = _mcSort.dir;
  const val = (x) =>
    k === "absdelta"
      ? Math.abs(x.delta)
      : k === "call"
      ? x.call
        ? (x.call.changed ? 2 : 0) + (x.call.a ? 1 : 0)
        : -1
      : x[k];
  rows.sort((a, b) => {
    const va = val(a),
      vb = val(b);
    if (va === vb) return 0;
    if (va === null || va === undefined) return 1;
    if (vb === null || vb === undefined) return -1;
    return (typeof va === "string" ? va.localeCompare(vb) : va - vb) * dir;
  });
  return rows;
}

function _mcCallTxt(c, removed) {
  if (!c) return removed ? "removed" : "—";
  const t = (p) => (p ? "pass" : "fail");
  return `${t(c.b)} → ${removed ? "removed" : t(c.a)}`;
}

function _mcRenderModal() {
  const ov = document.getElementById("mcmp-modal");
  if (!ov) return;
  // per-sample summary from the JSON metadata
  const meta = (typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {};
  const sRows = Object.keys(meta)
    .filter((s) => ((meta[s] || {}).masking || {}).compare)
    .sort()
    .map((s) => {
      const m = meta[s].masking,
        c = m.compare;
      const g = c.max_gain,
        l = c.max_loss;
      return (
        `<tr><td>${_esc(s)}</td><td>${_esc(m.mode || "exclude")}</td><td class="num">${c.strains_compared}</td>` +
        `<td class="num" style="color:#2e7d32">${c.strains_up}</td><td class="num" style="color:#c62828">${c.strains_down}</td>` +
        `<td class="num">${_mcSigned(c.mean_delta)}</td>` +
        `<td>${g ? `${_esc(g.name)} (${_mcSigned(g.delta)})` : "—"}</td>` +
        `<td>${l ? `${_esc(l.name)} (${_mcSigned(l.delta)})` : "—"}</td>` +
        `<td class="num">${c.n_dropped || 0}</td></tr>`
      );
    })
    .join("");
  ov.querySelector(".mcmp-samples").innerHTML = sRows
    ? `<table class="mcmp-table mcmp-stbl"><thead><tr><th>Sample</th><th>Mode</th><th>Strains</th><th>Up</th><th>Down</th>` +
      `<th>Mean change</th><th>Largest gain</th><th>Largest loss</th><th title="Detected without the mask only">Removed</th></tr></thead>` +
      `<tbody>${sRows}</tbody></table>`
    : "";

  const rows = _mcFiltered();
  const head = ov.querySelector(".mcmp-table:not(.mcmp-stbl) thead");
  head.innerHTML =
    "<tr>" +
    _MC_COLS
      .map(
        (c) =>
          `<th data-k="${c.key}" class="${c.num ? "num" : ""}"${c.tip ? ` title="${_esc(c.tip)}"` : ""}>${c.label}${
            _mcSort.key === c.key ? (_mcSort.dir > 0 ? " ▴" : " ▾") : ""
          }</th>`,
      )
      .join("") +
    "</tr>";
  head.querySelectorAll("th").forEach((th) =>
    th.addEventListener("click", () => {
      const k = th.dataset.k;
      _mcSort = { key: k, dir: _mcSort.key === k ? -_mcSort.dir : k === "sample" || k === "organism" ? 1 : -1 };
      _mcRenderModal();
    }),
  );
  const f1 = (v) => (v === null || v === undefined ? "—" : (+v).toFixed(1));
  ov.querySelector(".mcmp-table:not(.mcmp-stbl) tbody").innerHTML = rows.length
    ? rows
        .map((x) => {
          const cc = x.call && x.call.changed ? (x.call.a ? "mc-gain" : "mc-loss") : x.removed ? "mc-loss" : "";
          const dc = x.delta >= 0.5 ? "#2e7d32" : x.delta <= -0.5 ? "#c62828" : "#666";
          return (
            `<tr class="${x.removed ? "mc-removed" : ""}"><td>${_esc(x.sample)}</td><td><i>${_esc(
              x.organism,
            )}</i></td>` +
            `<td class="num">${f1(x.before)}</td><td class="num"><b>${x.removed ? "—" : f1(x.after)}</b></td>` +
            `<td class="num" style="color:${dc}">${_mcSigned(x.delta)}</td>` +
            `<td class="${cc}">${_mcCallTxt(x.call, x.removed)}</td>` +
            `<td class="num">${_fmtInt(x.readsB)}</td><td class="num">${_fmtInt(x.readsA)}</td>` +
            `<td class="num">${f1(x.brdB)}</td><td class="num">${f1(x.brdA)}</td>` +
            `<td class="num">${x.notScored === null ? "—" : x.notScored.toFixed(1)}</td></tr>`
          );
        })
        .join("")
    : `<tr><td colspan="${_MC_COLS.length}" style="color:#888;padding:10px">No rows.</td></tr>`;
}

function _mcCsv() {
  const rows = _mcFiltered();
  const q = (v) => {
    const s = v === null || v === undefined ? "" : String(v);
    return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  };
  const lines = [
    [
      "sample",
      "organism",
      "level",
      "tass_before",
      "tass_after",
      "change",
      "call_before",
      "call_after",
      "cutoff",
      "reads_before",
      "reads_after",
      "breadth_pct_before",
      "breadth_pct_after",
      "not_scored_pct",
    ].join(","),
  ];
  rows.forEach((x) =>
    lines.push(
      [
        x.sample,
        x.organism,
        x.level,
        x.before,
        x.removed ? "" : x.after,
        x.delta.toFixed(4),
        x.call ? (x.call.b ? "pass" : "fail") : "",
        x.call ? (x.removed ? "removed" : x.call.a ? "pass" : "fail") : x.removed ? "removed" : "",
        x.call ? x.call.thr : "",
        x.readsB,
        x.readsA,
        x.brdB,
        x.brdA,
        x.notScored,
      ]
        .map(q)
        .join(","),
    ),
  );
  const blob = new Blob([lines.join("\n") + "\n"], { type: "text/csv" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = "masking_before_after.csv";
  document.body.appendChild(a);
  a.click();
  setTimeout(() => {
    URL.revokeObjectURL(a.href);
    a.remove();
  }, 0);
}

/* ── wiring ────────────────────────────────────────────────────────────── */
(function () {
  const sel = document.getElementById("mcmp-mode");
  if (sel)
    sel.addEventListener("change", () => {
      _ttMaskCmpModeVal = sel.value;
      const hv = document.getElementById("hm-value-col");
      // the mode applies to the TASS value; switch the heatmap to it
      if (_ttMaskCmpModeVal !== "masked" && hv && hv.value !== "TASS Score") hv.value = "TASS Score";
      ttMaskCmpRefreshPanel();
      if (typeof redraw === "function") redraw();
    });
  const cb = document.getElementById("mcmp-show-tables");
  if (cb)
    cb.addEventListener("change", () => {
      _ttMaskCmpShowTables = cb.checked;
      if (typeof redraw === "function") redraw();
    });
  const btn = document.getElementById("mcmp-open-btn");
  if (btn) btn.addEventListener("click", _mcOpenModal);
  // SAMPLE_META / DATA arrive during init (or a later file drop): refresh then
  window.addEventListener("load", () => setTimeout(ttMaskCmpRefreshPanel, 0));
})();
