/* ═══════════════════════════════════════════════════════════════════════════
       -  §  ADAPTIVE TASS vs TASS (default)   (match_paths --adaptive_tass)
       -     TASS (default) is always computed and is the report's score unless
       -     the user switches the report to adaptive TASS. With --adaptive_tass
       -     (off by default) a second, adaptive score is computed per organism,
       -     with the weights shifted by the organism's context in the sample.
       -     Rows then carry "TASS (adaptive)" / "TASS Change (adaptive)" /
       -     "Adaptive Penalty" / "Adaptive Cutoff Shift" / "Adaptive Reasons",
       -     and SAMPLE_META[s].adaptive_tass holds the per-sample settings and
       -     change statistics.
       -
       -       • right panel "Adaptive TASS":
       -           - "TASS used in report": TASS (default) | Adaptive TASS.
       -             Adaptive swaps "TASS Score" (and the roll-up "Species TASS"
       -             / "Genus TASS") in every row that has an adaptive score, so
       -             every tab, filter, cutoff check, roll-up and export follows.
       -             Cutoffs themselves do not change. The default values are
       -             kept on each row (non-enumerable) and put back on switch
       -             back and while a session file is written (29_session_state).
       -           - heatmap change view (adaptive − default), optional inline
       -             other-score chips in the Summary table, a one-line summary
       -             and a comparison window (settings, scatter, table, CSV).
       -       • hover: Summary / Table TASS cells and heatmap cells show both.
       -       • banner pill while adaptive TASS is in use.
═══════════════════════════════════════════════════════════════════════════ */
let _ttAdaptScoreVal = "default"; // "default" | "adaptive": the TASS the whole report uses
let _ttAdaptHmDelta = false; // heatmap shows adaptive − default
let _ttAdaptShowTables = false;
let _ttAdaptAvail = null;
let _ttAdaptAvailN = -1;
let _ttAdaptAvailRef = null;
const _AD_EPS = 0.5; // TASS points: smaller changes count as "unchanged"
const _AD_SWAP = "__adDefault"; // per-row default TASS fields while adaptive is in use
const _AD_SWAP_COLS = ["TASS Score", "Species TASS", "Genus TASS"];

function ttAdaptHas(r) {
  const v = r ? r["TASS (adaptive)"] : null;
  return v !== null && v !== undefined && v !== "";
}
function ttAdaptAvailable() {
  const d = typeof DATA !== "undefined" ? DATA : null;
  const n = d ? d.length : 0;
  if (_ttAdaptAvailN !== n || _ttAdaptAvailRef !== d) {
    _ttAdaptAvailN = n;
    _ttAdaptAvailRef = d;
    _ttAdaptAvail = n > 0 && d.some(ttAdaptHas);
  }
  return _ttAdaptAvail;
}
/** True while adaptive TASS is the report's score. */
function ttAdaptScoreActive() {
  return _ttAdaptScoreVal === "adaptive" && ttAdaptAvailable();
}
/** Heatmap TASS view: "standard" (the report's score as is), "adaptive" (the
 *  report uses adaptive TASS) or "delta" (adaptive − default). */
function ttAdaptMode() {
  if (!ttAdaptAvailable()) return "standard";
  if (_ttAdaptHmDelta) return "delta";
  return _ttAdaptScoreVal === "adaptive" ? "adaptive" : "standard";
}
function ttAdaptShowInTables() {
  return _ttAdaptShowTables && ttAdaptAvailable();
}
/** TASS (default) of a row, whichever score the report is using. */
function ttTassDefault(r) {
  if (!r) return 0;
  const saved = r[_AD_SWAP];
  return saved && "TASS Score" in saved ? num(saved["TASS Score"]) : num(r["TASS Score"]);
}
function _adStd(r) {
  return ttTassDefault(r);
}
function _adAd(r) {
  return num(r["TASS (adaptive)"]);
}
// Value of a row under a mode. Rows from samples run without --adaptive_tass
// have no adaptive field: their TASS is the default one (change 0).
function ttAdaptValue(r, mode) {
  const std = _adStd(r);
  if (!ttAdaptHas(r)) return mode === "delta" ? 0 : std;
  const ad = _adAd(r);
  return mode === "adaptive" ? ad : mode === "delta" ? ad - std : std;
}

/* ── Report-wide switch: swap the score columns in DATA ─────────────────── */
function _adParentIndex(rows) {
  const sp = new Map(),
    ge = new Map();
  rows.forEach((r) => {
    if (!ttAdaptHas(r)) return;
    const lv = r["Level"] || "Strain",
      s = r["Specimen ID"];
    if (lv === "Species") {
      sp.set(`${s}\u0001${r["Subkey"] == null ? "" : r["Subkey"]}`, r);
      sp.set(`${s}\u0002${r["Detected Organism"] || ""}`, r);
    } else if (lv === "Genus") {
      ge.set(`${s}\u0002${r["Detected Organism"] || ""}`, r);
    }
  });
  return { sp, ge };
}
function _adRestoreRow(r) {
  const saved = r[_AD_SWAP];
  if (!saved) return;
  _AD_SWAP_COLS.forEach((c) => {
    if (c in saved) r[c] = saved[c];
  });
  delete r[_AD_SWAP];
}
// adaptive=true: each row with an adaptive score takes it as "TASS Score"; its
// roll-up Species / Genus TASS take the parent rows' adaptive scores (same
// sample). adaptive=false: every swapped row gets its default values back.
function _adApplyRows(rows, adaptive) {
  if (!Array.isArray(rows)) return;
  if (!adaptive) {
    rows.forEach((r) => r && _adRestoreRow(r));
    return;
  }
  const idx = _adParentIndex(rows);
  rows.forEach((r) => {
    if (!r || !ttAdaptHas(r)) return;
    let saved = r[_AD_SWAP];
    if (!saved) {
      saved = {};
      _AD_SWAP_COLS.forEach((c) => {
        if (c in r) saved[c] = r[c];
      });
      // Non-enumerable: never shown as a column, never serialised.
      Object.defineProperty(r, _AD_SWAP, { value: saved, writable: true, configurable: true, enumerable: false });
    }
    const s = r["Specimen ID"];
    r["TASS Score"] = _adAd(r);
    if ("Species TASS" in saved) {
      const p =
        idx.sp.get(`${s}\u0001${r["Subkey"] == null ? "" : r["Subkey"]}`) ||
        idx.sp.get(`${s}\u0002${r["Species Name"] || ""}`);
      r["Species TASS"] = p ? _adAd(p) : saved["Species TASS"];
    }
    if ("Genus TASS" in saved) {
      const p = idx.ge.get(`${s}\u0002${r["Genus Name"] || r["Genus"] || ""}`);
      r["Genus TASS"] = p ? _adAd(p) : saved["Genus TASS"];
    }
  });
}
let _adAppliedMode = "default";
let _adAppliedRef = null;
let _adAppliedLen = -1;
/** Make DATA match the chosen score. Cheap when nothing changed; called at the
 *  top of every redraw() (21_compare_coverages.js) so new or reloaded data
 *  follows the switch too. */
function ttAdaptEnsure() {
  if (typeof DATA === "undefined") return;
  const want = ttAdaptScoreActive() ? "adaptive" : "default";
  if (want === _adAppliedMode && DATA === _adAppliedRef && DATA.length === _adAppliedLen) return;
  _adApplyRows(DATA, want === "adaptive");
  _adAppliedMode = want;
  _adAppliedRef = DATA;
  _adAppliedLen = DATA.length;
  if (typeof _invalidateFilterCache === "function") _invalidateFilterCache();
}
/** Put the default values back (e.g. while a session file is written).
 *  ttAdaptEnsure() — or the next redraw — re-applies the chosen score. */
function ttAdaptSuspend() {
  if (typeof DATA === "undefined" || _adAppliedMode !== "adaptive") return false;
  _adApplyRows(DATA, false);
  _adAppliedMode = "default";
  return true;
}
function ttAdaptResume() {
  ttAdaptEnsure();
}
// Session state (29_session_state.js)
function ttAdaptCaptureConfig() {
  return { score: _ttAdaptScoreVal, heatmapChange: _ttAdaptHmDelta, showInTables: _ttAdaptShowTables };
}
function ttAdaptLoadConfig(cfg) {
  if (!cfg || typeof cfg !== "object") return;
  _ttAdaptScoreVal = cfg.score === "adaptive" ? "adaptive" : "default";
  _ttAdaptHmDelta = !!cfg.heatmapChange;
  _ttAdaptShowTables = !!cfg.showInTables;
  ttAdaptEnsure();
  ttAdaptRefreshPanel();
}

/* ── formatting helpers ─────────────────────────────────────────────────── */
function _adSigned(v, d) {
  if (v === null || v === undefined || !isFinite(v)) return "—";
  const a = Math.abs(v).toFixed(d === undefined ? 1 : d);
  return (v > 0 ? "+" : v < 0 ? "−" : "±") + a;
}
// Pass / fail of both scores against the current cutoff for the row's sample type.
function _adCall(r, std, ad) {
  const thr = typeof thresholdForRow === "function" ? thresholdForRow(r) : NaN;
  if (!isFinite(thr)) return null;
  const s = std >= thr,
    a = ad >= thr;
  return { thr, s, a, changed: s !== a };
}
function _adCallTxt(c) {
  if (!c) return "—";
  const t = (p) => (p ? "pass" : "fail");
  return `${t(c.s)} → ${t(c.a)}`;
}
function _adReasonList(r) {
  return String((r && r["Adaptive Reasons"]) || "")
    .split(/;\s*/)
    .map((x) => x.trim())
    .filter(Boolean);
}
function _adNumOrNull(v) {
  if (v === null || v === undefined || v === "") return null;
  const n = parseFloat(v);
  return isNaN(n) ? null : n;
}
const _AD_INUSE = " (in use)";

// Tooltip block (heatmap cells and other d3 tips).
function ttAdaptTipHTML(r) {
  if (!ttAdaptHas(r)) return "";
  const std = _adStd(r),
    ad = _adAd(r);
  const on = ttAdaptScoreActive();
  const c = _adCall(r, std, ad);
  const pen = _adNumOrNull(r["Adaptive Penalty"]);
  const rs = _adReasonList(r);
  return (
    `<div style="margin-top:4px;padding-top:3px;border-top:1px solid rgba(255,255,255,.25);color:#cfd8dc;font-size:.92em">` +
    `Adaptive TASS <b>${ad.toFixed(1)}</b>${on ? _AD_INUSE : ""} · TASS (default) <b>${std.toFixed(1)}</b>${
      on ? "" : _AD_INUSE
    } (${_adSigned(ad - std)})` +
    (pen !== null && pen < 0.9999 ? ` · ambiguity ×${pen.toFixed(2)}` : "") +
    (c && c.changed
      ? `<br><b style="color:${c.a ? "#a5d6a7" : "#ffab91"}">Call: default ${c.s ? "pass" : "fail"} → adaptive ${
          c.a ? "pass" : "fail"
        }</b> (cutoff ${c.thr})`
      : "") +
    (rs.length
      ? `<br><span style="opacity:.85">${rs
          .slice(0, 2)
          .map((x) => _esc(x))
          .join("<br>")}${rs.length > 2 ? ` <i>(+${rs.length - 2} more)</i>` : ""}</span>`
      : "") +
    `</div>`
  );
}
// Plain-text version for native title= tooltips (Summary / Table cells).
function ttAdaptTitle(r) {
  if (!ttAdaptHas(r)) return "";
  const std = _adStd(r),
    ad = _adAd(r);
  const on = ttAdaptScoreActive();
  const c = _adCall(r, std, ad);
  const pen = _adNumOrNull(r["Adaptive Penalty"]);
  const rs = _adReasonList(r);
  return (
    `Adaptive TASS ${ad.toFixed(1)}${on ? _AD_INUSE : ""}  ·  TASS (default) ${std.toFixed(1)}${
      on ? "" : _AD_INUSE
    }  (${_adSigned(ad - std)})` +
    (pen !== null && pen < 0.9999 ? `\nAmbiguity mark-down ×${pen.toFixed(2)}` : "") +
    (c && c.changed
      ? `\nCall: default ${c.s ? "PASS" : "FAIL"} → adaptive ${c.a ? "PASS" : "FAIL"} (cutoff ${c.thr})`
      : "") +
    (rs.length ? `\n• ${rs.slice(0, 3).join("\n• ")}${rs.length > 3 ? `\n(+${rs.length - 3} more)` : ""}` : "")
  );
}
// Join helper for cells that may also carry the masking title.
function ttAdaptJoinTitles(a, b) {
  return [a, b].filter(Boolean).join("\n\n");
}
// Small chip after a TASS with the score NOT in use (only with "Show in tables"
// on): A = adaptive while the report uses the default, D = default while it
// uses adaptive. Green / red when that score would flip the call.
function ttAdaptInlineHTML(r) {
  if (!ttAdaptShowInTables() || !ttAdaptHas(r)) return "";
  const std = _adStd(r),
    ad = _adAd(r);
  const on = ttAdaptScoreActive();
  const other = on ? std : ad,
    shown = on ? ad : std;
  const d = other - shown;
  const c = _adCall(r, std, ad);
  const otherPass = c ? (on ? c.s : c.a) : null;
  const col = c && c.changed ? (otherPass ? "#2e7d32" : "#c62828") : "#7b6fa8";
  const arrow = Math.abs(d) < 0.05 ? "=" : d > 0 ? "▴" : "▾";
  return (
    `<span class="adt-inline" style="color:${col}" title="${on ? "TASS (default)" : "Adaptive TASS"}">` +
    `${on ? "D" : "A"}${arrow}${other.toFixed(1)}</span>`
  );
}

/* ── Rows for the comparison ─────────────────────────────────────────── */
function _adRows(level) {
  const lvl = level || "Strain";
  const out = [];
  (typeof DATA !== "undefined" ? DATA : []).forEach((r) => {
    if (!ttAdaptHas(r)) return;
    if ((r["Level"] || "Strain") !== lvl) return;
    const std = _adStd(r),
      ad = _adAd(r);
    out.push({
      row: r,
      sample: r["Specimen ID"],
      organism: r["Detected Organism"],
      level: lvl,
      std,
      ad,
      delta: ad - std,
      penalty: _adNumOrNull(r["Adaptive Penalty"]),
      cutShift: _adNumOrNull(r["Adaptive Cutoff Shift"]),
      reasons: _adReasonList(r),
      call: _adCall(r, std, ad),
    });
  });
  return out;
}

function _adMetaSamples() {
  const meta = (typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {};
  return Object.keys(meta)
    .filter((s) => (meta[s] || {}).adaptive_tass)
    .sort()
    .map((s) => ({ sample: s, m: meta[s].adaptive_tass }));
}

function _adSummaryLine() {
  const rows = _adRows("Strain");
  if (!rows.length) return "";
  const on = ttAdaptScoreActive();
  const up = rows.filter((x) => x.delta >= _AD_EPS).length;
  const down = rows.filter((x) => x.delta <= -_AD_EPS).length;
  const mean = rows.reduce((s, x) => s + x.delta, 0) / rows.length;
  const gained = rows.filter((x) => x.call && x.call.changed && x.call.a).length;
  const lost = rows.filter((x) => x.call && x.call.changed && !x.call.a).length;
  const fallback = _adMetaSamples().filter((x) => x.m.applied === false);
  const pl = (n) => (n === 1 ? "" : "s");
  return (
    (on
      ? `<b>Adaptive TASS</b> is in use across the report (same cutoffs). Against TASS (default), `
      : `<b>TASS (default)</b> is the report's score. `) +
    `<b>${rows.length}</b> strain${pl(rows.length)}${on ? ":" : " also have an adaptive score:"} ` +
    `<span style="color:#2e7d32">${up} higher</span>, ` +
    `<span style="color:#c62828">${down} lower</span>, ${rows.length - up - down} within ±${_AD_EPS} ` +
    `(mean ${_adSigned(mean)}).` +
    `<br>At the current cutoffs adaptive TASS ${on ? "" : "would "}` +
    `<span style="color:#2e7d32">add${on ? "s" : ""} ${gained}</span> and ` +
    `<span style="color:#c62828">drop${on ? "s" : ""} ${lost}</span> strain call${pl(gained + lost)}.` +
    (fallback.length
      ? `<br><span style="color:#8d6e63">${fallback.length} sample${pl(fallback.length)} kept default weights (${_esc(
          fallback[0].m.reason || "not applied",
        )}${fallback.length > 1 ? ", …" : ""}).</span>`
      : "")
  );
}

function ttAdaptRefreshPanel() {
  const sec = document.getElementById("adt-section");
  const pill = document.getElementById("banner-adt");
  const on = ttAdaptAvailable();
  const active = ttAdaptScoreActive();
  if (document.body) document.body.classList.toggle("tt-adaptive-tass", active);
  if (pill) pill.hidden = !active;
  if (!sec) return;
  sec.hidden = !on;
  const badge = document.getElementById("hm-adt-badge");
  if (!on) {
    if (badge) badge.hidden = true;
    return;
  }
  const line = document.getElementById("adt-summary");
  if (line) line.innerHTML = _adSummaryLine();
  const sel = document.getElementById("adt-score");
  if (sel && sel.value !== _ttAdaptScoreVal) sel.value = _ttAdaptScoreVal;
  const hd = document.getElementById("adt-hm-delta");
  if (hd) hd.checked = _ttAdaptHmDelta;
  const st = document.getElementById("adt-show-tables");
  if (st) st.checked = _ttAdaptShowTables;
  const stl = document.getElementById("adt-show-tables-lbl");
  if (stl)
    stl.textContent = active ? "Show TASS (default) beside TASS in tables" : "Show adaptive TASS beside TASS in tables";
  if (badge) {
    const m = ttAdaptMode();
    badge.hidden = m === "standard";
    badge.textContent = m === "adaptive" ? "Adaptive TASS" : m === "delta" ? "TASS change (adaptive − default)" : "";
  }
}

/* ── Comparison window ──────────────────────────────────────────────── */
let _adSort = { key: "absdelta", dir: -1 };

const _AD_EXPLAIN =
  `<details class="adt-explain"><summary>How adaptive TASS differs from TASS (default)</summary>` +
  `<p><b>TASS (default)</b> scores every organism in a sample with one weight vector (for the sample type ` +
  `and platform, or from the CLI). It is always computed and is the report's score unless you switch ` +
  `<b>TASS used in report</b> (right panel, Adaptive TASS) to adaptive.</p>` +
  `<p><b>Adaptive TASS</b> (<code>--adaptive_tass</code>, off by default) starts from the same weights and ` +
  `shifts them per organism according to its context in the sample, then re-normalises them:</p><ul>` +
  `<li><b>complexity</b>: rich, even communities put more weight on MinHash and Gini.</li>` +
  `<li><b>ani</b>: another species with ANI at or above the cutoff is also detected: more MinHash and Gini, ` +
  `less breadth, and an ambiguity mark-down (needs <code>--enable_matrix</code>).</li>` +
  `<li><b>sibling</b>: other species of the same genus are detected: more MinHash, less breadth, and a ` +
  `smaller mark-down.</li>` +
  `<li><b>low_abundance</b>: under 3% of the reads at its level: more MinHash and Gini.</li>` +
  `<li><b>rollup</b>: species / genus built from two or more children: more breadth and Gini, less MinHash.</li>` +
  `</ul><p>The ambiguity mark-down is <code>1 − max_penalty × strength × ambiguity × (1 − MinHash)</code>, ` +
  `so a call is only marked down when its own unique (MinHash) support is weak. Strains of the same ` +
  `species never count against each other. <b>Cutoff shift</b> is the adaptive cutoff change the model ` +
  `suggests, in TASS points; it is reported only and never applied, so switching the report to adaptive ` +
  `TASS keeps the same cutoffs. The weights are heuristics: calibrate them on labelled data before ` +
  `relying on adaptive calls.</p></details>`;

function _adModalSub() {
  return ttAdaptScoreActive()
    ? `<b>Adaptive TASS</b> is in use across the report (same cutoffs). The call column shows what changes against <b>TASS (default)</b>.`
    : `<b>TASS (default)</b> is the report's score. <b>Adaptive TASS</b> is an extra score computed with ` +
        `<code>--adaptive_tass</code>; the call column shows what it would change at the current cutoffs.`;
}

function _adOpenModal() {
  let ov = document.getElementById("adt-modal");
  if (!ov) {
    ov = document.createElement("div");
    ov.id = "adt-modal";
    ov.className = "mcmp-modal";
    ov.innerHTML =
      `<div class="mcmp-dialog" role="dialog" aria-modal="true" aria-label="Adaptive TASS vs TASS (default)">` +
      `<div class="mcmp-head"><div><div class="mcmp-title">Adaptive TASS vs TASS (default)</div>` +
      `<div class="mcmp-sub adt-sub"></div></div>` +
      `<button type="button" class="mcmp-x" aria-label="Close">×</button></div>` +
      `<div class="adt-body">` +
      _AD_EXPLAIN +
      `<div class="adt-samples"></div>` +
      `<div class="mcmp-ctrls">` +
      `<label>Level <select class="adt-level"><option>Strain</option><option>Species</option><option>Genus</option></select></label>` +
      `<label>Sample <select class="adt-sample"><option value="">All</option></select></label>` +
      `<label><input type="checkbox" class="adt-only"> Only changed (|change| ≥ ${_AD_EPS} or call changed)</label>` +
      `<button type="button" class="mcmp-csv adt-csv"><i class="fas fa-download"></i> CSV</button></div>` +
      `<div class="adt-plotrow"><div class="adt-plot"></div><div class="adt-plotnote"></div></div>` +
      `<div class="mcmp-tablewrap"><table class="mcmp-table adt-orgs"><thead></thead><tbody></tbody></table></div>` +
      `</div></div>`;
    document.body.appendChild(ov);
    ov.addEventListener("click", (e) => {
      if (e.target === ov || e.target.closest(".mcmp-x")) ov.classList.remove("show");
    });
    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape") ov.classList.remove("show");
    });
    ov.querySelectorAll(".adt-level,.adt-sample,.adt-only").forEach((el) =>
      el.addEventListener("change", _adRenderModal),
    );
    ov.querySelector(".adt-csv").addEventListener("click", _adCsv);
  }
  const ssel = ov.querySelector(".adt-sample");
  const samples = uniq(["Strain", "Species", "Genus"].flatMap((l) => _adRows(l).map((x) => x.sample))).sort();
  ssel.innerHTML =
    `<option value="">All</option>` + samples.map((s) => `<option value="${_esc(s)}">${_esc(s)}</option>`).join("");
  _adRenderModal();
  ov.classList.add("show");
}

function _adCols() {
  const on = ttAdaptScoreActive();
  return [
    { key: "sample", label: "Sample" },
    { key: "organism", label: "Organism" },
    {
      key: "std",
      label: "TASS (default)",
      num: true,
      tip: on
        ? "Default weights (not in use: the report uses adaptive TASS)"
        : "The report's TASS (used for every cutoff)",
    },
    {
      key: "ad",
      label: "Adaptive TASS",
      num: true,
      tip: on ? "--adaptive_tass score (in use across the report)" : "--adaptive_tass score",
    },
    { key: "delta", label: "Change", num: true, tip: "Adaptive − default, in TASS points" },
    {
      key: "call",
      label: "Call (default → adaptive)",
      tip: "Pass / fail of each score against the current TASS cutoff for the sample type",
    },
    { key: "penalty", label: "Ambiguity ×", num: true, tip: "Multiplicative mark-down (1 = none)" },
    { key: "cutShift", label: "Cutoff shift", num: true, tip: "Suggested cutoff change in TASS points (not applied)" },
    { key: "reasons", label: "Why" },
  ];
}

function _adFiltered() {
  const ov = document.getElementById("adt-modal");
  const lvl = ov.querySelector(".adt-level").value;
  const smp = ov.querySelector(".adt-sample").value;
  const only = ov.querySelector(".adt-only").checked;
  let rows = _adRows(lvl).filter((x) => !smp || x.sample === smp);
  if (only) rows = rows.filter((x) => Math.abs(x.delta) >= _AD_EPS || (x.call && x.call.changed));
  const k = _adSort.key,
    dir = _adSort.dir;
  const val = (x) =>
    k === "absdelta"
      ? Math.abs(x.delta)
      : k === "call"
      ? x.call
        ? (x.call.changed ? 2 : 0) + (x.call.a ? 1 : 0)
        : -1
      : k === "reasons"
      ? x.reasons.join("; ")
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

function _adSampleTable() {
  const rows = _adMetaSamples()
    .map(({ sample, m }) => {
      const st = (m.levels || {}).strain || {};
      const smp = m.sample || {};
      const status =
        m.applied === false
          ? `<span style="color:#8d6e63" title="${_esc(m.reason || "")}">default weights kept</span>`
          : `applied${
              m.errors ? ` <span style="color:#c62828">(${m.errors} error${m.errors === 1 ? "" : "s"})</span>` : ""
            }`;
      const cx = _adNumOrNull(smp.complexity);
      return (
        `<tr><td>${_esc(sample)}</td><td>${status}</td>` +
        `<td class="num">${m.strength !== undefined ? m.strength : "—"}</td>` +
        `<td>${_esc((m.terms || []).join(", ") || "none")}</td><td>${_esc(m.base || "sampletype")}</td>` +
        `<td class="num">${cx === null ? "—" : cx.toFixed(2)}</td>` +
        `<td>${m.ani_available ? "yes" : '<span style="color:#888">no</span>'}</td>` +
        `<td class="num">${st.organisms !== undefined ? `${st.changed}/${st.organisms}` : "—"}</td>` +
        `<td class="num" style="color:#2e7d32">${st.raised !== undefined ? st.raised : "—"}</td>` +
        `<td class="num" style="color:#c62828">${st.lowered !== undefined ? st.lowered : "—"}</td>` +
        `<td class="num">${st.mean_delta !== undefined ? _adSigned(st.mean_delta, 2) : "—"}</td>` +
        `<td class="num">${st.max_abs_delta !== undefined ? (+st.max_abs_delta).toFixed(2) : "—"}</td></tr>`
      );
    })
    .join("");
  return rows
    ? `<table class="mcmp-table mcmp-stbl"><thead><tr><th>Sample</th><th>Status</th><th>Strength</th><th>Terms</th>` +
        `<th>Base</th><th title="Species-level richness × evenness, 0-1">Complexity</th><th title="--enable_matrix ANI available">ANI</th>` +
        `<th title="Strains whose adaptive TASS differs from TASS (default) by any amount">Strains changed</th>` +
        `<th title="Adaptive higher than default, by any amount">Higher (any)</th>` +
        `<th title="Adaptive lower than default, by any amount">Lower (any)</th>` +
        `<th title="Mean adaptive − default over all strains, TASS points">Mean change</th><th>Max |change|</th></tr></thead><tbody>${rows}</tbody></table>`
    : "";
}

// Default (x) vs adaptive (y) scatter: points off the diagonal moved; the
// shaded quadrants are calls that flip at the cutoff.
function _adDrawPlot(host, rows) {
  host.innerHTML = "";
  if (typeof d3 === "undefined" || !rows.length) return;
  const W = 330,
    H = 270,
    m = { l: 40, r: 10, t: 10, b: 34 };
  const iw = W - m.l - m.r,
    ih = H - m.t - m.b;
  const svg = d3.select(host).append("svg").attr("width", W).attr("height", H).attr("class", "adt-svg");
  const g = svg.append("g").attr("transform", `translate(${m.l},${m.t})`);
  const x = d3.scaleLinear().domain([0, 100]).range([0, iw]);
  const y = d3.scaleLinear().domain([0, 100]).range([ih, 0]);
  const thrs = uniq(rows.map((r) => (r.call ? r.call.thr : null)).filter((t) => t !== null && isFinite(t)));
  if (thrs.length === 1) {
    const t = thrs[0];
    g.append("rect")
      .attr("x", x(t))
      .attr("y", y(t))
      .attr("width", iw - x(t))
      .attr("height", ih - y(t))
      .attr("fill", "#fdecea")
      .attr("opacity", 0.7);
    g.append("rect")
      .attr("x", 0)
      .attr("y", 0)
      .attr("width", x(t))
      .attr("height", y(t))
      .attr("fill", "#e8f5e9")
      .attr("opacity", 0.7);
    g.append("line")
      .attr("x1", x(t))
      .attr("x2", x(t))
      .attr("y1", 0)
      .attr("y2", ih)
      .attr("stroke", "#9aa5b1")
      .attr("stroke-dasharray", "3,3");
    g.append("line")
      .attr("x1", 0)
      .attr("x2", iw)
      .attr("y1", y(t))
      .attr("y2", y(t))
      .attr("stroke", "#9aa5b1")
      .attr("stroke-dasharray", "3,3");
  }
  g.append("line")
    .attr("x1", x(0))
    .attr("y1", y(0))
    .attr("x2", x(100))
    .attr("y2", y(100))
    .attr("stroke", "#607d8b")
    .attr("stroke-width", 1);
  g.append("g").attr("transform", `translate(0,${ih})`).call(d3.axisBottom(x).ticks(5)).attr("font-size", 10);
  g.append("g").call(d3.axisLeft(y).ticks(5)).attr("font-size", 10);
  svg
    .append("text")
    .attr("x", m.l + iw / 2)
    .attr("y", H - 4)
    .attr("text-anchor", "middle")
    .attr("font-size", 11)
    .attr("fill", "#52606d")
    .text("TASS (default)");
  svg
    .append("text")
    .attr("transform", `translate(11,${m.t + ih / 2}) rotate(-90)`)
    .attr("text-anchor", "middle")
    .attr("font-size", 11)
    .attr("fill", "#52606d")
    .text("Adaptive TASS");
  const colr = (r) =>
    r.call && r.call.changed
      ? r.call.a
        ? "#2e7d32"
        : "#c62828"
      : Math.abs(r.delta) >= _AD_EPS
      ? "#7b6fa8"
      : "#b0bec5";
  g.selectAll("circle")
    .data(rows)
    .enter()
    .append("circle")
    .attr("cx", (r) => x(Math.max(0, Math.min(100, r.std))))
    .attr("cy", (r) => y(Math.max(0, Math.min(100, r.ad))))
    .attr("r", 3.5)
    .attr("fill", colr)
    .attr("fill-opacity", 0.8)
    .attr("stroke", "#fff")
    .attr("stroke-width", 0.6)
    .on("mouseover", (ev, r) => {
      if (typeof showTip !== "function") return;
      showTip(
        `<b>${_esc(r.organism)}</b><br>Sample: ${_esc(r.sample)}<br>TASS (default) ${r.std.toFixed(
          1,
        )} → adaptive ${r.ad.toFixed(1)} (${_adSigned(r.delta)})` +
          (r.call ? `<br>Call: ${_adCallTxt(r.call)} (cutoff ${r.call.thr})` : "") +
          (r.reasons.length
            ? `<br><span style="opacity:.85">${r.reasons.slice(0, 3).map(_esc).join("<br>")}</span>`
            : ""),
        ev,
      );
    })
    .on("mouseout", () => typeof hideTip === "function" && hideTip());
}

function _adRenderModal() {
  const ov = document.getElementById("adt-modal");
  if (!ov) return;
  const on = ttAdaptScoreActive();
  ov.querySelector(".adt-sub").innerHTML = _adModalSub();
  ov.querySelector(".adt-samples").innerHTML = _adSampleTable();
  const rows = _adFiltered();
  _adDrawPlot(ov.querySelector(".adt-plot"), rows);
  const up = rows.filter((x) => x.delta >= _AD_EPS).length,
    down = rows.filter((x) => x.delta <= -_AD_EPS).length;
  const gained = rows.filter((x) => x.call && x.call.changed && x.call.a).length,
    lost = rows.filter((x) => x.call && x.call.changed && !x.call.a).length;
  ov.querySelector(".adt-plotnote").innerHTML = rows.length
    ? `<b>${rows.length}</b> organism${rows.length === 1 ? "" : "s"} shown.<br>` +
      `Points on the diagonal score the same both ways; above it adaptive is higher ` +
      `(<span style="color:#2e7d32">${up}</span>), below it lower (<span style="color:#c62828">${down}</span>; ` +
      `changes of ${_AD_EPS} point or more).<br>` +
      `Dashed lines are the current cutoff. <span style="color:#2e7d32">Green</span> points pass only with ` +
      `adaptive TASS (${gained}); <span style="color:#c62828">red</span> points pass with TASS (default) but ` +
      `fail with adaptive (${lost}).<br><span style="color:#6b7480">${
        on
          ? "The report currently uses adaptive TASS for every call (same cutoffs)."
          : "The report uses TASS (default) for every call."
      }</span>`
    : `<span style="color:#888">No organisms with an adaptive score at this level.</span>`;

  const cols = _adCols();
  const head = ov.querySelector(".adt-orgs thead");
  head.innerHTML =
    "<tr>" +
    cols
      .map(
        (c) =>
          `<th data-k="${c.key}" class="${c.num ? "num" : ""}"${c.tip ? ` title="${_esc(c.tip)}"` : ""}>${c.label}${
            _adSort.key === c.key ? (_adSort.dir > 0 ? " ▴" : " ▾") : ""
          }</th>`,
      )
      .join("") +
    "</tr>";
  head.querySelectorAll("th").forEach((th) =>
    th.addEventListener("click", () => {
      const k = th.dataset.k;
      _adSort = {
        key: k,
        dir: _adSort.key === k ? -_adSort.dir : k === "sample" || k === "organism" || k === "reasons" ? 1 : -1,
      };
      _adRenderModal();
    }),
  );
  const f1 = (v) => (v === null || v === undefined ? "—" : (+v).toFixed(1));
  ov.querySelector(".adt-orgs tbody").innerHTML = rows.length
    ? rows
        .map((x) => {
          const cc = x.call && x.call.changed ? (x.call.a ? "mc-gain" : "mc-loss") : "";
          const dc = x.delta >= _AD_EPS ? "#2e7d32" : x.delta <= -_AD_EPS ? "#c62828" : "#666";
          const b = (v, inUse) => (inUse ? `<b>${f1(v)}</b>` : f1(v));
          return (
            `<tr><td>${_esc(x.sample)}</td><td><i>${_esc(x.organism)}</i></td>` +
            `<td class="num">${b(x.std, !on)}</td><td class="num">${b(x.ad, on)}</td>` +
            `<td class="num" style="color:${dc}">${_adSigned(x.delta)}</td>` +
            `<td class="${cc}">${_adCallTxt(x.call)}</td>` +
            `<td class="num">${x.penalty === null ? "—" : x.penalty.toFixed(3)}</td>` +
            `<td class="num">${x.cutShift === null ? "—" : _adSigned(x.cutShift)}</td>` +
            `<td class="adt-why">${x.reasons.map(_esc).join("<br>") || "—"}</td></tr>`
          );
        })
        .join("")
    : `<tr><td colspan="${cols.length}" style="color:#888;padding:10px">No rows.</td></tr>`;
}

function _adCsv() {
  const rows = _adFiltered();
  const q = (v) => {
    const s = v === null || v === undefined ? "" : String(v);
    return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  };
  const lines = [
    [
      "sample",
      "organism",
      "level",
      "tass_default",
      "tass_adaptive",
      "change",
      "call_default",
      "call_adaptive",
      "cutoff",
      "ambiguity_penalty",
      "cutoff_shift",
      "reasons",
    ].join(","),
  ];
  rows.forEach((x) =>
    lines.push(
      [
        x.sample,
        x.organism,
        x.level,
        x.std,
        x.ad,
        x.delta.toFixed(4),
        x.call ? (x.call.s ? "pass" : "fail") : "",
        x.call ? (x.call.a ? "pass" : "fail") : "",
        x.call ? x.call.thr : "",
        x.penalty,
        x.cutShift,
        x.reasons.join("; "),
      ]
        .map(q)
        .join(","),
    ),
  );
  const blob = new Blob([lines.join("\n") + "\n"], { type: "text/csv" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = "tass_default_vs_adaptive.csv";
  document.body.appendChild(a);
  a.click();
  setTimeout(() => {
    URL.revokeObjectURL(a.href);
    a.remove();
  }, 0);
}

/* ── wiring ────────────────────────────────────────────────────────────── */
(function () {
  const _redraw = () => {
    if (typeof redraw === "function") redraw();
  };
  // TASS used across the report
  const sel = document.getElementById("adt-score");
  if (sel)
    sel.addEventListener("change", () => {
      _ttAdaptScoreVal = sel.value === "adaptive" ? "adaptive" : "default";
      ttAdaptEnsure();
      ttAdaptRefreshPanel();
      _redraw();
    });
  // Heatmap change view (adaptive − default)
  const hd = document.getElementById("adt-hm-delta");
  if (hd)
    hd.addEventListener("change", () => {
      _ttAdaptHmDelta = hd.checked;
      if (_ttAdaptHmDelta) {
        // the view applies to the TASS value; switch the heatmap to it
        const hv = document.getElementById("hm-value-col");
        if (hv && hv.value !== "TASS Score") hv.value = "TASS Score";
        // one alternative heatmap TASS at a time: reset the masking view
        if (typeof _ttMaskCmpModeVal !== "undefined" && _ttMaskCmpModeVal !== "masked") {
          _ttMaskCmpModeVal = "masked";
          if (typeof ttMaskCmpRefreshPanel === "function") ttMaskCmpRefreshPanel();
        }
      }
      ttAdaptRefreshPanel();
      _redraw();
    });
  // choosing a masking heatmap view turns the adaptive change view off
  // (53_mask_compare.js); the report-wide score is left alone
  const msel = document.getElementById("mcmp-mode");
  if (msel)
    msel.addEventListener("change", () => {
      if (msel.value !== "masked" && _ttAdaptHmDelta) {
        _ttAdaptHmDelta = false;
        ttAdaptRefreshPanel();
        _redraw();
      }
    });
  const cb = document.getElementById("adt-show-tables");
  if (cb)
    cb.addEventListener("change", () => {
      _ttAdaptShowTables = cb.checked;
      _redraw();
    });
  const btn = document.getElementById("adt-open-btn");
  if (btn) btn.addEventListener("click", _adOpenModal);
  const pill = document.getElementById("banner-adt");
  if (pill) {
    pill.addEventListener("click", _adOpenModal);
    pill.addEventListener("keydown", (e) => {
      if (e.key === "Enter" || e.key === " ") {
        e.preventDefault();
        _adOpenModal();
      }
    });
  }
  // SAMPLE_META / DATA arrive during init (or a later file drop): refresh then
  window.addEventListener("load", () => setTimeout(ttAdaptRefreshPanel, 0));
})();
