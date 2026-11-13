/* ═══════════════════════════════════════════════════════════════════════════
       -  §  CONTROL SAMPLES  (show / hide negative + positive controls)
       -     A sample is a control when its metadata carries a control_type
       -     (match_paths.py writes "negative" / "positive" for samplesheet
       -     controls). The "Show controls" switch above the sample list hides
       -     or shows all of them at once.
       -
       -     Hiding drives the EXISTING sampleHidden map, exactly like the
       -     sample QC "hide" action (41_sample_flags.js) and the sidebar
       -     search auto-hide (03_sample_sidebar.js), so every tab, chart,
       -     table and export honours it with no further changes. Same
       -     contract as those two: only samples WE hid (_ctlAutoHidden) are
       -     ever given back, and none of the three gives back a sample that
       -     another one is still holding down.
       -
       -     The rows stay in DATA while hidden, so other features that look
       -     across samples (e.g. index hopping, 54_index_hopping.js) still see
       -     the controls. Default view: BOOT.controls_view
       -     (make_report.py --controls / params.report_controls), "hide".
       -     A report whose samples are ALL controls opens with them shown
       -     (unless a saved session asks otherwise), so it is never empty.
═══════════════════════════════════════════════════════════════════════════ */

let TT_CONTROLS_VIEW = "hide"; // "show" | "hide" (controls start hidden)
const _ctlAutoHidden = new Set(); // ids WE hid

function ttControlType(id) {
  const m = (typeof SAMPLE_META !== "undefined" && SAMPLE_META && SAMPLE_META[id]) || {};
  const ct = String(m.control_type == null ? "" : m.control_type)
    .trim()
    .toLowerCase();
  return ct && ct !== "none" && ct !== "null" && ct !== "false" ? ct : "";
}

function _ttCtlAllSamples() {
  const ids = new Set();
  if (typeof _allSampleIds === "function") {
    try {
      _allSampleIds().forEach((s) => s && ids.add(s));
    } catch (e) {}
  }
  if (typeof SAMPLE_META !== "undefined" && SAMPLE_META) Object.keys(SAMPLE_META).forEach((s) => ids.add(s));
  if (typeof DATA !== "undefined") DATA.forEach((r) => r["Specimen ID"] && ids.add(r["Specimen ID"]));
  return Array.from(ids);
}

/** True when `id` is a control the switch is currently hiding: the sample
 *  list (03_sample_sidebar.js) leaves it out instead of showing it greyed. */
function ttControlListHidden(id) {
  return TT_CONTROLS_VIEW === "hide" && !!ttControlType(id) && !!sampleHidden[id];
}

function ttControlSamples() {
  return _ttCtlAllSamples().filter((s) => ttControlType(s));
}

/** Apply the current view to sampleHidden. Returns true when anything changed. */
function ttControlsApply() {
  let changed = false;
  const hide = TT_CONTROLS_VIEW === "hide";
  _ttCtlAllSamples().forEach((id) => {
    const want = hide && !!ttControlType(id);
    if (want) {
      if (!sampleHidden[id]) {
        sampleHidden[id] = true;
        _ctlAutoHidden.add(id);
        changed = true;
      }
    } else if (_ctlAutoHidden.has(id)) {
      const heldBySearch = typeof _filterAutoHidden !== "undefined" && _filterAutoHidden.has(id);
      const heldByQc = typeof _flagAutoHidden !== "undefined" && _flagAutoHidden.has(id);
      if (!heldBySearch && !heldByQc) sampleHidden[id] = false;
      _ctlAutoHidden.delete(id);
      changed = true;
    }
  });
  if (changed && typeof _invalidateFilterCache === "function") _invalidateFilterCache();
  return changed;
}

/* ── Positive-control check per detection (Summary detections table) ──────
   For every detection, was the same organism also detected in the run's
   positive control(s)? Read from DATA, so it works while controls are hidden.
   Match: same Level + taxid; else the same species name (another strain or
   level). Result cached per DATA size + control list. */
let _ttPcCache = { key: null, byTax: null, bySp: null, pcs: [] };
function _ttPcIndex() {
  const pcs = ttControlSamples().filter((s) => ttControlType(s) === "positive");
  const key = (typeof DATA !== "undefined" ? DATA.length : 0) + "|" + pcs.join(",");
  if (_ttPcCache.key === key) return _ttPcCache;
  const byTax = new Map(),
    bySp = new Map();
  const pcSet = new Set(pcs);
  const truthy = typeof isTruthy === "function" ? isTruthy : (v) => !!v;
  const put = (map, k, hit) => {
    if (!k) return;
    const list = map.get(k) || [];
    const i = list.findIndex((h) => h.sample === hit.sample);
    if (i < 0) list.push(hit);
    else if (hit.tass > list[i].tass) list[i] = hit;
    map.set(k, list);
  };
  if (pcs.length && typeof DATA !== "undefined") {
    DATA.forEach((r) => {
      const s = r && r["Specimen ID"];
      if (!pcSet.has(s)) return;
      const hit = {
        sample: s,
        name: r["Detected Organism"] || "",
        tass: Number(r["TASS Score"]) || 0,
        reads: Number(r["# Reads Aligned"]) || 0,
        pass: truthy(r["Passes Threshold"]),
      };
      put(byTax, (r.Level || "Strain") + "|" + String(r["Taxonomic ID #"] || ""), hit);
      put(
        bySp,
        String(r["Species Name"] || "")
          .trim()
          .toLowerCase(),
        hit,
      );
    });
  }
  _ttPcCache = { key, byTax, bySp, pcs };
  return _ttPcCache;
}
function _ttPcEsc(s) {
  return String(s == null ? "" : s)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}
/** { state: "pass" | "below" | "absent", exact, hits, pcs } or null when the
 *  run has no positive control (or the row IS a positive control). */
function ttPosCtlStatus(r) {
  if (!r) return null;
  const idx = _ttPcIndex();
  if (!idx.pcs.length || idx.pcs.includes(r["Specimen ID"])) return null;
  let hits = idx.byTax.get((r.Level || "Strain") + "|" + String(r["Taxonomic ID #"] || "")) || [];
  let exact = true;
  if (!hits.length) {
    const sp = String(r["Species Name"] || "")
      .trim()
      .toLowerCase();
    hits = (sp && idx.bySp.get(sp)) || [];
    exact = false;
  }
  if (!hits.length) return { state: "absent", exact: false, hits, pcs: idx.pcs };
  return { state: hits.some((h) => h.pass) ? "pass" : "below", exact, hits, pcs: idx.pcs };
}
/** Small "PC" tag after the organism name: green = in the positive control
 *  (passing), amber = there but below the cutoff, grey = not in it. */
/* Hover cards for the PC tags. Native `title` tooltips do not show in every
   viewer (and the summary row swallows them), so the tags carry a key into
   _ttPcTips and one delegated listener shows the report's own card (showTip),
   like the index-hop icons (54_index_hopping.js). */
const _ttPcTips = new Map(); // key → card HTML (bounded: one per sample/level/taxid)
const _TT_PC_LABEL = {
  pass: ["Also in the positive control", "#00897b"],
  below: ["In the positive control, not passing there", "#b26a00"],
  absent: ["Not in the positive control", "#90a4ae"],
};
function _ttPcTipHTML(r, st) {
  const e = _ttPcEsc;
  const [lab, col] = _TT_PC_LABEL[st.state];
  let h =
    `<div style="max-width:420px"><b>${e(r["Detected Organism"] || "")}</b> · <span style="color:#ccc">${e(
      r["Specimen ID"],
    )}</span><br>` +
    `<span style="display:inline-block;margin:3px 0;padding:0 6px;border-radius:3px;background:${col};color:#fff;font-weight:700;font-size:0.85em">${lab}</span>`;
  if (st.state === "absent") {
    h += `<br><span style="color:#ccc">Positive control${st.pcs.length > 1 ? "s" : ""} checked: ${st.pcs
      .map(e)
      .join(", ")}. The organism was not detected there at any TASS.</span>`;
  } else {
    if (!st.exact)
      h += `<br><span style="color:#ffd580;font-size:0.85em">Matched by species (the control carries a different strain / level).</span>`;
    h +=
      `<table style="border-collapse:collapse;font-size:0.85em;margin:5px 0">` +
      `<tr style="color:#ffd580"><td style="padding-right:10px">Control</td>${
        st.exact ? "" : '<td style="padding-right:10px">As</td>'
      }<td style="padding-right:10px">TASS</td><td style="padding-right:10px">Reads</td><td>Passes</td></tr>` +
      st.hits
        .map(
          (x) =>
            `<tr><td style="padding-right:10px">${e(x.sample)}</td>${
              st.exact ? "" : `<td style="padding-right:10px"><i>${e(x.name)}</i></td>`
            }<td style="padding-right:10px">${x.tass.toFixed(1)}</td><td style="padding-right:10px">${Math.round(
              x.reads,
            ).toLocaleString()}</td><td>${
              x.pass ? '<span style="color:#80cbc4">yes</span>' : '<span style="color:#ffcc80">no</span>'
            }</td></tr>`,
        )
        .join("") +
      `</table><span style="color:#aaa;font-size:0.8em">A detection that matches what the positive control carries can come from carry-over or cross-contamination from it; compare reads and TASS.</span>`;
  }
  return h + `</div>`;
}
/** Small "PC" tag after the organism name: green = in the positive control
 *  (passing), amber = there but below the cutoff, grey = not in it. */
function ttPosCtlBadgeHTML(r) {
  const st = ttPosCtlStatus(r);
  if (!st) return "";
  const key = [r["Specimen ID"], r.Level || "Strain", r["Taxonomic ID #"], r["Detected Organism"]].join("|");
  _ttPcTips.set(key, _ttPcTipHTML(r, st));
  const sym = st.state === "absent" ? "–" : "✓";
  return `<span class="pc-mark pc-${st.state}" data-pc-tip="${_ttPcEsc(key)}">PC${sym}</span>`;
}
const _TT_PC_LEGEND_TIP =
  `<div style="max-width:380px"><b>Positive-control check</b><br>` +
  `<span style="color:#ccc">For each detection: was the same organism also found in the run's positive control? ` +
  `Matched on the same taxid and level, or else the same species. Works while controls are hidden.</span>` +
  `<table style="border-collapse:collapse;font-size:0.85em;margin:5px 0">` +
  `<tr><td style="padding-right:8px"><span class="pc-mark pc-pass">PC✓</span></td><td>also in the positive control, and passing there</td></tr>` +
  `<tr><td style="padding-right:8px"><span class="pc-mark pc-below">PC✓</span></td><td>detected in the positive control, but not passing there</td></tr>` +
  `<tr><td style="padding-right:8px"><span class="pc-mark pc-absent">PC–</span></td><td>not detected in the positive control</td></tr>` +
  `</table><span style="color:#aaa;font-size:0.8em">Hover a tag in the table for the control's TASS and reads.</span></div>`;
(function _ttPcWireTips() {
  if (typeof document === "undefined") return;
  const target = (ev) => (ev.target && ev.target.closest ? ev.target.closest("[data-pc-tip]") : null);
  const html = (el) => {
    const k = el.getAttribute("data-pc-tip");
    return k === "legend" ? _TT_PC_LEGEND_TIP : _ttPcTips.get(k) || "";
  };
  document.addEventListener("mouseover", (ev) => {
    const el = target(ev);
    if (!el || typeof showTip !== "function") return;
    ev.stopPropagation();
    const h = html(el);
    if (h) showTip(h, ev);
  });
  document.addEventListener("mousemove", (ev) => {
    if (target(ev) && typeof moveTip === "function") moveTip(ev);
  });
  document.addEventListener("mouseout", (ev) => {
    const el = target(ev);
    if (el && !el.contains(ev.relatedTarget) && typeof hideTip === "function") hideTip();
  });
})();

function ttControlsSyncUi() {
  const pcLeg = document.getElementById("sum-pc-legend");
  if (pcLeg) pcLeg.hidden = !_ttPcIndex().pcs.length;
  const row = document.getElementById("controls-toggle-row");
  const cb = document.getElementById("show-controls");
  const info = document.getElementById("controls-toggle-info");
  if (!row) return;
  const ctl = ttControlSamples();
  row.hidden = ctl.length === 0;
  row.style.display = ctl.length ? "flex" : "none"; // inline style, so [hidden] alone would not hide it
  if (cb) cb.checked = TT_CONTROLS_VIEW !== "hide";
  if (info) {
    const by = {};
    ctl.forEach((s) => {
      const t = ttControlType(s);
      by[t] = (by[t] || 0) + 1;
    });
    info.textContent = ctl.length
      ? `(${Object.keys(by)
          .sort()
          .map((t) => `${by[t]} ${t}`)
          .join(", ")})`
      : "";
    row.title = ctl.length ? `Control samples: ${ctl.join(", ")}` : "";
  }
}

/** Switch the view ("show" | "hide") and repaint. */
function ttControlsSetView(view) {
  TT_CONTROLS_VIEW = view === "hide" ? "hide" : "show";
  ttControlsApply();
  if (typeof buildSampleList === "function") buildSampleList();
  if (typeof redraw === "function") redraw();
  ttControlsSyncUi();
}

/** Called once from the report init, before the first sample list / redraw. */
function ttControlsInit(view) {
  const boot = typeof BOOT !== "undefined" && BOOT ? BOOT.controls_view : null;
  const want = view || boot || "hide";
  TT_CONTROLS_VIEW = want === "show" ? "show" : "hide";
  // Default-hidden must never leave an empty report: when every sample is a
  // control (a controls-only run), open with them shown. An explicit saved
  // session view (`view`) is honoured as-is.
  if (!view && TT_CONTROLS_VIEW === "hide") {
    const all = _ttCtlAllSamples();
    if (all.length && all.every((s) => ttControlType(s))) TT_CONTROLS_VIEW = "show";
  }
  ttControlsApply();
  ttControlsSyncUi();
}

/** After detections are added or removed (uploads, clear): newly loaded
 *  controls follow the current view. The caller redraws. */
function ttControlsOnDataChanged() {
  ttControlsApply();
  ttControlsSyncUi();
}

/** A user's explicit eye / Show all / remove action takes ownership of the
 *  sample's visibility away from the switch. */
function ttControlsRelease(id) {
  _ctlAutoHidden.delete(id);
}

(function () {
  if (typeof document === "undefined") return;
  const wire = () => {
    const cb = document.getElementById("show-controls");
    if (cb && !cb._ttWired) {
      cb._ttWired = true;
      cb.addEventListener("change", () => ttControlsSetView(cb.checked ? "show" : "hide"));
    }
  };
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", wire);
  else wire();
})();
