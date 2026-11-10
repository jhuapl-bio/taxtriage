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
       -     (make_report.py --controls / params.report_controls), "show".
═══════════════════════════════════════════════════════════════════════════ */

let TT_CONTROLS_VIEW = "show"; // "show" | "hide"
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

function ttControlsSyncUi() {
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
  TT_CONTROLS_VIEW = (view || boot) === "hide" ? "hide" : "show";
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
