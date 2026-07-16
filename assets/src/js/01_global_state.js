/* ═══════════════════════════════════════════════════════════════════════════
       -  §  GLOBAL STATE
       -     Shared mutable state used across every tab: sample colors, the
       -     hidden / rescale flags, custom sample order, and the global DATA
       -     array (set by the page-level template).
═══════════════════════════════════════════════════════════════════════════ */
const sampleColors = {};
const sampleHidden = {};
const sampleRescale = {}; // per-sample flag: scale legacy 0–1 TASS/Coverage ×100
const perTypeTass = {}; // per-sample-type TASS override: { "stool": 65, "blood": 80, ... }
const _loadingSampleIds = new Set(); // samples currently ingesting an attached file (drives the per-row spinner)
let _sampleOrder = []; // custom display order for samples (index → id)

/* ── Specimen grouping (merge DNA/RNA libraries of one specimen) ──────────
       A "specimen" groups one or more samples (e.g. a separate DNA and RNA
       library from the same swab) into a single biological unit. The grouping
       is supplied by the samplesheet `specimen` column, which flows through the
       pipeline into SAMPLE_META[sample].specimen, and can be overridden live in
       the report via SPECIMEN_OVERRIDE. When specimenMergeEnabled is on, the
       cross-sample views aggregate hits / TASS / coverage per specimen instead
       of per sample, and prevalence is computed over specimens. */
let specimenMergeEnabled = false;
const SPECIMEN_OVERRIDE = {}; // sample_name → specimen group (live UI edits win)

/** Resolve the specimen a sample belongs to, independent of the merge toggle.
 *  Priority: live override → samplesheet metadata → the sample name itself
 *  (a sample with no specimen assignment is its own specimen). */
function specimenOf(sample) {
  const s = String(sample == null ? "" : sample);
  if (!s) return s;
  if (Object.prototype.hasOwnProperty.call(SPECIMEN_OVERRIDE, s) && SPECIMEN_OVERRIDE[s]) {
    return String(SPECIMEN_OVERRIDE[s]);
  }
  const meta = (typeof SAMPLE_META !== "undefined" && SAMPLE_META[s]) || {};
  const sp = meta.specimen != null ? meta.specimen : meta.specimen_id != null ? meta.specimen_id : meta.specimen_group;
  return sp != null && String(sp).trim() ? String(sp).trim() : s;
}

/** The key cross-sample views should group rows by: the specimen when merging
 *  is enabled, otherwise the raw sample (Specimen ID) name. */
function specimenKey(sample) {
  return specimenMergeEnabled ? specimenOf(sample) : String(sample == null ? "" : sample);
}

/** Give every multi-sample specimen a stable color so legends, the heatmap,
 *  the coverage plots and the right-panel all render the merged unit in a
 *  consistent hue. A specimen inherits its first member's color (visually
 *  linking them); brand-new specimen names with no colored member fall back to
 *  the shared PALETTE. Idempotent — only fills gaps, never recolors. */
function _ensureSpecimenColors() {
  if (typeof specimenGroups !== "function") return;
  let i = Object.keys(sampleColors).length;
  specimenGroups().forEach((members, g) => {
    if (members.length > 1 && !sampleColors[g]) {
      sampleColors[g] = sampleColors[members[0]] || PALETTE[i++ % PALETTE.length];
    }
  });
}

/** True when merging is active AND this sample has been folded into a
 *  specimen that carries a different name (i.e. it no longer appears in the
 *  merged view under its own id). Per-sample section / indicator builders use
 *  this to avoid rendering an empty "no detections" row for a member sample
 *  whose rows were relabelled onto its specimen. */
function _isMergedAway(sample) {
  if (!specimenMergeEnabled) return false;
  const s = String(sample == null ? "" : sample);
  return !!s && specimenOf(s) !== s;
}

/** Whether any sample is actually assigned to a multi-sample specimen (via
 *  metadata or a live override). Lets the UI hide the merge control when the
 *  run has no specimen grouping to act on. */
function hasSpecimenGrouping() {
  const groups = specimenGroups();
  for (const members of groups.values()) if (members.length > 1) return true;
  return false;
}

/** Map every known specimen → sorted list of its member samples, using
 *  SAMPLE_META (all samples) plus any live overrides and anything seen in DATA. */
function specimenGroups() {
  const groups = new Map();
  const add = (sample) => {
    const s = String(sample == null ? "" : sample);
    if (!s) return;
    const g = specimenOf(s);
    if (!groups.has(g)) groups.set(g, new Set());
    groups.get(g).add(s);
  };
  Object.keys((typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {}).forEach(add);
  if (typeof DATA !== "undefined") DATA.forEach((r) => add(r["Specimen ID"]));
  const out = new Map();
  Array.from(groups.keys())
    .sort()
    .forEach((g) => out.set(g, Array.from(groups.get(g)).sort()));
  return out;
}

/** Small, consistent visual marker used beside sample/specimen names in the
 *  Summary, full Table and Metadata tables. Pass either a collapsed detection
 *  row (__mergedFrom) or a raw sample id. Returns an empty string unless the
 *  item is currently part of an active multi-sample specimen. */
function _mergedSampleBadgeHTML(item) {
  if (!(typeof specimenMergeEnabled !== "undefined" && specimenMergeEnabled)) return "";
  const esc = (v) =>
    String(v == null ? "" : v)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  let specimen = "",
    members = [];
  if (item && typeof item === "object") {
    specimen = String(item["Specimen ID"] || "");
    members = Array.isArray(item.__specimenMembers)
      ? item.__specimenMembers.slice()
      : Array.isArray(item.__mergedFrom)
      ? item.__mergedFrom.slice()
      : [];
    if (members.length < 2 && typeof specimenGroups === "function") {
      members = (specimenGroups().get(specimen) || members).slice();
    }
  } else {
    const sample = String(item == null ? "" : item);
    specimen = typeof specimenOf === "function" ? specimenOf(sample) : sample;
    const groups = typeof specimenGroups === "function" ? specimenGroups() : new Map();
    members = (groups.get(specimen) || []).slice();
  }
  if (members.length < 2) return "";
  const title = `Merged specimen: ${specimen} (${members.length} samples: ${members.join(", ")})`;
  return (
    `<span class="merged-sample-badge" title="${esc(title)}" ` +
    `style="display:inline-flex;align-items:center;gap:3px;margin-left:5px;padding:0 5px;` +
    `border:1px solid #90caf9;border-radius:8px;background:#e3f2fd;color:#0d47a1;` +
    `font-size:9px;font-weight:700;line-height:1.55;vertical-align:middle;white-space:nowrap">` +
    `<i class="fas fa-layer-group" aria-hidden="true"></i> merged ×${members.length}</span>`
  );
}

/** Fingerprint of everything that changes specimen aggregation, for the
 *  filteredData() cache key. Changes when the merge toggle flips or any live
 *  SPECIMEN_OVERRIDE assignment is edited. */
function _hashSpecimenMerge() {
  const ov = Object.keys(SPECIMEN_OVERRIDE)
    .sort()
    .map((k) => k + "=" + SPECIMEN_OVERRIDE[k])
    .join(",");
  return (specimenMergeEnabled ? "1" : "0") + "|" + ov;
}

/** User-chosen resolutions for metadata that conflicts across the samples of a
 *  merged specimen (e.g. differing lat/long). Keyed specimen → { field: value }.
 *  Populated by the merge modal's conflict chooser; read best-effort by views
 *  that surface per-specimen metadata. */
const SPECIMEN_META_RESOLVED = {};

/** Prevalence denominator: the number of distinct specimens (respecting the
 *  merge toggle) that have at least one positive hit (Passes Threshold)
 *  anywhere in the run. Deliberately decoupled from the live TASS slider so
 *  prevalence is measured against everything that tested positive, not just
 *  what clears the current threshold. */
function positiveHitSpecimenCount() {
  const seen = new Set();
  if (typeof DATA !== "undefined") {
    DATA.forEach((r) => {
      if (typeof isTruthy === "function" ? isTruthy(r["Passes Threshold"]) : r["Passes Threshold"]) {
        const k = specimenKey(r["Specimen ID"]);
        if (k) seen.add(k);
      }
    });
  }
  return seen.size;
}

/** Return arr sorted by _sampleOrder; unknowns go to the end. `arr` may
 *  contain raw sample ids or specimen/group names (callers pass whichever
 *  filteredData() rows carry under "Specimen ID" — a group name once
 *  specimen merge is on), so the lookup below must resolve both. */
function _orderedSamples(arr) {
  if (!_sampleOrder.length) return arr;
  const idx = _sampleOrGroupIndexMap();
  return [...arr].sort((a, b) => {
    const ia = idx[a] !== undefined ? idx[a] : 9999;
    const ib = idx[b] !== undefined ? idx[b] : 9999;
    return ia - ib;
  });
}

/** Build a lookup from "row group key" → right-panel display position, for
 *  tables/plots that sort by the row's "Specimen ID" field. That field holds
 *  the raw sample id normally, but holds the *specimen/group name* once
 *  specimen merge is enabled (see specimenKey()) — a string that never
 *  appears in _sampleOrder itself (which only lists individual sample ids).
 *  A plain `_sampleOrder.map((id,i)=>[id,i])` lookup therefore misses every
 *  merged group and silently sorts it to the end, which is why grouped views
 *  used to ignore the sidebar drag order. Here, each specimen/group name is
 *  additionally keyed to the earliest position any of its member samples
 *  occupies in _sampleOrder, so dragging a specimen's row block in the
 *  sidebar (or a whole merged-specimen header) is reflected everywhere rows
 *  are grouped by "Specimen ID". */
function _sampleOrGroupIndexMap() {
  const idx = {};
  if (!_sampleOrder.length) return idx;
  _sampleOrder.forEach((id, i) => {
    if (idx[id] === undefined) idx[id] = i;
    const g = typeof specimenOf === "function" ? specimenOf(id) : id;
    if (g && idx[g] === undefined) idx[g] = i;
  });
  return idx;
}
let sortCol = null;
let sortAsc = true;
let visibleCols = [];
let activeTab = "summary";

/* ── Follow-up / watchlist ─────────────────────────────────────────────
         A cross-tab set of organisms the user has starred for follow-up.
         Keyed by Taxonomic ID # (falls back to organism name) so a single
         star follows an organism across every sample and every tab.
         watchFilterMode controls how filteredData() treats the list:
           "all"  — no filtering (default)
           "only" — show only starred organisms everywhere
           "hide" — hide starred organisms everywhere (cull what you don't care about) */
const watchlist = new Set();
let watchFilterMode = "all";

function _watchKey(r) {
  const tid = String((r && r["Taxonomic ID #"]) || "").trim();
  if (tid) return "tax:" + tid;
  return (
    "name:" +
    String((r && r["Detected Organism"]) || "")
      .trim()
      .toLowerCase()
  );
}
function _isWatched(r) {
  return watchlist.has(_watchKey(r));
}
function _watchLabelForKey(key) {
  // Best-effort display name for a stored key (used by the chip list).
  if (key.startsWith("name:")) return key.slice(5);
  for (let i = 0; i < DATA.length; i++) {
    if (_watchKey(DATA[i]) === key) return DATA[i]["Detected Organism"] || key;
  }
  return key.replace(/^tax:/, "Taxon ");
}
// Star markup for organism name cells. interactive=true makes it a clickable
// toggle (.watch-star with data-watch-key). inline=true renders it in the normal
// flow right after the organism text instead of absolutely on the cell's right edge.
function _watchStarHTML(r, interactive, rightOffset, inline) {
  const key = _watchKey(r);
  const on = watchlist.has(key);
  const cls = "watch-star" + (on ? " on" : "");
  const icon = on ? "fas fa-star" : "far fa-star";
  const title = on ? "On the follow-up list — click to remove" : "Add to follow-up list";
  const ro = rightOffset != null ? rightOffset : 48;
  const style = inline
    ? `margin-left:5px;font-size:11px;cursor:pointer;vertical-align:middle;`
    : `position:absolute;right:${ro}px;top:50%;transform:translateY(-50%);font-size:11px;cursor:pointer;z-index:1;`;
  if (interactive) {
    return `<i class="${cls}" data-watch-key="${encodeURIComponent(
      key,
    )}" role="button" tabindex="0" title="${title}" style="${style}"><i class="${icon}"></i></i>`;
  }
  // Non-interactive marker: only render when starred.
  return on
    ? `<i class="watch-star on" title="On the follow-up list" style="${style}"><i class="${icon}"></i></i>`
    : "";
}
function _setWatch(key, on) {
  if (on) watchlist.add(key);
  else watchlist.delete(key);
  _updateWatchPanel();
  if (typeof redraw === "function") redraw();
}
function _toggleWatchKey(key) {
  _setWatch(key, !watchlist.has(key));
}
function _clearWatchlist() {
  watchlist.clear();
  _updateWatchPanel();
  if (typeof redraw === "function") redraw();
}
function _setWatchFilterMode(mode) {
  watchFilterMode = mode;
  document.querySelectorAll(".watch-filter-btn").forEach((b) => {
    b.classList.toggle("active", b.dataset.watchMode === mode);
  });
  if (typeof redraw === "function") redraw();
}
// Rebuild the follow-up panel: count badge + removable chips. Hidden when empty.
function _updateWatchPanel() {
  const panel = document.getElementById("watch-panel");
  if (!panel) return;
  const n = watchlist.size;
  const badge = document.getElementById("watch-count-badge");
  if (badge) badge.textContent = String(n);
  panel.style.display = n > 0 ? "block" : "none";
  const chips = document.getElementById("watch-chips");
  if (chips) {
    if (!n) {
      chips.innerHTML = "";
    } else {
      chips.innerHTML = Array.from(watchlist)
        .map((key) => {
          const name = _watchLabelForKey(key);
          const safe = encodeURIComponent(key);
          const disp = String(name).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
          return (
            `<span class="watch-chip"><i class="fas fa-star" style="color:#f5a623;font-size:0.8em"></i>` +
            `<span class="watch-chip-name" title="${disp}">${disp}</span>` +
            `<span class="watch-chip-x" data-watch-key="${safe}" role="button" tabindex="0" ` +
            `title="Remove from follow-up list">&times;</span></span>`
          );
        })
        .join("");
    }
  }
}
// Delegated wiring for chip removal + filter buttons + clear-all (panel is static markup).
(function _initWatchPanel() {
  const ready = () => {
    const chips = document.getElementById("watch-chips");
    if (chips && !chips._wired) {
      chips._wired = true;
      chips.addEventListener("click", (e) => {
        const x = e.target.closest(".watch-chip-x");
        if (!x) return;
        _setWatch(decodeURIComponent(x.dataset.watchKey), false);
      });
    }
    document.querySelectorAll(".watch-filter-btn").forEach((b) => {
      if (b._wired) return;
      b._wired = true;
      b.addEventListener("click", () => _setWatchFilterMode(b.dataset.watchMode));
    });
    const clr = document.getElementById("watch-clear-btn");
    if (clr && !clr._wired) {
      clr._wired = true;
      clr.addEventListener("click", _clearWatchlist);
    }
  };
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", ready);
  else ready();
})();
// Delegated handler so stars rendered inside any tab's HTML toggle the watchlist.
document.addEventListener("click", (e) => {
  const star = e.target.closest(".watch-star[data-watch-key]");
  if (!star) return;
  e.stopPropagation();
  e.preventDefault();
  _toggleWatchKey(decodeURIComponent(star.dataset.watchKey));
});

// Palette for samples
// Colorblind-safe categorical palette (Okabe-Ito + Paul Tol muted),
// ordered so the earliest-assigned sample colors stay maximally
// distinct under deuteranopia/protanopia (no adjacent red↔green).
const PALETTE = [
  "#0072B2", // blue
  "#E69F00", // orange
  "#009E73", // bluish green
  "#CC79A7", // reddish purple
  "#56B4E9", // sky blue
  "#D55E00", // vermillion
  "#117733", // dark green
  "#882255", // wine
  "#DDCC77", // sand
  "#AA4499", // purple
  "#44AA99", // teal
  "#999933", // olive
  "#332288", // indigo
  "#661100", // brown
  "#6699CC", // steel blue
];
