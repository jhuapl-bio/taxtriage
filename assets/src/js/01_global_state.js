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

/** Return arr sorted by _sampleOrder; unknowns go to the end. */
function _orderedSamples(arr) {
  if (!_sampleOrder.length) return arr;
  const idx = Object.fromEntries(_sampleOrder.map((id, i) => [id, i]));
  return [...arr].sort((a, b) => {
    const ia = idx[a] !== undefined ? idx[a] : 9999;
    const ib = idx[b] !== undefined ? idx[b] : 9999;
    return ia - ib;
  });
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
