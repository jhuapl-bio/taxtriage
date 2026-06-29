/* ═══════════════════════════════════════════════════════════════════════════
       -  §  FILTER LISTENERS
       -     Wires every filter input on the page to redraw(). Text + numeric
       -     inputs are debounced (180 ms via _redrawDebounced) so individual
       -     keystrokes / arrow taps don't kick off a redraw on each event.
       -     Add new filter UI element IDs to the IDs list here so they
       -     trigger a redraw on `change`.
═══════════════════════════════════════════════════════════════════════════ */
// Debounce text + numeric filter inputs so each keystroke / arrow tap
// doesn't trigger a full redraw on a 14k-row dataset.
let _redrawTimer = null;
function _redrawDebounced() {
  clearTimeout(_redrawTimer);
  _redrawTimer = setTimeout(redraw, 180);
}
["filter-text", "filter-min"].forEach((id) => document.getElementById(id).addEventListener("input", _redrawDebounced));

// Min TASS: keep the range slider and the number input synchronized.
// Either control can be the source; we mirror the value to the other,
// clamp into [0,100], and trigger a debounced redraw.
(function () {
  const num = document.getElementById("filter-min");
  const range = document.getElementById("filter-min-range");
  if (!num || !range) return;
  function _clamp(v) {
    const n = parseFloat(v);
    if (isNaN(n)) return null;
    return Math.max(0, Math.min(100, n));
  }
  // Slider drag → mirror into the number box (live, debounced redraw).
  range.addEventListener("input", () => {
    const v = _clamp(range.value);
    if (v == null) return;
    num.value = v;
    _redrawDebounced();
  });
  // Number box → mirror into the slider. We don't trigger redraw here
  // because the existing filter-min input listener above already does.
  num.addEventListener("input", () => {
    const v = _clamp(num.value);
    if (v == null) return;
    range.value = v;
  });
})();
[
  "filter-ic",
  "filter-scope",
  "filter-hc",
  "filter-pass",
  "filter-below-vfamr",
  "filter-novelty-sub",
  "hm-show-vals",
  "filter-mt-dna",
  "filter-mt-rna",
  "filter-mt-both",
].forEach((id) => {
  const el = document.getElementById(id);
  if (el) el.addEventListener("change", redraw);
});

/* ═══ Search box enhancements: live status, clear button, quick chips,
             keyboard shortcut. Reuses filteredData()'s validity flag. ═══ */
(function () {
  const input = document.getElementById("filter-text");
  const clearBtn = document.getElementById("filter-clear");
  const statusEl = document.getElementById("filter-status");
  const scopeEl = document.getElementById("filter-scope");
  const icEl = document.getElementById("filter-ic");
  const chipsEl = document.getElementById("filter-chips");
  if (!input) return;

  // Lightweight, debounced count of what the current pattern matches —
  // independent of the TASS threshold so it reflects the raw search.
  function _updateSearchStatus() {
    const txt = (input.value || "").trim();
    const hasTxt = txt.length > 0;
    document.body.classList.toggle("has-search", hasTxt);
    _renderChips();
    if (!hasTxt) {
      document.body.classList.remove("search-invalid");
      statusEl.textContent = "";
      statusEl.className = "";
      return;
    }
    const ic = icEl.checked;
    const scope = scopeEl.value || "both";
    let rx = null;
    try {
      rx = new RegExp(txt, ic ? "i" : "");
    } catch (e) {
      rx = null;
    }
    if (!rx) {
      document.body.classList.add("search-invalid");
      statusEl.textContent = "incomplete pattern";
      statusEl.className = "warn";
      return;
    }
    document.body.classList.remove("search-invalid");
    const orgs = new Set();
    const samps = new Set();
    for (let i = 0; i < DATA.length; i++) {
      const r = DATA[i];
      const sid = r["Specimen ID"] || "";
      if (sampleHidden[sid]) continue;
      const inS = rx.test(sid);
      const inO = rx.test(r["Detected Organism"] || "");
      const hit = scope === "sample" ? inS : scope === "organism" ? inO : inS || inO;
      if (hit) {
        orgs.add(r["Taxonomic ID #"] || r["Detected Organism"] || "");
        if (sid) samps.add(sid);
      }
    }
    const nO = orgs.size,
      nS = samps.size;
    if (nO === 0) {
      statusEl.textContent = "no matches";
      statusEl.className = "zero";
    } else {
      statusEl.textContent = `${nO} organism${nO === 1 ? "" : "s"} · ${nS} sample${nS === 1 ? "" : "s"}`;
      statusEl.className = "";
    }
  }

  let _statusTimer = null;
  function _statusDebounced() {
    clearTimeout(_statusTimer);
    _statusTimer = setTimeout(_updateSearchStatus, 120);
  }
  input.addEventListener("input", _statusDebounced);
  scopeEl.addEventListener("change", _updateSearchStatus);
  icEl.addEventListener("change", _updateSearchStatus);

  // Options panel toggle
  const optsToggle = document.getElementById("filter-opts-toggle");
  const metaPanel = document.getElementById("filter-search-meta");
  if (optsToggle && metaPanel) {
    optsToggle.addEventListener("click", () => {
      const open = metaPanel.classList.toggle("open");
      optsToggle.classList.toggle("open", open);
      optsToggle.setAttribute("aria-expanded", String(open));
      optsToggle.querySelector("i") &&
        (optsToggle.innerHTML = `<i class="fas fa-sliders-h" style="font-size:0.85em"></i> Options ${
          open ? "▴" : "▾"
        }`);
    });
  }

  // Clear button
  clearBtn.addEventListener("click", () => {
    input.value = "";
    input.dispatchEvent(new Event("input", { bubbles: true }));
    _updateSearchStatus();
    input.focus();
  });

  // Esc clears while focused; "/" focuses search from anywhere.
  input.addEventListener("keydown", (e) => {
    if (e.key === "Escape" && input.value) {
      e.preventDefault();
      input.value = "";
      input.dispatchEvent(new Event("input", { bubbles: true }));
      _updateSearchStatus();
    }
  });
  document.addEventListener("keydown", (e) => {
    if (e.key !== "/") return;
    const t = e.target;
    const typing =
      t && (t.tagName === "INPUT" || t.tagName === "TEXTAREA" || t.tagName === "SELECT" || t.isContentEditable);
    if (typing) return;
    if (document.body.classList.contains("help-active")) return;
    e.preventDefault();
    input.focus();
    input.select();
  });

  // ── Quick-filter chips: recent searches + a one-click "Clear" reset.
  // The bar stays hidden (CSS :empty) until there is something to show.
  const _recent = [];
  function _applyChip(pattern) {
    input.value = pattern;
    input.dispatchEvent(new Event("input", { bubbles: true }));
    _updateSearchStatus();
  }
  function _renderChips() {
    if (!chipsEl) return;
    const parts = [];
    if ((input.value || "").trim()) {
      parts.push(
        `<span class="filter-chip" data-pat="" title="Clear the search"><i class="fas fa-xmark" style="font-size:0.85em;margin-right:3px;opacity:.6"></i>Clear</span>`,
      );
    }
    _recent
      .slice(0, 4)
      .forEach((p) =>
        parts.push(
          `<span class="filter-chip recent" data-pat="${p.replace(
            /"/g,
            "&quot;",
          )}" title="Recent: ${p}"><i class="fas fa-clock-rotate-left" style="font-size:0.85em;margin-right:3px;opacity:.6"></i>${
            p.length > 16 ? p.slice(0, 15) + "…" : p
          }</span>`,
        ),
      );
    chipsEl.innerHTML = parts.join("");
  }
  chipsEl?.addEventListener("click", (e) => {
    const chip = e.target.closest(".filter-chip");
    if (!chip) return;
    _applyChip(chip.getAttribute("data-pat") || "");
  });
  // Remember a committed search (Enter) as a recent chip.
  input.addEventListener("keydown", (e) => {
    if (e.key !== "Enter") return;
    const v = (input.value || "").trim();
    if (v && !_recent.includes(v)) {
      _recent.unshift(v);
      if (_recent.length > 4) _recent.pop();
      _renderChips();
    }
  });
  _renderChips();
  _updateSearchStatus();
})();
// Threshold-rollup toggle: drive ROLLUP_PASS + invalidate filter cache.
(function () {
  const el = document.getElementById("filter-rollup");
  if (!el) return;
  ROLLUP_PASS = el.checked;
  el.addEventListener("change", () => {
    ROLLUP_PASS = el.checked;
    _invalidateFilterCache();
    redraw();
  });
})();
// ── Info-tip icon hover (sidebar tooltips) ──────────────────────────────
document.querySelectorAll(".info-tip-btn[data-tip]").forEach((btn) => {
  btn.addEventListener("mouseenter", (ev) => showTip(btn.dataset.tip, ev));
  btn.addEventListener("mousemove", moveTip);
  btn.addEventListener("mouseleave", hideTip);
  // Prevent label click from toggling the associated checkbox
  btn.addEventListener("click", (ev) => ev.preventDefault());
});

// ── Rollup row hover: highlight all visible tab buttons ─────────────────
(function () {
  const row = document.getElementById("rollup-row");
  if (!row) return;
  row.addEventListener("mouseenter", () => {
    document.querySelectorAll(".tab-btn:not(.hidden)").forEach((b) => b.classList.add("rollup-highlight"));
  });
  row.addEventListener("mouseleave", () => {
    document.querySelectorAll(".tab-btn").forEach((b) => b.classList.remove("rollup-highlight"));
  });
})();

// View-level dropdown: strain/species/genus granularity → invalidate the
// filtered-data cache (its result depends on Level) and redraw everything.
function _updateHistTabDisabled() {
  const level = document.getElementById("view-level")?.value || "Strain";
  const histBtn = document.getElementById("hist-tab-btn");
  const warnIcon = document.getElementById("hist-tab-warn");
  if (!histBtn) return;
  const shouldDisable = level === "Species" || level === "Genus";
  histBtn.classList.toggle("tab-disabled", shouldDisable);
  if (warnIcon) warnIcon.style.display = shouldDisable ? "" : "none";
  // Remove any previously attached tooltip listeners before re-attaching
  if (histBtn._disabledTipOver) histBtn.removeEventListener("mouseover", histBtn._disabledTipOver);
  if (histBtn._disabledTipMove) histBtn.removeEventListener("mousemove", histBtn._disabledTipMove);
  if (histBtn._disabledTipOut) histBtn.removeEventListener("mouseout", histBtn._disabledTipOut);
  if (shouldDisable) {
    const _msg = `<b><i class="fas fa-triangle-exclamation" style="color:#f59f00"></i> Histograms unavailable</b><br><span style="font-size:0.88em;color:#ccc">Depth histograms are only available at <b>Strain</b> view level — they require individual alignment data that is not present for ${level}-level aggregates.<br>Switch <b>View level</b> back to <b>Strain</b> to enable this tab.</span>`;
    histBtn._disabledTipOver = (e) => showTip(_msg, e);
    histBtn._disabledTipMove = (e) => moveTip(e);
    histBtn._disabledTipOut = () => hideTip();
    histBtn.addEventListener("mouseover", histBtn._disabledTipOver);
    histBtn.addEventListener("mousemove", histBtn._disabledTipMove);
    histBtn.addEventListener("mouseout", histBtn._disabledTipOut);
    // If the histogram tab is currently active, switch to summary
    if (activeTab === "histogram") {
      const summaryBtn = document.querySelector('.tab-btn[data-tab="summary"]');
      if (summaryBtn) summaryBtn.click();
    }
  }
}
(function () {
  const el = document.getElementById("view-level");
  if (!el) return;
  el.addEventListener("change", () => {
    _invalidateFilterCache();
    _updateHistTabDisabled();
    redraw();
  });
})();
// ── Client-side hierarchy synthesis ──────────────────────────────────
// Reports built from the flat TSV/XLSX path — and the bundled sample data
// in heatmap_boot.js / pages.js — have no taxonomic hierarchy: every row
// is an un-leveled strain with no Species/Genus TASS. So that the View
// level dropdown works EVERYWHERE (not just JSON-hierarchy reports), we
// synthesize Species (Subkey) and Genus rollup rows here when none exist.
//
// The synthesized parent TASS is a proxy: the MAX member TASS (a species
// is detected at least as strongly as its best strain). JSON-hierarchy
// reports already carry properly-computed Species/Genus TASS rows, so this
// is skipped for them. Idempotent — guarded by a flag.
let _HIERARCHY_SYNTHESIZED = false;
function _binomial(name) {
  const parts = String(name || "")
    .trim()
    .split(/\s+/);
  return parts.slice(0, 2).join(" ") || String(name || "");
}
function _synthesizeHierarchy() {
  if (_HIERARCHY_SYNTHESIZED) return;

  // Step 1: normalize — every row without a Level gets "Strain"
  DATA.forEach((r) => {
    if (!r["Level"]) r["Level"] = "Strain";
  });

  // Step 2: collect which (specimen × subkey) and (specimen × genus) pairs
  // already have pre-built Species / Genus rows (e.g. from BOOT hierarchy JSON).
  // We only synthesize proxy rows for pairs that are NOT already covered.
  const existingSpeciesKeys = new Set();
  const existingGenusKeys = new Set();
  for (const r of DATA) {
    const spec = r["Specimen ID"] || "";
    if (r["Level"] === "Species") {
      existingSpeciesKeys.add(spec + "\0" + (r["Subkey"] || r["Taxonomic ID #"] || ""));
    }
    if (r["Level"] === "Genus") {
      existingGenusKeys.add(spec + "\0" + (r["Genus"] || ""));
    }
  }

  const num_ = (v) => {
    const n = parseFloat(v);
    return isNaN(n) ? 0 : n;
  };
  const byGenus = new Map(); // (spec \0 genus)  → strain rows needing a proxy
  const bySpecies = new Map(); // (spec \0 subkey) → strain rows needing a proxy

  for (const r of DATA) {
    if (r["Level"] !== "Strain") continue;
    const spec = r["Specimen ID"] || "";
    const sub = r["Subkey"] || r["Taxonomic ID #"] || "";
    const gen = r["Genus"] || "";
    if (sub) {
      const k = spec + "\0" + sub;
      if (!existingSpeciesKeys.has(k)) {
        (bySpecies.get(k) || bySpecies.set(k, []).get(k)).push(r);
      }
    }
    if (gen) {
      const k = spec + "\0" + gen;
      if (!existingGenusKeys.has(k)) {
        (byGenus.get(k) || byGenus.set(k, []).get(k)).push(r);
      }
    }
  }

  // genus-level proxy TASS = max strain TASS within the genus group
  const genusTassByKey = new Map();
  for (const [k, rows] of byGenus) {
    genusTassByKey.set(k, Math.max(...rows.map((r) => num_(r["TASS Score"]))));
  }

  const newRows = [];

  // Species proxy rows + back-fill Species/Genus TASS onto contributing strains
  for (const [k, rows] of bySpecies) {
    const rep = rows.reduce((a, b) => (num_(b["TASS Score"]) > num_(a["TASS Score"]) ? b : a));
    const spTass = num_(rep["TASS Score"]);
    const gKey = (rep["Specimen ID"] || "") + "\0" + (rep["Genus"] || "");
    const gTass = genusTassByKey.has(gKey) ? genusTassByKey.get(gKey) : spTass;
    rows.forEach((r) => {
      if (r["Species TASS"] == null) r["Species TASS"] = spTass;
      if (r["Genus TASS"] == null) r["Genus TASS"] = gTass;
    });
    const sRow = Object.assign({}, rep);
    sRow["Level"] = "Species";
    sRow["Detected Organism"] = _binomial(rep["Detected Organism"]);
    sRow["Taxonomic ID #"] = rep["Subkey"] || rep["Taxonomic ID #"] || "";
    sRow["TASS Score"] = spTass;
    sRow["Species TASS"] = spTass;
    sRow["Genus TASS"] = gTass;
    sRow["# Reads Aligned"] = rows.reduce((s, r) => s + num_(r["# Reads Aligned"]), 0);
    sRow["Coverage"] = Math.max(...rows.map((r) => num_(r["Coverage"])));
    newRows.push(sRow);
  }

  // Genus proxy rows
  for (const [k, rows] of byGenus) {
    const rep = rows.reduce((a, b) => (num_(b["TASS Score"]) > num_(a["TASS Score"]) ? b : a));
    const gTass = genusTassByKey.get(k);
    const gRow = Object.assign({}, rep);
    gRow["Level"] = "Genus";
    gRow["Detected Organism"] = rep["Genus"] || _binomial(rep["Detected Organism"]);
    gRow["TASS Score"] = gTass;
    gRow["Species TASS"] = gTass;
    gRow["Genus TASS"] = gTass;
    gRow["# Reads Aligned"] = rows.reduce((s, r) => s + num_(r["# Reads Aligned"]), 0);
    gRow["Coverage"] = Math.max(...rows.map((r) => num_(r["Coverage"])));
    newRows.push(gRow);
  }

  if (newRows.length) DATA.push(...newRows);
  _HIERARCHY_SYNTHESIZED = true;
}

// Keep the View-level dropdown in sync with the data. We synthesize a
// hierarchy when needed (above), so all three levels are normally usable.
// Only options with genuinely zero rows are disabled (never the whole
// control), so the dropdown never looks dead.
function _syncViewLevelOptions() {
  _synthesizeHierarchy();
  const el = document.getElementById("view-level");
  if (!el) return;
  const have = new Set(DATA.map((r) => r["Level"] || "Strain"));
  let changed = false;
  for (const opt of el.options) {
    const present = opt.value === "Strain" ? true : have.has(opt.value);
    opt.disabled = !present;
    // Strip any stale "(no rows)" suffix, then re-add if needed.
    const base = opt.textContent.replace(/\s*\(no rows\)\s*$/, "");
    opt.textContent = present ? base : base + " (no rows)";
    if (!present && el.value === opt.value) {
      el.value = "Strain";
      changed = true;
    }
  }
  el.disabled = false;
  el.title =
    "Switch the granularity of the TABLE, heatmap and TASS/coverage plots " +
    "between strain (key), species (subkey) and genus. Rolls strains up " +
    "into their parent aggregation — most useful when a sample has lots of " +
    "strain discordance (many near-tied strains of one species). Note: it " +
    "updates the row-level tables and plots, not the sunburst/summary " +
    "genus charts, which always span all ranks.";
  if (changed) {
    _invalidateFilterCache();
    redraw();
  }
}
document.getElementById("filter-mc").addEventListener("change", redraw);

// ── Kingdom / Domain dropdown multiselect ──────────────────────────
(function () {
  const btn = document.getElementById("filter-kingdom-btn");
  const panel = document.getElementById("filter-kingdom-panel");
  const badge = document.getElementById("fk-badge");
  const clearBtn = document.getElementById("fk-clear-btn");
  const checkboxes = Array.from(document.querySelectorAll(".fk-cb"));

  function updateBtn() {
    const sel = checkboxes.filter((cb) => cb.checked);
    if (sel.length === 0) {
      btn.querySelector(".fk-label").textContent = "All kingdoms";
      btn.classList.remove("has-sel");
    } else {
      const names = sel.map((cb) => cb.closest(".fk-option").textContent.trim());
      const label = names.length <= 2 ? names.join(", ") : names.slice(0, 2).join(", ") + "\u2026";
      btn.querySelector(".fk-label").textContent = label;
      badge.textContent = sel.length;
      btn.classList.add("has-sel");
    }
    btn.setAttribute("aria-expanded", panel.classList.contains("open") ? "true" : "false");
  }

  btn.addEventListener("click", (e) => {
    e.stopPropagation();
    panel.classList.toggle("open");
    btn.setAttribute("aria-expanded", panel.classList.contains("open") ? "true" : "false");
  });

  document.addEventListener("click", (e) => {
    if (!btn.contains(e.target) && !panel.contains(e.target)) {
      panel.classList.remove("open");
      btn.setAttribute("aria-expanded", "false");
    }
  });

  checkboxes.forEach((cb) => {
    cb.addEventListener("change", () => {
      updateBtn();
      _invalidateFilterCache();
      redraw();
    });
  });

  clearBtn.addEventListener("click", () => {
    checkboxes.forEach((cb) => (cb.checked = false));
    updateBtn();
    _invalidateFilterCache();
    redraw();
  });
})();
[
  "hm-value-col",
  "hm-rank",
  "hm-scale",
  "tass-rank",
  "tass-scale",
  "tass-level",
  "tass-show-cutoff",
  "sun-metric",
  "cov-scale",
  "explore-bubble-scale",
].forEach((id) => {
  const el = document.getElementById(id);
  if (el) el.addEventListener("change", redraw);
});
// Summary tab genera-comparison controls
["summary-genus-metric", "summary-genus-topn", "summary-genus-values"].forEach((id) => {
  const el = document.getElementById(id);
  if (el)
    el.addEventListener("change", () => {
      if (typeof drawSummary === "function") drawSummary();
    });
});
// VF/AMR genus chart scale selector
const _protScaleEl = document.getElementById("prot-genus-scale");
if (_protScaleEl) _protScaleEl.addEventListener("change", drawProtGenus);
["hm-cell-color"].forEach((id) => {
  const el = document.getElementById(id);
  if (el) el.addEventListener("input", redraw);
});
