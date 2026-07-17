/* ═══════════════════════════════════════════════════════════════════════
         CROSS-SAMPLE ORGANISM ANALYSIS
         Rolls the per-(sample,organism) records up to one row per organism and
         offers four lenses over the cohort. Everything is computed client-side
         from the same filteredData() the rest of the UI uses, so it scales to
         hundreds of samples without any server round-trip.
      ═══════════════════════════════════════════════════════════════════════ */
const _XS = {
  view: "table",
  sortKey: "sampleCount",
  sortAsc: false,
  lastAgg: null,
  cmpMode: "matrix", // Feature Compare: "matrix" | "parallel"
  cmpMetrics: null, // selected metric keys (null ⇒ default set)
  cmpSort: "sampleCount", // organism ordering metric in the matrix
};
let _xsWired = false;

function _xsMean(a) {
  return a.length ? a.reduce((s, v) => s + v, 0) / a.length : 0;
}
function _xsMed(a) {
  if (!a.length) return 0;
  const s = [...a].sort((x, y) => x - y);
  const m = Math.floor(s.length / 2);
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}
// Combine a specimen's per-member values (one entry per member sample where the
// organism was seen) into a single score. `total` is the specimen's full member
// count in view, used only by the "detection" mode (fraction detected × 100).
function _xsReduceMembers(vals, method, total) {
  if (!vals.length) return 0;
  switch (method) {
    case "median":
      return _xsMed(vals);
    case "mean":
      return _xsMean(vals);
    case "min":
      return Math.min(...vals);
    case "detection":
      return total ? (vals.length / total) * 100 : 0;
    case "max":
    default:
      return Math.max(...vals);
  }
}
function _xsQuart(a, q) {
  // linear-interpolated quantile on a pre-sorted array
  if (!a.length) return 0;
  const pos = (a.length - 1) * q;
  const base = Math.floor(pos);
  const rest = pos - base;
  return a[base + 1] !== undefined ? a[base] + rest * (a[base + 1] - a[base]) : a[base];
}

// ── ANI helpers (shared by Feature Compare + capability detection) ───────
// A record "supports ANI" when match_paths.py annotated it. JSON records
// carry an explicit boolean "ANI Annotated"; flat XLSX/TSV records instead
// expose a "High ANI Matches" column (presence ⇒ ANI was computed). Older
// exports have neither, so they read as unsupported / out of date.
function _aniAnnotated(r) {
  if (!r) return false;
  if (r["ANI Annotated"] === true) return true;
  if (r["ANI Annotated"] === false) return false;
  return Object.prototype.hasOwnProperty.call(r, "High ANI Matches");
}
// Normalize a record's high-ANI partners to [{key, ani_pct}], parsing the
// JSON list form or the serialized "taxid(pct%);taxid(pct%)" string form.
function _aniMatchesFor(r) {
  if (!r) return [];
  const v = r["High ANI Matches"];
  if (Array.isArray(v)) {
    return v
      .map((m) => ({ key: String((m && m.key) || ""), ani_pct: parseFloat((m && m.ani_pct) || 0) || 0 }))
      .filter((m) => m.key);
  }
  if (typeof v === "string" && v.trim()) {
    return v
      .split(";")
      .map((s) => {
        const mm = s.trim().match(/^(.+?)\(([\d.]+)%?\)$/);
        if (mm) return { key: mm[1].trim(), ani_pct: parseFloat(mm[2]) || 0 };
        return s.trim() ? { key: s.trim(), ani_pct: 0 } : null;
      })
      .filter(Boolean);
  }
  return [];
}

// Aggregate filteredData() → one entry per organism (keyed by taxid).
// Per sample we keep the MAX TASS / coverage so a single organism that
// appears at several taxonomic levels in one sample counts once.
function _xsAggregate(fd) {
  // Distinct specimens present in the current (filtered) view — numerators
  // reference these. When specimen merge is on, DNA/RNA libraries of the same
  // specimen collapse to a single key here.
  const viewSpecimens = uniq(fd.map((r) => specimenKey(r["Specimen ID"])).filter(Boolean));
  // Prevalence denominator: specimens with ANY positive hit across the whole
  // run (Passes Threshold), independent of the live TASS slider. Falls back to
  // the in-view specimen count when the report carries no Passes-Threshold
  // flags (older data), so prevalence never divides by zero.
  const N = positiveHitSpecimenCount() || viewSpecimens.length || 1;
  const byOrg = new Map();
  // Per-specimen ANI capability: which specimens carry ANI annotation at all.
  const aniSamples = new Set();
  const noAniSamples = new Set();
  // Member samples seen (in view) per specimen — the denominator for the
  // within-specimen "detection %" aggregation.
  const specimenMembers = new Map();
  fd.forEach((r) => {
    const sn = specimenKey(r["Specimen ID"]);
    if (!sn) return;
    const member = String(r["Specimen ID"] == null ? "" : r["Specimen ID"]);
    if (!specimenMembers.has(sn)) specimenMembers.set(sn, new Set());
    if (member) specimenMembers.get(sn).add(member);
    (_aniAnnotated(r) ? aniSamples : noAniSamples).add(sn);
    const taxid = (r["Taxonomic ID #"] || "").toString();
    const key = taxid || r["Detected Organism"] || "";
    if (!key) return;
    let e = byOrg.get(key);
    if (!e) {
      e = {
        taxid: taxid || "—",
        name: r["Detected Organism"] || key,
        cat: r["Microbial Category"] || "Unknown",
        // Nested: specimen → (member sample → best TASS / coverage for that
        // member). Reduced to one value per specimen at row-build time using
        // the user-selected aggregation (max / median / mean / min / detection).
        tassMemb: new Map(),
        covMemb: new Map(),
        reads: 0,
        readsBySpec: new Map(), // specimen → summed reads (passing view)
        hc: false,
        aniMatches: new Map(), // partner taxid → max ANI %
        aniAnnotated: false,
      };
      byOrg.set(key, e);
    }
    const t = num(r["TASS Score"]);
    const c = num(r["Coverage"]);
    if (!e.tassMemb.has(sn)) e.tassMemb.set(sn, new Map());
    if (!e.covMemb.has(sn)) e.covMemb.set(sn, new Map());
    const _tm = e.tassMemb.get(sn);
    const _cm = e.covMemb.get(sn);
    if (!_tm.has(member) || t > _tm.get(member)) _tm.set(member, t);
    if (!_cm.has(member) || c > _cm.get(member)) _cm.set(member, c);
    const _rd = num(r["# Reads Aligned"]);
    e.reads += _rd;
    e.readsBySpec.set(sn, (e.readsBySpec.get(sn) || 0) + _rd);
    if (isTruthy(r["High Consequence"])) e.hc = true;
    if (_aniAnnotated(r)) e.aniAnnotated = true;
    _aniMatchesFor(r).forEach((m) => {
      if (!m.key || m.key === key) return;
      if (!e.aniMatches.has(m.key) || m.ani_pct > e.aniMatches.get(m.key)) e.aniMatches.set(m.key, m.ani_pct);
    });
  });
  const aggMethod = typeof specimenTassAgg !== "undefined" ? specimenTassAgg : "max";
  // ── Detection breakdown (pass / below / total), decoupled from the slider ──
  // fdAll re-runs the same base + specimen filters but WITHOUT the TASS gate,
  // so we can count specimens where an organism is detected yet sits below the
  // current cutoff. pass = specimens clearing the threshold (what's in fd),
  // below = detected-but-under, total = the prevalence denominator N.
  // detByOrg[key] = { tass, cov, reads } — each a Map(specimen → value) over the
  // un-thresholded data, so we can draw "all specimens (incl below)" distributions.
  const detByOrg = new Map();
  try {
    const fdAll = typeof filteredData === "function" ? filteredData({ ignoreThreshold: true }) : fd;
    fdAll.forEach((r) => {
      const sn = specimenKey(r["Specimen ID"]);
      if (!sn) return;
      const taxid = (r["Taxonomic ID #"] || "").toString();
      const dkey = taxid || r["Detected Organism"] || "";
      if (!dkey) return;
      let d = detByOrg.get(dkey);
      if (!d) {
        d = { tass: new Map(), cov: new Map(), reads: new Map() };
        detByOrg.set(dkey, d);
      }
      const _t = num(r["TASS Score"]);
      const _c = num(r["Coverage"]);
      const _rd = num(r["# Reads Aligned"]);
      if (!d.tass.has(sn) || _t > d.tass.get(sn)) d.tass.set(sn, _t);
      if (!d.cov.has(sn) || _c > d.cov.get(sn)) d.cov.set(sn, _c);
      d.reads.set(sn, (d.reads.get(sn) || 0) + _rd);
    });
  } catch (_e) {
    /* fall back to pass-only counts if the un-thresholded pass fails */
  }
  // Total specimen cohort size (all specimens, positive-hit or not) — the
  // denominator for the pass / below / total breakdown.
  const _totalSpecimens = typeof totalSpecimenCount === "function" ? totalSpecimenCount() || N : N;
  const rows = [];
  byOrg.forEach((e, _key) => {
    // Reduce each specimen's per-member values to a single score using the
    // selected aggregation. detection% divides by the specimen's total member
    // samples in view, so a hit in 1 of 2 libraries reads as 50%.
    const tassMap = new Map();
    const covMap = new Map();
    e.tassMemb.forEach((memb, sn) => {
      const total = (specimenMembers.get(sn) || new Set()).size || memb.size || 1;
      tassMap.set(sn, _xsReduceMembers([...memb.values()], aggMethod, total));
    });
    e.covMemb.forEach((memb, sn) => {
      const total = (specimenMembers.get(sn) || new Set()).size || memb.size || 1;
      // "detection" only makes sense for the TASS column; keep coverage on the
      // mean of member coverages in that mode so the number stays meaningful.
      covMap.set(sn, _xsReduceMembers([...memb.values()], aggMethod === "detection" ? "mean" : aggMethod, total));
    });
    const tass = [...tassMap.values()];
    const cov = [...covMap.values()];
    const readsVals = [...e.readsBySpec.values()];
    const _passCount = tassMap.size;
    const _detected = detByOrg.get(_key);
    const _detCount = _detected ? Math.max(_detected.tass.size, _passCount) : _passCount;
    const _belowCount = Math.max(0, _detCount - _passCount);
    // "All specimens (incl below/0)" vectors: detected values padded with zeros
    // for every specimen where the organism was not seen at all.
    const _pad = (arr, n) => {
      const a = arr.slice();
      while (a.length < n) a.push(0);
      return a;
    };
    const _detTass = _detected ? [..._detected.tass.values()] : tass.slice();
    const _detCov = _detected ? [..._detected.cov.values()] : cov.slice();
    const _detReads = _detected ? [..._detected.reads.values()] : readsVals.slice();
    const tassAll = _pad(_detTass, _totalSpecimens);
    const covAll = _pad(_detCov, _totalSpecimens);
    const readsAll = _pad(_detReads, _totalSpecimens);
    rows.push({
      taxid: e.taxid,
      name: e.name,
      cat: e.cat,
      hc: e.hc,
      sampleCount: _passCount,
      passCount: _passCount,
      belowCount: _belowCount,
      detCount: _detCount,
      total: _totalSpecimens,
      samplePct: Math.min(100, (_passCount / N) * 100),
      meanTass: _xsMean(tass),
      medianTass: _xsMed(tass),
      maxTass: tass.length ? Math.max(...tass) : 0,
      minTass: tass.length ? Math.min(...tass) : 0,
      meanCov: _xsMean(cov),
      medianCov: _xsMed(cov),
      maxCov: cov.length ? Math.max(...cov) : 0,
      minCov: cov.length ? Math.min(...cov) : 0,
      reads: e.reads,
      tassVals: tass,
      covVals: cov,
      readsVals: readsVals,
      tassAll: tassAll,
      covAll: covAll,
      readsAll: readsAll,
      tassMap: tassMap,
      covMap: covMap,
      samples: new Set(tassMap.keys()),
      aniMatches: e.aniMatches,
      aniAnnotated: e.aniAnnotated,
    });
  });
  // ── ANI grouping: union-find over organisms that share ANI ≥ threshold ──
  // Only edges where BOTH endpoints are present in this view are used, so
  // group sizes reflect what the analyst is actually looking at.
  const _present = new Set(rows.map((r) => r.taxid));
  const _parent = {};
  const _find = (x) => {
    while (_parent[x] !== undefined && _parent[x] !== x) {
      _parent[x] = _parent[_parent[x]] !== undefined ? _parent[_parent[x]] : _parent[x];
      x = _parent[x];
    }
    return x;
  };
  const _union = (a, b) => {
    _parent[a] = _parent[a] === undefined ? a : _parent[a];
    _parent[b] = _parent[b] === undefined ? b : _parent[b];
    _parent[_find(a)] = _find(b);
  };
  rows.forEach((r) => {
    if (_parent[r.taxid] === undefined) _parent[r.taxid] = r.taxid;
  });
  rows.forEach((r) => {
    r.aniMatches.forEach((pct, partner) => {
      if (_present.has(partner)) _union(r.taxid, partner);
    });
  });
  const _groupMembers = new Map();
  rows.forEach((r) => {
    const g = _find(r.taxid);
    r.aniGroup = g;
    if (!_groupMembers.has(g)) _groupMembers.set(g, []);
    _groupMembers.get(g).push(r);
  });
  rows.forEach((r) => {
    r.aniGroupSize = (_groupMembers.get(r.aniGroup) || []).length;
  });
  return {
    rows,
    totalSamples: N,
    totalSpecimensAll: _totalSpecimens,
    sampleList: viewSpecimens,
    aniSamples,
    noAniSamples,
    aniGroups: _groupMembers,
    aniSupported: aniSamples.size > 0,
  };
}

// Apply the category + min-sample controls to an aggregate.
function _xsFilterRows(agg) {
  const cat = (document.getElementById("xs-cat") || {}).value || "";
  const minP = parseInt((document.getElementById("xs-minprev") || {}).value || "1", 10);
  return agg.rows.filter((r) => (!cat || r.cat === cat) && r.sampleCount >= minP);
}

function _xsCatColor(cat) {
  return _CAT_COLORS[(cat || "Unknown").split(";")[0]] || "#607d8b";
}

// ── Specimen merge: toggle + membership editor ───────────────────────────
// Sync the toggle button label/style and the "sample(s)"/"specimen(s)" unit
// wording to the current merge state. Hides the controls when the run has no
// specimen grouping to act on (nothing would change).
function _xsUpdateSpecimenControls() {
  const has = typeof hasSpecimenGrouping === "function" ? hasSpecimenGrouping() : false;
  const tog = document.getElementById("xs-specimen-toggle");
  const edit = document.getElementById("xs-specimen-edit");
  const lbl = document.getElementById("xs-specimen-lbl");
  // The editor is always useful (you can create grouping live); the toggle is
  // only meaningful once at least one multi-sample specimen exists.
  if (tog) {
    tog.style.display = has ? "" : "none";
    tog.textContent = "Merge: " + (specimenMergeEnabled ? "On" : "Off");
    tog.style.background = specimenMergeEnabled ? "#1565c0" : "#fff";
    tog.style.color = specimenMergeEnabled ? "#fff" : "#1565c0";
  }
  if (lbl) lbl.style.display = "";
  if (edit) edit.style.display = "";
  // The combine-method dropdown only matters while merge is on (it governs how
  // multiple libraries of one specimen collapse into a single score).
  const agg = document.getElementById("xs-specimen-agg");
  if (agg) {
    agg.style.display = has && specimenMergeEnabled ? "" : "none";
    if (agg.value !== specimenTassAgg) agg.value = specimenTassAgg;
  }
  // Swap the unit noun in the summary note.
  const unit = document.getElementById("xs-unit");
  if (unit) unit.textContent = specimenMergeEnabled ? "specimen(s)" : "sample(s)";
}

function _xsToggleSpecimenMerge() {
  const next = !specimenMergeEnabled;
  if (window.specimenMerge && typeof window.specimenMerge.setEnabled === "function") {
    window.specimenMerge.setEnabled(next);
    return;
  }
  specimenMergeEnabled = next;
  if (typeof _invalidateFilterCache === "function") _invalidateFilterCache();
  if (typeof _rebuildMapMarkers === "function") _rebuildMapMarkers();
  const fd = typeof filteredData === "function" ? filteredData() : [];
  _drawCrossSample(fd, null);
}

// Build (once) and open a modal that lists every sample and the specimen it is
// assigned to. Editing a specimen name and applying re-groups the samples live
// via SPECIMEN_OVERRIDE, then re-renders with merge enabled.
function _xsOpenSpecimenEditor() {
  let back = document.getElementById("xs-specimen-modal");
  if (back) back.remove();
  back = document.createElement("div");
  back.id = "xs-specimen-modal";
  back.style.cssText =
    "position:fixed;inset:0;background:rgba(0,0,0,0.45);z-index:10000;display:flex;align-items:center;justify-content:center;";

  const groups = specimenGroups();
  // Flatten to [sample, currentSpecimen] rows, grouped for readability.
  const rowsHtml = [];
  groups.forEach((members, g) => {
    members.forEach((s) => {
      const safeS = String(s).replace(/"/g, "&quot;");
      const safeG = String(g).replace(/"/g, "&quot;");
      const isGrouped = members.length > 1;
      rowsHtml.push(
        `<tr>` +
          `<td style="padding:4px 8px;font-size:0.85em;white-space:nowrap;">${safeS}</td>` +
          `<td style="padding:4px 8px;">` +
          `<input class="xs-sp-input" data-sample="${safeS}" value="${safeG}" ` +
          `style="width:100%;box-sizing:border-box;padding:3px 6px;border:1px solid ${
            isGrouped ? "#1565c0" : "#ccc"
          };border-radius:4px;font-size:0.85em;" />` +
          `</td></tr>`,
      );
    });
  });

  back.innerHTML =
    `<div style="background:#fff;border-radius:8px;max-width:560px;width:92%;max-height:80vh;display:flex;flex-direction:column;box-shadow:0 8px 30px rgba(0,0,0,0.3);">` +
    `<div style="padding:14px 18px;border-bottom:1px solid #eee;display:flex;align-items:center;gap:8px;">` +
    `<b style="font-size:1.02em;color:#1565c0;">Edit specimens</b>` +
    `<span style="flex:1"></span>` +
    `<span style="cursor:pointer;font-size:1.3em;color:#888;" id="xs-sp-close">&times;</span>` +
    `</div>` +
    `<div style="padding:6px 18px;font-size:0.82em;color:#666;">Give samples the same <b>specimen</b> value to merge them (e.g. a DNA and RNA library from one swab). Prevalence and per-organism TASS / coverage / reads then aggregate per specimen.</div>` +
    `<div style="overflow:auto;padding:4px 18px 10px;">` +
    `<table style="width:100%;border-collapse:collapse;">` +
    `<thead><tr>` +
    `<th style="text-align:left;padding:4px 8px;font-size:0.78em;color:#888;border-bottom:1px solid #eee;">Sample</th>` +
    `<th style="text-align:left;padding:4px 8px;font-size:0.78em;color:#888;border-bottom:1px solid #eee;">Specimen</th>` +
    `</tr></thead><tbody>${rowsHtml.join("")}</tbody></table>` +
    `</div>` +
    `<div style="padding:12px 18px;border-top:1px solid #eee;display:flex;gap:8px;justify-content:flex-end;">` +
    `<button id="xs-sp-reset" style="padding:5px 12px;border:1px solid #ccc;border-radius:5px;background:#fff;color:#666;cursor:pointer;font-size:0.85em;">Reset to samplesheet</button>` +
    `<span style="flex:1"></span>` +
    `<button id="xs-sp-cancel" style="padding:5px 12px;border:1px solid #ccc;border-radius:5px;background:#fff;color:#666;cursor:pointer;font-size:0.85em;">Cancel</button>` +
    `<button id="xs-sp-apply" style="padding:5px 14px;border:1px solid #1565c0;border-radius:5px;background:#1565c0;color:#fff;cursor:pointer;font-size:0.85em;font-weight:600;">Apply</button>` +
    `</div></div>`;

  document.body.appendChild(back);

  const close = () => back.remove();
  back.addEventListener("click", (e) => {
    if (e.target === back) close();
  });
  document.getElementById("xs-sp-close").addEventListener("click", close);
  document.getElementById("xs-sp-cancel").addEventListener("click", close);

  document.getElementById("xs-sp-reset").addEventListener("click", () => {
    Object.keys(SPECIMEN_OVERRIDE).forEach((k) => delete SPECIMEN_OVERRIDE[k]);
    close();
    if (window.specimenMerge && typeof window.specimenMerge.setEnabled === "function") {
      window.specimenMerge.setEnabled(specimenMergeEnabled);
      return;
    }
    if (typeof _invalidateFilterCache === "function") _invalidateFilterCache();
    if (typeof _rebuildMapMarkers === "function") _rebuildMapMarkers();
    const fd = typeof filteredData === "function" ? filteredData() : [];
    _drawCrossSample(fd, null);
  });

  document.getElementById("xs-sp-apply").addEventListener("click", () => {
    back.querySelectorAll(".xs-sp-input").forEach((inp) => {
      const sample = inp.getAttribute("data-sample");
      const val = (inp.value || "").trim();
      // Determine the metadata default so we only keep genuine overrides.
      const meta = (typeof SAMPLE_META !== "undefined" && SAMPLE_META[sample]) || {};
      const metaSp =
        meta.specimen != null ? meta.specimen : meta.specimen_id != null ? meta.specimen_id : meta.specimen_group;
      const dflt = metaSp != null && String(metaSp).trim() ? String(metaSp).trim() : sample;
      if (!val || val === dflt) delete SPECIMEN_OVERRIDE[sample];
      else SPECIMEN_OVERRIDE[sample] = val;
    });
    // Turn merge on so the user immediately sees the effect of their grouping.
    close();
    if (window.specimenMerge && typeof window.specimenMerge.setEnabled === "function") {
      window.specimenMerge.setEnabled(true);
      return;
    }
    specimenMergeEnabled = true;
    if (typeof _invalidateFilterCache === "function") _invalidateFilterCache();
    if (typeof _rebuildMapMarkers === "function") _rebuildMapMarkers();
    const fd = typeof filteredData === "function" ? filteredData() : [];
    _drawCrossSample(fd, null);
  });
}

function _drawCrossSample(fd, samples) {
  const wrap = document.getElementById("xs-wrap");
  if (!wrap) return;
  const agg = _xsAggregate(fd);
  _XS.lastAgg = agg;
  const nsamp = document.getElementById("xs-nsamp");
  if (nsamp) nsamp.textContent = agg.totalSamples;
  // Reflect the current specimen-merge state in the toggle + swap the unit
  // wording between "sample(s)" and "specimen(s)".
  _xsUpdateSpecimenControls();

  // populate category dropdown once per render (preserve selection)
  const catSel = document.getElementById("xs-cat");
  if (catSel) {
    const cur = catSel.value;
    const cats = uniq(agg.rows.map((r) => r.cat)).sort();
    catSel.innerHTML =
      '<option value="">All</option>' +
      cats.map((c) => `<option value="${c}"${c === cur ? " selected" : ""}>${c}</option>`).join("");
  }

  if (!_xsWired) {
    _xsWired = true;
    document.querySelectorAll(".xs-subtab").forEach((b) =>
      b.addEventListener("click", () => {
        _XS.view = b.getAttribute("data-xs");
        document.querySelectorAll(".xs-subtab").forEach((x) => x.classList.toggle("active", x === b));
        _xsRenderBody();
      }),
    );
    ["xs-cat", "xs-minprev", "xs-metric", "xs-topn"].forEach((id) => {
      const el = document.getElementById(id);
      if (el) el.addEventListener("change", _xsRenderBody);
    });
    const exp = document.getElementById("xs-export");
    if (exp) exp.addEventListener("click", _xsExportCsv);
    // Specimen-merge toggle + membership editor.
    const spTog = document.getElementById("xs-specimen-toggle");
    if (spTog) spTog.addEventListener("click", _xsToggleSpecimenMerge);
    const spEdit = document.getElementById("xs-specimen-edit");
    if (spEdit) spEdit.addEventListener("click", _xsOpenSpecimenEditor);
    const spAgg = document.getElementById("xs-specimen-agg");
    if (spAgg)
      spAgg.addEventListener("change", () => {
        specimenTassAgg = spAgg.value || "max";
        // The combine method feeds the global specimen merge in filteredData(),
        // so refresh every tab (redraw invalidates the filter cache first).
        if (typeof redraw === "function") redraw();
        else _drawCrossSample(typeof filteredData === "function" ? filteredData() : [], null);
      });
    // sort handler (delegated) for the frequency table
    const body = document.getElementById("xs-body");
    if (body)
      body.addEventListener("click", (ev) => {
        const th = ev.target.closest("th[data-key]");
        if (!th) return;
        const k = th.getAttribute("data-key");
        if (_XS.sortKey === k) _XS.sortAsc = !_XS.sortAsc;
        else {
          _XS.sortKey = k;
          _XS.sortAsc = k === "name" || k === "cat";
        }
        _xsRenderBody();
      });
  }

  // toggle which controls are relevant to the active view
  _xsRenderBody();
}

// re-render only the body using the cached aggregate (controls already wired)
function _xsRenderBody() {
  const agg = _XS.lastAgg;
  if (!agg) return;
  const body = document.getElementById("xs-body");
  if (!body) return;
  // contextual control visibility
  const showTopN = _XS.view !== "table";
  const showMetric = _XS.view === "pca" || _XS.view === "cooc";
  // Feature Compare carries its own mode + metric controls inside the body.
  const tnLbl = document.getElementById("xs-topn-lbl"),
    tn = document.getElementById("xs-topn"),
    mLbl = document.getElementById("xs-metric-lbl"),
    m = document.getElementById("xs-metric"),
    exp = document.getElementById("xs-export");
  if (tnLbl) tnLbl.style.display = showTopN ? "" : "none";
  if (tn) tn.style.display = showTopN ? "" : "none";
  if (mLbl) mLbl.style.display = showMetric ? "" : "none";
  if (m) m.style.display = showMetric ? "" : "none";
  if (exp) exp.style.display = _XS.view === "table" ? "" : "none";

  const rows = _xsFilterRows(agg);
  if (!rows.length) {
    body.innerHTML = '<p style="color:#888;padding:1em">No organisms match the current filters.</p>';
    return;
  }
  if (_XS.view === "table") _xsRenderTable(rows, agg);
  else if (_XS.view === "dist") _xsRenderDist(rows, agg);
  else if (_XS.view === "pca") _xsRenderPCA(rows, agg);
  else if (_XS.view === "cooc") _xsRenderCooc(rows, agg);
  else if (_XS.view === "compare") _xsRenderCompare(rows, agg);
}

// ── View 5: multi-feature cross-comparison (matrix + parallel coords) ────
// Compares many organisms across many features at once. Default = a
// colour-coded organism × metric matrix; toggles to parallel coordinates.
// ANI grouping (from high_ani_matches) clusters near-identical references;
// when no sample in view carries ANI, the ANI feature degrades to a notice
// and the sub-tab shows a warning icon.
const _CMP_METRICS = [
  {
    key: "samplePct",
    label: "% Samples",
    fmt: (v) => v.toFixed(0) + "%",
    desc: "prevalence across samples in view",
  },
  {
    key: "sampleCount",
    label: "# Samples",
    fmt: (v) => String(Math.round(v)),
    desc: "number of samples detected in",
  },
  { key: "meanTass", label: "Mean TASS", fmt: (v) => v.toFixed(1), desc: "mean TASS across samples" },
  { key: "medianTass", label: "Median TASS", fmt: (v) => v.toFixed(1), desc: "median TASS across samples" },
  { key: "maxTass", label: "Max TASS", fmt: (v) => v.toFixed(1), desc: "best TASS in any sample" },
  { key: "meanCov", label: "Mean Cov", fmt: (v) => v.toFixed(1), desc: "mean coverage across samples" },
  { key: "maxCov", label: "Max Cov", fmt: (v) => v.toFixed(1), desc: "best coverage in any sample" },
  {
    key: "reads",
    label: "Total Reads",
    fmt: (v) => _fmtBig(v).short,
    desc: "total aligned reads across samples",
  },
  {
    key: "aniGroupSize",
    label: "ANI Group",
    fmt: (v) => String(Math.round(v)),
    desc: "# organisms sharing high ANI (in view)",
    ani: true,
  },
];
const _CMP_DEFAULT = ["samplePct", "meanTass", "maxTass", "meanCov", "reads"];

function _cmpActiveMetrics(aniSupported) {
  let keys = _XS.cmpMetrics && _XS.cmpMetrics.length ? _XS.cmpMetrics : _CMP_DEFAULT.slice();
  // Auto-include ANI Group when ANI is available and the user hasn't chosen.
  if (!_XS.cmpMetrics && aniSupported) keys = keys.concat(["aniGroupSize"]);
  return _CMP_METRICS.filter((m) => keys.includes(m.key) && (!m.ani || aniSupported));
}

// Stable, readable colour per ANI group (skips singletons → grey).
const _CMP_GROUP_PALETTE = [
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
];
function _cmpGroupColors(rows) {
  const sizes = new Map();
  rows.forEach((r) => sizes.set(r.aniGroup, (sizes.get(r.aniGroup) || 0) + 1));
  const multi = [...sizes.entries()].filter(([, n]) => n > 1).map(([g]) => g);
  const map = new Map();
  multi.forEach((g, i) => map.set(g, _CMP_GROUP_PALETTE[i % _CMP_GROUP_PALETTE.length]));
  return map; // group → colour (singletons absent)
}

function _xsRenderCompare(rows, agg) {
  const topN = parseInt((document.getElementById("xs-topn") || {}).value || "25", 10);
  const aniSupported = !!agg.aniSupported;
  const metrics = _cmpActiveMetrics(aniSupported);
  // order organisms by the chosen sort metric, then keep ANI-group members
  // adjacent so similar references read together.
  const sortKey = _XS.cmpSort in (rows[0] || {}) ? _XS.cmpSort : "sampleCount";
  const ordered = [...rows]
    .sort((a, b) => (b[sortKey] || 0) - (a[sortKey] || 0))
    .slice(0, Math.max(2, Math.min(60, topN)));
  const grpColors = _cmpGroupColors(ordered);
  // Group-aware ordering: primary by group's best sort value, then within group.
  const groupBest = new Map();
  ordered.forEach((r) => {
    const cur = groupBest.get(r.aniGroup);
    if (cur === undefined || (r[sortKey] || 0) > cur) groupBest.set(r.aniGroup, r[sortKey] || 0);
  });
  ordered.sort((a, b) => {
    const d = (groupBest.get(b.aniGroup) || 0) - (groupBest.get(a.aniGroup) || 0);
    if (d) return d;
    if (a.aniGroup !== b.aniGroup) return String(a.aniGroup).localeCompare(String(b.aniGroup));
    return (b[sortKey] || 0) - (a[sortKey] || 0);
  });

  // ── controls (mode toggle + metric chips + sort) ──────────────────────
  const aniNote = aniSupported
    ? agg.noAniSamples.size
      ? `<span style="color:#b26a00"><i class="fas fa-triangle-exclamation"></i> ${agg.noAniSamples.size} sample(s) lack ANI — excluded from ANI grouping.</span>`
      : `<span style="color:#2e7d32"><i class="fas fa-circle-check"></i> ANI grouping active.</span>`
    : `<span style="color:#b26a00"><i class="fas fa-triangle-exclamation"></i> No ANI data in view — ANI comparison unavailable. Re-run the pipeline with <code>--enable_matrix</code> to enable it. Other features still compare normally.</span>`;
  const modeBtns =
    `<div class="cmp-modebar">` +
    `<button class="cmp-mode${
      _XS.cmpMode === "matrix" ? " active" : ""
    }" data-cmpmode="matrix"><i class="fas fa-table-cells"></i> Matrix</button>` +
    `<button class="cmp-mode${
      _XS.cmpMode === "parallel" ? " active" : ""
    }" data-cmpmode="parallel"><i class="fas fa-chart-line"></i> Parallel</button>` +
    `</div>`;
  const chips = _CMP_METRICS
    .filter((m) => !m.ani || aniSupported)
    .map((m) => {
      const on = metrics.some((x) => x.key === m.key);
      return `<button class="cmp-chip${on ? " on" : ""}" data-cmpmetric="${m.key}" title="${m.desc}">${
        m.label
      }</button>`;
    })
    .join("");
  const controls =
    `<div class="cmp-controls">${modeBtns}` +
    `<span class="cmp-chip-label">Features:</span>${chips}</div>` +
    `<div class="xs-note" style="margin-top:.2em">${ordered.length} organism(s) across ${agg.totalSamples} sample(s). ${aniNote}</div>` +
    `<div id="cmp-body"></div>`;
  body_set(controls);

  if (!metrics.length) {
    document.getElementById("cmp-body").innerHTML =
      '<p style="color:#888;padding:1em">Select at least one feature to compare.</p>';
  } else if (_XS.cmpMode === "parallel") {
    _cmpRenderParallel(ordered, metrics, agg, grpColors);
  } else {
    _cmpRenderMatrix(ordered, metrics, agg, grpColors);
  }

  // wire mode + chip controls (re-render on change)
  document.querySelectorAll("#xs-body .cmp-mode").forEach((b) =>
    b.addEventListener("click", () => {
      _XS.cmpMode = b.getAttribute("data-cmpmode");
      _xsRenderBody();
    }),
  );
  document.querySelectorAll("#xs-body .cmp-chip").forEach((b) =>
    b.addEventListener("click", () => {
      const k = b.getAttribute("data-cmpmetric");
      // Materialize current selection then toggle the clicked metric.
      const cur = new Set(metrics.map((m) => m.key));
      if (cur.has(k)) cur.delete(k);
      else cur.add(k);
      _XS.cmpMetrics = [..._CMP_METRICS.map((m) => m.key).filter((mk) => cur.has(mk))];
      _xsRenderBody();
    }),
  );
}

// Matrix: rows = organisms (ANI-grouped), columns = features.
// Each column is min–max normalised independently for the colour ramp.
function _cmpRenderMatrix(orgs, metrics, agg, grpColors) {
  const host = document.getElementById("cmp-body");
  const colMin = {},
    colMax = {};
  metrics.forEach((m) => {
    const vals = orgs.map((o) => +o[m.key] || 0);
    colMin[m.key] = Math.min(...vals);
    colMax[m.key] = Math.max(...vals);
  });
  const color = d3.scaleSequential(d3.interpolateBlues);
  // Metrics on a fixed 0–100 scale — when their column has no spread (e.g. every
  // organism is detected in 100% of its samples), fall back to an ABSOLUTE
  // shade so a uniform high column reads as strong blue instead of a flat grey
  // that looks broken.
  const _BOUNDED_0_100 = new Set(["samplePct", "meanTass", "medianTass", "maxTass", "minTass", "meanCov", "maxCov"]);
  const _cellT = (m, v) => {
    const lo = colMin[m.key],
      hi = colMax[m.key];
    if (hi > lo) return (v - lo) / (hi - lo);
    return _BOUNDED_0_100.has(m.key) ? Math.max(0, Math.min(1, v / 100)) : null; // null ⇒ neutral
  };
  const cellBg = (m, v) => {
    const t = _cellT(m, v);
    return t == null ? "#eef3f8" : color(0.12 + 0.78 * t);
  };
  const txtColor = (m, v) => {
    const t = _cellT(m, v);
    return t != null && t > 0.6 ? "#fff" : "#1a2733";
  };
  let html =
    '<div class="cmp-matrix-wrap"><table class="cmp-matrix"><thead><tr>' +
    '<th class="cmp-org-h">Organism</th>' +
    metrics.map((m) => `<th title="${m.desc}">${m.label}</th>`).join("") +
    "</tr></thead><tbody>";
  let lastGroup = null;
  orgs.forEach((o, oi) => {
    const gc = grpColors.get(o.aniGroup);
    const newGroup = gc && o.aniGroup !== lastGroup;
    lastGroup = o.aniGroup;
    const band = gc
      ? `<span class="cmp-grp-band" style="background:${gc}" title="ANI group of ${o.aniGroupSize} similar reference(s)"></span>`
      : '<span class="cmp-grp-band cmp-grp-none"></span>';
    const star = o.hc ? '<span title="high-consequence" style="color:#cc0000">● </span>' : "";
    html += `<tr${newGroup ? ' class="cmp-grp-start"' : ""}>`;
    html +=
      `<td class="cmp-org cmp-org-click" data-oi="${oi}" title="Click to compare coverage of this organism across its samples">${band}<i>${star}${o.name}</i>` +
      (gc ? `<span class="cmp-grp-tag" style="border-color:${gc};color:${gc}">ANI×${o.aniGroupSize}</span>` : "") +
      ` <i class="fas fa-chart-column cmp-org-covicon"></i></td>`;
    metrics.forEach((m) => {
      const v = +o[m.key] || 0;
      // On-chart numbers removed — the colour ramp carries the relative value
      // and the exact figure (plus per-sample breakdown) is shown on hover.
      html += `<td class="cmp-cell" style="background:${cellBg(m, v)};cursor:pointer" data-oi="${oi}" data-mk="${
        m.key
      }"></td>`;
    });
    html += "</tr>";
  });
  html += "</tbody></table></div>";
  host.innerHTML = html;
  // Per-cell: hover shows the relevant sample list, click opens a full
  // per-sample table popup for that organism × metric.
  host.querySelectorAll(".cmp-cell").forEach((td) => {
    const o = orgs[+td.getAttribute("data-oi")];
    const m = _CMP_METRICS.find((x) => x.key === td.getAttribute("data-mk"));
    if (!o || !m) return;
    td.addEventListener("mouseover", (ev) => showTip(_cmpSampleTip(o, m, agg), ev));
    td.addEventListener("mousemove", moveTip);
    td.addEventListener("mouseout", hideTip);
    td.addEventListener("click", () => {
      hideTip();
      _xsOpenSamplePopup(o, m, agg);
    });
  });
  // Organism label: open the coverage comparison overlay for that organism
  // across every sample it was detected in.
  host.querySelectorAll(".cmp-org-click").forEach((td) => {
    const o = orgs[+td.getAttribute("data-oi")];
    if (!o) return;
    td.addEventListener("click", () => _openCoverageModalForOrg(o, agg));
  });
}

// Parallel coordinates: one polyline per organism across normalised axes.
function _cmpRenderParallel(orgs, metrics, agg, grpColors) {
  const host = document.getElementById("cmp-body");
  if (metrics.length < 2) {
    host.innerHTML = '<p style="color:#888;padding:1em">Select at least two features for parallel coordinates.</p>';
    return;
  }
  host.innerHTML = '<div id="cmp-par-svg" style="overflow-x:auto"></div>';
  const wrap = document.getElementById("cmp-par-svg");
  const W = Math.max(metrics.length * 130, wrap.clientWidth || 760);
  const H = 360,
    mT = 26,
    mB = 54,
    mL = 28,
    mR = 28;
  const innerH = H - mT - mB;
  const axX = d3
    .scalePoint()
    .domain(metrics.map((m) => m.key))
    .range([mL, W - mR]);
  const scales = {};
  metrics.forEach((m) => {
    const vals = orgs.map((o) => +o[m.key] || 0);
    const lo = Math.min(...vals),
      hi = Math.max(...vals);
    scales[m.key] = d3
      .scaleLinear()
      .domain(lo === hi ? [lo - 1, hi + 1] : [lo, hi])
      .range([mT + innerH, mT]);
  });
  const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  // axes
  metrics.forEach((m) => {
    const gx = axX(m.key);
    svg
      .append("line")
      .attr("x1", gx)
      .attr("x2", gx)
      .attr("y1", mT)
      .attr("y2", mT + innerH)
      .attr("stroke", "#cfd8dc");
    const sc = scales[m.key];
    sc.ticks(4).forEach((t) => {
      svg
        .append("text")
        .attr("x", gx - 5)
        .attr("y", sc(t))
        .attr("text-anchor", "end")
        .attr("dominant-baseline", "middle")
        .style("font-size", "8px")
        .style("fill", "#90a4ae")
        .text(m.fmt(t));
    });
    svg
      .append("text")
      .attr("x", gx)
      .attr("y", mT + innerH + 18)
      .attr("text-anchor", "middle")
      .style("font-size", "9px")
      .style("font-weight", "600")
      .style("fill", "#455a64")
      .text(m.label);
  });
  const line = (o) => metrics.map((m, i) => `${i ? "L" : "M"}${axX(m.key)},${scales[m.key](+o[m.key] || 0)}`).join(" ");
  orgs.forEach((o) => {
    const gc = grpColors.get(o.aniGroup) || _xsCatColor(o.cat);
    const path = svg
      .append("path")
      .attr("d", line(o))
      .attr("fill", "none")
      .attr("stroke", gc)
      .attr("stroke-width", 1.4)
      .attr("stroke-opacity", 0.55)
      .style("cursor", "pointer");
    const tip =
      `<b><i>${o.name}</i></b><br>` +
      metrics.map((m) => `${m.label}: <b>${m.fmt(+o[m.key] || 0)}</b>`).join("<br>") +
      (grpColors.get(o.aniGroup) ? `<br>ANI group of ${o.aniGroupSize}` : "");
    path
      .on("mouseover", (ev) => {
        path.attr("stroke-width", 3).attr("stroke-opacity", 1);
        showTip(tip, ev);
      })
      .on("mousemove", moveTip)
      .on("mouseout", () => {
        path.attr("stroke-width", 1.4).attr("stroke-opacity", 0.55);
        hideTip();
      });
  });
}
