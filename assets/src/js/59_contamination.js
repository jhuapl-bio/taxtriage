/* ═══════════════════════════════════════════════════════════════════════════
       -  §  CROSSOVER CHECK  (run-wide, per-organism meta-analysis; flag only)
       -     In a run with many samples, an organism that is abundant in one
       -     sample often turns up at a low level in others: splash or aerosol
       -     at extraction / library prep, carry-over on the instrument,
       -     barcode bleed. Reagent / kit background shows a different shape:
       -     present in many samples AND in the negative control. A low-level
       -     detection may still be real; this check scores how well
       -     contamination explains it, so it gets a second look.
       -
       -     Nothing to load or upload: it runs in the browser over the
       -     detections already in the report (DATA, every level, below-cutoff
       -     detections included), so it works on any report this viewer
       -     opens, older ones included. Missing columns are skipped (see
       -     "Data it needs" below). TASS, calls and the data are never changed.
       -
       -     For every organism (Level + Taxonomic ID #, else name) the samples
       -     that carry it are compared with each other. Per detection in a
       -     sample R, four evidence terms in [0, 1]:
       -       ab    abundance vs the strongest OTHER sample S (depth-normalised,
       -             reads per million input): r = rpm_R / rpm_S;
       -             ab = log10(1/r) / log10(fullFold), 0 when r ≥ 1
       -       prev  run prevalence: share of the other samples carrying it,
       -             reaching 1 at prevFull %
       -       cov   thin coverage: 1 − breadth / breadthRef        (if known)
       -       neg   present in a negative control at ≥ 10 % of R's level
       -             (depth-normalised)                             (if any NC)
       -     Two pattern scores share the weights (weighted means in [0, 1];
       -     terms without data drop out):
       -       crossover   = mean(ab, prev, cov)   — a much stronger sample nearby
       -       background  = mean(neg, prev, cov)  — kit / reagent / environment
       -     Flagged as crossover when its score ≥ cutoff, S has ≥ minSourceReads,
       -     r ≤ 1/minFold and the organism is in ≥ minSamples samples; else as
       -     background when its score ≥ cutoff and neg ≥ 0.5. The score shown
       -     is the flagged pattern's, or the higher of the two.
       -     Libraries of one specimen are never each other's source; controls
       -     are never flagged (positive controls can be a source, negative
       -     controls are evidence); in-silico samples are left out. "Compare
       -     within" limits the comparison to samples that share a metadata
       -     value (run, plate, batch …); negative controls count everywhere.
       -
       -     Data it needs: "# Reads Aligned" (required); depth from
       -     SAMPLE_META[s].total_reads, else inferred from the RPM column,
       -     else the sample's aligned reads; breadth from "Breadth %" or
       -     "Coverage"; controls from SAMPLE_META[s].control_type.
       -
       -     Engine: ttCtCompute (pure; scripts/test_contamination.js).
       -     Lookups: ttCtFor(row) / ttCtSourceFor(row) by unit + level + taxon,
       -     so row copies from filteredData() (rescue, specimen merge) resolve.
       -     Views: the Summary → Crossover sub-tab (ttCtRender), markers in
       -     the Summary detections / Table (ttCtBadgeHTML) and heatmap
       -     (ttCtHeatmapMarks) while marking is on, the filteredData() gate
       -     (ttCtFilterKey / ttCtRowHidden) for Only / Hide flagged, and the
       -     banner pill (#banner-ct) as the quick on / off switch.
═══════════════════════════════════════════════════════════════════════════ */

const TT_CT_DEFAULTS = {
  cutoff: 0.5, // flag at score ≥ this
  wAbund: 50, // weights (relative; shown as shares)
  wPrev: 15,
  wCov: 15,
  wNeg: 20,
  fullFold: 1000, // abundance term reaches 1 at 1:fullFold of the top sample
  minFold: 10, // crossover needs the top sample ≥ minFold × higher
  prevFull: 50, // prevalence term reaches 1 when this % of the other samples carry it
  breadthRef: 20, // breadth (%) treated as well covered
  minSourceReads: 50, // a source needs at least this many reads
  minSamples: 2, // crossover needs the organism in ≥ this many samples
};
const TT_CT_PRESETS = {
  sensitive: { cutoff: 0.4, fullFold: 100, minFold: 5, minSourceReads: 20 },
  balanced: { cutoff: 0.5, fullFold: 1000, minFold: 10, minSourceReads: 50 },
  strict: { cutoff: 0.65, fullFold: 10000, minFold: 50, minSourceReads: 200 },
};
const TT_CT_NC_REF = 0.1; // NC at ≥ 10 % of the sample's level → neg term = 1
const TT_CT_STORE_KEY = "taxtriage:contamination:v1";
const _TT_CT_SEP = "\u0001";

const TT_CT = {
  on: false, // mark flagged detections across the report
  view: "all", // "all" | "only" | "hide" (report filter; not saved in the browser)
  p: Object.assign({}, TT_CT_DEFAULTS),
  scope: "", // metadata column the comparison stays within ("" = whole run)
  excluded: new Set(), // units left out of the pool (session only)
  ui: { q: "", sort: "flags", only: true, topN: 30, focus: null },
  _cache: { key: null, value: null },
  _rev: 0,
};

/* ── small helpers ───────────────────────────────────────────────────────── */
function _ttCtEsc(v) {
  return String(v == null ? "" : v)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}
function _ttCtFmt(v, d) {
  const n = Number(v);
  return isFinite(n) ? n.toLocaleString(undefined, { maximumFractionDigits: d == null ? 1 : d }) : "—";
}
function _ttCtClamp01(v) {
  return v > 1 ? 1 : v > 0 ? v : 0;
}
/** "1:340", "1:12k", "1:2.5M" for a ratio r (< 1); "≥ top" otherwise. */
function _ttCtFold(r) {
  if (r == null || !isFinite(r)) return "—";
  if (r >= 1) return "top";
  const f = 1 / Math.max(r, 1e-12);
  const s =
    f >= 1e6
      ? _ttCtFmt(f / 1e6, f < 1e7 ? 1 : 0) + "M"
      : f >= 1e4
      ? _ttCtFmt(f / 1e3, 0) + "k"
      : _ttCtFmt(f, f < 10 ? 1 : 0);
  return "1:" + s;
}
function _ttCtLevel(r) {
  return String((r && r.Level) || "Strain").toLowerCase();
}
/** Taxon identity of a row: the taxid, else the organism name (older / TSV rows). */
function _ttCtTaxKey(r) {
  const t = String(r["Taxonomic ID #"] == null ? "" : r["Taxonomic ID #"]).trim();
  if (t) return t;
  const n = String(r["Detected Organism"] == null ? "" : r["Detected Organism"])
    .trim()
    .toLowerCase();
  return n ? "name:" + n : "";
}
function _ttCtBreadth(r) {
  for (const c of ["Breadth %", "Coverage"]) {
    const v = r[c];
    if (v === null || v === undefined || v === "") continue;
    const n = Number(v);
    if (isFinite(n) && n >= 0) return Math.min(100, n);
  }
  return null;
}
function _ttCtIsNegType(t) {
  return /^(neg|ntc|blank|nc$|no.?template|extraction.?blank)/.test(String(t || "").toLowerCase());
}
function _ttCtHash(s) {
  let h = 2166136261;
  s = String(s);
  for (let i = 0; i < s.length; i++) {
    h ^= s.charCodeAt(i);
    h = Math.imul(h, 16777619);
  }
  return ((h >>> 0) % 10007) / 10007; // [0, 1)
}

/* ── the analysis (pure) ─────────────────────────────────────────────────── */
/** rows: detection rows. opts: { params, unitOf(s), specimenOf(s), groupOf(s),
 *  skip(s), ctlType(s), depthOf(s) } — all optional.
 *  Returns { recs, orgs, units, sources, norm, hasNeg, hasBreadth, weights, params }:
 *    recs     "unit␁level␁tax" → one detection's record (score, terms, verdict)
 *    orgs     "level␁tax" → { name, recs[], nWith, nPool, nFlag, nCross, nBg, top … }
 *    units    unit → { members, ctl, group, nDet{lvl}, nFlag{lvl}, nSrc{lvl} }
 *    sources  "unit␁level␁tax" of a source → { n, recipients[] } */
function ttCtCompute(rows, opts) {
  const o = opts || {};
  const P = Object.assign({}, TT_CT_DEFAULTS, o.params || {});
  const unitOf = o.unitOf || ((s) => s);
  const specOf = o.specimenOf || ((s) => s);
  const groupOf = o.groupOf || (() => "");
  const skip = o.skip || (() => false);
  const ctlType = o.ctlType || (() => "");
  const depthOf = o.depthOf || (() => null);
  const SEP = _TT_CT_SEP;

  // 1. per sample: the strongest row of each (level, taxon)
  const perSample = new Map();
  (rows || []).forEach((r) => {
    if (!r) return;
    const s = String(r["Specimen ID"] == null ? "" : r["Specimen ID"]);
    if (!s || skip(s)) return;
    const tk = _ttCtTaxKey(r);
    if (!tk) return;
    const lvl = _ttCtLevel(r);
    const k = lvl + SEP + tk;
    const reads = Math.max(0, Math.trunc(Number(r["# Reads Aligned"]) || 0));
    let S = perSample.get(s);
    if (!S) perSample.set(s, (S = { taxa: new Map(), aligned: 0, rpmDepth: [] }));
    if (lvl === "strain") S.aligned += reads;
    const rpm = Number(r.RPM);
    if (reads > 0 && isFinite(rpm) && rpm > 0) S.rpmDepth.push((reads * 1e6) / rpm);
    const cur = S.taxa.get(k);
    if (!cur || reads > cur.reads)
      S.taxa.set(k, { reads, breadth: _ttCtBreadth(r), name: r["Detected Organism"] || "" });
  });

  // 2. sequencing depth per sample (for depth-normalised abundance)
  const RANK = { input: 0, rpm: 1, aligned: 2 };
  let norm = "input";
  perSample.forEach((S, s) => {
    let d = Number(depthOf(s));
    let how = "input";
    if (!(d > 0) && S.rpmDepth.length) {
      const a = S.rpmDepth.sort((x, y) => x - y);
      d = a[a.length >> 1];
      how = "rpm";
    }
    if (!(d > 0)) {
      d = S.aligned;
      if (!(d > 0)) S.taxa.forEach((t) => (d = Math.max(d || 0, t.reads)));
      how = "aligned";
    }
    S.depth = d > 0 ? d : 1;
    if (RANK[how] > RANK[norm]) norm = how;
  });

  // 3. units (a sample, or a merged specimen = the sum of its libraries)
  const units = new Map();
  perSample.forEach((S, s) => {
    const id = String(unitOf(s));
    let U = units.get(id);
    if (!U) {
      U = {
        id,
        specimen: String(specOf(s)),
        group: String(groupOf(s) || ""),
        members: [],
        depth: 0,
        taxa: new Map(),
        ctl: "",
        nDet: {},
        nFlag: {},
        nCross: {},
        nBg: {},
        nSrc: {},
      };
      units.set(id, U);
    }
    U.members.push(s);
    U.depth += S.depth;
    const ct = String(ctlType(s) || "").toLowerCase();
    if (ct) U.ctl = _ttCtIsNegType(ct) || U.ctl === "negative" ? "negative" : U.ctl || "positive";
    S.taxa.forEach((t, k) => {
      const c = U.taxa.get(k);
      if (!c) U.taxa.set(k, { reads: t.reads, breadth: t.breadth, name: t.name });
      else {
        c.reads += t.reads;
        if (t.breadth != null && (c.breadth == null || t.breadth > c.breadth)) c.breadth = t.breadth;
      }
    });
  });

  // 4. who carries what
  const nPool = new Map(); // group → regular (non-control) units
  let hasNeg = false;
  let hasBreadth = false;
  const lists = new Map(); // group␁key → { k, group, items }
  const ncMax = new Map(); // key → strongest negative control
  const topAll = new Map(); // key → strongest non-negative unit, run-wide
  units.forEach((U) => {
    if (U.ctl === "negative") hasNeg = true;
    if (!U.ctl) nPool.set(U.group, (nPool.get(U.group) || 0) + 1);
    U.taxa.forEach((t, k) => {
      t.rpm = (t.reads / U.depth) * 1e6;
      if (t.breadth != null) hasBreadth = true;
      const lvl = k.slice(0, k.indexOf(SEP));
      if (!U.ctl) U.nDet[lvl] = (U.nDet[lvl] || 0) + 1;
      if (t.reads <= 0) return;
      if (U.ctl === "negative") {
        const c = ncMax.get(k);
        if (!c || t.rpm > c.rpm) ncMax.set(k, { rpm: t.rpm, unit: U.id, reads: t.reads });
        return;
      }
      const a = topAll.get(k);
      if (!a || t.rpm > a.rpm) topAll.set(k, { rpm: t.rpm, unit: U.id, reads: t.reads });
      const gk = U.group + SEP + k;
      let L = lists.get(gk);
      if (!L) lists.set(gk, (L = { k, group: U.group, items: [] }));
      L.items.push({ U, t });
    });
  });

  const W = {
    ab: Math.max(0, Number(P.wAbund) || 0),
    prev: Math.max(0, Number(P.wPrev) || 0),
    cov: hasBreadth ? Math.max(0, Number(P.wCov) || 0) : 0,
    neg: hasNeg ? Math.max(0, Number(P.wNeg) || 0) : 0,
  };
  const lFull = Math.log10(Math.max(1.001, Number(P.fullFold) || TT_CT_DEFAULTS.fullFold));
  const minFold = Math.max(1, Number(P.minFold) || 1);
  const prevFull = Math.max(0.01, (Number(P.prevFull) || TT_CT_DEFAULTS.prevFull) / 100);
  const bRef = Math.max(0.01, Number(P.breadthRef) || TT_CT_DEFAULTS.breadthRef);
  const minSrc = Math.max(0, Number(P.minSourceReads) || 0);
  const minSamples = Math.max(2, Math.trunc(Number(P.minSamples) || 2));
  const cutoff = Number(P.cutoff);

  const recs = new Map();
  const orgs = new Map();
  const sent = new Map();
  const orgOf = (k, name) => {
    let O = orgs.get(k);
    if (!O) {
      const i = k.indexOf(SEP);
      const ta = topAll.get(k);
      const nc = ncMax.get(k);
      O = {
        key: k,
        level: k.slice(0, i),
        tax: k.slice(i + 1),
        name,
        recs: [],
        nWith: 0,
        nPool: 0,
        nFlag: 0,
        nCross: 0,
        nBg: 0,
        maxScore: 0,
        top: ta ? ta.unit : null,
        topRpm: ta ? ta.rpm : 0,
        topReads: ta ? ta.reads : 0,
        ncUnit: nc ? nc.unit : null,
        ncRpm: nc ? nc.rpm : 0,
      };
      orgs.set(k, O);
    }
    return O;
  };

  lists.forEach((L) => {
    const k = L.k;
    const lvl = k.slice(0, k.indexOf(SEP));
    const items = L.items;
    const O = orgOf(k, items[0].t.name);
    const pool = nPool.get(L.group) || 0;
    const nWith = items.reduce((n, it) => n + (it.U.ctl ? 0 : 1), 0);
    const nAll = items.length;
    O.nWith += nWith;
    O.nPool += pool;
    const srcs = items.filter((it) => it.t.reads >= minSrc).sort((a, b) => b.t.rpm - a.t.rpm);
    const nc = ncMax.get(k);
    items.forEach((it) => {
      const U = it.U,
        t = it.t;
      let top = null;
      for (const s of srcs) {
        if (s.U !== U && s.U.specimen !== U.specimen) {
          top = s;
          break;
        }
      }
      const r = top && top.t.rpm > 0 ? t.rpm / top.t.rpm : null;
      const rec = {
        unit: U.id,
        key: k,
        level: lvl,
        name: t.name,
        group: L.group,
        role: U.ctl === "positive" ? "pos" : "sample",
        reads: t.reads,
        rpm: t.rpm,
        breadth: t.breadth,
        ratio: r,
        top: top ? top.U.id : null,
        topReads: top ? top.t.reads : null,
        topRpm: top ? top.t.rpm : null,
        nWith,
        nAll,
        nPool: pool,
        ncUnit: nc ? nc.unit : null,
        ncRpm: nc ? nc.rpm : 0,
        terms: null,
        score: 0,
        suspect: false,
        kind: null,
      };
      if (rec.role === "sample") {
        const c = {
          ab: r != null && r < 1 ? _ttCtClamp01(Math.log10(1 / r) / lFull) : 0,
          prev: pool > 1 ? _ttCtClamp01((nWith - 1) / (pool - 1) / prevFull) : 0,
          cov: W.cov > 0 && t.breadth != null ? _ttCtClamp01(1 - t.breadth / bRef) : null,
          neg: hasNeg ? (nc && nc.rpm > 0 && t.rpm > 0 ? _ttCtClamp01(nc.rpm / (TT_CT_NC_REF * t.rpm)) : 0) : null,
        };
        const mean = (names) => {
          let sw = 0,
            s = 0;
          for (const n of names) {
            if (W[n] > 0 && c[n] != null) {
              sw += W[n];
              s += W[n] * c[n];
            }
          }
          return sw > 0 ? s / sw : 0;
        };
        rec.terms = c;
        rec.scores = { crossover: mean(["ab", "prev", "cov"]), background: hasNeg ? mean(["neg", "prev", "cov"]) : 0 };
        const crossOk =
          t.reads > 0 && rec.scores.crossover >= cutoff && r != null && r <= 1 / minFold && nAll >= minSamples;
        const bgOk = t.reads > 0 && hasNeg && W.neg > 0 && rec.scores.background >= cutoff && c.neg >= 0.5;
        rec.suspect = crossOk || bgOk;
        rec.kind = crossOk ? "crossover" : bgOk ? "background" : null;
        rec.score = rec.kind ? rec.scores[rec.kind] : Math.max(rec.scores.crossover, rec.scores.background);
        if (rec.score > O.maxScore) O.maxScore = rec.score;
        if (rec.suspect) {
          O.nFlag++;
          U.nFlag[lvl] = (U.nFlag[lvl] || 0) + 1;
          if (rec.kind === "crossover") {
            O.nCross++;
            U.nCross[lvl] = (U.nCross[lvl] || 0) + 1;
            const dk = top.U.id + SEP + k;
            let list = sent.get(dk);
            if (!list) {
              sent.set(dk, (list = []));
              top.U.nSrc[lvl] = (top.U.nSrc[lvl] || 0) + 1;
            }
            list.push({ unit: U.id, ratio: r, reads: t.reads, score: rec.score });
          } else {
            O.nBg++;
            U.nBg[lvl] = (U.nBg[lvl] || 0) + 1;
          }
        }
      }
      recs.set(U.id + SEP + k, rec);
      O.recs.push(rec);
    });
  });

  // negative controls: plotted as evidence, never scored
  units.forEach((U) => {
    if (U.ctl !== "negative") return;
    U.taxa.forEach((t, k) => {
      if (t.reads <= 0) return;
      const ta = topAll.get(k);
      const O = orgOf(k, t.name);
      const rec = {
        unit: U.id,
        key: k,
        level: k.slice(0, k.indexOf(SEP)),
        name: t.name,
        group: U.group,
        role: "neg",
        reads: t.reads,
        rpm: t.rpm,
        breadth: t.breadth,
        ratio: ta && ta.rpm > 0 ? t.rpm / ta.rpm : null,
        top: ta ? ta.unit : null,
        topReads: ta ? ta.reads : null,
        topRpm: ta ? ta.rpm : null,
        score: 0,
        suspect: false,
        kind: null,
      };
      recs.set(U.id + SEP + k, rec);
      O.recs.push(rec);
    });
  });

  const sources = new Map();
  sent.forEach((list, dk) => {
    list.sort((a, b) => b.reads - a.reads);
    sources.set(dk, { n: list.length, recipients: list.slice(0, 6) });
  });
  return { recs, orgs, units, sources, norm, hasNeg, hasBreadth, weights: W, params: P };
}

/* ── state: the report's rows → a cached result ──────────────────────────── */
function _ttCtMetaValue(sample, col) {
  if (!col) return "";
  let rec = null;
  if (typeof RUN_META !== "undefined" && Array.isArray(RUN_META)) {
    if (!_ttCtMetaValue._idx || _ttCtMetaValue._n !== RUN_META.length) {
      const idx = {};
      RUN_META.forEach((r) => {
        if (r && r.sample_name != null) idx[r.sample_name] = r;
      });
      _ttCtMetaValue._idx = idx;
      _ttCtMetaValue._n = RUN_META.length;
    }
    rec = _ttCtMetaValue._idx[sample] || null;
  }
  let v = rec ? rec[col] : null;
  if (v == null && typeof SAMPLE_META !== "undefined" && SAMPLE_META && SAMPLE_META[sample])
    v = SAMPLE_META[sample][col];
  if (v == null) return "";
  return (Array.isArray(v) ? v : [v])
    .map((x) => String(x).trim())
    .filter(Boolean)
    .sort()
    .join("; ");
}
function _ttCtUnitOf(sample) {
  return typeof specimenKey === "function" ? specimenKey(sample) : String(sample == null ? "" : sample);
}
function _ttCtIsInsilico(sample) {
  const m = (typeof SAMPLE_META !== "undefined" && SAMPLE_META && SAMPLE_META[sample]) || null;
  return !!(m && m.insilico === true);
}
function _ttCtCtlType(sample) {
  if (typeof ttControlType === "function") return ttControlType(sample);
  const m = (typeof SAMPLE_META !== "undefined" && SAMPLE_META && SAMPLE_META[sample]) || {};
  return m.control_type || "";
}
function _ttCtDepthOf(sample) {
  const m = (typeof SAMPLE_META !== "undefined" && SAMPLE_META && SAMPLE_META[sample]) || null;
  if (!m) return null;
  const v = Number(m.total_reads != null ? m.total_reads : m.input_data && m.input_data.reads);
  return isFinite(v) && v > 0 ? v : null;
}
/** Cheap: called for every row by the markers and the filteredData() gate.
 *  Every settings change bumps TT_CT._rev (_ttCtChanged / ttCtInvalidate). */
function _ttCtStateKey() {
  const merge = typeof specimenMergeEnabled !== "undefined" && specimenMergeEnabled;
  return (
    TT_CT._rev +
    "|" +
    (typeof DATA !== "undefined" && DATA ? DATA.length : 0) +
    "|" +
    (typeof SAMPLE_META_EPOCH !== "undefined" ? SAMPLE_META_EPOCH : 0) +
    "|" +
    (merge && typeof _hashSpecimenMerge === "function"
      ? _hashSpecimenMerge()
      : typeof SPECIMEN_OVERRIDE !== "undefined" && SPECIMEN_OVERRIDE
      ? "0:" + Object.keys(SPECIMEN_OVERRIDE).length
      : "0")
  );
}
/** The analysis for the current settings (cached), or null when there is no data. */
function ttCtResult() {
  if (typeof DATA === "undefined" || !DATA || !DATA.length) return null;
  const key = _ttCtStateKey();
  if (TT_CT._cache.key === key) return TT_CT._cache.value;
  let value = null;
  try {
    const scope = TT_CT.scope;
    const ex = TT_CT.excluded;
    value = ttCtCompute(DATA, {
      params: TT_CT.p,
      unitOf: _ttCtUnitOf,
      specimenOf: typeof specimenOf === "function" ? specimenOf : null,
      groupOf: scope ? (s) => _ttCtMetaValue(s, scope) : null,
      skip: (s) => _ttCtIsInsilico(s) || (ex.size > 0 && (ex.has(s) || ex.has(_ttCtUnitOf(s)))),
      ctlType: _ttCtCtlType,
      depthOf: _ttCtDepthOf,
    });
  } catch (e) {
    if (typeof console !== "undefined" && console.error) console.error("[crossover] compute failed", e);
    value = null;
  }
  TT_CT._cache = { key, value };
  return value;
}
function ttCtInvalidate() {
  TT_CT._rev++;
  TT_CT._cache = { key: null, value: null };
}
function _ttCtLookup(map, r) {
  if (!r || !map) return null;
  const sid = String(r["Specimen ID"] == null ? "" : r["Specimen ID"]);
  const k = _ttCtLevel(r) + _TT_CT_SEP + _ttCtTaxKey(r);
  const unit = _ttCtUnitOf(sid);
  return map.get(unit + _TT_CT_SEP + k) || (sid !== unit ? map.get(sid + _TT_CT_SEP + k) : null) || null;
}
/** The contamination record of a row (any copy of it), or null. */
function ttCtFor(r) {
  const res = ttCtResult();
  return res ? _ttCtLookup(res.recs, r) : null;
}
function ttCtIsSuspect(r) {
  const x = ttCtFor(r);
  return !!(x && x.suspect);
}
/** Source side: samples this row's organism may have crossed into, or null. */
function ttCtSourceFor(r) {
  const res = ttCtResult();
  return res ? _ttCtLookup(res.sources, r) : null;
}

/* ── filteredData() gate: Only flagged / Hide flagged ────────────────────── */
function ttCtFilterKey() {
  if (TT_CT.view === "all") return "";
  return "ct:" + TT_CT.view + "|" + _ttCtStateKey();
}
function ttCtRowHidden(r) {
  if (TT_CT.view === "all") return false;
  const sus = ttCtIsSuspect(r);
  return TT_CT.view === "only" ? !sus : sus;
}

/* ── markers + hover cards (Summary detections, Table) ───────────────────── */
function _ttCtMarkAttrs(r, kind) {
  const a = (v) => _ttCtEsc(v == null ? "" : v);
  return (
    `data-ct-tip="${kind}" data-ct-s="${a(r["Specimen ID"])}" data-ct-l="${a(r.Level || "Strain")}" ` +
    `data-ct-k="${a(r["Taxonomic ID #"])}" data-ct-n="${a(r["Detected Organism"])}"`
  );
}
/** Markers after the organism name. "" while marking is off. */
function ttCtBadgeHTML(r) {
  if (!TT_CT.on || !r) return "";
  let html = "";
  const x = ttCtFor(r);
  if (x && x.suspect) {
    const bg = x.kind === "background";
    html += ` <span class="ct-mark${bg ? " ct-bg" : ""}" ${_ttCtMarkAttrs(r, "sus")} tabindex="0" aria-label="${
      bg ? "possible background (also in the negative control)" : "possible crossover from another sample"
    }"><i class="fas ${bg ? "fa-flask" : "fa-right-left"}" aria-hidden="true"></i> ${
      bg ? "background?" : "crossover?"
    }</span>`;
  }
  if (ttCtSourceFor(r))
    html += ` <i class="fas fa-arrow-right-from-bracket ct-src-mark" ${_ttCtMarkAttrs(
      r,
      "src",
    )} tabindex="0" aria-label="possible crossover source"></i>`;
  return html;
}
function _ttCtTermsTable(x) {
  if (!x || !x.terms) return "";
  const res = ttCtResult();
  const W = (res && res.weights) || {};
  // Shares within the pattern shown (crossover: ab/prev/cov; background: neg/prev/cov).
  const pat = x.kind || (x.scores && x.scores.background > x.scores.crossover ? "background" : "crossover");
  const inPat = (n) => (pat === "background" ? n !== "ab" : n !== "neg");
  let sw = 0;
  ["ab", "prev", "cov", "neg"].forEach((n) => {
    if (inPat(n) && W[n] > 0 && x.terms[n] != null) sw += W[n];
  });
  const lab = {
    ab: "Abundance vs top",
    prev: "Run prevalence",
    cov: "Thin coverage",
    neg: "In negative control",
  };
  const why = {
    ab: x.ratio != null ? _ttCtFold(x.ratio) + (x.top ? " of " + _ttCtEsc(x.top) : "") : "no other source",
    prev: `${x.nWith}/${x.nPool} samples`,
    cov: x.breadth != null ? _ttCtFmt(x.breadth, 1) + "% breadth" : "",
    neg:
      x.ncRpm > 0
        ? `${_ttCtEsc(x.ncUnit)} at ${_ttCtFmt((100 * x.ncRpm) / Math.max(x.rpm, 1e-9), 0)}% of here`
        : "not in NC",
  };
  const rows = ["ab", "prev", "cov", "neg"]
    .filter((n) => x.terms[n] != null && W[n] > 0)
    .map((n) => {
      const share = inPat(n) && sw > 0 ? W[n] / sw : 0;
      const v = x.terms[n];
      return (
        `<tr><td style="color:#90caf9;padding-right:8px;white-space:nowrap">${lab[n]}</td>` +
        `<td style="padding-right:8px"><span class="ct-tbar"><span style="width:${Math.round(
          v * 100,
        )}%"></span></span> ${v.toFixed(2)}</td>` +
        `<td style="color:#b0bec5;padding-right:8px">${inPat(n) ? "× " + Math.round(share * 100) + "%" : "—"}</td>` +
        `<td style="color:#b0bec5">${why[n]}</td></tr>`
      );
    })
    .join("");
  const sc = x.scores || {};
  return (
    `<table style="border-collapse:collapse;font-size:0.85em;margin:5px 0">${rows}</table>` +
    `<div style="font-size:0.8em;color:#b0bec5">Pattern scores: crossover ${(sc.crossover || 0).toFixed(2)}` +
    (res && res.hasNeg ? ` · background ${(sc.background || 0).toFixed(2)}` : "") +
    ` (weights shown for ${pat})</div>`
  );
}
function ttCtSuspectTipHTML(r, x) {
  x = x || ttCtFor(r);
  if (!x) return "";
  const esc = _ttCtEsc;
  const bg = x.kind === "background";
  const head = x.suspect
    ? bg
      ? `<b style="color:#b39ddb"><i class="fas fa-flask"></i> Possible background (negative control)</b>`
      : `<b style="color:#90caf9"><i class="fas fa-right-left"></i> Possible crossover from another sample</b>`
    : `<b style="color:#90caf9"><i class="fas fa-vials"></i> Crossover check</b>`;
  const lead = x.suspect
    ? bg
      ? `${esc(
          x.name,
        )} is also in the negative control at a comparable level, the pattern of reagent / kit or environmental background.`
      : `${esc(x.name)} is ${_ttCtFold(x.ratio).replace("1:", "")}× more abundant in ${esc(
          x.top,
        )} (depth-normalised), and is seen in ${x.nWith} of ${
          x.nPool
        } samples. Reads crossing over from there could explain this level.`
    : `Not flagged at the current settings.`;
  return (
    head +
    `<br><span style="font-size:0.85em;color:#ccc">${lead}</span>` +
    `<div style="font-size:0.85em;color:#eee;margin-top:4px">${esc(r["Specimen ID"] || x.unit)}: ${_ttCtFmt(
      x.reads,
      0,
    )} reads · ${_ttCtFmt(x.rpm, 1)} RPM${x.breadth != null ? ` · ${_ttCtFmt(x.breadth, 1)}% breadth` : ""}${
      x.top
        ? `<br>Top other sample: ${esc(x.top)} · ${_ttCtFmt(x.topReads, 0)} reads · ${_ttCtFmt(x.topRpm, 1)} RPM`
        : ""
    }</div>` +
    _ttCtTermsTable(x) +
    `<span style="font-size:0.85em;color:#ccc">Score <b style="color:#fff">${x.score.toFixed(2)}</b> (flag at ≥ ${
      TT_CT.p.cutoff
    }). Flag only: TASS and the call are unchanged. Summary → Crossover shows this organism across every sample.</span>`
  );
}
function ttCtSourceTipHTML(r, src) {
  src = src || ttCtSourceFor(r);
  if (!src) return "";
  const esc = _ttCtEsc;
  const rows = src.recipients
    .map(
      (x) =>
        `<tr><td style="padding-right:10px">${esc(x.unit)}</td><td style="padding-right:10px">${_ttCtFmt(
          x.reads,
          0,
        )}</td><td style="padding-right:10px">${_ttCtFold(x.ratio)}</td><td>${x.score.toFixed(2)}</td></tr>`,
    )
    .join("");
  return (
    `<b style="color:#64b5f6"><i class="fas fa-arrow-right-from-bracket"></i> Possible crossover source</b><br>` +
    `<span style="font-size:0.85em;color:#ccc">${esc(r["Detected Organism"])} is strongest here (${_ttCtFmt(
      r["# Reads Aligned"],
      0,
    )} reads in ${esc(r["Specimen ID"])}). A low-level detection of it is flagged as possible crossover in ${
      src.n
    } other sample${src.n === 1 ? "" : "s"}${
      src.n > src.recipients.length ? " (top " + src.recipients.length + ")" : ""
    }.</span>` +
    `<table style="border-collapse:collapse;font-size:0.85em;margin:5px 0">` +
    `<tr style="color:#90caf9"><td style="padding-right:10px">Flagged in</td><td style="padding-right:10px">Reads</td><td style="padding-right:10px">vs here</td><td>Score</td></tr>${rows}</table>` +
    `<span style="font-size:0.85em;color:#ccc">Informational only: this detection, its TASS and its call are not affected.</span>`
  );
}
/** One line for other tooltips (heatmap cells). "" when not flagged or marking is off. */
function ttCtTipLineHTML(r) {
  const x = TT_CT.on ? ttCtFor(r) : null;
  if (!x || !x.suspect) return "";
  return x.kind === "background"
    ? `<br><span style="color:#b39ddb"><i class="fas fa-flask"></i> Possible background (also in ${_ttCtEsc(
        x.ncUnit,
      )}) · score ${x.score.toFixed(2)}</span>`
    : `<br><span style="color:#90caf9"><i class="fas fa-right-left"></i> Possible crossover: ${_ttCtFold(
        x.ratio,
      )} of ${_ttCtEsc(x.top)} · score ${x.score.toFixed(2)}</span>`;
}
function _ttCtFindRow(sample, level, key, name) {
  if (typeof DATA === "undefined" || !DATA) return null;
  const unit = _ttCtUnitOf(sample);
  let hit = null;
  for (let i = 0; i < DATA.length; i++) {
    const r = DATA[i];
    if (String(r.Level || "Strain") !== level) continue;
    if (key ? String(r["Taxonomic ID #"]) !== key : String(r["Detected Organism"]) !== name) continue;
    const s = String(r["Specimen ID"]);
    if (s === sample) return r;
    if (!hit && _ttCtUnitOf(s) === unit) hit = r;
  }
  return hit;
}
(function () {
  if (typeof document === "undefined") return;
  const target = (ev) => (ev.target && ev.target.closest ? ev.target.closest("[data-ct-tip]") : null);
  const html = (el) => {
    const d = el.dataset;
    const base = _ttCtFindRow(d.ctS, d.ctL, d.ctK, d.ctN);
    if (!base) return "";
    const r = Object.assign({}, base, { "Specimen ID": d.ctS });
    if (d.ctTip === "sus") return ttCtSuspectTipHTML(r, ttCtFor(base));
    if (d.ctTip === "src") return ttCtSourceTipHTML(r, ttCtSourceFor(base));
    return "";
  };
  const show = (ev) => {
    const el = target(ev);
    if (!el || typeof showTip !== "function") return;
    ev.stopPropagation();
    const h = html(el);
    if (!h) return;
    if (ev.clientX != null) showTip(h, ev);
    else {
      const b = el.getBoundingClientRect();
      showTip(h, { clientX: b.right, clientY: b.top });
    }
  };
  document.addEventListener("mouseover", show);
  document.addEventListener("focusin", show);
  document.addEventListener("mousemove", (ev) => {
    if (target(ev) && typeof moveTip === "function") moveTip(ev);
  });
  const hide = (ev) => {
    const el = target(ev);
    if (el && !el.contains(ev.relatedTarget) && typeof hideTip === "function") hideTip();
  };
  document.addEventListener("mouseout", hide);
  document.addEventListener("focusout", hide);
})();

/* ── heatmap: a corner (top-left) on flagged cells ───────────────────────── */
/** Called by drawHeatmap after the cells. lookup: "sample|||organism" → row. */
function ttCtHeatmapMarks(svg, organisms, samples, lookup, xScale, yScale) {
  const note = typeof document !== "undefined" ? document.getElementById("hm-ct-note") : null;
  if (!TT_CT.on) {
    if (note) note.hidden = true;
    return 0;
  }
  let n = 0;
  const g = svg.append("g").attr("class", "hm-ct").style("pointer-events", "none");
  organisms.forEach((org) => {
    samples.forEach((sp) => {
      const r = lookup[`${sp}|||${org}`];
      const x = r ? ttCtFor(r) : null;
      if (!x || !x.suspect) return;
      n++;
      const x0 = xScale(sp),
        y0 = yScale(org),
        s = Math.min(9, xScale.bandwidth() * 0.36);
      g.append("path")
        .attr("class", "hm-ct-corner" + (x.kind === "background" ? " bg" : ""))
        .attr("d", `M${x0},${y0} L${x0 + s},${y0} L${x0},${y0 + s} Z`);
    });
  });
  if (note) {
    note.hidden = false;
    note.innerHTML =
      `<span class="hm-ct-key" aria-hidden="true"></span> crossover? · <span class="hm-ct-key bg" aria-hidden="true"></span> background? ` +
      `— possible crossover (${n} cell${n === 1 ? "" : "s"} in view; Summary → Crossover for details)`;
  }
  return n;
}

/* ── saved settings ──────────────────────────────────────────────────────── */
function _ttCtSave() {
  try {
    localStorage.setItem(TT_CT_STORE_KEY, JSON.stringify({ on: TT_CT.on, p: TT_CT.p, scope: TT_CT.scope }));
  } catch (e) {
    /* storage blocked: the setting still holds for this page view */
  }
}
function _ttCtApplyParams(p) {
  if (!p || typeof p !== "object") return;
  Object.keys(TT_CT_DEFAULTS).forEach((k) => {
    const v = Number(p[k]);
    if (p[k] !== undefined && p[k] !== null && p[k] !== "" && isFinite(v) && v >= 0) TT_CT.p[k] = v;
  });
}
function _ttCtApplyConfig(cfg) {
  if (!cfg || typeof cfg !== "object") return;
  if (typeof cfg.on === "boolean") TT_CT.on = cfg.on;
  _ttCtApplyParams(cfg.p);
  if (typeof cfg.scope === "string") TT_CT.scope = cfg.scope;
  if (cfg.view === "all" || cfg.view === "only" || cfg.view === "hide") TT_CT.view = cfg.view;
  if (Array.isArray(cfg.excluded)) TT_CT.excluded = new Set(cfg.excluded.map(String));
}
function _ttCtLoad() {
  try {
    const raw = typeof localStorage !== "undefined" ? localStorage.getItem(TT_CT_STORE_KEY) : null;
    if (raw) _ttCtApplyConfig(JSON.parse(raw));
  } catch (e) {
    /* ignore */
  }
  TT_CT.view = "all";
  TT_CT.excluded = new Set();
}
/** Session export / import (29_session_state.js). */
function ttCtCaptureConfig() {
  return {
    on: TT_CT.on,
    p: Object.assign({}, TT_CT.p),
    scope: TT_CT.scope,
    view: TT_CT.view,
    excluded: Array.from(TT_CT.excluded),
  };
}
function ttCtLoadConfig(cfg) {
  if (!cfg) return; // a session saved before this check existed: keep the current settings
  _ttCtApplyConfig(cfg);
  ttCtInvalidate();
  ttCtSyncUi();
}

/* ── "Compare within" choices ────────────────────────────────────────────── */
function _ttCtScopeChoices() {
  let cands = [];
  try {
    if (typeof window !== "undefined" && window.metaGrouping && typeof window.metaGrouping.candidates === "function")
      cands = window.metaGrouping.candidates() || [];
  } catch (e) {
    cands = [];
  }
  return cands.map((c) => ({ key: c.key, label: c.label || c.key }));
}
function _ttCtCheckScope() {
  if (TT_CT.scope && !_ttCtScopeChoices().some((c) => c.key === TT_CT.scope)) {
    TT_CT.scope = "";
    ttCtInvalidate();
  }
}

/* ── applying a change ───────────────────────────────────────────────────── */
let _ttCtViewTimer = null;
let _ttCtReportTimer = null;
function _ttCtRedrawReport() {
  const run = () => {
    try {
      if (typeof redraw === "function") redraw();
    } catch (e) {
      if (typeof console !== "undefined" && console.warn) console.warn("[crossover] redraw skipped:", e);
    }
    ttCtSyncUi();
  };
  if (typeof ttBusyRun === "function") {
    try {
      const p = ttBusyRun("Crossover check…", run);
      if (p && p.catch) p.catch(() => {});
      return;
    } catch (e) {
      /* fall through */
    }
  }
  run();
}
/** Settings changed. The view redraws soon; the rest of the report too when it
 *  shows the result (marking on, or Only / Hide flagged). */
function _ttCtChanged(opts) {
  const o = Object.assign({ save: true, reportNow: false }, opts || {});
  ttCtInvalidate();
  if (o.save) _ttCtSave();
  clearTimeout(_ttCtViewTimer);
  _ttCtViewTimer = setTimeout(() => {
    ttCtRenderResults();
    ttCtSyncUi();
  }, 90);
  if (TT_CT.on || TT_CT.view !== "all" || o.reportNow) {
    clearTimeout(_ttCtReportTimer);
    _ttCtReportTimer = setTimeout(_ttCtRedrawReport, o.reportNow ? 0 : 450);
  }
}
function ttCtSet(patch, opts) {
  const before = TT_CT.on + "|" + TT_CT.view;
  _ttCtApplyConfig(patch || {});
  if (patch && patch.scope === "") TT_CT.scope = "";
  if (!TT_CT.on && patch && patch.on === false) TT_CT.view = "all"; // switching off never leaves rows filtered
  const after = TT_CT.on + "|" + TT_CT.view;
  _ttCtChanged(Object.assign({ reportNow: before !== after }, opts || {}));
  _ttCtSyncControls();
}
function ttCtToggle(on) {
  ttCtSet({ on: on === undefined ? !TT_CT.on : !!on });
}

/* ── banner pill + popover ───────────────────────────────────────────────── */
function _ttCtComparable() {
  if (typeof DATA === "undefined" || !DATA) return false;
  const seen = new Set();
  for (let i = 0; i < DATA.length && seen.size < 2; i++) {
    const s = DATA[i] && DATA[i]["Specimen ID"];
    if (s && !_ttCtIsInsilico(s)) seen.add(s);
  }
  return seen.size >= 2;
}
function _ttCtViewLevel() {
  const el = typeof document !== "undefined" ? document.getElementById("view-level") : null;
  return String((el && el.value) || "Strain").toLowerCase();
}
function _ttCtCounts() {
  const res = ttCtResult();
  const lvl = _ttCtViewLevel();
  const out = { level: lvl, flagged: 0, cross: 0, bg: 0, units: 0, compared: 0, orgs: 0 };
  if (!res) return out;
  res.units.forEach((U) => {
    if (U.ctl) return;
    out.compared++;
    const f = U.nFlag[lvl] || 0;
    if (f) out.units++;
    out.flagged += f;
    out.cross += U.nCross[lvl] || 0;
    out.bg += U.nBg[lvl] || 0;
  });
  res.orgs.forEach((O) => {
    if (O.level === lvl && O.nFlag) out.orgs++;
  });
  return out;
}
/** Banner text for the PDF cover; also keeps the pill current. */
function ttCtSyncUi() {
  if (typeof document === "undefined") return { text: "" };
  const pill = document.getElementById("banner-ct");
  const label = document.getElementById("banner-ct-label");
  const legend = document.getElementById("sum-ct-legend");
  const ok = _ttCtComparable();
  let text = "";
  let c = null;
  if (ok && TT_CT.on) {
    _ttCtCheckScope();
    c = _ttCtCounts();
    text =
      `${c.flagged} ${c.level}-level detection${c.flagged === 1 ? "" : "s"} flagged as possible crossover ` +
      `(${c.cross} crossover, ${c.bg} background; ${c.units} of ${c.compared} sample${c.compared === 1 ? "" : "s"})`;
  }
  if (legend) legend.hidden = !(ok && TT_CT.on);
  if (pill) {
    pill.hidden = !ok;
    pill.classList.toggle("on", !!(ok && TT_CT.on));
    pill.classList.toggle("has-flags", !!(c && c.flagged > 0));
    pill.classList.toggle("filtering", TT_CT.view !== "all");
    pill.setAttribute("aria-pressed", TT_CT.on ? "true" : "false");
    if (label)
      label.textContent = TT_CT.on
        ? `Crossover ${c ? c.flagged : 0}${
            TT_CT.view === "only" ? " · only" : TT_CT.view === "hide" ? " · hidden" : ""
          }`
        : "Crossover check";
    pill.setAttribute(
      "aria-label",
      TT_CT.on ? `Crossover check is on: ${text}. Click to turn off.` : "Crossover check is off. Click to turn on.",
    );
    if (!document.getElementById("banner-ct-pop")) _ttCtInitPopover(pill);
    _ttCtRenderPopover();
  }
  return { text };
}
function _ttCtRenderPopover() {
  const pop = document.getElementById("banner-ct-pop");
  if (!pop) return;
  const on = TT_CT.on;
  const c = on ? _ttCtCounts() : null;
  const seg = (attr, cur, items) =>
    items
      .map(
        ([v, t, icon]) =>
          `<button type="button" data-ct-act="${attr}" data-ct-v="${v}" class="${cur === v ? "active" : ""}">${
            icon ? `<i class="fas ${icon}"></i> ` : ""
          }${t}</button>`,
      )
      .join("");
  const preset = Object.keys(TT_CT_PRESETS).find((k) =>
    Object.keys(TT_CT_PRESETS[k]).every((p) => Number(TT_CT.p[p]) === TT_CT_PRESETS[k][p]),
  );
  let html =
    `<div class="bmp-title"><i class="fas fa-vials"></i> Crossover check ` +
    `<button type="button" class="ct-switch${on ? " on" : ""}" data-ct-act="toggle" role="switch" aria-checked="${
      on ? "true" : "false"
    }"><span></span>${on ? "On" : "Off"}</button></div>` +
    `<div class="bmp-note">Compares every organism across all samples of the report: low-level detections of an organism that is far more abundant in another sample (crossover), or that is also in the negative control (background), are flagged for a second look. ` +
    `Computed here from the report itself; TASS, calls and exports are unchanged.</div>`;
  if (on && c)
    html += `<div class="bmp-note"><b style="color:#eceff1">${c.flagged}</b> ${c.level}-level detection${
      c.flagged === 1 ? "" : "s"
    } flagged (${c.cross} crossover · ${c.bg} background) in <b style="color:#eceff1">${c.units}</b> of ${
      c.compared
    } samples, across ${c.orgs} organism${c.orgs === 1 ? "" : "s"}.</div>`;
  html +=
    `<div class="bhp-controls"><span class="ct-pop-lbl">Sensitivity</span>` +
    seg("preset", preset || "", [
      ["sensitive", "Sensitive"],
      ["balanced", "Balanced"],
      ["strict", "Strict"],
    ]) +
    (preset ? "" : `<span class="bmp-note" style="margin:0">custom</span>`) +
    `</div>`;
  if (on)
    html +=
      `<div class="bhp-controls"><span class="ct-pop-lbl">In the report</span>` +
      seg("view", TT_CT.view, [
        ["all", "Show all", "fa-list"],
        ["only", "Only flagged", "fa-filter"],
        ["hide", "Hide flagged", "fa-eye-slash"],
      ]) +
      `</div>`;
  html +=
    `<div class="bhp-controls"><button type="button" data-ct-act="open"><i class="fas fa-chart-simple"></i> Open the analysis (weights, cutoffs, per organism)</button></div>` +
    `<div class="bmp-note">Click the pill to switch marking on or off. A flag is a prompt to look closer, not a call: a shared true positive (an outbreak strain, a common commensal) can look the same in read counts alone.</div>`;
  pop.innerHTML = html;
}
let _ttCtPopPinned = false;
let _ttCtPopHideTimer = null;
function _ttCtPopSet(show, pin) {
  const pill = document.getElementById("banner-ct");
  const pop = document.getElementById("banner-ct-pop");
  if (!pill || !pop) return;
  if (pin !== undefined) _ttCtPopPinned = pin;
  const visible = show || _ttCtPopPinned;
  clearTimeout(_ttCtPopHideTimer);
  pop.classList.toggle("show", visible);
  pill.classList.toggle("open", visible);
  if (visible) {
    const r = pill.getBoundingClientRect();
    const w = pop.offsetWidth;
    pop.style.top = `${Math.round(r.bottom + 8)}px`;
    pop.style.left = `${Math.round(Math.max(8, Math.min(r.left, window.innerWidth - w - 8)))}px`;
  }
}
function _ttCtInitPopover(pill) {
  const pop = document.createElement("div");
  pop.id = "banner-ct-pop";
  pop.className = "banner-mask-pop banner-hop-pop banner-ct-pop";
  pop.setAttribute("role", "dialog");
  pop.setAttribute("aria-label", "Crossover check");
  document.body.appendChild(pop);
  pill.setAttribute("aria-describedby", "banner-ct-pop");
  const hideSoon = () => {
    clearTimeout(_ttCtPopHideTimer);
    _ttCtPopHideTimer = setTimeout(() => {
      if (!pop.contains(document.activeElement)) _ttCtPopSet(false);
    }, 180);
  };
  pill.addEventListener("mouseenter", () => _ttCtPopSet(true));
  pill.addEventListener("mouseleave", hideSoon);
  pill.addEventListener("focus", () => _ttCtPopSet(true));
  pill.addEventListener("blur", hideSoon);
  pop.addEventListener("mouseenter", () => clearTimeout(_ttCtPopHideTimer));
  pop.addEventListener("mouseleave", hideSoon);
  pop.addEventListener("focusout", hideSoon);
  pill.addEventListener("click", (ev) => {
    ev.stopPropagation();
    ttCtToggle();
    _ttCtPopSet(true);
  });
  pill.addEventListener("keydown", (ev) => {
    if (ev.key === "Enter" || ev.key === " ") {
      ev.preventDefault();
      ttCtToggle();
      _ttCtPopSet(true);
    }
    if (ev.key === "Escape") _ttCtPopSet(false, false);
  });
  pop.addEventListener("keydown", (ev) => {
    if (ev.key === "Escape") {
      _ttCtPopSet(false, false);
      pill.focus();
    }
  });
  pop.addEventListener("click", (ev) => {
    ev.stopPropagation();
    const b = ev.target.closest ? ev.target.closest("[data-ct-act]") : null;
    if (!b || b.disabled) return;
    const act = b.dataset.ctAct;
    if (act === "toggle") ttCtToggle();
    else if (act === "preset") ttCtSet({ p: TT_CT_PRESETS[b.dataset.ctV] });
    else if (act === "view") ttCtSet({ view: b.dataset.ctV }, { save: false });
    else if (act === "open") {
      _ttCtPopSet(false, false);
      if (typeof ttOpenSummaryInner === "function") ttOpenSummaryInner("crossover");
    }
  });
  document.addEventListener("click", () => _ttCtPopSet(false, false));
  window.addEventListener("resize", () => {
    if (pop.classList.contains("show")) _ttCtPopSet(true);
  });
  return pop;
}

/* ═══════════════════════════════════════════════════════════════════════════
   The analysis view: Summary → Crossover
   ═══════════════════════════════════════════════════════════════════════════ */
const _TT_CT_WEIGHTS = [
  [
    "wAbund",
    "Abundance vs top sample",
    "How far below the strongest other sample this detection sits (depth-normalised).",
  ],
  ["wPrev", "Run prevalence", "How many of the other samples carry the organism."],
  ["wCov", "Thin coverage", "Little of the genome covered here (low breadth)."],
  ["wNeg", "In negative control", "The organism is also in a negative control at a comparable level."],
];
const _TT_CT_COLORS = {
  crossover: "#1565c0", // flagged: from another sample
  background: "#5e35b1", // flagged: also in the negative control
  clean: "#b0bec5", // not flagged
  source: "#0b2545", // the strongest sample
  pos: "#2e7d32",
  neg: "#c62828",
};

function _ttCtBuildSkeleton(root) {
  const num = (id, label, tip, min, max, step) =>
    `<label class="ct-num" title="${_ttCtEsc(
      tip,
    )}">${label}<input type="number" id="${id}" min="${min}" max="${max}" step="${step}"></label>`;
  root.innerHTML = `
  <div class="chart-wrap ct-wrap">
    <div class="chart-title"><i class="fas fa-vials"></i> Crossover check — every organism across the run</div>
    <div class="xs-note">Each organism is compared across <b>all</b> samples of the report (below-cutoff detections included, sidebar filters not applied):
      a low-level detection of something far more abundant in another sample scores as possible <b style="color:${
        _TT_CT_COLORS.crossover
      }">crossover</b>; one that is also in the negative control as possible <b style="color:${
        _TT_CT_COLORS.background
      }">background</b>. Adjust the weights and cutoffs and the plot, the flags and (when marking is on) every tab follow. Flag only: TASS and calls never change.</div>
    <div class="ct-controls">
      <div class="ct-row">
        <span class="ct-lbl">Sensitivity</span>
        <div class="ct-seg" id="ct-preset" role="group" aria-label="Sensitivity preset">
          <button type="button" data-v="sensitive" title="More flags: lower score cutoff, 1:100 counts fully">Sensitive</button>
          <button type="button" data-v="balanced" title="The defaults">Balanced</button>
          <button type="button" data-v="strict" title="Fewer flags: higher cutoff, needs a 1:10,000 gap for full weight">Strict</button>
        </div>
        <label class="ct-slider" title="A detection is flagged when its crossover score reaches this">Flag at score ≥
          <input type="range" id="ct-cutoff" min="0.2" max="0.95" step="0.01"><b id="ct-cutoff-v"></b></label>
        <span class="ct-spacer"></span>
        <label class="ct-marktoggle" title="Show crossover? / background? markers in the Detections table, the Table tab and the heatmap"><input type="checkbox" id="ct-on"> Mark flagged across the report</label>
        <div class="ct-seg" id="ct-view" role="group" aria-label="Flagged detections in the report">
          <button type="button" data-v="all"><i class="fas fa-list"></i> Show all</button>
          <button type="button" data-v="only"><i class="fas fa-filter"></i> Only flagged</button>
          <button type="button" data-v="hide"><i class="fas fa-eye-slash"></i> Hide flagged</button>
        </div>
        <button type="button" class="ct-btn" id="ct-tsv" title="Every detection's score and terms at this level (TSV)"><i class="fas fa-download"></i> Scores</button>
      </div>
      <details class="ct-adv" id="ct-adv" open>
        <summary>Weights &amp; thresholds <button type="button" class="ct-link" id="ct-reset">reset to defaults</button></summary>
        <div class="ct-grid">
          <div class="ct-box">
            <div class="ct-box-h">Evidence weights <span class="ct-dim">(relative, 0–100)</span></div>
            ${_TT_CT_WEIGHTS
              .map(
                ([k, l, t]) =>
                  `<label class="ct-w" title="${_ttCtEsc(
                    t,
                  )}"><span class="ct-w-l">${l}</span><input type="range" id="ct-${k}" data-p="${k}" min="0" max="100" step="5"><b id="ct-${k}-v"></b></label>`,
              )
              .join("")}
            <div class="ct-dim" id="ct-w-note"></div>
          </div>
          <div class="ct-box">
            <div class="ct-box-h">Cutoffs</div>
            <label class="ct-num" title="The abundance term reaches its full value when a detection is this many times below the strongest other sample">Full weight at 1:
              <select id="ct-fullFold"><option value="100">100</option><option value="1000">1,000</option><option value="10000">10,000</option><option value="100000">100,000</option></select></label>
            ${num(
              "ct-minFold",
              "Crossover: top sample ≥ N× higher",
              "A crossover flag needs the strongest other sample to be at least this many times higher (depth-normalised)",
              1,
              1e6,
              1,
            )}
            ${num(
              "ct-minSourceReads",
              "Source needs ≥ reads",
              "A sample only counts as a source with at least this many reads of the organism",
              0,
              1e9,
              10,
            )}
            ${num(
              "ct-minSamples",
              "Organism in ≥ samples",
              "A crossover flag needs the organism in at least this many samples (source included)",
              2,
              10000,
              1,
            )}
            ${num(
              "ct-prevFull",
              "Prevalence full weight at %",
              "The prevalence term reaches its full value when this share of the other samples carry the organism",
              1,
              100,
              5,
            )}
            ${num(
              "ct-breadthRef",
              "Well covered at breadth %",
              "Breadth at or above this counts as well covered (thin-coverage term = 0)",
              0.1,
              100,
              1,
            )}
            <label class="ct-num" title="Only compare samples that share this metadata value (e.g. the same run, plate or batch). Negative controls count everywhere.">Compare within
              <select id="ct-scope"><option value="">whole report</option></select></label>
          </div>
        </div>
      </details>
      <div class="ct-row">
        <input type="search" id="ct-q" placeholder="Find organism…" aria-label="Find organism">
        <label>Sort <select id="ct-sort">
          <option value="flags">Most flagged</option><option value="score">Highest score</option>
          <option value="prev">Prevalence</option><option value="top">Top abundance</option><option value="name">Name</option>
        </select></label>
        <label title="List only organisms with at least one flagged detection"><input type="checkbox" id="ct-only"> Only organisms with flags</label>
        <label>Show <select id="ct-topn"><option value="15">15</option><option value="30">30</option><option value="60">60</option><option value="120">120</option><option value="0">all</option></select></label>
        <span class="ct-spacer"></span>
        <span class="ct-dim" id="ct-level-note"></span>
      </div>
    </div>
    <div id="ct-summary" class="ct-summary"></div>
    <div id="ct-legend" class="xs-note ct-legend"></div>
    <div id="ct-plot" class="ct-plot"></div>
    <details id="ct-samples-wrap" class="ct-samples-wrap" open>
      <summary>Per-sample crossover burden <span class="ct-dim">— untick a sample to leave it out of the comparison (e.g. a known spiked or failed sample); click a name to highlight it in the plot</span></summary>
      <div id="ct-samples"></div>
    </details>
  </div>`;
}

function _ttCtSyncControls() {
  const root = typeof document !== "undefined" ? document.getElementById("ct-root") : null;
  if (!root || !root.dataset.built) return;
  const $ = (id) => document.getElementById(id);
  const p = TT_CT.p;
  const act = document.activeElement;
  const setv = (el, v) => {
    if (el && el !== act) el.value = v;
  };
  setv($("ct-cutoff"), p.cutoff);
  if ($("ct-cutoff-v")) $("ct-cutoff-v").textContent = Number(p.cutoff).toFixed(2);
  const res = ttCtResult();
  const W = (res && res.weights) || {};
  const avail = { wAbund: true, wPrev: true, wCov: !res || res.hasBreadth, wNeg: !res || res.hasNeg };
  const pct = (names) => {
    const t = names.reduce((s, n) => s + (W[n] || 0), 0);
    return names.map((n) => (t > 0 ? Math.round((100 * (W[n] || 0)) / t) : 0) + "%");
  };
  _TT_CT_WEIGHTS.forEach(([k]) => {
    const el = $("ct-" + k);
    setv(el, p[k]);
    if (el) el.disabled = !avail[k];
    const lab = $("ct-" + k + "-v");
    if (lab) lab.textContent = !avail[k] ? "n/a" : String(Math.round(Number(p[k]) || 0));
    const row = el && el.closest(".ct-w");
    if (row) row.classList.toggle("ct-na", !avail[k]);
  });
  const cx = pct(["ab", "prev", "cov"]),
    bx = pct(["neg", "prev", "cov"]);
  const notes = [
    `Crossover score = abundance ${cx[0]} · prevalence ${cx[1]} · coverage ${cx[2]}` +
      (res && res.hasNeg
        ? `; background score = negative control ${bx[0]} · prevalence ${bx[1]} · coverage ${bx[2]}`
        : ""),
  ];
  if (res && !res.hasNeg) notes.push("no negative control in this report: the control term is off");
  if (res && !res.hasBreadth) notes.push("no breadth / coverage column: the coverage term is off");
  if ($("ct-w-note")) $("ct-w-note").textContent = notes.join(" · ");
  const ff = $("ct-fullFold");
  if (ff && ff !== act) {
    if (!Array.from(ff.options).some((o) => Number(o.value) === Number(p.fullFold))) {
      const o = document.createElement("option");
      o.value = p.fullFold;
      o.textContent = _ttCtFmt(p.fullFold, 0);
      ff.appendChild(o);
    }
    ff.value = String(p.fullFold);
  }
  ["minFold", "minSourceReads", "minSamples", "prevFull", "breadthRef"].forEach((k) => setv($("ct-" + k), p[k]));
  const sc = $("ct-scope");
  if (sc && sc !== act) {
    _ttCtCheckScope();
    const ch = _ttCtScopeChoices();
    sc.innerHTML =
      `<option value="">whole report</option>` +
      ch
        .map(
          (c) =>
            `<option value="${_ttCtEsc(c.key)}"${c.key === TT_CT.scope ? " selected" : ""}>${_ttCtEsc(
              c.label,
            )}</option>`,
        )
        .join("");
    sc.disabled = !ch.length;
    sc.title = ch.length ? "" : "No groupable metadata column in this report";
  }
  const preset = Object.keys(TT_CT_PRESETS).find((k) =>
    Object.keys(TT_CT_PRESETS[k]).every((q) => Number(p[q]) === TT_CT_PRESETS[k][q]),
  );
  root.querySelectorAll("#ct-preset button").forEach((b) => b.classList.toggle("active", b.dataset.v === preset));
  root.querySelectorAll("#ct-view button").forEach((b) => b.classList.toggle("active", b.dataset.v === TT_CT.view));
  if ($("ct-on")) $("ct-on").checked = TT_CT.on;
  const u = TT_CT.ui;
  setv($("ct-q"), u.q);
  setv($("ct-sort"), u.sort);
  if ($("ct-only")) $("ct-only").checked = u.only;
  setv($("ct-topn"), String(u.topN));
}

function _ttCtWire(root) {
  const $ = (id) => document.getElementById(id);
  const setP = (k, v) => {
    TT_CT.p[k] = v;
    _ttCtSyncControls();
    _ttCtChanged();
  };
  $("ct-cutoff").addEventListener("input", (e) => setP("cutoff", Number(e.target.value)));
  _TT_CT_WEIGHTS.forEach(([k]) => $("ct-" + k).addEventListener("input", (e) => setP(k, Number(e.target.value))));
  $("ct-fullFold").addEventListener("change", (e) => setP("fullFold", Number(e.target.value)));
  ["minFold", "minSourceReads", "minSamples", "prevFull", "breadthRef"].forEach((k) =>
    $("ct-" + k).addEventListener("change", (e) => {
      const v = Number(e.target.value);
      const lo = Number(e.target.min),
        hi = Number(e.target.max);
      if (!isFinite(v) || v < lo || v > hi) {
        e.target.value = TT_CT.p[k];
        return;
      }
      setP(k, k === "minSamples" ? Math.round(v) : v);
    }),
  );
  $("ct-scope").addEventListener("change", (e) => {
    TT_CT.scope = e.target.value;
    _ttCtChanged();
  });
  root.querySelector("#ct-preset").addEventListener("click", (e) => {
    const b = e.target.closest("button");
    if (b) ttCtSet({ p: TT_CT_PRESETS[b.dataset.v] });
  });
  root.querySelector("#ct-view").addEventListener("click", (e) => {
    const b = e.target.closest("button");
    if (!b) return;
    // Filtering the report only makes sense with the markers on.
    ttCtSet(b.dataset.v === "all" ? { view: "all" } : { view: b.dataset.v, on: true }, {
      save: b.dataset.v === "all" ? false : true,
    });
  });
  $("ct-on").addEventListener("change", (e) => ttCtSet({ on: e.target.checked }));
  $("ct-reset").addEventListener("click", (e) => {
    e.preventDefault();
    e.stopPropagation();
    TT_CT.p = Object.assign({}, TT_CT_DEFAULTS);
    TT_CT.scope = "";
    _ttCtSyncControls();
    _ttCtChanged();
  });
  $("ct-tsv").addEventListener("click", () => ttCtDownloadTsv());
  let qTimer = null;
  $("ct-q").addEventListener("input", (e) => {
    TT_CT.ui.q = e.target.value;
    clearTimeout(qTimer);
    qTimer = setTimeout(ttCtRenderResults, 150);
  });
  $("ct-sort").addEventListener("change", (e) => {
    TT_CT.ui.sort = e.target.value;
    ttCtRenderResults();
  });
  $("ct-only").addEventListener("change", (e) => {
    TT_CT.ui.only = e.target.checked;
    ttCtRenderResults();
  });
  $("ct-topn").addEventListener("change", (e) => {
    TT_CT.ui.topN = Number(e.target.value) || 0;
    ttCtRenderResults();
  });
  $("ct-samples").addEventListener("change", (e) => {
    const cb = e.target.closest("input[data-ct-unit]");
    if (!cb) return;
    const u = cb.dataset.ctUnit;
    if (cb.checked) TT_CT.excluded.delete(u);
    else TT_CT.excluded.add(u);
    _ttCtChanged({ save: false });
  });
  $("ct-samples").addEventListener("click", (e) => {
    const b = e.target.closest("[data-ct-focus]");
    if (!b) return;
    const u = b.dataset.ctFocus;
    TT_CT.ui.focus = TT_CT.ui.focus === u ? null : u;
    ttCtRenderResults();
  });
  $("ct-summary").addEventListener("click", (e) => {
    const b = e.target.closest("[data-ct-clear-focus]");
    if (!b) return;
    TT_CT.ui.focus = null;
    ttCtRenderResults();
  });
}

/** Draw (or refresh) the view. Safe to call when its pane is not in the page. */
function ttCtRender() {
  const root = typeof document !== "undefined" ? document.getElementById("ct-root") : null;
  if (!root) return;
  if (!root.dataset.built) {
    _ttCtBuildSkeleton(root);
    _ttCtWire(root);
    root.dataset.built = "1";
  }
  _ttCtSyncControls();
  ttCtRenderResults();
}

function _ttCtOrgRows(res) {
  const lvl = _ttCtViewLevel();
  const u = TT_CT.ui;
  let rx = null;
  if (u.q && u.q.trim()) {
    try {
      rx = new RegExp(u.q.trim(), "i");
    } catch (e) {
      const q = u.q.trim().toLowerCase();
      rx = { test: (s) => String(s).toLowerCase().includes(q) };
    }
  }
  let rows = [];
  res.orgs.forEach((O) => {
    if (O.level !== lvl) return;
    if (u.only && !O.nFlag) return;
    if (!u.only && O.recs.filter((x) => x.role !== "neg").length < 2 && !O.ncRpm) return; // nothing to compare
    if (rx && !rx.test(O.name) && !rx.test(O.tax)) return;
    if (u.focus && !O.recs.some((x) => x.unit === u.focus)) return;
    rows.push(O);
  });
  const by = {
    flags: (a, b) => b.nFlag - a.nFlag || b.maxScore - a.maxScore || b.nWith - a.nWith,
    score: (a, b) => b.maxScore - a.maxScore || b.nFlag - a.nFlag,
    prev: (a, b) => b.nWith / Math.max(1, b.nPool) - a.nWith / Math.max(1, a.nPool) || b.nFlag - a.nFlag,
    top: (a, b) => b.topRpm - a.topRpm,
    name: (a, b) => String(a.name).localeCompare(String(b.name)),
  };
  rows.sort(by[u.sort] || by.flags);
  const total = rows.length;
  if (u.topN > 0) rows = rows.slice(0, u.topN);
  return { rows, total, lvl };
}

function ttCtRenderResults() {
  if (typeof document === "undefined") return;
  const root = document.getElementById("ct-root");
  if (!root || !root.dataset.built) return;
  const sumEl = document.getElementById("ct-summary");
  const legEl = document.getElementById("ct-legend");
  const plot = document.getElementById("ct-plot");
  const sampEl = document.getElementById("ct-samples");
  if (typeof hideTip === "function") hideTip();
  const res = _ttCtComparable() ? ttCtResult() : null;
  if (!res) {
    sumEl.innerHTML = `<span class="ct-dim">Needs detections from at least two samples.</span>`;
    legEl.innerHTML = "";
    plot.innerHTML = "";
    sampEl.innerHTML = "";
    return;
  }
  const c = _ttCtCounts();
  const { rows, total, lvl } = _ttCtOrgRows(res);
  const lvlNote = document.getElementById("ct-level-note");
  if (lvlNote)
    lvlNote.textContent = `Level: ${lvl} (follows View level) · abundance from ${
      res.norm === "input"
        ? "reads per million input reads"
        : res.norm === "rpm"
        ? "the RPM column"
        : "share of aligned reads (no input read counts in this report)"
    }`;
  const esc = _ttCtEsc;
  const focus = TT_CT.ui.focus;
  sumEl.innerHTML =
    `<span class="ct-kpi"><b>${c.flagged}</b> flagged detection${c.flagged === 1 ? "" : "s"}</span>` +
    `<span class="ct-kpi" style="color:${_TT_CT_COLORS.crossover}"><b>${c.cross}</b> crossover</span>` +
    `<span class="ct-kpi" style="color:${_TT_CT_COLORS.background}"><b>${c.bg}</b> background</span>` +
    `<span class="ct-kpi"><b>${c.units}</b> / ${c.compared} samples affected</span>` +
    `<span class="ct-kpi"><b>${c.orgs}</b> organism${c.orgs === 1 ? "" : "s"}</span>` +
    (TT_CT.excluded.size ? `<span class="ct-kpi ct-dim">${TT_CT.excluded.size} left out of the pool</span>` : "") +
    (TT_CT.scope ? `<span class="ct-kpi ct-dim">compared within ${esc(TT_CT.scope)}</span>` : "") +
    (focus
      ? `<span class="ct-kpi ct-focus">highlighting <b>${esc(
          focus,
        )}</b> <button type="button" class="ct-link" data-ct-clear-focus="1">clear</button></span>`
      : "") +
    (TT_CT.on
      ? ""
      : `<span class="ct-kpi ct-dim">marking off — tick “Mark flagged across the report” to see these in every tab</span>`);
  const dot = (col, shape, hollow) =>
    shape === "d"
      ? `<span class="ct-lg-d" style="background:${col}"></span>`
      : `<span class="ct-lg-dot" style="${
          hollow ? `border:1.5px solid ${col};background:#fff` : `background:${col}`
        }"></span>`;
  legEl.innerHTML =
    `<b>Prevalence</b>: samples carrying the organism, <span style="color:${_TT_CT_COLORS.crossover}">flagged</span> part first. ` +
    `<b>Abundance vs top sample</b>: one dot per sample, depth-normalised and placed relative to the strongest other sample (right edge = top); ` +
    `shaded = at least ${_ttCtFmt(
      TT_CT.p.minFold,
      0,
    )}× below (crossover range), dotted line = full weight (1:${_ttCtFmt(TT_CT.p.fullFold, 0)}). ` +
    `${dot(_TT_CT_COLORS.source)} top source ${dot(_TT_CT_COLORS.crossover)} crossover? ${dot(
      _TT_CT_COLORS.background,
    )} background? ${dot(_TT_CT_COLORS.clean)} not flagged ${dot(
      "#78909c",
      "",
      true,
    )} hollow = breadth below ${_ttCtFmt(TT_CT.p.breadthRef, 1)}% ` +
    `${dot(_TT_CT_COLORS.neg, "d")} negative / ${dot(_TT_CT_COLORS.pos, "d")} positive control. ` +
    `<b>Crossover score</b>: violin of the samples' scores (the crossover or background pattern, whichever fits better), red line = flag cutoff. Hover anything for the numbers.`;
  if (!rows.length) {
    plot.innerHTML = `<div class="ct-empty">${
      TT_CT.ui.only
        ? "No organism has a flagged detection at these settings. Untick “Only organisms with flags” to see every organism, or try Sensitive."
        : "No organism matches."
    }</div>`;
  } else {
    _ttCtDrawPlot(plot, rows, res);
    if (total > rows.length) {
      const more = document.createElement("div");
      more.className = "ct-dim";
      more.style.margin = "4px 0 0 4px";
      more.textContent = `Showing ${rows.length} of ${total} organisms (Show → all for every one).`;
      plot.appendChild(more);
    }
  }
  _ttCtRenderSamples(sampEl, res, lvl);
}

function _ttCtDrawPlot(wrap, rows, res) {
  wrap.innerHTML = "";
  if (typeof d3 === "undefined") {
    wrap.innerHTML = `<div class="ct-empty">The plotting library (d3) did not load.</div>`;
    return;
  }
  const P = TT_CT.p;
  const focus = TT_CT.ui.focus;
  const Wv = Math.max(860, wrap.clientWidth || 1000);
  const rowH = 34,
    mT = 30,
    mB = 8,
    mL = 230,
    prevW = 96,
    gap = 26,
    scW = 170,
    verdW = 200;
  const abW = Math.max(260, Wv - (mL + prevW + gap + gap + scW + gap + verdW + 10));
  const W = mL + prevW + gap + abW + gap + scW + gap + verdW + 10;
  const H = rows.length * rowH + mT + mB;
  const prevX0 = mL,
    abX0 = prevX0 + prevW + gap,
    scX0 = abX0 + abW + gap,
    verdX0 = scX0 + scW + gap;
  const D = Math.min(7, Math.max(4, Math.ceil(Math.log10(Math.max(10, P.fullFold))) + 2));
  const xa = (r) => abX0 + (1 + Math.max(-D, Math.min(0, Math.log10(Math.max(r, 1e-12))) / 1) / D) * abW;
  const xs = (v) => scX0 + Math.max(0, Math.min(1, v)) * scW;
  const svg = d3
    .select(wrap)
    .append("svg")
    .attr("class", "ct-svg")
    .attr("width", W)
    .attr("height", H)
    .attr("viewBox", `0 0 ${W} ${H}`);
  const hdr = (cx, t) =>
    svg
      .append("text")
      .attr("x", cx)
      .attr("y", 11)
      .attr("text-anchor", "middle")
      .style("font-size", "10px")
      .style("font-weight", "600")
      .style("fill", "#455a64")
      .text(t);
  hdr(prevX0 + prevW / 2, "Prevalence");
  hdr(abX0 + abW / 2, "Abundance vs top sample (depth-normalised)");
  hdr(scX0 + scW / 2, "Crossover score");
  hdr(verdX0 + verdW / 2, "What we see");
  // abundance axis
  const plotTop = mT,
    plotBot = H - mB;
  const xMin = xa(1 / P.minFold);
  svg
    .append("rect")
    .attr("x", abX0)
    .attr("y", plotTop - 4)
    .attr("width", Math.max(0, xMin - abX0))
    .attr("height", plotBot - plotTop + 4)
    .attr("fill", "#e3f2fd");
  for (let d = 0; d <= D; d++) {
    const x = xa(Math.pow(10, -d));
    svg
      .append("line")
      .attr("x1", x)
      .attr("x2", x)
      .attr("y1", plotTop - 4)
      .attr("y2", plotBot)
      .attr("stroke", "#eceff1");
    svg
      .append("text")
      .attr("x", x)
      .attr("y", 24)
      .attr("text-anchor", d === 0 ? "end" : d === D ? "start" : "middle")
      .style("font-size", "9px")
      .style("fill", "#90a4ae")
      .text(d === 0 ? "top" : d === D ? "≤ " + _ttCtFold(Math.pow(10, -d)) : _ttCtFold(Math.pow(10, -d)));
  }
  svg
    .append("line")
    .attr("x1", xMin)
    .attr("x2", xMin)
    .attr("y1", plotTop - 4)
    .attr("y2", plotBot)
    .attr("stroke", _TT_CT_COLORS.crossover)
    .attr("stroke-dasharray", "4,3");
  const xF = xa(1 / P.fullFold);
  svg
    .append("line")
    .attr("x1", xF)
    .attr("x2", xF)
    .attr("y1", plotTop - 4)
    .attr("y2", plotBot)
    .attr("stroke", _TT_CT_COLORS.crossover)
    .attr("stroke-dasharray", "1,3");
  // score axis
  [0, 0.25, 0.5, 0.75, 1].forEach((v) => {
    svg
      .append("line")
      .attr("x1", xs(v))
      .attr("x2", xs(v))
      .attr("y1", plotTop - 4)
      .attr("y2", plotBot)
      .attr("stroke", "#eceff1");
    svg
      .append("text")
      .attr("x", xs(v))
      .attr("y", 24)
      .attr("text-anchor", "middle")
      .style("font-size", "9px")
      .style("fill", "#90a4ae")
      .text(v);
  });
  svg
    .append("line")
    .attr("x1", xs(P.cutoff))
    .attr("x2", xs(P.cutoff))
    .attr("y1", plotTop - 4)
    .attr("y2", plotBot)
    .attr("stroke", "#c62828")
    .attr("stroke-dasharray", "4,3");

  const tip = (html) => (ev) => {
    if (typeof showTip === "function") showTip(html(), ev);
  };
  const move = (ev) => typeof moveTip === "function" && moveTip(ev);
  const hide = () => typeof hideTip === "function" && hideTip();

  rows.forEach((O, i) => {
    const y0 = mT + i * rowH;
    const cy = y0 + rowH / 2;
    const g = svg.append("g").attr("class", "ct-rowg");
    if (i % 2 === 0)
      g.append("rect")
        .attr("x", 0)
        .attr("y", y0)
        .attr("width", mL - 6)
        .attr("height", rowH)
        .attr("fill", "#fafbfc");
    // name
    const nm = String(O.name || O.tax);
    const lab = nm.length > 34 ? nm.slice(0, 33) + "…" : nm;
    g.append("text")
      .attr("x", mL - 10)
      .attr("y", cy + 4)
      .attr("text-anchor", "end")
      .style("font-size", "11px")
      .style("font-style", "italic")
      .style("fill", O.nFlag ? "#4e342e" : "#37474f")
      .style("font-weight", O.nFlag ? 600 : 400)
      .style("cursor", "default")
      .text(lab)
      .on(
        "mouseover",
        tip(() => _ttCtOrgTipHTML(O)),
      )
      .on("mousemove", move)
      .on("mouseout", hide);
    // prevalence bar
    const pool = Math.max(1, O.nPool);
    const fr = Math.min(1, O.nWith / pool),
      ff = Math.min(fr, O.nFlag / pool);
    g.append("rect")
      .attr("x", prevX0)
      .attr("y", cy - 6)
      .attr("width", prevW)
      .attr("height", 12)
      .attr("rx", 2)
      .attr("fill", "#eceff1");
    g.append("rect")
      .attr("x", prevX0 + ff * prevW)
      .attr("y", cy - 6)
      .attr("width", Math.max(0, (fr - ff) * prevW))
      .attr("height", 12)
      .attr("fill", "#bbdefb");
    g.append("rect")
      .attr("x", prevX0)
      .attr("y", cy - 6)
      .attr("width", ff * prevW)
      .attr("height", 12)
      .attr("fill", _TT_CT_COLORS.crossover)
      .attr("opacity", 0.85);
    g.append("text")
      .attr("x", prevX0 + prevW - 3)
      .attr("y", cy + 4)
      .attr("text-anchor", "end")
      .style("font-size", "9px")
      .style("fill", "#263238")
      .text(`${O.nWith}/${O.nPool}`);
    g.append("rect")
      .attr("x", prevX0)
      .attr("y", cy - 8)
      .attr("width", prevW)
      .attr("height", 16)
      .attr("fill", "transparent")
      .on(
        "mouseover",
        tip(() => _ttCtOrgTipHTML(O)),
      )
      .on("mousemove", move)
      .on("mouseout", hide);
    // abundance dots
    const sorted = O.recs.slice().sort((a, b) => (a.suspect ? 1 : 0) - (b.suspect ? 1 : 0));
    sorted.forEach((x) => {
      const isTop = x.role !== "neg" && O.top === x.unit;
      const r = isTop ? 1 : x.ratio == null ? 1 : x.ratio;
      const px = xa(r);
      const jit = (_ttCtHash(x.unit) - 0.5) * (rowH - 12);
      const py = cy + (isTop ? 0 : jit);
      const col =
        x.role === "neg"
          ? _TT_CT_COLORS.neg
          : x.role === "pos"
          ? _TT_CT_COLORS.pos
          : isTop
          ? _TT_CT_COLORS.source
          : x.suspect
          ? _TT_CT_COLORS[x.kind]
          : _TT_CT_COLORS.clean;
      const thin = x.breadth != null && x.breadth < P.breadthRef;
      const faded = focus && x.unit !== focus;
      const isFocus = focus && x.unit === focus;
      let el;
      if (x.role === "neg" || x.role === "pos") {
        const s = 4.5;
        el = g
          .append("path")
          .attr("d", `M${px},${py - s} L${px + s},${py} L${px},${py + s} L${px - s},${py} Z`)
          .attr("fill", col);
      } else {
        el = g
          .append("circle")
          .attr("cx", px)
          .attr("cy", py)
          .attr("r", isTop ? 5 : x.suspect ? 4 : 3.2)
          .attr("fill", thin && !isTop ? "#fff" : col)
          .attr("stroke", col)
          .attr("stroke-width", thin ? 1.5 : 0.6)
          .attr("fill-opacity", x.suspect || isTop ? 0.95 : 0.7);
      }
      el.attr("opacity", faded ? 0.18 : 1)
        .style("cursor", "help")
        .on(
          "mouseover",
          tip(() => _ttCtDotTipHTML(x, O)),
        )
        .on("mousemove", move)
        .on("mouseout", hide);
      if (isFocus)
        g.append("circle")
          .attr("cx", px)
          .attr("cy", py)
          .attr("r", 7.5)
          .attr("fill", "none")
          .attr("stroke", "#212121")
          .attr("stroke-width", 1.5)
          .style("pointer-events", "none");
    });
    // score violin + dots
    const scores = O.recs.filter((x) => x.role === "sample" && O.top !== x.unit).map((x) => x.score);
    if (scores.length >= 3) {
      const bw = 0.07,
        steps = 40,
        pts = [];
      let mx = 0;
      for (let s = 0; s <= steps; s++) {
        const v = s / steps;
        let k = 0;
        scores.forEach((q) => (k += Math.exp(-0.5 * ((v - q) / bw) ** 2)));
        pts.push([v, k]);
        if (k > mx) mx = k;
      }
      const h = rowH / 2 - 3;
      const area = d3
        .area()
        .x((d) => xs(d[0]))
        .y0((d) => cy + (d[1] / mx) * h)
        .y1((d) => cy - (d[1] / mx) * h)
        .curve(d3.curveBasis);
      g.append("path").attr("d", area(pts)).attr("fill", "#cfd8dc").attr("opacity", 0.75);
    }
    O.recs.forEach((x) => {
      if (x.role !== "sample" || O.top === x.unit) return;
      const faded = focus && x.unit !== focus;
      g.append("circle")
        .attr("cx", xs(x.score))
        .attr("cy", cy + (_ttCtHash(x.unit + "s") - 0.5) * (rowH - 14))
        .attr("r", x.suspect ? 3.2 : 2.4)
        .attr("fill", x.suspect ? _TT_CT_COLORS[x.kind] : "#607d8b")
        .attr("opacity", faded ? 0.15 : 0.9)
        .style("cursor", "help")
        .on(
          "mouseover",
          tip(() => _ttCtDotTipHTML(x, O)),
        )
        .on("mousemove", move)
        .on("mouseout", hide);
    });
    // verdict
    const v1 = O.nCross
      ? `${O.nCross} likely crossover`
      : O.nBg
      ? `${O.nBg} background (in NC)`
      : O.nWith > 1
      ? "no flags"
      : O.ncRpm
      ? "only in the negative control"
      : "single sample";
    const v2 = O.nCross
      ? `from ${O.top}` + (O.nBg ? ` · ${O.nBg} background` : "")
      : O.top
      ? `top: ${O.top} (${_ttCtFmt(O.topReads, 0)} reads)`
      : "";
    g.append("text")
      .attr("x", verdX0)
      .attr("y", cy - 2)
      .style("font-size", "10.5px")
      .style("font-weight", O.nFlag ? 700 : 400)
      .style("fill", O.nCross ? _TT_CT_COLORS.crossover : O.nBg ? _TT_CT_COLORS.background : "#607d8b")
      .text(v1);
    g.append("text")
      .attr("x", verdX0)
      .attr("y", cy + 10)
      .style("font-size", "9.5px")
      .style("fill", "#78909c")
      .text(v2.length > 38 ? v2.slice(0, 37) + "…" : v2);
    g.append("line")
      .attr("x1", 0)
      .attr("x2", W)
      .attr("y1", y0 + rowH)
      .attr("y2", y0 + rowH)
      .attr("stroke", "#f1f3f5");
  });
}

function _ttCtOrgTipHTML(O) {
  const esc = _ttCtEsc;
  return (
    `<b>${esc(O.name)}</b> <span style="color:#90a4ae">(${esc(O.level)}${
      O.tax.indexOf("name:") === 0 ? "" : " · taxid " + esc(O.tax)
    })</span><br>` +
    `<span style="font-size:0.85em;color:#ddd">In ${O.nWith} of ${O.nPool} samples · strongest in <b>${esc(
      O.top || "—",
    )}</b> (${_ttCtFmt(O.topReads, 0)} reads, ${_ttCtFmt(O.topRpm, 1)} RPM)` +
    (O.ncRpm ? `<br>Also in negative control ${esc(O.ncUnit)} (${_ttCtFmt(O.ncRpm, 1)} RPM)` : "") +
    `<br>Flagged: <b style="color:#90caf9">${O.nCross}</b> crossover · <b style="color:#b39ddb">${
      O.nBg
    }</b> background · highest score ${O.maxScore.toFixed(2)}</span>`
  );
}
function _ttCtDotTipHTML(x, O) {
  const esc = _ttCtEsc;
  if (x.role === "neg")
    return (
      `<b style="color:#ef9a9a">◆ Negative control ${esc(
        x.unit,
      )}</b><br><span style="font-size:0.85em;color:#ddd">${esc(O.name)}: ${_ttCtFmt(x.reads, 0)} reads · ${_ttCtFmt(
        x.rpm,
        1,
      )} RPM` +
      (x.ratio != null ? ` · ${_ttCtFold(x.ratio)} of ${esc(x.top)}` : "") +
      `<br>Used as evidence for background (reagent / kit) signal; never flagged itself.</span>`
    );
  if (x.role === "pos")
    return `<b style="color:#a5d6a7">◆ Positive control ${esc(
      x.unit,
    )}</b><br><span style="font-size:0.85em;color:#ddd">${esc(O.name)}: ${_ttCtFmt(x.reads, 0)} reads · ${_ttCtFmt(
      x.rpm,
      1,
    )} RPM<br>Can be a crossover source; never flagged itself.</span>`;
  if (O.top === x.unit)
    return (
      `<b style="color:#90caf9">● Top sample ${esc(x.unit)}</b><br><span style="font-size:0.85em;color:#ddd">${esc(
        O.name,
      )}: ${_ttCtFmt(x.reads, 0)} reads · ${_ttCtFmt(x.rpm, 1)} RPM${
        x.breadth != null ? ` · ${_ttCtFmt(x.breadth, 1)}% breadth` : ""
      }<br>The strongest sample for this organism, so the reference the others are compared with.` +
      (O.nCross ? ` Possible source of ${O.nCross} flagged crossover detection${O.nCross === 1 ? "" : "s"}.` : "") +
      `</span>`
    );
  return ttCtSuspectTipHTML({ "Specimen ID": x.unit, "Detected Organism": O.name }, x);
}

function _ttCtAllUnits() {
  const seen = new Map();
  (typeof DATA !== "undefined" && DATA ? DATA : []).forEach((r) => {
    const s = r && r["Specimen ID"];
    if (!s || _ttCtIsInsilico(s)) return;
    const u = _ttCtUnitOf(s);
    if (!seen.has(u)) seen.set(u, s);
  });
  return seen;
}
function _ttCtRenderSamples(el, res, lvl) {
  const esc = _ttCtEsc;
  const all = _ttCtAllUnits();
  const rows = [];
  all.forEach((s, u) => {
    const U = res.units.get(u);
    const ctl = U
      ? U.ctl
      : String(_ttCtCtlType(s) || "")
      ? _ttCtIsNegType(_ttCtCtlType(s))
        ? "negative"
        : "positive"
      : "";
    const n = U ? U.nDet[lvl] || 0 : 0;
    const f = U ? U.nFlag[lvl] || 0 : 0;
    rows.push({
      u,
      ctl,
      inPool: !TT_CT.excluded.has(u),
      n,
      f,
      cross: U ? U.nCross[lvl] || 0 : 0,
      bg: U ? U.nBg[lvl] || 0 : 0,
      src: U ? U.nSrc[lvl] || 0 : 0,
      share: n ? f / n : 0,
    });
  });
  rows.sort((a, b) => b.inPool - a.inPool || b.share - a.share || b.f - a.f || (a.u < b.u ? -1 : 1));
  const focus = TT_CT.ui.focus;
  el.innerHTML =
    `<div class="ct-samples-scroll"><table class="ct-samples"><thead><tr>` +
    `<th title="In the comparison pool">Use</th><th>Sample</th><th class="num" title="Detections at this level">Detections</th>` +
    `<th class="num">Flagged</th><th style="min-width:150px">Share flagged</th><th class="num" title="crossover / background">Crossover · bg</th>` +
    `<th class="num" title="Detections here that are the source of a flagged crossover elsewhere">Source of</th></tr></thead><tbody>` +
    rows
      .map(
        (r) =>
          `<tr class="${r.inPool ? "" : "ct-out"}${focus === r.u ? " ct-focus-row" : ""}">` +
          `<td><input type="checkbox" data-ct-unit="${esc(r.u)}" ${r.inPool ? "checked" : ""} aria-label="Use ${esc(
            r.u,
          )} in the comparison"></td>` +
          `<td><button type="button" class="ct-link" data-ct-focus="${esc(r.u)}">${esc(r.u)}</button>${
            r.ctl
              ? ` <span class="ct-ctl ct-ctl-${r.ctl === "negative" ? "neg" : "pos"}">${
                  r.ctl === "negative" ? "NC" : "PC"
                }</span>`
              : ""
          }</td>` +
          `<td class="num">${r.inPool ? r.n : "—"}</td><td class="num">${r.inPool && !r.ctl ? r.f : "—"}</td>` +
          `<td>${
            r.inPool && !r.ctl
              ? `<span class="ct-sbar"><span class="x" style="width:${Math.round(
                  (100 * r.cross) / Math.max(1, r.n),
                )}%"></span><span class="b" style="width:${Math.round(
                  (100 * r.bg) / Math.max(1, r.n),
                )}%"></span></span> ${Math.round(100 * r.share)}%`
              : ""
          }</td>` +
          `<td class="num">${r.inPool && !r.ctl ? `${r.cross} · ${r.bg}` : ""}</td><td class="num">${
            r.inPool && r.src ? r.src : ""
          }</td></tr>`,
      )
      .join("") +
    `</tbody></table></div>`;
}

function ttCtDownloadTsv() {
  const res = ttCtResult();
  if (!res) return;
  const lvl = _ttCtViewLevel();
  const head = [
    "Sample",
    "Level",
    "Taxonomic ID #",
    "Detected Organism",
    "Reads",
    "RPM",
    "Breadth %",
    "Top other sample",
    "Top reads",
    "Ratio to top",
    "Samples with organism",
    "Samples in pool",
    "Negative control RPM",
    "Term abundance",
    "Term prevalence",
    "Term coverage",
    "Term negative control",
    "Score",
    "Flagged",
    "Kind",
  ];
  const cell = (v) => String(v == null ? "" : v).replace(/[\t\r\n]+/g, " ");
  const f = (v, d) => (v == null || !isFinite(v) ? "" : Number(v).toFixed(d));
  const lines = [head.join("\t")];
  res.recs.forEach((x) => {
    if (x.level !== lvl || x.role !== "sample") return;
    const t = x.terms || {};
    lines.push(
      [
        x.unit,
        x.level,
        x.key.split(_TT_CT_SEP)[1].indexOf("name:") === 0 ? "" : x.key.split(_TT_CT_SEP)[1],
        x.name,
        x.reads,
        f(x.rpm, 3),
        f(x.breadth, 2),
        x.top,
        x.topReads,
        x.ratio == null ? "" : x.ratio.toExponential(3),
        x.nWith,
        x.nPool,
        f(x.ncRpm, 3),
        f(t.ab, 3),
        f(t.prev, 3),
        f(t.cov, 3),
        f(t.neg, 3),
        f(x.score, 3),
        x.suspect ? "true" : "false",
        x.kind || "",
      ]
        .map(cell)
        .join("\t"),
    );
  });
  const text =
    `# TaxTriage crossover check · ${new Date().toISOString()} · settings ${JSON.stringify(TT_CT.p)}` +
    (TT_CT.scope ? ` · within ${TT_CT.scope}` : "") +
    (TT_CT.excluded.size ? ` · excluded ${Array.from(TT_CT.excluded).join(",")}` : "") +
    "\n" +
    lines.join("\n") +
    "\n";
  if (typeof _downloadText === "function") _downloadText(text, "crossover_scores.tsv", "text/tab-separated-values");
}

/** After detections are added / removed (uploads, clear). The caller redraws. */
function ttCtOnDataChanged() {
  // Drop pool exclusions for samples no longer loaded.
  const all = _ttCtAllUnits();
  Array.from(TT_CT.excluded).forEach((u) => all.has(u) || TT_CT.excluded.delete(u));
  ttCtInvalidate();
  ttCtSyncUi();
}

/* ── start-up ────────────────────────────────────────────────────────────── */
_ttCtLoad();
(function () {
  if (typeof document === "undefined") return;
  const start = () => {
    ttCtSyncUi();
    // Marking was left on: draw once more so the markers are there.
    if (TT_CT.on && _ttCtComparable()) {
      try {
        if (typeof redraw === "function") redraw();
      } catch (e) {
        /* the report is not built yet; its first draw picks the markers up */
      }
      ttCtSyncUi();
    }
  };
  if (document.readyState === "complete") setTimeout(start, 0);
  else window.addEventListener("load", () => setTimeout(start, 0), { once: true });
})();
