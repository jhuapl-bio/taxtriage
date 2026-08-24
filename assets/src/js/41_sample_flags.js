/* ═══════════════════════════════════════════════════════════════════════════
       -  §  SAMPLE FLAGS  (whole-sample QC rules)
       -     Everything else in the report filters DETECTIONS (rows). This
       -     section filters SAMPLES: a small, generic rule engine that answers
       -     "does this whole sample look usable?" — enough reads, enough
       -     distinct organisms above a TASS cutoff, the right metadata — and
       -     marks the ones that fail.
       -
       -     Model:  one rule = { field on some source } { operator } { value },
       -     with an action of "flag" (mark it, keep it visible) or "hide"
       -     (also drop it from every chart/table). Rules combine with ANY / ALL.
       -
       -     Sources a rule can read:
       -       meta     SAMPLE_META  — total_reads, aligned_reads, platform, …
       -       derived  computed from DATA — unique organisms above TASS,
       -                detections passing threshold, max TASS, …
       -       runmeta  RUN_META / uploaded metadata CSV — any column
       -       data     any numeric detection column, aggregated per sample
       -                (max / min / mean / sum / count > 0)
       -
       -     Defaults ship from the pipeline: make_report.py bakes the
       -     params.report_flag_* settings into BOOT.sample_flags, which seeds
       -     TT_FLAG_RULES here. Users edit them live in the Sample QC dialog
       -     (see 42_sample_flags_ui.js).
       -
       -     "hide" is implemented by driving the EXISTING sampleHidden map
       -     (like the sidebar search auto-hide does) rather than adding a
       -     second gate to filteredData() — so every tab, export and
       -     aggregation honours it with no further changes.
═══════════════════════════════════════════════════════════════════════════ */

/* ── Rule storage ────────────────────────────────────────────────────────── */
let TT_FLAG_RULES = []; // [{ id, on, source, field, op, value, agg, tass, action }]
let TT_FLAG_LOGIC = "any"; // "any" = flag if ANY rule trips, "all" = only if every rule trips
let TT_FLAG_ENABLED = true; // master switch for the whole engine
let TT_FLAG_MISSING_FAILS = false; // treat a missing/blank value as a trip
let TT_FLAG_HIDE_ALL = false; // sidebar "Hide flagged" — promotes every flag to a hide
let _TT_FLAG_SEQ = 0; // rule id counter
const _flagAutoHidden = new Set(); // ids WE hid (mirrors _filterAutoHidden's contract)

/* ── Operators ───────────────────────────────────────────────────────────── */
/*  num:  true  → numeric only        false → text only        null → both
    novalue: operator takes no right-hand value (empty / not empty).        */
const TT_FLAG_OPS = [
  { op: "<", label: "is less than", num: true },
  { op: "<=", label: "is at most", num: true },
  { op: ">", label: "is greater than", num: true },
  { op: ">=", label: "is at least", num: true },
  { op: "==", label: "equals", num: null },
  { op: "!=", label: "does not equal", num: null },
  { op: "contains", label: "contains", num: false },
  { op: "!contains", label: "does not contain", num: false },
  { op: "regex", label: "matches regex", num: false },
  { op: "empty", label: "is empty / missing", num: null, novalue: true },
  { op: "!empty", label: "has any value", num: null, novalue: true },
];

/* ── Field catalogues ────────────────────────────────────────────────────── */
/*  SAMPLE_META keys the pipeline always writes. Anything else present in
    SAMPLE_META is discovered at runtime and appended (see ttFlagFields).   */
const TT_FLAG_META_FIELDS = [
  { key: "total_reads", label: "Total reads (sample)", num: true },
  { key: "aligned_reads", label: "Aligned reads", num: true },
  { key: "total_organism_reads", label: "Organism-assigned reads", num: true },
  { key: "num_species_groups", label: "Species groups (pipeline)", num: true },
  { key: "num_subkeys", label: "Species keys (pipeline)", num: true },
  { key: "num_keys", label: "Strain keys (pipeline)", num: true },
  { key: "num_toplevelkeys", label: "Genus keys (pipeline)", num: true },
  { key: "sample_type", label: "Sample type", num: false },
  { key: "platform", label: "Sequencing platform", num: false },
  { key: "control_type", label: "Control type", num: false },
  { key: "specimen", label: "Specimen group", num: false },
];

/*  Counts computed from DATA. `needsTass` fields take an extra TASS cutoff
    so "> 3 organisms above TASS 75" is one rule rather than two.          */
const TT_FLAG_DERIVED_FIELDS = [
  { key: "unique_taxids_above_tass", label: "Distinct organisms above TASS", num: true, needsTass: true },
  { key: "unique_taxids", label: "Distinct organisms (any TASS)", num: true },
  { key: "passing_detections", label: "Detections passing threshold", num: true },
  { key: "detections", label: "Detection rows", num: true },
  { key: "high_consequence", label: "High-consequence organisms", num: true },
  { key: "unique_genera", label: "Distinct genera", num: true },
  { key: "max_tass", label: "Highest TASS score", num: true },
  { key: "reads_aligned_sum", label: "Sum of # Reads Aligned", num: true },
];

/** The sources a rule may name. Anything else is rejected on load. */
const _TT_FLAG_SOURCE_KEYS = new Set(["meta", "derived", "runmeta", "data"]);

const TT_FLAG_AGGS = [
  { key: "max", label: "max" },
  { key: "min", label: "min" },
  { key: "mean", label: "mean" },
  { key: "sum", label: "sum" },
  { key: "count", label: "count > 0" },
];

/*  Keys in SAMPLE_META that are structural rather than measurements — never
    worth offering as a rule field. */
const _TT_FLAG_META_SKIP = new Set([
  "sample_name",
  "weights",
  "best_cutoffs",
  "best_cutoffs_by_domain",
  "commit_id",
  "workflow_revision",
]);

/** Every field a rule can be built on, grouped by source. Recomputed on demand
 *  (cheap) so metadata uploaded after load shows up in the dialog. */
function ttFlagFields() {
  const meta = TT_FLAG_META_FIELDS.slice();
  const seen = new Set(meta.map((f) => f.key));
  const sm = (typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {};
  Object.keys(sm).forEach((s) => {
    const rec = sm[s] || {};
    Object.keys(rec).forEach((k) => {
      if (seen.has(k) || _TT_FLAG_META_SKIP.has(k)) return;
      const v = rec[k];
      if (v != null && typeof v === "object" && !Array.isArray(v)) return; // nested payloads
      seen.add(k);
      meta.push({ key: k, label: _ttFlagPrettyKey(k), num: _ttFlagLooksNumeric(sm, k) });
    });
  });

  //  Metadata columns: RUN_META is the primary source, but _ttFlagValueFor
  //  falls back to SAMPLE_META, so descriptive sample fields belong here too —
  //  a run whose samplesheet columns all landed in SAMPLE_META would otherwise
  //  show an empty "Metadata column" list.
  const runmeta = [];
  const rseen = new Set();
  ((typeof RUN_META !== "undefined" && RUN_META) || []).forEach((rec) => {
    Object.keys(rec || {}).forEach((k) => {
      if (k === "sample_name" || rseen.has(k)) return;
      rseen.add(k);
      runmeta.push({ key: k, label: _ttFlagPrettyKey(k), num: false });
    });
  });
  ["sample_type", "platform", "control_type", "specimen"].forEach((k) => {
    if (rseen.has(k)) return;
    const present = Object.keys(sm).some((s) => {
      const v = (sm[s] || {})[k];
      return v != null && v !== "";
    });
    if (!present) return;
    rseen.add(k);
    runmeta.push({ key: k, label: _ttFlagPrettyKey(k), num: false });
  });
  runmeta.sort((a, b) => a.label.localeCompare(b.label));

  //  Numeric detection columns only — aggregating a text column per sample
  //  has no meaning, and the metadata sources above already cover text.
  const cols = ((typeof ALL_COLS !== "undefined" && ALL_COLS) || []).filter(
    (c) => typeof NUMERIC !== "undefined" && NUMERIC && NUMERIC.has(c),
  );
  const data = cols.map((c) => ({ key: c, label: c, num: true }));

  return { meta, derived: TT_FLAG_DERIVED_FIELDS.slice(), runmeta, data };
}

function _ttFlagPrettyKey(k) {
  if (typeof _metaKeyLabel === "function") return _metaKeyLabel(k);
  return String(k)
    .replace(/_/g, " ")
    .replace(/\b\w/g, (c) => c.toUpperCase());
}

/** Does this SAMPLE_META key hold numbers in most samples? Decides which
 *  operators the dialog offers; the comparison itself always re-checks. */
function _ttFlagLooksNumeric(sm, key) {
  let n = 0,
    numeric = 0;
  for (const s in sm) {
    const v = (sm[s] || {})[key];
    if (v == null || v === "") continue;
    n++;
    if (typeof v === "number" || (typeof v === "string" && v.trim() !== "" && !isNaN(Number(v)))) numeric++;
    if (n >= 12) break;
  }
  return n > 0 && numeric === n;
}

/** Look up one field descriptor (for labels + numeric hints). */
function ttFlagFieldDef(rule) {
  if (!rule) return null;
  const all = ttFlagFields();
  const list = all[rule.source] || [];
  return list.find((f) => f.key === rule.field) || { key: rule.field, label: rule.field, num: true };
}

/* ── Rule helpers ────────────────────────────────────────────────────────── */
function ttFlagNewRule(seed) {
  const r = Object.assign(
    {
      id: "r" + ++_TT_FLAG_SEQ,
      on: true,
      source: "meta",
      field: "total_reads",
      op: "<",
      value: "10000",
      agg: "max",
      tass: 75,
      action: "flag",
    },
    seed || {},
  );
  //  Ids must be unique even when a seed carries one (pipeline defaults do).
  if (!r.id || TT_FLAG_RULES.some((x) => x.id === r.id)) r.id = "r" + ++_TT_FLAG_SEQ;
  return r;
}

/** Human-readable rendering of one rule — used by badges, tooltips and the
 *  flagged-sample list so the wording is identical everywhere. */
function ttFlagRuleLabel(rule) {
  if (!rule) return "";
  const def = ttFlagFieldDef(rule);
  const opDef = TT_FLAG_OPS.find((o) => o.op === rule.op) || { label: rule.op };
  let name = def.label;
  if (rule.source === "data") name = `${rule.agg || "max"}(${name})`;
  if (rule.source === "derived" && def.needsTass) name = `${name} ${_ttFlagNum(rule.tass)}`;
  return opDef.novalue ? `${name} ${opDef.label}` : `${name} ${opDef.label} ${_ttFlagFmt(rule.value)}`;
}

function _ttFlagNum(v) {
  const n = Number(v);
  return isFinite(n) ? n : 0;
}
function _ttFlagFmt(v) {
  const n = Number(v);
  if (v !== "" && v != null && isFinite(n)) return n.toLocaleString();
  return String(v == null ? "" : v);
}
function _ttFlagEsc(v) {
  return String(v == null ? "" : v)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

/* ── Per-sample inputs ───────────────────────────────────────────────────── */
/*  Strain-level rows only. DATA also carries Species / Genus roll-up rows
    (real ones from the JSON hierarchy, synthesized ones from
    _synthesizeHierarchy) — counting those would inflate every "distinct
    organism" figure two- or three-fold.                                    */
let _TTF_ROWS_IDX = { key: null, map: null };
function _ttFlagRowsBySample() {
  const key =
    (typeof DATA !== "undefined" && DATA ? DATA.length : 0) +
    "|" +
    (typeof SAMPLE_META_EPOCH !== "undefined" ? SAMPLE_META_EPOCH : 0);
  if (_TTF_ROWS_IDX.key === key && _TTF_ROWS_IDX.map) return _TTF_ROWS_IDX.map;
  const m = new Map();
  const src = (typeof DATA !== "undefined" && DATA) || [];
  for (let i = 0; i < src.length; i++) {
    const r = src[i];
    if ((r["Level"] || "Strain") !== "Strain") continue;
    const s = String(r["Specimen ID"] || "");
    if (!s) continue;
    let arr = m.get(s);
    if (!arr) m.set(s, (arr = []));
    arr.push(r);
  }
  _TTF_ROWS_IDX = { key, map: m };
  return m;
}

let _TTF_RUNMETA_IDX = { n: -1, map: null };
function _ttFlagRunMetaBySample() {
  const list = (typeof RUN_META !== "undefined" && RUN_META) || [];
  if (_TTF_RUNMETA_IDX.map && _TTF_RUNMETA_IDX.n === list.length) return _TTF_RUNMETA_IDX.map;
  const m = new Map();
  list.forEach((rec) => {
    if (rec && rec.sample_name != null) m.set(String(rec.sample_name), rec);
  });
  _TTF_RUNMETA_IDX = { n: list.length, map: m };
  return m;
}

/** Every sample the report knows about, whether or not it has detections. */
function _ttFlagAllSamples() {
  if (typeof _allSampleIds === "function") {
    try {
      return Array.from(_allSampleIds()).filter(Boolean);
    } catch (e) {
      /* fall through */
    }
  }
  const s = new Set();
  Object.keys((typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {}).forEach((k) => s.add(k));
  _ttFlagRowsBySample().forEach((_, k) => s.add(k));
  return Array.from(s).filter(Boolean);
}

/* ── Value resolution ────────────────────────────────────────────────────── */
function _ttFlagDerivedValue(field, rows, rule) {
  const n = typeof num === "function" ? num : (v) => (isNaN(parseFloat(v)) ? 0 : parseFloat(v));
  const truthy = typeof isTruthy === "function" ? isTruthy : (v) => !!v;
  const orgKey = (r) => String(r["Taxonomic ID #"] || r["Detected Organism"] || "");
  switch (field) {
    case "detections":
      return rows.length;
    case "unique_taxids":
      return new Set(rows.map(orgKey).filter(Boolean)).size;
    case "unique_taxids_above_tass": {
      const thr = _ttFlagNum(rule && rule.tass);
      const s = new Set();
      rows.forEach((r) => {
        if (n(r["TASS Score"]) >= thr) {
          const k = orgKey(r);
          if (k) s.add(k);
        }
      });
      return s.size;
    }
    case "passing_detections":
      return rows.filter((r) => truthy(r["Passes Threshold"])).length;
    case "high_consequence":
      return new Set(
        rows
          .filter((r) => truthy(r["High Consequence"]))
          .map(orgKey)
          .filter(Boolean),
      ).size;
    case "unique_genera":
      return new Set(rows.map((r) => r["Genus"]).filter(Boolean)).size;
    case "max_tass":
      return rows.length ? _ttMax(rows.map((r) => n(r["TASS Score"]))) : 0;
    case "reads_aligned_sum":
      return rows.reduce((s, r) => s + n(r["# Reads Aligned"]), 0);
  }
  return null;
}

function _ttFlagAggValue(rows, col, agg) {
  const vals = [];
  for (let i = 0; i < rows.length; i++) {
    const v = parseFloat(rows[i][col]);
    if (!isNaN(v)) vals.push(v);
  }
  if (!vals.length) return null; // no data → "missing", handled by the caller
  switch (agg) {
    case "min":
      return _ttMin(vals);
    case "mean":
      return vals.reduce((a, b) => a + b, 0) / vals.length;
    case "sum":
      return vals.reduce((a, b) => a + b, 0);
    case "count":
      return vals.filter((v) => v > 0).length;
    default:
      return _ttMax(vals);
  }
}

/** The value a rule should be tested against for one sample.
 *  Returns undefined when the sample carries nothing for that field — the
 *  caller decides whether "missing" counts as a trip. */
function _ttFlagValueFor(rule, sample, rows) {
  switch (rule.source) {
    case "meta": {
      const rec = ((typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {})[sample] || {};
      return rec[rule.field];
    }
    case "runmeta": {
      //  RUN_META first (that is what the Metadata & Mapping tab edits and
      //  what an uploaded CSV populates), then SAMPLE_META as a fallback —
      //  a samplesheet column can land in either depending on the run.
      const rec = _ttFlagRunMetaBySample().get(sample) || {};
      let v = rec[rule.field];
      if (_ttFlagIsMissing(v)) {
        const sm = ((typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {})[sample] || {};
        v = sm[rule.field];
      }
      return Array.isArray(v) ? v.join(", ") : v;
    }
    case "derived":
      return _ttFlagDerivedValue(rule.field, rows, rule);
    case "data": {
      const v = _ttFlagAggValue(rows, rule.field, rule.agg || "max");
      return v == null ? undefined : v;
    }
  }
  return undefined;
}

function _ttFlagIsMissing(v) {
  if (v === undefined || v === null) return true;
  if (typeof v === "string" && v.trim() === "") return true;
  if (typeof v === "number" && !isFinite(v)) return true;
  return false;
}

/** Does `actual` satisfy the rule's comparison? Missing values never satisfy a
 *  positive comparison — TT_FLAG_MISSING_FAILS decides whether they trip. */
function _ttFlagCompare(actual, rule) {
  const op = rule.op;
  if (op === "empty") return _ttFlagIsMissing(actual);
  if (op === "!empty") return !_ttFlagIsMissing(actual);
  if (_ttFlagIsMissing(actual)) return !!TT_FLAG_MISSING_FAILS;

  const raw = rule.value == null ? "" : String(rule.value);
  if (op === "<" || op === "<=" || op === ">" || op === ">=") {
    const a = Number(actual),
      b = Number(raw);
    if (!isFinite(a) || !isFinite(b)) return !!TT_FLAG_MISSING_FAILS;
    if (op === "<") return a < b;
    if (op === "<=") return a <= b;
    if (op === ">") return a > b;
    return a >= b;
  }
  if (op === "==" || op === "!=") {
    //  Numeric on both sides → numeric equality (so "5" == 5.0). Otherwise a
    //  case-insensitive string match, which is what metadata comparison wants.
    const a = Number(actual),
      b = Number(raw);
    let eq;
    if (isFinite(a) && isFinite(b) && String(actual).trim() !== "") eq = a === b;
    else eq = String(actual).trim().toLowerCase() === raw.trim().toLowerCase();
    return op === "==" ? eq : !eq;
  }
  if (op === "contains" || op === "!contains") {
    const hit = String(actual).toLowerCase().includes(raw.trim().toLowerCase());
    return op === "contains" ? hit : !hit;
  }
  if (op === "regex") {
    try {
      return new RegExp(raw, "i").test(String(actual));
    } catch (e) {
      return false; // invalid pattern mid-typing — never trips
    }
  }
  return false;
}

/* ── Evaluation ──────────────────────────────────────────────────────────── */
/*  Memoized on a fingerprint of the rules plus every input value they read.
    Hashing the referenced values (rather than just DATA.length) is what makes
    a metadata cell edit re-flag the sample without any explicit invalidation
    hook — the fields in play are a handful, so the walk is trivial.        */
let _TTF_EVAL = { key: null, value: null };

function ttFlagInvalidate() {
  _TTF_EVAL = { key: null, value: null };
  _TTF_ROWS_IDX = { key: null, map: null };
  _TTF_RUNMETA_IDX = { n: -1, map: null };
}

function _ttFlagSignature(samples) {
  const parts = [
    TT_FLAG_ENABLED ? 1 : 0,
    TT_FLAG_LOGIC,
    TT_FLAG_MISSING_FAILS ? 1 : 0,
    TT_FLAG_HIDE_ALL ? 1 : 0,
    (typeof DATA !== "undefined" && DATA ? DATA.length : 0) + "",
    samples.length + "",
  ];
  TT_FLAG_RULES.forEach((r) => {
    parts.push([r.on ? 1 : 0, r.source, r.field, r.op, r.value, r.agg, r.tass, r.action].join("~"));
  });
  //  Only metadata sources can change without DATA.length changing.
  const watched = TT_FLAG_RULES.filter((r) => r.on && (r.source === "meta" || r.source === "runmeta"));
  if (watched.length) {
    const sm = (typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {};
    const rm = _ttFlagRunMetaBySample();
    samples.forEach((s) => {
      watched.forEach((r) => {
        const rec = r.source === "meta" ? sm[s] || {} : rm.get(s) || {};
        const v = rec[r.field];
        parts.push(Array.isArray(v) ? v.join(",") : v == null ? "" : String(v));
      });
    });
  }
  return parts.join("");
}

/** Map(sample → { flagged, hide, hits: [{rule, actual, text}] }).
 *  Every known sample gets an entry, so callers can look up unconditionally. */
function ttFlagEvaluate() {
  const samples = _ttFlagAllSamples();
  const key = _ttFlagSignature(samples);
  if (_TTF_EVAL.key === key && _TTF_EVAL.value) return _TTF_EVAL.value;

  const active = TT_FLAG_ENABLED ? TT_FLAG_RULES.filter((r) => r.on) : [];
  const rowsBy = _ttFlagRowsBySample();
  const out = new Map();

  samples.forEach((s) => {
    const rows = rowsBy.get(s) || [];
    const hits = [];
    let tripped = 0;
    active.forEach((rule) => {
      const actual = _ttFlagValueFor(rule, s, rows);
      if (!_ttFlagCompare(actual, rule)) return;
      tripped++;
      hits.push({
        rule,
        actual,
        text: `${ttFlagRuleLabel(rule)} — actual: ${_ttFlagIsMissing(actual) ? "missing" : _ttFlagFmt(actual)}`,
      });
    });
    //  ANY: one trip is enough. ALL: every enabled rule must trip (a sample
    //  with no rules at all is never flagged).
    const flagged = active.length > 0 && (TT_FLAG_LOGIC === "all" ? tripped === active.length : tripped > 0);
    const hide = flagged && (TT_FLAG_HIDE_ALL || hits.some((h) => h.rule.action === "hide"));
    out.set(s, { flagged, hide, hits: flagged ? hits : [] });
  });

  _TTF_EVAL = { key, value: out };
  return out;
}

/* ── Lookups used by the tabs ────────────────────────────────────────────── */
/** Flag state for a sample id, a specimen/group name, or a detection row.
 *  With specimen merge on, a row's "Specimen ID" is a GROUP name that has no
 *  entry of its own — the group inherits the flags of its member libraries. */
function ttFlagStateFor(item) {
  const st = ttFlagEvaluate();
  let name = "";
  if (item && typeof item === "object") name = String(item["Specimen ID"] || "");
  else name = String(item == null ? "" : item);
  if (!name) return null;
  const direct = st.get(name);
  if (direct) return direct;
  if (typeof specimenGroups === "function") {
    const members = specimenGroups().get(name);
    if (members && members.length) {
      const hits = [];
      let flagged = false,
        hide = false;
      members.forEach((m) => {
        const s = st.get(m);
        if (!s || !s.flagged) return;
        flagged = true;
        hide = hide || s.hide;
        s.hits.forEach((h) => hits.push(Object.assign({ sample: m }, h)));
      });
      return { flagged, hide, hits, merged: true, members: members.slice() };
    }
  }
  return null;
}

function ttFlagIsFlagged(item) {
  const s = ttFlagStateFor(item);
  return !!(s && s.flagged);
}

/** Ids of every flagged sample, in sidebar order where one exists. */
function ttFlagFlaggedSamples() {
  const st = ttFlagEvaluate();
  const out = [];
  st.forEach((v, k) => {
    if (v.flagged) out.push(k);
  });
  return typeof _orderedSamples === "function" ? _orderedSamples(out) : out.sort();
}

function ttFlagCounts() {
  const st = ttFlagEvaluate();
  let flagged = 0,
    hidden = 0,
    total = 0;
  st.forEach((v) => {
    total++;
    if (v.flagged) flagged++;
    if (v.hide) hidden++;
  });
  return { flagged, hidden, total, rules: TT_FLAG_RULES.filter((r) => r.on).length };
}

/* ── Shared visual: the badge + its tooltip ──────────────────────────────── */
/*  Deliberately mirrors _mergedSampleBadgeHTML() so a sample name can carry
    both markers without them looking like two different design languages.  */
function _flagBadgeHTML(item, opts) {
  const st = ttFlagStateFor(item);
  if (!st || !st.flagged) return "";
  const o = opts || {};
  const n = st.hits.length;
  const name = item && typeof item === "object" ? String(item["Specimen ID"] || "") : String(item || "");
  const title = ttFlagPlainReasons(name);
  const icon = st.hide ? "fa-eye-slash" : "fa-flag";
  const label = o.compact ? "" : st.hide ? " hidden" : " flagged";
  return (
    `<span class="tt-flag-badge${st.hide ? " tt-flag-hidden" : ""}" data-flag-sample="${_ttFlagEsc(name)}" ` +
    `title="${_ttFlagEsc(title)}">` +
    `<i class="fas ${icon}" aria-hidden="true"></i>${label}${n > 1 ? " ×" + n : ""}</span>`
  );
}

/** Plain-text reason list (for `title=` attributes and SVG <title> nodes). */
function ttFlagPlainReasons(sample) {
  const st = ttFlagStateFor(sample);
  if (!st || !st.flagged) return "";
  const head = st.hide ? `${sample} — flagged and hidden by a sample QC rule:` : `${sample} — flagged by sample QC:`;
  return [head].concat(st.hits.map((h) => "• " + (h.sample ? h.sample + ": " : "") + h.text)).join("\n");
}

/** Rich tooltip body for showTip(). */
function ttFlagTipHTML(sample) {
  const st = ttFlagStateFor(sample);
  if (!st || !st.flagged) return "";
  return (
    `<b><i class="fas ${st.hide ? "fa-eye-slash" : "fa-flag"}"></i> ${_ttFlagEsc(sample)}</b>` +
    `<br><span style="color:#ccc;font-size:0.9em">${
      st.hide ? "Hidden from every view by a sample QC rule." : "Flagged by sample QC — still shown everywhere."
    }</span><ul style="margin:.35em 0 0;padding-left:1.1em;color:#eee;font-size:0.9em">` +
    st.hits.map((h) => `<li>${_ttFlagEsc((h.sample ? h.sample + ": " : "") + h.text)}</li>`).join("") +
    `</ul>`
  );
}

/* ── Applying the "hide" action ──────────────────────────────────────────── */
/*  Same contract as _applyFilterAutoHide(): only ever touch samples WE hid,
    so a sample the user hid by hand (or the sidebar search hid) is left alone.
    Returns true when visibility actually changed, so the caller can redraw. */
function ttFlagApplyHide() {
  const st = ttFlagEvaluate();
  let changed = false;
  _ttFlagAllSamples().forEach((id) => {
    const want = !!(st.get(id) || {}).hide;
    if (want) {
      if (!sampleHidden[id]) {
        sampleHidden[id] = true;
        _flagAutoHidden.add(id);
        changed = true;
      } else {
        //  Already hidden for another reason — remember nothing, so we never
        //  "restore" something we did not hide.
      }
    } else if (_flagAutoHidden.has(id)) {
      //  Only give visibility back if the search auto-hide isn't also holding
      //  it down; otherwise clearing a QC rule would undo the search filter.
      if (typeof _filterAutoHidden === "undefined" || !_filterAutoHidden.has(id)) sampleHidden[id] = false;
      _flagAutoHidden.delete(id);
      changed = true;
    }
  });
  if (changed && typeof _invalidateFilterCache === "function") _invalidateFilterCache();
  return changed;
}

/** Re-evaluate, re-apply hiding and repaint everything that shows flags.
 *  The single entry point the UI calls after any rule edit. */
function ttFlagRefresh(opts) {
  const o = opts || {};
  ttFlagInvalidate();
  const changed = ttFlagApplyHide();
  if (typeof ttFlagRenderSummary === "function") ttFlagRenderSummary();
  if (o.rebuildList !== false && typeof buildSampleList === "function") buildSampleList();
  if (o.redraw !== false && typeof redraw === "function") redraw();
  else if (changed && typeof redraw === "function") redraw();
  if (typeof _buildRunMetaTable === "function" && activeTab === "runmeta") _buildRunMetaTable();
}

/* ── Defaults from the pipeline ──────────────────────────────────────────── */
/*  make_report.py serializes params.report_flag_* into BOOT.sample_flags:
      { enabled, logic, missing_fails, rules: [{source, field, op, value, …}] }
    Anything malformed is skipped rather than throwing — a bad param must not
    take the whole report down.                                             */
function ttFlagLoadConfig(cfg) {
  TT_FLAG_RULES = [];
  if (!cfg || typeof cfg !== "object") {
    TT_FLAG_ENABLED = true;
    TT_FLAG_LOGIC = "any";
    TT_FLAG_MISSING_FAILS = false;
    TT_FLAG_HIDE_ALL = false;
    return;
  }
  TT_FLAG_ENABLED = cfg.enabled !== false;
  TT_FLAG_LOGIC = cfg.logic === "all" ? "all" : "any";
  TT_FLAG_MISSING_FAILS = !!cfg.missing_fails;
  TT_FLAG_HIDE_ALL = !!cfg.hide_all;
  (Array.isArray(cfg.rules) ? cfg.rules : []).forEach((r) => {
    //  Reject anything that could never match rather than installing a rule
    //  that silently does nothing: an unknown source resolves to undefined and
    //  an unknown operator always compares false, so either one would sit in
    //  the dialog looking active while catching zero samples.
    if (!r || !r.source || !r.field || !r.op) return;
    if (!_TT_FLAG_SOURCE_KEYS.has(String(r.source))) {
      console.warn("[taxtriage] ignoring sample QC rule with unknown source:", r);
      return;
    }
    if (!TT_FLAG_OPS.some((o) => o.op === String(r.op))) {
      console.warn("[taxtriage] ignoring sample QC rule with unknown operator:", r);
      return;
    }
    TT_FLAG_RULES.push(
      ttFlagNewRule({
        on: r.on !== false,
        source: String(r.source),
        field: String(r.field),
        op: String(r.op),
        value: r.value == null ? "" : String(r.value),
        agg: r.agg || "max",
        tass: r.tass != null ? Number(r.tass) : 75,
        action: r.action === "hide" ? "hide" : "flag",
      }),
    );
  });
}

/** Serialize for session export. */
function ttFlagCaptureConfig() {
  return {
    enabled: TT_FLAG_ENABLED,
    logic: TT_FLAG_LOGIC,
    missing_fails: TT_FLAG_MISSING_FAILS,
    hide_all: TT_FLAG_HIDE_ALL,
    rules: TT_FLAG_RULES.map((r) => ({
      on: r.on,
      source: r.source,
      field: r.field,
      op: r.op,
      value: r.value,
      agg: r.agg,
      tass: r.tass,
      action: r.action,
    })),
  };
}

/** Called once from __ttRunInit(), before the first buildSampleList/redraw. */
function ttFlagsInit() {
  const boot = (typeof BOOT !== "undefined" && BOOT && BOOT.sample_flags) || null;
  ttFlagLoadConfig(boot);
  ttFlagInvalidate();
  ttFlagApplyHide();
  if (typeof ttFlagRenderSummary === "function") ttFlagRenderSummary();
}
