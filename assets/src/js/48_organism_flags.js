/* ═══════════════════════════════════════════════════════════════════════════
       -  §  ORGANISM FLAGS  (per-detection QC rules)
       -     The sample rules in 41_sample_flags.js judge WHOLE samples. This
       -     section does the same job one level down: it judges individual
       -     detections (organism rows) and either HIGHLIGHTS them (flag) or
       -     REMOVES them from every view (hide).
       -
       -     Model:  one rule = one or more CONDITIONS that must all hold
       -             (AND), plus an action. Rules are independent of each other
       -             (OR): a row is flagged when any rule trips. That lets a
       -             single rule say "genus is Streptococcus AND fewer than 50
       -             reads" while another says "≥ 99 % ANI to a stronger hit".
       -
       -     Condition sources:
       -       col      any detection column on the row itself (reads, TASS,
       -                coverage, Genus, Family, Microbial Category …)
       -       derived  computed from the row's sample context — ANI to another
       -                hit in the same sample, the genus's reads in the
       -                sample, this organism's share / rank within its genus,
       -                aligned-vs-classifier support, prevalence, lineage text
       -
       -     Everything is evaluated against the FULL dataset (DATA), not the
       -     filtered view, so a flag never flickers as the TASS slider moves.
       -     Peers for the in-sample figures are rows of the same sample at the
       -     same taxonomic level (strain vs strain, species vs species).
       -
       -     "hide" is a gate in filteredData()'s base predicate, so every tab,
       -     export and aggregation honours it. Defaults ship from the pipeline
       -     (params.report_org_flag_* → BOOT.organism_flags).
═══════════════════════════════════════════════════════════════════════════ */

/* ── State ───────────────────────────────────────────────────────────────── */
let TT_OFLAG_RULES = []; // [{ id, on, action, conds: [{ id, source, field, op, value, peer }] }]
let TT_OFLAG_ENABLED = true;
/*  View mode — the sidebar's Highlight / Hide flagged / Only flagged control.
      "all"   highlight flagged rows; a rule whose action is "hide" still hides
      "hide"  every flagged row is removed from the views
      "only"  show nothing but the flagged rows (per-rule hides are ignored)   */
let TT_OFLAG_VIEW = "all";
let _TT_OFLAG_SEQ = 0;

const TT_OFLAG_VIEWS = [
  { key: "all", label: "Highlight flagged" },
  { key: "hide", label: "Hide flagged" },
  { key: "only", label: "Only flagged" },
];

function ttOFlagNormalizeView(v) {
  const s = String(v == null ? "" : v).toLowerCase();
  return TT_OFLAG_VIEWS.some((m) => m.key === s) ? s : "all";
}

/* ── Operators (same vocabulary as the sample rules) ─────────────────────── */
const TT_OFLAG_OPS = [
  { op: "<", label: "is less than", num: true },
  { op: "<=", label: "is at most", num: true },
  { op: ">", label: "is greater than", num: true },
  { op: ">=", label: "is at least", num: true },
  { op: "==", label: "equals", num: null },
  { op: "!=", label: "does not equal", num: null },
  { op: "contains", label: "contains", num: false },
  { op: "!contains", label: "does not contain", num: false },
  { op: "in", label: "is one of", num: false },
  { op: "!in", label: "is not one of", num: false },
  { op: "regex", label: "matches regex", num: false },
  { op: "empty", label: "is empty / missing", num: null, novalue: true },
  { op: "!empty", label: "has any value", num: null, novalue: true },
];
const _TT_OFLAG_OP_KEYS = new Set(TT_OFLAG_OPS.map((o) => o.op));

/* ── Field catalogue ─────────────────────────────────────────────────────── */
/*  `needsPeer`: the ANI fields take a partner qualifier.
      stronger  only count a partner that out-ranks this row (more reads, then
                higher TASS, then a stable tie-break) — so of two near-identical
                references only the WEAKER one trips, and hiding it never
                removes both.
      any       any other hit in the sample.                                  */
const TT_OFLAG_DERIVED_FIELDS = [
  {
    key: "ani_max",
    label: "Shared ANI % with another hit in this sample",
    num: true,
    needsPeer: true,
    ani: true,
  },
  { key: "ani_partners", label: "# hits in this sample sharing high ANI", num: true, needsPeer: true, ani: true },
  { key: "genus_reads", label: "Genus reads in this sample", num: true },
  { key: "genus_share", label: "% of its genus's reads in this sample", num: true },
  { key: "genus_rank", label: "Rank within its genus (by reads, 1 = top)", num: true },
  { key: "genus_members", label: "# hits in the same genus in this sample", num: true },
  { key: "sample_read_share", label: "% of the sample's organism-aligned reads", num: true },
  { key: "k2_ratio", label: "Aligned ÷ classifier (K2) reads", num: true },
  { key: "prevalence", label: "# samples this organism is detected in", num: true },
  { key: "lineage", label: "Taxonomy lineage (any rank)", num: false },
];

const TT_OFLAG_PEERS = [
  { key: "stronger", label: "stronger hit" },
  { key: "any", label: "any hit" },
];

/*  Row fields that are structural or nested, never worth a condition. */
const _TT_OFLAG_COL_SKIP = new Set(["High ANI Matches", "ANI Annotated"]);

/*  Columns offered first because they are what people actually filter on. */
const _TT_OFLAG_COL_FIRST = [
  "# Reads Aligned",
  "TASS Score",
  "Genus",
  "Detected Organism",
  "Species Name",
  "Family",
  "Microbial Category",
  "Breadth %",
  "Coverage",
  "K2 Reads",
  "% Reads",
  "RPM",
  "Taxonomic ID #",
];

function ttOFlagFields() {
  const all = ((typeof ALL_COLS !== "undefined" && ALL_COLS) || []).filter((c) => !_TT_OFLAG_COL_SKIP.has(c));
  const isNum = (c) => typeof NUMERIC !== "undefined" && NUMERIC && NUMERIC.has(c);
  const ordered = _TT_OFLAG_COL_FIRST
    .filter((c) => all.includes(c))
    .concat(all.filter((c) => !_TT_OFLAG_COL_FIRST.includes(c)));
  //  A handful of columns are numeric in every run even when NUMERIC was built
  //  from a sample that lacked them; don't offer "contains" on read counts.
  const knownNum = new Set(["# Reads Aligned", "TASS Score", "K2 Reads", "Breadth %", "Coverage", "% Reads", "RPM"]);
  const col = ordered.map((c) => ({ key: c, label: c, num: isNum(c) || knownNum.has(c) }));
  return { col, derived: TT_OFLAG_DERIVED_FIELDS.slice() };
}

function ttOFlagFieldDef(cond) {
  if (!cond) return null;
  const list = ttOFlagFields()[cond.source] || [];
  return (
    list.find((f) => f.key === cond.field) || {
      key: cond.field,
      label: cond.field,
      num: cond.source === "derived" ? true : null,
    }
  );
}

/* ── Rule helpers ────────────────────────────────────────────────────────── */
function ttOFlagNewCond(seed) {
  const c = Object.assign(
    { id: "c" + ++_TT_OFLAG_SEQ, source: "col", field: "# Reads Aligned", op: "<", value: "10", peer: "stronger" },
    seed || {},
  );
  c.id = "c" + ++_TT_OFLAG_SEQ;
  return c;
}

function ttOFlagNewRule(seed) {
  const s = seed || {};
  const conds = (Array.isArray(s.conds) && s.conds.length ? s.conds : [{}]).map((c) => ttOFlagNewCond(c));
  return {
    id: "o" + ++_TT_OFLAG_SEQ,
    on: s.on !== false,
    action: s.action === "hide" ? "hide" : "flag",
    conds,
  };
}

function _ttOFlagEsc(v) {
  return String(v == null ? "" : v)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}
function _ttOFlagFmt(v) {
  if (typeof v === "number") {
    if (!isFinite(v)) return String(v);
    return Math.abs(v) >= 1000 || Number.isInteger(v) ? Math.round(v).toLocaleString() : String(+v.toFixed(3));
  }
  const n = Number(v);
  if (v !== "" && v != null && isFinite(n) && String(v).trim() !== "") return n.toLocaleString();
  return String(v == null ? "" : v);
}

function ttOFlagCondLabel(cond) {
  const def = ttOFlagFieldDef(cond) || { label: cond.field };
  const opDef = TT_OFLAG_OPS.find((o) => o.op === cond.op) || { label: cond.op };
  let name = def.label;
  if (def.needsPeer) name += cond.peer === "any" ? " (any hit)" : " (stronger hit)";
  return opDef.novalue ? `${name} ${opDef.label}` : `${name} ${opDef.label} ${_ttOFlagFmt(cond.value)}`;
}

function ttOFlagRuleLabel(rule) {
  if (!rule) return "";
  return (rule.conds || []).map(ttOFlagCondLabel).join(" and ");
}

/* ── Row identity ────────────────────────────────────────────────────────── */
function _ttOFlagOrgId(r) {
  return String(r["Taxonomic ID #"] || r["Detected Organism"] || "");
}
function _ttOFlagKey(r) {
  return [String(r["Specimen ID"] || ""), r["Level"] || "Strain", _ttOFlagOrgId(r), String(r["Subkey"] || "")].join(
    "\u0001",
  );
}

/* ── ANI helpers ─────────────────────────────────────────────────────────── */
function _ttOFlagAniAnnotated(r) {
  if (r["ANI Annotated"] === true) return true;
  if (r["ANI Annotated"] === false) return false;
  return Object.prototype.hasOwnProperty.call(r, "High ANI Matches");
}
function _ttOFlagAniList(r) {
  const v = r["High ANI Matches"];
  if (Array.isArray(v)) {
    return v
      .map((m) => ({ key: String((m && m.key) || ""), pct: parseFloat((m && m.ani_pct) || 0) || 0 }))
      .filter((m) => m.key);
  }
  if (typeof v === "string" && v.trim()) {
    return v
      .split(";")
      .map((s) => {
        const mm = s.trim().match(/^(.+?)\(([\d.]+)%?\)$/);
        if (mm) return { key: mm[1].trim(), pct: parseFloat(mm[2]) || 0 };
        return null;
      })
      .filter(Boolean);
  }
  return [];
}

/* ── Context (per sample × level) ────────────────────────────────────────── */
function _ttOFlagNum(v) {
  const n = parseFloat(v);
  return isNaN(n) ? 0 : n;
}

/*  "Does `b` out-rank `a`?" — reads first, then TASS, then a stable tie-break
    on the organism id so exactly one of an identical pair counts as stronger. */
function _ttOFlagStronger(b, a) {
  const rb = _ttOFlagNum(b["# Reads Aligned"]),
    ra = _ttOFlagNum(a["# Reads Aligned"]);
  if (rb !== ra) return rb > ra;
  const tb = _ttOFlagNum(b["TASS Score"]),
    ta = _ttOFlagNum(a["TASS Score"]);
  if (tb !== ta) return tb > ta;
  return _ttOFlagOrgId(b) < _ttOFlagOrgId(a);
}

let _TT_OFLAG_CTX = { key: null, value: null };

function _ttOFlagContext() {
  const src = (typeof DATA !== "undefined" && DATA) || [];
  const key = src.length + "|" + (typeof SAMPLE_META_EPOCH !== "undefined" ? SAMPLE_META_EPOCH : 0);
  if (_TT_OFLAG_CTX.key === key && _TT_OFLAG_CTX.value) return _TT_OFLAG_CTX.value;
  const groups = new Map(); // sample|level → { rows, byId: Map, genus: Map(genus → {reads, rows}) , reads }
  const prevalence = new Map(); // level|orgId → Set(sample)
  for (let i = 0; i < src.length; i++) {
    const r = src[i];
    const s = String(r["Specimen ID"] || "");
    if (!s) continue;
    const lvl = r["Level"] || "Strain";
    const gk = s + "\u0001" + lvl;
    let g = groups.get(gk);
    if (!g) groups.set(gk, (g = { rows: [], byId: new Map(), genus: new Map(), reads: 0 }));
    g.rows.push(r);
    const oid = _ttOFlagOrgId(r);
    if (oid && !g.byId.has(oid)) g.byId.set(oid, r);
    const reads = _ttOFlagNum(r["# Reads Aligned"]);
    g.reads += reads;
    const gen = String(r["Genus"] || r["Genus Name"] || "").trim();
    if (gen) {
      let ge = g.genus.get(gen);
      if (!ge) g.genus.set(gen, (ge = { reads: 0, rows: [] }));
      ge.reads += reads;
      ge.rows.push(r);
    }
    const pk = lvl + "\u0001" + oid;
    let ps = prevalence.get(pk);
    if (!ps) prevalence.set(pk, (ps = new Set()));
    ps.add(s);
  }
  //  Genus rank is a sort per genus; do it once here rather than per rule.
  groups.forEach((g) => {
    g.genus.forEach((ge) => {
      ge.rows.sort((a, b) => (_ttOFlagStronger(a, b) ? -1 : _ttOFlagStronger(b, a) ? 1 : 0));
      ge.rank = new Map();
      ge.rows.forEach((r, i) => ge.rank.set(r, i + 1));
    });
  });
  _TT_OFLAG_CTX = { key, value: { groups, prevalence } };
  return _TT_OFLAG_CTX.value;
}

/*  Returns { v, detail } — `detail` is extra wording for the reason text
    (e.g. which partner the ANI came from). v === undefined means missing. */
function _ttOFlagDerived(cond, r, ctx) {
  const s = String(r["Specimen ID"] || "");
  const lvl = r["Level"] || "Strain";
  const g = ctx.groups.get(s + "\u0001" + lvl);
  const gen = String(r["Genus"] || r["Genus Name"] || "").trim();
  const ge = g && gen ? g.genus.get(gen) : null;
  const reads = _ttOFlagNum(r["# Reads Aligned"]);
  switch (cond.field) {
    case "ani_max":
    case "ani_partners": {
      if (!_ttOFlagAniAnnotated(r)) return { v: undefined, detail: "no ANI data (run with --enable_matrix)" };
      const me = _ttOFlagOrgId(r);
      let best = 0,
        bestRow = null,
        n = 0;
      _ttOFlagAniList(r).forEach((m) => {
        if (!g || m.key === me) return;
        const p = g.byId.get(m.key);
        if (!p || p === r) return;
        if (cond.peer !== "any" && !_ttOFlagStronger(p, r)) return;
        if (m.pct > best) {
          best = m.pct;
          bestRow = p;
        }
        n++; // every partner the pipeline recorded is already ≥ --ani_threshold
      });
      if (cond.field === "ani_partners") return { v: n };
      return {
        v: best,
        detail: bestRow ? `vs ${bestRow["Detected Organism"] || _ttOFlagOrgId(bestRow)}` : "",
      };
    }
    case "genus_reads":
      return { v: ge ? ge.reads : undefined };
    case "genus_share":
      return { v: ge && ge.reads > 0 ? (reads / ge.reads) * 100 : undefined, detail: gen ? `of ${gen}` : "" };
    case "genus_rank":
      return { v: ge ? ge.rank.get(r) : undefined, detail: ge ? `of ${ge.rows.length} in ${gen}` : "" };
    case "genus_members":
      return { v: ge ? ge.rows.length : undefined };
    case "sample_read_share":
      return { v: g && g.reads > 0 ? (reads / g.reads) * 100 : undefined };
    case "k2_ratio": {
      const k2 = _ttOFlagNum(r["K2 Reads"]);
      return { v: k2 > 0 ? reads / k2 : undefined };
    }
    case "prevalence": {
      const ps = ctx.prevalence.get(lvl + "\u0001" + _ttOFlagOrgId(r));
      return { v: ps ? ps.size : 0 };
    }
    case "lineage": {
      const parts = [
        "Domain",
        "Superkingdom",
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species Name",
        "Detected Organism",
      ]
        .map((k) => String(r[k] == null ? "" : r[k]).trim())
        .filter(Boolean);
      return { v: parts.join("; ") };
    }
  }
  return { v: undefined };
}

function _ttOFlagValue(cond, r, ctx) {
  if (cond.source === "derived") return _ttOFlagDerived(cond, r, ctx);
  const v = r[cond.field];
  return { v: Array.isArray(v) ? v.join(", ") : v };
}

function _ttOFlagIsMissing(v) {
  if (v === undefined || v === null) return true;
  if (typeof v === "string" && v.trim() === "") return true;
  if (typeof v === "number" && !isFinite(v)) return true;
  return false;
}

function _ttOFlagList(raw) {
  return String(raw == null ? "" : raw)
    .split(/[,;|\n]+/)
    .map((s) => s.trim().toLowerCase())
    .filter(Boolean);
}

/** Missing values never satisfy a comparison — only "is empty" matches them. */
function _ttOFlagCompare(actual, cond) {
  const op = cond.op;
  if (op === "empty") return _ttOFlagIsMissing(actual);
  if (op === "!empty") return !_ttOFlagIsMissing(actual);
  if (_ttOFlagIsMissing(actual)) return false;
  const raw = cond.value == null ? "" : String(cond.value);
  if (op === "<" || op === "<=" || op === ">" || op === ">=") {
    const a = Number(actual),
      b = Number(raw);
    if (!isFinite(a) || !isFinite(b) || raw.trim() === "") return false;
    if (op === "<") return a < b;
    if (op === "<=") return a <= b;
    if (op === ">") return a > b;
    return a >= b;
  }
  if (op === "==" || op === "!=") {
    const a = Number(actual),
      b = Number(raw);
    let eq;
    if (isFinite(a) && isFinite(b) && String(actual).trim() !== "" && raw.trim() !== "") eq = a === b;
    else eq = String(actual).trim().toLowerCase() === raw.trim().toLowerCase();
    return op === "==" ? eq : !eq;
  }
  if (op === "contains" || op === "!contains") {
    const needle = raw.trim().toLowerCase();
    if (!needle) return false; // an empty box mid-edit must not match everything
    const hit = String(actual).toLowerCase().includes(needle);
    return op === "contains" ? hit : !hit;
  }
  if (op === "in" || op === "!in") {
    const list = _ttOFlagList(raw);
    if (!list.length) return false;
    const hit = list.includes(String(actual).trim().toLowerCase());
    return op === "in" ? hit : !hit;
  }
  if (op === "regex") {
    if (!raw) return false;
    try {
      return new RegExp(raw, "i").test(String(actual));
    } catch (e) {
      return false;
    }
  }
  return false;
}

/* ── Evaluation ──────────────────────────────────────────────────────────── */
let _TT_OFLAG_EVAL = { key: null, value: null };

function ttOFlagInvalidate() {
  _TT_OFLAG_EVAL = { key: null, value: null };
  _TT_OFLAG_CTX = { key: null, value: null };
}

function _ttOFlagRulesSig() {
  return TT_OFLAG_RULES.map(
    (r) =>
      (r.on ? 1 : 0) +
      "~" +
      r.action +
      "~" +
      (r.conds || []).map((c) => [c.source, c.field, c.op, c.value, c.peer].join("^")).join("&"),
  ).join("|");
}

/** Cheap fingerprint for filteredData()'s cache key. Only states that can
 *  HIDE something need to be part of it; pure highlighting never changes rows. */
function ttOFlagFilterKey() {
  if (!TT_OFLAG_ENABLED || !TT_OFLAG_RULES.some((r) => r.on)) return "";
  return TT_OFLAG_VIEW + "#" + _ttOFlagRulesSig();
}

function _ttOFlagActiveRules() {
  if (!TT_OFLAG_ENABLED) return [];
  return TT_OFLAG_RULES.filter((r) => r.on && r.conds && r.conds.length);
}

/** Map(rowKey → { flagged, hide, hits: [{ rule, text }] }) for FLAGGED rows
 *  only, plus `.anyFlagged` / `.onlyActive` on the map itself. */
function ttOFlagEvaluate() {
  const src = (typeof DATA !== "undefined" && DATA) || [];
  const key = [TT_OFLAG_ENABLED ? 1 : 0, TT_OFLAG_VIEW, src.length, _ttOFlagRulesSig()].join("#");
  if (_TT_OFLAG_EVAL.key === key && _TT_OFLAG_EVAL.value) return _TT_OFLAG_EVAL.value;
  const out = new Map();
  const active = _ttOFlagActiveRules();
  if (active.length) {
    const ctx = _ttOFlagContext();
    for (let i = 0; i < src.length; i++) {
      const r = src[i];
      let hits = null;
      for (let j = 0; j < active.length; j++) {
        const rule = active[j];
        const parts = [];
        let ok = true;
        for (let k = 0; k < rule.conds.length; k++) {
          const c = rule.conds[k];
          const res = _ttOFlagValue(c, r, ctx);
          if (!_ttOFlagCompare(res.v, c)) {
            ok = false;
            break;
          }
          const opDef = TT_OFLAG_OPS.find((o) => o.op === c.op) || {};
          parts.push(
            ttOFlagCondLabel(c) +
              (opDef.novalue ? "" : ` (actual: ${_ttOFlagFmt(res.v)}${res.detail ? " " + res.detail : ""})`),
          );
        }
        if (!ok) continue;
        (hits || (hits = [])).push({ rule, text: parts.join(" and ") });
      }
      if (!hits) continue;
      const hide =
        TT_OFLAG_VIEW === "only" ? false : TT_OFLAG_VIEW === "hide" || hits.some((h) => h.rule.action === "hide");
      out.set(_ttOFlagKey(r), { flagged: true, hide, hits, row: r });
    }
  }
  out.anyFlagged = out.size > 0;
  //  "Only flagged" with nothing flagged would blank the report; fall back to
  //  showing everything, exactly as the sample rules do.
  out.onlyActive = TT_OFLAG_VIEW === "only" && out.anyFlagged;
  _TT_OFLAG_EVAL = { key, value: out };
  return out;
}

/** Gate used by filteredData(): should this RAW data row be dropped? */
function ttOFlagRowHidden(r) {
  if (!TT_OFLAG_ENABLED) return false;
  const st = ttOFlagEvaluate();
  if (!st.anyFlagged) return false;
  const e = st.get(_ttOFlagKey(r));
  if (st.onlyActive) return !e;
  return !!(e && e.hide);
}

/** Flag state for any row shown in a view — including a specimen-merged row,
 *  whose "Specimen ID" is a group name: it inherits its member libraries'. */
function ttOFlagStateFor(r) {
  if (!r || typeof r !== "object") return null;
  const st = ttOFlagEvaluate();
  if (!st.anyFlagged) return null;
  const direct = st.get(_ttOFlagKey(r));
  if (direct) return direct;
  if (typeof specimenGroups === "function") {
    let members = null;
    try {
      members = specimenGroups().get(String(r["Specimen ID"] || ""));
    } catch (e) {
      members = null;
    }
    if (members && members.length) {
      const hits = [];
      let hide = false;
      members.forEach((m) => {
        const k = _ttOFlagKey(Object.assign({}, r, { "Specimen ID": m }));
        const e = st.get(k);
        if (!e) return;
        hide = hide || e.hide;
        e.hits.forEach((h) => hits.push(Object.assign({ sample: m }, h)));
      });
      if (hits.length) return { flagged: true, hide, hits, merged: true };
    }
  }
  return null;
}

/* ── Counts / reasons ────────────────────────────────────────────────────── */
function _ttOFlagViewLevel() {
  try {
    const el = typeof document !== "undefined" && document.getElementById("view-level");
    return (el && el.value) || "Strain";
  } catch (e) {
    return "Strain";
  }
}

/** Counts at the taxonomic level currently on screen, so "12 flagged" matches
 *  the rows the reader can actually see. */
function ttOFlagCounts() {
  const st = ttOFlagEvaluate();
  const lvl = _ttOFlagViewLevel();
  const src = (typeof DATA !== "undefined" && DATA) || [];
  let total = 0,
    flagged = 0,
    hidden = 0;
  const orgs = new Set();
  for (let i = 0; i < src.length; i++) {
    const r = src[i];
    if ((r["Level"] || "Strain") !== lvl) continue;
    total++;
    const e = st.get(_ttOFlagKey(r));
    if (e) {
      flagged++;
      orgs.add(_ttOFlagOrgId(r));
      if (e.hide) hidden++;
    }
  }
  if (st.onlyActive) hidden = total - flagged;
  return {
    total,
    flagged,
    hidden,
    organisms: orgs.size,
    level: lvl,
    view: TT_OFLAG_VIEW,
    onlyActive: st.onlyActive,
    rules: _ttOFlagActiveRules().length,
  };
}

function ttOFlagPlainReasons(r) {
  const st = ttOFlagStateFor(r);
  if (!st) return "";
  const head = `${r["Detected Organism"] || ""} — ${st.hide ? "flagged and hidden" : "flagged"} by organism QC:`;
  return [head].concat(st.hits.map((h) => "• " + (h.sample ? h.sample + ": " : "") + h.text)).join("\n");
}

function ttOFlagTipHTML(r) {
  const st = ttOFlagStateFor(r);
  if (!st) return "";
  return (
    `<b><i class="fas fa-flag"></i> ${_ttOFlagEsc(r["Detected Organism"] || "")}</b>` +
    `<br><span style="color:#ccc;font-size:0.9em">${
      st.hide ? "Flagged and hidden by an organism QC rule." : "Flagged by organism QC — still shown."
    }</span><ul style="margin:.35em 0 0;padding-left:1.1em;color:#eee;font-size:0.9em">` +
    st.hits.map((h) => `<li>${_ttOFlagEsc((h.sample ? h.sample + ": " : "") + h.text)}</li>`).join("") +
    `</ul>`
  );
}

/** Inline badge next to an organism name (Table / Summary tables). */
function _orgFlagBadgeHTML(r) {
  const st = ttOFlagStateFor(r);
  if (!st) return "";
  const n = st.hits.length;
  return (
    `<span class="tt-oflag-badge${st.hide ? " tt-flag-hidden" : ""}" title="${_ttOFlagEsc(ttOFlagPlainReasons(r))}">` +
    `<i class="fas fa-flag" aria-hidden="true"></i> QC${n > 1 ? " ×" + n : ""}</span>`
  );
}

/** Text for the export column. */
function ttOFlagExportText(r) {
  const st = ttOFlagStateFor(r);
  if (!st) return "";
  const t = st.hits.map((h) => h.text).join(" | ");
  return st.hide ? "hidden: " + t : t;
}

/* ── Refresh ─────────────────────────────────────────────────────────────── */
function ttOFlagRefresh() {
  ttOFlagInvalidate();
  if (typeof _invalidateFilterCache === "function") _invalidateFilterCache();
  if (typeof ttOFlagRenderSummary === "function") ttOFlagRenderSummary();
  if (typeof redraw === "function") redraw();
}

/* ── Config load / capture ───────────────────────────────────────────────── */
const _TT_OFLAG_SOURCES = new Set(["col", "derived"]);

function _ttOFlagCleanCond(c) {
  if (!c || typeof c !== "object") return null;
  const source = String(c.source || "col");
  const field = String(c.field || "").trim();
  const op = String(c.op || "").trim();
  if (!_TT_OFLAG_SOURCES.has(source) || !field || !_TT_OFLAG_OP_KEYS.has(op)) {
    if (typeof console !== "undefined") console.warn("[taxtriage] ignoring organism QC condition:", c);
    return null;
  }
  if (source === "derived" && !TT_OFLAG_DERIVED_FIELDS.some((f) => f.key === field)) {
    if (typeof console !== "undefined") console.warn("[taxtriage] unknown derived organism field:", c);
    return null;
  }
  return {
    source,
    field,
    op,
    value: c.value == null ? "" : String(c.value),
    peer: c.peer === "any" ? "any" : "stronger",
  };
}

function ttOFlagLoadConfig(cfg) {
  TT_OFLAG_RULES = [];
  TT_OFLAG_ENABLED = true;
  TT_OFLAG_VIEW = "all";
  if (!cfg || typeof cfg !== "object") return;
  TT_OFLAG_ENABLED = cfg.enabled !== false;
  TT_OFLAG_VIEW = ttOFlagNormalizeView(cfg.view);
  (Array.isArray(cfg.rules) ? cfg.rules : []).forEach((r) => {
    if (!r || typeof r !== "object") return;
    //  A flat rule ({source, field, op, value}) is a one-condition rule.
    const rawConds = Array.isArray(r.conds) ? r.conds : Array.isArray(r.conditions) ? r.conditions : [r];
    const conds = rawConds.map(_ttOFlagCleanCond).filter(Boolean);
    if (!conds.length || conds.length !== rawConds.length) return; // never half-install a rule
    TT_OFLAG_RULES.push(ttOFlagNewRule({ on: r.on !== false, action: r.action, conds }));
  });
}

function ttOFlagCaptureConfig() {
  return {
    enabled: TT_OFLAG_ENABLED,
    view: TT_OFLAG_VIEW,
    rules: TT_OFLAG_RULES.map((r) => ({
      on: r.on,
      action: r.action,
      conds: r.conds.map((c) => {
        const o = { source: c.source, field: c.field, op: c.op, value: c.value };
        //  Only the ANI fields read `peer`; leave it off everything else.
        if (c.source === "derived" && (c.field === "ani_max" || c.field === "ani_partners")) o.peer = c.peer;
        return o;
      }),
    })),
  };
}

/** Called once from __ttRunInit(), before the first redraw. */
function ttOFlagsInit() {
  const boot = (typeof BOOT !== "undefined" && BOOT && BOOT.organism_flags) || null;
  ttOFlagLoadConfig(boot);
  ttOFlagInvalidate();
  if (typeof _invalidateFilterCache === "function") _invalidateFilterCache();
  if (typeof ttOFlagRenderSummary === "function") ttOFlagRenderSummary();
}
