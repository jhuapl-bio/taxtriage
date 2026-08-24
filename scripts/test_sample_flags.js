// Unit-test the sample QC rule engine (assets/src/js/41_sample_flags.js).
//
// 41_sample_flags.js is pure declarations + functions with no top-level side
// effects, so it can be evaluated whole inside a vm context with stubs for the
// report globals it reads. Everything below exercises the ENGINE only; the
// dialog (42_sample_flags_ui.js) needs a DOM and is covered by the browser
// checks instead.
//
//   node scripts/test_sample_flags.js
const fs = require("fs"),
  vm = require("vm");

let fails = 0;
const ck = (n, c, e) => {
  console.log((c ? "PASS  " : "FAIL  ") + n + (e !== undefined ? "   " + JSON.stringify(e) : ""));
  if (!c) fails++;
};

// ── fixture ───────────────────────────────────────────────────────────────
// S1  shallow (5k reads), 2 strains at TASS 90 / 80, one Species roll-up row
//     that must NOT be counted as an organism
// S2  deep (50k reads), 3 strains at TASS 99 / 98 / 50
// S3  metadata only: no detections, no total_reads
function makeCtx() {
  const row = (s, tax, tass, extra) =>
    Object.assign(
      {
        "Specimen ID": s,
        "Taxonomic ID #": String(tax),
        "Detected Organism": "Org" + tax,
        "TASS Score": tass,
        Level: "Strain",
        Genus: "G" + String(tax)[0],
        "Passes Threshold": tass >= 90 ? "true" : "false",
        "High Consequence": tass >= 99 ? "true" : "false",
        "# Reads Aligned": tass,
        Coverage: tass / 2,
      },
      extra || {},
    );

  const DATA = [
    row("S1", 11, 90),
    row("S1", 12, 80),
    // A species roll-up row for the same organisms. Counting it would inflate
    // every "distinct organism" figure — the engine must skip non-Strain rows.
    row("S1", 19, 95, { Level: "Species" }),
    row("S2", 21, 99),
    row("S2", 22, 98),
    row("S2", 23, 50),
  ];

  const ctx = {
    console,
    DATA,
    ALL_COLS: ["TASS Score", "Coverage", "# Reads Aligned", "Detected Organism"],
    NUMERIC: new Set(["TASS Score", "Coverage", "# Reads Aligned"]),
    SAMPLE_META: {
      S1: { sample_name: "S1", total_reads: 5000, aligned_reads: 100, sample_type: "nasal" },
      S2: { sample_name: "S2", total_reads: 50000, aligned_reads: 40000, sample_type: "stool" },
      S3: { sample_name: "S3", aligned_reads: 0, sample_type: "nasal" },
    },
    RUN_META: [
      { sample_name: "S1", site: "clinic-alpha" },
      { sample_name: "S2", site: "clinic-beta" },
      { sample_name: "S3", site: "" },
    ],
    SAMPLE_META_EPOCH: 0,
    sampleHidden: {},
    BOOT: {},
    // helpers the engine borrows from the rest of the bundle
    _ttMax: (a) => a.reduce((m, v) => (Number(v) > m ? Number(v) : m), -Infinity),
    _ttMin: (a) => a.reduce((m, v) => (Number(v) < m ? Number(v) : m), Infinity),
    num: (v) => (isNaN(parseFloat(v)) ? 0 : parseFloat(v)),
    isTruthy: (v) => ["true", "yes", "1", "t"].includes(String(v).trim().toLowerCase()),
    _invalidateFilterCache: () => {},
  };
  vm.createContext(ctx);
  //  The engine's top-level `let` bindings (TT_FLAG_RULES, TT_FLAG_LOGIC, …)
  //  live in the vm script's lexical scope; only `function` declarations become
  //  properties of the context object. These accessors are declared in the same
  //  scope, so they can read and write those bindings on the test's behalf.
  const ACCESSORS = `
    function __ttSet(o) {
      if (o.rules !== undefined) TT_FLAG_RULES = o.rules;
      if (o.logic !== undefined) TT_FLAG_LOGIC = o.logic;
      if (o.enabled !== undefined) TT_FLAG_ENABLED = o.enabled;
      if (o.missingFails !== undefined) TT_FLAG_MISSING_FAILS = o.missingFails;
      if (o.hideAll !== undefined) TT_FLAG_VIEW = o.hideAll ? "hide" : "all";
      if (o.view !== undefined) TT_FLAG_VIEW = ttFlagNormalizeView(o.view);
      if (o.exclude !== undefined) TT_FLAG_EXCLUDE = ttFlagParseExclude(o.exclude);
    }
    function __ttGet() {
      return { rules: TT_FLAG_RULES, logic: TT_FLAG_LOGIC, enabled: TT_FLAG_ENABLED,
               missingFails: TT_FLAG_MISSING_FAILS, view: TT_FLAG_VIEW,
               hideAll: TT_FLAG_VIEW === "hide", exclude: TT_FLAG_EXCLUDE };
    }`;
  vm.runInContext(fs.readFileSync("assets/src/js/41_sample_flags.js", "utf8") + ACCESSORS, ctx, {
    filename: "41_sample_flags.js",
  });
  return ctx;
}

// Replace the rule set and read back the verdict for every sample.
function verdict(ctx, rules, opts) {
  const o = opts || {};
  ctx.__ttSet({
    rules: rules.map((r) => ctx.ttFlagNewRule(r)),
    logic: o.logic || "any",
    enabled: o.enabled !== false,
    missingFails: !!o.missingFails,
    view: o.view !== undefined ? o.view : o.hideAll ? "hide" : "all",
    exclude: o.exclude,
  });
  ctx.ttFlagInvalidate();
  const out = {};
  ctx.ttFlagEvaluate().forEach((v, k) => (out[k] = v));
  return out;
}
const flaggedNames = (v) =>
  Object.keys(v)
    .filter((k) => v[k].flagged)
    .sort();

// ── 1. sources ────────────────────────────────────────────────────────────
let ctx = makeCtx();

let v = verdict(ctx, [{ source: "meta", field: "total_reads", op: "<", value: "10000" }]);
ck(
  "meta: total_reads < 10000 catches only the shallow sample",
  JSON.stringify(flaggedNames(v)) === '["S1"]',
  flaggedNames(v),
);
ck("meta: a sample with NO total_reads does not trip by default", v.S3.flagged === false);

v = verdict(ctx, [{ source: "derived", field: "unique_taxids_above_tass", op: "<", value: "2", tass: 95 }]);
ck(
  "derived: organisms above TASS 95 — S1 has 0, S2 has 2",
  JSON.stringify(flaggedNames(v)) === '["S1","S3"]',
  flaggedNames(v),
);
ck(
  "derived: the Species roll-up row is not counted as an organism",
  ctx._ttFlagDerivedValue("unique_taxids", ctx._ttFlagRowsBySample().get("S1"), {}) === 2,
  ctx._ttFlagDerivedValue("unique_taxids", ctx._ttFlagRowsBySample().get("S1"), {}),
);
ck(
  "derived: passing_detections counts only rows over their own threshold",
  ctx._ttFlagDerivedValue("passing_detections", ctx._ttFlagRowsBySample().get("S2"), {}) === 2,
);
ck("derived: max_tass", ctx._ttFlagDerivedValue("max_tass", ctx._ttFlagRowsBySample().get("S2"), {}) === 99);

v = verdict(ctx, [{ source: "runmeta", field: "site", op: "contains", value: "alpha" }]);
ck("runmeta: contains matches the right sample", JSON.stringify(flaggedNames(v)) === '["S1"]', flaggedNames(v));

v = verdict(ctx, [{ source: "runmeta", field: "site", op: "empty", value: "" }]);
ck("runmeta: is-empty catches the blank cell", JSON.stringify(flaggedNames(v)) === '["S3"]', flaggedNames(v));

// sample_type lives in SAMPLE_META, never in RUN_META — the runmeta source
// must fall back to it rather than reporting the column as missing.
v = verdict(ctx, [{ source: "runmeta", field: "sample_type", op: "==", value: "nasal" }]);
ck(
  "runmeta: falls back to SAMPLE_META for a column RUN_META lacks",
  JSON.stringify(flaggedNames(v)) === '["S1","S3"]',
  flaggedNames(v),
);

v = verdict(ctx, [{ source: "data", field: "Coverage", agg: "max", op: "<", value: "46" }]);
ck("data: max(Coverage) — S1 max is 45, S2 max is 49.5", JSON.stringify(flaggedNames(v)) === '["S1"]', flaggedNames(v));
// S1 means 42.5 (45, 40); S2 means 41.17 (49.5, 49, 25). max picks S1, mean
// picks S2 — the aggregation genuinely changes the answer.
v = verdict(ctx, [{ source: "data", field: "Coverage", agg: "mean", op: "<", value: "42" }]);
ck(
  "data: mean(Coverage) selects a different sample than max",
  JSON.stringify(flaggedNames(v)) === '["S2"]',
  flaggedNames(v),
);

// ── 2. combining rules ────────────────────────────────────────────────────
const twoRules = [
  { source: "meta", field: "total_reads", op: "<", value: "10000" }, // S1
  { source: "derived", field: "max_tass", op: ">=", value: "99" }, // S2
];
v = verdict(ctx, twoRules, { logic: "any" });
ck("logic any: union of both rules", JSON.stringify(flaggedNames(v)) === '["S1","S2"]', flaggedNames(v));
v = verdict(ctx, twoRules, { logic: "all" });
ck("logic all: no sample trips both", flaggedNames(v).length === 0, flaggedNames(v));

v = verdict(ctx, twoRules, { enabled: false });
ck("master switch off: nothing is flagged", flaggedNames(v).length === 0);

v = verdict(ctx, [Object.assign({ on: false }, twoRules[0])]);
ck("a disabled rule is ignored", flaggedNames(v).length === 0);

// ── 3. missing values ─────────────────────────────────────────────────────
const missingRule = [{ source: "meta", field: "total_reads", op: "<", value: "10000" }];
v = verdict(ctx, missingRule, { missingFails: false });
ck("missing value does not trip by default", v.S3.flagged === false);
v = verdict(ctx, missingRule, { missingFails: true });
ck("missing value trips when TT_FLAG_MISSING_FAILS is on", v.S3.flagged === true);

// ── 4. hide action drives sampleHidden ────────────────────────────────────
ctx = makeCtx();
verdict(ctx, [{ source: "meta", field: "total_reads", op: "<", value: "10000", action: "hide" }]);
let changed = ctx.ttFlagApplyHide();
ck("hide action: sets sampleHidden", changed === true && ctx.sampleHidden.S1 === true);
ck("hide action: leaves other samples alone", ctx.sampleHidden.S2 !== true);

// Clearing the rule must give visibility back.
verdict(ctx, []);
ctx.ttFlagApplyHide();
ck("clearing the rule restores visibility", ctx.sampleHidden.S1 === false);

// A sample the USER hid must never be un-hidden by the engine.
ctx = makeCtx();
ctx.sampleHidden.S2 = true; // as if the eye icon were clicked
verdict(ctx, [{ source: "meta", field: "total_reads", op: ">", value: "1000", action: "hide" }]);
ctx.ttFlagApplyHide();
verdict(ctx, []);
ctx.ttFlagApplyHide();
ck("a manually hidden sample stays hidden after the rules are cleared", ctx.sampleHidden.S2 === true);

// TT_FLAG_HIDE_ALL promotes every flag to a hide.
ctx = makeCtx();
v = verdict(ctx, [{ source: "meta", field: "total_reads", op: "<", value: "10000", action: "flag" }], {
  hideAll: true,
});
ck("hide-all promotes a flag-only rule to a hide", v.S1.hide === true);
v = verdict(ctx, [{ source: "meta", field: "total_reads", op: "<", value: "10000", action: "flag" }]);
ck("without hide-all a flag-only rule does not hide", v.S1.hide === false);

// ── 5. config load / capture ──────────────────────────────────────────────
ctx = makeCtx();
ctx.ttFlagLoadConfig({
  logic: "all",
  missing_fails: true,
  rules: [
    { source: "meta", field: "total_reads", op: "<", value: 10000, action: "hide" },
    { source: "nope", field: "x", op: "<", value: 1 }, // bad source
    { source: "meta", field: "", op: "<", value: 1 }, // no field
    { source: "meta", field: "aligned_reads", op: "!!", value: 1 }, // bad operator
  ],
});
ck("loadConfig keeps only the valid rule", ctx.__ttGet().rules.length === 1, ctx.__ttGet().rules.length);
ck("loadConfig carries logic and missing_fails", ctx.__ttGet().logic === "all" && ctx.__ttGet().missingFails === true);
const captured = ctx.ttFlagCaptureConfig();
ck(
  "captureConfig round-trips",
  captured.rules.length === 1 && captured.logic === "all" && captured.missing_fails === true,
);
ctx.ttFlagLoadConfig(null);
ck("loadConfig(null) clears the rule set", ctx.__ttGet().rules.length === 0);
ctx.ttFlagLoadConfig(captured);
ck("re-loading a captured config restores it", ctx.__ttGet().rules.length === 1 && ctx.__ttGet().logic === "all");

// Rule ids must stay unique even when a seeded config supplies its own.
ctx.ttFlagLoadConfig({
  rules: [
    { id: "dup", source: "meta", field: "total_reads", op: "<", value: 1 },
    { id: "dup", source: "meta", field: "aligned_reads", op: "<", value: 1 },
  ],
});
ck("rule ids are unique even with duplicate seeds", ctx.__ttGet().rules[0].id !== ctx.__ttGet().rules[1].id);

// ── 6. specimen groups inherit their members' flags ───────────────────────
ctx = makeCtx();
ctx.specimenGroups = () => new Map([["SPEC1", ["S1", "S2"]]]);
verdict(ctx, [{ source: "meta", field: "total_reads", op: "<", value: "10000" }]);
const grp = ctx.ttFlagStateFor("SPEC1");
ck("a merged specimen inherits a flagged member", grp && grp.flagged === true && grp.merged === true);
ck("the inherited reason names the member sample", grp.hits.length === 1 && grp.hits[0].sample === "S1", grp.hits);
ck("a detection row resolves through its Specimen ID", ctx.ttFlagIsFlagged({ "Specimen ID": "S1" }) === true);
ck("an unflagged sample resolves to not-flagged", ctx.ttFlagIsFlagged("S2") === false);

// ── 7. counts, labels and the badge ───────────────────────────────────────
ctx = makeCtx();
verdict(ctx, [{ source: "meta", field: "total_reads", op: "<", value: "10000", action: "hide" }]);
const counts = ctx.ttFlagCounts();
ck(
  "counts: flagged / hidden / total / rules",
  counts.flagged === 1 && counts.hidden === 1 && counts.total === 3 && counts.rules === 1,
  counts,
);
ck(
  "ruleLabel reads as a sentence",
  ctx.ttFlagRuleLabel(ctx.__ttGet().rules[0]) === "Total reads (sample) is less than 10,000",
  ctx.ttFlagRuleLabel(ctx.__ttGet().rules[0]),
);
ck(
  "derived labels carry the TASS cutoff",
  ctx.ttFlagRuleLabel(
    ctx.ttFlagNewRule({ source: "derived", field: "unique_taxids_above_tass", op: "<", value: "2", tass: 95 }),
  ) === "Distinct organisms above TASS 95 is less than 2",
);
ck(
  "data labels carry the aggregation",
  ctx.ttFlagRuleLabel(ctx.ttFlagNewRule({ source: "data", field: "Coverage", agg: "mean", op: "<", value: "1" })) ===
    "mean(Coverage) is less than 1",
);
ck("badge renders for a flagged sample", ctx._flagBadgeHTML("S1").includes("tt-flag-badge"));
ck("badge is empty for a clean sample", ctx._flagBadgeHTML("S2") === "");
ck(
  "reasons list the actual value",
  ctx.ttFlagPlainReasons("S1").includes("actual: 5,000"),
  ctx.ttFlagPlainReasons("S1"),
);

// ── 8. field catalogue ────────────────────────────────────────────────────
const fields = ctx.ttFlagFields();
ck(
  "catalogue exposes all four sources",
  ["meta", "derived", "runmeta", "data"].every((k) => Array.isArray(fields[k])),
);
ck(
  "data fields are numeric detection columns only",
  JSON.stringify(fields.data.map((f) => f.key)) === '["TASS Score","Coverage","# Reads Aligned"]',
  fields.data.map((f) => f.key),
);
ck(
  "runmeta catalogue includes the RUN_META column",
  fields.runmeta.some((f) => f.key === "site"),
);
ck(
  "runmeta catalogue also offers sample_type from SAMPLE_META",
  fields.runmeta.some((f) => f.key === "sample_type"),
);
ck(
  "meta catalogue keeps structural keys out",
  !fields.meta.some((f) => f.key === "weights" || f.key === "best_cutoffs"),
);

// ── 9. host exclusion from the detection counts ───────────────────────────
// Without dehosting, human sits in DATA like any other detection. A
// low-complexity rule ("fewer than 2 organisms") exists precisely to catch a
// sample that is nothing but host, so host must not be one of the organisms it
// counts.
ctx = makeCtx();
const hostRow = (s, tass) => ({
  "Specimen ID": s,
  "Taxonomic ID #": "9606",
  "Detected Organism": "Homo sapiens",
  Species: "Homo sapiens",
  Genus: "Homo",
  "TASS Score": tass,
  Level: "Strain",
  "Passes Threshold": "true",
  "High Consequence": "false",
  "# Reads Aligned": 500000,
  Coverage: 99,
});
// S4: host plus ONE real organism — the classic "nothing grew" mNGS sample.
ctx.DATA.push(hostRow("S4", 99));
ctx.DATA.push({
  "Specimen ID": "S4",
  "Taxonomic ID #": "41",
  "Detected Organism": "Org41",
  "TASS Score": 96,
  Level: "Strain",
  Genus: "G4",
  "Passes Threshold": "true",
  "High Consequence": "false",
  "# Reads Aligned": 40,
  Coverage: 48,
});
ctx.SAMPLE_META.S4 = {
  sample_name: "S4",
  total_reads: 60000,
  aligned_reads: 50000,
};
ctx.ttFlagInvalidate();

const s4rows = () => ctx._ttFlagRowsBySample().get("S4");
ck(
  "host is excluded from unique_taxids by default",
  ctx._ttFlagDerivedValue("unique_taxids", s4rows(), {}) === 1,
  ctx._ttFlagDerivedValue("unique_taxids", s4rows(), {}),
);
v = verdict(ctx, [
  {
    source: "derived",
    field: "unique_taxids_above_tass",
    op: "<",
    value: "2",
    tass: 90,
  },
]);
ck("a host + 1-organism sample still trips the low-complexity rule", flaggedNames(v).includes("S4"), flaggedNames(v));
// Turning the exclusion off is what "count host like anything else" means.
v = verdict(
  ctx,
  [
    {
      source: "derived",
      field: "unique_taxids_above_tass",
      op: "<",
      value: "2",
      tass: 90,
    },
  ],
  {
    exclude: "",
  },
);
ck(
  "with the exclusion cleared, host counts as an organism and S4 passes",
  !flaggedNames(v).includes("S4"),
  flaggedNames(v),
);

// Name matching, for a run whose rows carry no taxid.
ctx.__ttSet({ exclude: "Homo sapiens" });
ctx.ttFlagInvalidate();
ck(
  "an exclusion given as a name matches the organism text",
  ctx._ttFlagDerivedValue("unique_taxids", s4rows(), {}) === 1,
);
// A bare taxid must never match text — "41" is a substring of nothing here, but
// a name-style match on digits would be a silent, very confusing filter.
ctx.__ttSet({ exclude: "9606" });
ctx.ttFlagInvalidate();
ck(
  "aggregated columns skip host too",
  ctx._ttFlagAggValue(s4rows(), "# Reads Aligned", "max") === 40,
  ctx._ttFlagAggValue(s4rows(), "# Reads Aligned", "max"),
);
ck(
  "exclusion parsing accepts commas, spaces and duplicates",
  JSON.stringify(ctx.ttFlagParseExclude("9606, 9606  10090")) === '["9606","10090"]',
  ctx.ttFlagParseExclude("9606, 9606  10090"),
);

// ── 10. view modes: all / hide / only ─────────────────────────────────────
ctx = makeCtx();
const lowReads = [{ source: "meta", field: "total_reads", op: "<", value: "10000", action: "flag" }];

verdict(ctx, lowReads, { view: "all" });
ctx.ttFlagApplyHide();
ck("view 'all': a flag-only rule hides nothing", ctx.sampleHidden.S1 !== true && ctx.sampleHidden.S2 !== true);

v = verdict(ctx, lowReads, { view: "hide" });
ctx.ttFlagApplyHide();
ck("view 'hide': the flagged sample is hidden", ctx.sampleHidden.S1 === true && v.S1.hide === true);
ck("view 'hide': an unflagged sample stays visible", ctx.sampleHidden.S2 !== true);

v = verdict(ctx, lowReads, { view: "only" });
ctx.ttFlagApplyHide();
ck("view 'only': the flagged sample is visible", ctx.sampleHidden.S1 === false || ctx.sampleHidden.S1 === undefined);
ck("view 'only': unflagged samples are hidden", ctx.sampleHidden.S2 === true && ctx.sampleHidden.S3 === true);
ck("view 'only': a flagged sample never reports hide", v.S1.hide === false);
let c10 = ctx.ttFlagCounts();
ck("view 'only': counts report the unflagged ones as hidden", c10.hidden === 2 && c10.onlyActive === true, c10);

// A per-rule "hide" action must not fight the mode: in "only" the flagged
// sample is what the user asked to see.
v = verdict(ctx, [Object.assign({}, lowReads[0], { action: "hide" })], {
  view: "only",
});
ctx.ttFlagApplyHide();
ck("view 'only' overrides a rule's own hide action", ctx.sampleHidden.S1 !== true && v.S1.hide === false);

// Nothing flagged: "only" must not blank the report.
ctx = makeCtx();
verdict(ctx, [{ source: "meta", field: "total_reads", op: "<", value: "1" }], {
  view: "only",
});
ctx.ttFlagApplyHide();
ck(
  "view 'only' with nothing flagged leaves every sample visible",
  ctx.sampleHidden.S1 !== true && ctx.sampleHidden.S2 !== true && ctx.sampleHidden.S3 !== true,
);
ck("counts say 'only' is inert when nothing is flagged", ctx.ttFlagCounts().onlyActive === false);

// Switching back restores what the mode hid.
ctx = makeCtx();
verdict(ctx, lowReads, { view: "only" });
ctx.ttFlagApplyHide();
verdict(ctx, lowReads, { view: "all" });
ctx.ttFlagApplyHide();
ck("leaving 'only' gives the unflagged samples back", ctx.sampleHidden.S2 === false && ctx.sampleHidden.S3 === false);

// ── 11. view / exclusion survive a config round-trip ──────────────────────
ctx = makeCtx();
ctx.ttFlagLoadConfig({
  view: "only",
  exclude_taxids: ["9606", "10090"],
  rules: [],
});
ck(
  "loadConfig reads view and exclude_taxids",
  ctx.__ttGet().view === "only" && JSON.stringify(ctx.__ttGet().exclude) === '["9606","10090"]',
  ctx.__ttGet(),
);
const cap2 = ctx.ttFlagCaptureConfig();
ck(
  "captureConfig writes view, the legacy hide_all and the exclusion",
  cap2.view === "only" && cap2.hide_all === false && JSON.stringify(cap2.exclude_taxids) === '["9606","10090"]',
  cap2,
);
// Sessions saved before the tri-state existed only carry hide_all.
ctx.ttFlagLoadConfig({ hide_all: true, rules: [] });
ck("a legacy hide_all:true config loads as view 'hide'", ctx.__ttGet().view === "hide");
ctx.ttFlagLoadConfig({ hide_all: false, rules: [] });
ck("a legacy hide_all:false config loads as view 'all'", ctx.__ttGet().view === "all");
ctx.ttFlagLoadConfig({ view: "nonsense", rules: [] });
ck("an unknown view falls back to 'all'", ctx.__ttGet().view === "all");
ctx.ttFlagLoadConfig({ rules: [] });
ck(
  "a config with no exclusion key keeps the host default",
  JSON.stringify(ctx.__ttGet().exclude) === '["9606"]',
  ctx.__ttGet().exclude,
);
ctx.ttFlagLoadConfig({ exclude_taxids: "", rules: [] });
ck("an explicitly empty exclusion means count everything", ctx.__ttGet().exclude.length === 0);

// ── 12. classifier (K2) vs aligned reads ─────────────────────────────────
// K2 Reads is what Kraken2/Centrifuge assigned to the organism; # Reads Aligned
// is what the aligner corroborated. A big gap is the classic database-artefact
// signature, and it is invisible to either column on its own.
ctx = makeCtx();
const k2row = (s, tax, k2, aligned) => ({
  "Specimen ID": s,
  "Taxonomic ID #": String(tax),
  "Detected Organism": "Org" + tax,
  "TASS Score": 80,
  Level: "Strain",
  Genus: "G" + String(tax)[0],
  "Passes Threshold": "true",
  "High Consequence": "false",
  "K2 Reads": k2,
  "# Reads Aligned": aligned,
});
// S5: 900 classifier reads, 12 aligned — one organism the aligner does not back.
ctx.DATA.push(k2row("S5", 51, 900, 12));
ctx.DATA.push(k2row("S5", 52, 400, 380));
// S6: both organisms align in line with what the classifier called.
ctx.DATA.push(k2row("S6", 61, 500, 460));
ctx.DATA.push(k2row("S6", 62, 300, 240));
ctx.SAMPLE_META.S5 = { sample_name: "S5", total_reads: 90000, aligned_reads: 400 };
ctx.SAMPLE_META.S6 = { sample_name: "S6", total_reads: 90000, aligned_reads: 700 };
ctx.ttFlagInvalidate();
const rowsOf = (s) => ctx._ttFlagRowsBySample().get(s);

ck(
  "k2_reads_sum adds the classifier column",
  ctx._ttFlagDerivedValue("k2_reads_sum", rowsOf("S5"), {}) === 1300,
  ctx._ttFlagDerivedValue("k2_reads_sum", rowsOf("S5"), {}),
);
ck(
  "aligned_to_k2_ratio is aligned ÷ classified",
  Math.abs(ctx._ttFlagDerivedValue("aligned_to_k2_ratio", rowsOf("S5"), {}) - 392 / 1300) < 1e-9,
  ctx._ttFlagDerivedValue("aligned_to_k2_ratio", rowsOf("S5"), {}),
);
// The fixture's original samples carry no K2 Reads at all: a ratio of 0/0 is
// undefined, and must read as MISSING rather than as a very suspicious zero.
ck(
  "aligned_to_k2_ratio is missing, not 0, when nothing was classified",
  ctx._ttFlagDerivedValue("aligned_to_k2_ratio", rowsOf("S1"), {}) === null,
);
v = verdict(ctx, [{ source: "derived", field: "aligned_to_k2_ratio", op: "<", value: "0.5" }]);
ck("a sample with no classifier reads does not trip a ratio rule", !flaggedNames(v).includes("S1"), flaggedNames(v));
ck("the ratio rule catches the classifier-carried sample", flaggedNames(v).includes("S5"), flaggedNames(v));
ck("the ratio rule leaves the corroborated sample alone", !flaggedNames(v).includes("S6"), flaggedNames(v));

ck(
  "unsupported_k2_organisms counts the unaligned organism only",
  ctx._ttFlagDerivedValue("unsupported_k2_organisms", rowsOf("S5"), { k2min: 50 }) === 1,
  ctx._ttFlagDerivedValue("unsupported_k2_organisms", rowsOf("S5"), { k2min: 50 }),
);
ck(
  "unsupported_k2_organisms is 0 where alignment backs the classifier",
  ctx._ttFlagDerivedValue("unsupported_k2_organisms", rowsOf("S6"), { k2min: 50 }) === 0,
);
// The floor is the point: below it the comparison is noise, so a 900-read
// organism stops counting once the floor is above it.
ck(
  "the K2 floor excludes organisms below it",
  ctx._ttFlagDerivedValue("unsupported_k2_organisms", rowsOf("S5"), { k2min: 1000 }) === 0,
);
v = verdict(ctx, [{ source: "derived", field: "unsupported_k2_organisms", op: ">=", value: "1", k2min: 50 }]);
ck("unsupported_k2_organisms flags the right sample", JSON.stringify(flaggedNames(v)) === '["S5"]', flaggedNames(v));
ck(
  "the rule label names the K2 floor",
  ctx.ttFlagRuleLabel(
    ctx.ttFlagNewRule({ source: "derived", field: "unsupported_k2_organisms", op: ">=", value: "1", k2min: 50 }),
  ) === "Classifier-only organisms (≥ 50 K2 reads) is at least 1",
  ctx.ttFlagRuleLabel(
    ctx.ttFlagNewRule({ source: "derived", field: "unsupported_k2_organisms", op: ">=", value: "1", k2min: 50 }),
  ),
);
// Host exclusion applies here too — human is often the loudest classifier call.
ctx.DATA.push({
  "Specimen ID": "S6",
  "Taxonomic ID #": "9606",
  "Detected Organism": "Homo sapiens",
  "TASS Score": 99,
  Level: "Strain",
  "K2 Reads": 500000,
  "# Reads Aligned": 3,
  "Passes Threshold": "true",
  "High Consequence": "false",
});
ctx.ttFlagInvalidate();
ck(
  "host never counts as a classifier-only organism",
  ctx._ttFlagDerivedValue("unsupported_k2_organisms", rowsOf("S6"), { k2min: 50 }) === 0,
);
ck(
  "host is kept out of the classifier read total",
  ctx._ttFlagDerivedValue("k2_reads_sum", rowsOf("S6"), {}) === 800,
  ctx._ttFlagDerivedValue("k2_reads_sum", rowsOf("S6"), {}),
);
// k2min must survive the config round-trip like tass does.
ctx.ttFlagLoadConfig({
  rules: [{ source: "derived", field: "unsupported_k2_organisms", op: ">=", value: 2, k2min: 250 }],
});
ck("loadConfig keeps k2min", ctx.__ttGet().rules[0].k2min === 250, ctx.__ttGet().rules[0]);
ck("captureConfig writes k2min", ctx.ttFlagCaptureConfig().rules[0].k2min === 250);

console.log(fails ? `\n${fails} check(s) FAILED` : "\nAll checks passed");
process.exit(fails ? 1 : 0);
