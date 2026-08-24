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
      if (o.hideAll !== undefined) TT_FLAG_HIDE_ALL = o.hideAll;
    }
    function __ttGet() {
      return { rules: TT_FLAG_RULES, logic: TT_FLAG_LOGIC, enabled: TT_FLAG_ENABLED,
               missingFails: TT_FLAG_MISSING_FAILS, hideAll: TT_FLAG_HIDE_ALL };
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
    hideAll: !!o.hideAll,
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
ck("meta: total_reads < 10000 catches only the shallow sample", JSON.stringify(flaggedNames(v)) === '["S1"]', flaggedNames(v));
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
ck("data: mean(Coverage) selects a different sample than max", JSON.stringify(flaggedNames(v)) === '["S2"]', flaggedNames(v));

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
v = verdict(ctx, [{ source: "meta", field: "total_reads", op: "<", value: "10000", action: "flag" }], { hideAll: true });
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
ck("captureConfig round-trips", captured.rules.length === 1 && captured.logic === "all" && captured.missing_fails === true);
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
ck("reasons list the actual value", ctx.ttFlagPlainReasons("S1").includes("actual: 5,000"), ctx.ttFlagPlainReasons("S1"));

// ── 8. field catalogue ────────────────────────────────────────────────────
const fields = ctx.ttFlagFields();
ck("catalogue exposes all four sources", ["meta", "derived", "runmeta", "data"].every((k) => Array.isArray(fields[k])));
ck("data fields are numeric detection columns only", JSON.stringify(fields.data.map((f) => f.key)) ===
  '["TASS Score","Coverage","# Reads Aligned"]', fields.data.map((f) => f.key));
ck("runmeta catalogue includes the RUN_META column", fields.runmeta.some((f) => f.key === "site"));
ck("runmeta catalogue also offers sample_type from SAMPLE_META", fields.runmeta.some((f) => f.key === "sample_type"));
ck("meta catalogue keeps structural keys out", !fields.meta.some((f) => f.key === "weights" || f.key === "best_cutoffs"));

console.log(fails ? `\n${fails} check(s) FAILED` : "\nAll checks passed");
process.exit(fails ? 1 : 0);
