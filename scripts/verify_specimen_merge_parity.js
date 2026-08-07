/* Parity + speed check for the specimen-merge optimisation.
 *
 *   node scripts/verify_specimen_merge_parity.js
 *
 * Loads the ORIGINAL 01_global_state.js / 02_utilities.js (from git HEAD) and
 * the CURRENT working-tree versions into two separate VMs, runs both over the
 * same synthetic 164-specimen run under a matrix of filter settings, and
 * asserts the outputs are byte-identical. Then reports the speedup.
 *
 * The optimisation is pure memoisation + removing a redundant O(n^2) rescan,
 * so any difference in output is a bug.
 */
const fs = require("fs");
const path = require("path");
const vm = require("vm");
const { execSync } = require("child_process");

const ROOT = path.join(__dirname, "..");
const SRC = path.join(ROOT, "assets", "src", "js");
const FILES = ["01_global_state.js", "02_utilities.js"];

const headSource = (f) => execSync(`git show HEAD:assets/src/js/${f}`, { cwd: ROOT, maxBuffer: 1 << 28 }).toString();
const workingSource = (f) => fs.readFileSync(path.join(SRC, f), "utf8");

// ── Synthetic run: 164 specimens x 2 libraries ───────────────────────────
// The HEAD implementation is O(groups x rows), so a full 164x2 run takes it
// minutes. Parity is a correctness check, not a scale test — default to a
// smaller but structurally identical run and let the size be raised via env.
//   SPECIMENS=164 ORGS=45 node scripts/verify_specimen_merge_parity.js
function makeRun() {
  const SPECIMENS = parseInt(process.env.SPECIMENS || "25", 10);
  const ORGS = parseInt(process.env.ORGS || "16", 10);
  const LEVELS = ["Strain", "Species", "Genus"];
  const SAMPLE_META = {};
  const DATA = [];
  let seed = 42;
  const rnd = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
  for (let s = 1; s <= SPECIMENS; s++) {
    const spec = "SPEC" + s;
    ["DNA", "RNA"].forEach((mt) => {
      const sample = `${spec}_${mt}`;
      SAMPLE_META[sample] = {
        specimen: spec,
        mol_type: mt,
        total_reads: 1e6 + Math.floor(rnd() * 1e6),
        total_organism_reads: 5e5 + Math.floor(rnd() * 1e5),
      };
      for (let o = 0; o < ORGS; o++) {
        const taxid = 1000 + Math.floor(rnd() * 400);
        const tass = Math.round(rnd() * 100);
        DATA.push({
          "Specimen ID": sample,
          "Detected Organism": "Organism " + taxid,
          "Taxonomic ID #": String(taxid),
          Subkey: "sk" + (taxid % 50),
          Level: LEVELS[o % 3],
          "TASS Score": tass,
          "Species TASS": Math.round(rnd() * 100),
          "Genus TASS": Math.round(rnd() * 100),
          Coverage: Math.round(rnd() * 100),
          "# Reads Aligned": Math.floor(rnd() * 5000),
          "# Unique Reads Aligned": Math.floor(rnd() * 3000),
          "Mean Depth": Math.round(rnd() * 40),
          "Mean MapQ": Math.round(rnd() * 60),
          "Mean BaseQ": Math.round(rnd() * 40),
          "Microbial Category": o % 4 === 0 ? "Bacteria" : o % 4 === 1 ? "Viruses" : "Fungi",
          "High Consequence": o % 17 === 0 ? "true" : "",
          "Passes Threshold": tass >= 50 ? "true" : "",
          Domain: o % 4 === 1 ? "Viruses" : "Bacteria",
          Kingdom: "",
          "Mol Type": mt,
        });
      }
    });
  }
  return { SAMPLE_META, DATA };
}

function makeCtx(sources, run) {
  const CONTROLS = {
    "filter-text": { value: "" },
    "filter-ic": { checked: true },
    "filter-scope": { value: "both" },
    "view-level": { value: "Strain" },
    "filter-min": { value: "0" },
    "filter-hc": { checked: false },
    "filter-pass": { checked: false },
    "filter-mc": { selectedOptions: [] },
    "filter-mt-dna": { checked: true },
    "filter-mt-rna": { checked: true },
    "filter-mt-both": { checked: true },
  };
  const ctx = {
    console,
    CONTROLS,
    document: {
      getElementById: (id) => CONTROLS[id] || null,
      querySelectorAll: () => [],
      querySelector: () => null,
      addEventListener: () => {},
      createElement: () => ({ style: {}, classList: { add() {}, remove() {} }, appendChild() {}, setAttribute() {} }),
      body: { appendChild: () => {}, classList: { add() {}, remove() {} } },
    },
    window: {},
    DATA: run.DATA,
    SAMPLE_META: run.SAMPLE_META,
    sampleHidden: {},
    sampleRescale: {},
    sampleColors: {},
    PALETTE: ["#111"],
    perTypeTass: {},
    ROLLUP_PASS: false,
    watchFilterMode: "all",
    watchlist: new Set(),
    _watchKey: (r) => r["Specimen ID"] + "|" + r["Taxonomic ID #"],
    applyRescale: (r) => r,
    rowPassInfo: (r) => {
      const t = parseFloat(r["TASS Score"]);
      return { strain: t, thr: 50, strainPass: t >= 50, effectivePass: t >= 50 };
    },
    _xsReduceMembers: (vals, method, total) => {
      const s = vals.slice().sort((a, b) => a - b);
      if (method === "min") return s[0];
      if (method === "mean") return s.reduce((a, b) => a + b, 0) / s.length;
      if (method === "median") return s[Math.floor((s.length - 1) / 2)];
      if (method === "detection") return (s.length / (total || 1)) * 100;
      return s[s.length - 1];
    },
    showTip: () => {},
    moveTip: () => {},
    hideTip: () => {},
  };
  vm.createContext(ctx);
  sources.forEach(({ name, code }) => vm.runInContext(code, ctx, { filename: name }));
  return ctx;
}

const run = makeRun();
const OLD = makeCtx(
  FILES.map((f) => ({ name: "HEAD:" + f, code: headSource(f) })),
  run,
);
const NEW = makeCtx(
  FILES.map((f) => ({ name: f, code: workingSource(f) })),
  run,
);

// ── Scenario matrix ──────────────────────────────────────────────────────
const scenarios = [
  { name: "merge off, no filters", merge: false, controls: {} },
  { name: "merge on,  no filters", merge: true, controls: {} },
  { name: "merge on,  pass-only", merge: true, controls: { "filter-pass": { checked: true } } },
  { name: "merge on,  high-consequence", merge: true, controls: { "filter-hc": { checked: true } } },
  { name: "merge on,  text=Organism 12", merge: true, controls: { "filter-text": { value: "Organism 12" } } },
  { name: "merge on,  Species level", merge: true, controls: { "view-level": { value: "Species" } } },
  { name: "merge on,  DNA only", merge: true, controls: { "filter-mt-rna": { checked: false } } },
  { name: "merge on,  agg=mean", merge: true, controls: {}, agg: "mean" },
  { name: "merge on,  agg=median", merge: true, controls: {}, agg: "median" },
  { name: "merge on,  agg=detection", merge: true, controls: {}, agg: "detection" },
  { name: "merge on,  override regroup", merge: true, controls: {}, override: true },
  { name: "merge on,  one sample hidden", merge: true, controls: {}, hide: "SPEC3_RNA" },
];

const DEFAULTS = JSON.parse(
  JSON.stringify({
    "filter-text": { value: "" },
    "filter-ic": { checked: true },
    "filter-scope": { value: "both" },
    "view-level": { value: "Strain" },
    "filter-min": { value: "0" },
    "filter-hc": { checked: false },
    "filter-pass": { checked: false },
    "filter-mc": { selectedOptions: [] },
    "filter-mt-dna": { checked: true },
    "filter-mt-rna": { checked: true },
    "filter-mt-both": { checked: true },
  }),
);

function applyScenario(ctx, sc) {
  Object.keys(DEFAULTS).forEach((k) => {
    ctx.CONTROLS[k] = JSON.parse(JSON.stringify(DEFAULTS[k]));
  });
  Object.entries(sc.controls).forEach(([k, v]) => Object.assign(ctx.CONTROLS[k], v));
  vm.runInContext(`specimenMergeEnabled = ${sc.merge ? "true" : "false"};`, ctx);
  vm.runInContext(`specimenTassAgg = ${JSON.stringify(sc.agg || "max")};`, ctx);
  vm.runInContext(
    `Object.keys(SPECIMEN_OVERRIDE).forEach(k => delete SPECIMEN_OVERRIDE[k]);` +
      (sc.override ? `SPECIMEN_OVERRIDE["SPEC1_DNA"]="MERGED_A";SPECIMEN_OVERRIDE["SPEC2_RNA"]="MERGED_A";` : ""),
    ctx,
  );
  ctx.sampleHidden = {};
  if (sc.hide) ctx.sampleHidden[sc.hide] = true;
  ctx._invalidateFilterCache();
}

// Stable serialisation: sort rows and drop key order differences.
const norm = (rows) =>
  rows
    .map((r) => {
      const o = {};
      Object.keys(r)
        .sort()
        .forEach((k) => {
          const v = r[k];
          o[k] = typeof v === "number" && !Number.isInteger(v) ? Number(v.toFixed(9)) : v;
        });
      return JSON.stringify(o);
    })
    .sort()
    .join("\n");

let failures = 0;
let oldTotal = 0,
  newTotal = 0;

console.log(
  `Parity: HEAD vs working tree — ${Object.keys(run.SAMPLE_META).length / 2} specimens x 2 libraries, ` +
    `${run.DATA.length.toLocaleString()} rows\n`,
);
for (const sc of scenarios) {
  applyScenario(OLD, sc);
  applyScenario(NEW, sc);

  let t = process.hrtime.bigint();
  const a = OLD.filteredData();
  const aAll = OLD.filteredData({ ignoreThreshold: true });
  const aTot = OLD.totalSpecimenCount();
  const aPos = OLD.positiveHitSpecimenCount();
  const oldMs = Number(process.hrtime.bigint() - t) / 1e6;

  t = process.hrtime.bigint();
  const b = NEW.filteredData();
  const bAll = NEW.filteredData({ ignoreThreshold: true });
  const bTot = NEW.totalSpecimenCount();
  const bPos = NEW.positiveHitSpecimenCount();
  const newMs = Number(process.hrtime.bigint() - t) / 1e6;

  oldTotal += oldMs;
  newTotal += newMs;

  const same =
    norm(a) === norm(b) && norm(aAll) === norm(bAll) && aTot === bTot && aPos === bPos && a.length === b.length;
  if (!same) failures++;
  console.log(
    `  ${same ? "OK  " : "DIFF"}  ${sc.name.padEnd(32)} rows ${String(a.length).padStart(6)}  ` +
      `tot ${aTot} pos ${aPos}   ${oldMs.toFixed(0).padStart(6)} ms -> ${newMs.toFixed(0).padStart(5)} ms  ` +
      `(${(oldMs / Math.max(newMs, 0.001)).toFixed(1)}x)`,
  );
  if (!same) {
    console.log(`        rows ${a.length} vs ${b.length}, all ${aAll.length} vs ${bAll.length}`);
    console.log(`        totals ${aTot}/${aPos} vs ${bTot}/${bPos}`);
  }
}

console.log(
  `\nTotal: ${oldTotal.toFixed(0)} ms -> ${newTotal.toFixed(0)} ms  ` +
    `(${(oldTotal / Math.max(newTotal, 0.001)).toFixed(1)}x faster)`,
);
console.log(failures ? `\n${failures} SCENARIO(S) DIFFER` : "\nALL SCENARIOS IDENTICAL");
process.exit(failures ? 1 : 0);
