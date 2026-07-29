/* Benchmark the specimen-merge hot path.
 *
 *   node scripts/bench_specimen_merge.js
 *
 * Loads the real report sources (01_global_state.js + 02_utilities.js) into a
 * VM with a stub DOM, synthesises a run the size of the one that's slow
 * (164 specimens x 2 libraries), and times the calls a single redraw makes.
 * Run it before and after a change to see the effect.
 */
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const SRC = path.join(__dirname, "..", "assets", "src", "js");

// ── Synthetic run ────────────────────────────────────────────────────────
// 164 specimens, 2 libraries each (DNA + RNA) = 328 samples.
const SPECIMENS = 164;
const ORGS_PER_SAMPLE = 45; // ~14.7k detection rows, typical for this report
const LEVELS = ["Strain", "Species", "Genus"];

const SAMPLE_META = {};
const DATA = [];
let seed = 42;
const rnd = () => ((seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff);

for (let s = 1; s <= SPECIMENS; s++) {
  const spec = "SPEC" + s;
  ["DNA", "RNA"].forEach((mt) => {
    const sample = `${spec}_${mt}`;
    SAMPLE_META[sample] = { specimen: spec, mol_type: mt };
    for (let o = 0; o < ORGS_PER_SAMPLE; o++) {
      const taxid = 1000 + Math.floor(rnd() * 400);
      const tass = Math.round(rnd() * 100);
      DATA.push({
        "Specimen ID": sample,
        "Detected Organism": "Organism " + taxid,
        "Taxonomic ID #": String(taxid),
        Subkey: "sk" + (taxid % 50),
        Level: LEVELS[o % 3],
        "TASS Score": tass,
        Coverage: Math.round(rnd() * 100),
        "# Reads Aligned": Math.floor(rnd() * 5000),
        "# Unique Reads Aligned": Math.floor(rnd() * 3000),
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

// ── Stub DOM: only the controls filteredData() reads ─────────────────────
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
const document = {
  getElementById: (id) => CONTROLS[id] || null,
  querySelectorAll: () => [],
  querySelector: () => null,
  addEventListener: () => {},
  createElement: () => ({ style: {}, classList: { add() {}, remove() {} }, appendChild() {}, setAttribute() {} }),
  body: { appendChild: () => {}, classList: { add() {}, remove() {} } },
};

const ctx = {
  console,
  document,
  window: {},
  DATA,
  SAMPLE_META,
  sampleHidden: {},
  sampleRescale: {},
  sampleColors: {},
  PALETTE: ["#111"],
  perTypeTass: {},
  ROLLUP_PASS: false,
  watchFilterMode: "all",
  watchlist: new Set(),
  TASS_THRESHOLD: 50,
  _watchKey: (r) => r["Specimen ID"] + "|" + r["Taxonomic ID #"],
  applyRescale: (r) => r,
  rowPassInfo: (r) => {
    const t = parseFloat(r["TASS Score"]);
    return { strain: t, thr: 50, strainPass: t >= 50, effectivePass: t >= 50 };
  },
  showTip: () => {},
  moveTip: () => {},
  hideTip: () => {},
};
vm.createContext(ctx);
for (const f of ["01_global_state.js", "02_utilities.js"]) {
  vm.runInContext(fs.readFileSync(path.join(SRC, f), "utf8"), ctx, { filename: f });
}

// ── Timing helpers ───────────────────────────────────────────────────────
const time = (label, n, fn) => {
  fn(); // warm
  const t0 = process.hrtime.bigint();
  for (let i = 0; i < n; i++) fn();
  const ms = Number(process.hrtime.bigint() - t0) / 1e6;
  console.log(`  ${label.padEnd(46)} ${ms.toFixed(1).padStart(9)} ms  (${n}x, ${(ms / n).toFixed(2)} ms each)`);
  return ms;
};

console.log(`Synthetic run: ${SPECIMENS} specimens x 2 libraries, ${DATA.length.toLocaleString()} rows\n`);

// `let specimenMergeEnabled` is a lexical binding inside the VM's global
// scope — assigning ctx.specimenMergeEnabled would create a shadowed property
// the sources never see. Flip it with an in-context assignment instead.
vm.runInContext("specimenMergeEnabled = true;", ctx);
ctx._invalidateFilterCache();
console.log("merge actually on:", vm.runInContext("specimenMergeEnabled", ctx), "\n");

console.log("MERGE ON");
const total = {};
total.groups = time("specimenGroups()", 200, () => ctx.specimenGroups());
total.hasGrouping = time("hasSpecimenGrouping()", 200, () => ctx.hasSpecimenGrouping());
total.specimenOf = time("specimenOf() x DATA.length", 20, () => {
  for (let i = 0; i < DATA.length; i++) ctx.specimenOf(DATA[i]["Specimen ID"]);
});
total.badge = time("_mergedSampleBadgeHTML() x 200 rows", 1, () => {
  for (let i = 0; i < 200; i++) ctx._mergedSampleBadgeHTML(DATA[i]["Specimen ID"]);
});
total.totalCount = time("totalSpecimenCount()", 50, () => ctx.totalSpecimenCount());
total.posCount = time("positiveHitSpecimenCount()", 50, () => ctx.positiveHitSpecimenCount());
total.fdCold = time("filteredData() cold (cache cleared)", 5, () => {
  ctx._invalidateFilterCache();
  ctx.filteredData();
});
total.fdWarm = time("filteredData() warm (cached)", 200, () => ctx.filteredData());
total.fdIgnore = time("filteredData({ignoreThreshold:true})", 5, () => ctx.filteredData({ ignoreThreshold: true }));

// A cross-sample redraw: ~10 filteredData() calls + the un-thresholded pass +
// the two specimen counts, all after an invalidation (what a merge toggle does).
const redraw = () => {
  ctx._invalidateFilterCache();
  for (let i = 0; i < 10; i++) ctx.filteredData();
  ctx.filteredData({ ignoreThreshold: true });
  ctx.totalSpecimenCount();
  ctx.positiveHitSpecimenCount();
  for (let i = 0; i < 200; i++) ctx._mergedSampleBadgeHTML(DATA[i]["Specimen ID"]);
};
console.log("");
time("FULL REDRAW after merge toggle", 1, redraw);
