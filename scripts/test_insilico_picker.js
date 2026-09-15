/* In-Silico tab series picker.
 *
 *   node scripts/test_insilico_picker.js
 *
 * The tab renders one group per parent sample, each a screen and a half of
 * tables and charts. At 20+ samples that is a scroll marathon, so the groups
 * are chosen from a searchable dropdown (the shared ttMultiSelect) instead of
 * all being poured onto the page.
 *
 * These checks pin down the behaviour that makes that safe to rely on:
 *   - a large suite renders ONE group by default, not all of them
 *   - the dropdown lists every series and filters as you type
 *   - ticking a second series renders it alongside the first
 *   - All / None do what they say
 *   - a small suite (<= 3 series) still shows everything, and a single series
 *     spends no control at all
 *   - selection survives a redraw, and a stale key from a previously loaded
 *     file cannot leave the tab empty
 *
 * 43_multiselect.js and 45_tab_insilico.js are loaded as-is from source, so
 * edits to either are exercised here automatically.
 */
const fs = require("fs");
const path = require("path");
const { JSDOM } = require("jsdom");

const SRC = path.join(__dirname, "..", "assets", "src", "js");

let failures = 0;
function check(label, cond, extra) {
  if (cond) {
    console.log("  ok   " + label);
  } else {
    failures++;
    console.log("  FAIL " + label + (extra ? "\n       " + extra : ""));
  }
}

// The fixture is cloned from the real committed payload (assets/pages.js) rather
// than hand-built: the group renderer reads dozens of per-dataset fields, and a
// synthetic stub drifts out of date the moment one is added.
function loadRealSuite() {
  const js = fs.readFileSync(path.join(__dirname, "..", "assets", "pages.js"), "utf8");
  const sandbox = { window: {}, globalThis: {} };
  const fn = new Function("window", "globalThis", js + "\nreturn window.HEATMAP_BOOT;");
  const boot = fn(sandbox.window, sandbox.globalThis);
  const suite = boot && boot.insilico_suite;
  if (!suite || !(suite.groups || []).length) {
    console.error("assets/pages.js carries no in-silico suite — cannot run this test.");
    process.exit(1);
  }
  return suite;
}

const REAL_SUITE = loadRealSuite();

// n series, each a deep copy of a real group under a distinct parent name.
function makeSuite(n) {
  const tmpl = REAL_SUITE.groups[0];
  const groups = [];
  for (let i = 0; i < n; i++) {
    const g = JSON.parse(JSON.stringify(tmpl));
    g.parent = "Sample_" + String(i + 1).padStart(2, "0");
    (g.datasets || []).forEach((d) => {
      if (d.id) d.id = String(d.id).replace(/^[^_]+(?:_[A-Z])?/, g.parent);
    });
    groups.push(g);
  }
  return Object.assign({}, REAL_SUITE, { enabled: true, groups });
}

function boot(suite) {
  const dom = new JSDOM(
    `<!doctype html><html><body>
       <div id="insilico-params-panel"></div>
       <div id="insilico-groups"></div>
       <div id="insilico-empty" style="display:none"></div>
       <button id="insilico-tab-btn"></button>
     </body></html>`,
    { pretendToBeVisual: true },
  );
  const { window } = dom;
  global.window = window;
  global.document = window.document;
  global.navigator = window.navigator;
  window.HEATMAP_BOOT = { insilico_suite: suite, has_insilico_suite: true };

  for (const f of ["43_multiselect.js", "45_tab_insilico.js"]) {
    const code = fs.readFileSync(path.join(SRC, f), "utf8");
    window.eval(code);
  }
  window.drawInsilico();
  return window;
}

const groupCount = (w) => w.document.getElementById("insilico-groups").querySelectorAll(":scope > div").length;
const picker = (w) => w.document.getElementById("insilico-group-picker");
const options = (w) => Array.from(w.document.querySelectorAll("#insilico-series-ms .tt-ms-opt"));

// ── a large suite: one series by default, everything listed ────────────────
console.log("large suite (24 series)");
let w = boot(makeSuite(24));
check("only one group rendered by default", groupCount(w) === 1, "got " + groupCount(w));
check("picker is present", !!picker(w) && picker(w).style.display !== "none");

const btn = w.document.querySelector("#insilico-series-ms .tt-ms-btn");
check("collapsed button shows the selected sample", /Sample_01/.test(btn.textContent), btn.textContent);

btn.dispatchEvent(new w.MouseEvent("click", { bubbles: true }));
check("all 24 series are listed", options(w).length === 24, "got " + options(w).length);

// ── search ─────────────────────────────────────────────────────────────────
const search = w.document.querySelector("#insilico-series-ms .tt-ms-search");
check("the list has a search box", !!search);
search.value = "Sample_1"; // matches Sample_10 … Sample_19
search.dispatchEvent(new w.window.Event("input", { bubbles: true }));
const filtered = options(w).map((o) => o.textContent.trim());
check(
  "typing filters the list (Sample_10 … Sample_19 = 10)",
  filtered.length === 10,
  "got " + filtered.length + ": " + filtered.slice(0, 4).join(", "),
);
check(
  "every filtered row matches the query",
  filtered.every((t) => t.includes("Sample_1")),
  filtered.join(", "),
);

// ── ticking a second series ────────────────────────────────────────────────
const cb = w.document.querySelector('#insilico-series-ms .tt-ms-cb[data-key*="Sample_12"]');
check("the filtered row has a checkbox", !!cb);
cb.checked = true;
cb.dispatchEvent(new w.window.Event("change", { bubbles: true }));
check("a second series renders alongside the first", groupCount(w) === 2, "got " + groupCount(w));

// ── All / None ─────────────────────────────────────────────────────────────
// The panel re-renders on every change, so these must be re-queried, never cached.
check(
  "All / None buttons exist",
  !!w.document.querySelector("#insilico-series-ms .tt-ms-all") &&
    !!w.document.querySelector("#insilico-series-ms .tt-ms-none"),
);
w.document.querySelector("#insilico-series-ms .tt-ms-all").dispatchEvent(new w.MouseEvent("click", { bubbles: true }));
check("All renders every series", groupCount(w) === 24, "got " + groupCount(w));
w.document.querySelector("#insilico-series-ms .tt-ms-none").dispatchEvent(new w.MouseEvent("click", { bubbles: true }));
check(
  "None renders no groups but explains itself",
  groupCount(w) === 1 && /No series selected/.test(w.document.getElementById("insilico-groups").textContent),
);

// ── selection survives a redraw ────────────────────────────────────────────
console.log("redraw");
w = boot(makeSuite(24));
const cb2 = w.document.querySelector('#insilico-series-ms .tt-ms-cb[data-key*="Sample_07"]');
cb2.checked = true;
cb2.dispatchEvent(new w.window.Event("change", { bubbles: true }));
const before = groupCount(w);
w.drawInsilico();
check("the same groups are shown after a redraw", groupCount(w) === before, before + " -> " + groupCount(w));

// ── a stale selection from a previously loaded file ───────────────────────
// Dropping a new dataset on the report re-renders the tab with different group
// keys. The old selection must not survive as an empty tab.
console.log("stale selection");
w = boot(makeSuite(24));
const stale = w.document.querySelector('#insilico-series-ms .tt-ms-cb[data-key*="Sample_09"]');
stale.checked = true;
stale.dispatchEvent(new w.window.Event("change", { bubbles: true }));
const renamed = makeSuite(24);
renamed.groups.forEach((g, i) => (g.parent = "Other_" + String(i + 1).padStart(2, "0")));
w.HEATMAP_BOOT.insilico_suite = renamed;
w.drawInsilico();
check("a stale selection falls back to the default, not an empty tab", groupCount(w) === 1, "got " + groupCount(w));
check(
  "the new series are listed",
  /Other_01/.test(w.document.querySelector("#insilico-series-ms .tt-ms-btn").textContent),
);

// ── small suites ───────────────────────────────────────────────────────────
console.log("small suites");
w = boot(makeSuite(3));
check("3 series: all shown by default", groupCount(w) === 3, "got " + groupCount(w));
check("3 series: picker still offered", picker(w) && picker(w).style.display !== "none");

w = boot(makeSuite(1));
check("1 series: rendered", groupCount(w) === 1);
check("1 series: no picker (nothing to choose)", picker(w) === null || picker(w).style.display === "none");

// ── an empty suite hides the tab content ───────────────────────────────────
console.log("empty suite");
w = boot({ enabled: false, groups: [] });
check("no groups rendered", groupCount(w) === 0);
check("empty-state block is shown", w.document.getElementById("insilico-empty").style.display === "");

console.log(failures ? `\n${failures} check(s) FAILED` : "\nAll checks passed.");
process.exit(failures ? 1 : 0);
