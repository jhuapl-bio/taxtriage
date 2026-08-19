/* Headless test for docs/javascripts/pathogen_table.js.
 *
 * Boots the real widget in jsdom against the real pathogen CSV and asserts
 * on rendered DOM — search, facet AND/OR semantics, excluding-self counts,
 * taxonomy drilldown, sorting, pagination, detail drawer and CSV export.
 *
 * Run:  node scripts/test_pathogen_table.js
 */

const fs = require("fs");
const path = require("path");
const { JSDOM } = require("jsdom");

const ROOT = path.resolve(__dirname, "..");
// The page fetches this from raw.githubusercontent at runtime; for tests we
// use a local checkout of the pipeline repo so the suite stays offline.
const CSV = process.env.PATHOGEN_CSV || path.resolve(ROOT, "..", "taxtriage", "assets", "pathogen_sheet.csv");
const SCRIPT = path.join(ROOT, "docs", "javascripts", "pathogen_table.js");

let failures = 0;
let checks = 0;

function ok(label, cond, extra) {
  checks++;
  if (cond) {
    console.log(`  ✓ ${label}`);
  } else {
    failures++;
    console.log(`  ✗ ${label}${extra ? "  → " + extra : ""}`);
  }
}

function eq(label, actual, expected) {
  ok(label, actual === expected, `got ${JSON.stringify(actual)}, want ${JSON.stringify(expected)}`);
}

const csvText = fs.readFileSync(CSV, "utf8");

// Independent CSV reader for ground truth — deliberately not the one under
// test, so a parser bug cannot make the assertions agree with themselves.
function refParse(text) {
  const out = [];
  let row = [],
    field = "",
    q = false;
  text = text.replace(/^\uFEFF/, "").replace(/\r\n?/g, "\n");
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (q) {
      if (c === '"' && text[i + 1] === '"') {
        field += '"';
        i++;
      } else if (c === '"') q = false;
      else field += c;
    } else if (c === '"') q = true;
    else if (c === ",") {
      row.push(field);
      field = "";
    } else if (c === "\n") {
      row.push(field);
      field = "";
      if (row.length > 1 || row[0] !== "") out.push(row);
      row = [];
    } else field += c;
  }
  if (field !== "" || row.length) {
    row.push(field);
    if (row.length > 1 || row[0] !== "") out.push(row);
  }
  return { cols: out[0].map((h) => h.trim()), rows: out.slice(1) };
}

const parsed = refParse(csvText);
const col = (n) => parsed.cols.indexOf(n);
const payload = {
  cols: parsed.cols,
  rows: parsed.rows
    .map((cells) =>
      parsed.cols.map((c, i) => {
        const v = (cells[i] || "").trim();
        if (["pathogenic_sites", "commensal_sites", "alternative_names"].includes(c))
          return v
            ? v
                .split(",")
                .map((s) => s.trim())
                .filter(Boolean)
            : [];
        if (c === "high_consequence") return v.toUpperCase() === "TRUE";
        if (c === "status") return { estbalished: "established", etsablished: "established" }[v] || v;
        return v;
      }),
    )
    .filter((r) => r[parsed.cols.indexOf("name")]),
};

/* ── ground truth computed independently from the raw payload ─────────── */

const truth = {
  total: payload.rows.length,
  primary: payload.rows.filter((r) => r[col("general_classification")] === "primary").length,
  highConsequence: payload.rows.filter((r) => r[col("high_consequence")] === true).length,
  bloodSite: payload.rows.filter((r) => r[col("pathogenic_sites")].includes("blood")).length,
  rnaPrimary: payload.rows.filter((r) => r[col("mol_type")] === "rna" && r[col("general_classification")] === "primary")
    .length,
  primaryOrCommensal: payload.rows.filter((r) => ["primary", "commensal"].includes(r[col("general_classification")]))
    .length,
  established: payload.rows.filter((r) => r[col("status")] === "established").length,
};

async function mount(text, opts) {
  opts = opts || {};
  const dom = new JSDOM(`<!doctype html><html><body><div class="pt-app" id="pathogen-app"></div></body></html>`, {
    url: "https://example.org/pathogen-sheet/",
    runScripts: "outside-only",
  });
  const { window } = dom;
  const fetched = [];
  if (opts.docsRef) window.TAXTRIAGE_DOCS = opts.docsRef;
  window.fetch = (url) => {
    fetched.push(url);
    return opts.fail
      ? Promise.resolve({ ok: false, status: 404, text: () => Promise.resolve("") })
      : Promise.resolve({ ok: true, status: 200, text: () => Promise.resolve(text) });
  };
  window.__fetched = fetched;
  window.URL.createObjectURL = () => "blob:stub";
  window.URL.revokeObjectURL = () => {};
  window.HTMLAnchorElement.prototype.click = function () {};
  window.eval(fs.readFileSync(SCRIPT, "utf8"));
  window.document.dispatchEvent(new window.Event("DOMContentLoaded"));
  await new Promise((r) => setTimeout(r, 80));
  return window;
}

async function main() {
  const dom = new JSDOM(`<!doctype html><html><body><div class="pt-app" id="pathogen-app"></div></body></html>`, {
    url: "https://example.org/pathogen-sheet/",
    runScripts: "outside-only",
  });
  const { window } = dom;

  // Minimal shims for the browser APIs the widget touches.
  window.fetch = () => Promise.resolve({ ok: true, status: 200, text: () => Promise.resolve(csvText) });
  const downloads = [];
  window.URL.createObjectURL = (blob) => {
    downloads.push(blob);
    return "blob:stub";
  };
  window.URL.revokeObjectURL = () => {};
  window.HTMLAnchorElement.prototype.click = function () {};

  window.eval(fs.readFileSync(SCRIPT, "utf8"));
  window.document.dispatchEvent(new window.Event("DOMContentLoaded"));
  await new Promise((r) => setTimeout(r, 60)); // let the fetch promise settle

  const doc = window.document;
  const $ = (s) => doc.querySelector(s);
  const $$ = (s) => Array.from(doc.querySelectorAll(s));

  const countText = () => $(".pt-count").textContent;
  const shownTotal = () =>
    parseInt(
      countText()
        .replace(/,/g, "")
        .match(/(\d+) of/)[1],
      10,
    );
  const bodyRows = () => $$(".pt-table tbody tr[data-name]").length;

  const fire = (el, type) => el.dispatchEvent(new window.Event(type, { bubbles: true }));

  function facetOption(key, value) {
    return $$(`.pt-facet[data-key="${key}"] input[type=checkbox]`).find((i) => i.value === value);
  }
  function optionCount(key, value) {
    const input = facetOption(key, value);
    if (!input) return null;
    return parseInt(input.closest(".pt-opt").querySelector(".pt-opt-n").textContent, 10);
  }
  function toggle(key, value) {
    const input = facetOption(key, value);
    if (!input) throw new Error(`no option ${key}=${value}`);
    input.checked = !input.checked;
    fire(input, "change");
  }
  function reset() {
    $(".pt-reset").click();
  }
  function search(q) {
    $(".pt-search input").value = q;
    fire($(".pt-search input"), "input");
    return new Promise((r) => setTimeout(r, 200)); // debounce is 140ms
  }

  console.log("\nrender");
  ok("widget mounted", !!$(".pt-table"), "no table found");
  eq("all rows counted", shownTotal(), truth.total);
  eq("first page is 50 rows", bodyRows(), 50);
  eq("column headers", $$(".pt-table thead th").length, 11);

  console.log("\nfacet counts (unfiltered)");
  eq("classification: primary", optionCount("general_classification", "primary"), truth.primary);
  eq("high_consequence: TRUE", optionCount("high_consequence", "TRUE"), truth.highConsequence);
  eq("pathogenic_sites: blood", optionCount("pathogenic_sites", "blood"), truth.bloodSite);
  eq("status: established (typos folded in)", optionCount("status", "established"), truth.established);
  ok("typo 'estbalished' is gone", !facetOption("status", "estbalished"));

  console.log("\nsingle facet");
  toggle("general_classification", "primary");
  eq("filters to primary", shownTotal(), truth.primary);
  ok("chip appears", $$(".pt-chip").length === 1, `${$$(".pt-chip").length} chips`);

  console.log("\nOR within a group");
  toggle("general_classification", "commensal");
  eq("primary OR commensal", shownTotal(), truth.primaryOrCommensal);
  eq(
    "excluding-self count unchanged by own selection",
    optionCount("general_classification", "primary"),
    truth.primary,
  );
  reset();

  console.log("\nAND across groups");
  toggle("mol_type", "rna");
  toggle("general_classification", "primary");
  eq("rna AND primary", shownTotal(), truth.rnaPrimary);
  eq(
    "cross-facet count narrows",
    optionCount("mol_type", "rna"),
    payload.rows.filter((r) => r[col("general_classification")] === "primary" && r[col("mol_type")] === "rna").length,
  );
  reset();
  eq("reset restores everything", shownTotal(), truth.total);

  console.log("\nmulti-value facet");
  toggle("pathogenic_sites", "blood");
  eq("blood site membership", shownTotal(), truth.bloodSite);
  reset();

  console.log("\nsearch");
  await search("mycobacterium");
  const mycoExpected = payload.rows.filter((r) =>
    [
      "name",
      "alternative_names",
      "taxid",
      "genus",
      "family",
      "order",
      "class",
      "phylum",
      "kingdom",
      "pathology",
      "assembly_accession",
    ]
      .map((c) => {
        const v = r[col(c)];
        return Array.isArray(v) ? v.join(", ") : String(v == null ? "" : v);
      })
      .join(" ")
      .toLowerCase()
      .includes("mycobacterium"),
  ).length;
  eq("free-text search", shownTotal(), mycoExpected);
  ok("matches are highlighted", $$(".pt-mark").length > 0);

  console.log("\nsearch + facet compose");
  toggle("general_classification", "primary");
  const both = payload.rows.filter((r) => {
    const blob = [
      "name",
      "alternative_names",
      "taxid",
      "genus",
      "family",
      "order",
      "class",
      "phylum",
      "kingdom",
      "pathology",
      "assembly_accession",
    ]
      .map((c) => {
        const v = r[col(c)];
        return Array.isArray(v) ? v.join(", ") : String(v == null ? "" : v);
      })
      .join(" ")
      .toLowerCase();
    return blob.includes("mycobacterium") && r[col("general_classification")] === "primary";
  }).length;
  eq("search AND facet", shownTotal(), both);
  reset();
  await search("");

  console.log("\ntaxonomy drilldown");
  const phylumSel = $('.pt-tax select[data-rank="phylum"]');
  ok("phylum enabled", !phylumSel.disabled);
  ok("class starts disabled", $('.pt-tax select[data-rank="class"]').disabled);
  phylumSel.value = "Bacillota";
  fire(phylumSel, "change");
  eq("phylum filter applied", shownTotal(), payload.rows.filter((r) => r[col("phylum")] === "Bacillota").length);
  ok("class now enabled", !$('.pt-tax select[data-rank="class"]').disabled);
  const classSel = $('.pt-tax select[data-rank="class"]');
  const classOpts = Array.from(classSel.options)
    .slice(1)
    .map((o) => o.value);
  ok(
    "class options constrained to phylum",
    classOpts.every((c) => payload.rows.some((r) => r[col("phylum")] === "Bacillota" && r[col("class")] === c)),
    classOpts.join(","),
  );
  classSel.value = classOpts[0];
  fire(classSel, "change");
  eq(
    "class narrows further",
    shownTotal(),
    payload.rows.filter((r) => r[col("phylum")] === "Bacillota" && r[col("class")] === classOpts[0]).length,
  );
  // Changing the parent must clear the child. Re-query first: every filter
  // change re-renders the facet rail, so earlier references are detached.
  const phylumLive = $('.pt-tax select[data-rank="phylum"]');
  phylumLive.value = "";
  fire(phylumLive, "change");
  eq("clearing parent clears child", shownTotal(), truth.total);
  eq("child select reset", $('.pt-tax select[data-rank="class"]').value, "");
  ok("child select disabled again", $('.pt-tax select[data-rank="class"]').disabled);
  reset();

  console.log("\nsorting");
  const firstBefore = $(".pt-table tbody tr .pt-name").textContent;
  const nameTh = $$(".pt-table thead th").find((t) => t.dataset.key === "name");
  nameTh.click();
  const firstDesc = $(".pt-table tbody tr .pt-name").textContent;
  ok("clicking sorted column reverses it", firstBefore !== firstDesc, `${firstBefore} vs ${firstDesc}`);
  eq("aria-sort reflects direction", nameTh.getAttribute("aria-sort"), "descending");
  const taxTh = $$(".pt-table thead th").find((t) => t.dataset.key === "taxid");
  taxTh.click();
  const ids = $$(".pt-table tbody tr").map((tr) => Number(tr.children[1].textContent.trim()));
  ok(
    "taxid sorts numerically",
    ids.every((v, i) => i === 0 || ids[i - 1] <= v),
    ids.slice(0, 5).join(","),
  );

  console.log("\npagination");
  reset();
  const per = $(".pt-per");
  per.value = "25";
  fire(per, "change");
  eq("page size honoured", bodyRows(), 25);
  const firstPage = $(".pt-table tbody tr").dataset.name;
  $(".pt-next").click();
  ok("next page changes rows", $(".pt-table tbody tr").dataset.name !== firstPage);
  ok("prev is enabled on page 2", !$(".pt-prev").disabled);
  $(".pt-prev").click();
  eq("prev returns to page 1", $(".pt-table tbody tr").dataset.name, firstPage);

  console.log("\ndetail drawer");
  $(".pt-table tbody tr").dispatchEvent(new window.Event("click", { bubbles: true }));
  ok("drawer opens", $(".pt-drawer").classList.contains("is-open"));
  ok("drawer has a title", !!$(".pt-drawer h3").textContent.trim());
  ok("drawer links tax ID to NCBI", !!$('.pt-drawer a[href*="ncbi.nlm.nih.gov"]'));
  $(".pt-drawer-close").dispatchEvent(new window.Event("click", { bubbles: true }));
  ok("drawer closes", !$(".pt-drawer").classList.contains("is-open"));

  console.log("\nCSV export");
  reset();
  toggle("general_classification", "primary");
  $(".pt-export").click();
  ok("a blob was produced", downloads.length === 1, `${downloads.length} downloads`);
  const csv = await downloads[0].text();
  const lines = csv.split("\n");
  eq("export row count matches filter", lines.length - 1, truth.primary);
  ok("export has the full header", lines[0].startsWith("name,taxid,general_classification"), lines[0].slice(0, 60));
  ok("quoted fields are escaped", csv.includes('""') || !csv.includes('"'), "unbalanced quoting");

  console.log("\nempty state");
  reset();
  await search("zzzzzznotarealorganism");
  eq("no matches", shownTotal(), 0);
  ok("empty message shown", $(".pt-table tbody .pt-empty") !== null);

  /* ── schema tolerance ──────────────────────────────────────────────────
   * The live sheet may not carry every column. A missing column must read as
   * blank and drop its facet, not break the page.
   */
  console.log("\nschema tolerance: CSV with no `status` column");
  {
    const lines = csvText.replace(/\r\n?/g, "\n").split("\n");
    const header = refParse(csvText).cols;
    const si = header.indexOf("status");
    // Rebuild the CSV without the status column, quoting preserved.
    const drop = (row) => row.filter((_, i) => i !== si);
    const q = (v) => (/[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v);
    const reparsed = refParse(csvText);
    const noStatus =
      drop(reparsed.cols).map(q).join(",") + "\n" + reparsed.rows.map((r) => drop(r).map(q).join(",")).join("\n");

    const w = await mount(noStatus);
    const d = w.document;
    const total = parseInt(
      d
        .querySelector(".pt-count")
        .textContent.replace(/,/g, "")
        .match(/(\d+) of/)[1],
      10,
    );
    eq("all rows still load without `status`", total, truth.total);
    ok("status facet is hidden", !d.querySelector('.pt-facet[data-key="status"]'), "facet still rendered");
    ok(
      "other facets still render",
      !!d.querySelector('.pt-facet[data-key="general_classification"]'),
      "classification facet missing",
    );
    const statusIdx = 3; // Organism, Tax ID, Class., Status
    const firstRowCells = d.querySelectorAll(".pt-table tbody tr td");
    eq("status cell renders blank", firstRowCells[statusIdx].textContent.trim(), "");
    ok("table still has all 11 columns", d.querySelectorAll(".pt-table thead th").length === 11);
  }

  console.log("\nCSV parsing edge cases");
  {
    const tricky =
      "name,taxid,general_classification,pathogenic_sites,reference\n" +
      '"Smith, John et al.",99,primary,"blood, gut","He said ""hello"", loudly"\n' +
      'Plain organism,100,commensal,skin,"line one\nline two"\n';
    const w = await mount(tricky);
    const d = w.document;
    const total = parseInt(
      d
        .querySelector(".pt-count")
        .textContent.replace(/,/g, "")
        .match(/(\d+) of/)[1],
      10,
    );
    eq("two records despite embedded newline", total, 2);
    const names = Array.from(d.querySelectorAll(".pt-table tbody tr .pt-name")).map((e) => e.textContent.trim());
    ok("comma inside quotes kept in one field", names.includes("Smith, John et al."), names.join(" | "));
    const sitesCell = d.querySelectorAll(".pt-table tbody tr")[names.indexOf("Smith, John et al.")].children[5];
    eq("multi-value split on comma", sitesCell.textContent.trim(), "blood, gut");
  }

  console.log("\nfetch failure");
  {
    const w = await mount("", { fail: true });
    const html = w.document.getElementById("pathogen-app").innerHTML;
    ok("shows an error message", /Could not load/.test(html));
    ok("offers a direct GitHub link", /github\.com\/jhuapl-bio\/taxtriage/.test(html));
  }

  /* ── version pinning ───────────────────────────────────────────────────
   * Versioned docs must read the pathogen sheet from the ref that shipped with
   * that release, not from main. scripts/write_docs_ref.py sets the global.
   */
  console.log("\nversion-pinned data ref");
  {
    const w = await mount(csvText, { docsRef: { ref: "v3.3.8", version: "3.3" } });
    const url = w.__fetched[0];
    ok(
      "fetches from the pinned tag",
      url === "https://raw.githubusercontent.com/jhuapl-bio/taxtriage/v3.3.8/assets/pathogen_sheet.csv",
      url,
    );
    const total = parseInt(
      w.document
        .querySelector(".pt-count")
        .textContent.replace(/,/g, "")
        .match(/(\d+) of/)[1],
      10,
    );
    eq("still loads normally when pinned", total, truth.total);
  }
  {
    const w = await mount(csvText); // no global set
    ok(
      "falls back to main when unpinned",
      w.__fetched[0] === "https://raw.githubusercontent.com/jhuapl-bio/taxtriage/main/assets/pathogen_sheet.csv",
      w.__fetched[0],
    );
  }
  {
    const w = await mount("", { fail: true, docsRef: { ref: "v3.1.0", version: "3.1" } });
    const html = w.document.getElementById("pathogen-app").innerHTML;
    ok("error link points at the pinned ref", /blob\/v3\.1\.0\/assets/.test(html), html.slice(0, 200));
  }

  console.log(`\n${checks - failures}/${checks} checks passed`);
  if (failures) {
    console.log(`${failures} FAILED`);
    process.exit(1);
  }
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
