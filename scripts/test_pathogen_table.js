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
  const files = [];
  let pendingBlob = null;
  window.URL.createObjectURL = (b) => {
    pendingBlob = b;
    return "blob:stub";
  };
  window.URL.revokeObjectURL = () => {};
  window.HTMLAnchorElement.prototype.click = function () {
    if (pendingBlob) {
      files.push({ name: this.download, blob: pendingBlob });
      pendingBlob = null;
    }
  };
  window.__files = files;
  const opened = [];
  window.open = (u) => {
    opened.push(u);
    return null;
  };
  window.__opened = opened;
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
  eq("column headers", $$(".pt-table thead th").length, 14);
  {
    const heads = $$(".pt-table thead th").map((t) => t.dataset.key);
    ok("assembly column present", heads.includes("assembly_accession"), heads.join(","));
    ok("reference column present", heads.includes("reference"), heads.join(","));
  }

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
    ok("table still has all 14 columns", d.querySelectorAll(".pt-table thead th").length === 14);
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

  /* ── request-an-organism CTA ────────────────────────────────────────────
   * Only current docs (`latest`, PR previews, local serve) may invite edits;
   * frozen release docs must not, since the request would apply to a sheet
   * that version does not track.
   */
  console.log("\nrequest entry points");
  {
    const w = await mount(csvText, { docsRef: { ref: "main", version: "latest" } });
    const d = w.document;
    const btn = d.querySelector(".pt-toolbar button.pt-request");
    ok("toolbar button shown in latest", !!btn, "no request button");
    btn.dispatchEvent(new w.Event("click", { bubbles: true }));
    const alt = d.querySelector(".pt-request-drawer .pt-rq-alt a");
    ok(
      "builder still offers the guided single-organism form",
      alt && /\/issues\/new\?.*template=add_organism\.yml/.test(alt.getAttribute("href")),
      alt && alt.getAttribute("href"),
    );
    ok("guided form opens in a new tab safely", alt && alt.rel === "noopener" && alt.target === "_blank");
  }
  {
    const w = await mount(csvText, { docsRef: { ref: "v3.3.9", version: "3.3.9" } });
    ok(
      "hidden in a frozen release version",
      !w.document.querySelector(".pt-request"),
      "request entry point leaked into a release build",
    );
  }
  {
    const w = await mount(csvText, { docsRef: { ref: "dev", version: "dev" } });
    ok("shown for a branch-named build", !!w.document.querySelector(".pt-toolbar button.pt-request"));
  }

  console.log("\nrequest builder: staging multiple organisms");
  {
    const w = await mount(csvText, { docsRef: { ref: "main", version: "latest" } });
    const d = w.document;
    const fire = (el, type) => el.dispatchEvent(new w.Event(type, { bubbles: true }));
    const drawer = () => d.querySelector(".pt-request-drawer");

    d.querySelector(".pt-toolbar .pt-request").dispatchEvent(new w.Event("click", { bubbles: true }));
    ok("builder opens", drawer().classList.contains("is-open"));
    ok("has an Add another control", !!drawer().querySelector(".pt-rq-add"));

    const set = (id, val) => {
      const n = drawer().querySelector(`[data-f="${id}"]`);
      n.value = val;
      fire(n, "change");
    };
    const selectSites = (id, vals) => {
      const n = drawer().querySelector(`[data-f="${id}"]`);
      Array.from(n.options).forEach((o) => {
        o.selected = vals.includes(o.value);
      });
      fire(n, "change");
    };
    const click = (sel) =>
      drawer()
        .querySelector(sel)
        .dispatchEvent(new w.Event("click", { bubbles: true }));

    // Required-field validation must block, not silently proceed.
    set("name", "Testium unum");
    click(".pt-rq-add");
    const err = drawer().querySelector(".pt-rq-error");
    ok("blocks on missing required fields", !err.hidden, "no error shown");
    eq("nothing staged yet", drawer().querySelectorAll(".pt-rq-list li").length, 0);

    set("taxid", "not-a-number");
    set("general_classification", "primary");
    set("reference", "Doe J. Journal 2020. doi:10/x");
    click(".pt-rq-add");
    ok("rejects a non-numeric tax ID", /numeric/i.test(drawer().querySelector(".pt-rq-error").textContent));

    set("taxid", "1313");
    selectSites("pathogenic_sites", ["blood", "resp"]);
    click(".pt-rq-add");
    eq("first organism staged", drawer().querySelectorAll(".pt-rq-list li").length, 1);
    eq("form cleared for the next", drawer().querySelector('[data-f="name"]').value, "");

    set("name", "Testium duo");
    set("taxid", "999");
    set("general_classification", "commensal");
    set("reference", "Roe R. Journal 2021.");
    click(".pt-rq-add");
    eq("second organism staged", drawer().querySelectorAll(".pt-rq-list li").length, 2);

    // Remove one, then re-add a third.
    drawer()
      .querySelector(".pt-rq-del")
      .dispatchEvent(new w.Event("click", { bubbles: true }));
    eq("staged entry removable", drawer().querySelectorAll(".pt-rq-list li").length, 1);

    set("name", "Testium tria");
    set("taxid", "555");
    set("general_classification", "potential");
    set("reference", "Poe P. Journal 2022.");
    click(".pt-rq-add");
    eq("back to two staged", drawer().querySelectorAll(".pt-rq-list li").length, 2);

    click(".pt-rq-open");
    eq("one issue opened", w.__opened.length, 1);
    const url = w.__opened[0];
    const body = decodeURIComponent(new URL(url).searchParams.get("body"));
    const title = new URL(url).searchParams.get("title");
    ok("title reflects the count", /Add 2 organisms/.test(decodeURIComponent(title)), title);
    ok("body lists both organisms", body.includes("Testium duo") && body.includes("Testium tria"), body.slice(0, 200));
    ok("removed organism is absent", !body.includes("Testium unum"));
    ok(
      "body carries a paste-ready CSV block",
      body.includes("```csv") && body.includes("name,taxid,general_classification"),
    );
    ok("issue is labelled", url.includes("labels=pathogen-sheet"));

    const csvLines = body.split("```csv")[1].split("```")[0].trim().split("\n");
    eq("CSV block has header + 2 rows", csvLines.length, 3);
    eq("CSV columns match the sheet", csvLines[0].split(",").length, 21);
  }

  console.log("\nrequest builder: single organism without staging");
  {
    const w = await mount(csvText, { docsRef: { ref: "main", version: "latest" } });
    const d = w.document;
    const drawer = () => d.querySelector(".pt-request-drawer");
    const fire = (el, type) => el.dispatchEvent(new w.Event(type, { bubbles: true }));
    d.querySelector(".pt-toolbar .pt-request").dispatchEvent(new w.Event("click", { bubbles: true }));
    const set = (id, val) => {
      const n = drawer().querySelector(`[data-f="${id}"]`);
      n.value = val;
      fire(n, "change");
    };
    set("name", "Solo organism");
    set("taxid", "42");
    set("general_classification", "primary");
    set("reference", "Solo S. Journal 2023.");
    drawer()
      .querySelector(".pt-rq-open")
      .dispatchEvent(new w.Event("click", { bubbles: true }));
    eq("opens without pressing Add another", w.__opened.length, 1);
    const t = decodeURIComponent(new URL(w.__opened[0]).searchParams.get("title"));
    eq("singular title", t, "Add organism: Solo organism");
  }

  console.log("\nrequest builder: gating and prefill");
  {
    const w = await mount(csvText, { docsRef: { ref: "v3.3.9", version: "3.3.9" } });
    ok("no builder button in a frozen release", !w.document.querySelector(".pt-toolbar .pt-request"));
  }
  {
    const w = await mount(csvText, { docsRef: { ref: "main", version: "latest" } });
    const d = w.document;
    d.querySelector(".pt-search input").value = "Nocardia zzz";
    d.querySelector(".pt-search input").dispatchEvent(new w.Event("input", { bubbles: true }));
    await new Promise((r) => setTimeout(r, 220));
    const cta = d.querySelector(".pt-empty .pt-request");
    ok("empty state offers the builder", !!cta);
    cta.dispatchEvent(new w.Event("click", { bubbles: true }));
    const drawer = d.querySelector(".pt-request-drawer");
    ok("builder opens from the empty state", drawer.classList.contains("is-open"));
    eq("organism name prefilled from the query", drawer.querySelector('[data-f="name"]').value, "Nocardia zzz");
  }

  console.log("\nassembly and reference columns");
  {
    const w = await mount(csvText);
    const d = w.document;
    const heads = Array.from(d.querySelectorAll(".pt-table thead th")).map((t) => t.dataset.key);
    const accIdx = heads.indexOf("assembly_accession");
    const refIdx = heads.indexOf("reference");

    // Find a row that actually carries both, so the assertions are meaningful.
    const rows = Array.from(d.querySelectorAll(".pt-table tbody tr"));
    const withAcc = rows.find((r) => r.children[accIdx].textContent.trim());
    ok("some row shows an assembly accession", !!withAcc);
    const link = withAcc && withAcc.children[accIdx].querySelector("a");
    ok(
      "assembly links to NCBI datasets",
      link && /ncbi\.nlm\.nih\.gov\/datasets\/genome\//.test(link.getAttribute("href")),
      link && link.getAttribute("href"),
    );

    const withRef = rows.find((r) => r.children[refIdx].textContent.trim());
    ok("some row shows a reference", !!withRef);
    ok(
      "reference cell keeps the full text in a title",
      withRef && withRef.children[refIdx].getAttribute("title").length > 0,
    );
    ok("reference cell is clamped for layout", withRef && withRef.children[refIdx].classList.contains("pt-ref"));

    // Export must still emit the canonical schema regardless of table columns.
    d.querySelector(".pt-export").dispatchEvent(new w.Event("click", { bubbles: true }));
  }

  console.log("\nfull-width page");
  {
    const w = await mount(csvText);
    ok("body marked full width on this page", w.document.body.classList.contains("pt-fullwidth"));
  }

  console.log("\nrequest builder: download for private sharing");
  {
    const w = await mount(csvText, { docsRef: { ref: "main", version: "latest" } });
    const d = w.document;
    const fire = (el, type) => el.dispatchEvent(new w.Event(type, { bubbles: true }));
    const drawer = () => d.querySelector(".pt-request-drawer");
    d.querySelector(".pt-toolbar .pt-request").dispatchEvent(new w.Event("click", { bubbles: true }));
    ok("download button offered alongside the issue", !!drawer().querySelector(".pt-rq-download"));

    const set = (id, val) => {
      const n = drawer().querySelector(`[data-f="${id}"]`);
      n.value = val;
      fire(n, "change");
    };
    set("name", "Private organism");
    set("taxid", "77");
    set("general_classification", "primary");
    set("reference", "Confidential source 2024.");
    fire(drawer().querySelector(".pt-rq-add"), "click");
    set("name", "Second private");
    set("taxid", "78");
    set("general_classification", "potential");
    set("reference", "Another source 2024.");
    fire(drawer().querySelector(".pt-rq-download"), "click");

    eq("a file was produced", w.__files.length, 1);
    eq("sensible filename", w.__files[0].name, "taxtriage-organism-request.md");
    const text = await w.__files[0].blob.text();
    ok(
      "contains both organisms",
      text.includes("Private organism") && text.includes("Second private"),
      text.slice(0, 120),
    );
    ok("contains the paste-ready CSV block", text.includes("```csv"));
    eq("no issue was opened", w.__opened.length, 0);
  }
  {
    // Downloading with nothing entered must report, not emit an empty file.
    const w = await mount(csvText, { docsRef: { ref: "main", version: "latest" } });
    const d = w.document;
    d.querySelector(".pt-toolbar .pt-request").dispatchEvent(new w.Event("click", { bubbles: true }));
    d.querySelector(".pt-rq-download").dispatchEvent(new w.Event("click", { bubbles: true }));
    eq("no empty file written", w.__files.length, 0);
    ok("explains why", !d.querySelector(".pt-rq-error").hidden);
  }

  /* ── provenance (request_type) ──────────────────────────────────────────
   * Every existing row is APL-derived. New requests are stamped from the route
   * actually used, and external requesters may not claim APL-derived.
   */
  console.log("\nrequest_type provenance");
  {
    const w = await mount(csvText);
    const d = w.document;
    const heads = Array.from(d.querySelectorAll(".pt-table thead th")).map((t) => t.dataset.key);
    ok("request_type is a table column", heads.includes("request_type"), heads.join(","));
    ok("request_type is a facet", !!d.querySelector('.pt-facet[data-key="request_type"]'));
    const opts = Array.from(d.querySelectorAll('.pt-facet[data-key="request_type"] input')).map((i) => i.value);
    eq("existing sheet is uniformly APL-derived", JSON.stringify(opts), JSON.stringify(["APL-derived"]));
  }
  {
    const w = await mount(csvText, { docsRef: { ref: "main", version: "latest" } });
    const d = w.document;
    const fire = (el, type) => el.dispatchEvent(new w.Event(type, { bubbles: true }));
    const drawer = () => d.querySelector(".pt-request-drawer");
    d.querySelector(".pt-toolbar .pt-request").dispatchEvent(new w.Event("click", { bubbles: true }));

    const sel = drawer().querySelector('[data-f="request_type"]');
    ok("builder offers a request type picker", !!sel);
    const apl = Array.from(sel.options).find((o) => o.value === "APL-derived");
    ok("APL-derived is visible in the list", !!apl);
    ok("APL-derived cannot be chosen", apl && apl.disabled, "option is selectable");
    ok(
      "the other two are selectable",
      Array.from(sel.options)
        .filter((o) => !o.disabled)
        .map((o) => o.value)
        .join(",") === "git-tracked,external-local",
      Array.from(sel.options)
        .map((o) => o.value + (o.disabled ? "(disabled)" : ""))
        .join(","),
    );

    const set = (id, val) => {
      const n = drawer().querySelector(`[data-f="${id}"]`);
      n.value = val;
      fire(n, "change");
    };
    set("name", "Route test");
    set("taxid", "12");
    set("general_classification", "primary");
    set("reference", "Ref 2024.");
    // Deliberately pick the wrong route, then download: the file must win.
    set("request_type", "git-tracked");
    fire(drawer().querySelector(".pt-rq-download"), "click");
    const text = await w.__files[0].blob.text();
    ok("download stamps external-local", /external-local/.test(text), text.slice(0, 300));
    ok("download does not record git-tracked", !/git-tracked/.test(text));
  }
  {
    const w = await mount(csvText, { docsRef: { ref: "main", version: "latest" } });
    const d = w.document;
    const fire = (el, type) => el.dispatchEvent(new w.Event(type, { bubbles: true }));
    const drawer = () => d.querySelector(".pt-request-drawer");
    d.querySelector(".pt-toolbar .pt-request").dispatchEvent(new w.Event("click", { bubbles: true }));
    const set = (id, val) => {
      const n = drawer().querySelector(`[data-f="${id}"]`);
      n.value = val;
      fire(n, "change");
    };
    set("name", "Issue route");
    set("taxid", "13");
    set("general_classification", "primary");
    set("reference", "Ref 2024.");
    set("request_type", "external-local");
    fire(drawer().querySelector(".pt-rq-open"), "click");
    const body = decodeURIComponent(new w.URL(w.__opened[0]).searchParams.get("body"));
    ok("issue stamps git-tracked", /git-tracked/.test(body), body.slice(0, 300));
    ok("issue does not record external-local", !/external-local/.test(body));
    const csvLines = body.split("```csv")[1].split("```")[0].trim().split("\n");
    eq("CSV block carries the new column", csvLines[0].split(",").length, 21);
    ok("last CSV field is the provenance", csvLines[1].trim().endsWith("git-tracked"), csvLines[1]);
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
