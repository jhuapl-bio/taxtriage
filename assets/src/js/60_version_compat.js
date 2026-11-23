/* ═══════════════════════════════════════════════════════════════════════════
       -  §  VERSION COMPATIBILITY  (older TaxTriage outputs in this report)
       -     Everything the report can be handed — per-sample .paths.json,
       -     combined all.odr.json, .organisms.report.txt / .tsv / .xlsx tables,
       -     an older interactive report (.html) and saved state files — is
       -     checked here for the TaxTriage version that produced it and for
       -     what it is missing, so older files load with the right scale and
       -     the user is told, in one dialog, what will and will not work.
       -
       -       ttCompatFixScales(rows)    1.x tables stored TASS (and breadth)
       -                                  as 0–1; today's report is 0–100.
       -                                  Rescales in place; returns what changed.
       -       ttCompatAssess(kind, …)    one entry per file: version, era,
       -                                  status (ok / limited / failed), notes.
       -       ttCompatDialog(entries)    the "older version" dialog.
       -       ttCompatAfterUpload(…)     status-line note + dialog when needed.
       -       ttExtractReportBoot(html)  the data block of an older report.
       -
       -     TT_STATE_FORMAT is the version of the exported state envelope
       -     (29_session_state.js): 1 = data + UI, 2 = + app state,
       -     3 = + report_format / saved_by stamps. Older envelopes load;
       -     a newer one is attempted and the user is warned.
═══════════════════════════════════════════════════════════════════════════ */

const TT_STATE_FORMAT = 3;
// Outputs from before this pipeline release are labelled "older" in the dialog.
const TT_COMPAT_BASELINE = "v3.0.0";
// Columns the report cannot work without.
const TT_COMPAT_REQUIRED = ["Detected Organism", "Specimen ID"];

function _ttCompatEsc(v) {
  return String(v == null ? "" : v)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

/** "v3.3.11" / "3.3.11" / "v3.3.11-dev" → [3, 3, 11]; null when not a version. */
function ttVersionParse(v) {
  const m = String(v == null ? "" : v)
    .trim()
    .match(/^v?(\d+)(?:\.(\d+))?(?:\.(\d+))?/i);
  return m ? [Number(m[1]), Number(m[2] || 0), Number(m[3] || 0)] : null;
}
/** < 0, 0, > 0 like a comparator; null when either is not a version. */
function ttVersionCmp(a, b) {
  const x = ttVersionParse(a),
    y = ttVersionParse(b);
  if (!x || !y) return null;
  for (let i = 0; i < 3; i++) if (x[i] !== y[i]) return x[i] - y[i];
  return 0;
}

/* ── legacy scales ─────────────────────────────────────────────────────── */
function _ttCompatNums(rows, col) {
  const out = [];
  for (const r of rows) {
    const v = r ? r[col] : null;
    if (v === null || v === undefined || v === "") continue;
    const n = Number(String(v).replace(/%/g, ""));
    if (isFinite(n)) out.push(n);
  }
  return out;
}
/** 1.x tables stored TASS as 0–1 (and Coverage / breadth as a 0–1 fraction).
 *  Rescale in place to today's 0–100 when the whole file is on the old scale.
 *  Coverage is only touched when TASS was, so a modern file whose detections
 *  all happen to have < 1 % breadth is never inflated. */
function ttCompatFixScales(rows) {
  const info = { tass: false, coverage: false };
  if (!Array.isArray(rows) || !rows.length) return info;
  const tass = _ttCompatNums(rows, "TASS Score");
  if (!tass.length) return info;
  const max = Math.max.apply(null, tass);
  if (!(max <= 1 && tass.some((v) => v > 0 && v < 1))) return info;
  info.tass = true;
  const scale = (r, c, d) => {
    const n = Number(String(r[c]).replace(/%/g, ""));
    if (r[c] !== null && r[c] !== undefined && r[c] !== "" && isFinite(n)) r[c] = Math.round(n * 100 * d) / d;
  };
  rows.forEach((r) => {
    scale(r, "TASS Score", 10);
    ["Species TASS", "Genus TASS"].forEach((c) => c in r && scale(r, c, 10));
  });
  const cov = _ttCompatNums(rows, "Coverage");
  if (cov.length && Math.max.apply(null, cov) <= 1) {
    info.coverage = true;
    rows.forEach((r) => {
      scale(r, "Coverage", 100);
      if ("Breadth %" in r) scale(r, "Breadth %", 100);
    });
  }
  return info;
}

/* ── per-file assessment ───────────────────────────────────────────────── */
/** kind: "paths" | "combined" | "table" | "report" | "state" | "unknown".
 *  ctx: { name, rows, revisions[], scale, fatal, formatVersion }. */
function ttCompatAssess(kind, ctx) {
  const c = ctx || {};
  const rows = Array.isArray(c.rows) ? c.rows : [];
  const revs = Array.from(new Set((c.revisions || []).filter((v) => v && !/^(na|null|none)$/i.test(String(v)))));
  const notes = [];
  let status = "ok";
  let version = revs.length ? revs.join(", ") : "";
  let era = "";

  if (c.fatal) {
    return { name: c.name, kind, version: version || "unknown", era: c.era || "", status: "failed", notes: [c.fatal] };
  }
  const cols = new Set();
  rows.slice(0, 200).forEach((r) => r && Object.keys(r).forEach((k) => cols.add(k)));
  if (rows.length) {
    const missing = TT_COMPAT_REQUIRED.filter((k) => !cols.has(k));
    if (missing.length)
      return {
        name: c.name,
        kind,
        version: version || "unknown",
        era: "",
        status: "failed",
        notes: [
          `Not a TaxTriage detections file this report can read: no ${missing
            .map((m) => `“${m}”`)
            .join(" / ")} column.`,
        ],
      };
  }
  if (!version) {
    if (c.scale && c.scale.tass) era = "TaxTriage 1.x table (0–1 TASS)";
    else if (kind === "table" && cols.has("Passes Threshold")) era = "TaxTriage 3.x table";
    else if (kind === "table" && cols.has("TASS Score")) era = "TaxTriage 2.x / 3.x table";
  } else {
    const old = revs.some(
      (v) => ttVersionCmp(v, TT_COMPAT_BASELINE) !== null && ttVersionCmp(v, TT_COMPAT_BASELINE) < 0,
    );
    if (old) era = `older than ${TT_COMPAT_BASELINE}`;
  }

  if (c.scale && c.scale.tass) {
    notes.push("TASS was on the old 0–1 scale; rescaled to 0–100 so the cutoffs apply.");
    status = "limited";
  }
  if (c.scale && c.scale.coverage) notes.push("Coverage was a 0–1 fraction; rescaled to %.");
  if (rows.length) {
    if (!cols.has("TASS Score")) {
      notes.push("No TASS Score: detections cannot be ranked or filtered by confidence.");
      status = "limited";
    }
    if (!cols.has("Taxonomic ID #")) notes.push("No taxonomic IDs: organisms are matched across samples by name.");
    if (kind === "table" || kind === "report") {
      if (!cols.has("Level") && kind === "table")
        notes.push("Species / genus rows are rebuilt in the browser (their TASS is the best strain's).");
    }
    if (kind === "table")
      notes.push(
        "Tables carry no coverage histograms or sample metadata: Histograms, Metadata and Mapping stay empty for these samples.",
      );
    if (!cols.has("Breadth %") && !cols.has("Coverage"))
      notes.push("No breadth / coverage column: coverage views and the crossover coverage term are off.");
    if (!cols.has("RPM") && !c.hasDepth)
      notes.push("No RPM or input read counts: abundance is compared as a share of aligned reads.");
    const lineage = ["Genus", "Family", "Order", "Phylum"].filter((k) => cols.has(k));
    if (lineage.length < 2) notes.push("Little or no taxonomic lineage: grouping by rank is limited.");
  }
  if (era && /older|1\.x/.test(era) && status === "ok") status = "limited";
  return { name: c.name, kind, version: version || "not recorded", era, status, notes };
}

/** The data block of an interactive report (any version that embeds one):
 *  `window.HEATMAP_BOOT = {...};` or `<script id="BOOTSTRAP" type="application/json">`.
 *  Returns the parsed object, or null. */
function ttExtractReportBoot(html) {
  const s = String(html || "");
  const tag = s.match(/<script[^>]*id=["']BOOTSTRAP["'][^>]*type=["']application\/json["'][^>]*>([\s\S]*?)<\/script>/i);
  if (tag) {
    try {
      return JSON.parse(tag[1]);
    } catch (e) {
      /* fall through */
    }
  }
  const at = s.indexOf("window.HEATMAP_BOOT");
  if (at < 0) return null;
  const eq = s.indexOf("=", at);
  const start = s.indexOf("{", eq);
  const end = s.indexOf("<" + "/script>", start); // split: this file is inlined into a <script>
  if (start < 0 || end < 0) return null;
  let body = s.slice(start, end).trim();
  body = body.replace(/;\s*$/, "");
  try {
    return JSON.parse(body);
  } catch (e) {
    // A later statement in the same <script>: cut at the last "};".
    const cut = body.lastIndexOf("}");
    try {
      return JSON.parse(body.slice(0, cut + 1));
    } catch (e2) {
      return null;
    }
  }
}

/* ── the dialog ────────────────────────────────────────────────────────── */
const _TT_COMPAT_STATUS = {
  ok: { icon: "fa-circle-check", label: "Loaded", cls: "ok" },
  limited: { icon: "fa-circle-info", label: "Loaded with limits", cls: "limited" },
  failed: { icon: "fa-circle-xmark", label: "Not loaded", cls: "failed" },
  newer: { icon: "fa-triangle-exclamation", label: "Newer format", cls: "limited" },
};
let _ttCompatLast = [];
/** entries: from ttCompatAssess. opts.title / opts.lead override the heading. */
function ttCompatDialog(entries, opts) {
  if (typeof document === "undefined") return;
  const list = (entries || []).filter(Boolean);
  if (!list.length) return;
  _ttCompatLast = list;
  const o = opts || {};
  const esc = _ttCompatEsc;
  const anyFailed = list.some((e) => e.status === "failed");
  const anyNewer = list.some((e) => e.status === "newer");
  const title =
    o.title ||
    (anyFailed
      ? "Some files could not be loaded"
      : anyNewer
      ? "Saved by a newer TaxTriage report"
      : "Files from older TaxTriage versions");
  const lead =
    o.lead ||
    (anyFailed
      ? "These files are from a TaxTriage version (or in a format) this report cannot read. Anything else you loaded is unaffected."
      : anyNewer
      ? "This file was saved by a newer version of the report. Its data loaded; settings this version does not know about were skipped."
      : "These files were produced by an older TaxTriage version. They loaded, but some views have less to show:");
  let ov = document.getElementById("tt-compat-overlay");
  if (ov) ov.remove();
  ov = document.createElement("div");
  ov.id = "tt-compat-overlay";
  ov.className = "tt-compat-overlay";
  const rows = list
    .map((e) => {
      const st = _TT_COMPAT_STATUS[e.status] || _TT_COMPAT_STATUS.ok;
      return (
        `<tr class="tt-compat-${st.cls}"><td class="tt-compat-file">${esc(e.name)}</td>` +
        `<td>${esc(e.version && e.version !== "unknown" ? e.version : "—")}${
          e.era ? `<div class="tt-compat-era">${esc(e.era)}</div>` : ""
        }</td>` +
        `<td class="tt-compat-st"><i class="fas ${st.icon}"></i> ${st.label}</td>` +
        `<td>${
          e.notes && e.notes.length ? `<ul>${e.notes.map((n) => `<li>${esc(n)}</li>`).join("")}</ul>` : "—"
        }</td></tr>`
      );
    })
    .join("");
  const help =
    anyNewer && !anyFailed
      ? `<div class="tt-compat-help">Open it in a report built by the newer pipeline to get everything back.</div>`
      : anyFailed
      ? `<div class="tt-compat-help"><b>What you can do</b><ul>` +
        `<li>Load the run's <code>*.paths.json</code> or <code>all.odr.json</code> instead (any 3.x version), or its <code>*.organisms.report.txt</code> tables (any version).</li>` +
        `<li>Reports older than the interactive format (PDF / static HTML) carry no data to import: use the run's report tables.</li>` +
        `<li>Re-run the report step with the current pipeline (<code>-resume</code> reuses the alignments) for every feature.</li></ul></div>`
      : `<div class="tt-compat-help">Re-running the report step with the current pipeline (<code>-resume</code> reuses the alignments) restores every feature for these samples.</div>`;
  ov.innerHTML =
    `<div class="tt-compat-modal" role="dialog" aria-modal="true" aria-labelledby="tt-compat-title">` +
    `<header><i class="fas ${anyFailed ? "fa-triangle-exclamation" : "fa-clock-rotate-left"}"></i>` +
    `<span id="tt-compat-title">${esc(title)}</span>` +
    `<button type="button" class="tt-compat-x" aria-label="Close">&times;</button></header>` +
    `<div class="tt-compat-body"><p>${esc(lead)}</p>` +
    `<div class="tt-compat-scroll"><table class="tt-compat-table"><thead><tr><th>File</th><th>Produced by</th><th>Status</th><th>Notes</th></tr></thead>` +
    `<tbody>${rows}</tbody></table></div>${help}</div>` +
    `<footer><span class="tt-compat-this">This report: state format ${TT_STATE_FORMAT}${
      typeof BOOT !== "undefined" && BOOT && BOOT.pipeline_revision
        ? ` · built by TaxTriage ${esc(BOOT.pipeline_revision)}`
        : ""
    }</span><button type="button" class="tt-compat-okbtn">OK</button></footer></div>`;
  document.body.appendChild(ov);
  const close = () => {
    ov.remove();
    document.removeEventListener("keydown", onKey, true);
  };
  const onKey = (e) => {
    if (e.key === "Escape") close();
  };
  document.addEventListener("keydown", onKey, true);
  ov.addEventListener("click", (e) => {
    if (e.target === ov) close();
  });
  ov.querySelector(".tt-compat-x").addEventListener("click", close);
  const ok = ov.querySelector(".tt-compat-okbtn");
  ok.addEventListener("click", close);
  setTimeout(() => ok.focus(), 20);
}
/** Re-open the last dialog (status-line "details" link). */
function ttCompatShowLast() {
  if (_ttCompatLast.length) ttCompatDialog(_ttCompatLast);
}

/** After an upload batch: a status line, and the dialog when it matters
 *  (a file failed, an old 0–1 scale was rescaled, or the producer is older
 *  than TT_COMPAT_BASELINE). Plain notes on current files stay in the line. */
function ttCompatAfterUpload(entries, statusEl) {
  const list = (entries || []).filter(Boolean);
  const flagged = list.filter((e) => e.status !== "ok");
  if (!flagged.length) return;
  if (statusEl) {
    const line = document.createElement("div");
    line.className = "tt-compat-line";
    const nFail = flagged.filter((e) => e.status === "failed").length;
    const nOld = flagged.length - nFail;
    const parts = [];
    if (nOld) parts.push(`${nOld} file${nOld === 1 ? "" : "s"} from an older TaxTriage version`);
    if (nFail) parts.push(`${nFail} file${nFail === 1 ? "" : "s"} could not be loaded`);
    line.classList.toggle("failed", nFail > 0);
    line.innerHTML =
      `<i class="fas ${nFail ? "fa-triangle-exclamation" : "fa-clock-rotate-left"}"></i> ${parts.join(" · ")}` +
      ` · <a href="#" class="tt-compat-more">details</a>`;
    line.querySelector("a").addEventListener("click", (e) => {
      e.preventDefault();
      ttCompatDialog(flagged);
    });
    statusEl.appendChild(line);
  }
  ttCompatDialog(flagged);
}
