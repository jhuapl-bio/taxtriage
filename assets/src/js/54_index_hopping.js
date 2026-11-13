/* ═══════════════════════════════════════════════════════════════════════════
       -  §  INDEX HOPPING  (Illumina cross-sample bleed-through; flag only)
       -     Two sources, one display:
       -       pipeline  --index_hopping ran bin/index_hopping.py; the rows
       -                 already carry the "Index Hop …" columns and each
       -                 sample's SAMPLE_META[s].index_hopping says how.
       -       browser   an index-hopping table (CSV / TSV / XLSX, bcl-convert
       -                 Index_Hopping_Counts layout) dropped on the upload
       -                 zone or loaded from the banner badge. The same test as
       -                 bin/index_hopping.py runs here over the loaded rows
       -                 (keep the two in step: same parser, same rates, same
       -                 Poisson tail, taxa matched by Level + Taxonomic ID #).
       -
       -     For a detection with N reads in sample R, every other sample D of
       -     the same lane with MORE reads of that taxon is a donor:
       -         expected = Σ reads_D × rate(D, R)
       -         flagged  = expected > 0 and P(X ≥ N | Poisson(factor × expected)) ≥ alpha
       -     Rates: "pair" (hopped rows of the two samples), "lane" (lane-wide
       -     hop fraction), "default" (no hopped rows: params.rate).
       -
       -     Results land on the rows as the "Index Hop …" columns, a default
       -     Organism QC rule (Index Hop Suspect == true → flag) highlights
       -     them, and the banner badge (#banner-hop) summarises per sample.
       -     TASS and calls are never changed.
═══════════════════════════════════════════════════════════════════════════ */

const TT_HOP_COLS = [
  "Index Hop Suspect",
  "Index Hop Expected Reads",
  "Index Hop Explained %",
  "Index Hop p",
  "Index Hop Source",
  "Index Hop Source Reads",
  "Index Hop Rate %",
  "Index Hop Rate Source",
  "Index Hop Source Of",
  "Index Hop Reads Out",
];
const TT_HOP_NUMERIC = [
  "Index Hop Expected Reads",
  "Index Hop Explained %",
  "Index Hop p",
  "Index Hop Source Reads",
  "Index Hop Rate %",
  "Index Hop Reads Out",
];
const TT_HOP_DEFAULTS = { rate: 0.001, factor: 2, alpha: 0.01 };
// Donor-side notes only for a donor expected to have sent at least this many
// reads into a flagged detection (MIN_SOURCE_READS in bin/index_hopping.py).
const TT_HOP_MIN_SOURCE_READS = 0.5;

const TT_HOP = {
  tables: [], // tables loaded in this browser: [{ name, lanes, hasIndex2 }]
  params: Object.assign({}, TT_HOP_DEFAULTS),
  result: null, // browser run: { samples: {sid: block like metadata.index_hopping}, files, at }
};

/* ── small helpers ───────────────────────────────────────────────────────── */
function _ttHopEsc(v) {
  return String(v == null ? "" : v)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}
function _ttHopRound(v, d) {
  // toFixed rounds the exact binary value, as Python's round() does
  // (13.645 -> 13.64), so browser and pipeline columns agree to the digit.
  const n = Number(v);
  return isFinite(n) ? Number(n.toFixed(d)) : 0;
}
function _ttHopNormId(s) {
  return String(s == null ? "" : s)
    .trim()
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "_")
    .replace(/^_+|_+$/g, "");
}
const _TT_HOP_BCL_SUFFIX = /_s\d+(_l\d{3})?(_r[12](_001)?)?$/;
function _ttHopIdCandidates(name, override) {
  if (override) return [_ttHopNormId(override)];
  const base = _ttHopNormId(name);
  const out = [base];
  const stripped = base.replace(_TT_HOP_BCL_SUFFIX, "");
  if (stripped && stripped !== base) out.push(stripped);
  return out;
}
function _ttHopPoissonSf(k, lam) {
  // P(X >= k), X ~ Poisson(lam) -- same arithmetic as bin/index_hopping.py
  if (k <= 0) return 1.0;
  if (lam <= 0) return 0.0;
  if (k > lam + 12.0 * Math.sqrt(lam) + 30.0) return 0.0;
  let term = Math.exp(-lam);
  let cdf = term;
  for (let i = 1; i < Math.trunc(k); i++) {
    term *= lam / i;
    cdf += term;
  }
  return Math.max(0, Math.min(1, 1 - cdf));
}

/* ── table parsing (mirrors parse_hopping_file) ──────────────────────────── */
function _ttHopNormHeader(h) {
  return String(h == null ? "" : h)
    .trim()
    .toLowerCase()
    .replace(/#/g, " ")
    .replace(/%/g, " pct ")
    .replace(/[^a-z0-9]+/g, "_")
    .replace(/^_+|_+$/g, "");
}
const _TT_HOP_ALIASES = {
  lane: ["lane"],
  sample: ["sampleid", "sample_id", "sample", "sample_name", "samplename", "sample_id_s"],
  project: ["sample_project", "project", "sampleproject"],
  i7: ["index", "index1", "index_1", "i7", "i7_index", "index_i7", "i7_index_id_sequence"],
  i5: ["index2", "index_2", "i5", "i5_index", "index_i5"],
  reads: ["reads", "read_count", "readcount", "num_reads", "clusters", "count", "number_of_reads", "reads_count"],
};
function _ttHopMatchHeader(row) {
  const cols = {};
  (row || []).forEach((raw, i) => {
    const h = _ttHopNormHeader(raw);
    Object.keys(_TT_HOP_ALIASES).forEach((field) => {
      if (_TT_HOP_ALIASES[field].includes(h) && !(field in cols)) cols[field] = i;
    });
  });
  return "reads" in cols && ("i7" in cols || "sample" in cols) ? cols : null;
}
function _ttHopNum(v) {
  const s = String(v == null ? "" : v)
    .trim()
    .replace(/[,\s]/g, "");
  if (!s) return null;
  const n = Number(s);
  return isFinite(n) ? n : null;
}
function _ttHopSeq(v) {
  return String(v == null ? "" : v)
    .toUpperCase()
    .replace(/[^ACGTN]/g, "");
}

/** rows (arrays of cells) → table, or null when no hopping header is found.
 *  `strict` (used when sniffing dropped files) also requires an index column,
 *  so an ordinary detections TSV is never mistaken for a hopping table. */
function ttHopParseRows(tables, name, strict) {
  let header = null;
  let body = null;
  for (const rows of tables) {
    for (let k = 0; k < rows.length; k++) {
      const cols = _ttHopMatchHeader(rows[k]);
      if (cols && (!strict || "i7" in cols)) {
        header = cols;
        body = rows.slice(k + 1);
        break;
      }
    }
    if (header) break;
  }
  if (!header) return null;
  const cell = (r, f) => {
    const i = header[f];
    return i != null && i < r.length && r[i] != null ? String(r[i]).trim() : "";
  };
  const records = [];
  for (const r of body) {
    if (!r || !r.some((x) => String(x == null ? "" : x).trim())) continue;
    const reads = _ttHopNum(cell(r, "reads"));
    if (reads == null) continue; // repeated header, footer or section title
    records.push({
      lane: cell(r, "lane") || "all",
      sample: cell(r, "sample"),
      project: cell(r, "project"),
      i7: _ttHopSeq(cell(r, "i7")),
      i5: _ttHopSeq(cell(r, "i5")),
      reads,
    });
  }
  const own = new Map(); // "lane\u0001sample" -> record with most reads
  records.forEach((rec) => {
    if (!rec.sample) return;
    const k = rec.lane + "\u0001" + rec.sample;
    const cur = own.get(k);
    if (!cur || rec.reads > cur.reads) own.set(k, rec);
  });
  const lanes = {};
  const laneOf = (lane) => lanes[lane] || (lanes[lane] = { samples: {}, hopped: [] });
  own.forEach((rec) => {
    laneOf(rec.lane).samples[rec.sample] = { i7: rec.i7, i5: rec.i5, reads: rec.reads, project: rec.project };
  });
  records.forEach((rec) => {
    if (rec.sample && own.get(rec.lane + "\u0001" + rec.sample) === rec) return;
    laneOf(rec.lane).hopped.push([rec.i7, rec.i5, rec.reads]);
  });
  const hasIndex2 = Object.values(lanes).some((L) => Object.values(L.samples).some((s) => s.i5));
  return { name, lanes, hasIndex2 };
}

/** Quote-aware delimited text → rows; delimiter chosen like read_delimited. */
function _ttHopSplitText(text, name) {
  text = String(text || "").replace(/^﻿/, "");
  const head = text
    .split(/\r?\n/)
    .filter((l) => l.trim())
    .slice(0, 20);
  const count = (ch) => head.reduce((n, l) => n + l.split(ch).length - 1, 0);
  const tabs = count("\t"),
    commas = count(","),
    semis = count(";");
  const lower = String(name || "").toLowerCase();
  let delim = "\t";
  if (tabs && (tabs >= commas || /\.(tsv|txt)$/.test(lower))) delim = "\t";
  else if (commas) delim = ",";
  else if (semis) delim = ";";
  const rows = [];
  let row = [],
    field = "",
    inQ = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (inQ) {
      if (ch === '"') {
        if (text[i + 1] === '"') {
          field += '"';
          i++;
        } else inQ = false;
      } else field += ch;
    } else if (ch === '"' && field === "") inQ = true;
    else if (ch === delim) {
      row.push(field);
      field = "";
    } else if (ch === "\n" || ch === "\r") {
      if (ch === "\r" && text[i + 1] === "\n") i++;
      row.push(field);
      rows.push(row);
      row = [];
      field = "";
    } else field += ch;
  }
  if (field !== "" || row.length) {
    row.push(field);
    rows.push(row);
  }
  return rows;
}

function _ttHopReadFile(file, asBuffer) {
  return new Promise((resolve, reject) => {
    const fr = new FileReader();
    fr.onload = (e) => resolve(e.target.result);
    fr.onerror = () => reject(new Error(`${file.name}: could not be read`));
    if (asBuffer) fr.readAsArrayBuffer(file);
    else fr.readAsText(file);
  });
}

/** File → parsed hopping table, or null when it is not one (so the regular
 *  upload path can take it). Throws only for files that ARE hopping tables
 *  but cannot be used. */
async function ttHopTryParseFile(file, strict) {
  const name = file.name || "table";
  const lower = name.toLowerCase();
  if (/\.(json|gz|bam|fa|fasta|fastq)$/.test(lower)) return null;
  if (/\.(xlsx|xlsm|xls)$/.test(lower)) {
    if (typeof XLSX === "undefined") return null;
    const buf = await _ttHopReadFile(file, true);
    const wb = XLSX.read(new Uint8Array(buf), { type: "array" });
    const sheets = wb.SheetNames.map((sn) =>
      XLSX.utils.sheet_to_json(wb.Sheets[sn], { header: 1, defval: "", raw: false, blankrows: false }),
    );
    return ttHopParseRows(sheets, name, strict !== false);
  }
  if (!/\.(csv|tsv|txt)$/.test(lower)) return null;
  const text = await _ttHopReadFile(file, false);
  return ttHopParseRows([_ttHopSplitText(text, name)], name, strict !== false);
}

/* ── rates (mirrors LaneRates) ───────────────────────────────────────────── */
function _ttHopLaneRates(L) {
  const byI7 = {},
    byI5 = {};
  Object.keys(L.samples).forEach((sid) => {
    const s = L.samples[sid];
    if (s.i7) (byI7[s.i7] = byI7[s.i7] || []).push(sid);
    if (s.i5) (byI5[s.i5] = byI5[s.i5] || []).push(sid);
  });
  const pairHops = new Map();
  let hoppedTotal = 0;
  L.hopped.forEach(([i7, i5, reads]) => {
    const src = byI7[i7] || [],
      dst = byI5[i5] || [];
    const combos = [];
    src.forEach((a) => dst.forEach((b) => a !== b && combos.push([a, b])));
    if (!combos.length) return;
    hoppedTotal += reads;
    const share = reads / combos.length;
    combos.forEach(([a, b]) => {
      const k = [a, b].sort().join("\u0001");
      pairHops.set(k, (pairHops.get(k) || 0) + share);
    });
  });
  const proper = Object.values(L.samples).reduce((n, s) => n + s.reads, 0);
  const denom = hoppedTotal + proper;
  const laneRate = L.hopped.length && denom > 0 ? hoppedTotal / denom : null;
  return {
    laneRate,
    rate(a, b, def) {
      if (laneRate == null) return [def, "default"];
      const k = [a, b].sort().join("\u0001");
      if (pairHops.has(k)) {
        const n = L.samples[a].reads + L.samples[b].reads;
        if (n > 0) return [(pairHops.get(k) + 1) / n, "pair"];
      }
      return [laneRate, "lane"];
    },
  };
}

/* ── the test, over the loaded rows (mirrors analyse) ────────────────────── */
function _ttHopLevel(r) {
  return String(r.Level || "Strain").toLowerCase();
}
function _ttHopOverride(sample, tableName) {
  // A pipeline run records which SampleID each sample matched (incl. the
  // samplesheet index_sample_id); reuse it for the same table.
  const ih = ((typeof SAMPLE_META !== "undefined" && SAMPLE_META[sample]) || {}).index_hopping;
  const f = ih && Array.isArray(ih.files) ? ih.files.find((x) => x && x.file === tableName) : null;
  return f ? f.sample_id : null;
}

function ttHopCompute(rows, tables, params) {
  const P = Object.assign({}, TT_HOP_DEFAULTS, params || {});
  const bySample = new Map();
  rows.forEach((r) => {
    const s = String(r["Specimen ID"] || "");
    if (!s) return;
    if (!bySample.has(s)) bySample.set(s, []);
    bySample.get(s).push(r);
  });
  const names = Array.from(bySample.keys()).sort();

  const lookup = tables.map((t) => {
    const ids = {};
    Object.values(t.lanes).forEach((L) =>
      Object.keys(L.samples).forEach((sid) => {
        const k = _ttHopNormId(sid);
        if (!(k in ids)) ids[k] = sid;
      }),
    );
    return ids;
  });
  const membership = {};
  names.forEach((s) => {
    membership[s] = {};
    tables.forEach((t, fi) => {
      const cands = _ttHopIdCandidates(s, _ttHopOverride(s, t.name));
      const hit = cands.find((c) => c in lookup[fi]);
      if (hit) membership[s][fi] = lookup[fi][hit];
    });
  });
  const laneRates = tables.map((t) => {
    const o = {};
    Object.keys(t.lanes).forEach((l) => (o[l] = _ttHopLaneRates(t.lanes[l])));
    return o;
  });

  const pairRate = new Map(); // "d\u0001r" -> [rate, source, lane]
  names.forEach((d) =>
    names.forEach((r) => {
      if (d === r) return;
      const md = membership[d],
        mr = membership[r];
      let best = null;
      Object.keys(md).forEach((fi) => {
        if (!(fi in mr)) return;
        const sd = md[fi],
          sr = mr[fi];
        Object.keys(tables[fi].lanes).forEach((lane) => {
          const L = tables[fi].lanes[lane];
          if (!(sd in L.samples) || !(sr in L.samples)) return;
          const [rate, src] = laneRates[fi][lane].rate(sd, sr, P.rate);
          if (!best || rate > best[0]) best = [rate, src, lane];
        });
      });
      if (best) pairRate.set(d + "\u0001" + r, best);
    }),
  );

  const readsAt = {};
  const rowAt = {}; // the row each readsAt entry came from (donor-side notes)
  names.forEach((s) => {
    const m = new Map();
    const rm = new Map();
    bySample.get(s).forEach((r) => {
      const key = String(r["Taxonomic ID #"] || "");
      if (!key) return;
      const k = _ttHopLevel(r) + "\u0001" + key;
      const n = Math.trunc(Number(r["# Reads Aligned"]) || 0);
      if (!m.has(k) || n > m.get(k)) {
        m.set(k, n);
        rm.set(k, r);
      }
    });
    readsAt[s] = m;
    rowAt[s] = rm;
  });
  const sent = new Map(); // donor row -> [{sample, expected, observed}]

  const perRow = new Map();
  const samples = {};
  names.forEach((r) => {
    const fileIdx = Object.keys(membership[r]);
    if (!fileIdx.length) {
      samples[r] = {
        applied: false,
        reason: "sample not found in the index hopping file(s)",
        files: tables.map((t) => t.name),
      };
      return;
    }
    const donors = names.filter((d) => pairRate.has(d + "\u0001" + r));
    const counts = {},
      nNodes = {};
    bySample.get(r).forEach((row) => {
      const key = String(row["Taxonomic ID #"] || "");
      const lvl = _ttHopLevel(row);
      const obs = Math.trunc(Number(row["# Reads Aligned"]) || 0);
      nNodes[lvl] = (nNodes[lvl] || 0) + 1;
      const contrib = [];
      donors.forEach((d) => {
        const nd = readsAt[d].get(lvl + "\u0001" + key) || 0;
        if (nd <= obs || nd <= 0) return;
        const [rate, src, lane] = pairRate.get(d + "\u0001" + r);
        contrib.push({ sample: d, reads: nd, rate, rate_source: src, lane, expected: nd * rate });
      });
      const expected = contrib.reduce((n, c) => n + c.expected, 0);
      const p = expected > 0 ? _ttHopPoissonSf(obs, P.factor * expected) : 0;
      const suspect = expected > 0 && obs > 0 && p >= P.alpha;
      contrib.sort((a, b) => b.expected - a.expected);
      const top = contrib[0] || null;
      perRow.set(row, {
        suspect,
        observed_reads: obs,
        expected_reads: _ttHopRound(expected, 4),
        p_value: _ttHopRound(p, 6),
        explained_fraction: obs > 0 ? _ttHopRound(Math.min(1, expected / obs), 4) : 0,
        top_donor: top ? top.sample : null,
        top_donor_reads: top ? top.reads : null,
        top_donor_rate: top ? top.rate : null,
        rate_source: top ? top.rate_source : null,
      });
      if (suspect) {
        counts[lvl] = (counts[lvl] || 0) + 1;
        contrib.forEach((c) => {
          if (c.expected < TT_HOP_MIN_SOURCE_READS) return;
          const drow = rowAt[c.sample].get(lvl + "\u0001" + key);
          if (!drow) return;
          if (!sent.has(drow)) sent.set(drow, []);
          sent.get(drow).push({ sample: r, expected: c.expected, observed: obs });
        });
      }
    });
    samples[r] = {
      applied: true,
      files: fileIdx.map((fi) => {
        const t = tables[fi],
          sid = membership[r][fi];
        const lanes = Object.keys(t.lanes)
          .filter((l) => sid in t.lanes[l].samples)
          .sort();
        const lane_rates = {};
        lanes.forEach((l) => (lane_rates[l] = laneRates[fi][l].laneRate));
        return {
          file: t.name,
          sample_id: sid,
          lanes,
          lane_rates,
          reads_in_file: lanes.reduce((n, l) => n + t.lanes[l].samples[sid].reads, 0),
        };
      }),
      donor_samples: donors,
      default_rate: P.rate,
      factor: P.factor,
      alpha: P.alpha,
      rate_sources: Array.from(new Set(donors.map((d) => pairRate.get(d + "\u0001" + r)[1]))).sort(),
      n_suspect: counts,
      n_nodes: nNodes,
      computed_in: "browser",
    };
  });
  // Donor side: a light note on the organism the hopped reads likely came FROM.
  const perSource = new Map();
  const nSource = {};
  sent.forEach((list, drow) => {
    list.sort((a, b) => b.expected - a.expected);
    perSource.set(drow, {
      n_recipients: list.length,
      expected_out: _ttHopRound(
        list.reduce((n, x) => n + x.expected, 0),
        4,
      ),
      recipients: list.slice(0, 5).map((x) => ({
        sample: x.sample,
        expected_reads: _ttHopRound(x.expected, 4),
        observed_reads: x.observed,
      })),
    });
    const s = String(drow["Specimen ID"] || "");
    const lvl = _ttHopLevel(drow);
    nSource[s] = nSource[s] || {};
    nSource[s][lvl] = (nSource[s][lvl] || 0) + 1;
  });
  Object.keys(nSource).forEach((s) => {
    if (samples[s] && samples[s].applied) samples[s].n_source = nSource[s];
  });
  return { perRow, perSource, samples };
}

/** "Negative_Control_Miseq (~16.6), S2 (~0.8)" (index_hop_source_text in make_report.py). */
function _ttHopSourceText(src) {
  return ((src && src.recipients) || [])
    .map((r) => `${r.sample} (~${(Number(r.expected_reads) || 0).toFixed(1)})`)
    .join(", ");
}
function _ttHopSourceCols(src) {
  if (!src || !src.recipients || !src.recipients.length) return {};
  return {
    "Index Hop Source Of": _ttHopSourceText(src),
    "Index Hop Reads Out": _ttHopRound(Number(src.expected_out) || 0, 2),
  };
}

/* ── Row markers + hover cards ──────────────────────────────────────────── */
/*  Two small icons after the organism name (Table tab + Summary detections):
      ⇄ (fa-shuffle, orange)  this detection is flagged as a possible index hop
      ↪ (fa-share, faint)     this organism may be the SOURCE of a flagged hop
    Each carries data-hop-* attributes; one delegated listener below shows a
    hover card with the details (showTip / moveTip / hideTip, 02_utilities). */
function _ttHopMarkAttrs(r, kind) {
  const a = (v) => _ttHopEsc(v == null ? "" : v);
  return (
    `data-hop-tip="${kind}" data-hop-s="${a(r["Specimen ID"])}" data-hop-l="${a(r.Level || "Strain")}" ` +
    `data-hop-k="${a(r["Taxonomic ID #"])}" data-hop-n="${a(r["Detected Organism"])}"`
  );
}
function ttHopSuspectBadgeHTML(r) {
  if (!r || r["Index Hop Suspect"] !== true) return "";
  return ` <i class="fas fa-shuffle hop-sus-mark" ${_ttHopMarkAttrs(
    r,
    "sus",
  )} aria-label="possible index hop" tabindex="0"></i>`;
}
/** Small faint marker next to the organism name of a donor row. */
function ttHopSourceBadgeHTML(r) {
  if (!r || !r["Index Hop Source Of"]) return "";
  return ` <i class="fas fa-share hop-src-mark" ${_ttHopMarkAttrs(
    r,
    "src",
  )} aria-label="possible index-hop source" tabindex="0"></i>`;
}
/** Both markers (whichever apply). */
function ttHopBadgesHTML(r) {
  return ttHopSuspectBadgeHTML(r) + ttHopSourceBadgeHTML(r);
}

function _ttHopFindRow(sample, level, key) {
  if (typeof DATA === "undefined") return null;
  return (
    DATA.find(
      (r) =>
        String(r["Specimen ID"]) === sample &&
        String(r.Level || "Strain") === level &&
        String(r["Taxonomic ID #"]) === key,
    ) || null
  );
}
function _ttHopFmt(v, d) {
  const n = Number(v);
  return isFinite(n) ? n.toLocaleString(undefined, { maximumFractionDigits: d == null ? 1 : d }) : "—";
}
const _TT_HOP_RATE_WHY = {
  pair: "measured for this pair of samples from the hopped index combinations in the table",
  lane: "the lane-wide hop fraction from the table (an upper bound for one sample)",
  default: "the default rate, because the table has no hopped rows",
};

/** Hover card for a flagged (⇄) row. */
function ttHopSuspectTipHTML(r) {
  const esc = _ttHopEsc;
  const obs = Math.trunc(Number(r["# Reads Aligned"]) || 0);
  const exp = Number(r["Index Hop Expected Reads"]) || 0;
  const src = r["Index Hop Source"] || "another sample";
  const srcRow = r["Index Hop Source"]
    ? _ttHopFindRow(r["Index Hop Source"], r.Level || "Strain", r["Taxonomic ID #"])
    : null;
  const rs = r["Index Hop Rate Source"] || "";
  const P = TT_HOP.params;
  return (
    `<b style="color:#ffb74d"><i class="fas fa-shuffle"></i> Possible index hop</b><br>` +
    `<span style="font-size:0.85em;color:#ccc">Illumina libraries in the same lane can swap index reads, so a few ` +
    `reads of an organism that is abundant in one sample can show up in another. The reads seen here are about ` +
    `what that would produce.</span>` +
    `<table style="border-collapse:collapse;font-size:0.85em;margin:5px 0">` +
    `<tr><td style="color:#90caf9;padding-right:10px">Reads here</td><td>${_ttHopFmt(obs, 0)} in ${esc(
      r["Specimen ID"],
    )}</td></tr>` +
    `<tr><td style="color:#90caf9;padding-right:10px">Expected from hopping</td><td>~${_ttHopFmt(
      exp,
      1,
    )} reads (${_ttHopFmt(r["Index Hop Explained %"], 1)}% of what is seen)</td></tr>` +
    `<tr><td style="color:#90caf9;padding-right:10px">Main source</td><td>${esc(src)}${
      r["Index Hop Source Reads"] != null
        ? ` — ${_ttHopFmt(r["Index Hop Source Reads"], 0)} reads of this organism`
        : ""
    }${srcRow ? ` (TASS ${_ttHopFmt(srcRow["TASS Score"], 1)})` : ""}</td></tr>` +
    `<tr><td style="color:#90caf9;padding-right:10px">Hop rate</td><td>${_ttHopFmt(r["Index Hop Rate %"], 4)}% — ${esc(
      _TT_HOP_RATE_WHY[rs] || rs,
    )}</td></tr>` +
    `<tr><td style="color:#90caf9;padding-right:10px">P(≥ ${_ttHopFmt(obs, 0)} reads by hopping)</td><td>${_ttHopFmt(
      r["Index Hop p"],
      4,
    )} (flagged when ≥ ${esc(P.alpha)}; expected reads × ${esc(P.factor)} allowed for rate uncertainty)</td></tr>` +
    `</table>` +
    `<span style="font-size:0.85em;color:#ccc">This may be bleed-through rather than an organism present in this ` +
    `specimen. Check whether it is expected here and compare with ${esc(src)}; a real detection usually has far ` +
    `more reads than hopping can explain. Flag only: TASS and the call are unchanged.</span>`
  );
}

/** Hover card for a source (↪) row. */
function ttHopSourceTipHTML(r) {
  const esc = _ttHopEsc;
  const text = String(r["Index Hop Source Of"] || "");
  const recips = [];
  const rx = /(.+?) \(~([0-9.,]+)\)(?:, |$)/g;
  let m;
  while ((m = rx.exec(text))) recips.push({ sample: m[1], expected: Number(m[2].replace(/,/g, "")) });
  const rows = recips
    .map((x) => {
      const rr = _ttHopFindRow(x.sample, r.Level || "Strain", r["Taxonomic ID #"]);
      return (
        `<tr><td style="padding-right:10px">${esc(x.sample)}</td>` +
        `<td style="text-align:right;padding-right:10px">${rr ? _ttHopFmt(rr["# Reads Aligned"], 0) : "—"}</td>` +
        `<td style="text-align:right;padding-right:10px">~${_ttHopFmt(x.expected, 1)}</td>` +
        `<td style="text-align:right">${rr ? _ttHopFmt(rr["Index Hop p"], 4) : "—"}</td></tr>`
      );
    })
    .join("");
  return (
    `<b style="color:#b39ddb"><i class="fas fa-share"></i> Possible index-hop source</b><br>` +
    `<span style="font-size:0.85em;color:#ccc">${esc(r["Detected Organism"])} has ${_ttHopFmt(
      r["# Reads Aligned"],
      0,
    )} reads in ${esc(r["Specimen ID"])}. About ${_ttHopFmt(
      r["Index Hop Reads Out"],
      1,
    )} of them may have hopped into the sample(s) below, where a low-level detection of the same organism is flagged as a possible index hop.</span>` +
    `<table style="border-collapse:collapse;font-size:0.85em;margin:5px 0">` +
    `<tr style="color:#90caf9"><td style="padding-right:10px">Flagged in</td><td style="padding-right:10px">Reads there</td>` +
    `<td style="padding-right:10px">From this sample</td><td>P (hopping)</td></tr>${rows}</table>` +
    `<span style="font-size:0.85em;color:#ccc">Informational only: this detection, its TASS and its call are not affected.</span>`
  );
}

/** Header help for the Index Hop columns (11_tab_table.js renderTableHeaders). */
const TT_HOP_COL_HELP = {
  "Index Hop Suspect":
    "True when index hopping from another sample of the same lane could explain this detection (flag only).",
  "Index Hop Expected Reads":
    "Reads expected to hop in from other samples of the lane that carry more of this organism.",
  "Index Hop Explained %": "Expected hopped reads as a % of the reads seen (capped at 100).",
  "Index Hop p": "P(reads ≥ observed) if all of them came from hopping. Flagged when ≥ alpha (default 0.01).",
  "Index Hop Source": "Sample contributing most of the expected hopped reads.",
  "Index Hop Source Reads": "Reads of this organism in that source sample.",
  "Index Hop Rate %": "Hop rate used for that source: per pair, lane-wide, or the default (see Rate Source).",
  "Index Hop Rate Source":
    "pair = measured for the two samples; lane = lane-wide fraction; default = no hopped rows in the table.",
  "Index Hop Source Of": "Samples this organism's reads may have hopped into (flagged there), with the expected reads.",
  "Index Hop Reads Out": "Total reads of this organism expected to have hopped into flagged detections elsewhere.",
};

(function () {
  if (typeof document === "undefined") return;
  const target = (ev) => (ev.target && ev.target.closest ? ev.target.closest("[data-hop-tip]") : null);
  const html = (el) => {
    const r = _ttHopFindRow(el.dataset.hopS, el.dataset.hopL, el.dataset.hopK);
    if (!r) return "";
    return el.dataset.hopTip === "sus" ? ttHopSuspectTipHTML(r) : ttHopSourceTipHTML(r);
  };
  document.addEventListener("mouseover", (ev) => {
    const el = target(ev);
    if (!el || typeof showTip !== "function") return;
    ev.stopPropagation();
    const h = html(el);
    if (h) showTip(h, ev);
  });
  document.addEventListener("mousemove", (ev) => {
    if (target(ev) && typeof moveTip === "function") moveTip(ev);
  });
  document.addEventListener("mouseout", (ev) => {
    const el = target(ev);
    if (el && !el.contains(ev.relatedTarget) && typeof hideTip === "function") hideTip();
  });
})();

/** A row's "Index Hop …" columns from its block (same rounding as make_report). */
function _ttHopRowCols(b) {
  return {
    "Index Hop Suspect": !!b.suspect,
    "Index Hop Expected Reads": _ttHopRound(Number(b.expected_reads) || 0, 2),
    "Index Hop Explained %": _ttHopRound(100 * (Number(b.explained_fraction) || 0), 1),
    "Index Hop p": _ttHopRound(Number(b.p_value) || 0, 4),
    "Index Hop Source": b.top_donor || "",
    "Index Hop Source Reads": b.top_donor_reads != null ? Math.trunc(b.top_donor_reads) : null,
    "Index Hop Rate %": b.top_donor_rate != null ? _ttHopRound(100 * b.top_donor_rate, 5) : null,
    "Index Hop Rate Source": b.rate_source || "",
  };
}

/* ── apply to DATA + UI ──────────────────────────────────────────────────── */
function _ttHopEnsureColumns() {
  if (typeof ALL_COLS !== "undefined") {
    TT_HOP_COLS.forEach((c) => {
      if (!ALL_COLS.includes(c)) ALL_COLS.push(c);
    });
  }
  if (typeof NUMERIC !== "undefined" && NUMERIC && NUMERIC.add) TT_HOP_NUMERIC.forEach((c) => NUMERIC.add(c));
}

function _ttHopEnsureFlagRule() {
  if (typeof TT_OFLAG_RULES === "undefined" || typeof ttOFlagNewRule !== "function") return false;
  const has = TT_OFLAG_RULES.some((r) => (r.conds || []).some((c) => c.field === "Index Hop Suspect"));
  if (has) return false;
  TT_OFLAG_RULES.push(
    ttOFlagNewRule({
      on: true,
      action: "flag",
      conds: [{ source: "col", field: "Index Hop Suspect", op: "==", value: "true" }],
    }),
  );
  return true;
}

/** Run the test over every loaded row with the tables loaded in this browser.
 *  opts.redraw (default true) refreshes the table, flags and banner. */
function ttHopRecompute(opts) {
  const o = Object.assign({ redraw: true }, opts || {});
  if (!TT_HOP.tables.length || typeof DATA === "undefined") return null;
  const res = ttHopCompute(DATA, TT_HOP.tables, TT_HOP.params);
  DATA.forEach((row) => {
    TT_HOP_COLS.forEach((c) => (row[c] = null));
    const b = res.perRow.get(row);
    if (b) Object.assign(row, _ttHopRowCols(b));
    const src = res.perSource.get(row);
    if (src) Object.assign(row, _ttHopSourceCols(src));
  });
  TT_HOP.result = { samples: res.samples, files: TT_HOP.tables.map((t) => t.name), at: new Date() };
  _ttHopEnsureColumns();
  const anySuspect = DATA.some((r) => r["Index Hop Suspect"] === true);
  if (anySuspect) _ttHopEnsureFlagRule();
  if (typeof ttOFlagInvalidate === "function") ttOFlagInvalidate();
  if (o.redraw) _ttHopRefreshUi();
  return res;
}

function _ttHopRefreshUi() {
  const safe = (fn) => {
    try {
      fn();
    } catch (e) {
      if (window.console && console.error) console.error("index hopping refresh failed", e);
    }
  };
  safe(() => typeof _invalidateFilterCache === "function" && _invalidateFilterCache());
  safe(() => {
    const el = document.getElementById("banner-sub");
    if (el && typeof _buildBannerSub === "function") el.textContent = _buildBannerSub();
  });
  safe(() => {
    if (typeof buildTable === "function") buildTable();
    if (typeof renderTableHeaders === "function") renderTableHeaders();
  });
  safe(() => typeof ttOFlagRefresh === "function" && ttOFlagRefresh());
  safe(() => _ttHopRenderPopover());
}

/** Called by the upload code: add parsed tables and recompute. */
function ttHopAddTables(tables) {
  (tables || []).forEach((t) => {
    if (!t) return;
    const i = TT_HOP.tables.findIndex((x) => x.name === t.name);
    if (i >= 0) TT_HOP.tables[i] = t;
    else TT_HOP.tables.push(t);
  });
  return ttHopRecompute();
}

/** Called after detections are added/removed (uploads, clear): keep a
 *  browser-side result current. The caller redraws. */
function ttHopOnDataChanged() {
  if (TT_HOP.tables.length) {
    ttHopRecompute({ redraw: false });
    return;
  }
  // Pipeline-annotated JSONs dropped on the report carry the columns already:
  // make sure they are listed and the default flag rule is in place.
  if (typeof DATA === "undefined") return;
  if (DATA.some((r) => "Index Hop Suspect" in r)) _ttHopEnsureColumns();
  if (
    DATA.some((r) => r["Index Hop Suspect"] === true) &&
    _ttHopEnsureFlagRule() &&
    typeof ttOFlagInvalidate === "function"
  )
    ttOFlagInvalidate();
}

function ttHopClearTables() {
  TT_HOP.tables = [];
  TT_HOP.result = null;
  // Back to whatever the pipeline wrote (or nothing): drop browser values
  // from rows that came from a pipeline without index hopping.
  const pipelineHas = _ttHopPipelineSamples().length > 0;
  if (!pipelineHas && typeof DATA !== "undefined") DATA.forEach((r) => TT_HOP_COLS.forEach((c) => delete r[c]));
  if (!pipelineHas && typeof ALL_COLS !== "undefined") {
    for (let i = ALL_COLS.length - 1; i >= 0; i--) if (TT_HOP_COLS.includes(ALL_COLS[i])) ALL_COLS.splice(i, 1);
  }
  _ttHopRefreshUi();
}

/* ── banner badge ────────────────────────────────────────────────────────── */
function _ttHopPipelineSamples() {
  const meta = (typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {};
  return Object.keys(meta).filter((s) => meta[s] && meta[s].index_hopping);
}

/** Per-sample status to display: the browser run when there is one, else the
 *  pipeline's metadata. */
function _ttHopStatus() {
  if (TT_HOP.result) return { source: "browser", samples: TT_HOP.result.samples };
  const meta = (typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {};
  const out = {};
  _ttHopPipelineSamples().forEach((s) => (out[s] = meta[s].index_hopping));
  return Object.keys(out).length ? { source: "pipeline", samples: out } : null;
}

function _ttHopFlaggedCounts() {
  const by = {};
  const src = {};
  let strain = 0;
  (typeof DATA !== "undefined" ? DATA : []).forEach((r) => {
    const lvl = _ttHopLevel(r);
    const s = String(r["Specimen ID"] || "");
    if (r["Index Hop Source Of"] && lvl === "strain") src[s] = (src[s] || 0) + 1;
    if (r["Index Hop Suspect"] !== true) return;
    by[s] = by[s] || { strain: 0, species: 0, genus: 0 };
    by[s][lvl] = (by[s][lvl] || 0) + 1;
    if (lvl === "strain") strain++;
  });
  return { by, strain, src };
}

/** Banner pill (#banner-hop). Returns { text } for the long (PDF cover) line. */
function _hopBannerNote(samples) {
  const pill = document.getElementById("banner-hop");
  const label = document.getElementById("banner-hop-label");
  const st = _ttHopStatus();
  const legend = document.getElementById("sum-hop-legend");
  if (legend) {
    const any =
      typeof DATA !== "undefined" && DATA.some((r) => r["Index Hop Suspect"] === true || r["Index Hop Source Of"]);
    legend.hidden = !any;
  }
  if (!st) {
    if (pill) pill.hidden = true;
    _ttHopPopSet(false, false);
    return { text: "" };
  }
  // A pipeline result: start the settings boxes at the run's own values.
  if (st.source === "pipeline" && !TT_HOP.tables.length && !TT_HOP._seeded) {
    const m = Object.values(st.samples).find((x) => x && x.applied) || {};
    ["rate", "factor", "alpha"].forEach((k) => {
      const v = Number(m[k === "rate" ? "default_rate" : k]);
      if (isFinite(v) && v > 0) TT_HOP.params[k] = v;
    });
    TT_HOP._seeded = true;
  }
  const shown = (samples && samples.length ? samples : Object.keys(st.samples)).filter((s) => st.samples[s]);
  const matched = shown.filter((s) => st.samples[s].applied).length;
  const total = samples && samples.length ? samples.length : shown.length;
  const { strain } = _ttHopFlaggedCounts();
  const text = `${strain} detection${
    strain === 1 ? "" : "s"
  } flagged as possible index hops (${matched}/${total} sample${total === 1 ? "" : "s"} in the hopping table)`;
  if (pill) {
    pill.hidden = false;
    pill.classList.toggle("has-flags", strain > 0);
    if (label) label.textContent = `Index hop ${strain} · ${matched}/${total}`;
    pill.setAttribute("aria-label", `Index hopping: ${text}`);
    if (!document.getElementById("banner-hop-pop")) _ttHopInitPopover(pill);
    _ttHopRenderPopover();
  }
  return { text };
}

function _ttHopRenderPopover() {
  const pop = document.getElementById("banner-hop-pop");
  if (!pop) return;
  const st = _ttHopStatus();
  const esc = _ttHopEsc;
  const { by, strain, src } = _ttHopFlaggedCounts();
  const P = TT_HOP.params;
  let html = "";
  if (st) {
    const names = Object.keys(st.samples).sort();
    const matched = names.filter((s) => st.samples[s].applied);
    const unmatched = names.filter((s) => !st.samples[s].applied);
    const fmtRate = (v) => (v == null ? "—" : `${_ttHopRound(100 * v, 4)}%`);
    html +=
      `<div class="bmp-title">Index hopping: ${strain} strain-level detection${strain === 1 ? "" : "s"} flagged</div>` +
      `<div class="bmp-note">${
        st.source === "browser"
          ? `Computed in this browser from <b>${esc(TT_HOP.result.files.join(", "))}</b> over the ${
              names.length
            } loaded sample${names.length === 1 ? "" : "s"}.`
          : "From the pipeline (--index_hopping). Load the table here to change the settings."
      } Flag only: TASS and calls are unchanged.</div>`;
    if (matched.length) {
      const multi = new Set(matched.flatMap((s) => (st.samples[s].files || []).map((f) => f.file))).size > 1;
      const rateTxt = (m) => {
        const srcs = m.rate_sources || [];
        const lr = (m.files || []).flatMap((x) => Object.values(x.lane_rates || {})).filter((v) => v != null);
        return srcs
          .map((src) =>
            src === "default"
              ? `default ${fmtRate(m.default_rate)}`
              : src === "lane"
              ? `lane ${lr.map(fmtRate).join(", ")}`
              : "per pair",
          )
          .join("; ");
      };
      html +=
        `<table><thead><tr><th>Sample</th><th>In table as</th><th>Lane</th><th>Hop rate</th>` +
        `<th class="num">Donors</th><th class="num" title="strain / species / genus rows">Flagged</th>` +
        `<th class="num" title="strain rows whose reads may have hopped into a flagged detection of another sample">Source of</th></tr></thead><tbody>` +
        matched
          .map((s) => {
            const m = st.samples[s];
            const ids = (m.files || []).map((f) => (multi ? `${f.file}: ` : "") + f.sample_id);
            const lanes = (m.files || []).flatMap((x) => x.lanes || []);
            const c = by[s] || { strain: 0, species: 0, genus: 0 };
            return (
              `<tr><td>${esc(s)}</td><td>${esc(ids.join(", "))}</td><td>${esc(lanes.join(", "))}</td>` +
              `<td>${esc(rateTxt(m) || "—")}</td><td class="num">${(m.donor_samples || []).length}</td>` +
              `<td class="num">${c.strain} / ${c.species || 0} / ${c.genus || 0}</td>` +
              `<td class="num">${src[s] || 0}</td></tr>`
            );
          })
          .join("") +
        `</tbody></table>`;
      let pass = null;
      try {
        if (strain && typeof filteredData === "function" && TT_OFLAG_VIEW === "all")
          pass = filteredData().filter((r) => r["Index Hop Suspect"] === true && _ttHopLevel(r) === "strain").length;
      } catch (e) {
        pass = null; // filters not built yet (first paint)
      }
      if (pass != null) {
        html += `<div class="bmp-note">${pass} of the ${strain} flagged strain row${
          strain === 1 ? "" : "s"
        } pass the current filters (TASS cutoff etc.)${
          pass ? "" : ", so <b>Only flagged</b> shows none until you lower the cutoff"
        }.</div>`;
      }
    }
    if (unmatched.length) {
      html += `<div class="bmp-note">Not in the table (not compared): ${unmatched.map(esc).join(", ")}</div>`;
    }
  } else {
    html += `<div class="bmp-title">Index hopping</div>`;
  }
  html +=
    `<div class="bhp-controls">` +
    `<label title="Hop rate used when a table has no hopped rows">Rate % <input type="number" id="bhp-rate" min="0" max="100" step="0.01" value="${_ttHopRound(
      100 * P.rate,
      6,
    )}"></label>` +
    `<label title="Multiplier on the expected hopped reads">Factor <input type="number" id="bhp-factor" min="0.01" step="0.5" value="${P.factor}"></label>` +
    `<label title="Flag when P(reads ≥ observed | hopping) ≥ alpha">Alpha <input type="number" id="bhp-alpha" min="0.0001" max="1" step="0.005" value="${P.alpha}"></label>` +
    `<button type="button" id="bhp-recompute" ${TT_HOP.tables.length ? "" : "disabled"} title="${
      TT_HOP.tables.length ? "Re-run with these settings" : "Load a hopping table first"
    }"><i class="fas fa-rotate"></i> Recompute</button>` +
    `</div><div class="bhp-controls">` +
    `<button type="button" id="bhp-load"><i class="fas fa-file-import"></i> Load hopping table…</button>` +
    `<button type="button" id="bhp-only" ${
      strain ? "" : "disabled"
    }><i class="fas fa-filter"></i> Only flagged</button>` +
    `<button type="button" id="bhp-hide" ${
      strain ? "" : "disabled"
    }><i class="fas fa-eye-slash"></i> Hide flagged</button>` +
    (TT_HOP.tables.length
      ? `<button type="button" id="bhp-clear"><i class="fas fa-xmark"></i> Clear table</button>`
      : "") +
    `<span id="bhp-msg" class="bmp-note"></span></div>` +
    `<div class="bmp-note">CSV, TSV or XLSX in the bcl-convert Index_Hopping_Counts layout (Lane, SampleID, index, index2, # Reads). ` +
    `Single-hop rates are an upper bound for unique dual indexes.</div>`;
  pop.innerHTML = html;
  _ttHopBindPopover(pop);
}

function _ttHopBindPopover(pop) {
  const num = (id, scale) => {
    const v = Number((pop.querySelector("#" + id) || {}).value);
    return isFinite(v) ? v * (scale || 1) : NaN;
  };
  const msg = (t, bad) => {
    const el = pop.querySelector("#bhp-msg");
    if (el) {
      el.textContent = t;
      el.style.color = bad ? "#ef9a9a" : "";
    }
  };
  const rc = pop.querySelector("#bhp-recompute");
  if (rc)
    rc.onclick = () => {
      const rate = num("bhp-rate", 0.01),
        factor = num("bhp-factor"),
        alpha = num("bhp-alpha");
      if (!(rate >= 0 && rate <= 1) || !(factor > 0) || !(alpha > 0 && alpha <= 1)) {
        msg("Rate 0–100 %, factor > 0, alpha in (0, 1].", true);
        return;
      }
      TT_HOP.params = { rate, factor, alpha };
      ttHopRecompute();
    };
  const ld = pop.querySelector("#bhp-load");
  if (ld) ld.onclick = () => document.getElementById("hop-upload-input")?.click();
  const view = (v) => {
    if (typeof TT_OFLAG_VIEW === "undefined") return;
    _ttHopEnsureFlagRule();
    TT_OFLAG_ENABLED = true;
    TT_OFLAG_VIEW = v;
    if (typeof ttOFlagRefresh === "function") ttOFlagRefresh();
  };
  const only = pop.querySelector("#bhp-only");
  if (only) only.onclick = () => view("only");
  const hide = pop.querySelector("#bhp-hide");
  if (hide) hide.onclick = () => view("hide");
  const clr = pop.querySelector("#bhp-clear");
  if (clr) clr.onclick = () => ttHopClearTables();
}

/* Popover behaviour mirrors the masking pill (26_resizable_right_sidebar.js). */
let _ttHopPopPinned = false;
let _ttHopPopHideTimer = null;
function _ttHopPopSet(show, pin) {
  const pill = document.getElementById("banner-hop");
  const pop = document.getElementById("banner-hop-pop");
  if (!pill || !pop) return;
  if (pin !== undefined) _ttHopPopPinned = pin;
  const visible = show || _ttHopPopPinned;
  clearTimeout(_ttHopPopHideTimer);
  pop.classList.toggle("show", visible);
  pill.classList.toggle("open", visible);
  pill.setAttribute("aria-expanded", visible ? "true" : "false");
  if (visible) {
    const r = pill.getBoundingClientRect();
    const w = pop.offsetWidth;
    pop.style.top = `${Math.round(r.bottom + 8)}px`;
    pop.style.left = `${Math.round(Math.max(8, Math.min(r.left, window.innerWidth - w - 8)))}px`;
  }
}
function _ttHopInitPopover(pill) {
  const pop = document.createElement("div");
  pop.id = "banner-hop-pop";
  pop.className = "banner-mask-pop banner-hop-pop";
  pop.setAttribute("role", "dialog");
  document.body.appendChild(pop);
  pill.setAttribute("aria-describedby", "banner-hop-pop");
  pill.setAttribute("role", "button");
  const hideSoon = () => {
    clearTimeout(_ttHopPopHideTimer);
    _ttHopPopHideTimer = setTimeout(() => _ttHopPopSet(false), 150);
  };
  pill.addEventListener("mouseenter", () => _ttHopPopSet(true));
  pill.addEventListener("mouseleave", hideSoon);
  pop.addEventListener("mouseenter", () => clearTimeout(_ttHopPopHideTimer));
  pop.addEventListener("mouseleave", () => {
    // keep it open while an input inside has focus
    if (!pop.contains(document.activeElement)) hideSoon();
  });
  pill.addEventListener("click", (ev) => {
    ev.stopPropagation();
    _ttHopPopSet(!_ttHopPopPinned, !_ttHopPopPinned);
  });
  pill.addEventListener("keydown", (ev) => {
    if (ev.key === "Escape") _ttHopPopSet(false, false);
    if (ev.key === "Enter" || ev.key === " ") {
      ev.preventDefault();
      _ttHopPopSet(!_ttHopPopPinned, !_ttHopPopPinned);
    }
  });
  pop.addEventListener("click", (ev) => ev.stopPropagation());
  document.addEventListener("click", () => _ttHopPopSet(false, false));
  window.addEventListener("resize", () => {
    if (pop.classList.contains("show")) _ttHopPopSet(true);
  });
  return pop;
}

/* ── dedicated file input (badge button + upload-zone link) ──────────────── */
async function ttHopLoadFiles(files) {
  const status = document.getElementById("upload-status");
  const tables = [];
  const errors = [];
  for (const f of Array.from(files || [])) {
    try {
      const t = await ttHopTryParseFile(f, false);
      if (t) tables.push(t);
      else errors.push(`${f.name}: no index-hopping header (need '# Reads' and 'index' or 'SampleID')`);
    } catch (e) {
      errors.push(`${f.name}: ${e && e.message ? e.message : e}`);
    }
  }
  if (tables.length) ttHopAddTables(tables);
  const note = tables.length ? ttHopSummaryText() : "";
  if (status) {
    const div = document.createElement("div");
    div.style.cssText = `margin-top:.3em;font-weight:600;color:${errors.length ? "#b33" : "#6a1b9a"}`;
    div.textContent = [note, ...errors].filter(Boolean).join(" · ");
    status.appendChild(div);
  }
  const msg = document.getElementById("bhp-msg");
  if (msg) msg.textContent = errors.join("; ");
  const inp = document.getElementById("hop-upload-input");
  if (inp) inp.value = "";
}

function ttHopSummaryText() {
  const st = _ttHopStatus();
  if (!st) return "";
  const names = Object.keys(st.samples);
  const matched = names.filter((s) => st.samples[s].applied).length;
  const { strain } = _ttHopFlaggedCounts();
  return `Index hopping: ${matched}/${names.length} sample(s) matched; ${strain} strain-level detection(s) flagged`;
}

(function () {
  if (typeof document === "undefined") return;
  const add = () => {
    if (document.getElementById("hop-upload-input")) return;
    const inp = document.createElement("input");
    inp.type = "file";
    inp.id = "hop-upload-input";
    inp.multiple = true;
    inp.accept = ".csv,.tsv,.txt,.xlsx,.xlsm,.xls";
    inp.style.display = "none";
    inp.addEventListener("change", () => ttHopLoadFiles(inp.files));
    document.body.appendChild(inp);
  };
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", add);
  else add();
})();
