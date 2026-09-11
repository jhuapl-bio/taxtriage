/* ═══════════════════════════════════════════════════════════════════════════
       -  §  EXPORT DATA  (multi-tab combined export)
       -
       -     The per-table "Export Table" control (02_utilities.js) scrapes ONE
       -     rendered <table> at a time, so pulling the Summary detections, the
       -     cross-sample rollup and the VF/AMR hits into one spreadsheet meant
       -     three separate downloads and a manual merge.
       -
       -     This module registers every tab's underlying dataset in ONE catalog
       -     (TT_EXPORT_DATASETS) built from the in-memory model — DATA,
       -     CONTIG_DATA, PROT, NOVELTY, RUN_META, INSILICO_SUITE — not from the
       -     DOM, so a dataset exports in full even when its tab was never opened
       -     and its table is paginated to 25 rows.
       -
       -     The same catalog is mirrored server-side in bin/export_data.py so
       -     `--export_data` produces byte-comparable tables with no browser in
       -     the loop. KEEP THE TWO IN SYNC: dataset ids and column headers are
       -     the contract between them.
       -
       -     Output shapes:
       -       xlsx-multi  one workbook, one sheet per selected dataset
       -       xlsx-wide   one sheet: datasets joined on Specimen ID x Organism
       -       csv-wide    the same wide join as a single CSV
       -       csv-stacked every selected dataset stacked, with a `Dataset` column
═══════════════════════════════════════════════════════════════════════════ */

/* Join keys for the wide shape. A dataset that can be keyed by sample and/or
   organism declares `joinKey` so _ttWideJoin() knows how to line it up against
   the detections backbone. */
const TT_JOIN = {
  NONE: null,
  SAMPLE: "sample", // one row per specimen
  ORG: "organism", // one row per taxid (run-wide rollup)
  SAMPLE_ORG: "sample+organism", // one row per specimen x taxid
};

function _ttNum(v) {
  const n = typeof v === "number" ? v : parseFloat(v);
  return Number.isFinite(n) ? n : null;
}
function _ttRound(v, dp) {
  const n = _ttNum(v);
  if (n === null) return "";
  const f = Math.pow(10, dp == null ? 3 : dp);
  return Math.round(n * f) / f;
}
/** Scalar-only view of a metadata object: nested objects / arrays are dropped
 *  (they have no cell representation) except short primitive arrays, which are
 *  joined with "; ". Mirrors _flatten_scalars() in bin/export_data.py. */
function _ttScalars(obj, prefix) {
  const out = {};
  Object.keys(obj || {}).forEach((k) => {
    const v = obj[k];
    if (v === null || v === undefined) return;
    const key = (prefix || "") + k;
    if (Array.isArray(v)) {
      if (v.length && v.every((x) => x === null || typeof x !== "object")) out[key] = v.join("; ");
      return;
    }
    if (typeof v === "object") return;
    out[key] = typeof v === "boolean" ? (v ? "Yes" : "No") : v;
  });
  return out;
}
/** Union of keys across rows, in first-seen order — the column list for a
 *  dataset whose shape comes from upstream data rather than a fixed schema. */
function _ttKeyUnion(rows) {
  const seen = [];
  const set = new Set();
  (rows || []).forEach((r) => {
    Object.keys(r || {}).forEach((k) => {
      if (!set.has(k)) {
        set.add(k);
        seen.push(k);
      }
    });
  });
  return seen;
}
function _ttHasRows(fn) {
  try {
    const v = fn();
    return Array.isArray(v) ? v.length > 0 : !!v;
  } catch (e) {
    return false;
  }
}

/** Does this row clear its sample's TASS cutoff? Follows the report's own
 *  rollup setting: with ROLLUP_PASS on, a strain rescued by its species or genus
 *  aggregate counts as passing, exactly as filteredData() treats it. */
function _ttRowPasses(r, info) {
  const i = info || (typeof rowPassInfo === "function" ? rowPassInfo(r) : null);
  if (!i) return isTruthy(r["Passes Threshold"]);
  return typeof ROLLUP_PASS !== "undefined" && ROLLUP_PASS ? i.effectivePass : i.strainPass;
}

/* ── Shared derivations ──────────────────────────────────────────────────── */

/** Per-specimen rollup used by the Sample Summary dataset and by the wide
 *  join's sample-level columns. */
function _ttSampleRollup(fd) {
  const bySample = new Map();
  (fd || []).forEach((r) => {
    const s = String(r["Specimen ID"] == null ? "" : r["Specimen ID"]);
    if (!s) return;
    let e = bySample.get(s);
    if (!e) {
      e = { detections: 0, passing: 0, reads: 0, maxTass: 0, topOrg: "", hc: 0, orgs: new Set(), cutoff: null };
      bySample.set(s, e);
    }
    e.detections++;
    e.orgs.add(r["Detected Organism"] || "");
    const t = _ttNum(r["TASS Score"]) || 0;
    if (t > e.maxTass) {
      e.maxTass = t;
      e.topOrg = r["Detected Organism"] || "";
    }
    const info = typeof rowPassInfo === "function" ? rowPassInfo(r) : null;
    if (_ttRowPasses(r, info)) e.passing++;
    if (e.cutoff == null && info && info.thr != null) e.cutoff = info.thr;
    if (isTruthy(r["High Consequence"])) e.hc++;
    e.reads += _ttNum(r["# Reads Aligned"]) || 0;
  });
  return bySample;
}

/** Per sample x organism coverage/contig rollup derived from CONTIG_DATA. */
function _ttCoverageRollup() {
  const out = [];
  (typeof CONTIG_DATA !== "undefined" ? CONTIG_DATA : []).forEach((cd) => {
    const contigs = Array.isArray(cd.contigs) ? cd.contigs : [];
    let len = 0;
    let covered = 0;
    let reads = 0;
    let depthSum = 0;
    contigs.forEach((c) => {
      const l = _ttNum(c.length) || 0;
      len += l;
      covered += _ttNum(c.covered_bases) || 0;
      reads += _ttNum(c.reads) || 0;
      depthSum += (_ttNum(c.mean_depth) || 0) * l;
    });
    const dh = cd.depth_histogram || {};
    out.push({
      "Specimen ID": cd.sample || "",
      "Detected Organism": cd.organism || "",
      "Taxonomic ID #": String(cd.taxon_id == null ? "" : cd.taxon_id),
      "Contigs": contigs.length,
      "Genome Length (bp)": len || "",
      "Covered Bases": covered || "",
      "Breadth %": len ? _ttRound((covered / len) * 100, 3) : "",
      "Mean Depth": len ? _ttRound(depthSum / len, 3) : "",
      "# Reads Aligned": reads || "",
      "Bases 0x": _ttNum(dh["0x"]) == null ? "" : dh["0x"],
      "Bases 1-5x": _ttNum(dh["1-5x"]) == null ? "" : dh["1-5x"],
      "Bases 5-10x": _ttNum(dh["5-10x"]) == null ? "" : dh["5-10x"],
      "Bases 10-50x": _ttNum(dh["10-50x"]) == null ? "" : dh["10-50x"],
      "Bases >50x": _ttNum(dh[">50x"]) == null ? "" : dh[">50x"],
    });
  });
  return out;
}

function _ttNoveltyRows(kind) {
  const rows = [];
  const samples = (typeof NOVELTY !== "undefined" && NOVELTY && NOVELTY.samples) || {};
  Object.keys(samples).forEach((sname) => {
    const sdata = samples[sname] || {};
    if (kind === "summary") {
      const s = sdata.summary || {};
      if (!Object.keys(s).length) return;
      rows.push(Object.assign({ "Specimen ID": sname }, _ttScalars(s)));
    } else {
      (sdata.candidates || []).forEach((c) => {
        const row = Object.assign({ "Specimen ID": sname }, _ttScalars(c));
        const p = c.pathogen;
        if (p && typeof p === "object") {
          row["Pathogen Match"] = p.name || p.organism || "Yes";
          if (p.taxid != null) row["Pathogen Taxid"] = p.taxid;
        }
        rows.push(row);
      });
    }
  });
  return rows;
}

function _ttInsilicoRows(kind) {
  const suite = typeof INSILICO_SUITE !== "undefined" ? INSILICO_SUITE : null;
  if (!suite || !Array.isArray(suite.groups)) return [];
  const rows = [];
  suite.groups.forEach((g) => {
    const head = {
      "Parent Sample": g.parent || "",
      "Platform": g.platform || "",
      "Series Kind": g.series_kind || "depth",
      "Level": g.level || "",
      "Read Unit": g.read_unit || "reads",
    };
    if (kind === "datasets") {
      (g.datasets || []).forEach((d) => {
        rows.push(
          Object.assign({}, head, {
            "Dataset ID": d.id || "",
            "Replicate": d.replicate == null ? "" : d.replicate,
            "Target Count": d.target_count == null ? "" : d.target_count,
            "Actual Count": d.actual_count == null ? "" : d.actual_count,
            "Total Master Reads": d.total_master_reads == null ? "" : d.total_master_reads,
            "Seed": d.seed == null ? "" : d.seed,
            "Observed Total Reads": d.observed_total_reads == null ? "" : d.observed_total_reads,
            "# Detected": d.n_detected == null ? "" : d.n_detected,
            "TP": d.tp == null ? "" : d.tp,
            "FP": d.fp == null ? "" : d.fp,
            "FN": d.fn == null ? "" : d.fn,
            "Precision": d.precision == null ? "" : d.precision,
            "Recall": d.recall == null ? "" : d.recall,
            "F1": d.f1 == null ? "" : d.f1,
          }),
        );
      });
    } else {
      (g.organisms || []).forEach((o) => {
        (o.series || []).forEach((s) => {
          rows.push(
            Object.assign({}, head, {
              "Taxonomic ID #": String(o.taxid == null ? "" : o.taxid),
              "Detected Organism": o.name || "",
              "Microbial Category": o.category || "",
              "Species Name": o.species || "",
              "Genus Name": o.genus || "",
              "Expected Fraction": o.expected_fraction == null ? "" : o.expected_fraction,
              "LoD Count": o.lod_count == null ? "" : o.lod_count,
              "Series Count": s.count == null ? "" : s.count,
              "Expected Reads": s.expected_reads == null ? "" : s.expected_reads,
              "Observed Reads": s.observed_reads == null ? "" : s.observed_reads,
              "TASS Score": s.tass == null ? "" : s.tass,
              "Detection Rate": s.detection_rate == null ? "" : s.detection_rate,
              "Detected": s.detected ? "Yes" : "No",
              "# Replicates": s.n_reps == null ? "" : s.n_reps,
            }),
          );
        });
      });
    }
  });
  return rows;
}

/* ── Metadata access ─────────────────────────────────────────────────────────
   Run metadata lives in two places: RUN_META (the per-sample collection /
   geography record, including anything the user uploaded on the Metadata tab)
   and SAMPLE_META (pipeline provenance). Merge them per sample so a pivot or a
   long export can key on "host_disease" or "location" without caring which feed
   carried it. */

/** Pipeline bookkeeping that would only add noise to a metadata pivot. */
const _TT_META_SKIP = new Set([
  "sample_name",
  "weights",
  "best_cutoffs",
  "best_cutoffs_by_domain",
  "best_cutoffs_source",
  "missing_positive_controls",
  "missing_insilico_controls",
  "missing_insilico_by_type",
  "negative_controls_used",
  "positive_controls_used",
  "insilico_controls_used",
  "insilico_simulator_types",
  "num_keys",
  "num_subkeys",
  "num_toplevelkeys",
  "num_species_groups",
  "commit_id",
  "workflow_revision",
]);

let _TT_META_INDEX = null;
let _TT_META_INDEX_SIG = null;

function _ttMetaIndex() {
  const sig =
    (typeof RUN_META !== "undefined" ? RUN_META.length : 0) +
    "|" +
    Object.keys((typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {}).length +
    "|" +
    (typeof SAMPLE_META_EPOCH !== "undefined" ? SAMPLE_META_EPOCH : 0);
  if (_TT_META_INDEX && _TT_META_INDEX_SIG === sig) return _TT_META_INDEX;
  const idx = new Map();
  const put = (sample, obj) => {
    if (!sample) return;
    const cur = idx.get(sample) || {};
    Object.keys(obj || {}).forEach((k) => {
      if (_TT_META_SKIP.has(k)) return;
      const v = obj[k];
      if (v === null || v === undefined || v === "") return;
      if (typeof v === "object" && !Array.isArray(v)) return;
      if (Array.isArray(v) && v.some((x) => x && typeof x === "object")) return;
      // RUN_META wins over SAMPLE_META: it is what the Metadata tab edits.
      if (cur[k] === undefined) cur[k] = Array.isArray(v) ? v.join("; ") : v;
    });
    idx.set(sample, cur);
  };
  (typeof RUN_META !== "undefined" ? RUN_META : []).forEach((m) => put(m.sample_name || m.sample_id, m));
  const sm = (typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {};
  Object.keys(sm).forEach((s) => put(s, sm[s]));
  _TT_META_INDEX = idx;
  _TT_META_INDEX_SIG = sig;
  return idx;
}

function _ttMetaFor(sample) {
  return _ttMetaIndex().get(String(sample == null ? "" : sample)) || {};
}

/** Metadata fields worth offering as a pivot axis, with how many distinct
 *  non-empty values each has — a field with one value makes a pointless column
 *  axis, so the picker shows the count and sorts the useful ones first. */
/** Numeric measurements that make fine metadata COLUMNS but useless pivot axes
 *  (a crosstab keyed on total_reads has one column per sample). */
const _TT_PIVOT_FIELD_SKIP = new Set([
  "total_reads",
  "aligned_reads",
  "total_organism_reads",
  "latitude",
  "longitude",
  "minmapq",
  "mapq_breadth_power",
  "mapq_gini_power",
  "control_fold_threshold",
  "min_conf_applied",
  "depth",
  "salinity",
]);

function _ttMetaFields(forPivot) {
  const counts = new Map();
  _ttMetaIndex().forEach((meta) => {
    Object.keys(meta).forEach((k) => {
      const v = meta[k];
      if (v === null || v === undefined || String(v).trim() === "") return;
      if (!counts.has(k)) counts.set(k, new Set());
      counts.get(k).add(String(v));
    });
  });
  const nSamples = Math.max(1, _ttMetaIndex().size);
  return [...counts.entries()]
    .map(([field, vals]) => ({ field: field, n: vals.size, vals: vals, label: _ttPrettyField(field) }))
    .filter((f) => {
      if (!f.n) return false;
      if (!forPivot) return true;
      if (_TT_PIVOT_FIELD_SKIP.has(f.field)) return false;
      // Continuous numerics give a crosstab one column per sample.
      const numeric = [...f.vals].every((v) => _ttNum(v) !== null);
      return !(numeric && f.n > 8 && f.n > 0.6 * nSamples);
    })
    .sort((a, b) => (b.n > 1) - (a.n > 1) || a.label.localeCompare(b.label));
}

function _ttPrettyField(k) {
  return String(k)
    .replace(/[_\-]+/g, " ")
    .replace(/\b\w/g, (c) => c.toUpperCase())
    .replace(/\bId\b/, "ID")
    .replace(/\bTass\b/, "TASS");
}

/* ── Pivot / crosstab ────────────────────────────────────────────────────────
   "How many hits to Influenza A came from each collection site / host type?" is
   a question no single tab answers: the counts live in the detections, the site
   lives in the run metadata. The pivot joins the two and aggregates, in either
   layout:

     wide  one row per organism, one column per metadata value (a crosstab)
     long  one row per organism x metadata value (tidy, for R / pandas / plots)
*/

const TT_PIVOT_ROWS = [
  { id: "organism", label: "Detected Organism", of: (r) => r["Detected Organism"] || "" },
  {
    id: "organism_taxid",
    label: "Organism + Taxid",
    of: (r) => (r["Detected Organism"] || "") + (r["Taxonomic ID #"] ? " (" + r["Taxonomic ID #"] + ")" : ""),
  },
  { id: "genus", label: "Genus", of: (r) => r["Genus Name"] || r["Genus"] || "" },
  { id: "category", label: "Microbial Category", of: (r) => r["Microbial Category"] || "Unknown" },
  { id: "domain", label: "Domain", of: (r) => r["Domain"] || r["Kingdom"] || "Unknown" },
  { id: "sample", label: "Specimen ID", of: (r) => r["Specimen ID"] || "" },
  { id: "sample_type", label: "Sample Type", of: (r) => r["Sample Type"] || "" },
];

const TT_PIVOT_MEASURES = [
  { id: "detections", label: "# Detections" },
  { id: "specimens", label: "# Specimens" },
  { id: "organisms", label: "# Distinct Organisms" },
  { id: "reads", label: "Total Reads Aligned" },
  { id: "mean_tass", label: "Mean TASS" },
  { id: "max_tass", label: "Max TASS" },
];

function _ttPivotRowDim(id) {
  return TT_PIVOT_ROWS.find((d) => d.id === id) || TT_PIVOT_ROWS[0];
}
function _ttPivotMeasure(id) {
  return TT_PIVOT_MEASURES.find((m) => m.id === id) || TT_PIVOT_MEASURES[0];
}

function _ttPivotCellInit() {
  return { n: 0, specimens: new Set(), organisms: new Set(), reads: 0, tassSum: 0, tassN: 0, tassMax: null };
}

function _ttPivotAccumulate(cell, r) {
  cell.n++;
  if (r["Specimen ID"]) cell.specimens.add(String(r["Specimen ID"]));
  const org = String(r["Taxonomic ID #"] || r["Detected Organism"] || "");
  if (org) cell.organisms.add(org);
  cell.reads += _ttNum(r["# Reads Aligned"]) || 0;
  const t = _ttNum(r["TASS Score"]);
  if (t !== null) {
    cell.tassSum += t;
    cell.tassN++;
    cell.tassMax = cell.tassMax === null ? t : Math.max(cell.tassMax, t);
  }
}

function _ttPivotValue(cell, measure) {
  if (!cell) return measure === "mean_tass" || measure === "max_tass" ? "" : 0;
  switch (measure) {
    case "specimens":
      return cell.specimens.size;
    case "organisms":
      return cell.organisms.size;
    case "reads":
      return cell.reads;
    case "mean_tass":
      return cell.tassN ? _ttRound(cell.tassSum / cell.tassN, 2) : "";
    case "max_tass":
      return cell.tassMax === null ? "" : _ttRound(cell.tassMax, 2);
    default:
      return cell.n;
  }
}

/** Build the pivot. opts: {rowDim, field, measure, shape, dropEmpty}
 *  `field` "" means no column axis — a plain per-row rollup. */
function _ttBuildPivot(ctx, opts) {
  const dim = _ttPivotRowDim(opts.rowDim);
  const measure = _ttPivotMeasure(opts.measure);
  const field = opts.field || "";
  const fieldLabel = field ? _ttPrettyField(field) : "";
  const rows = ctx.filtered ? ctx.fd : ctx.rowsAll;

  const grid = new Map(); // rowKey -> Map(colKey -> cell)
  const rowTotals = new Map(); // rowKey -> cell (the "All" column)
  const colKeys = new Map(); // colKey -> count, for ordering

  rows.forEach((r) => {
    const rk = String(dim.of(r) || "").trim();
    if (!rk) return;
    const meta = field ? _ttMetaFor(r["Specimen ID"]) : null;
    let ck = "";
    if (field) {
      const v = meta[field];
      ck = v === null || v === undefined || String(v).trim() === "" ? "(not recorded)" : String(v).trim();
    }
    if (!grid.has(rk)) grid.set(rk, new Map());
    const byCol = grid.get(rk);
    if (!byCol.has(ck)) byCol.set(ck, _ttPivotCellInit());
    _ttPivotAccumulate(byCol.get(ck), r);
    if (!rowTotals.has(rk)) rowTotals.set(rk, _ttPivotCellInit());
    _ttPivotAccumulate(rowTotals.get(rk), r);
    colKeys.set(ck, (colKeys.get(ck) || 0) + 1);
  });

  // Columns: most-populated first, with "(not recorded)" always last so a
  // sparsely annotated run does not lead with its gaps.
  const cols = [...colKeys.keys()].sort((a, b) => {
    if (a === "(not recorded)") return 1;
    if (b === "(not recorded)") return -1;
    return colKeys.get(b) - colKeys.get(a) || a.localeCompare(b, undefined, { numeric: true });
  });

  // Rows ordered by the measure's own total, so the biggest signal is on top.
  const rowKeys = [...grid.keys()].sort((a, b) => {
    const va = _ttPivotValue(rowTotals.get(a), measure.id);
    const vb = _ttPivotValue(rowTotals.get(b), measure.id);
    const na = typeof va === "number" ? va : -1;
    const nb = typeof vb === "number" ? vb : -1;
    return nb - na || String(a).localeCompare(String(b));
  });

  if (opts.shape === "long") {
    const columns = [dim.label, fieldLabel || "Group", "Measure", "Value"];
    const out = [];
    rowKeys.forEach((rk) => {
      const byCol = grid.get(rk);
      (field ? cols : [""]).forEach((ck) => {
        const cell = byCol.get(ck);
        if (!cell && opts.dropEmpty !== false) return;
        out.push({
          [dim.label]: rk,
          [fieldLabel || "Group"]: field ? ck : "All",
          Measure: measure.label,
          Value: _ttPivotValue(cell, measure.id),
        });
      });
    });
    return { columns: columns, rows: out, meta: { dim: dim, measure: measure, field: field, nCols: cols.length } };
  }

  // wide
  const columns = [dim.label].concat(field ? cols : []).concat(["Total (" + measure.label + ")"]);
  const out = rowKeys.map((rk) => {
    const byCol = grid.get(rk);
    const row = {};
    row[dim.label] = rk;
    if (field) cols.forEach((ck) => (row[ck] = _ttPivotValue(byCol.get(ck), measure.id)));
    row["Total (" + measure.label + ")"] = _ttPivotValue(rowTotals.get(rk), measure.id);
    return row;
  });
  return { columns: columns, rows: out, meta: { dim: dim, measure: measure, field: field, nCols: cols.length } };
}

/* ── The catalog ─────────────────────────────────────────────────────────────
   Each entry:
     id        stable identifier — also the --export_data_datasets token and the
               sheet-name seed. Must match bin/export_data.py.
     label     sheet / section name as the user sees it
     tab       which report tab the data backs (grouping in the picker)
     desc      one-line explanation under the checkbox
     joinKey   how the wide shape lines this dataset up (see TT_JOIN)
     filterable  true when the rows honour the sidebar filters
     available()  hide the entry when the run carries no such data
     build(ctx)   -> { columns: [...], rows: [{col: value}] }
*/
const TT_EXPORT_DATASETS = [
  {
    id: "detections",
    label: "Detections",
    tab: "Summary / Table",
    desc: "Every organism detection with all TASS, coverage and taxonomy columns.",
    joinKey: TT_JOIN.SAMPLE_ORG,
    filterable: true,
    defaultOn: true,
    available: () => _ttHasRows(() => DATA),
    build(ctx) {
      const cols = (typeof ALL_COLS !== "undefined" && ALL_COLS.length ? ALL_COLS : _ttKeyUnion(ctx.rowsAll)).slice();
      const src = ctx.filtered ? ctx.fd : ctx.rowsAll;
      // `Passes Threshold` is computed in the browser (the JSON carries False for
      // every record), so spell out the verdict and the cutoff it used — the
      // same two columns bin/export_data.py appends.
      const extra = ["TASS Cutoff", "Passes Cutoff"];
      const rows = src.map((r) => {
        const o = {};
        cols.forEach((c) => {
          const v = r[c];
          o[c] = typeof v === "boolean" ? (v ? "Yes" : "No") : v == null ? "" : v;
        });
        const info = typeof rowPassInfo === "function" ? rowPassInfo(r) : null;
        o["TASS Cutoff"] = info && info.thr != null ? _ttRound(info.thr, 2) : "";
        o["Passes Cutoff"] = _ttRowPasses(r, info) ? "Yes" : "No";
        return o;
      });
      return { columns: cols.concat(extra), rows: rows };
    },
  },
  {
    id: "detections_meta",
    label: "Detections + Metadata",
    tab: "Summary / Table",
    desc: "Every detection with its sample's run metadata joined on — long format, ready to pivot.",
    joinKey: TT_JOIN.SAMPLE_ORG,
    filterable: true,
    defaultOn: false,
    available: () => _ttHasRows(() => DATA) && _ttMetaFields().length > 0,
    build(ctx) {
      // The detections table plus every metadata column the run carries, one row
      // per detection. This is the shape you hand to a pivot table, R or pandas
      // to ask "how many hits to this organism came from each site / host?".
      const base = _ttDatasetById("detections").build(ctx);
      const fields = _ttMetaFields().map((f) => f.field);
      const cols = base.columns.concat(fields.map((f) => _ttPrettyField(f)));
      const src = ctx.filtered ? ctx.fd : ctx.rowsAll;
      const rows = base.rows.map((row, i) => {
        const meta = _ttMetaFor(src[i] && src[i]["Specimen ID"]);
        const out = Object.assign({}, row);
        fields.forEach((f) => {
          const v = meta[f];
          out[_ttPrettyField(f)] = v === null || v === undefined ? "" : v;
        });
        return out;
      });
      return { columns: cols, rows: rows };
    },
  },
  {
    id: "sample_summary",
    label: "Sample Summary",
    tab: "Summary",
    desc: "One row per specimen: detection counts, aligned reads, strongest hit.",
    joinKey: TT_JOIN.SAMPLE,
    filterable: true,
    defaultOn: true,
    available: () => _ttHasRows(() => DATA),
    build(ctx) {
      const roll = _ttSampleRollup(ctx.filtered ? ctx.fd : ctx.rowsAll);
      const cols = [
        "Specimen ID",
        "Specimen Group",
        "Sample Type",
        "Platform",
        "Total Reads",
        "Aligned Reads",
        "TASS Cutoff",
        "# Detections",
        "# Passing Cutoff",
        "# Distinct Organisms",
        "# High Consequence",
        "Max TASS Score",
        "Top Organism",
        "QC Flag",
      ];
      const rows = [];
      roll.forEach((e, sample) => {
        const meta = (typeof SAMPLE_META !== "undefined" && SAMPLE_META[sample]) || {};
        // Whole-sample QC verdict (41_sample_flags.js). Absent on reports built
        // before sample QC existed, hence the capability check.
        let flag = "";
        try {
          if (typeof ttFlagStateFor === "function") {
            const st = ttFlagStateFor(sample);
            if (st && st.flagged) {
              flag = (st.hits || []).map((h) => h.text).join(" | ") || "flagged";
              if (st.hide) flag = "hidden: " + flag;
            }
          }
        } catch (e) {
          flag = "";
        }
        rows.push({
          "Specimen ID": sample,
          "Specimen Group": typeof specimenOf === "function" ? specimenOf(sample) : sample,
          "Sample Type": meta.sample_type || "",
          "Platform": meta.platform || "",
          "Total Reads": meta.total_reads == null ? "" : meta.total_reads,
          "Aligned Reads": meta.aligned_reads == null ? "" : meta.aligned_reads,
          "TASS Cutoff": e.cutoff == null ? "" : _ttRound(e.cutoff, 2),
          "# Detections": e.detections,
          "# Passing Cutoff": e.passing,
          "# Distinct Organisms": e.orgs.size,
          "# High Consequence": e.hc,
          "Max TASS Score": _ttRound(e.maxTass, 3),
          "Top Organism": e.topOrg,
          "QC Flag": flag,
        });
      });
      rows.sort((a, b) => String(a["Specimen ID"]).localeCompare(String(b["Specimen ID"])));
      return { columns: cols, rows: rows };
    },
  },
  {
    id: "organism_summary",
    label: "Cross-Sample Organisms",
    tab: "Explore",
    desc: "Run-wide rollup per organism: prevalence plus TASS / coverage spread.",
    joinKey: TT_JOIN.ORG,
    filterable: true,
    defaultOn: true,
    available: () => typeof _xsAggregate === "function" && _ttHasRows(() => DATA),
    build(ctx) {
      const agg = _xsAggregate(ctx.filtered ? ctx.fd : ctx.rowsAll);
      const cols = [
        "Taxonomic ID #",
        "Detected Organism",
        "Microbial Category",
        "High Consequence",
        "# Specimens Passing",
        "# Specimens Below Cutoff",
        "# Specimens Detected",
        "# Specimens Total",
        "Prevalence %",
        "Detected Prevalence %",
        "Mean TASS",
        "Median TASS",
        "Min TASS",
        "Max TASS",
        "Mean Coverage",
        "Median Coverage",
        "Min Coverage",
        "Max Coverage",
        "Total Reads Aligned",
        "ANI Group",
        "ANI Group Size",
        "Specimens",
      ];
      const rows = agg.rows.map((r) => ({
        "Taxonomic ID #": r.taxid,
        "Detected Organism": r.name,
        "Microbial Category": r.cat,
        "High Consequence": r.hc ? "Yes" : "No",
        "# Specimens Passing": r.passCount,
        "# Specimens Below Cutoff": r.belowCount,
        "# Specimens Detected": r.detCount,
        "# Specimens Total": r.total,
        "Prevalence %": _ttRound(r.samplePct, 2),
        "Detected Prevalence %": _ttRound(r.detPct, 2),
        "Mean TASS": _ttRound(r.meanTass, 3),
        "Median TASS": _ttRound(r.medianTass, 3),
        "Min TASS": _ttRound(r.minTass, 3),
        "Max TASS": _ttRound(r.maxTass, 3),
        "Mean Coverage": _ttRound(r.meanCov, 3),
        "Median Coverage": _ttRound(r.medianCov, 3),
        "Min Coverage": _ttRound(r.minCov, 3),
        "Max Coverage": _ttRound(r.maxCov, 3),
        "Total Reads Aligned": r.reads,
        "ANI Group": r.aniGroup == null ? "" : r.aniGroup,
        "ANI Group Size": r.aniGroupSize == null ? "" : r.aniGroupSize,
        "Specimens": [...r.samples].sort().join("; "),
      }));
      rows.sort((a, b) => b["# Specimens Passing"] - a["# Specimens Passing"]);
      return { columns: cols, rows: rows };
    },
  },
  {
    id: "coverage",
    label: "Coverage Summary",
    tab: "Coverage / Histogram",
    desc: "Per specimen x organism breadth, mean depth and depth-bin base counts.",
    joinKey: TT_JOIN.SAMPLE_ORG,
    filterable: true,
    defaultOn: false,
    available: () => _ttHasRows(() => (typeof CONTIG_DATA !== "undefined" ? CONTIG_DATA : [])),
    build(ctx) {
      let rows = _ttCoverageRollup();
      if (ctx.filtered) {
        const keep = new Set(ctx.fd.map((r) => r["Specimen ID"] + "\u0001" + (r["Taxonomic ID #"] || "")));
        const samples = new Set(ctx.fd.map((r) => r["Specimen ID"]));
        rows = rows.filter(
          (r) => keep.has(r["Specimen ID"] + "\u0001" + r["Taxonomic ID #"]) || samples.has(r["Specimen ID"]),
        );
      }
      return { columns: _ttKeyUnion(rows), rows: rows };
    },
  },
  {
    id: "contigs",
    label: "Coverage by Contig",
    tab: "Coverage / Histogram",
    desc: "One row per aligned contig / reference sequence.",
    joinKey: TT_JOIN.NONE,
    filterable: true,
    defaultOn: false,
    available: () =>
      _ttHasRows(() =>
        (typeof CONTIG_DATA !== "undefined" ? CONTIG_DATA : []).some((c) => (c.contigs || []).length),
      ),
    build(ctx) {
      const samples = ctx.filtered ? new Set(ctx.fd.map((r) => r["Specimen ID"])) : null;
      const rows = [];
      (typeof CONTIG_DATA !== "undefined" ? CONTIG_DATA : []).forEach((cd) => {
        if (samples && !samples.has(cd.sample)) return;
        (cd.contigs || []).forEach((c) => {
          const dh = c.depth_histogram || {};
          rows.push({
            "Specimen ID": cd.sample || "",
            "Detected Organism": cd.organism || "",
            "Taxonomic ID #": String(cd.taxon_id == null ? "" : cd.taxon_id),
            "Contig": c.name || "",
            "Length (bp)": c.length == null ? "" : c.length,
            "# Reads Aligned": c.reads == null ? "" : c.reads,
            "Mean Depth": c.mean_depth == null ? "" : c.mean_depth,
            "Covered Bases": c.covered_bases == null ? "" : c.covered_bases,
            "Coverage": c.coverage == null ? "" : c.coverage,
            "Bases 0x": dh["0x"] == null ? "" : dh["0x"],
            "Bases 1-5x": dh["1-5x"] == null ? "" : dh["1-5x"],
            "Bases 5-10x": dh["5-10x"] == null ? "" : dh["5-10x"],
            "Bases 10-50x": dh["10-50x"] == null ? "" : dh["10-50x"],
            "Bases >50x": dh[">50x"] == null ? "" : dh[">50x"],
          });
        });
      });
      return { columns: _ttKeyUnion(rows), rows: rows };
    },
  },
  {
    id: "vfamr_hits",
    label: "VF/AMR Per-Gene Hits",
    tab: "VF/AMR",
    desc: "Every virulence-factor / resistance-gene hit with its source and identity.",
    joinKey: TT_JOIN.NONE,
    filterable: false,
    defaultOn: true,
    available: () => _ttHasRows(() => (typeof PROT !== "undefined" && PROT.per_gene_hits) || []),
    build() {
      const rows = ((typeof PROT !== "undefined" && PROT.per_gene_hits) || []).slice();
      return { columns: _ttKeyUnion(rows), rows: rows };
    },
  },
  {
    id: "vfamr_genus",
    label: "VF/AMR Genus Summary",
    tab: "VF/AMR",
    desc: "Per-genus counts of annotated virulence / resistance properties.",
    joinKey: TT_JOIN.NONE,
    filterable: false,
    defaultOn: true,
    available: () => _ttHasRows(() => (typeof PROT !== "undefined" && PROT.genus_summary) || []),
    build() {
      const rows = ((typeof PROT !== "undefined" && PROT.genus_summary) || []).slice();
      return { columns: _ttKeyUnion(rows), rows: rows };
    },
  },
  {
    id: "vfamr_amr",
    label: "AMR Genes",
    tab: "VF/AMR",
    desc: "Resistance genes with antibiotic class annotation.",
    joinKey: TT_JOIN.NONE,
    filterable: false,
    defaultOn: true,
    available: () => _ttHasRows(() => (typeof PROT !== "undefined" && PROT.amr_genes) || []),
    build() {
      const rows = ((typeof PROT !== "undefined" && PROT.amr_genes) || []).slice();
      return { columns: _ttKeyUnion(rows), rows: rows };
    },
  },
  {
    id: "novelty_summary",
    label: "Novelty Summary",
    tab: "Novelty",
    desc: "Per-sample reference-free novelty statistics.",
    joinKey: TT_JOIN.SAMPLE,
    filterable: false,
    defaultOn: true,
    available: () => _ttHasRows(() => _ttNoveltyRows("summary")),
    build() {
      const rows = _ttNoveltyRows("summary");
      return { columns: _ttKeyUnion(rows), rows: rows };
    },
  },
  {
    id: "novelty_candidates",
    label: "Novelty Candidates",
    tab: "Novelty",
    desc: "Candidate novel taxa per sample, with pathogen cross-references.",
    joinKey: TT_JOIN.NONE,
    filterable: false,
    defaultOn: true,
    available: () => _ttHasRows(() => _ttNoveltyRows("candidates")),
    build() {
      const rows = _ttNoveltyRows("candidates");
      return { columns: _ttKeyUnion(rows), rows: rows };
    },
  },
  {
    id: "run_metadata",
    label: "Run Metadata",
    tab: "Run Metadata",
    desc: "Per-sample run / collection metadata, including any uploaded columns.",
    joinKey: TT_JOIN.SAMPLE,
    filterable: false,
    defaultOn: true,
    available: () => _ttHasRows(() => (typeof RUN_META !== "undefined" ? RUN_META : [])),
    build() {
      const rows = (typeof RUN_META !== "undefined" ? RUN_META : []).map((m) => _ttScalars(m));
      return { columns: _ttKeyUnion(rows), rows: rows };
    },
  },
  {
    id: "sample_metadata",
    label: "Pipeline Sample Metadata",
    tab: "Run Metadata",
    desc: "Per-sample pipeline provenance: read counts, thresholds, controls.",
    joinKey: TT_JOIN.SAMPLE,
    filterable: false,
    defaultOn: false,
    available: () => _ttHasRows(() => Object.keys((typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {})),
    build() {
      const meta = (typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {};
      const rows = Object.keys(meta)
        .sort()
        .map((s) => Object.assign({ "Specimen ID": s }, _ttScalars(meta[s])));
      return { columns: _ttKeyUnion(rows), rows: rows };
    },
  },
  {
    id: "geo",
    label: "Sample Geography",
    tab: "Map",
    desc: "Coordinates and place fields for samples that carry them.",
    joinKey: TT_JOIN.SAMPLE,
    filterable: false,
    defaultOn: false,
    available: () =>
      _ttHasRows(() =>
        (typeof RUN_META !== "undefined" ? RUN_META : []).filter(
          (m) => _ttNum(m.latitude) !== null && _ttNum(m.longitude) !== null,
        ),
      ),
    build() {
      const rows = (typeof RUN_META !== "undefined" ? RUN_META : [])
        .filter((m) => _ttNum(m.latitude) !== null && _ttNum(m.longitude) !== null)
        .map((m) => ({
          "Specimen ID": m.sample_name || m.sample_id || "",
          "Latitude": _ttNum(m.latitude),
          "Longitude": _ttNum(m.longitude),
          "Location": m.location || "",
          "Country": m.sample_origin_country || "",
          "State/Province": m.sample_origin_state_province_territory || "",
          "Environmental Site": m.environmental_site || "",
          "Collection Time": m.collection_time || "",
          "Run ID": m.run_id || "",
        }));
      return { columns: _ttKeyUnion(rows), rows: rows };
    },
  },
  {
    id: "insilico_datasets",
    label: "In-Silico Datasets",
    tab: "In-Silico",
    desc: "Per subsample dataset: target vs actual reads and detection scoring.",
    joinKey: TT_JOIN.NONE,
    filterable: false,
    defaultOn: true,
    available: () => _ttHasRows(() => _ttInsilicoRows("datasets")),
    build() {
      const rows = _ttInsilicoRows("datasets");
      return { columns: _ttKeyUnion(rows), rows: rows };
    },
  },
  {
    id: "insilico_lod",
    label: "In-Silico LoD Series",
    tab: "In-Silico",
    desc: "Per-organism dilution series with limit-of-detection counts.",
    joinKey: TT_JOIN.NONE,
    filterable: false,
    defaultOn: true,
    available: () => _ttHasRows(() => _ttInsilicoRows("organisms")),
    build() {
      const rows = _ttInsilicoRows("organisms");
      return { columns: _ttKeyUnion(rows), rows: rows };
    },
  },
];

function _ttDatasetById(id) {
  return TT_EXPORT_DATASETS.find((d) => d.id === id) || null;
}

/* ── Per-dataset column selection ────────────────────────────────────────────
   A detections export is 53 columns wide; most people want a handful. Each
   dataset in the picker carries its own column multi-select — the same
   show-all / hide-all pattern as the tables' "Table Display" dialog — and the
   choice applies to every output shape, since they all build through
   _ttBuildDataset().

   Only `columns` is narrowed, never the row objects: the wide join keys on
   Specimen ID / Taxonomic ID # by name, and it must keep working even when the
   user has dropped those columns from what gets printed.                      */

const _TT_COL_SEL = {}; // dataset id -> [kept column names]  (absent = all)
const _TT_COL_CACHE = {}; // dataset id -> full column list, from the last build

/** Full (unfiltered) column list for a dataset, cached from the picker build. */
function _ttAllColumns(dsId, ctx) {
  if (_TT_COL_CACHE[dsId]) return _TT_COL_CACHE[dsId];
  const ds = _ttDatasetById(dsId);
  if (!ds) return [];
  const built = _ttBuildDataset(ds, ctx || _ttExportContext(true), true);
  _TT_COL_CACHE[dsId] = built.columns;
  return built.columns;
}

/** Columns actually exported for a dataset, in the dataset's own order. */
function _ttSelectedColumns(dsId, allCols) {
  const sel = _TT_COL_SEL[dsId];
  if (!sel || !sel.length) return allCols;
  const keep = new Set(sel);
  const out = allCols.filter((c) => keep.has(c));
  // A stale selection (the run changed shape under it) must not export nothing.
  return out.length ? out : allCols;
}

function _ttColSelActive(dsId, allCols) {
  const sel = _TT_COL_SEL[dsId];
  if (!sel) return false;
  return _ttSelectedColumns(dsId, allCols).length < allCols.length;
}

/** Build one dataset, guarding against a tab-specific helper throwing.
 *  `raw` skips the user's column selection — used to populate the picker. */
function _ttBuildDataset(ds, ctx, raw) {
  let out;
  try {
    out = ds.build(ctx) || { columns: [], rows: [] };
  } catch (e) {
    console.warn("[export-data] dataset failed:", ds.id, e);
    return { columns: [], rows: [], error: String((e && e.message) || e) };
  }
  const columns = out.columns || _ttKeyUnion(out.rows);
  _TT_COL_CACHE[ds.id] = columns;
  return {
    columns: raw ? columns : _ttSelectedColumns(ds.id, columns),
    rows: out.rows || [],
  };
}

function _ttExportContext(filtered) {
  const rowsAll = typeof DATA !== "undefined" ? DATA : [];
  let fd = rowsAll;
  try {
    if (filtered && typeof filteredData === "function") fd = filteredData();
  } catch (e) {
    fd = rowsAll;
  }
  return { filtered: !!filtered, rowsAll: rowsAll, fd: fd };
}

/* ── Wide join ───────────────────────────────────────────────────────────────
   Backbone: the detections rows (specimen x organism). Every other selected
   dataset contributes its columns to the matching row, namespaced by dataset
   label so "Mean Depth" from Coverage never collides with the detections
   column of the same name. Datasets with joinKey NONE (per-gene hits, in-silico
   series, novelty candidates — all many-per-key) cannot be widened and are
   reported back to the caller so the UI can say so. */
function _ttWideJoin(selectedIds, ctx) {
  const built = new Map();
  selectedIds.forEach((id) => {
    const ds = _ttDatasetById(id);
    if (ds) built.set(id, { ds: ds, data: _ttBuildDataset(ds, ctx) });
  });

  const base = built.get("detections");
  const backbone = base
    ? base.data
    : _ttBuildDataset(_ttDatasetById("detections"), ctx); // always needed as the spine

  const sampleKey = (r) => String(r["Specimen ID"] == null ? "" : r["Specimen ID"]);
  const orgKey = (r) => String(r["Taxonomic ID #"] || r["Detected Organism"] || "");

  const columns = backbone.columns.slice();
  const rows = backbone.rows.map((r) => Object.assign({}, r));
  const skipped = [];

  built.forEach((entry, id) => {
    if (id === "detections") return;
    const ds = entry.ds;
    const data = entry.data;
    if (!ds.joinKey || !data.rows.length) {
      if (!ds.joinKey) skipped.push(ds.label);
      return;
    }
    const index = new Map();
    data.rows.forEach((r) => {
      let k;
      if (ds.joinKey === TT_JOIN.SAMPLE) k = sampleKey(r);
      else if (ds.joinKey === TT_JOIN.ORG) k = orgKey(r);
      else k = sampleKey(r) + "\u0001" + orgKey(r);
      if (!index.has(k)) index.set(k, r);
    });
    // Namespace this dataset's columns, minus the keys already on the spine.
    const keyCols = new Set(["Specimen ID", "Taxonomic ID #", "Detected Organism"]);
    const outCols = data.columns.filter((c) => !keyCols.has(c));
    outCols.forEach((c) => columns.push(ds.label + " · " + c));
    rows.forEach((row) => {
      let k;
      if (ds.joinKey === TT_JOIN.SAMPLE) k = sampleKey(row);
      else if (ds.joinKey === TT_JOIN.ORG) k = orgKey(row);
      else k = sampleKey(row) + "\u0001" + orgKey(row);
      const hit = index.get(k);
      outCols.forEach((c) => {
        row[ds.label + " · " + c] = hit ? (hit[c] == null ? "" : hit[c]) : "";
      });
    });
  });

  return { columns: columns, rows: rows, skipped: skipped };
}

/** Stack every selected dataset into one table with a leading Dataset column
 *  and the union of all columns. The "single CSV" shape for datasets that do
 *  not share a join key. */
function _ttStack(selectedIds, ctx) {
  const columns = ["Dataset"];
  const seen = new Set(columns);
  const rows = [];
  selectedIds.forEach((id) => {
    const ds = _ttDatasetById(id);
    if (!ds) return;
    const data = _ttBuildDataset(ds, ctx);
    data.columns.forEach((c) => {
      if (!seen.has(c)) {
        seen.add(c);
        columns.push(c);
      }
    });
    data.rows.forEach((r) => rows.push(Object.assign({ Dataset: ds.label }, r)));
  });
  return { columns: columns, rows: rows };
}

/* ── Writers ─────────────────────────────────────────────────────────────── */

function _ttAoa(table) {
  const out = [table.columns.slice()];
  table.rows.forEach((r) => out.push(table.columns.map((c) => (r[c] == null ? "" : r[c]))));
  return out;
}

function _ttCsvText(table, delimiter) {
  const d = delimiter || ",";
  return _ttAoa(table)
    .map((row) => row.map((v) => _delimitedEscape(v, d)).join(d))
    .join("\r\n");
}

/** Excel sheet names: <=31 chars, no []:*?/\ — and unique within the book. */
function _ttSheetName(label, used) {
  let base = String(label || "Sheet").replace(/[\\/\?\*\[\]:]/g, "-").slice(0, 31);
  if (!base) base = "Sheet";
  let name = base;
  let n = 2;
  while (used.has(name)) {
    const suffix = "_" + n++;
    name = base.slice(0, 31 - suffix.length) + suffix;
  }
  used.add(name);
  return name;
}

function _ttExportRun(opts) {
  const ids = opts.datasets || [];
  if (!ids.length) {
    alert("Select at least one dataset to export.");
    return;
  }
  const ctx = _ttExportContext(opts.filtered);
  const stamp = new Date().toISOString().slice(0, 10);
  const base = _slug(opts.filename || "taxtriage-export") + "-" + stamp;
  const delimiter = opts.delimiter === "\\t" ? "\t" : opts.delimiter || ",";

  // ── Pivot / crosstab ──────────────────────────────────────────────────────
  // Built from the detections joined to run metadata, not from the dataset
  // checkboxes — it answers "how many hits to X per site / host / run", which is
  // a different question from "give me these tables".
  if (opts.format === "pivot") {
    const pOpts = opts.pivot || { rowDim: "organism", field: "", measure: "detections", shape: "wide" };
    const table = _ttBuildPivot(ctx, pOpts);
    if (!table.rows.length) {
      alert("Nothing to pivot with the current filters.");
      return;
    }
    const name = base + "-pivot-" + _slug(table.meta.dim.label) + (pOpts.field ? "-by-" + _slug(pOpts.field) : "");
    if (pOpts.filetype === "xlsx") {
      if (typeof XLSX === "undefined") {
        alert("XLSX export needs the SheetJS library, which this report could not load. Use CSV instead.");
        return;
      }
      const wb = XLSX.utils.book_new();
      const used = new Set();
      XLSX.utils.book_append_sheet(
        wb,
        XLSX.utils.aoa_to_sheet(_ttManifestAoa([], opts, ctx, table)),
        _ttSheetName("Export Info", used),
      );
      XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(_ttAoa(table)), _ttSheetName("Pivot", used));
      XLSX.writeFile(wb, name + ".xlsx");
    } else {
      _downloadText(
        _ttCsvText(table, delimiter),
        name + "." + (delimiter === "\t" ? "tsv" : "csv"),
        "text/plain;charset=utf-8",
      );
    }
    _ttExportToast(
      table.rows.length.toLocaleString() +
        " rows × " +
        table.columns.length +
        " columns — " +
        table.meta.measure.label +
        (pOpts.field ? " by " + _ttPrettyField(pOpts.field) : "") +
        " (" +
        pOpts.shape +
        " format)",
    );
    return;
  }

  if (opts.format === "xlsx-multi") {
    if (typeof XLSX === "undefined") {
      alert("XLSX export needs the SheetJS library, which this report could not load. Use one of the CSV formats.");
      return;
    }
    const wb = XLSX.utils.book_new();
    const used = new Set();
    let wrote = 0;
    // Sheet 1 is a manifest: what was exported, under which filters, when.
    XLSX.utils.book_append_sheet(
      wb,
      XLSX.utils.aoa_to_sheet(_ttManifestAoa(ids, opts, ctx)),
      _ttSheetName("Export Info", used),
    );
    ids.forEach((id) => {
      const ds = _ttDatasetById(id);
      if (!ds) return;
      const data = _ttBuildDataset(ds, ctx);
      XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(_ttAoa(data)), _ttSheetName(ds.label, used));
      wrote += data.rows.length;
    });
    XLSX.writeFile(wb, base + ".xlsx");
    _ttExportToast(ids.length + " sheet(s), " + wrote.toLocaleString() + " rows → " + base + ".xlsx");
    return;
  }

  if (opts.format === "xlsx-wide" || opts.format === "csv-wide") {
    const table = _ttWideJoin(ids, ctx);
    if (table.skipped.length) {
      _ttExportToast("Not joinable, left out of the wide sheet: " + table.skipped.join(", "), true);
    }
    if (opts.format === "csv-wide") {
      _downloadText(_ttCsvText(table, delimiter), base + "-wide." + (delimiter === "\t" ? "tsv" : "csv"), "text/plain;charset=utf-8");
    } else {
      if (typeof XLSX === "undefined") {
        alert("XLSX export needs the SheetJS library, which this report could not load. Use one of the CSV formats.");
        return;
      }
      const wb = XLSX.utils.book_new();
      const used = new Set();
      XLSX.utils.book_append_sheet(
        wb,
        XLSX.utils.aoa_to_sheet(_ttManifestAoa(ids, opts, ctx)),
        _ttSheetName("Export Info", used),
      );
      XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(_ttAoa(table)), _ttSheetName("Combined", used));
      XLSX.writeFile(wb, base + "-wide.xlsx");
    }
    _ttExportToast(table.rows.length.toLocaleString() + " rows × " + table.columns.length + " columns exported");
    return;
  }

  // csv-stacked
  const table = _ttStack(ids, ctx);
  _downloadText(
    _ttCsvText(table, delimiter),
    base + "." + (delimiter === "\t" ? "tsv" : "csv"),
    "text/plain;charset=utf-8",
  );
  _ttExportToast(table.rows.length.toLocaleString() + " rows from " + ids.length + " dataset(s) exported");
}

/** Provenance sheet — what the numbers in this workbook actually are. */
function _ttManifestAoa(ids, opts, ctx, pivotTable) {
  const rows = [
    ["TaxTriage combined data export"],
    ["Generated", new Date().toISOString()],
    ["Report built", (typeof BOOT !== "undefined" && BOOT.report_generated_at) || ""],
    ["Pipeline revision", (typeof BOOT !== "undefined" && BOOT.pipeline_revision) || ""],
    ["Pipeline commit", (typeof BOOT !== "undefined" && BOOT.pipeline_commit) || ""],
    ["Shape", opts.format],
    ["Rows", opts.filtered ? "Current sidebar filters applied" : "All rows (filters ignored)"],
    [],
  ];
  if (opts.filtered) {
    const val = (id) => {
      const el = document.getElementById(id);
      return el ? el.value : "";
    };
    const chk = (id) => {
      const el = document.getElementById(id);
      return el ? (el.checked ? "Yes" : "No") : "";
    };
    rows.push(["Active filters"]);
    rows.push(["Search", val("filter-text")]);
    rows.push(["Search scope", val("filter-scope")]);
    rows.push(["Min TASS score", val("filter-min")]);
    rows.push(["High consequence only", chk("filter-hc")]);
    rows.push(["Passing only", chk("filter-pass")]);
    rows.push(["Taxonomic level", val("view-level")]);
    const hidden = Object.keys((typeof sampleHidden !== "undefined" && sampleHidden) || {}).filter(
      (k) => sampleHidden[k],
    );
    rows.push(["Hidden samples", hidden.length ? hidden.join("; ") : "none"]);
    rows.push([
      "Specimen merge",
      typeof specimenMergeEnabled !== "undefined" && specimenMergeEnabled ? "On (" + specimenTassAgg + ")" : "Off",
    ]);
    rows.push([]);
  }
  if (pivotTable) {
    const p = opts.pivot || {};
    rows.push(["Pivot"]);
    rows.push(["Rows", pivotTable.meta.dim.label]);
    rows.push(["Columns", p.field ? _ttPrettyField(p.field) + " (run metadata)" : "(none — totals only)"]);
    rows.push(["Measure", pivotTable.meta.measure.label]);
    rows.push(["Layout", p.shape === "long" ? "Long (one row per pair)" : "Wide (one column per value)"]);
    rows.push(["Size", pivotTable.rows.length + " rows x " + pivotTable.columns.length + " columns"]);
    rows.push([]);
    rows.push(["Source", "Detections joined to each sample's run metadata"]);
    return rows;
  }
  rows.push(["Datasets included", "Tab", "Rows", "Columns"]);
  ids.forEach((id) => {
    const ds = _ttDatasetById(id);
    if (!ds) return;
    const built = _ttBuildDataset(ds, ctx);
    const all = _TT_COL_CACHE[id] || built.columns;
    rows.push([
      ds.label,
      ds.tab,
      built.rows.length,
      built.columns.length < all.length ? built.columns.length + " of " + all.length + " (narrowed)" : built.columns.length,
    ]);
  });
  const narrowed = ids.filter((id) => _ttColSelActive(id, _TT_COL_CACHE[id] || []));
  if (narrowed.length) {
    rows.push([]);
    rows.push(["Columns kept"]);
    narrowed.forEach((id) => {
      const ds = _ttDatasetById(id);
      rows.push([ds ? ds.label : id, (_TT_COL_SEL[id] || []).join(", ")]);
    });
  }
  return rows;
}

let _ttToastTimer = null;
function _ttExportToast(msg, warn) {
  let el = document.getElementById("data-export-toast");
  if (!el) {
    el = document.createElement("div");
    el.id = "data-export-toast";
    document.body.appendChild(el);
  }
  el.textContent = msg;
  el.classList.toggle("warn", !!warn);
  el.classList.add("show");
  clearTimeout(_ttToastTimer);
  _ttToastTimer = setTimeout(() => el.classList.remove("show"), warn ? 7000 : 4000);
}

/* ── The Export popup ────────────────────────────────────────────────────────
   ONE dialog behind the sidebar's single Export button. The type selector at the
   top switches between the two things worth taking out of a report:

     Data tables  the dataset picker + output shape (this module)
     Report PDF   the printable layout (_exportReportPdf in 02_utilities.js)

   The PDF pane owns #report-pdf-legend-place, which _exportReportPdf() reads —
   that element used to live in a separate modal of its own.                   */

function _ttEnsureExportModal() {
  let overlay = document.getElementById("data-export-overlay");
  if (overlay) return overlay;
  overlay = document.createElement("div");
  overlay.id = "data-export-overlay";
  overlay.className = "export-modal-overlay";
  overlay.innerHTML = `
    <div class="export-modal data-export-modal" role="dialog" aria-modal="true" aria-labelledby="data-export-title">
      <header>
        <i class="fas fa-file-export"></i>
        <span id="data-export-title">Export</span>
        <button type="button" id="data-export-close" title="Close" aria-label="Close">&times;</button>
      </header>

      <div class="tt-export-types" role="tablist" aria-label="What to export">
        <button type="button" class="tt-export-type active" data-extype="data" role="tab" aria-selected="true">
          <i class="fas fa-table"></i>
          <span class="tt-export-type-main">Data tables</span>
          <span class="tt-export-type-sub">Excel / CSV from any tabs</span>
        </button>
        <button type="button" class="tt-export-type" data-extype="pdf" role="tab" aria-selected="false">
          <i class="fas fa-file-pdf"></i>
          <span class="tt-export-type-main">Report PDF</span>
          <span class="tt-export-type-sub">Printable layout of this view</span>
        </button>
      </div>

      <div class="tt-export-pane" id="tt-export-pane-data">
        <div class="data-export-intro">
          Pick the tables you want — from as many tabs as you like — and get them back as one file.
        </div>
        <div class="data-export-toolbar">
          <button type="button" id="data-export-all">Select all</button>
          <button type="button" id="data-export-none">Clear</button>
          <label id="data-export-filtered-lbl" title="Export only what the current sidebar filters, TASS cutoff and sample visibility leave in view">
            <input type="checkbox" id="data-export-filtered" checked /> Apply current filters
          </label>
        </div>
        <div class="data-export-list" id="data-export-list"></div>
        <div class="data-export-opts">
          <label>Output
            <select id="data-export-format">
              <option value="xlsx-multi">Excel workbook — one sheet per table</option>
              <option value="xlsx-wide">Excel — single joined sheet (Sample × Organism)</option>
              <option value="csv-wide">CSV — single joined table (Sample × Organism)</option>
              <option value="csv-stacked">CSV — all tables stacked, with a Dataset column</option>
              <option value="pivot">Pivot — organisms × a metadata field (counts)</option>
            </select>
          </label>
          <label id="data-export-delim-lbl">Delimiter
            <select id="data-export-delim">
              <option value=",">Comma (,)</option>
              <option value="\\t">Tab</option>
              <option value=";">Semicolon (;)</option>
              <option value="|">Pipe (|)</option>
            </select>
          </label>
          <label>File name
            <input id="data-export-name" value="taxtriage-export" maxlength="60" />
          </label>
        </div>
        <div class="tt-pivot-opts" id="tt-pivot-opts" hidden>
          <div class="tt-pivot-hd"><i class="fas fa-table-cells"></i> Pivot — count detections against a metadata field</div>
          <div class="tt-pivot-grid">
            <label>Rows
              <select id="tt-pivot-row"></select>
            </label>
            <label>Columns (metadata)
              <select id="tt-pivot-field"></select>
            </label>
            <label>Count
              <select id="tt-pivot-measure"></select>
            </label>
            <label>Layout
              <select id="tt-pivot-shape">
                <option value="wide">Wide — one column per value</option>
                <option value="long">Long — one row per pair (tidy)</option>
              </select>
            </label>
            <label>File
              <select id="tt-pivot-filetype">
                <option value="csv">CSV</option>
                <option value="xlsx">Excel</option>
              </select>
            </label>
          </div>
          <div class="tt-pivot-preview" id="tt-pivot-preview"></div>
        </div>
        <div class="data-export-note" id="data-export-note"></div>
      </div>

      <div class="tt-export-pane" id="tt-export-pane-pdf" hidden>
        <div class="tt-pdf-body">
          <p>
            This prepares every report tab using the current filters, sample visibility, chart controls and table state.
          </p>
          <p class="tt-pdf-sub">
            Your browser's print dialog opens next — set the destination to <b>Save as PDF</b> to write the file.
            Landscape orientation, with <i>Headers and footers</i> off, gives the best result.
          </p>
          <label class="tt-pdf-legend">
            Sample-colour &amp; TASS-cutoff legend
            <select id="report-pdf-legend-place">
              <option value="cover" selected>On cover page</option>
              <option value="footer">Repeat on every page</option>
              <option value="off">Don't include</option>
            </select>
          </label>
          <p class="tt-pdf-foot">
            The legend maps each plot colour to its sample and records the applied TASS cutoffs.
          </p>
        </div>
      </div>

      <div class="export-modal-actions">
        <button type="button" id="data-export-cancel">Cancel</button>
        <button type="button" class="primary" id="data-export-save"><i class="fas fa-download"></i> Export</button>
      </div>
    </div>`;
  document.body.appendChild(overlay);

  const close = () => (overlay.style.display = "none");
  overlay.querySelector("#data-export-close").addEventListener("click", close);
  overlay.querySelector("#data-export-cancel").addEventListener("click", close);
  overlay.addEventListener("click", (e) => {
    if (e.target === overlay) close();
  });
  overlay.querySelectorAll(".tt-export-type").forEach((btn) => {
    btn.addEventListener("click", () => _ttSetExportType(btn.dataset.extype));
  });
  overlay.querySelector("#data-export-all").addEventListener("click", () => {
    overlay.querySelectorAll(".data-export-cb:not(:disabled)").forEach((cb) => (cb.checked = true));
    _ttSyncExportModal();
  });
  overlay.querySelector("#data-export-none").addEventListener("click", () => {
    overlay.querySelectorAll(".data-export-cb").forEach((cb) => (cb.checked = false));
    _ttSyncExportModal();
  });
  overlay.querySelector("#data-export-filtered").addEventListener("change", () => _ttFillExportDataList());
  overlay.querySelector("#data-export-format").addEventListener("change", () => _ttSyncExportModal());
  // #tt-pivot-filetype is in here too: it decides whether the delimiter control
  // is relevant, so the dialog has to re-sync when it changes.
  ["#tt-pivot-row", "#tt-pivot-field", "#tt-pivot-measure", "#tt-pivot-shape", "#tt-pivot-filetype"].forEach((sel) => {
    overlay.querySelector(sel).addEventListener("change", () => _ttSyncExportModal());
  });
  overlay.addEventListener("change", (e) => {
    if (e.target.classList && e.target.classList.contains("data-export-cb")) _ttSyncExportModal();
  });
  // The columns chip sits inside the item's <label>, so a plain click would
  // also toggle that table's checkbox. Swallow it and open the picker instead.
  overlay.addEventListener("click", (e) => {
    const chip = e.target.closest && e.target.closest(".data-export-cols");
    if (!chip) return;
    e.preventDefault();
    e.stopPropagation();
    _ttOpenColumnPicker(chip.dataset.colsFor);
  });
  overlay.querySelector("#data-export-save").addEventListener("click", () => {
    if (_TT_EXPORT_TYPE === "pdf") {
      close();
      if (typeof _exportReportPdf === "function") {
        _exportReportPdf().catch((err) => {
          if (typeof _hidePdfProgress === "function") _hidePdfProgress();
          alert("PDF export failed: " + (err && err.message ? err.message : err));
        });
      }
      return;
    }
    const fmt = overlay.querySelector("#data-export-format").value;
    const opts = {
      datasets: Array.from(overlay.querySelectorAll(".data-export-cb:checked")).map((cb) => cb.value),
      filtered: overlay.querySelector("#data-export-filtered").checked,
      format: fmt,
      delimiter: overlay.querySelector("#data-export-delim").value,
      filename: overlay.querySelector("#data-export-name").value,
      pivot: _ttPivotOptsFromModal(),
    };
    if (fmt !== "pivot" && !opts.datasets.length) {
      alert("Select at least one dataset to export.");
      return;
    }
    close();
    // Let the modal paint away before a large build blocks the thread.
    setTimeout(() => _ttExportRun(opts), 30);
  });
  _ttFillPivotControls();
  return overlay;
}

let _TT_EXPORT_TYPE = "data";

function _ttSetExportType(type) {
  const overlay = document.getElementById("data-export-overlay");
  if (!overlay) return;
  _TT_EXPORT_TYPE = type === "pdf" ? "pdf" : "data";
  overlay.querySelectorAll(".tt-export-type").forEach((b) => {
    const on = b.dataset.extype === _TT_EXPORT_TYPE;
    b.classList.toggle("active", on);
    b.setAttribute("aria-selected", on ? "true" : "false");
  });
  overlay.querySelector("#tt-export-pane-data").hidden = _TT_EXPORT_TYPE !== "data";
  overlay.querySelector("#tt-export-pane-pdf").hidden = _TT_EXPORT_TYPE !== "pdf";
  overlay.querySelector("#data-export-title").textContent =
    _TT_EXPORT_TYPE === "pdf" ? "Export — Report PDF" : "Export — Data";
  const save = overlay.querySelector("#data-export-save");
  save.innerHTML =
    _TT_EXPORT_TYPE === "pdf"
      ? '<i class="fas fa-print"></i> Prepare PDF'
      : '<i class="fas fa-download"></i> Export';
  overlay.querySelector(".data-export-modal").classList.toggle("pdf-mode", _TT_EXPORT_TYPE === "pdf");
  if (_TT_EXPORT_TYPE === "data") _ttSyncExportModal();
  else save.disabled = false;
}

/** Populate the pivot dropdowns. Metadata fields come from whatever the run
 *  actually carries (including columns uploaded on the Metadata tab), so the
 *  list is rebuilt each time the dialog opens. */
function _ttFillPivotControls() {
  const overlay = document.getElementById("data-export-overlay");
  if (!overlay) return;
  const rowSel = overlay.querySelector("#tt-pivot-row");
  const fieldSel = overlay.querySelector("#tt-pivot-field");
  const measSel = overlay.querySelector("#tt-pivot-measure");
  const prevRow = rowSel.value;
  const prevField = fieldSel.value;
  const prevMeas = measSel.value;

  rowSel.innerHTML = TT_PIVOT_ROWS.map((d) => '<option value="' + d.id + '">' + _ttEsc(d.label) + "</option>").join("");
  measSel.innerHTML = TT_PIVOT_MEASURES.map(
    (m) => '<option value="' + m.id + '">' + _ttEsc(m.label) + "</option>",
  ).join("");

  const fields = _ttMetaFields(true);
  fieldSel.innerHTML =
    '<option value="">(none — one row per group, totals only)</option>' +
    fields
      .map(
        (f) =>
          '<option value="' +
          _ttEsc(f.field) +
          '">' +
          _ttEsc(f.label) +
          " — " +
          f.n +
          (f.n === 1 ? " value" : " values") +
          "</option>",
      )
      .join("");

  rowSel.value = prevRow || "organism";
  measSel.value = prevMeas || "detections";
  if (prevField && fields.some((f) => f.field === prevField)) fieldSel.value = prevField;
  else {
    // Default to the first field that actually varies — a single-valued column
    // axis produces a one-column crosstab, which helps nobody.
    const first = fields.find((f) => f.n > 1) || fields[0];
    fieldSel.value = first ? first.field : "";
  }
}

function _ttPivotOptsFromModal() {
  const overlay = document.getElementById("data-export-overlay");
  if (!overlay) return null;
  return {
    rowDim: overlay.querySelector("#tt-pivot-row").value,
    field: overlay.querySelector("#tt-pivot-field").value,
    measure: overlay.querySelector("#tt-pivot-measure").value,
    shape: overlay.querySelector("#tt-pivot-shape").value,
    filetype: overlay.querySelector("#tt-pivot-filetype").value,
  };
}

/** (Re)build the checkbox list, with live row counts for the current filter
 *  state so the user can see what they are about to download. */
function _ttFillExportDataList() {
  const overlay = document.getElementById("data-export-overlay");
  if (!overlay) return;
  const list = overlay.querySelector("#data-export-list");
  const filtered = overlay.querySelector("#data-export-filtered").checked;
  const prev = new Map();
  overlay.querySelectorAll(".data-export-cb").forEach((cb) => prev.set(cb.value, cb.checked));
  const ctx = _ttExportContext(filtered);

  const byTab = new Map();
  TT_EXPORT_DATASETS.forEach((ds) => {
    let ok = false;
    try {
      ok = ds.available();
    } catch (e) {
      ok = false;
    }
    if (!ok) return;
    if (!byTab.has(ds.tab)) byTab.set(ds.tab, []);
    byTab.get(ds.tab).push(ds);
  });

  let html = "";
  byTab.forEach((items, tab) => {
    html += '<div class="data-export-group"><div class="data-export-group-hd">' + _ttEsc(tab) + "</div>";
    items.forEach((ds) => {
      const built = _ttBuildDataset(ds, ctx);
      const n = built.rows.length;
      const checked = prev.has(ds.id) ? prev.get(ds.id) : ds.defaultOn;
      const allCols = _TT_COL_CACHE[ds.id] || built.columns;
      const kept = _ttSelectedColumns(ds.id, allCols);
      const narrowed = kept.length < allCols.length;
      html +=
        '<label class="data-export-item' +
        (n ? "" : " empty") +
        '"><input type="checkbox" class="data-export-cb" value="' +
        ds.id +
        '"' +
        (checked && n ? " checked" : "") +
        (n ? "" : " disabled") +
        ' /><span class="data-export-item-main"><span class="data-export-item-label">' +
        _ttEsc(ds.label) +
        '</span><span class="data-export-item-desc">' +
        _ttEsc(ds.desc) +
        "</span></span>" +
        '<span class="data-export-meta">' +
        '<span class="data-export-count">' +
        (n ? n.toLocaleString() + " rows" : "no data") +
        "</span>" +
        (n
          ? '<button type="button" class="data-export-cols' +
            (narrowed ? " narrowed" : "") +
            '" data-cols-for="' +
            ds.id +
            '" title="Choose which columns to export from this table">' +
            '<i class="fas fa-table-columns"></i> ' +
            (narrowed ? kept.length + " of " + allCols.length : "all " + allCols.length) +
            " cols</button>"
          : "") +
        "</span></label>";
    });
    html += "</div>";
  });
  list.innerHTML = html || '<div class="data-export-empty">This report carries no exportable datasets.</div>';
  _ttSyncExportModal();
}

function _ttEsc(s) {
  return String(s == null ? "" : s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}

function _ttSyncExportModal() {
  const overlay = document.getElementById("data-export-overlay");
  if (!overlay || _TT_EXPORT_TYPE === "pdf") return;
  const fmt = overlay.querySelector("#data-export-format").value;
  const isPivot = fmt === "pivot";
  const isCsv = fmt.indexOf("csv") === 0 || (isPivot && overlay.querySelector("#tt-pivot-filetype").value === "csv");
  overlay.querySelector("#data-export-delim-lbl").style.display = isCsv ? "" : "none";
  overlay.querySelector("#tt-pivot-opts").hidden = !isPivot;
  // The pivot is built from the detections + metadata, not from the checkboxes,
  // so hide the table picker entirely in that mode rather than leaving a dimmed
  // list the user might think they still have to set.
  overlay.querySelector(".data-export-modal").classList.toggle("pivot-mode", isPivot);

  const note = overlay.querySelector("#data-export-note");
  const save = overlay.querySelector("#data-export-save");

  if (isPivot) {
    save.disabled = false;
    const opts = _ttPivotOptsFromModal();
    const ctx = _ttExportContext(overlay.querySelector("#data-export-filtered").checked);
    let table;
    try {
      table = _ttBuildPivot(ctx, opts);
    } catch (e) {
      table = null;
    }
    const prev = overlay.querySelector("#tt-pivot-preview");
    if (!table || !table.rows.length) {
      prev.innerHTML = '<span class="warn">Nothing to pivot with the current filters.</span>';
      save.disabled = true;
    } else {
      prev.innerHTML = _ttPivotPreviewHtml(table, opts);
    }
    note.textContent =
      "The pivot is built from the detections joined to each sample's run metadata — it ignores the table selection above.";
    note.className = "data-export-note";
    return;
  }

  const ids = Array.from(overlay.querySelectorAll(".data-export-cb:checked")).map((cb) => cb.value);
  save.disabled = !ids.length;
  if (!ids.length) {
    note.textContent = "Nothing selected.";
    note.className = "data-export-note warn";
    return;
  }
  if (fmt === "xlsx-wide" || fmt === "csv-wide") {
    const unjoinable = ids.map(_ttDatasetById).filter((d) => d && !d.joinKey);
    if (unjoinable.length) {
      note.innerHTML =
        '<i class="fas fa-triangle-exclamation"></i> ' +
        _ttEsc(unjoinable.map((d) => d.label).join(", ")) +
        " have several rows per organism, so they cannot be folded into one joined row. " +
        "They will be left out of this shape — use the workbook or stacked CSV to include them.";
      note.className = "data-export-note warn";
      return;
    }
    note.textContent = ids.length + " table(s) joined on Specimen ID × Organism into a single sheet.";
  } else if (fmt === "xlsx-multi") {
    note.textContent = ids.length + " table(s), one sheet each, plus an Export Info provenance sheet.";
  } else {
    note.textContent = ids.length + " table(s) stacked into one file with a leading Dataset column.";
  }
  note.className = "data-export-note";
}

/** A few cells of the pivot, so the user can see the shape before downloading. */
function _ttPivotPreviewHtml(table, opts) {
  const maxCols = 6;
  const maxRows = 4;
  const cols = table.columns.slice(0, maxCols);
  const more = table.columns.length - cols.length;
  let html =
    '<div class="tt-pivot-preview-hd">' +
    table.rows.length.toLocaleString() +
    " rows × " +
    table.columns.length +
    " columns · " +
    _ttEsc(table.meta.measure.label) +
    (opts.field ? " by " + _ttEsc(_ttPrettyField(opts.field)) : "") +
    "</div><table><thead><tr>";
  cols.forEach((c) => (html += "<th>" + _ttEsc(c) + "</th>"));
  if (more > 0) html += '<th class="tt-pivot-more">+' + more + "</th>";
  html += "</tr></thead><tbody>";
  table.rows.slice(0, maxRows).forEach((r) => {
    html += "<tr>";
    cols.forEach((c) => (html += "<td>" + _ttEsc(r[c] === "" || r[c] == null ? "–" : r[c]) + "</td>"));
    if (more > 0) html += '<td class="tt-pivot-more">…</td>';
    html += "</tr>";
  });
  html += "</tbody></table>";
  if (table.rows.length > maxRows)
    html += '<div class="tt-pivot-more-rows">…and ' + (table.rows.length - maxRows).toLocaleString() + " more rows</div>";
  return html;
}

/* ── Column picker ───────────────────────────────────────────────────────────
   One multi-select per dataset, layered over the Export popup. Mirrors the
   tables' "Table Display" dialog (multi-select + show-all / hide-all), with a
   filter box because the detections list runs to 53 columns.                  */

function _ttEnsureColumnModal() {
  let overlay = document.getElementById("tt-cols-overlay");
  if (overlay) return overlay;
  overlay = document.createElement("div");
  overlay.id = "tt-cols-overlay";
  overlay.className = "export-modal-overlay tt-cols-overlay";
  overlay.innerHTML = `
    <div class="export-modal tt-cols-modal" role="dialog" aria-modal="true" aria-labelledby="tt-cols-title">
      <header>
        <i class="fas fa-table-columns"></i>
        <span id="tt-cols-title">Columns</span>
        <button type="button" id="tt-cols-close" title="Close" aria-label="Close">&times;</button>
      </header>
      <div class="tt-cols-body">
        <div class="tt-cols-intro">
          Ctrl/\u2318-click or drag to select. Only the highlighted columns are written for this table;
          the rest are left out of every output shape.
        </div>
        <input id="tt-cols-filter" type="text" placeholder="Filter columns\u2026" autocomplete="off" spellcheck="false" />
        <select id="tt-cols-list" multiple size="14"></select>
        <div class="tt-cols-tools">
          <button type="button" id="tt-cols-all">Select all</button>
          <button type="button" id="tt-cols-none">Select none</button>
          <button type="button" id="tt-cols-visible">Select filtered</button>
          <span id="tt-cols-count"></span>
        </div>
      </div>
      <div class="export-modal-actions">
        <button type="button" id="tt-cols-cancel">Cancel</button>
        <button type="button" class="primary" id="tt-cols-apply">Apply</button>
      </div>
    </div>`;
  document.body.appendChild(overlay);

  const close = () => (overlay.style.display = "none");
  overlay.querySelector("#tt-cols-close").addEventListener("click", close);
  overlay.querySelector("#tt-cols-cancel").addEventListener("click", close);
  overlay.addEventListener("click", (e) => {
    if (e.target === overlay) close();
  });
  const list = overlay.querySelector("#tt-cols-list");
  const setAll = (on, visibleOnly) => {
    Array.from(list.options).forEach((o) => {
      if (visibleOnly && o.hidden) return;
      o.selected = on;
    });
    _ttColCount();
  };
  overlay.querySelector("#tt-cols-all").addEventListener("click", () => setAll(true));
  overlay.querySelector("#tt-cols-none").addEventListener("click", () => setAll(false));
  overlay.querySelector("#tt-cols-visible").addEventListener("click", () => setAll(true, true));
  list.addEventListener("change", _ttColCount);
  overlay.querySelector("#tt-cols-filter").addEventListener("input", (e) => {
    const q = e.target.value.trim().toLowerCase();
    Array.from(list.options).forEach((o) => {
      o.hidden = !!q && o.value.toLowerCase().indexOf(q) === -1;
    });
  });
  overlay.querySelector("#tt-cols-apply").addEventListener("click", () => {
    const dsId = overlay._ttColsFor;
    const all = Array.from(list.options).map((o) => o.value);
    const picked = Array.from(list.selectedOptions).map((o) => o.value);
    if (!picked.length) {
      alert("Keep at least one column — deselecting every column would export an empty table.");
      return;
    }
    // Remember "everything" as no selection at all, so a dataset that later
    // gains a column still exports it.
    if (picked.length === all.length) delete _TT_COL_SEL[dsId];
    else _TT_COL_SEL[dsId] = picked;
    close();
    _ttFillExportDataList();
  });
  return overlay;
}

function _ttColCount() {
  const overlay = document.getElementById("tt-cols-overlay");
  if (!overlay) return;
  const list = overlay.querySelector("#tt-cols-list");
  const n = list.selectedOptions.length;
  const total = list.options.length;
  const el = overlay.querySelector("#tt-cols-count");
  el.textContent = n + " of " + total + " selected";
  el.classList.toggle("warn", n === 0);
}

function _ttOpenColumnPicker(dsId) {
  const ds = _ttDatasetById(dsId);
  if (!ds) return;
  const overlay = _ttEnsureColumnModal();
  overlay._ttColsFor = dsId;
  const filtered = (document.getElementById("data-export-filtered") || {}).checked;
  const all = _ttAllColumns(dsId, _ttExportContext(!!filtered));
  const kept = new Set(_ttSelectedColumns(dsId, all));
  const list = overlay.querySelector("#tt-cols-list");
  list.innerHTML = "";
  all.forEach((c) => {
    const opt = document.createElement("option");
    opt.value = c;
    opt.textContent = c;
    opt.selected = kept.has(c);
    list.appendChild(opt);
  });
  overlay.querySelector("#tt-cols-filter").value = "";
  overlay.querySelector("#tt-cols-title").textContent = "Columns: " + ds.label;
  _ttColCount();
  overlay.style.display = "flex";
}

/** Open the single Export popup. `type` is "data" (default) or "pdf". */
function openExportModal(type) {
  const overlay = _ttEnsureExportModal();
  _ttFillPivotControls();
  _ttFillExportDataList();
  _ttSetExportType(type || "data");
  overlay.style.display = "flex";
}

/** Back-compat entry points: the sidebar tour and 02_utilities.js still speak
 *  in terms of the two old dialogs. */
function openExportDataModal() {
  openExportModal("data");
}
function _openReportPdfModal() {
  openExportModal("pdf");
}

function _ttInitExportData() {
  const btn = document.getElementById("export-btn");
  // 28_meta_csv.js binds the click (it owns the sidebar's tooltip wiring); bind
  // here too only if that never ran, so the button works either way.
  if (btn && !btn.dataset.ttExportBound) {
    btn.dataset.ttExportBound = "1";
    btn.addEventListener("click", () => openExportModal("data"));
  }
  // Ctrl/Cmd+Shift+E from anywhere in the report.
  document.addEventListener("keydown", (e) => {
    if ((e.ctrlKey || e.metaKey) && e.shiftKey && (e.key === "E" || e.key === "e")) {
      e.preventDefault();
      openExportModal("data");
    }
  });
}

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", _ttInitExportData);
} else {
  _ttInitExportData();
}
