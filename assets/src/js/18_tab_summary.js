/* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB: SUMMARY            (data-tab="summary"  —  the default tab)
       -     drawSummary() — overall run overview. Mirrors the high-level
       -     features of the ODR PDF: run-level KPIs, the highest-TASS / high-
       -     consequence detections, a cross-sample genera comparison bar chart,
       -     and (only when protein/VF-AMR annotation is enabled) an annotation
       -     summary. All derived client-side from the same filteredData() view
       -     every other tab uses, so it respects the sidebar filters.
═══════════════════════════════════════════════════════════════════════════ */
function _fmtInt(n) {
  n = Math.round(num(n));
  return n.toLocaleString();
}
// Scientific-notation for large counts (>= 1e5); returns {short, full}.
function _fmtBig(n) {
  n = num(n);
  const full = Math.round(n).toLocaleString();
  if (n >= 1e5) {
    // e.g. 1.23e7 — superscript-free, compact
    const exp = Math.floor(Math.log10(n));
    const mant = n / Math.pow(10, exp);
    return { short: `${mant.toFixed(2)}×10${_sup(exp)}`, full };
  }
  return { short: full, full };
}
function _sup(n) {
  const map = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
  return String(n)
    .split("")
    .map((c) => map[c] || c)
    .join("");
}
function _orgBadges(r) {
  // ● high consequence, Ⓓ/Ⓡ molecule type
  let s = "";
  if (isTruthy(r["High Consequence"])) s += '<span title="High consequence" style="color:#cc0000">●</span> ';
  const mt = String(r["Mol Type"] || "").toLowerCase();
  if (mt === "dna") s += '<span title="DNA" style="color:#1565c0;font-weight:700">Ⓓ</span> ';
  else if (mt === "rna") s += '<span title="RNA" style="color:#6a1b9a;font-weight:700">Ⓡ</span> ';
  return s;
}

// External link to the NCBI Taxonomy page for a detection row's taxon id.
function _ncbiLink(r) {
  const tid = String(r["Taxonomic ID #"] || "").trim();
  if (!tid || !/^\d+$/.test(tid)) return "";
  const url = `https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=${tid}`;
  return (
    ` <a href="${url}" target="_blank" rel="noopener" class="ncbi-link" title="View taxon ${tid} on NCBI"` +
    ` onclick="event.stopPropagation()"><i class="fas fa-up-right-from-square"></i></a>`
  );
}

// Keyword chips for a detection row, drawn from category/status, flags,
// taxonomy, and BSL level.
function _rowKeywords(r) {
  const chips = [];
  const cat = r["Microbial Category"];
  if (cat && cat !== "Unknown") chips.push({ t: cat, c: "kw-cat" });
  const status = r["Status"];
  if (status) chips.push({ t: status, c: "" });
  const annClass = r["Ann Class"];
  if (annClass) chips.push({ t: annClass, c: "" });
  if (isTruthy(r["High Consequence"])) chips.push({ t: "high-consequence", c: "kw-hc" });
  if (isTruthy(r["IsAnnotated"])) chips.push({ t: "annotated", c: "kw-annot" });
  if (_orgHasAmr(r)) chips.push({ t: "AMR", c: "kw-amr" });
  const mt = String(r["Mol Type"] || "").toLowerCase();
  if (mt === "dna" || mt === "rna") chips.push({ t: mt.toUpperCase(), c: "" });
  const genus = r["Genus"];
  if (genus) chips.push({ t: genus, c: "" });
  const fam = r["Family"];
  if (fam && fam !== genus) chips.push({ t: fam, c: "" });
  const bsl = r["BSL Level"];
  if (bsl) chips.push({ t: bsl, c: "kw-bsl" });
  return chips.map((ch) => `<span class="kw-chip ${ch.c}">${String(ch.t)}</span>`).join("");
}

// Build a Map of full CONTIG_DATA entries keyed by sample||organism||taxon.
let _SUM_HIST_MAP = null;
function _histMap() {
  if (_SUM_HIST_MAP) return _SUM_HIST_MAP;
  _SUM_HIST_MAP = new Map();
  (CONTIG_DATA || []).forEach((cd) => {
    const k = `${cd.sample}||${cd.organism}||${cd.taxon_id}`;
    if (cd.depth_histogram || cd.breadth_histogram) _SUM_HIST_MAP.set(k, cd);
  });
  return _SUM_HIST_MAP;
}
function _invalidateSummaryHistMap() {
  _SUM_HIST_MAP = null;
  // Drop the VF/AMR index (and its attached per-sample summary) so uploads /
  // state loads rebuild it against the new protein data instead of serving
  // a stale cache.
  try {
    _VFAMR_INDEX = null;
  } catch (e) {}
}
const _DEPTH_BINS = ["0x", "1-5x", "5-10x", "10-50x", ">50x"];
const _DEPTH_COLORS = ["#dee2e6", "#a5d8ff", "#74c0fc", "#4dabf7", "#1c7ed6"];
function _histFor(r) {
  return _histMap().get(`${r["Specimen ID"]}||${r["Detected Organism"]}||${r["Taxonomic ID #"]}`);
}
// Per-position coverage profile from a CONTIG_DATA entry's breadth histogram.
// Returns an array of 0–100 % values (one per genomic bin), or null.
function _posBins(cd) {
  const bh = cd && cd.breadth_histogram;
  if (bh && Array.isArray(bh.bins) && bh.bins.length) return bh.bins.map((v) => num(v));
  return null;
}
// Sparkline: per-position coverage profile when available, else depth-bin bars.
function _sparkHtml(cd) {
  const pos = _posBins(cd);
  if (pos && pos.length) {
    // Downsample to ~24 columns for a compact inline profile.
    const N = 24;
    const step = Math.max(1, Math.ceil(pos.length / N));
    const cols = [];
    for (let i = 0; i < pos.length; i += step) {
      const seg = pos.slice(i, i + step);
      cols.push(seg.reduce((s, v) => s + v, 0) / seg.length);
    }
    return (
      '<span class="spark spark-pos">' +
      cols
        .map((v) => {
          const h = Math.max(1, Math.round((Math.min(100, v) / 100) * 16));
          const c = v >= 50 ? "#1c7ed6" : v > 0 ? "#74c0fc" : "#dee2e6";
          return `<i style="height:${h}px;background:${c}"></i>`;
        })
        .join("") +
      "</span>"
    );
  }
  const hist = (cd && cd.depth_histogram) || {};
  const vals = _DEPTH_BINS.map((b) => num(hist[b]));
  const max = Math.max(1, ...vals);
  return (
    '<span class="spark">' +
    vals
      .map(
        (v, i) => `<i style="height:${Math.max(2, Math.round((v / max) * 16))}px;background:${_DEPTH_COLORS[i]}"></i>`,
      )
      .join("") +
    "</span>"
  );
}
function _histPreviewHtml(r, cd) {
  const pos = _posBins(cd);
  const head = `<b><i>${r["Detected Organism"]}</i></b><br><span style="color:#aaa">${r["Specimen ID"]}</span><br>`;
  if (pos && pos.length) {
    const W = 180,
      H = 56;
    const n = pos.length;
    const path = pos
      .map((v, i) => {
        const x = (i / Math.max(1, n - 1)) * W;
        const y = H - (Math.min(100, v) / 100) * H;
        return `${i === 0 ? "M" : "L"}${x.toFixed(1)},${y.toFixed(1)}`;
      })
      .join(" ");
    const area = `${path} L${W},${H} L0,${H} Z`;
    const bh = cd && cd.breadth_histogram;
    const binSzBp = (bh && bh.bin_size) || 0;
    const binSzFmt =
      binSzBp >= 1e6
        ? (binSzBp / 1e6).toFixed(2) + " Mb"
        : binSzBp >= 1e3
        ? (binSzBp / 1e3).toFixed(1) + " kb"
        : binSzBp + " bp";
    return (
      head +
      `<span style="color:#aaa;font-size:0.85em">% of each bin with ≥1× coverage · Y: 0–100%</span><br>` +
      `<svg width="${W}" height="${H + 12}">` +
      `<line x1="0" y1="0" x2="0" y2="${H}" stroke="#ccc" stroke-width="0.8"></line>` +
      `<text x="2" y="7" font-size="6" fill="#999">100%</text>` +
      `<text x="2" y="${H - 1}" font-size="6" fill="#999">0%</text>` +
      `<path d="${area}" fill="#d0ebff"></path>` +
      `<path d="${path}" fill="none" stroke="#1c7ed6" stroke-width="1"></path>` +
      `<text x="0" y="${H + 10}" font-size="7" fill="#666">start</text>` +
      `<text x="${W}" y="${H + 10}" font-size="7" text-anchor="end" fill="#666">end</text>` +
      `</svg>` +
      `<div style="color:#777;font-size:0.72em;margin-top:3px;line-height:1.4">Genomic coverage based on the percent of each bin with at least 1× coverage across the reference genome.${
        binSzBp ? ` Bin size for this genome: <b>${binSzFmt}</b>.` : ""
      } Bin sizes vary by genome (~0.3 kb – ~70 kb); viral genomes typically have much smaller bins than bacterial, which can significantly affect histogram appearance.</div>` +
      `<span style="color:#888;font-size:0.78em">↗ click to open in Histogram tab</span>`
    );
  }
  const hist = (cd && cd.depth_histogram) || {};
  const vals = _DEPTH_BINS.map((b) => num(hist[b]));
  const max = Math.max(1, ...vals);
  const bw = 24,
    gap = 4,
    H = 60;
  let bars = "";
  vals.forEach((v, i) => {
    const h = Math.max(1, Math.round((v / max) * H));
    const x = i * (bw + gap);
    bars +=
      `<rect x="${x}" y="${H - h}" width="${bw}" height="${h}" fill="${_DEPTH_COLORS[i]}"></rect>` +
      `<text x="${x + bw / 2}" y="${H + 10}" font-size="7" text-anchor="middle" fill="#666">${_DEPTH_BINS[i]}</text>`;
  });
  return (
    head +
    `<span style="color:#aaa;font-size:0.85em">coverage-depth distribution (positions)</span><br>` +
    `<svg width="${(bw + gap) * 5}" height="${H + 14}">${bars}</svg>` +
    `<br><span style="color:#888;font-size:0.78em">↗ click to open in Histogram tab</span>`
  );
}

// Jump to the Histogram tab and pre-select a given organism + sample.
function _jumpToHistogram(organism, taxonId, sample) {
  const btn = document.getElementById("hist-tab-btn") || document.querySelector('.tab-btn[data-tab="histogram"]');
  if (!btn || btn.classList.contains("hidden")) return;
  btn.click();
  setTimeout(() => {
    const orgSel = document.getElementById("hist-org-sel");
    const sampleSel = document.getElementById("hist-sample-sel");
    const orgKey = `${organism}||${taxonId}`;
    if (orgSel && Array.from(orgSel.options).some((o) => o.value === orgKey)) {
      orgSel.value = orgKey;
      orgSel.dispatchEvent(new Event("change", { bubbles: true }));
    }
    setTimeout(() => {
      if (sampleSel && sample && Array.from(sampleSel.options).some((o) => o.value === sample)) {
        sampleSel.value = sample;
        sampleSel.dispatchEvent(new Event("change", { bubbles: true }));
      }
      if (window.drawHistogram) window.drawHistogram();
    }, 40);
  }, 60);
}

// Summary detection-table state
// Group-by-Sample is on by default, so default-sort by sample for clean grouping.
let _sumSortCol = "Specimen ID";
let _sumSortAsc = true;
let _sumPage = 0;
const _sumPinned = new Set();
function _sumPageSize() {
  const sel = document.getElementById("summary-page-size");
  return sel ? parseInt(sel.value) || 0 : 25;
}
function _sumRowKey(r) {
  return `${r["Specimen ID"]}||${r["Detected Organism"]}||${r["Taxonomic ID #"]}`;
}

function drawSummary() {
  const fd = filteredData();
  const samples = uniq(fd.map((r) => r["Specimen ID"] || "")).filter(Boolean);

  // ── KPI cards ──────────────────────────────────────────────────────
  const totalInput = samples.reduce((s, sn) => s + (parseFloat((SAMPLE_META[sn] || {}).total_reads) || 0), 0);
  const totalOrg =
    samples.reduce((s, sn) => s + (parseFloat((SAMPLE_META[sn] || {}).total_organism_reads) || 0), 0) ||
    fd.reduce((s, r) => s + num(r["# Reads Aligned"]), 0);
  let pctClass = totalInput > 0 ? ((totalOrg / totalInput) * 100).toFixed(1) + "%" : "N/A";
  if (pctClass === "0.0%" && totalOrg > 0) pctClass = "<0.1%";

  const uniqOrgs = new Set(fd.map((r) => r["Taxonomic ID #"] || "").filter(Boolean)).size;
  const uniqGenera = new Set(fd.map((r) => r["Genus"]).filter(Boolean)).size;
  const hcCount = new Set(
    fd.filter((r) => isTruthy(r["High Consequence"])).map((r) => r["Taxonomic ID #"] || r["Detected Organism"]),
  ).size;
  const platforms = uniq(samples.map((sn) => (SAMPLE_META[sn] || {}).platform).filter((p) => p && p !== "unknown"));
  // Applied TASS cutoffs — per sample type when available, otherwise global fallback.
  const _kpiTypes = Array.from(
    new Set(DATA.map((r) => (r["Sample Type"] || "").trim().toLowerCase()).filter((t) => t && t !== "unknown")),
  ).sort();
  const _globalFallback = parseFloat((document.getElementById("filter-min") || {}).value) || 0;
  const _recCut = BEST_TASS_THRESH != null && !isNaN(BEST_TASS_THRESH) ? Number(BEST_TASS_THRESH) : null;
  let cutVal, cutSub, _typeVals;
  if (_kpiTypes.length > 0) {
    _typeVals = _kpiTypes.map((t) => ({
      type: t,
      applied: perTypeTass[t] != null ? perTypeTass[t] : _globalFallback,
      def: _defaultTassForType(t),
    }));
    const _appliedVals = _typeVals.map((tv) => tv.applied);
    const _allSame = _appliedVals.every((v) => v === _appliedVals[0]);
    if (_allSame) {
      cutVal = _appliedVals[0].toFixed(1);
      cutSub = _kpiTypes.length === 1 ? _kpiTypes[0] : "all types";
    } else {
      const _mn = Math.min(..._appliedVals),
        _mx = Math.max(..._appliedVals);
      cutVal = `${_mn.toFixed(0)}\u2013${_mx.toFixed(0)}`;
      cutSub = "by sample type";
    }
  } else {
    _typeVals = null;
    const appliedCut = _globalFallback;
    cutVal = appliedCut > 0 ? appliedCut.toFixed(1) : _recCut != null ? _recCut.toFixed(1) : "—";
    cutSub = _recCut != null ? `recommended ${_recCut.toFixed(1)}` : "applied filter";
  }

  const tReads = _fmtBig(totalInput);
  const aReads = _fmtBig(totalOrg);
  const cards = [
    {
      label: "Samples",
      value: samples.length,
      sub: platforms.length ? platforms.join(", ") : "",
      tipId: "kpi-tip-samples",
    },
    {
      label: "Total Reads",
      value: totalInput > 0 ? tReads.short : "N/A",
      full: tReads.full,
      sub: "input reads",
      tipId: "kpi-tip-reads",
    },
    {
      label: "% Classified",
      value: pctClass,
      sub: aReads.short + " aligned",
      subFull: aReads.full,
      tipId: "kpi-tip-pct",
    },
    { label: "Organisms", value: _fmtInt(uniqOrgs), sub: "unique taxa" },
    { label: "Genera", value: _fmtInt(uniqGenera), sub: "unique genera" },
    { label: "High Consequence", value: _fmtInt(hcCount), sub: "flagged pathogens" },
    { label: "TASS Cutoff", value: cutVal, sub: cutSub, tipId: "kpi-tip-tass" },
  ];
  const hcColor = hcCount > 0 ? "#c62828" : "#1565c0";
  const kpiRow = document.getElementById("summary-kpi-row");
  const infoIcon = `<i class="fas fa-circle-info" style="font-size:0.7em;opacity:0.45;margin-left:4px;vertical-align:middle"></i>`;
  kpiRow.innerHTML = cards
    .map(
      (c) =>
        `<div class="kpi-card${c.label === "High Consequence" && hcCount > 0 ? " kpi-hc" : ""}"${
          c.label === "High Consequence" ? ` style="border-left-color:${hcColor}"` : ""
        }${c.tipId ? ` id="${c.tipId}" style="cursor:help"` : ""}>` +
        `<div class="kpi-label">${c.label}${c.tipId || c.label === "High Consequence" ? infoIcon : ""}</div>` +
        `<div class="kpi-value"${c.full ? ` title="${c.full}"` : ""}>${c.value}</div>` +
        (c.sub ? `<div class="kpi-sub"${c.subFull ? ` title="${c.subFull}"` : ""}>${c.sub}</div>` : "") +
        `</div>`,
    )
    .join("");

  // High-consequence KPI hover: keep this compact so it cannot cover the
  // top of the summary page when many samples have flagged detections.
  const hcCard = kpiRow.querySelector(".kpi-hc");
  if (hcCard) {
    const byOrg = new Map();
    const hcSamples = new Set();
    fd.filter((r) => isTruthy(r["High Consequence"])).forEach((r) => {
      const org = r["Detected Organism"] || r["Taxonomic ID #"] || "Unknown";
      if (!byOrg.has(org)) byOrg.set(org, { samples: new Set(), maxTass: 0 });
      const entry = byOrg.get(org);
      if (r["Specimen ID"]) {
        entry.samples.add(r["Specimen ID"]);
        hcSamples.add(r["Specimen ID"]);
      }
      const tass = parseFloat(r["TASS Score"]) || 0;
      if (tass > entry.maxTass) entry.maxTass = tass;
    });
    // Sort organisms by descending sample count
    const orgRows = Array.from(byOrg.entries())
      .sort((a, b) => b[1].samples.size - a[1].samples.size)
      .map(
        ([org, entry]) =>
          `<tr>` +
          `<td style="text-align:right;padding-right:6px;font-weight:600;color:#ffb3b3">${entry.samples.size}</td>` +
          `<td style="text-align:left;padding-right:10px">${org}</td>` +
          `<td style="text-align:right;color:#ffd580;font-weight:600;white-space:nowrap">${
            entry.maxTass > 0 ? entry.maxTass.toFixed(1) : "—"
          }</td>` +
          `</tr>`,
      )
      .join("");
    const tip =
      `<b style="color:#ffb3b3">High-consequence detections</b><br>` +
      `<span style="font-size:0.85em;color:#ccc">${_fmtInt(byOrg.size)} unique organism${
        byOrg.size === 1 ? "" : "s"
      } &nbsp;·&nbsp; ${_fmtInt(hcSamples.size)} sample${hcSamples.size === 1 ? "" : "s"}</span>` +
      `<table style="margin-top:6px;border-collapse:collapse;width:100%">` +
      `<tr style="color:#aaa;font-size:0.8em"><td style="text-align:right;padding-right:6px"># Samples</td><td style="padding-right:10px">Organism</td><td style="text-align:right">Max TASS</td></tr>` +
      orgRows +
      `</table>`;
    hcCard.style.cursor = "help";
    hcCard.addEventListener("mouseover", (ev) => showTip(tip, ev));
    hcCard.addEventListener("mousemove", moveTip);
    hcCard.addEventListener("mouseout", hideTip);
  }

  // ── KPI tooltips: Samples ──────────────────────────────────────────
  const sampCard = document.getElementById("kpi-tip-samples");
  if (sampCard) {
    const sampRows = samples
      .map((sn) => {
        const m = SAMPLE_META[sn] || {};
        const reads = m.total_reads ? _fmtBig(m.total_reads).short : "—";
        const plat = m.platform && m.platform !== "unknown" ? m.platform : "";
        return (
          `<tr><td style="padding:2px 8px 2px 0;white-space:nowrap;color:#e0e0e0">${sn}</td>` +
          `<td style="padding:2px 8px 2px 0;text-align:right;color:#90caf9">${reads}</td>` +
          `<td style="padding:2px 0;color:#aaa;font-size:0.85em">${plat}</td></tr>`
        );
      })
      .join("");
    const sampTip =
      `<b style="color:#90caf9">Samples in view</b> <span style="color:#aaa;font-size:0.85em">(${samples.length})</span><br>` +
      `<table style="margin-top:5px;border-collapse:collapse">` +
      `<tr style="color:#888;font-size:0.78em"><td style="padding-right:8px">Name</td><td style="padding-right:8px;text-align:right">Reads</td><td>Platform</td></tr>` +
      sampRows +
      `</table>`;
    sampCard.addEventListener("mouseover", (ev) => showTip(sampTip, ev));
    sampCard.addEventListener("mousemove", moveTip);
    sampCard.addEventListener("mouseout", hideTip);
  }

  // ── KPI tooltips: Total Reads ──────────────────────────────────────
  const readsCard = document.getElementById("kpi-tip-reads");
  if (readsCard && totalInput > 0) {
    const rRows = samples
      .map((sn) => {
        const r = parseFloat((SAMPLE_META[sn] || {}).total_reads) || 0;
        const pct = totalInput > 0 ? ((r / totalInput) * 100).toFixed(1) + "%" : "—";
        return (
          `<tr><td style="padding:2px 8px 2px 0;white-space:nowrap;color:#e0e0e0">${sn}</td>` +
          `<td style="padding:2px 4px;text-align:right;color:#90caf9">${_fmtBig(r).short}</td>` +
          `<td style="padding:2px 0;color:#aaa;text-align:right">${pct}</td></tr>`
        );
      })
      .join("");
    const readsTip =
      `<b style="color:#90caf9">Total input reads</b><br>` +
      `<span style="color:#aaa;font-size:0.85em">Full count: ${tReads.full}</span><br>` +
      `<table style="margin-top:5px;border-collapse:collapse">` +
      `<tr style="color:#888;font-size:0.78em"><td style="padding-right:8px">Sample</td><td style="padding-right:4px;text-align:right">Reads</td><td style="text-align:right">Share</td></tr>` +
      rRows +
      `</table>`;
    readsCard.addEventListener("mouseover", (ev) => showTip(readsTip, ev));
    readsCard.addEventListener("mousemove", moveTip);
    readsCard.addEventListener("mouseout", hideTip);
  }

  // ── KPI tooltips: % Classified ─────────────────────────────────────
  const pctCard = document.getElementById("kpi-tip-pct");
  if (pctCard) {
    const pctRows = samples
      .map((sn) => {
        const m = SAMPLE_META[sn] || {};
        const inp = parseFloat(m.total_reads) || 0;
        const aln = parseFloat(m.total_organism_reads) || 0;
        const p = inp > 0 ? ((aln / inp) * 100).toFixed(1) + "%" : "—";
        return (
          `<tr><td style="padding:2px 8px 2px 0;white-space:nowrap;color:#e0e0e0">${sn}</td>` +
          `<td style="padding:2px 4px;text-align:right;color:#90caf9">${p}</td>` +
          `<td style="padding:2px 0;color:#aaa;text-align:right;font-size:0.85em">${_fmtBig(aln).short} / ${
            _fmtBig(inp).short
          }</td></tr>`
        );
      })
      .join("");
    const pctTip =
      `<b style="color:#90caf9">% Classified</b><br>` +
      `<span style="color:#aaa;font-size:0.85em">Reads assigned to any organism ÷ total input reads</span><br>` +
      `<table style="margin-top:5px;border-collapse:collapse">` +
      `<tr style="color:#888;font-size:0.78em"><td style="padding-right:8px">Sample</td><td style="padding-right:4px;text-align:right">% Class.</td><td style="text-align:right">Aligned / Total</td></tr>` +
      pctRows +
      `</table>`;
    pctCard.addEventListener("mouseover", (ev) => showTip(pctTip, ev));
    pctCard.addEventListener("mousemove", moveTip);
    pctCard.addEventListener("mouseout", hideTip);
  }

  // ── KPI tooltips: TASS Cutoff ──────────────────────────────────────
  const tassCard = document.getElementById("kpi-tip-tass");
  if (tassCard) {
    let tassTip =
      `<b style="color:#90caf9">TASS Cutoff</b><br>` +
      `<span style="color:#ccc;font-size:0.88em">Taxonomic Assignment Scoring System threshold.<br>` +
      `Detections below this score are hidden across all tabs.</span><br><br>`;
    if (_typeVals && _typeVals.length > 0) {
      tassTip +=
        `<table style="border-collapse:collapse;font-size:0.88em;width:100%">` +
        `<tr>` +
        `<th style="color:#90caf9;font-weight:500;text-align:left;padding-bottom:4px;padding-right:12px">Type</th>` +
        `<th style="color:#90caf9;font-weight:500;text-align:right;padding-bottom:4px;padding-right:12px">Applied</th>` +
        `<th style="color:#90caf9;font-weight:500;text-align:right;padding-bottom:4px">Default</th>` +
        `</tr>` +
        _typeVals
          .map((tv) => {
            const changed = Math.round(tv.applied) !== Math.round(tv.def);
            return (
              `<tr>` +
              `<td style="color:#fff;text-transform:capitalize;padding-right:12px;padding-bottom:2px">${tv.type}</td>` +
              `<td style="color:${
                changed ? "#ffb74d" : "#a5d6a7"
              };font-weight:700;text-align:right;padding-right:12px">${tv.applied.toFixed(0)}</td>` +
              `<td style="color:#78909c;text-align:right">${tv.def.toFixed(0)}</td>` +
              `</tr>`
            );
          })
          .join("") +
        `</table>` +
        `<div style="margin-top:6px;font-size:0.78em;color:#607d8b">Orange = modified from default · Green = at default</div>`;
    } else {
      tassTip +=
        `<table style="border-collapse:collapse;font-size:0.88em">` +
        `<tr><td style="color:#aaa;padding-right:8px">Applied cutoff</td><td style="color:#fff;font-weight:700">${cutVal}</td></tr>` +
        (_recCut != null
          ? `<tr><td style="color:#aaa;padding-right:8px">Recommended</td><td style="color:#90caf9;font-weight:700">${_recCut.toFixed(
              1,
            )}</td></tr>`
          : "") +
        `</table>`;
    }
    tassCard.addEventListener("mouseover", (ev) => showTip(tassTip, ev));
    tassCard.addEventListener("mousemove", moveTip);
    tassCard.addEventListener("mouseout", hideTip);
  }

  // ── Detections table (sortable / paginated / groupable) ────────────
  _renderSummaryTable(fd);

  // ── Genera comparison stacked bar chart ────────────────────────────
  _drawSummaryGenera(fd, samples);

  // ── Annotation summary (only if protein/VF-AMR enabled) ────────────
  _drawSummaryAnnotation(fd);

  // ── Cross-sample organism analysis (aggregate / PCA / co-occurrence)
  _drawCrossSample(fd, samples);
}

const _SUM_COLS = [
  {
    key: "Detected Organism",
    label: "Organism",
    align: "left",
    tip: "<b>Organism</b><br>Detected organism name for this row at the current view level (Strain/Species/Genus).<br><br><span style='color:#ccc'>Includes badges for risk/level/rescue context. Click row to pin for comparison.</span>",
  },
  {
    key: "_spark",
    label: "Breadth of Coverage",
    align: "center",
    nosort: true,
    tip: "<b>Breadth of Coverage</b><br>Mini coverage sparkline for this detection. It summarizes where reads cover the reference genome (position-based breadth).<br><br><span style='color:#ccc'>Hover the sparkline for quick distribution context. Click it to jump to the Histogram tab for full detail.</span>",
  },
  {
    key: "Specimen ID",
    label: "Sample",
    align: "left",
    tip: "<b>Sample</b><br>Specimen identifier this detection belongs to.<br><br><span style='color:#ccc'>Use with grouping and filters to compare organisms within or across samples.</span>",
  },
  {
    key: "Microbial Category",
    label: "Category",
    align: "left",
    tip: "<b>Category</b><br>Microbial category assigned to the detection (for triage/risk grouping).<br><br><span style='color:#ccc'>Useful for narrowing views to clinically relevant classes.</span>",
  },
  {
    key: "TASS Score",
    label: "TASS",
    align: "right",
    num: true,
    tip: "<b>TASS Score</b><br>Taxonomic Assignment Scoring System confidence score for the row (higher is stronger evidence).<br><br><span style='color:#ccc'>This is the primary score used by the cutoff filter. Rows below cutoff are hidden unless rescued by rollup logic or VF/AMR exceptions.</span>",
  },
  {
    key: "_rescue",
    label: "Rescue",
    align: "center",
    nosort: true,
    tip: "<b>Rescue</b><br>Shows whether a row passes on its own or is kept via rollup.<br><br><span style='color:#ccc'><b>pass</b>: row TASS passes cutoff.<br><b>↑ species / ↑ genus</b>: row TASS is below cutoff, but parent Species/Genus TASS passes.<br><b>below</b>: no rescue at current threshold.</span>",
  },
  {
    key: "_vfamr",
    label: "VF / AMR",
    align: "center",
    nosort: true,
    tip: "<b>VF / AMR</b><br>Presence/count summary of virulence-factor and antimicrobial-resistance gene hits linked to this detection's genus/sample context.<br><br><span style='color:#ccc'>Hover cells to see gene details. Click to jump to the VF/AMR view filtered to the relevant genus.</span>",
  },
  {
    key: "_novelty",
    label: "Novelty",
    align: "center",
    nosort: true,
    tip: "<b>Novelty</b><br>Indicates whether this organism's genus or species has reference-free novelty signal in the same sample, based on the Novelty tab classifier.<br><br><span style='color:#ccc'><b>● gen</b>: genus-level novelty candidate detected.<br><b>● sp</b>: species-level novelty candidate detected.<br><b>—</b>: no novelty signal for this lineage.<br><br>Matched by name (taxid when available). Open the Novelty tab for full details.</span>",
  },
  {
    key: "% Reads",
    label: "% Reads",
    align: "right",
    num: true,
    tip: "<b>% Reads</b><br>Percent of reads in the sample attributed to this detection row.<br><br><span style='color:#ccc'>Higher values indicate larger representation in that sample.</span>",
  },
  {
    key: "# Reads Aligned",
    label: "Reads",
    align: "right",
    num: true,
    tip: "<b>Reads</b><br>Absolute read count aligned/assigned to this row.<br><br><span style='color:#ccc'>Pair with % Reads to compare abundance across different sequencing depths.</span>",
  },
  {
    key: "Breadth %",
    label: "Breadth %",
    align: "right",
    num: true,
    tip: "<b>Breadth %</b><br>Percent of the reference breadth covered by reads for this row.<br><br><span style='color:#ccc'>Higher breadth generally increases confidence that the organism signal is real.</span>",
  },
  {
    key: "Sample Type",
    label: "Sample Type",
    align: "left",
    tip: "<b>Sample Type</b><br>Body site / specimen-type metadata for the sample.<br><br><span style='color:#ccc'>Can change applied TASS cutoff when per-type thresholds are enabled.</span>",
  },
  {
    key: "_kw",
    label: "Keywords",
    align: "left",
    nosort: true,
    tip: "<b>Keywords</b><br>Quick tags extracted from organism metadata and annotations (e.g., AMR, virulence, respiratory).<br><br><span style='color:#ccc'>A rapid scan aid for triage and communication.</span>",
  },
];

function _summaryRescueCellHTML(row, rowKey) {
  const pi = rowPassInfo(row);
  if (!pi || isNaN(pi.strain)) return '<span style="color:#adb5bd">-</span>';

  const strainTxt = isNaN(pi.strain) ? "?" : pi.strain.toFixed(1);
  const speciesTxt = isNaN(pi.species) ? "?" : pi.species.toFixed(1);
  const genusTxt = isNaN(pi.genus) ? "?" : pi.genus.toFixed(1);
  const rkAttr = rowKey ? ` data-rk="${String(rowKey).replace(/"/g, "&quot;")}"` : "";

  if (pi.rescued && pi.rescueLevel === "genus") {
    return `<span class="sum-rescue-chip"${rkAttr} style="display:inline-block;font-size:9px;font-weight:700;padding:0 5px;border-radius:10px;line-height:1.5;background:#f3e5f5;color:#6a1b9a;cursor:help"
            title="Rescued at genus level: row TASS ${strainTxt} is below cutoff ${pi.thr}, but Genus TASS ${genusTxt} passes. Species TASS ${speciesTxt}.">&#x2191; genus</span>`;
  }
  if (pi.rescued) {
    return `<span class="sum-rescue-chip"${rkAttr} style="display:inline-block;font-size:9px;font-weight:700;padding:0 5px;border-radius:10px;line-height:1.5;background:#ffe0b2;color:#c2410c;cursor:help"
            title="Rescued at species level: row TASS ${strainTxt} is below cutoff ${pi.thr}, but Species TASS ${speciesTxt} passes. Genus TASS ${genusTxt}.">&#x2191; species</span>`;
  }
  if (pi.strainPass) {
    return `<span class="sum-rescue-chip"${rkAttr} style="display:inline-block;font-size:9px;font-weight:700;padding:0 5px;border-radius:10px;line-height:1.5;background:#e8f5e9;color:#2e7d32;cursor:help"
            title="Passes on own evidence: row TASS ${strainTxt} is at or above cutoff ${pi.thr}.">pass</span>`;
  }
  return `<span class="sum-rescue-chip"${rkAttr} style="display:inline-block;font-size:9px;font-weight:700;padding:0 5px;border-radius:10px;line-height:1.5;background:#ffebee;color:#c62828;cursor:help"
          title="Below cutoff: row TASS ${strainTxt}, cutoff ${pi.thr}. Species TASS ${speciesTxt}; Genus TASS ${genusTxt}.">below</span>`;
}

function _summaryRescueTipHTML(row) {
  const pi = rowPassInfo(row);
  if (!pi || isNaN(pi.strain)) return null;

  const sName = row["Specimen ID"] || "";
  const org = row["Detected Organism"] || "Unknown";
  const lvl = row["Level"] || "Strain";
  const taxid = row["Taxonomic ID #"] || "";
  const strain = isNaN(pi.strain) ? 0 : pi.strain;
  const species = isNaN(pi.species) ? 0 : pi.species;
  const genus = isNaN(pi.genus) ? 0 : pi.genus;
  const status = pi.rescued
    ? `Rescued by ${pi.rescueLevel}`
    : pi.strainPass
    ? "Passes on own evidence"
    : "Below cutoff";
  const statusColor = pi.rescued ? "#c2410c" : pi.strainPass ? "#2e7d32" : "#c62828";
  const maxDen = Math.max(1, pi.thr, strain, species, genus);
  const bar = (v, color) =>
    `<span style="display:inline-block;width:64px;height:6px;background:#eceff1;border-radius:3px;overflow:hidden;vertical-align:middle;margin-left:6px"><span style="display:inline-block;height:100%;width:${Math.min(
      100,
      (Math.max(0, v) / maxDen) * 100,
    ).toFixed(1)}%;background:${color}"></span></span>`;

  let html =
    `<div style="font-weight:700;margin-bottom:4px;color:#90caf9">Rescue Summary</div>` +
    `<div style="margin-bottom:5px"><i>${_novEsc(org)}</i><br><span style="color:#b0bec5">${_novEsc(
      sName,
    )}</span><span style="color:#90a4ae"> · ${_novEsc(lvl)}${
      taxid ? ` · taxid ${_novEsc(String(taxid))}` : ""
    }</span></div>` +
    `<div style="font-size:0.9em;margin-bottom:6px"><b style="color:${statusColor}">${status}</b> · cutoff <b>${pi.thr.toFixed(
      1,
    )}</b></div>` +
    `<table style="border-collapse:collapse;font-size:0.82em;width:100%">` +
    `<tr><td style="padding:1px 8px 1px 0;color:#ffd580">Strain</td><td style="padding:1px 0;text-align:right">${strain.toFixed(
      1,
    )}${bar(strain, "#1971c2")}</td></tr>` +
    `<tr><td style="padding:1px 8px 1px 0;color:#ffd580">Species</td><td style="padding:1px 0;text-align:right">${species.toFixed(
      1,
    )}${bar(species, "#c2410c")}</td></tr>` +
    `<tr><td style="padding:1px 8px 1px 0;color:#ffd580">Genus</td><td style="padding:1px 0;text-align:right">${genus.toFixed(
      1,
    )}${bar(genus, "#6a1b9a")}</td></tr>` +
    `</table>`;

  // Row-linked novelty context (per relevant genus), so values differ per row.
  const genusNameRaw = (row["Genus Name"] || row["Genus"] || "").trim();
  const genusKey = genusNameRaw.toLowerCase();
  if (genusKey) {
    const closed = _novClosedGenus(sName) || {};
    const c = closed[genusKey] || null;
    const cand = (
      ((_novSamples()[sName] || {}).candidates || []).filter((x) => (x.rank || "").toLowerCase() === "genus") || []
    ).find((x) => ((x.name || "").trim().toLowerCase() || "") === genusKey);
    const bucket =
      c && cand ? "Both" : c ? "Aligned only" : cand ? `${_novClsShort()} only` : "No genus novelty evidence";
    const bucketColor = c && cand ? "#1864ab" : c ? "#d9480f" : cand ? "#2b8a3e" : "#9e9e9e";
    const cReads = c ? +c.reads || 0 : 0;
    const cPct = c ? +c.pct || 0 : 0;
    const mmReads = cand ? +cand.reads || 0 : 0;
    const mmFracS = cand ? (+cand.frac_of_sample || 0) * 100 : 0;
    const mmFracHr = cand ? (+cand.frac_of_highrank || 0) * 100 : 0;
    html +=
      `<div style="margin-top:6px;padding-top:6px;border-top:1px solid rgba(255,255,255,0.16)">` +
      `<div style="font-weight:700;color:#90caf9;margin-bottom:3px">Row-linked Genus Context</div>` +
      `<div style="font-size:0.8em;color:#cfd8dc;margin-bottom:3px"><i>${_novEsc(
        genusNameRaw,
      )}</i> · <b style="color:${bucketColor}">${_novEsc(bucket)}</b></div>` +
      `<table style="border-collapse:collapse;font-size:0.8em;width:100%">` +
      `<tr><td style="padding:1px 8px 1px 0;color:#90caf9">Aligned reads</td><td style="text-align:right">${
        c ? cReads.toLocaleString() : "-"
      }</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#90caf9">Aligned % reads</td><td style="text-align:right">${
        c ? cPct.toFixed(3) + "%" : "-"
      }</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#f59f00">${_novEsc(
        _novClsShort(),
      )} ${_novUnit()}</td><td style="text-align:right">${cand ? mmReads.toLocaleString() : "-"}</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#f59f00">${_novEsc(
        _novClsShort(),
      )} % of ${_novUnit()}</td><td style="text-align:right">${cand ? mmFracS.toFixed(3) + "%" : "-"}</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#f59f00">${_novEsc(
        _novClsShort(),
      )} % genus+ pool</td><td style="text-align:right">${cand ? mmFracHr.toFixed(1) + "%" : "-"}</td></tr>` +
      `</table>` +
      `</div>`;
  }

  const ns = _novSamples();
  const nSm = (ns && ns[sName] && ns[sName].summary) || null;
  if (nSm && typeof _novMcFractions === "function") {
    const fr = _novMcFractions(nSm);
    const pct = (x) => `${((x || 0) * 100).toFixed(1)}%`;
    html +=
      `<div style="margin-top:6px;padding-top:6px;border-top:1px solid rgba(255,255,255,0.16)">` +
      `<div style="font-weight:700;color:#90caf9;margin-bottom:3px">Sample-level Method Split</div>` +
      `<table style="border-collapse:collapse;font-size:0.8em;width:100%">` +
      `<tr><td style="padding:1px 8px 1px 0;color:#90caf9">Aligned</td><td style="text-align:right">${pct(
        fr.align,
      )}</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#2f9e44">${_novEsc(
        _novClsShort(),
      )} (all)</td><td style="text-align:right">${pct((fr.mmSp || 0) + (fr.mmHr || 0))}</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#f59f00">${_novEsc(
        _novClsShort(),
      )} genus+</td><td style="text-align:right">${pct(fr.mmHr)}</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#9e9e9e">Kraken2 only</td><td style="text-align:right">${pct(
        fr.k2only,
      )}</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#b71c1c">Dark matter</td><td style="text-align:right">${pct(
        fr.dark,
      )}</td></tr>` +
      `</table>` +
      `</div>`;
  }
  return html;
}

function _renderSummaryTable(fd) {
  const wrap = document.getElementById("summary-top-wrap");
  // Append faded below-cutoff rows that carry VF/AMR hits (toggle-gated).
  const _belowExtra = _belowCutoffExtraRows();
  if (_belowExtra.length) fd = fd.concat(_belowExtra);
  // Append faded sub-threshold (aligned) + no-alignment novelty rows (toggle-gated).
  const _novSubExtra = _noveltySubThresholdExtraRows();
  if (_novSubExtra.length) fd = fd.concat(_novSubExtra);
  // Append indicator rows for samples with VF/AMR hits but no detection rows.
  const _vfOnlySum = _vfamrOnlySampleRows(fd);
  if (_vfOnlySum.length) fd = fd.concat(_vfOnlySum);
  // Append novelty-only indicator rows for samples with novelty genus evidence.
  const _novOnlySum = _noveltyOnlySampleRows(fd);
  if (_novOnlySum.length) fd = fd.concat(_novOnlySum);
  // Append empty-sample indicator rows for samples with no detections at all.
  const _emptyOnlySum = _emptyOnlySampleRows(fd);
  if (_emptyOnlySum.length) fd = fd.concat(_emptyOnlySum);
  // Show rollup hint only when rescued rows are actually present
  const _sumRollupHint = document.getElementById("sum-rollup-hint");
  if (_sumRollupHint) {
    const _hasRescued =
      ROLLUP_PASS &&
      fd.some((r) => {
        const p = rowPassInfo(r);
        return p && p.rescued;
      });
    _sumRollupHint.style.display = _hasRescued ? "" : "none";
  }
  // Legend: show all three badges when rollup on; single current-level badge when off
  const _legendRollup = document.getElementById("sum-level-legend-rollup");
  const _legendSimple = document.getElementById("sum-level-legend-simple");
  if (_legendRollup) _legendRollup.style.display = ROLLUP_PASS ? "" : "none";
  if (_legendSimple) {
    _legendSimple.style.display = ROLLUP_PASS ? "none" : "";
    if (!ROLLUP_PASS) {
      const _vl = document.getElementById("view-level")?.value || "Strain";
      const _lvlStyles = {
        Strain: ["#e3f2fd", "#1565c0"],
        Species: ["#f3e5f5", "#6a1b9a"],
        Genus: ["#e8f5e9", "#2e7d32"],
      };
      const _ls = _lvlStyles[_vl] || ["#f5f5f5", "#555"];
      _legendSimple.innerHTML = `<span style="display:inline-block;font-size:9px;font-weight:700;padding:0 4px;border-radius:3px;vertical-align:middle;line-height:1.5;background:${_ls[0]};color:${_ls[1]};border:1px solid ${_ls[1]}44">${_vl}</span> = taxonomic level of each row ·`;
    }
  }
  if (!fd.length) {
    wrap.innerHTML = '<p style="color:#888;padding:1em">No organisms match the current filters.</p>';
    _updateSummaryPager(0, 0, 0, 0);
    return;
  }
  const grouped = document.getElementById("summary-group-sample")?.checked;
  const _sumViewLevel = document.getElementById("view-level")?.value || "Strain";
  const _visSumCols = _SUM_COLS.filter(
    (c) =>
      !(c.key === "_spark" && _sumViewLevel !== "Strain") &&
      !(c.key === "_vfamr" && !HAS_PROT) &&
      !(c.key === "_novelty" && !HAS_NOVELTY),
  );

  // Sort
  const rows = [...fd];
  const col = _sumSortCol;
  const isNum = (_SUM_COLS.find((c) => c.key === col) || {}).num;
  const cmp = (a, b) => {
    if (isNum) {
      const d = num(a[col]) - num(b[col]);
      return _sumSortAsc ? d : -d;
    }
    const d = String(a[col] || "").localeCompare(String(b[col] || ""));
    return _sumSortAsc ? d : -d;
  };
  const _sumGrpIdx = _sampleOrder.length ? Object.fromEntries(_sampleOrder.map((id, i) => [id, i])) : {};
  rows.sort((a, b) => {
    // When grouping, order groups by _sampleOrder (right-panel order), then sort within.
    if (grouped) {
      const sia = String(a["Specimen ID"] || "");
      const sib = String(b["Specimen ID"] || "");
      const ia = _sumGrpIdx[sia] !== undefined ? _sumGrpIdx[sia] : 9999;
      const ib = _sumGrpIdx[sib] !== undefined ? _sumGrpIdx[sib] : 9999;
      if (ia !== ib) return ia - ib;
      if (col !== "Specimen ID") return cmp(a, b);
      return num(b["TASS Score"]) - num(a["TASS Score"]);
    }
    return cmp(a, b);
  });

  // Pagination over the (optionally grouped) row list
  const pageSize = _sumPageSize();
  const total = rows.length;
  const pages = pageSize > 0 ? Math.max(1, Math.ceil(total / pageSize)) : 1;
  _sumPage = Math.max(0, Math.min(_sumPage, pages - 1));
  const start = pageSize > 0 ? _sumPage * pageSize : 0;
  const end = pageSize > 0 ? Math.min(start + pageSize, total) : total;
  const slice = rows.slice(start, end);

  // Header
  let html = "<table><thead><tr>";
  _visSumCols.forEach((c) => {
    const sortable = !c.nosort;
    const ind = c.key === _sumSortCol ? `<span class="sort-ind">${_sumSortAsc ? "▲" : "▼"}</span>` : "";
    const tipAttr = c.tip ? ` data-sum-tip="${String(c.tip).replace(/"/g, "&quot;")}"` : "";
    const titleAttr = c.tip ? ` title="${String(c.label).replace(/"/g, "&quot;")} column details"` : "";
    html += `<th class="${sortable ? "sortable" : ""}" data-col="${c.key}"${tipAttr}${titleAttr} style="text-align:${
      c.align
    }">${c.label}${ind}</th>`;
  });
  html += "</tr></thead><tbody>";

  let lastSample = null;
  slice.forEach((r) => {
    if (grouped && r["Specimen ID"] !== lastSample) {
      lastSample = r["Specimen ID"];
      const _grpSwatchColor = sampleColors[lastSample] || "#90a4ae";
      html += `<tr class="grp-row"><td colspan="${_visSumCols.length}"> <span style="display:inline-block;width:10px;height:10px;border-radius:2px;background:${_grpSwatchColor};vertical-align:middle;margin:0 5px 1px 0;border:1px solid rgba(0,0,0,0.18);flex-shrink:0"></span>${lastSample}</td></tr>`;
    }
    // VF/AMR-only indicator row (sample has no passing detections).
    if (r.__vfamrOnly) {
      html +=
        `<tr class="vfamr-only-row" data-vfamr-sample="${String(r["Specimen ID"]).replace(/"/g, "&quot;")}">` +
        `<td colspan="${_visSumCols.length}">${_vfamrOnlyMessageHTML(r["Specimen ID"], r.__vfamr)}</td></tr>`;
      return;
    }
    // Novelty-only indicator row (sample has no passing detections).
    if (r.__noveltyOnly) {
      html +=
        `<tr class="novelty-only-row" data-novelty-sample="${String(r["Specimen ID"]).replace(/"/g, "&quot;")}">` +
        `<td colspan="${_visSumCols.length}">${_noveltyOnlyMessageHTML(r["Specimen ID"], r.__novelty)}</td></tr>`;
      return;
    }
    // Empty-sample indicator row (no detections and no VF/AMR or novelty evidence).
    if (r.__emptyOnly) {
      html +=
        `<tr class="empty-only-row" data-empty-sample="${String(r["Specimen ID"]).replace(/"/g, "&quot;")}">` +
        `<td colspan="${_visSumCols.length}">${_emptyOnlyMessageHTML(r["Specimen ID"], r.__emptyMeta)}</td></tr>`;
      return;
    }
    const cat = r["Microbial Category"] || "Unknown";
    const cc = _CAT_COLORS[cat] || "#555";
    const hc = isTruthy(r["High Consequence"]);
    const key = _sumRowKey(r);
    const hist = _histFor(r);
    const _pi_s = rowPassInfo(r);
    const _rescued_s = _pi_s && _pi_s.rescued;
    const _rescueBadge = _rescued_s
      ? `<span style="display:inline-block;margin-left:5px;font-size:9px;font-weight:700;` +
        `padding:0 4px;border-radius:3px;vertical-align:middle;line-height:1.5;` +
        `background:#ffe0b2;color:#c2410c;" ` +
        `title="Strain TASS ${_pi_s.strain.toFixed(1)} is below the cutoff (${
          _pi_s.thr
        }) but kept visible because its ${_pi_s.rescueLevel} aggregation passes. ` +
        `Turn off &quot;Roll up threshold&quot; to hide it.">&#x2191; ${_pi_s.rescueLevel}</span>`
      : "";
    html +=
      `<tr data-key="${key}" class="${hc ? "hc-row" : ""}${_rescued_s ? " rescued-row" : ""}${
        r.__belowCutoffVFAMR ? " below-cutoff-row" : ""
      }${r.__belowCutoffAligned ? " below-cutoff-aligned-row" : ""}${r.__noveltyNoAlign ? " novelty-noalign-row" : ""}${
        _sumPinned.has(key) ? " row-pinned" : ""
      }">` +
      `<td style="position:relative;padding-right:70px;padding-left:4px;white-space:nowrap;overflow:hidden;">` +
      `<span style="display:inline-block;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;vertical-align:middle;max-width:100%">${_orgBadges(
        r,
      )}<i>${r["Detected Organism"] || ""}</i>${_ncbiLink(r)}${_rescueBadge}${_belowCutoffBadgeHTML(
        r,
      )}${_subThresholdBadgeHTML(r)}</span>` +
      // Star sits inline, right after the organism text.
      `${_watchStarHTML(r, true, null, true)}` +
      // Right-anchored group: level badge (flush right), then the pin hint
      // which only appears on hover/pin and shifts the badge left.
      `<span style="position:absolute;right:6px;top:50%;transform:translateY(-50%);display:inline-flex;align-items:center;white-space:nowrap;">${
        ROLLUP_PASS
          ? (function () {
              const _lvl = r["Level"] || "Strain";
              const _s = {
                Strain: ["#e3f2fd", "#1565c0"],
                Species: ["#f3e5f5", "#6a1b9a"],
                Genus: ["#e8f5e9", "#2e7d32"],
              }[_lvl] || ["#f5f5f5", "#555"];
              return `<span style="display:inline-block;font-size:9px;font-weight:700;padding:0 4px;border-radius:3px;line-height:1.5;background:${
                _s[0]
              };color:${_s[1]};border:1px solid ${_s[1]}44" title="${_lvl}-level ${
                _lvl === "Strain" ? "detection" : "aggregate"
              }">${_lvl}</span>`;
            })()
          : ""
      }<span class="row-pin-hint"></span></span></td>` +
      (_visSumCols.some((c) => c.key === "_spark")
        ? `<td style="text-align:center">${
            hist
              ? `<span class="spark-cell" data-spark="${key}">${_sparkHtml(hist)}</span>`
              : '<span style="color:#ccc">—</span>'
          }</td>`
        : "") +
      `<td>${r["Specimen ID"] || ""}</td>` +
      `<td><span style="color:${cc};font-weight:600">${cat}</span></td>` +
      `<td style="text-align:right">${num(r["TASS Score"]).toFixed(1)}</td>` +
      `<td style="text-align:center">${_summaryRescueCellHTML(r, key)}</td>` +
      (HAS_PROT ? `<td style="text-align:center">${_vfamrCellHTML(r)}</td>` : "") +
      (HAS_NOVELTY ? `<td style="text-align:center">${_noveltyCellHTML(r)}</td>` : "") +
      `<td style="text-align:right">${num(r["% Reads"]).toFixed(2)}</td>` +
      `<td style="text-align:right" title="${_fmtInt(r["# Reads Aligned"])}">${
        _fmtBig(r["# Reads Aligned"]).short
      }</td>` +
      `<td style="text-align:right">${num(r["Breadth %"]).toFixed(1)}</td>` +
      `<td>${
        r["Sample Type"] && r["Sample Type"] !== "unknown" ? r["Sample Type"] : '<span style="color:#ccc">—</span>'
      }</td>` +
      `<td>${_rowKeywords(r)}</td>` +
      "</tr>";
  });
  html += "</tbody></table>";
  wrap.innerHTML = html;

  // Wire header sort
  wrap.querySelectorAll("th.sortable").forEach((th) => {
    th.addEventListener("click", () => {
      const c = th.dataset.col;
      if (_sumSortCol === c) _sumSortAsc = !_sumSortAsc;
      else {
        _sumSortCol = c;
        _sumSortAsc = false;
      }
      _sumPage = 0;
      _renderSummaryTable(filteredData());
    });
  });
  // Column help tooltips for key summary headers.
  wrap.querySelectorAll("thead th[data-sum-tip]").forEach((th) => {
    th.style.cursor = "help";
    th.addEventListener("mouseover", (ev) => showTip(th.dataset.sumTip, ev));
    th.addEventListener("mousemove", moveTip);
    th.addEventListener("mouseout", hideTip);
  });
  // Novelty chip hover tooltips in the summary table body.
  wrap.querySelectorAll("tbody span[data-sum-tip]").forEach((el) => {
    el.addEventListener("mouseover", (ev) => showTip(el.dataset.sumTip, ev));
    el.addEventListener("mousemove", moveTip);
    el.addEventListener("mouseout", hideTip);
  });
  // Wire row pin (click) — ignore clicks on the spark, group rows, links, and star icon
  wrap.querySelectorAll("tbody tr[data-key]").forEach((tr) => {
    tr.addEventListener("click", (e) => {
      if (
        e.target.closest(".spark") ||
        e.target.closest("[data-spark]") ||
        e.target.closest("a") ||
        e.target.closest(".watch-star") ||
        e.target.closest(".vfamr-cell")
      )
        return;
      const k = tr.dataset.key;
      if (_sumPinned.has(k)) {
        _sumPinned.delete(k);
        tr.classList.remove("row-pinned");
      } else {
        _sumPinned.add(k);
        tr.classList.add("row-pinned");
      }
      _updateSumPinBar();
    });
  });
  // Wire spark hover preview + click-to-open in Histogram tab
  wrap.querySelectorAll("[data-spark]").forEach((sp) => {
    const key = sp.dataset.spark;
    const cd = _histMap().get(key);
    if (!cd) return;
    const parts = key.split("||");
    const r = { "Detected Organism": parts[1], "Specimen ID": parts[0], "Taxonomic ID #": parts[2] };
    sp.style.cursor = "pointer";
    sp.addEventListener("mouseover", (ev) => showTip(_histPreviewHtml(r, cd), ev));
    sp.addEventListener("mousemove", moveTip);
    sp.addEventListener("mouseout", hideTip);
    sp.addEventListener("click", (e) => {
      e.stopPropagation();
      hideTip();
      _jumpToHistogram(parts[1], parts[2], parts[0]);
    });
  });
  // Strain-breakdown hover tooltip for Species / Genus rows in the summary table.
  // Build a key→row lookup from the current page slice so we can resolve data on hover.
  const _sumSliceMap = new Map();
  slice.forEach(function (r) {
    _sumSliceMap.set(_sumRowKey(r), r);
  });
  wrap.querySelectorAll("tbody tr[data-key]").forEach(function (tr) {
    // skip group-header rows (they have no data-key — already filtered by selector)
    const _sumR = _sumSliceMap.get(tr.dataset.key);
    if (!_sumR) return;
    const _sumLvl = _sumR["Level"] || "Strain";
    if (_sumLvl !== "Species" && _sumLvl !== "Genus") return;
    const _sumOrgTd = tr.querySelector("td:first-child");
    if (!_sumOrgTd) return;
    _sumOrgTd.style.cursor = "help";
    _sumOrgTd.addEventListener("mouseenter", function (e) {
      const _tipStrains = _getTopStrainsForRow(_sumR);
      if (_tipStrains.length) {
        const _tipHtml = _strainTooltipHtml(_tipStrains, _sumLvl, _sumR["Detected Organism"]);
        if (_tipHtml) showTip(_tipHtml, e);
      }
    });
    _sumOrgTd.addEventListener("mousemove", moveTip);
    _sumOrgTd.addEventListener("mouseleave", hideTip);
  });

  // Rescue badge hover: compact rescue diagnostics + trimmed novelty-style split.
  wrap.querySelectorAll(".sum-rescue-chip[data-rk]").forEach(function (chip) {
    const _r = _sumSliceMap.get(chip.getAttribute("data-rk"));
    if (!_r) return;
    chip.addEventListener("mouseover", function (ev) {
      const tip = _summaryRescueTipHTML(_r);
      if (tip) showTip(tip, ev);
    });
    chip.addEventListener("mousemove", moveTip);
    chip.addEventListener("mouseout", hideTip);
  });

  // VF / AMR cell: hover → gene-list tooltip, click → jump to VF/AMR tab
  // filtered by the row's genus.
  wrap.querySelectorAll("tbody tr[data-key] .vfamr-cell[data-vfamr]").forEach(function (cell) {
    const tr = cell.closest("tr[data-key]");
    const _vfR = tr ? _sumSliceMap.get(tr.dataset.key) : null;
    if (!_vfR) return;
    cell.addEventListener("mouseover", function (ev) {
      ev.stopPropagation();
      const tip = _vfamrTip(_vfR);
      if (tip) showTip(tip, ev);
    });
    cell.addEventListener("mousemove", moveTip);
    cell.addEventListener("mouseout", hideTip);
    cell.addEventListener("click", function (e) {
      e.stopPropagation();
      hideTip();
      const g = cell.dataset.vfamrGenus || _vfR["Genus"] || "";
      if (g) _jumpToProteins(g);
    });
  });

  // VF/AMR-only indicator rows: hover → gene list, "View VF/AMR" → tab.
  const _vfOnlyMap = new Map();
  slice.forEach((r) => {
    if (r.__vfamrOnly) _vfOnlyMap.set(String(r["Specimen ID"]), r.__vfamr);
  });
  wrap.querySelectorAll("tbody tr.vfamr-only-row[data-vfamr-sample]").forEach((tr) => {
    const smp = tr.getAttribute("data-vfamr-sample");
    const s = _vfOnlyMap.get(smp);
    const msg = tr.querySelector(".vfamr-only-msg");
    if (msg) {
      msg.addEventListener("mouseover", (ev) => showTip(_vfamrOnlyTip(smp, s), ev));
      msg.addEventListener("mousemove", moveTip);
      msg.addEventListener("mouseout", hideTip);
    }
    const lnk = tr.querySelector(".vfamr-only-link");
    if (lnk)
      lnk.addEventListener("click", (e) => {
        e.stopPropagation();
        hideTip();
        _vfamrOnlyOpen(s);
      });
  });

  // Novelty-only indicator rows: hover → genus status summary, open novelty.
  const _novOnlyMap = new Map();
  slice.forEach((r) => {
    if (r.__noveltyOnly) _novOnlyMap.set(String(r["Specimen ID"]), r.__novelty);
  });
  wrap.querySelectorAll("tbody tr.novelty-only-row[data-novelty-sample]").forEach((tr) => {
    const smp = tr.getAttribute("data-novelty-sample");
    const s = _novOnlyMap.get(smp);
    const msg = tr.querySelector(".novelty-only-msg");
    if (msg) {
      msg.addEventListener("mouseover", (ev) => showTip(_noveltyOnlyTip(smp, s), ev));
      msg.addEventListener("mousemove", moveTip);
      msg.addEventListener("mouseout", hideTip);
    }
    const lnk = tr.querySelector(".novelty-only-link");
    if (lnk)
      lnk.addEventListener("click", (e) => {
        e.stopPropagation();
        hideTip();
        _noveltyOnlyOpen(smp);
      });
  });

  // Empty-only indicator rows have no interactive tip — nothing to wire.

  _updateSummaryPager(total, start, end, pages);
}

function _updateSummaryPager(total, start, end, pages) {
  const countEl = document.getElementById("summary-row-count");
  if (countEl) {
    if (!total) countEl.textContent = "No rows";
    else if (pages > 1) countEl.textContent = `Rows ${start + 1}–${end} of ${total.toLocaleString()}`;
    else countEl.textContent = `${total.toLocaleString()} row${total !== 1 ? "s" : ""}`;
  }
  const pager = document.getElementById("summary-pager");
  const show = pages > 1;
  if (pager) pager.style.display = show ? "flex" : "none";
  if (show) {
    const lbl = document.getElementById("summary-page-label");
    if (lbl) lbl.textContent = `Page ${_sumPage + 1} of ${pages}`;
    document.getElementById("summary-first").disabled = _sumPage === 0;
    document.getElementById("summary-prev").disabled = _sumPage === 0;
    document.getElementById("summary-next").disabled = _sumPage >= pages - 1;
    document.getElementById("summary-last").disabled = _sumPage >= pages - 1;
  }
}

/* ─── VF / AMR per-organism index (summary-table column) ────────────────
         Builds a memoized lookup of Virulence-Factor (per_gene_hits) and AMR
         (amr_genes) gene hits keyed by BOTH species and genus, so a detection
         row in the summary table can show a compact "n VF · n AMR" cell that
         hovers to the gene list and clicks through to the VF / AMR tab. Mirrors
         the aggregation used by _drawSummaryAnnotation but keyed for fast row
         lookup. Memoized the same way as _AMR_GENERA (built once after load). */
let _VFAMR_INDEX = null;
function _vfamrIndex() {
  if (_VFAMR_INDEX) return _VFAMR_INDEX;
  const byGenus = new Map();
  const bySpecies = new Map();
  const _ensure = (map, key) => {
    if (!map.has(key)) map.set(key, { vf: new Map(), amr: new Map(), samples: new Set() });
    return map.get(key);
  };
  const _lbl = (r) =>
    r["Gene"] || r["Gene Name"] || r.gene_name || r["Product"] || r.product || r["Name"] || r["Antibiotics"] || "";
  const _add = (rows, kind) => {
    (rows || []).forEach((r) => {
      const gn = String(r.Genus || r.genus || "")
        .trim()
        .toLowerCase();
      const sp = String(r.Species || r.species || "")
        .trim()
        .toLowerCase();
      const lbl = _lbl(r);
      const smp = r["Specimen ID"] || r.Sample || r.sample || null;
      const stamp = (e) => {
        const m = kind === "amr" ? e.amr : e.vf;
        if (lbl) {
          if (!m.has(lbl)) m.set(lbl, new Set());
          if (smp) m.get(lbl).add(smp);
        }
        if (smp) e.samples.add(smp);
      };
      if (gn) stamp(_ensure(byGenus, gn));
      if (sp) stamp(_ensure(bySpecies, sp));
    });
  };
  _add(PROT.per_gene_hits, "vf");
  _add(PROT.amr_genes, "amr");
  _VFAMR_INDEX = { byGenus, bySpecies };
  return _VFAMR_INDEX;
}

// Resolve the best VF/AMR annotation entry for a detection row: prefer a
// species-level match (more specific), fall back to genus-level.
function _vfamrForRow(r) {
  if (!HAS_PROT) return null;
  const idx = _vfamrIndex();
  const org = String(r["Detected Organism"] || "")
    .trim()
    .toLowerCase();
  const gn = String(r["Genus"] || "")
    .trim()
    .toLowerCase();
  // Resolve species-level (exact or "Genus species" first-two-token) and
  // genus-level annotation entries separately, so the tooltip can report
  // how specific this row's match is (species hits out of the genus total).
  let speciesInfo = null;
  if (org && idx.bySpecies.size) {
    if (idx.bySpecies.has(org)) {
      speciesInfo = idx.bySpecies.get(org);
    } else {
      const two = org.split(/\s+/).slice(0, 2).join(" ");
      if (two && idx.bySpecies.has(two)) speciesInfo = idx.bySpecies.get(two);
    }
  }
  const genusInfo = gn && idx.byGenus.has(gn) ? idx.byGenus.get(gn) : null;
  // Prefer the more specific species match; fall back to genus.
  const info = speciesInfo || genusInfo;
  if (!info) return null;
  const level = speciesInfo ? "species" : "genus";
  // Restrict the VF/AMR hits to THIS row's sample so a detection only
  // surfaces genes actually found in that sample — not every sample in
  // which the organism was detected. The index stamps each gene with the
  // set of samples it was seen in (keyed by "Specimen ID" for VF hits and
  // "Sample" for AMR hits, both equal to the detection row's Specimen ID).
  const smp = r["Specimen ID"] || r.Sample || r.sample || null;
  const _filterBySample = (m) => {
    if (!m) return new Map();
    if (smp == null) return m;
    const out = new Map();
    m.forEach((sset, gene) => {
      if (sset.has(smp)) out.set(gene, sset);
    });
    return out;
  };
  const vfMap = _filterBySample(info.vf);
  const amrMap = _filterBySample(info.amr);
  // Sample-scoped VF+AMR totals for the species and genus entries, used to
  // compute the species-specificity percentage shown atop the tooltip.
  const _total = (i) => (i ? _filterBySample(i.vf).size + _filterBySample(i.amr).size : 0);
  const speciesTotal = _total(speciesInfo);
  const genusTotal = _total(genusInfo);
  return {
    info,
    level,
    vf: vfMap.size,
    amr: amrMap.size,
    vfMap,
    amrMap,
    speciesTotal,
    genusTotal,
  };
}

// Compact cell shown in the summary table's "VF / AMR" column.
function _vfamrCellHTML(r) {
  const v = _vfamrForRow(r);
  if (!v || (!v.vf && !v.amr)) return '<span style="color:#ccc">—</span>';
  const parts = [];
  if (v.vf) parts.push(`<span class="vfamr-chip vfamr-vf">${v.vf} VF</span>`);
  if (v.amr) parts.push(`<span class="vfamr-chip vfamr-amr">${v.amr} AMR</span>`);
  const lvl = v.level === "species" ? "sp" : "gen";
  return (
    `<span class="vfamr-cell" data-vfamr="1" data-vfamr-genus="${r["Genus"] || ""}" ` +
    `style="cursor:pointer">${parts.join("")}<span class="vfamr-lvl">${lvl}</span></span>`
  );
}

function _noveltyCellHTML(r) {
  if (!HAS_NOVELTY) return '<span style="color:#ccc">—</span>';
  const sName = r["Specimen ID"] || "";
  const ns = _novSamples();
  const sampleNov = ns[sName] || {};
  const candidates = sampleNov.candidates || [];
  if (!candidates.length) return '<span style="color:#ccc">—</span>';

  // Primary: taxid match (candidates have taxid field).
  const rowTaxid = String(r["Taxonomic ID #"] || "").trim();
  const subkeyTaxid = String(r["Subkey"] || "").trim();

  // Secondary: name match (lowercased).
  const genusNameRaw = (r["Genus Name"] || r["Genus"] || "").trim();
  const speciesNameRaw = (r["Species Name"] || "").trim();
  const genusKey = genusNameRaw.toLowerCase();
  const speciesKey = speciesNameRaw.toLowerCase();

  // Check for species-level novelty candidate.
  const speciesCand = candidates.find((c) => {
    if ((c.rank || "").toLowerCase() !== "species") return false;
    if (rowTaxid && c.taxid && String(c.taxid) === rowTaxid) return true;
    if (subkeyTaxid && c.taxid && String(c.taxid) === subkeyTaxid) return true;
    return speciesKey && (c.name || "").trim().toLowerCase() === speciesKey;
  });

  // Check for genus-level novelty candidate.
  const genusCand = candidates.find((c) => {
    if ((c.rank || "").toLowerCase() !== "genus") return false;
    return genusKey && (c.name || "").trim().toLowerCase() === genusKey;
  });

  if (!speciesCand && !genusCand) return '<span style="color:#ccc">—</span>';

  const lvl = speciesCand ? "sp" : "gen";
  const cand = speciesCand || genusCand;
  const fracPct =
    cand.frac_of_highrank != null
      ? (parseFloat(cand.frac_of_highrank) * 100).toFixed(1) + "% of genus+ pool"
      : cand.frac_of_sample != null
      ? (parseFloat(cand.frac_of_sample) * 100).toFixed(3) + "% of sample"
      : "";
  const reads = cand.reads ? parseInt(cand.reads).toLocaleString() + ` ${_novUnit()}` : "";
  const tipParts = [
    `<b style="color:#f76707">Novelty signal — ${lvl === "sp" ? "species" : "genus"} level</b>`,
    `<i>${_novEsc(cand.name || genusNameRaw)}</i>`,
    reads && fracPct ? `${reads} · ${fracPct}` : reads || fracPct,
    `<span style="color:#aaa;font-size:0.85em">Classifier: ${_novEsc(_novClsShort())}</span>`,
  ]
    .filter(Boolean)
    .join("<br>");

  return (
    `<span style="display:inline-block;font-size:9px;font-weight:700;padding:0 5px;border-radius:10px;` +
    `line-height:1.5;background:#fff4e6;color:#e8590c;cursor:help" ` +
    `data-sum-tip="${String(tipParts).replace(/"/g, "&quot;")}" ` +
    `title="Novelty ${lvl === "sp" ? "species" : "genus"}-level signal detected">● ${lvl}</span>`
  );
}

// Hover tooltip for a VF/AMR summary cell: lists unique gene names with
// the number of samples each was seen in, split into VF and AMR sections.
function _vfamrTip(r) {
  const v = _vfamrForRow(r);
  if (!v) return "";
  const _section = (heading, m, color) => {
    if (!m || !m.size) return "";
    const items = Array.from(m.entries())
      .map(([n, s]) => [n, s.size])
      .sort((a, b) => b[1] - a[1]);
    let t =
      `<div style="margin-top:5px"><b style="color:${color}">${heading}</b> ` +
      `<span style="color:#aaa;font-size:0.85em">(${items.length})</span>`;
    t += `<table style="border-collapse:collapse;margin-top:2px;font-size:0.86em">`;
    items.slice(0, 10).forEach(([n, c]) => {
      t +=
        `<tr><td style="padding-right:12px;color:#e0e0e0"><i>${n}</i></td>` +
        `<td style="text-align:right;color:#90caf9;white-space:nowrap">${c} sample${c === 1 ? "" : "s"}</td></tr>`;
    });
    if (items.length > 10)
      t += `<tr><td colspan="2" style="color:#aaa;padding-top:2px">…and ${items.length - 10} more</td></tr>`;
    t += `</table></div>`;
    return t;
  };
  const lvlStyle = v.level === "species" ? "background:#f3e5f5;color:#6a1b9a" : "background:#e8f5e9;color:#2e7d32";
  let t =
    `<b>${r["Detected Organism"] || r["Genus"] || ""} — VF / AMR hits</b> ` +
    `<span style="${lvlStyle};padding:0 4px;border-radius:3px;font-size:0.78em">${v.level}-level</span>`;
  // Species-specificity line: what fraction of this genus's VF/AMR hits
  // (in this sample) are assigned at the species level. At genus level
  // there is no species match, so state that explicitly.
  if (v.level === "species" && v.genusTotal > 0) {
    const pct = Math.round((v.speciesTotal / v.genusTotal) * 100);
    t +=
      `<div style="margin-top:3px;font-size:0.86em;color:#90caf9">` +
      `<b>${pct}%</b> species match ` +
      `<span style="color:#aaa">(${v.speciesTotal} of ${v.genusTotal} genus-level hit${
        v.genusTotal === 1 ? "" : "s"
      })</span></div>`;
  } else if (v.level === "genus") {
    t += `<div style="margin-top:3px;font-size:0.86em;color:#ffb74d">No species-level hits found</div>`;
  }
  t += _section("Virulence Factors", v.vfMap || v.info.vf, "#e53935");
  t += _section("AMR Genes", v.amrMap || v.info.amr, "#fb8c00");
  if (!v.vf && !v.amr) return "";
  t += `<div style="margin-top:6px;color:#90caf9;font-size:0.82em">Click to open in the VF / AMR tab →</div>`;
  return t;
}

// Memoized set of genera (lowercased) that have AMR gene hits anywhere.
let _AMR_GENERA = null;
function _orgHasAmr(r) {
  if (!HAS_PROT) return false;
  if (!_AMR_GENERA) {
    _AMR_GENERA = new Set(
      (PROT.amr_genes || [])
        .map((a) =>
          String(a.Genus || a.genus || "")
            .trim()
            .toLowerCase(),
        )
        .filter(Boolean),
    );
  }
  const g = String(r["Genus"] || "")
    .trim()
    .toLowerCase();
  return g && _AMR_GENERA.has(g);
}

/* ─── Below-cutoff organisms with VF/AMR hits ───────────────────────────
         Returns flagged copies of DATA rows that the TASS cutoff would normally
         hide (they fail the threshold gate at every rolled-up level) BUT for
         which one or more Virulence-Factor / AMR genes were detected for the
         organism's genus IN THE SAME SAMPLE. These are surfaced as faded rows in
         the two detection tables only (Summary → Detections and the Table tab) —
         never injected into filteredData(), so KPIs, charts and the heatmap stay
         restricted to detections that pass the cutoff.

         Requires:  --annotate data present (HAS_PROT) AND the sidebar toggle
         (#filter-below-vfamr) on. Mirrors the base predicates of filteredData()
         so the faded rows honour every active sidebar filter. */
// ── VF/AMR-only samples ──────────────────────────────────────────────────
// Some samples produce NO detection rows at all (nothing passes, and nothing
// is rescued) yet still carry Virulence-Factor / AMR gene hits. Surface that
// as an indicator row so the signal isn't silently lost.
// Summarize VF/AMR gene hits per sample, independent of the detection table.
function _vfamrSampleSummary() {
  if (!HAS_PROT) return {};
  // Memoize on the (already-cached) VF/AMR index object so this O(n) scan of
  // per_gene_hits + amr_genes runs once per dataset instead of on every table
  // render — important for large --annotate reports (e.g. 25 samples / 180 MB).
  const _idx = typeof _vfamrIndex === "function" ? _vfamrIndex() : null;
  if (_idx && _idx._sampleSummary) return _idx._sampleSummary;
  const out = {};
  const lblOf = (r) =>
    String(
      r["Gene"] || r["Gene Name"] || r.gene_name || r["Product"] || r.product || r["Name"] || r["Antibiotics"] || "",
    ).trim();
  const add = (rows, kind) =>
    (rows || []).forEach((r) => {
      const smp = r["Specimen ID"] || r.Sample || r.sample;
      if (!smp) return;
      const o = (out[smp] = out[smp] || { vf: 0, amr: 0, genes: new Set(), genera: new Set() });
      if (kind === "amr") o.amr++;
      else o.vf++;
      const lbl = lblOf(r);
      if (lbl) o.genes.add(lbl);
      const gn = String(r.Genus || r.genus || "").trim();
      if (gn) o.genera.add(gn);
    });
  add(PROT.per_gene_hits, "vf");
  add(PROT.amr_genes, "amr");
  if (_idx) _idx._sampleSummary = out;
  return out;
}
// Synthetic indicator rows for non-hidden samples that have VF/AMR hits but
// contribute zero rows to the supplied detection list (fdRows).
function _vfamrOnlySampleRows(fdRows) {
  if (!HAS_PROT) return [];
  const summary = _vfamrSampleSummary();
  const ids = Object.keys(summary);
  if (!ids.length) return [];
  const present = new Set((fdRows || []).map((r) => r["Specimen ID"]));
  const txt = (document.getElementById("filter-text") || { value: "" }).value.trim();
  const ic = (document.getElementById("filter-ic") || { checked: true }).checked;
  const scope = (document.getElementById("filter-scope") || {}).value || "both";
  let rx = null;
  if (txt) {
    try {
      rx = new RegExp(txt, ic ? "i" : "");
    } catch (e) {
      rx = null;
    }
  }
  const out = [];
  ids.forEach((smp) => {
    if (sampleHidden[smp]) return;
    if (present.has(smp)) return;
    // When the user is text-searching organisms, don't surface sample-level
    // indicators (would be noise); for sample/both scope require a name match.
    if (rx && scope === "organism") return;
    if (rx && scope !== "organism" && !rx.test(smp)) return;
    out.push({ "Specimen ID": smp, "Detected Organism": "", __vfamrOnly: true, __vfamr: summary[smp] });
  });
  return out;
}
// Plain-text + HTML for the indicator cell.
function _vfamrOnlyMessageHTML(sample, s) {
  s = s || { vf: 0, amr: 0, genera: new Set() };
  const esc = (x) =>
    String(x == null ? "" : x)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  const gset = s.genera || new Set();
  const generaList = [...gset]
    .slice(0, 4)
    .map((g) => `<i>${esc(g)}</i>`)
    .join(", ");
  const generaTxt = generaList ? ` · genera: ${generaList}${gset.size > 4 ? ` +${gset.size - 4}` : ""}` : "";
  return (
    `<span class="vfamr-only-msg" style="display:inline-flex;align-items:center;gap:.5em;flex-wrap:wrap;cursor:help">` +
    `<i class="fas fa-dna" style="color:#b45309"></i>` +
    `<span><b>No passing detections for this sample</b> — but ${s.vf + s.amr} VF/AMR gene hit(s) detected ` +
    `(<span style="color:#6a1b9a;font-weight:600">${s.vf} VF</span>, ` +
    `<span style="color:#c2410c;font-weight:600">${s.amr} AMR</span>)${generaTxt}.</span>` +
    `<button type="button" class="vfamr-only-link" style="margin-left:.3em;border:1px solid #f0c36d;background:#fff;color:#b45309;border-radius:4px;padding:1px 8px;cursor:pointer;font-size:.92em;white-space:nowrap">View VF/AMR →</button>` +
    `</span>` +
    _pathEvidenceHTML(sample)
  );
}
// Hover tooltip listing the genes found for a VF/AMR-only sample.
function _vfamrOnlyTip(sample, s) {
  s = s || { genes: new Set() };
  const esc = (x) =>
    String(x == null ? "" : x)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  const genes = [...(s.genes || [])];
  let html = `<b>${esc(sample)}</b><br><span style="color:#9fb0c3">VF/AMR gene hits (${genes.length}):</span>`;
  if (genes.length) {
    html += "<br>" + genes.slice(0, 25).map(esc).join(", ");
    if (genes.length > 25) html += `<br>…and ${genes.length - 25} more`;
  }
  return html;
}
// Open the VF/AMR tab, optionally filtered to one of the sample's genera.
function _vfamrOnlyOpen(s) {
  const g = s && s.genera && s.genera.size ? [...s.genera][0] : "";
  if (typeof _jumpToProteins === "function") _jumpToProteins(g || "");
}

function _sampleCutoffById(sample) {
  const r = (DATA || []).find((x) => (x["Specimen ID"] || x.sample) === sample);
  if (r) return thresholdForRow(r);
  return parseFloat((document.getElementById("filter-min") || { value: 0 }).value) || 0;
}

function _noveltySampleSummary() {
  if (!HAS_NOVELTY) return {};
  const out = {};
  const samples = _novSamples();
  Object.keys(samples || {}).forEach((smp) => {
    const candRows = ((samples[smp] || {}).candidates || []).filter((c) => (c.rank || "").toLowerCase() === "genus");
    const mm = {};
    candRows.forEach((c) => {
      mm[(c.name || "").trim().toLowerCase()] = c;
    });
    const closed = _novClosedGenus(smp) || {};
    const keys = Array.from(new Set([...Object.keys(closed), ...Object.keys(mm)]));
    if (!keys.length) return;

    const thr = _sampleCutoffById(smp);
    let pass = 0,
      below = 0,
      mmOnly = 0;
    const genera = keys
      .map((k) => {
        const c = closed[k] || null;
        const m = mm[k] || null;
        const tass = c ? +c.tass || 0 : NaN;
        const src = c && m ? "both" : c ? "aligned_only" : "mmseqs_only";
        let status = "novelty-only";
        if (c) {
          if (!isNaN(tass) && tass >= thr) {
            status = "pass";
            pass++;
          } else {
            status = "below";
            below++;
          }
        } else {
          mmOnly++;
        }
        return {
          name: (c && c.name) || (m && m.name) || k,
          src,
          status,
          tass: c ? tass : null,
          readsAligned: c ? +c.reads || 0 : 0,
          pctAligned: c ? +c.pct || 0 : 0,
          mmHits: m ? +m.reads || 0 : 0,
          mmFracSample: m ? (+m.frac_of_sample || 0) * 100 : 0,
        };
      })
      .sort((a, b) => {
        const rank = { pass: 3, below: 2, "novelty-only": 1 };
        const d = (rank[b.status] || 0) - (rank[a.status] || 0);
        if (d) return d;
        return (b.tass || 0) - (a.tass || 0);
      });

    const sm = (samples[smp] || {}).summary || {};
    out[smp] = {
      cutoff: thr,
      total: genera.length,
      pass,
      below,
      mmOnly,
      darkPct: (+sm.dark_fraction || 0) * 100,
      genera,
    };
  });
  return out;
}

function _noveltyOnlySampleRows(fdRows) {
  if (!HAS_NOVELTY) return [];
  const summary = _noveltySampleSummary();
  const ids = Object.keys(summary);
  if (!ids.length) return [];
  // Only count real detection rows as "present"; allow synthetic indicator
  // rows (VF/AMR-only, novelty-only) to coexist for the same sample.
  const present = new Set(
    (fdRows || []).filter((r) => !r.__vfamrOnly && !r.__noveltyOnly).map((r) => r["Specimen ID"]),
  );
  const txt = (document.getElementById("filter-text") || { value: "" }).value.trim();
  const ic = (document.getElementById("filter-ic") || { checked: true }).checked;
  const scope = (document.getElementById("filter-scope") || {}).value || "both";
  let rx = null;
  if (txt) {
    try {
      rx = new RegExp(txt, ic ? "i" : "");
    } catch (e) {
      rx = null;
    }
  }
  const out = [];
  ids.forEach((smp) => {
    if (sampleHidden[smp]) return;
    if (present.has(smp)) return;
    if (rx && scope === "organism") return;
    if (rx && scope !== "organism" && !rx.test(smp)) return;
    out.push({ "Specimen ID": smp, "Detected Organism": "", __noveltyOnly: true, __novelty: summary[smp] });
  });
  return out;
}

function _noveltyOnlyMessageHTML(sample, s) {
  s = s || { total: 0, pass: 0, below: 0, mmOnly: 0, darkPct: 0, genera: [] };
  const esc = (x) =>
    String(x == null ? "" : x)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  const top = (s.genera || [])
    .slice(0, 4)
    .map((g) => `<i>${esc(g.name)}</i>`)
    .join(", ");
  const topTxt = top
    ? ` · top genera: ${top}${(s.genera || []).length > 4 ? ` +${(s.genera || []).length - 4}` : ""}`
    : "";
  return (
    `<span class="novelty-only-msg" style="display:inline-flex;align-items:center;gap:.5em;flex-wrap:wrap;cursor:help">` +
    `<i class="fas fa-lightbulb" style="color:#3949ab"></i>` +
    `<span><b>No passing detections for this sample</b> — novelty genus evidence in ${s.total} genus(es)` +
    ` (<span style="color:#2e7d32;font-weight:600">${s.pass} pass</span>, ` +
    `<span style="color:#c62828;font-weight:600">${s.below} below</span>, ` +
    `<span style="color:#ef6c00;font-weight:600">${s.mmOnly} ${_novClsLc()}-only</span>)` +
    `${topTxt} · dark matter ${(+s.darkPct || 0).toFixed(1)}%.</span>` +
    `<button type="button" class="novelty-only-link" style="margin-left:.3em;border:1px solid #9fa8da;background:#fff;color:#303f9f;border-radius:4px;padding:1px 8px;cursor:pointer;font-size:.92em;white-space:nowrap">View Novelty →</button>` +
    `</span>` +
    _pathEvidenceHTML(sample)
  );
}

function _noveltyOnlyTip(sample, s) {
  s = s || { cutoff: 0, genera: [] };
  const esc = (x) =>
    String(x == null ? "" : x)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  const rows = (s.genera || []).slice(0, 14);
  let html =
    `<b>${esc(sample)}</b><br>` +
    `<span style="color:#9fb0c3">Novelty genus status (cutoff ${(+s.cutoff || 0).toFixed(1)})</span>` +
    `<table style="border-collapse:collapse;margin-top:4px;font-size:0.84em">` +
    `<tr><th style="text-align:left;padding-right:10px;color:#90caf9">Genus</th>` +
    `<th style="text-align:right;padding-right:10px;color:#90caf9">Status</th>` +
    `<th style="text-align:right;padding-right:10px;color:#90caf9">Genus TASS</th>` +
    `<th style="text-align:right;color:#90caf9">${esc(_novClsShort())} hits</th></tr>`;
  rows.forEach((g) => {
    const stCol = g.status === "pass" ? "#2e7d32" : g.status === "below" ? "#c62828" : "#ef6c00";
    html +=
      `<tr><td style="padding-right:10px"><i>${esc(g.name)}</i></td>` +
      `<td style="text-align:right;padding-right:10px;color:${stCol}">${esc(g.status)}</td>` +
      `<td style="text-align:right;padding-right:10px">${g.tass == null ? "-" : (+g.tass).toFixed(1)}</td>` +
      `<td style="text-align:right">${g.mmHits ? Number(g.mmHits).toLocaleString() : "-"}</td></tr>`;
  });
  if ((s.genera || []).length > rows.length)
    html += `<tr><td colspan="4" style="padding-top:2px;color:#aaa">…and ${
      (s.genera || []).length - rows.length
    } more</td></tr>`;
  html += `</table>`;
  return html;
}

function _noveltyOnlyOpen(sample) {
  const btn = document.getElementById("novelty-tab-btn") || document.querySelector('.tab-btn[data-tab="novelty"]');
  if (!btn || btn.classList.contains("hidden") || btn.classList.contains("tab-disabled")) return;
  btn.click();
  _novState.sample = sample;
  _novCmpMode = "genus";
  setTimeout(() => {
    const sel = document.getElementById("novelty-sample-sel");
    if (sel) sel.value = sample;
    const gBtn = document.querySelector('#nov-cmp-toggle .nov-seg-btn[data-cmp="genus"]');
    if (gBtn) gBtn.click();
    else {
      _drawNoveltyCandidates();
      _drawNoveltyGenusCompare();
    }
  }, 0);
}

function _emptyOnlySampleRows(fdRows) {
  const allSamples = Object.keys(SAMPLE_META || {});
  if (!allSamples.length) return [];
  // Exclude samples already represented — including other indicator rows so we
  // don't double-render a sample that has VF/AMR-only or novelty-only evidence.
  const present = new Set((fdRows || []).filter((r) => !r.__emptyOnly).map((r) => r["Specimen ID"]));
  const txt = (document.getElementById("filter-text") || { value: "" }).value.trim();
  const ic = (document.getElementById("filter-ic") || { checked: true }).checked;
  const scope = (document.getElementById("filter-scope") || {}).value || "both";
  let rx = null;
  if (txt) {
    try {
      rx = new RegExp(txt, ic ? "i" : "");
    } catch (e) {
      rx = null;
    }
  }
  const out = [];
  allSamples.forEach((smp) => {
    if (sampleHidden[smp]) return;
    if (present.has(smp)) return;
    if (rx && scope === "organism") return;
    if (rx && scope !== "organism" && !rx.test(smp)) return;
    out.push({
      "Specimen ID": smp,
      "Detected Organism": "",
      __emptyOnly: true,
      __emptyMeta: SAMPLE_META[smp] || {},
    });
  });
  return out;
}

// ── Pathogen cross-reference (no-alignment evidence) ─────────────────────
// Listed pathogens can surface with ZERO reference alignment when they appear
// only in reference-free novelty candidates or in VF/AMR gene hits. These
// helpers flag those cases on the sample-level indicator rows.
const _PATH_SEV = { primary: 4, opportunistic: 3, potential: 2, commensal: 1 };
function _pathClassColor(c) {
  return (
    { primary: "#c62828", opportunistic: "#ef6c00", potential: "#1565c0", commensal: "#2e7d32" }[
      String(c || "").toLowerCase()
    ] || "#6a1b9a"
  );
}
function _pathByName(name) {
  if (!name) return null;
  const r = PATHO.by_name[String(name).trim().toLowerCase()];
  return r ? { name: r.n, class: r.c, hc: !!r.hc, status: r.s, match: "name", genus: r.g } : null;
}
function _pathByGenus(genus) {
  if (!genus) return null;
  const g = PATHO.by_genus[String(genus).trim().toLowerCase()];
  return g ? { name: String(genus), class: g.c, hc: !!g.hc, match: "genus", genusTotal: g.n } : null;
}
// Client-side pathogen matcher for a (taxid, name, rank) triple — mirrors the
// server's bin/make_report.py _match_pathogen so the report flags listed pathogens
// even on reports generated without the server pre-stamp. Order: exact taxid →
// exact name → genus rollup. Genus rollup uses the candidate name (genus rank) or
// its first token (species/strain), and is skipped for phage/host-virus/satellite
// names so they are not mis-attributed to a host genus. (Note: a virus species like
// "Bundibugyo virus" whose NCBI genus is "Orthoebolavirus" can only match by exact
// taxid or via a genus-rank "Orthoebolavirus" candidate — its name does not contain
// the genus, so first-token rollup intentionally won't reach it.)
const _PATH_ROLLUP_RANKS = new Set(["genus", "species", "subspecies", "strain", "serotype", "serovar"]);
function _matchPathTaxon(taxid, name, rank) {
  if (!HAS_PATHO) return null;
  taxid = String(taxid == null ? "" : taxid).trim();
  name = String(name == null ? "" : name).trim();
  rank = String(rank == null ? "" : rank)
    .trim()
    .toLowerCase();
  if (taxid && PATHO.by_taxid[taxid]) {
    const h = PATHO.by_taxid[taxid];
    return { match: "taxid", name: h.n, cls: h.c, status: h.s, hc: !!h.hc, genus: h.g };
  }
  const nl = name.toLowerCase();
  if (nl && PATHO.by_name[nl]) {
    const h = PATHO.by_name[nl];
    return { match: "name", name: h.n, cls: h.c, status: h.s, hc: !!h.hc, genus: h.g };
  }
  if (name && _PATH_ROLLUP_RANKS.has(rank)) {
    let gkey = "";
    if (rank === "genus") gkey = nl;
    else if (!/ phage|prophage| virus| satellite/.test(nl)) gkey = name.split(/\s+/)[0].toLowerCase();
    if (gkey && PATHO.by_genus[gkey]) {
      const g = PATHO.by_genus[gkey];
      return {
        match: "genus",
        name: gkey.charAt(0).toUpperCase() + gkey.slice(1),
        cls: g.c,
        hc: !!g.hc,
        genus: gkey,
        genusTotal: g.n,
      };
    }
  }
  return null;
}
// Novelty candidates for a sample that correspond to listed pathogens. Matched
// client-side (falls back to the server-stamped candidate.pathogen). Returns a
// list sorted high-consequence/severity/reads first.
function _sampleNoveltyPathogens(sample) {
  if (!HAS_PATHO || !HAS_NOVELTY) return [];
  const cands = (_novSamples()[sample] || {}).candidates || [];
  const out = [];
  const seen = new Set();
  cands.forEach((c) => {
    let p = _matchPathTaxon(c.taxid, c.name, c.rank);
    if (!p && c.pathogen) {
      p = { match: c.pathogen.match, name: c.pathogen.name, cls: c.pathogen.class, hc: !!c.pathogen.hc };
    }
    if (!p) return;
    const key = String(c.taxid || c.name || "");
    if (seen.has(key)) return;
    seen.add(key);
    out.push({
      name: c.name || p.name,
      taxid: c.taxid,
      rank: c.rank,
      cls: p.cls,
      hc: !!p.hc,
      match: p.match,
      genus: p.genus || "",
      genusTotal: p.genusTotal || 0,
      reads: +c.reads || 0,
    });
  });
  out.sort((a, b) => b.hc - a.hc || (_PATH_SEV[b.cls] || 0) - (_PATH_SEV[a.cls] || 0) || b.reads - a.reads);
  return out;
}
// VF/AMR gene hits in this sample whose organism is a listed pathogen. Prefers the
// server-stamped row.pathogen (matched by canonical taxid resolved from the hit's
// Source ID via the bvbrc reference — robust to mis-parsed Genus/Species text), and
// falls back to client-side species/genus name matching. By construction these rows
// only appear for samples with no passing/aligned detections, so every hit here is
// "no-alignment" evidence.
function _sampleVfamrPathogens(sample) {
  if (!HAS_PATHO || !HAS_PROT) return [];
  const byOrg = new Map();
  const scan = (rows, kind) => {
    (rows || []).forEach((r) => {
      const s = r["Specimen ID"] || r.Sample || r.sample || null;
      if (sample && s !== sample) return;
      const sp = String(r.Species || r.species || "").trim();
      const gn = String(r.Genus || r.genus || "").trim();
      let hit = null;
      let level = "species";
      let key = "s:" + sp.toLowerCase();
      // 1) server-stamped pathogen match (by canonical taxid / name / genus).
      if (r.pathogen) {
        const p = r.pathogen;
        hit = { class: p.class, hc: !!p.hc, name: p.name, genus: p.genus };
        level = p.match === "genus" ? "genus" : "species";
        key = p.taxid ? "t:" + p.taxid : p.genus ? "g:" + p.genus : "s:" + sp.toLowerCase();
      }
      // 2) fall back to client-side name/genus matching.
      if (!hit) {
        hit = _pathByName(sp);
        if (!hit && sp) {
          const two = sp.split(/\s+/).slice(0, 2).join(" ");
          hit = _pathByName(two);
        }
        if (hit) {
          hit = { class: hit.class, hc: hit.hc, name: hit.name, genus: hit.genus };
        }
      }
      if (!hit) {
        const g = _pathByGenus(gn);
        if (g) {
          hit = { class: g.class, hc: g.hc, name: g.name, genus: gn.toLowerCase() };
          level = "genus";
          key = "g:" + gn.toLowerCase();
        }
      }
      if (!hit) return;
      const label = level === "species" ? sp || hit.name : gn || hit.name;
      let e = byOrg.get(key);
      if (!e) {
        e = { org: label, level, cls: hit.class, hc: !!hit.hc, genus: hit.genus || "", vf: 0, amr: 0 };
        byOrg.set(key, e);
      }
      if (kind === "amr") e.amr += 1;
      else e.vf += 1;
    });
  };
  scan(PROT.per_gene_hits, "vf");
  scan(PROT.amr_genes, "amr");
  const out = [...byOrg.values()];
  out.sort((a, b) => b.hc - a.hc || (_PATH_SEV[b.cls] || 0) - (_PATH_SEV[a.cls] || 0));
  return out;
}
// Kraken2-classified taxa for a sample that are listed pathogens but had NO
// reference alignment (K2 Reads > 0, # Reads Aligned == 0). These rows already
// exist in DATA (carried from paths.json per organism) but may sit below the TASS
// cutoff; this surfaces them as no-alignment pathogen evidence regardless.
function _sampleKraken2Pathogens(sample) {
  if (!HAS_PATHO) return [];
  const out = [];
  const seen = new Set();
  (DATA || []).forEach((r) => {
    if ((r["Specimen ID"] || r.sample) !== sample) return;
    const k2 = +r["K2 Reads"] || 0;
    const aligned = +r["# Reads Aligned"] || 0;
    if (k2 <= 0 || aligned > 0) return; // Kraken2-classified but not aligned
    const taxid = r["Taxonomic ID #"] || r["Subkey"] || "";
    const name = r["Detected Organism"] || r["Species Name"] || "";
    const genus = r["Genus"] || r["Genus Name"] || "";
    let p = _matchPathTaxon(taxid, name, "species");
    if (!p) p = _matchPathTaxon("", genus, "genus");
    if (!p) return;
    const key = String(taxid || name);
    if (seen.has(key)) return;
    seen.add(key);
    out.push({
      name: name || p.name,
      taxid,
      cls: p.cls,
      hc: !!p.hc,
      match: p.match,
      genus: p.genus || genus.toLowerCase(),
      reads: k2,
    });
  });
  out.sort((a, b) => b.hc - a.hc || (_PATH_SEV[b.cls] || 0) - (_PATH_SEV[a.cls] || 0) || b.reads - a.reads);
  return out;
}
// Per-genus rollup: "x / y listed-pathogen taxa under genus G seen with no
// alignment". y = number of listed pathogen species the sheet has under that genus;
// x = how many distinct ones this sample surfaced (via novelty / Kraken2).
function _pathGenusRollup(items) {
  const byGenus = new Map();
  (items || []).forEach((p) => {
    const g = (p.genus || "").toLowerCase();
    if (!g) return;
    let e = byGenus.get(g);
    if (!e) {
      e = { genus: g, seen: new Set(), total: (PATHO.by_genus[g] || {}).n || 0, cls: p.cls, hc: !!p.hc };
      byGenus.set(g, e);
    }
    e.seen.add(String(p.taxid || p.name));
    if ((_PATH_SEV[p.cls] || 0) > (_PATH_SEV[e.cls] || 0)) e.cls = p.cls;
    e.hc = e.hc || !!p.hc;
  });
  return [...byGenus.values()]
    .filter((e) => e.total > 0)
    .sort((a, b) => b.hc - a.hc || (_PATH_SEV[b.cls] || 0) - (_PATH_SEV[a.cls] || 0));
}
function _pathChip(label, sub, cls, hc, esc) {
  const col = _pathClassColor(cls);
  return (
    `<span style="display:inline-flex;align-items:center;gap:.3em;border:1px solid ${col};` +
    `border-radius:4px;padding:0 6px;font-size:.9em;color:${col};background:${col}14">` +
    (hc ? `<i class="fas fa-triangle-exclamation" title="high-consequence agent"></i> ` : "") +
    `<b>${esc(label)}</b>` +
    (sub ? `<span style="color:#777;font-weight:400">${esc(sub)}</span>` : "") +
    (cls
      ? `<span style="color:${col};font-weight:700;text-transform:uppercase;font-size:.8em">${esc(cls)}</span>`
      : "") +
    `</span>`
  );
}
// Combined "listed pathogen — no alignment" evidence block for a sample, drawn from
// three reference-free / no-alignment sources: novelty LCA candidates, Kraken2
// classifications that never aligned, and VF/AMR gene hits. Adds a per-genus rollup
// (x / y listed-pathogen taxa under a genus). Returns "" when there is no evidence.
function _pathEvidenceHTML(sample) {
  if (!HAS_PATHO) return "";
  const esc = (x) =>
    String(x == null ? "" : x)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  const nov = _sampleNoveltyPathogens(sample);
  const k2 = _sampleKraken2Pathogens(sample);
  const vf = _sampleVfamrPathogens(sample);
  if (!nov.length && !k2.length && !vf.length) return "";
  const N = 6; // chips shown per source before "+N more"
  let html =
    `<div class="path-evidence" style="flex-basis:100%;display:flex;flex-wrap:wrap;gap:.4em;align-items:center;margin-top:5px">` +
    `<span style="color:#b71c1c;font-weight:700"><i class="fas fa-biohazard"></i> Listed pathogen — no alignment:</span>`;
  const _section = (label, items, render) => {
    if (!items.length) return;
    html += `<span style="color:#555;font-weight:600;font-size:.9em">${label}:</span>`;
    items.slice(0, N).forEach(render);
    if (items.length > N) html += `<span style="color:#777">+${items.length - N} more</span>`;
  };
  _section("novelty", nov, (p) => {
    const sub = `tax ${p.taxid || "?"} · ${p.rank || "?"}` + (p.match === "genus" ? " · genus rollup" : "");
    html += _pathChip(p.name, sub, p.cls, p.hc, esc);
  });
  _section("Kraken2", k2, (p) => {
    html += _pathChip(p.name, `tax ${p.taxid || "?"} · ${(+p.reads).toLocaleString()} K2 reads`, p.cls, p.hc, esc);
  });
  _section("VF/AMR", vf, (p) => {
    const parts = [];
    if (p.vf) parts.push(`${p.vf} VF`);
    if (p.amr) parts.push(`${p.amr} AMR`);
    html += _pathChip(p.org, `${parts.join(" / ")} · ${p.level}`, p.cls, p.hc, esc);
  });
  // Per-genus rollup across novelty + Kraken2 + VF/AMR evidence (VF/AMR items map to
  // their canonical sheet genus via the taxid match, so e.g. Zaire/Orthoebolavirus
  // hits roll up correctly even when the merged sheet's Genus text is mis-parsed).
  const vfForRoll = vf.map((p) => ({ taxid: p.org, name: p.org, genus: p.genus, cls: p.cls, hc: p.hc }));
  const roll = _pathGenusRollup(nov.concat(k2).concat(vfForRoll));
  if (roll.length) {
    html += `<span style="flex-basis:100%;height:0"></span>`; // line break within the flex row
    html += `<span style="color:#555;font-weight:600;font-size:.9em"><i class="fas fa-sitemap"></i> genus rollup:</span>`;
    roll.slice(0, 8).forEach((g) => {
      const col = _pathClassColor(g.cls);
      const gname = g.genus.charAt(0).toUpperCase() + g.genus.slice(1);
      html +=
        `<span title="listed pathogen taxa seen with no alignment / total listed pathogen species in this genus" ` +
        `style="display:inline-flex;align-items:center;gap:.3em;border:1px dashed ${col};border-radius:4px;padding:0 6px;font-size:.88em;color:${col}">` +
        (g.hc ? `<i class="fas fa-triangle-exclamation"></i> ` : "") +
        `<i>${esc(gname)}</i> <b>${g.seen.size}/${g.total}</b> pathogenic</span>`;
    });
  }
  html += `</div>`;
  return html;
}

function _emptyOnlyMessageHTML(sample, meta) {
  meta = meta || {};
  const esc = (x) =>
    String(x == null ? "" : x)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  const totalReads = parseInt(meta.total_reads) || 0;
  const alignedReads = parseInt(meta.aligned_reads) || 0;
  const numKeys = parseInt(meta.num_keys) || 0;
  const numSubkeys = parseInt(meta.num_subkeys) || 0;
  const numToplevel = parseInt(meta.num_toplevelkeys) || 0;
  const hadOrgs = numKeys > 0 || numSubkeys > 0 || numToplevel > 0;
  const alignPct = totalReads > 0 ? ((alignedReads / totalReads) * 100).toFixed(1) : null;

  let reason, icon, color;
  if (!hadOrgs && alignedReads === 0) {
    reason = "No Kraken2 hits · no alignment hits";
    icon = "fa-circle-xmark";
    color = "#9ca3af";
  } else if (!hadOrgs) {
    const aLabel =
      alignPct !== null
        ? `${alignedReads.toLocaleString()} reads aligned (${alignPct}%)`
        : `${alignedReads.toLocaleString()} reads aligned`;
    reason = `No organisms detected · ${aLabel}`;
    icon = "fa-circle-minus";
    color = "#6b7280";
  } else {
    const orgLabel = numKeys
      ? `${numKeys} strain${numKeys !== 1 ? "s" : ""}`
      : numSubkeys
      ? `${numSubkeys} species`
      : `${numToplevel} genus-level`;
    reason = `All detections below TASS cutoff · ${orgLabel} in report, none passing`;
    icon = "fa-circle-arrow-down";
    color = "#b45309";
  }

  const parts = [];
  if (totalReads) parts.push(`${totalReads.toLocaleString()} total reads`);
  if (hadOrgs && alignedReads && alignPct !== null)
    parts.push(`${alignedReads.toLocaleString()} aligned (${alignPct}%)`);
  const platform = meta.platform && meta.platform !== "unknown" ? esc(meta.platform) : null;
  if (platform) parts.push(platform);
  const ct = meta.control_type && meta.control_type !== "none" ? esc(meta.control_type) + " control" : null;
  if (ct) parts.push(ct);

  // Pathogen evidence with no alignment (novelty candidates / VF-AMR hits). When
  // present, suppress the "no VF/AMR · no novelty signal" tail — it would
  // contradict the evidence shown below.
  const pathEv = _pathEvidenceHTML(sample);
  return (
    `<span class="empty-only-msg" style="display:inline-flex;align-items:center;gap:.5em;flex-wrap:wrap">` +
    `<i class="fas ${icon}" style="color:${color}"></i>` +
    `<span><b>${esc(sample)}</b>: ${reason}` +
    (parts.length ? ` · <span style="color:#6b7280;font-size:.92em">${parts.join(" · ")}</span>` : "") +
    (!pathEv && (HAS_PROT || HAS_NOVELTY)
      ? ` — ${[HAS_PROT ? "no VF/AMR" : "", HAS_NOVELTY ? "no novelty signal" : ""].filter(Boolean).join(" · ")}`
      : "") +
    `</span>` +
    `</span>` +
    pathEv
  );
}

function _belowCutoffExtraRows() {
  const tgl = document.getElementById("filter-below-vfamr");
  if (!HAS_PROT || !tgl || !tgl.checked) return [];
  const onlyP = document.getElementById("filter-pass").checked;
  // "Passes threshold only" is incompatible with showing failing rows.
  if (onlyP) return [];
  const txt = document.getElementById("filter-text").value.trim();
  const ic = document.getElementById("filter-ic").checked;
  const scope = (document.getElementById("filter-scope") || {}).value || "both";
  const viewLevel = (document.getElementById("view-level") || {}).value || "Strain";
  const onlyHC = document.getElementById("filter-hc").checked;
  const mcSel = Array.from(document.getElementById("filter-mc").selectedOptions).map((o) => o.value);
  const mcFilter = mcSel.length ? new Set(mcSel) : null;
  const kingdomSel2 = Array.from(document.querySelectorAll(".fk-cb:checked")).map((cb) => cb.value);
  const kingdomFilter2 = kingdomSel2.length ? new Set(kingdomSel2) : null;
  const mtDNA = document.getElementById("filter-mt-dna") ? document.getElementById("filter-mt-dna").checked : true;
  const mtRNA = document.getElementById("filter-mt-rna") ? document.getElementById("filter-mt-rna").checked : true;
  const mtBoth = document.getElementById("filter-mt-both") ? document.getElementById("filter-mt-both").checked : true;
  let rx = null;
  if (txt) {
    try {
      rx = new RegExp(txt, ic ? "i" : "");
    } catch (e) {
      rx = null;
    }
  }
  const anyRescale = _hasAnyRescale();
  const out = [];
  for (const r of DATA) {
    if ((r["Level"] || "Strain") !== viewLevel) continue;
    if (sampleHidden[r["Specimen ID"]]) continue;
    if (rx) {
      const inS = rx.test(r["Specimen ID"] || "");
      const inO = rx.test(r["Detected Organism"] || "");
      if (scope === "sample") {
        if (!inS) continue;
      } else if (scope === "organism") {
        if (!inO) continue;
      } else if (!inS && !inO) {
        continue;
      }
    }
    if (onlyHC && !isTruthy(r["High Consequence"])) continue;
    if (mcFilter && !mcFilter.has(r["Microbial Category"] || "Unknown")) continue;
    if (kingdomFilter2) {
      const b2 = _rowKingdomBucket(r);
      if (!b2 || !kingdomFilter2.has(b2)) continue;
    }
    const mt = (r["Mol Type"] || "").toLowerCase();
    if (mt === "dna" && !mtDNA) continue;
    if (mt === "rna" && !mtRNA) continue;
    if (mt !== "dna" && mt !== "rna" && !mtBoth) continue;
    if (watchFilterMode !== "all" && watchlist.size) {
      const watched = watchlist.has(_watchKey(r));
      if (watchFilterMode === "only" && !watched) continue;
      if (watchFilterMode === "hide" && watched) continue;
    }
    // Must be EXCLUDED by the threshold gate — otherwise it's already
    // visible in filteredData() and must not be duplicated here.
    const pi = rowPassInfo(r);
    if (isNaN(pi.strain)) continue;
    const excluded = ROLLUP_PASS ? !pi.effectivePass : pi.strain < pi.thr;
    if (!excluded) continue;
    // Must carry at least one VF or AMR gene hit in its own sample.
    const v = _vfamrForRow(r);
    if (!v || (!v.vf && !v.amr)) continue;
    const inj = anyRescale ? applyRescale(Object.assign({}, r)) : Object.assign({}, r);
    inj.__belowCutoffVFAMR = true;
    out.push(inj);
  }
  return out;
}

// Small badge appended to the organism cell of a faded below-cutoff row.
function _belowCutoffBadgeHTML(r) {
  if (!r || !r.__belowCutoffVFAMR) return "";
  const pi = rowPassInfo(r);
  const sc = !isNaN(pi.strain) ? pi.strain.toFixed(1) : "—";
  const thr = pi.thr;
  return (
    `<span class="below-cutoff-badge" title="TASS ${sc} is below the cutoff (${thr}) — hidden by score, ` +
    `but shown because VF/AMR genes were detected for this organism's genus in this sample. ` +
    `Uncheck &quot;Show below-cutoff organisms with VF/AMR hits&quot; to hide.">&#x2193; below cutoff</span>`
  );
}

// Does this detection row have a matching novelty candidate (species or genus)
// in its own sample? Mirrors the matching used by _noveltyCellHTML so a row's
// novelty signal and its "surfaced because…" reason stay consistent. Returns
// {cand, lvl} ("sp" | "gen") or null.
function _rowNoveltyMatch(r) {
  if (!HAS_NOVELTY) return null;
  const candidates = (_novSamples()[r["Specimen ID"] || ""] || {}).candidates || [];
  if (!candidates.length) return null;
  const rowTaxid = String(r["Taxonomic ID #"] || "").trim();
  const subkeyTaxid = String(r["Subkey"] || "").trim();
  const genusKey = (r["Genus Name"] || r["Genus"] || "").trim().toLowerCase();
  const speciesKey = (r["Species Name"] || "").trim().toLowerCase();
  const orgKey = (r["Detected Organism"] || "").trim().toLowerCase();
  const speciesCand = candidates.find((c) => {
    if ((c.rank || "").toLowerCase() !== "species") return false;
    if (rowTaxid && c.taxid && String(c.taxid) === rowTaxid) return true;
    if (subkeyTaxid && c.taxid && String(c.taxid) === subkeyTaxid) return true;
    const cn = (c.name || "").trim().toLowerCase();
    return cn && (cn === speciesKey || cn === orgKey);
  });
  const genusCand = candidates.find((c) => {
    if ((c.rank || "").toLowerCase() !== "genus") return false;
    return genusKey && (c.name || "").trim().toLowerCase() === genusKey;
  });
  if (speciesCand) return { cand: speciesCand, lvl: "sp" };
  if (genusCand) return { cand: genusCand, lvl: "gen" };
  return null;
}

// Extra rows surfaced by the "sub-threshold & novelty-supported" toggle:
//   • below-cutoff rows that still have reference alignments  (__belowCutoffAligned)
//   • below-cutoff rows with NO alignment but a novelty genus/species match
//     in the same sample                                       (__noveltyNoAlign)
// Mirrors the sidebar-filter predicates of _belowCutoffExtraRows() so the same
// text/category/HC/mol-type filters apply. Rows already surfaced by the VF/AMR
// toggle are skipped to avoid duplicates.
function _noveltySubThresholdExtraRows() {
  const tgl = document.getElementById("filter-novelty-sub");
  if (!tgl || !tgl.checked) return [];
  const onlyP = document.getElementById("filter-pass").checked;
  if (onlyP) return []; // "passes only" is incompatible with showing failing rows
  const txt = document.getElementById("filter-text").value.trim();
  const ic = document.getElementById("filter-ic").checked;
  const scope = (document.getElementById("filter-scope") || {}).value || "both";
  const viewLevel = (document.getElementById("view-level") || {}).value || "Strain";
  const onlyHC = document.getElementById("filter-hc").checked;
  const mcSel = Array.from(document.getElementById("filter-mc").selectedOptions).map((o) => o.value);
  const mcFilter = mcSel.length ? new Set(mcSel) : null;
  const kingdomSel3 = Array.from(document.querySelectorAll(".fk-cb:checked")).map((cb) => cb.value);
  const kingdomFilter3 = kingdomSel3.length ? new Set(kingdomSel3) : null;
  const mtDNA = document.getElementById("filter-mt-dna") ? document.getElementById("filter-mt-dna").checked : true;
  const mtRNA = document.getElementById("filter-mt-rna") ? document.getElementById("filter-mt-rna").checked : true;
  const mtBoth = document.getElementById("filter-mt-both") ? document.getElementById("filter-mt-both").checked : true;
  let rx = null;
  if (txt) {
    try {
      rx = new RegExp(txt, ic ? "i" : "");
    } catch (e) {
      rx = null;
    }
  }
  // Skip rows already shown by the VF/AMR below-cutoff toggle.
  const vfTgl = document.getElementById("filter-below-vfamr");
  const vfOn = HAS_PROT && vfTgl && vfTgl.checked;
  const anyRescale = _hasAnyRescale();
  const out = [];
  for (const r of DATA) {
    if ((r["Level"] || "Strain") !== viewLevel) continue;
    if (sampleHidden[r["Specimen ID"]]) continue;
    if (rx) {
      const inS = rx.test(r["Specimen ID"] || "");
      const inO = rx.test(r["Detected Organism"] || "");
      if (scope === "sample") {
        if (!inS) continue;
      } else if (scope === "organism") {
        if (!inO) continue;
      } else if (!inS && !inO) {
        continue;
      }
    }
    if (onlyHC && !isTruthy(r["High Consequence"])) continue;
    if (mcFilter && !mcFilter.has(r["Microbial Category"] || "Unknown")) continue;
    if (kingdomFilter3) {
      const b3 = _rowKingdomBucket(r);
      if (!b3 || !kingdomFilter3.has(b3)) continue;
    }
    const mt = (r["Mol Type"] || "").toLowerCase();
    if (mt === "dna" && !mtDNA) continue;
    if (mt === "rna" && !mtRNA) continue;
    if (mt !== "dna" && mt !== "rna" && !mtBoth) continue;
    if (watchFilterMode !== "all" && watchlist.size) {
      const watched = watchlist.has(_watchKey(r));
      if (watchFilterMode === "only" && !watched) continue;
      if (watchFilterMode === "hide" && watched) continue;
    }
    // Must be EXCLUDED by the threshold gate — otherwise it's already visible.
    const pi = rowPassInfo(r);
    if (isNaN(pi.strain)) continue;
    const excluded = ROLLUP_PASS ? !pi.effectivePass : pi.strain < pi.thr;
    if (!excluded) continue;
    const aligned = (+r["# Reads Aligned"] || 0) > 0;
    const nov = _rowNoveltyMatch(r);
    // Surface only rows that carry signal: aligned (below cutoff) or novelty-backed.
    if (!aligned && !nov) continue;
    // Don't duplicate rows already surfaced by the VF/AMR toggle.
    if (vfOn) {
      const v = _vfamrForRow(r);
      if (v && (v.vf || v.amr)) continue;
    }
    const inj = anyRescale ? applyRescale(Object.assign({}, r)) : Object.assign({}, r);
    if (aligned) inj.__belowCutoffAligned = true;
    else inj.__noveltyNoAlign = true;
    if (nov) {
      inj.__noveltyMatch = true;
      inj.__noveltyMatchLvl = nov.lvl;
    }
    out.push(inj);
  }
  return out;
}

// Badge(s) for a sub-threshold / novelty-supported row, appended to the
// organism cell so the reason it's visible is obvious.
function _subThresholdBadgeHTML(r) {
  if (!r || (!r.__belowCutoffAligned && !r.__noveltyNoAlign)) return "";
  const pi = rowPassInfo(r);
  const sc = !isNaN(pi.strain) ? pi.strain.toFixed(1) : "—";
  let html = "";
  if (r.__belowCutoffAligned) {
    const reads = _fmtInt(r["# Reads Aligned"]);
    html +=
      `<span class="subthr-aligned-badge" title="TASS ${sc} is below the cutoff (${pi.thr}) — hidden by score, ` +
      `but ${reads} reference read(s) aligned in this sample.">&#x2193; sub-threshold (aligned)</span>`;
  }
  if (r.__noveltyNoAlign) {
    html +=
      `<span class="novelty-noalign-badge" title="No reference alignment, but the novelty classifier placed this ` +
      `organism's ${r.__noveltyMatchLvl === "gen" ? "genus" : "species"} in this sample (TASS ${sc}).">` +
      `&#10022; novelty, no alignment</span>`;
  }
  return html;
}

// Wire summary detection-table controls (pager + group + page-size).
(function _initSummaryTableControls() {
  function _go(n) {
    _sumPage = n;
    _renderSummaryTable(filteredData());
  }
  const ids = {
    first: () => _go(0),
    prev: () => _go(_sumPage - 1),
    next: () => _go(_sumPage + 1),
    last: () => {
      const total = filteredData().length;
      const ps = _sumPageSize();
      _go(ps > 0 ? Math.max(0, Math.ceil(total / ps) - 1) : 0);
    },
  };
  Object.entries(ids).forEach(([k, fn]) => {
    const el = document.getElementById("summary-" + k);
    if (el) el.addEventListener("click", fn);
  });
  const ps = document.getElementById("summary-page-size");
  if (ps)
    ps.addEventListener("change", () => {
      _sumPage = 0;
      _renderSummaryTable(filteredData());
    });
  const grp = document.getElementById("summary-group-sample");
  if (grp)
    grp.addEventListener("change", () => {
      // grouping reads best when sorted by sample; force that when enabled
      if (grp.checked) {
        _sumSortCol = "Specimen ID";
        _sumSortAsc = true;
      }
      _sumPage = 0;
      _renderSummaryTable(filteredData());
    });
})();

/* ── Summary tab pin-bar ─────────────────────────────────────────────── */
function _updateSumPinBar() {
  const bar = document.getElementById("sum-pin-bar");
  const countEl = document.getElementById("sum-pin-bar-count");
  if (!bar) return;
  const n = _sumPinned.size;
  bar.style.display = n > 0 ? "flex" : "none";
  if (countEl) countEl.textContent = `${n} row${n === 1 ? "" : "s"} pinned`;
}

(function _initSumPinBar() {
  const compareBtn = document.getElementById("sum-compare-btn");
  if (compareBtn)
    compareBtn.addEventListener("click", () => {
      // Gather rows matching sumPinned keys from full DATA
      const pinnedRows = DATA.filter((r) => _sumPinned.has(_sumRowKey(r)));
      if (!pinnedRows.length) return;
      _openCompareModalForRows(pinnedRows);
    });
  const unpinBtn = document.getElementById("sum-unpin-all-btn");
  if (unpinBtn)
    unpinBtn.addEventListener("click", () => {
      _sumPinned.clear();
      _updateSumPinBar();
      _renderSummaryTable(filteredData());
    });
})();

// Persistent per-sample hide state for the genera comparison legend.
// Survives redraws triggered by metric/topN changes.
let _genusLegendHidden = {};

function _drawSummaryGenera(fd, samples) {
  const wrap = document.getElementById("summary-genus-wrap");
  if (!wrap) return;
  // If the pane hasn't laid out yet getBoundingClientRect returns 0 — defer
  // one animation frame so the browser can measure the real container width.
  if (wrap.getBoundingClientRect().width === 0) {
    requestAnimationFrame(() => _drawSummaryGenera(fd, samples));
    return;
  }
  wrap.innerHTML = "";
  const metric = (document.getElementById("summary-genus-metric") || {}).value || "# Reads Aligned";
  const topN = parseInt((document.getElementById("summary-genus-topn") || {}).value || "15", 10);
  const showValues = !!(document.getElementById("summary-genus-values") || {}).checked;

  const rows = fd.filter((r) => r["Genus"]);
  if (!rows.length) {
    wrap.innerHTML = '<p style="color:#888;padding:1em">No genus-level data for the current filters.</p>';
    return;
  }

  // genus -> {total, perSample: {sample: val}}
  const genera = {};
  rows.forEach((r) => {
    const g = r["Genus"];
    const sn = r["Specimen ID"] || "Unknown";
    let v;
    if (metric === "count") v = 1;
    else if (metric === "TASS Score") v = num(r["TASS Score"]);
    else v = num(r["# Reads Aligned"]);
    if (!genera[g]) genera[g] = { total: 0, perSample: {} };
    if (metric === "TASS Score") {
      genera[g].perSample[sn] = Math.max(genera[g].perSample[sn] || 0, v);
      genera[g].total = Math.max(genera[g].total, v);
    } else {
      genera[g].perSample[sn] = (genera[g].perSample[sn] || 0) + v;
      genera[g].total += v;
    }
  });

  const sampleList = _orderedSamples(samples.filter((sn) => sn));
  const color = d3.scaleOrdinal().domain(sampleList).range(PALETTE);

  // Only count visible samples toward totals for bar scaling
  const visibleSamples = sampleList.filter((sn) => !_genusLegendHidden[sn]);

  // Recompute totals using only visible samples for ordering + x-scale
  const visGenera = {};
  rows.forEach((r) => {
    const g = r["Genus"];
    const sn = r["Specimen ID"] || "Unknown";
    if (_genusLegendHidden[sn]) return;
    let v;
    if (metric === "count") v = 1;
    else if (metric === "TASS Score") v = num(r["TASS Score"]);
    else v = num(r["# Reads Aligned"]);
    if (!visGenera[g]) visGenera[g] = 0;
    if (metric === "TASS Score") visGenera[g] = Math.max(visGenera[g], v);
    else visGenera[g] += v;
  });

  const ordered = Object.entries(genera)
    .sort((a, b) => {
      const va = visGenera[a[0]] || 0;
      const vb = visGenera[b[0]] || 0;
      return vb - va;
    })
    .slice(0, topN);

  // Legend dimensions
  const legendRowH = 17;
  const legendCtrlH = 46; // header (10) + hide/show row (24) + gap + divider
  const legendW = 140;
  const legendPad = 12;

  // getBoundingClientRect gives the real rendered width even before clientWidth settles.
  const W = Math.max(
    600,
    wrap.getBoundingClientRect().width ||
      (wrap.parentElement && wrap.parentElement.getBoundingClientRect().width) ||
      wrap.clientWidth ||
      900,
  );
  const rowH = 26,
    mT = 10,
    mB = 10,
    mL = 150,
    mR = legendW + legendPad + 50; // +50 leaves ~42 px for value labels before the legend
  const H = Math.max(ordered.length * rowH + mT + mB, legendCtrlH + sampleList.length * legendRowH + mT + 10);
  const iW = Math.max(120, W - mL - mR);
  const maxVisTotal = d3.max(ordered, (d) => visGenera[d[0]] || 0) || 1;
  const x = d3.scaleLinear().domain([0, maxVisTotal]).range([0, iW]);

  const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  const g = svg.append("g").attr("transform", `translate(${mL},${mT})`);

  // ── Bars ──────────────────────────────────────────────────────────────
  ordered.forEach((entry, i) => {
    const [genus, info] = entry;
    const y = i * rowH;

    // Build a per-genus hover tip (total + per-sample breakdown) so the value
    // detail lives on hover instead of as an on-chart label.
    const _metricLbl = metric === "count" ? "# organisms" : metric;
    const _totLbl = metric === "TASS Score" ? "max" : "total";
    const _fmtV = (val) =>
      metric === "count" ? val + " organism(s)" : metric === "TASS Score" ? (+val).toFixed(1) : _fmtInt(val);
    const _visTot = visGenera[genus] || 0;
    const _bd = visibleSamples
      .map((sn) => [sn, info.perSample[sn] || 0])
      .filter((s) => s[1] > 0)
      .sort((a, b) => b[1] - a[1]);
    const _genusTip =
      `<b>${genus}</b><br><span style="color:#9bb">${_totLbl} ${_metricLbl}:</span> <b>${_fmtV(_visTot)}</b>` +
      (_bd.length
        ? `<br><br><b>Per sample:</b><br>` + _bd.map(([sn, val]) => `• ${sn}: <b>${_fmtV(val)}</b>`).join("<br>")
        : "");

    // genus label (hover shows the full breakdown)
    svg
      .append("text")
      .attr("x", mL - 8)
      .attr("y", mT + y + rowH / 2)
      .attr("text-anchor", "end")
      .attr("dominant-baseline", "middle")
      .style("font-size", "11px")
      .style("font-style", "italic")
      .style("fill", "#333")
      .style("cursor", "pointer")
      .text(genus.length > 22 ? genus.slice(0, 21) + "…" : genus)
      .on("mouseover", (ev) => showTip(_genusTip, ev))
      .on("mousemove", moveTip)
      .on("mouseout", hideTip);

    let xacc = 0;
    const segs = visibleSamples.map((sn) => [sn, info.perSample[sn] || 0]).filter((s) => s[1] > 0);
    // Short label for a segment/total value (compact so it fits on the bar).
    const _segLbl = (val) =>
      metric === "count" ? String(val) : metric === "TASS Score" ? (+val).toFixed(1) : _fmtBig(val).short;
    segs.forEach(([sn, val]) => {
      const segX = x(xacc);
      const segW = Math.max(0, x(xacc + val) - x(xacc));
      g.append("rect")
        .attr("x", segX)
        .attr("y", y + 3)
        .attr("width", segW)
        .attr("height", rowH - 6)
        .attr("fill", color(sn))
        .attr("stroke", "#fff")
        .on("mouseover", (ev) =>
          showTip(`<b>${genus}</b><br>${sn}<br>${metric === "count" ? val + " organism(s)" : _fmtInt(val)}`, ev),
        )
        .on("mousemove", moveTip)
        .on("mouseout", hideTip);
      // Per-segment value — only drawn when "Show values" is on AND the segment
      // is wide enough to hold the text, so labels never overlap each other.
      // Narrow segments stay legible via hover. ~6 px per char + padding.
      if (showValues) {
        const _txt = _segLbl(val);
        if (segW >= _txt.length * 6 + 6) {
          g.append("text")
            .attr("x", segX + segW / 2)
            .attr("y", y + rowH / 2)
            .attr("text-anchor", "middle")
            .attr("dominant-baseline", "middle")
            .style("font-size", "9px")
            .style("font-weight", "600")
            .style("fill", "#fff")
            .style("pointer-events", "none")
            .text(_txt);
        }
      }
    });

    // Per-genus total at the bar end — one per row, so these never overlap
    // each other. Only shown when "Show values" is on. Clamped left of the
    // legend so the number can't collide with it.
    if (showValues) {
      const _visTotG = visGenera[genus] || 0;
      if (_visTotG > 0) {
        const _tx = Math.min(x(_visTotG) + 6, iW + 44);
        g.append("text")
          .attr("x", _tx)
          .attr("y", y + rowH / 2)
          .attr("dominant-baseline", "middle")
          .style("font-size", "10px")
          .style("font-weight", "700")
          .style("fill", "#455a64")
          .style("pointer-events", "none")
          .text(_segLbl(_visTotG));
      }
    }
  });

  // ── Vertical right-aligned legend ─────────────────────────────────────
  const legendX = W - legendW - 4;
  const legendG = svg.append("g").attr("transform", `translate(${legendX},${mT})`);

  // "Samples" header
  legendG
    .append("text")
    .attr("x", 0)
    .attr("y", 11)
    .style("font-size", "10px")
    .style("font-weight", "600")
    .style("fill", "#444")
    .text("Samples");

  // Hide All / Show All links (y=26 gives clear gap below header)
  const allHidden = sampleList.every((sn) => _genusLegendHidden[sn]);
  const allVisible = sampleList.every((sn) => !_genusLegendHidden[sn]);

  legendG
    .append("text")
    .attr("x", 0)
    .attr("y", 26)
    .style("font-size", "9px")
    .style("fill", allHidden ? "#aaa" : "#1c7ed6")
    .style("cursor", allHidden ? "default" : "pointer")
    .style("text-decoration", allHidden ? "none" : "underline")
    .text("Hide all")
    .on("click", () => {
      if (allHidden) return;
      sampleList.forEach((sn) => {
        _genusLegendHidden[sn] = true;
      });
      _drawSummaryGenera(fd, samples);
    });

  legendG
    .append("text")
    .attr("x", 54)
    .attr("y", 26)
    .style("font-size", "9px")
    .style("fill", allVisible ? "#aaa" : "#1c7ed6")
    .style("cursor", allVisible ? "default" : "pointer")
    .style("text-decoration", allVisible ? "none" : "underline")
    .text("Show all")
    .on("click", () => {
      if (allVisible) return;
      _genusLegendHidden = {};
      _drawSummaryGenera(fd, samples);
    });

  // Divider sits below the hide/show row with a small gap
  legendG
    .append("line")
    .attr("x1", 0)
    .attr("x2", legendW - 8)
    .attr("y1", 34)
    .attr("y2", 34)
    .attr("stroke", "#ddd")
    .attr("stroke-width", 1);

  // Per-sample items — each row is legendRowH tall, swatch+text vertically centred
  sampleList.forEach((sn, si) => {
    const isHidden = !!_genusLegendHidden[sn];
    // Items start at legendCtrlH (46) so first item top = 46, well below divider at 34
    const iy = legendCtrlH + si * legendRowH;

    const itemG = legendG
      .append("g")
      .attr("transform", `translate(0,${iy})`)
      .style("cursor", "pointer")
      .on("click", () => {
        _genusLegendHidden[sn] = !_genusLegendHidden[sn];
        _drawSummaryGenera(fd, samples);
      })
      .on("mouseover", function () {
        d3.select(this).select("rect").attr("stroke", "#555");
      })
      .on("mouseout", function () {
        d3.select(this).select("rect").attr("stroke", "none");
      });

    // Color swatch — top-aligned at y=0, height=11
    itemG
      .append("rect")
      .attr("width", 11)
      .attr("height", 11)
      .attr("y", 0)
      .attr("rx", 2)
      .attr("fill", color(sn))
      .attr("opacity", isHidden ? 0.25 : 1)
      .attr("stroke", "none");

    // Sample name — baseline vertically centred with the swatch
    itemG
      .append("text")
      .attr("x", 15)
      .attr("y", 9)
      .style("font-size", "9px")
      .style("fill", isHidden ? "#bbb" : "#444")
      .style("text-decoration", isHidden ? "line-through" : "none")
      .text(sn.length > 17 ? sn.slice(0, 16) + "…" : sn);
  });
}
