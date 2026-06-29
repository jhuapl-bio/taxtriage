/* ═══════════════════════════════════════════════════════════════════════════
       -  §  CROSS-ENTRY COMPARISON  (organism-profile similarity / enrichment)
       -     Compares the organism-hit profiles of either individual samples or
       -     aggregated metadata groups (country, host disease, …) within a run.
       -     Analyses: cosine similarity, Jaccard overlap, enrichment (log2 fold
       -     vs the rest), and a raw profile heatmap. Profiles are built at genus
       -     or species level from presence/abundance, TASS, or coverage.
       -     Functions: _cmpAvailableFields, _cmpBuildProfiles, _cmpCosine,
       -     _cmpJaccard, _cmpEnrichment, _cmpRenderSquare, _cmpRenderRect,
       -     _buildComparison.
═══════════════════════════════════════════════════════════════════════════ */
let _cmpAnalysis = "profile";
let _cmpWired = false;

// Candidate metadata fields (ordered) offered in the field selector.
// NOTE: exposed via a *hoisted function declaration* (not a `const`) so it
// is reachable even when _cmpAvailableFields() runs during init() — before
// execution reaches this point in the script. A plain `const` here threw a
// temporal-dead-zone ReferenceError ("Cannot access '_CMP_FIELD_ORDER'
// before initialization") that aborted init() and left the report blank.
function _cmpFieldOrder() {
  return [
    "sample_origin_country",
    "sample_origin_state_province_territory",
    "host_disease",
    "host_scientific_name",
    "environmental_site",
    "sequencing_instrument",
    "sequencing_platform",
    "library_preparation_kit",
    "submitter_organization_name",
    "organism",
    "run_id",
    "location",
  ];
}

// Metadata fields that actually have ≥1 non-empty value across RUN_META.
// Array values (e.g. host_disease) are considered populated if any element is non-empty.
function _cmpAvailableFields() {
  return _cmpFieldOrder().filter((f) =>
    (RUN_META || []).some((r) => {
      const _v = r[f];
      if (_v == null) return false;
      if (Array.isArray(_v)) return _v.some((s) => String(s).trim() !== "");
      return String(_v).trim() !== "";
    }),
  );
}

// Build per-entry organism profiles.
// Returns { entries:[{id,label,group}], taxa:[...], profile:{id:{taxon:value}}, unit, value }
function _cmpBuildProfiles() {
  const unit = (document.getElementById("cmp-unit") || {}).value || "sample";
  const field = (document.getElementById("cmp-field") || {}).value || "";
  const level = (document.getElementById("cmp-level") || {}).value || "genus";
  const value = (document.getElementById("cmp-value") || {}).value || "presence";

  // sample_name → metadata field value
  const fieldBy = {};
  (RUN_META || []).forEach((r) => {
    const v = r[field];
    if (v != null && String(v).trim() !== "") fieldBy[r.sample_name] = String(v).trim();
  });

  const taxonOf = (r) =>
    level === "genus"
      ? (r["Genus Name"] || r["Genus"] || "").trim()
      : (r["Species Name"] || r["Detected Organism"] || "").trim();

  // The contribution of a single detection row, given the value basis.
  const valOf = (r) => {
    if (value === "tass") {
      const t = level === "genus" ? parseFloat(r["Genus TASS"]) : parseFloat(r["Species TASS"]);
      const fallback = parseFloat(r["TASS Score"]);
      return Math.max(isNaN(t) ? 0 : t, isNaN(fallback) ? 0 : fallback);
    }
    if (value === "coverage") {
      const b = parseFloat(r["Breadth %"]);
      const c = parseFloat(r["Coverage"]);
      return Math.max(isNaN(b) ? 0 : b, isNaN(c) ? 0 : c);
    }
    // presence / abundance → read count (presence implied by >0)
    const reads = parseFloat(r["# Reads Aligned"]);
    return isNaN(reads) ? 1 : Math.max(reads, 1);
  };

  // How to combine multiple rows of the same taxon within an entry.
  const combine = value === "presence" ? (a, b) => a + b : (a, b) => Math.max(a, b);

  const profile = {}; // entryId → {taxon → value}
  const entryGroup = {}; // entryId → metadata group value (for sample mode coloring)
  const taxaSet = new Set();

  filteredData().forEach((r) => {
    const sample = r["Specimen ID"] || "";
    if (!sample || sampleHidden[sample]) return;
    const taxon = taxonOf(r);
    if (!taxon) return;
    let entryId;
    if (unit === "group") {
      entryId = fieldBy[sample];
      if (entryId == null) return; // sample has no value for this field → excluded
    } else {
      entryId = sample;
      if (fieldBy[sample] != null) entryGroup[entryId] = fieldBy[sample];
    }
    const v = valOf(r);
    if (!profile[entryId]) profile[entryId] = {};
    profile[entryId][taxon] = profile[entryId][taxon] == null ? v : combine(profile[entryId][taxon], v);
    taxaSet.add(taxon);
  });

  const entries = Object.keys(profile)
    .sort((a, b) => a.localeCompare(b))
    .map((id) => ({ id, label: id, group: unit === "group" ? id : entryGroup[id] || null }));

  return {
    unit,
    value,
    field,
    level,
    entries,
    taxa: [...taxaSet],
    profile,
  };
}

// Cosine similarity between two {taxon:value} maps.
function _cmpCosine(a, b) {
  let dot = 0,
    na = 0,
    nb = 0;
  for (const k in a) {
    na += a[k] * a[k];
    if (b[k] != null) dot += a[k] * b[k];
  }
  for (const k in b) nb += b[k] * b[k];
  if (na === 0 || nb === 0) return 0;
  return dot / (Math.sqrt(na) * Math.sqrt(nb));
}

// Jaccard overlap of the presence sets of two {taxon:value} maps.
function _cmpJaccard(a, b) {
  const sa = new Set(Object.keys(a).filter((k) => a[k] > 0));
  const sb = new Set(Object.keys(b).filter((k) => b[k] > 0));
  if (!sa.size && !sb.size) return 0;
  let inter = 0;
  sa.forEach((k) => {
    if (sb.has(k)) inter++;
  });
  const uni = sa.size + sb.size - inter;
  return uni === 0 ? 0 : inter / uni;
}

// Enrichment: per entry × taxon log2 fold-change of mean value vs all OTHER entries.
// Returns { taxa:[topTaxa], matrix:{entryId:{taxon:log2fc}} }
function _cmpEnrichment(pf) {
  const { entries, taxa, profile } = pf;
  const eps = 1e-6;
  const meanOver = (ids, taxon) => {
    if (!ids.length) return 0;
    let s = 0;
    ids.forEach((id) => (s += (profile[id] && profile[id][taxon]) || 0));
    return s / ids.length;
  };
  const allIds = entries.map((e) => e.id);
  const matrix = {};
  const score = {}; // taxon → max |log2fc| (for ranking)
  entries.forEach((e) => {
    matrix[e.id] = {};
    const others = allIds.filter((id) => id !== e.id);
    taxa.forEach((t) => {
      const mIn = (profile[e.id] && profile[e.id][t]) || 0;
      const mOut = meanOver(others, t);
      const fc = Math.log2((mIn + eps) / (mOut + eps));
      matrix[e.id][t] = fc;
      score[t] = Math.max(score[t] || 0, Math.abs(fc));
    });
  });
  const topTaxa = [...taxa].sort((a, b) => (score[b] || 0) - (score[a] || 0)).slice(0, 25);
  return { taxa: topTaxa, matrix };
}

// ── Square (entry × entry) similarity heatmap ──────────────────────────
function _cmpRenderSquare(wrap, entries, simFn, fmt) {
  wrap.innerHTML = "";
  const ids = entries.map((e) => e.id);
  const n = ids.length;
  const cell = Math.max(16, Math.min(40, Math.floor(560 / n)));
  const labelLen = entries.reduce((m, e) => Math.max(m, String(e.label).length), 0);
  const marginL = Math.min(280, Math.max(90, labelLen * 7.5));
  const marginT = Math.min(220, Math.max(70, labelLen * 5.6));
  const W = marginL + n * cell + 20; // no right-side legend
  const H = marginT + n * cell + 76; // extra bottom space for horizontal legend
  const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H).style("overflow", "visible");
  const g = svg.append("g").attr("transform", `translate(${marginL},${marginT})`);
  const color = d3.scaleSequential(d3.interpolateBlues).domain([0, 1]);

  // Precompute matrix
  const M = ids.map((a) =>
    ids.map((b) =>
      simFn(
        entries.find((e) => e.id === a),
        entries.find((e) => e.id === b),
      ),
    ),
  );

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const v = M[i][j];
      g.append("rect")
        .attr("x", j * cell)
        .attr("y", i * cell)
        .attr("width", cell - 1)
        .attr("height", cell - 1)
        .attr("rx", 2)
        .attr("fill", color(v))
        .on("mouseover", (ev) =>
          showTip(`<b>${entries[i].label}</b> vs <b>${entries[j].label}</b><br>Similarity: ${fmt(v)}`, ev),
        )
        .on("mousemove", moveTip)
        .on("mouseout", hideTip);
      if (cell >= 28) {
        g.append("text")
          .attr("x", j * cell + (cell - 1) / 2)
          .attr("y", i * cell + (cell - 1) / 2)
          .attr("dy", "0.35em")
          .attr("text-anchor", "middle")
          .style("font-size", "9px")
          .attr("fill", v > 0.6 ? "#fff" : "#333")
          .text(fmt(v));
      }
    }
  }
  // Row labels
  entries.forEach((e, i) => {
    g.append("text")
      .attr("x", -6)
      .attr("y", i * cell + (cell - 1) / 2)
      .attr("dy", "0.35em")
      .attr("text-anchor", "end")
      .style("font-size", "10px")
      .attr("fill", "#333")
      .text(e.label.length > 28 ? e.label.slice(0, 27) + "…" : e.label);
  });
  // Column labels (rotated)
  entries.forEach((e, j) => {
    g.append("text")
      .attr("transform", `translate(${j * cell + (cell - 1) / 2}, -6) rotate(-45)`)
      .attr("text-anchor", "start")
      .style("font-size", "10px")
      .attr("fill", "#333")
      .text(e.label.length > 24 ? e.label.slice(0, 23) + "…" : e.label);
  });
  // Legend — horizontal bar below the matrix
  _cmpColorLegend(svg, marginL, marginT + n * cell + 22, color, 0, 1, fmt);
}

// ── Rectangular (entry × taxon) heatmap (profile or enrichment) ─────────
function _cmpRenderRect(wrap, entries, taxa, getVal, colorFn, fmt, domain) {
  wrap.innerHTML = "";
  const nR = entries.length,
    nC = taxa.length;
  if (!nC) return;
  const cellW = Math.max(16, Math.min(34, Math.floor(620 / nC)));
  const cellH = Math.max(16, Math.min(34, Math.floor(420 / Math.max(nR, 1))));
  const taxLen = taxa.reduce((m, t) => Math.max(m, String(t).length), 0);
  const labelLen = entries.reduce((m, e) => Math.max(m, String(e.label).length), 0);
  const marginL = Math.min(280, Math.max(100, labelLen * 7.5));
  const marginT = Math.min(240, Math.max(80, taxLen * 5.8));
  const W = marginL + nC * cellW + 20; // no right-side legend; minimal right padding
  const H = marginT + nR * cellH + 76; // extra bottom room for horizontal legend
  const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H).style("overflow", "visible");
  const g = svg.append("g").attr("transform", `translate(${marginL},${marginT})`);

  entries.forEach((e, i) => {
    taxa.forEach((t, j) => {
      const v = getVal(e.id, t);
      g.append("rect")
        .attr("x", j * cellW)
        .attr("y", i * cellH)
        .attr("width", cellW - 1)
        .attr("height", cellH - 1)
        .attr("rx", 2)
        .attr("fill", colorFn(v))
        .on("mouseover", (ev) => showTip(`<b>${e.label}</b><br>${t}: ${fmt(v)}`, ev))
        .on("mousemove", moveTip)
        .on("mouseout", hideTip);
    });
    g.append("text")
      .attr("x", -6)
      .attr("y", i * cellH + (cellH - 1) / 2)
      .attr("dy", "0.35em")
      .attr("text-anchor", "end")
      .style("font-size", "10px")
      .attr("fill", "#333")
      .text(e.label.length > 28 ? e.label.slice(0, 27) + "…" : e.label);
  });
  taxa.forEach((t, j) => {
    g.append("text")
      .attr("transform", `translate(${j * cellW + (cellW - 1) / 2}, -6) rotate(-45)`)
      .attr("text-anchor", "start")
      .style("font-size", "9.5px")
      .attr("fill", "#333")
      .text(t.length > 26 ? t.slice(0, 25) + "…" : t);
  });
  // Place the color legend horizontally below the heatmap cells
  _cmpColorLegend(svg, marginL, marginT + nR * cellH + 22, colorFn, domain[0], domain[1], fmt);
}

// Horizontal color-scale legend (below the heatmap cells).
function _cmpColorLegend(svg, x, y, colorFn, lo, hi, fmt) {
  const w = 140,
    h = 10,
    steps = 40;
  for (let i = 0; i < steps; i++) {
    const t = i / (steps - 1);
    const v = lo + (hi - lo) * t;
    svg
      .append("rect")
      .attr("x", x + t * w)
      .attr("y", y)
      .attr("width", w / steps + 0.5)
      .attr("height", h)
      .attr("fill", colorFn(v));
  }
  svg
    .append("text")
    .attr("x", x - 3)
    .attr("y", y + 8)
    .style("font-size", "9px")
    .attr("fill", "#555")
    .attr("text-anchor", "end")
    .text(fmt(lo));
  svg
    .append("text")
    .attr("x", x + w + 3)
    .attr("y", y + 8)
    .style("font-size", "9px")
    .attr("fill", "#555")
    .attr("text-anchor", "start")
    .text(fmt(hi));
}

function _buildComparison() {
  const section = document.getElementById("cmp-section");
  if (!section) return;
  const fields = _cmpAvailableFields();

  const unitSel = document.getElementById("cmp-unit");
  const fieldSel = document.getElementById("cmp-field");
  const levelSel = document.getElementById("cmp-level");
  const valueSel = document.getElementById("cmp-value");
  const fieldLabel = document.getElementById("cmp-field-label");
  const desc = document.getElementById("cmp-desc");
  const noData = document.getElementById("cmp-no-data");
  const chart = document.getElementById("cmp-chart");

  // Populate field selector. Show EVERY known metadata field, but grey out
  // (disable) the ones with no values in the current run and explain why on
  // hover — rather than hiding them — so users can see what's missing.
  if (fieldSel) {
    const prev = fieldSel.value;
    const avail = new Set(fields);
    fieldSel.innerHTML = _cmpFieldOrder()
      .map((f) => {
        const ok = avail.has(f);
        const attrs = ok
          ? ""
          : ` disabled title="No values for “${_metaKeyLabel(
              f,
            )}” in this run — add this column in the Run Metadata tab to enable it." style="color:#aaa"`;
        return `<option value="${f}"${attrs}>${_metaKeyLabel(f)}${ok ? "" : " — no data"}</option>`;
      })
      .join("");
    // Selection priority: keep a still-valid previous choice; otherwise
    // default the group/label field to State / Province / Territory; else
    // fall back to the first field that actually has data.
    let target = null;
    if (prev && avail.has(prev)) target = prev;
    else if (avail.has("sample_origin_state_province_territory")) target = "sample_origin_state_province_territory";
    else target = fields[0] || null;
    if (target) fieldSel.value = target;
  }

  // Group comparison needs at least one populated metadata field — disable
  // the "Metadata group" unit option (greyed, with a reason) when none exist.
  if (unitSel) {
    const grpOpt = unitSel.querySelector('option[value="group"]');
    if (grpOpt) {
      const anyField = fields.length > 0;
      grpOpt.disabled = !anyField;
      grpOpt.title = anyField
        ? ""
        : "No metadata fields with values — group comparison needs at least one populated metadata column.";
      grpOpt.style.color = anyField ? "" : "#aaa";
      if (!anyField && unitSel.value === "group") unitSel.value = "sample";
    }
  }

  // Wire controls once
  if (!_cmpWired) {
    _cmpWired = true;
    [unitSel, fieldSel, levelSel, valueSel].forEach((el) => {
      if (el) el.addEventListener("change", _buildComparison);
    });
    document.querySelectorAll("#cmp-analysis-tabs .cmp-atab").forEach((btn) => {
      btn.addEventListener("click", () => {
        _cmpAnalysis = btn.getAttribute("data-analysis");
        _buildComparison();
      });
    });
  }

  // Reflect active analysis on the buttons + relevance of controls
  document.querySelectorAll("#cmp-analysis-tabs .cmp-atab").forEach((btn) => {
    btn.classList.toggle("active", btn.getAttribute("data-analysis") === _cmpAnalysis);
  });
  // Jaccard ignores the value basis (presence only)
  if (valueSel) valueSel.disabled = _cmpAnalysis === "jaccard";
  // Group mode needs a field; sample mode uses field only as a label
  if (fieldLabel) {
    const u = (unitSel || {}).value || "sample";
    fieldLabel.textContent = u === "group" ? "Group by field" : "Label / color field";
    if (fieldSel) fieldSel.style.opacity = u === "group" && !fields.length ? 0.5 : 1;
  }

  const pf = _cmpBuildProfiles();

  if (pf.entries.length < 2) {
    if (chart) chart.innerHTML = "";
    if (noData) noData.style.display = "block";
    if (desc) desc.textContent = "";
    return;
  }
  if (noData) noData.style.display = "none";

  const lvl = pf.level === "genus" ? "genus" : "species";
  const valLbl =
    pf.value === "tass" ? "TASS score" : pf.value === "coverage" ? "coverage (breadth %)" : "read abundance";
  const unitLbl = pf.unit === "group" ? `${_metaKeyLabel(pf.field)} groups` : "samples";

  if (_cmpAnalysis === "cosine") {
    if (desc)
      desc.textContent = `Cosine similarity between ${unitLbl} based on ${lvl}-level ${valLbl} profiles (1 = identical, 0 = no shared signal). ${pf.entries.length} entries × ${pf.taxa.length} taxa.`;
    _cmpRenderSquare(
      chart,
      pf.entries,
      (a, b) => _cmpCosine(pf.profile[a.id], pf.profile[b.id]),
      (v) => v.toFixed(2),
    );
  } else if (_cmpAnalysis === "jaccard") {
    if (desc)
      desc.textContent = `Jaccard overlap between ${unitLbl} based on shared ${lvl}-level hits (presence/absence). ${pf.entries.length} entries × ${pf.taxa.length} taxa.`;
    _cmpRenderSquare(
      chart,
      pf.entries,
      (a, b) => _cmpJaccard(pf.profile[a.id], pf.profile[b.id]),
      (v) => v.toFixed(2),
    );
  } else if (_cmpAnalysis === "enrichment") {
    const enr = _cmpEnrichment(pf);
    const maxAbs = Math.max(
      1,
      ...enr.taxa.map((t) => Math.max(...pf.entries.map((e) => Math.abs(enr.matrix[e.id][t] || 0)))),
    );
    const color = d3.scaleSequential(d3.interpolateRdBu).domain([maxAbs, -maxAbs]);
    if (desc)
      desc.textContent = `Enrichment (log₂ fold-change of ${lvl}-level ${valLbl}, each entry vs. the rest). Red = enriched, blue = depleted. Top ${enr.taxa.length} most differential taxa.`;
    _cmpRenderRect(
      chart,
      pf.entries,
      enr.taxa,
      (id, t) => enr.matrix[id][t] || 0,
      (v) => color(v),
      (v) => (v >= 0 ? "+" : "") + v.toFixed(2),
      [-maxAbs, maxAbs],
    );
  } else {
    // profile heatmap
    const totals = {};
    pf.taxa.forEach((t) => {
      totals[t] = pf.entries.reduce((s, e) => s + ((pf.profile[e.id] && pf.profile[e.id][t]) || 0), 0);
    });
    const topTaxa = [...pf.taxa].sort((a, b) => totals[b] - totals[a]).slice(0, 30);
    const vmax = Math.max(
      1,
      ...topTaxa.map((t) => Math.max(...pf.entries.map((e) => (pf.profile[e.id] && pf.profile[e.id][t]) || 0))),
    );
    const color = d3.scaleSequential(d3.interpolateYlGnBu).domain([0, vmax]);
    const isInt = pf.value === "presence";
    const fmt = (v) => (isInt ? Math.round(v).toLocaleString() : v.toFixed(1));
    if (desc)
      desc.textContent = `${lvl.charAt(0).toUpperCase() + lvl.slice(1)}-level ${valLbl} per ${
        pf.unit === "group" ? "group" : "sample"
      } (top ${topTaxa.length} taxa by total).`;
    _cmpRenderRect(
      chart,
      pf.entries,
      topTaxa,
      (id, t) => (pf.profile[id] && pf.profile[id][t]) || 0,
      (v) => (v > 0 ? color(v) : "#f3f5f8"),
      fmt,
      [0, vmax],
    );
  }
}
