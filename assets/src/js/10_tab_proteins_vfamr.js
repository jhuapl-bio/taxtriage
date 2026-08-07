/* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB: VF / AMR           (data-tab="proteins"  —  hidden if no PROT)
       -     Largest section in the file. Renders three side-by-side panels:
       -        drawProtGenus()    — VF hits per genus (count rows)
       -        drawProtProperty() — group-by-property bar chart
       -        drawProteins()     — top-level dispatcher; calls both above
       -     Click an organism row → _showProtDetail() opens a per-organism
       -     detail view; that detail view internally builds a comparison
       -     bubble plot via _drawProtCompare().
       -     Hot path: _renderProtRows() / _filterProt() — both use the
       -     memoized _getTassPrefixSet() Set for O(1) "is this organism in
       -     TASS?" lookups (replaces a former O(N²) .some() scan).
       -     _dedupRows() de-duplicates incoming protein rows with a fast
       -     non-JSON.stringify key.
═══════════════════════════════════════════════════════════════════════════ */
/* ── Memoized "is this species present in the TASS report?" lookup.
         The VF/AMR table did this per row by scanning all of DATA, which
         with 14k rows × thousands of protein hits exploded to billions of
         ops. Now we build (once per render) a Set of all whitespace- /
         comma-bounded prefixes of every Detected Organism in DATA. Lookup
         per row is then a single Set.has() call. Cache is invalidated when
         DATA.length changes. */
let _TASS_PREFIX_CACHE = { len: -1, set: null };
function _getTassPrefixSet() {
  if (_TASS_PREFIX_CACHE.len === DATA.length && _TASS_PREFIX_CACHE.set) {
    return _TASS_PREFIX_CACHE.set;
  }
  const set = new Set();
  for (let i = 0; i < DATA.length; i++) {
    const o = (DATA[i]["Detected Organism"] || "").trim().toLowerCase();
    if (!o) continue;
    set.add(o);
    // Add every whitespace- or comma-bounded leading prefix so that
    // a species stem like "escherichia coli" matches a longer org
    // string like "escherichia coli k-12" without needing startsWith.
    for (let j = 0; j < o.length; j++) {
      const ch = o.charCodeAt(j);
      if (ch === 32 /* space */ || ch === 44 /* comma */) set.add(o.slice(0, j));
    }
  }
  _TASS_PREFIX_CACHE = { len: DATA.length, set };
  return set;
}

/* ── Microbial category of the organism a VF/AMR hit is linked to ───────────
       A hit carries the annotation database's Species / Genus, not a detection
       row, so the microbial category lives on the DATA side. Build the same kind
       of memoized prefix index _getTassPrefixSet() uses — organism name (and each
       whitespace/comma-bounded leading prefix, so "escherichia coli" resolves a
       detection recorded as "Escherichia coli O157:H7") → Set of categories seen
       under that name. Genus names are indexed too, as a coarser fallback.
       Cache is keyed on DATA.length, matching the TASS set's invalidation. */
let _PROT_CAT_CACHE = { len: -1, map: null };
function _getOrgCategoryMap() {
  if (_PROT_CAT_CACHE.len === DATA.length && _PROT_CAT_CACHE.map) return _PROT_CAT_CACHE.map;
  const map = new Map();
  const add = (name, cat) => {
    if (!name) return;
    let s = map.get(name);
    if (!s) map.set(name, (s = new Set()));
    s.add(cat);
  };
  for (let i = 0; i < DATA.length; i++) {
    const r = DATA[i];
    const cat = r["Microbial Category"] || "Unknown";
    const o = (r["Detected Organism"] || "").trim().toLowerCase();
    if (o) {
      add(o, cat);
      for (let j = 0; j < o.length; j++) {
        const ch = o.charCodeAt(j);
        if (ch === 32 || ch === 44) add(o.slice(0, j), cat);
      }
    }
    const g = (r["Genus"] || "").trim().toLowerCase();
    if (g) add(g, cat);
  }
  _PROT_CAT_CACHE = { len: DATA.length, map };
  return map;
}
/* make_report.py stamps every VF/AMR row with a `pathogen` descriptor resolved
     from the row's Source ID against the pathogen sheet, by canonical taxid:
       {match, name, class, status, hc, genus, taxid}
     `class` is exactly the microbial-category vocabulary the filter-mc control
     uses (primary / commensal / opportunistic / potential). This is authoritative
     and, crucially, independent of whether the organism ever produced a detection
     row — so a primary pathogen that never aligned, or that fell below the TASS
     threshold, still classifies as Primary instead of collapsing to Unknown. */
const _PROT_CLASS_LABEL = {
  primary: "Primary",
  commensal: "Commensal",
  opportunistic: "Opportunistic",
  potential: "Potential",
};
function _protPathogenClass(r) {
  if (!r) return null;
  // Charts read the raw PROT rows (`pathogen`); the table reads rows that have
  // been through _protExpandRow (`_pathogen` + a flattened "Pathogen Class").
  const p = r._pathogen || r.pathogen;
  let c = p && typeof p === "object" ? p.class : null;
  if (!c) c = r["Pathogen Class"];
  c = String(c || "")
    .trim()
    .toLowerCase();
  if (!c) return null;
  return _PROT_CLASS_LABEL[c] || c.charAt(0).toUpperCase() + c.slice(1);
}
/** True when this hit is a known pathogen of a class we must never filter away
       on the grounds of "not detected" — it has a canonical pathogen-sheet entry. */
function _protIsPathogenHit(r) {
  return !!_protPathogenClass(r);
}
/** Categories attributable to a hit row. null → nothing could be resolved.
       Order matters: the server-side pathogen stamp wins over guessing from the
       detections table, because the detection may not exist at all. */
function _protRowCategories(r) {
  const stamped = _protPathogenClass(r);
  if (stamped) return new Set([stamped]);
  const map = _getOrgCategoryMap();
  const sp = String(r["Species"] || r["species"] || "")
    .trim()
    .toLowerCase();
  if (sp && map.has(sp)) return map.get(sp);
  const org = String(r["Organism"] || r["Reference Organism"] || "")
    .trim()
    .toLowerCase();
  if (org && map.has(org)) return map.get(org);
  const gn = String(r["Genus"] || r["genus"] || "")
    .trim()
    .toLowerCase();
  if (gn && map.has(gn)) return map.get(gn);
  return null;
}
/** The current global Microbial Category selection, or null when unrestricted. */
function _protCategoryFilter() {
  const el = document.getElementById("filter-mc");
  if (!el || !el.selectedOptions) return null;
  const sel = Array.from(el.selectedOptions).map((o) => o.value);
  return sel.length ? new Set(sel) : null; // nothing selected = show everything
}
/** Does this hit's linked organism fall in the selected categories?
       A hit with no matching detection (an "Ext" row) has no category of its own,
       so it is treated as "Unknown" — exactly how filteredData() treats a
       detection row with a blank Microbial Category. */
function _protRowInCategories(r, catFilter) {
  if (!catFilter) return true;
  const cats = _protRowCategories(r);
  if (!cats || !cats.size) return catFilter.has("Unknown");
  for (const c of cats) if (catFilter.has(c)) return true;
  return false;
}

/* ── Sample identity for a hit row ──────────────────────────────────────────
       The merged XLSX names the sample column "Specimen ID" in the Per-Gene Hits
       sheet but "Sample" in the AMR Genes sheet (create_report.py writes them from
       different frames). The charts read PROT.* raw — only the table runs rows
       through _applyProtColRemap — so reading "Specimen ID" alone silently gave
       every AMR row an empty sample id, and every AMR row was then dropped from
       both charts. Accept every spelling here. */
function _protRowSample(r) {
  return r["Specimen ID"] || r["Sample"] || r["sample"] || r["sample_name"] || "";
}
function _protRowGenus(r) {
  return String(r["Genus"] || r["genus"] || "").trim() || "Unknown";
}
// Clean category labels used when a row carries no Property of its own. The AMR
// Genes sheet has no Property column at all (it has Classification, which is a
// long semicolon-joined string), so those rows used to bucket as "Unknown" and
// sat in the legend as a separate grey series next to the CARD hits they belong
// with. Fall back to the row's category label instead.
const _PROT_CATEGORY_LABEL = {
  amr: "Antibiotic Resistance",
  vf: "Virulence Factor",
  transporter: "Transporter",
  drug_target: "Drug Target",
  other: "Unknown",
};
function _protRowProp(r) {
  const explicit = r["Property"] || r["Class"] || r["property"];
  if (explicit) return explicit;
  return _PROT_CATEGORY_LABEL[protPropCategory("", r)] || "Unknown";
}

/* ── Shared visibility gate for the two VF/AMR charts ───────────────────────
       Previously both charts intersected hits against the set of
       `Specimen ID||Genus` pairs in filteredData(). Any hit whose genus had no
       matching detection row was discarded without a trace, so the charts read
       "No genus annotation data." while the table below them listed thousands of
       hits. That happens routinely:
         • de-novo / unaligned annotation (no reference row exists at all)
         • taxonomy renames — the detections table says "Orthoebolavirus" while
           the annotation DB still says "Ebolavirus"
         • any hit attributed to a genus below the TASS threshold
       The charts now follow the same rule the Annotation Table already uses: the
       `prot-tass-filter` control ("All hits" by default). Sample visibility is
       always honoured; the in-TASS restriction is applied only when the user
       explicitly asks for it, via the same _getTassPrefixSet() lookup the table
       uses, so all three views agree. */
function _protVisibleSamples() {
  const s = new Set();
  filteredData().forEach((r) => {
    const raw = r["Specimen ID"];
    s.add(typeof specimenKey === "function" ? specimenKey(raw) : String(raw == null ? "" : raw));
  });
  return s;
}
function _protHitFilter() {
  const visibleSamples = _protVisibleSamples();
  const tassMode = (document.getElementById("prot-tass-filter") || {}).value || "all";
  const tassSet = tassMode === "all" ? null : _getTassPrefixSet();
  const catFilter = _protCategoryFilter();
  return (r) => {
    const raw = _protRowSample(r);
    if (sampleHidden[raw]) return false;
    const key = typeof specimenKey === "function" ? specimenKey(raw) : String(raw);
    if (!visibleSamples.has(key)) return false;
    if (!_protRowInCategories(r, catFilter)) return false;
    if (!tassSet) return true;
    const sp = String(r["Species"] || r["species"] || "")
      .trim()
      .toLowerCase();
    const gn = _protRowGenus(r).toLowerCase();
    const inTass = (sp && tassSet.has(sp)) || (gn !== "unknown" && tassSet.has(gn));
    return tassMode === "in" ? inTass : !inTass;
  };
}
/** Every VF/AMR hit row that passes the shared gate, as one flat array. */
function _protVisibleHits() {
  const keep = _protHitFilter();
  const out = [];
  (PROT.per_gene_hits || []).forEach((r) => {
    if (keep(r)) out.push(r);
  });
  (PROT.amr_genes || []).forEach((r) => {
    if (keep(r)) out.push(r);
  });
  return out;
}
/** Total hit rows present before any filtering — used for the empty-state copy. */
function _protTotalHits() {
  return (PROT.per_gene_hits || []).length + (PROT.amr_genes || []).length;
}

/* ── Shared "does this hit survive the UI filters?" predicate ───────────────
     The Summary tab reports VF/AMR counts of its own (the per-row chips and the
     "no passing detections but N VF/AMR hits" indicator). Those used to count
     every raw row in PROT, so they disagreed with the VF/AMR tab as soon as any
     filter was applied — reporting Drug Target hits against host organisms
     (DrugBank targets are human/bovine proteins, genus "Homo"/"Bos") that the
     VF/AMR tab hides by default. Both surfaces now share this predicate.
     Only the category-level filters are applied here — sample visibility is
     handled by the caller, which knows its own sample context. */
function _protHitPassesUiFilters(r, opts) {
  const o = opts || {};
  const prop = _protRowProp(r);
  if (prop && PROT_HIDDEN_PROPS.has(prop)) return false;
  if (!_protRowInCategories(r, o.catFilter === undefined ? _protCategoryFilter() : o.catFilter)) return false;
  const minPid = o.minPid === undefined ? _protMinPid() : o.minPid;
  if (minPid > 0) {
    const pid = parseFloat(r["%id"] ?? r["pident"] ?? r["%ID"]);
    if (!isNaN(pid) && pid < minPid) return false;
  }
  return true;
}
function _protMinPid() {
  const v = parseFloat((document.getElementById("prot-pid-min") || {}).value);
  return isNaN(v) ? 0 : v;
}
/* Cache key covering everything _protHitPassesUiFilters depends on, so callers
     that memoize their own rollups know when to recompute. */
function _protUiFilterKey() {
  const cats = _protCategoryFilter();
  return [[...PROT_HIDDEN_PROPS].sort().join(","), cats ? [...cats].sort().join(",") : "*", _protMinPid()].join("|");
}

/* On first draw, hide every property that is not a Virulence Factor or an AMR
     determinant (Drug Target, Transporter, unclassified). Runs once per dataset so
     that redraws never clobber the user's own legend clicks; _resetProtHiddenDefaults()
     re-arms it when a new dataset is uploaded. */
// property -> one representative hit row, so the legend tooltip can re-derive the
// property's category (which needs the row's Source column, not just its name).
let _protPropRow = {};
function _applyProtHiddenDefaults(hits) {
  const seen = new Map(); // property -> representative row (for the Source signal)
  (hits || []).forEach((r) => {
    const p = _protRowProp(r);
    if (p && !seen.has(p)) seen.set(p, r);
  });
  (PROT.genus_summary || []).forEach((r) => {
    const p = r["Property"] || "Unknown";
    if (p && !seen.has(p)) seen.set(p, r);
  });
  _protPropRow = Object.fromEntries(seen);
  if (PROT_HIDDEN_DEFAULTS_APPLIED) return;
  PROT_HIDDEN_DEFAULTS_APPLIED = true;
  seen.forEach((row, p) => {
    if (!protPropVisibleByDefault(p, row)) PROT_HIDDEN_PROPS.add(p);
  });
  // Never let the default blank the chart out: if nothing at all was classified
  // as VF/AMR, show everything rather than an empty panel.
  if (seen.size && PROT_HIDDEN_PROPS.size >= seen.size) PROT_HIDDEN_PROPS.clear();
}

function drawProteins() {
  drawProtGenus();
  drawProtProperty();
  if (window._initProtTable) window._initProtTable();
}

function drawProtGenus() {
  const wrap = document.getElementById("prot-genus-svg");
  wrap.innerHTML = "";

  // Build aggregation from per-row hit data so sampleHidden filtering works.
  // Fall back to pre-aggregated genus_summary only when no per-row data exists.
  const hasByRow = _protTotalHits() > 0;
  const agg = {};
  const propSet = new Set();

  if (hasByRow) {
    const visibleHits = _protVisibleHits();
    _applyProtHiddenDefaults(visibleHits);
    visibleHits.forEach((r) => {
      const genus = _protRowGenus(r);
      const prop = _protRowProp(r);
      propSet.add(prop);
      if (!agg[genus]) agg[genus] = {};
      agg[genus][prop] = (agg[genus][prop] || 0) + 1;
    });
  } else {
    (PROT.genus_summary || []).forEach((r) => {
      const genus = r["Genus"] || "Unknown";
      const prop = r["Property"] || "Unknown";
      propSet.add(prop);
      if (!agg[genus]) agg[genus] = {};
      agg[genus][prop] = (agg[genus][prop] || 0) + (parseInt(r["# Hits"]) || 0);
    });
    _applyProtHiddenDefaults([]);
  }

  if (!Object.keys(agg).length) {
    // Distinguish "nothing was ever loaded" from "the current filters hid it all"
    // — the old message said the former in both cases, which is what made this
    // failure impossible to diagnose from the report itself.
    const total = _protTotalHits();
    wrap.innerHTML = total
      ? `<p style="color:#999">No genus annotation data for the current filters ` +
        `— all ${total.toLocaleString()} VF/AMR hit${total === 1 ? "" : "s"} are filtered out. ` +
        `Check the sample toggles, the detections filters, and the “${
          (document.getElementById("prot-tass-filter") || {}).value === "in" ? "In TASS report" : "Not in TASS report"
        }” setting on the Annotation Table.</p>`
      : '<p style="color:#999">No genus annotation data.</p>';
    return;
  }

  const genera = Object.keys(agg).sort();
  const props = [...propSet].filter(Boolean).sort();
  const activeProps = props.filter((p) => !PROT_HIDDEN_PROPS.has(p));

  const marginL = 130,
    marginT = 36,
    marginR = 20,
    marginB = 80;
  const W = Math.max(600, wrap.clientWidth || 900);
  const H = 340;
  const iW = W - marginL - marginR;
  const iH = H - marginT - marginB;

  const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H);

  // Clip bars to plot area
  svg
    .append("defs")
    .append("clipPath")
    .attr("id", "pg-clip")
    .append("rect")
    .attr("width", iW)
    .attr("height", iH + 4);

  const g = svg.append("g").attr("transform", `translate(${marginL},${marginT})`);

  const scaleType = (document.getElementById("prot-genus-scale") || {}).value || "log";

  const x0 = d3.scaleBand().domain(genera).range([0, iW]).paddingInner(0.15).paddingOuter(0.1);
  const x1 = d3.scaleBand().domain(activeProps).range([0, x0.bandwidth()]).padding(0.05);
  const yMaxRaw = d3.max(genera, (gn) => d3.sum(activeProps, (p) => agg[gn][p] || 0)) || 0;
  const yMax = yMaxRaw || 10;

  let y;
  if (scaleType === "log") {
    y = d3
      .scaleLog()
      .domain([1, Math.max(yMax, 2)])
      .range([iH, 0])
      .clamp(true);
  } else if (scaleType === "sqrt") {
    y = d3.scaleSqrt().domain([0, yMax]).range([iH, 0]).nice();
  } else {
    y = d3.scaleLinear().domain([0, yMax]).range([iH, 0]).nice();
  }
  const color = d3.scaleOrdinal(d3.schemePastel1).domain(props);

  // Y axis (static)
  const yAxisCall = scaleType === "log" ? d3.axisLeft(y).ticks(5, "~s") : d3.axisLeft(y).ticks(5);
  g.append("g").attr("class", "axis").call(yAxisCall);

  // Y axis label (rotated)
  g.append("text")
    .attr("transform", "rotate(-90)")
    .attr("x", -iH / 2)
    .attr("y", -marginL + 14)
    .attr("text-anchor", "middle")
    .style("font-size", "0.74em")
    .style("fill", "#555")
    .text("Hit Count");

  // X axis group — updated on zoom
  const xAxisG = g.append("g").attr("class", "axis").attr("transform", `translate(0,${iH})`);

  // Zoom-capture rect (behind bars)
  const zoomRect = g
    .append("rect")
    .attr("width", iW)
    .attr("height", iH)
    .attr("fill", "none")
    .attr("pointer-events", "all")
    .attr("cursor", "grab");

  // Bars group (clipped)
  const barsG = g.append("g").attr("clip-path", "url(#pg-clip)");

  function _drawBars(x0c, x1c) {
    barsG.selectAll("g.pg-grp").remove();
    genera.forEach((gn) => {
      const gg = barsG
        .append("g")
        .attr("class", "pg-grp")
        .attr("transform", `translate(${x0c(gn)},0)`)
        .style("cursor", "pointer")
        .on("click", () => _showProtDetail(gn, null));
      activeProps.forEach((p) => {
        const v = agg[gn][p] || 0;
        if (!v) return;
        // For log scale clamp value to domain minimum (1) so bar is always visible
        const yVal = scaleType === "log" ? Math.max(1, v) : v;
        const barTop = y(yVal);
        const barH = Math.max(1, iH - barTop);
        gg.append("rect")
          .attr("x", x1c(p))
          .attr("width", Math.max(0, x1c.bandwidth()))
          .attr("y", barTop)
          .attr("height", barH)
          .attr("fill", color(p))
          .attr("rx", 2)
          .on("mouseover", (ev) =>
            showTip(`<b>${gn}</b> — ${p}<br>Hits: <b>${v}</b><br><small>Click for details</small>`, ev),
          )
          .on("mousemove", moveTip)
          .on("mouseout", hideTip);
      });
    });
  }

  function _updateXAxis(x0c) {
    xAxisG
      .call(d3.axisBottom(x0c).tickSizeOuter(0))
      .selectAll(".tick")
      .each(function (d) {
        const xPos = x0c(d) + x0c.bandwidth() / 2;
        const outOfView = xPos < 0 || xPos > iW;
        const hasHits = activeProps.some((p) => (agg[d] || {})[p] > 0);
        const hide = outOfView || !hasHits || yMaxRaw <= 0;
        d3.select(this).style("display", hide ? "none" : "");
      })
      .selectAll("text")
      .attr("transform", "rotate(-30)")
      .style("text-anchor", "end");
  }

  _updateXAxis(x0);
  if (activeProps.length) {
    _drawBars(x0, x1);
  } else {
    svg
      .append("text")
      .attr("x", marginL + iW / 2)
      .attr("y", marginT + iH / 2)
      .attr("text-anchor", "middle")
      .attr("font-size", 12)
      .attr("fill", "#888")
      .text("All properties hidden — use legend to re-enable");
  }

  // D3 zoom — rescale band x0 range
  const zoom = d3
    .zoom()
    .scaleExtent([1, 10])
    .translateExtent([
      [0, 0],
      [iW, iH],
    ])
    .on("zoom", (event) => {
      const t = event.transform;
      x0.range([t.applyX(0), t.applyX(iW)]);
      x1.range([0, x0.bandwidth()]);
      _updateXAxis(x0);
      if (activeProps.length) {
        _drawBars(x0, x1);
      }
    });

  zoomRect.call(zoom);

  // Reset zoom button
  svg
    .append("text")
    .attr("x", marginL + iW - 2)
    .attr("y", marginT - 8)
    .attr("text-anchor", "end")
    .attr("font-size", 10)
    .attr("fill", "#1565c0")
    .style("cursor", "pointer")
    .text("⟳ Reset zoom")
    .on("click", () => zoomRect.call(zoom.transform, d3.zoomIdentity));

  // One-click escape from the VF/AMR-only default. Without this the only way to
  // get Drug Target / Transporter back is to find and click each legend square.
  const _anyHidden = props.some((p) => PROT_HIDDEN_PROPS.has(p));
  svg
    .append("text")
    .attr("x", marginL + iW - 84)
    .attr("y", marginT - 8)
    .attr("text-anchor", "end")
    .attr("font-size", 10)
    .attr("fill", "#1565c0")
    .style("cursor", "pointer")
    .text(_anyHidden ? "▣ Show all categories" : "▢ VF/AMR only")
    .on("mouseover", (ev) =>
      showTip(
        _anyHidden
          ? "Show every annotation category<br><small>Drug Target, Transporter and unclassified hits are hidden by default</small>"
          : "Show only Virulence Factor and AMR categories",
        ev,
      ),
    )
    .on("mousemove", moveTip)
    .on("mouseout", hideTip)
    .on("click", () => {
      hideTip();
      if (_anyHidden) PROT_HIDDEN_PROPS.clear();
      else
        props.forEach((p) => {
          if (!protPropVisibleByDefault(p, _protPropRow[p])) PROT_HIDDEN_PROPS.add(p);
        });
      drawProtGenus();
      drawProtProperty();
      if (window._renderProtCatLegend) window._renderProtCatLegend();
      if (window._filterProtExternal) window._filterProtExternal();
    });

  // Legend — above bars (toggleable)
  const legendG = svg.append("g").attr("class", "pg-legend");
  props.forEach((p, i) => {
    const lx = marginL + i * 160;
    const isHidden = PROT_HIDDEN_PROPS.has(p);
    const legItem = legendG
      .append("g")
      .attr("transform", `translate(${lx},0)`)
      .style("cursor", "pointer")
      .on("click", () => {
        if (PROT_HIDDEN_PROPS.has(p)) PROT_HIDDEN_PROPS.delete(p);
        else PROT_HIDDEN_PROPS.add(p);
        drawProtGenus();
        drawProtProperty();
        // Keep the table's category chips in step with this legend.
        if (window._renderProtCatLegend) window._renderProtCatLegend();
        if (window._filterProtExternal) window._filterProtExternal();
      })
      .on("mouseover", (ev) => {
        const action = PROT_HIDDEN_PROPS.has(p) ? "Enable" : "Disable";
        const cat = protPropCategory(p, (_protPropRow && _protPropRow[p]) || null);
        const why =
          PROT_HIDDEN_PROPS.has(p) && !PROT_DEFAULT_VISIBLE_CATEGORIES.has(cat)
            ? "<br><small>Hidden by default — only Virulence Factor and AMR categories are shown on load.</small>"
            : "";
        showTip(`${action} ${p}<br><small>Click to ${action.toLowerCase()} in chart and table</small>${why}`, ev);
      })
      .on("mousemove", moveTip)
      .on("mouseout", hideTip);

    legItem
      .append("rect")
      .attr("x", 0)
      .attr("y", 4)
      .attr("width", 12)
      .attr("height", 12)
      .attr("fill", color(p))
      .attr("rx", 2)
      .attr("opacity", isHidden ? 0.25 : 1);
    legItem
      .append("text")
      .attr("x", 16)
      .attr("y", 14)
      .attr("font-size", 11)
      .attr("fill", isHidden ? "#aaa" : "#222")
      .text(p + (isHidden ? " (hidden)" : ""));
  });
}

function drawProtProperty() {
  const wrap = document.getElementById("prot-prop-svg");
  wrap.innerHTML = "";

  // Same shared gate as drawProtGenus — see _protHitFilter() above.
  const _allHits = _protVisibleHits();
  const _propGenus = {};
  _allHits.forEach((r) => {
    const prop = _protRowProp(r);
    const gen = _protRowGenus(r);
    if (!_propGenus[prop]) _propGenus[prop] = {};
    _propGenus[prop][gen] = (_propGenus[prop][gen] || 0) + 1;
  });

  // Always build from per-row hits for accurate sample-visibility filtering.
  // Only use pre-aggregated metadata_counts when no per-row data exists.
  const hasByRow = _protTotalHits() > 0;
  let propRows;
  if (hasByRow) {
    const agg = {};
    _allHits.forEach((r) => {
      const prop = _protRowProp(r);
      agg[prop] = (agg[prop] || 0) + 1;
    });
    propRows = Object.entries(agg).map(([value, count]) => ({ field: "property", value, count }));
  } else {
    const counts = PROT.metadata_counts || [];
    propRows = counts.filter((r) => r["field"] === "property");
  }
  if (!propRows.length) {
    const total = _protTotalHits();
    wrap.innerHTML = total
      ? `<p style="color:#999">No property metadata for the current filters — all ${total.toLocaleString()} hits are filtered out.</p>`
      : '<p style="color:#999">No property metadata.</p>';
    return;
  }

  // Respect the legend's hidden-category set so this panel, the genus chart and
  // the table all show the same categories (VF/AMR only until the user says otherwise).
  const data = propRows
    .filter((r) => !PROT_HIDDEN_PROPS.has(r["value"]))
    .map((r) => ({ label: r["value"] || "", value: parseInt(r["count"]) || 0 }))
    .sort((a, b) => b.value - a.value)
    .slice(0, 20);
  if (!data.length) {
    wrap.innerHTML =
      '<p style="color:#999">All annotation categories are hidden — click “Show all categories” on the genus chart above.</p>';
    return;
  }

  // Dynamic row height — ensure labels never overlap
  const rowH = Math.max(24, Math.min(40, Math.floor(500 / Math.max(data.length, 1))));
  const marginL = Math.max(140, data.reduce((m, d) => Math.max(m, d.label.length * 6.2), 0) + 10);
  const marginT = 20,
    marginR = 60,
    marginB = 36;
  const W = Math.max(420, wrap.clientWidth || 600);
  const chartH = data.length * rowH;
  const H = chartH + marginT + marginB;
  const iW = W - marginL - marginR;
  const xMax = d3.max(data, (d) => d.value) || 1;

  const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H).style("overflow", "visible");

  // Clip path so bars don't overflow on zoom
  svg.append("defs").append("clipPath").attr("id", "prop-clip").append("rect").attr("width", iW).attr("height", chartH);

  const g = svg.append("g").attr("transform", `translate(${marginL},${marginT})`);

  // Zoom capture rect — inserted BEFORE bars so bars remain on top in z-order
  // and receive mouse events; wheel/scroll events bubble through to the rect anyway
  const zoomRect = g
    .append("rect")
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", iW)
    .attr("height", chartH)
    .attr("fill", "none")
    .attr("pointer-events", "all")
    .attr("cursor", "grab");

  // Clipped group for bars + value labels
  const barG = g.append("g").attr("clip-path", "url(#prop-clip)");

  const x = d3.scaleLinear().domain([0, xMax]).range([0, iW]).nice();
  const y = d3
    .scaleBand()
    .domain(data.map((d) => d.label))
    .range([0, chartH])
    .padding(0.2);

  const xAxisG = g
    .append("g")
    .attr("class", "axis")
    .attr("transform", `translate(0,${chartH})`)
    .call(d3.axisBottom(x).ticks(5));
  g.append("g").attr("class", "axis").call(d3.axisLeft(y).tickSize(0)).select(".domain").remove();
  svg
    .append("text")
    .attr("x", marginL + iW / 2)
    .attr("y", H - 4)
    .attr("text-anchor", "middle")
    .style("font-size", "10px")
    .attr("fill", "#666")
    .text("Hit count");

  const color = d3.scaleOrdinal(d3.schemeTableau10);
  data.forEach((d, i) => {
    const genusCounts = _propGenus[d.label] || {};
    const topGenus = Object.entries(genusCounts)
      .sort((a, b) => b[1] - a[1])
      .map(([g, c]) => `${g}: ${c}`)
      .join("<br>");
    const genusTip = topGenus ? `<br><small>${topGenus}</small>` : "";
    barG
      .append("rect")
      .datum(d)
      .attr("width", x(d.value))
      .attr("y", y(d.label))
      .attr("height", y.bandwidth())
      .attr("fill", color(i))
      .attr("rx", 2)
      .on("mouseover", (ev) => showTip(`<b>${d.label}</b>: ${d.value} hits${genusTip}`, ev))
      .on("mousemove", moveTip)
      .on("mouseout", hideTip);
  });

  // D3 zoom on X axis — scroll/pinch to zoom, double-click to reset
  const zoom = d3
    .zoom()
    .scaleExtent([1, 10])
    .translateExtent([
      [0, 0],
      [iW, chartH],
    ])
    .extent([
      [0, 0],
      [iW, chartH],
    ])
    .on("zoom", (event) => {
      const newX = event.transform.rescaleX(x);
      const dom = newX.domain();
      if (dom[0] < 0) {
        const span = dom[1] - dom[0];
        newX.domain([0, span]);
      }
      xAxisG.call(d3.axisBottom(newX).ticks(5));
      // Check for label overlap → rotate if needed
      const ticks = xAxisG.selectAll("text");
      let overlap = false;
      const tickNodes = ticks.nodes();
      for (let i = 1; i < tickNodes.length; i++) {
        const a = tickNodes[i - 1].getBoundingClientRect();
        const b = tickNodes[i].getBoundingClientRect();
        if (a.right > b.left - 2) {
          overlap = true;
          break;
        }
      }
      ticks.attr("transform", overlap ? "rotate(-30)" : null).style("text-anchor", overlap ? "end" : "middle");
      // Redraw bars
      barG
        .selectAll("rect")
        .attr("x", (d) => newX(0))
        .attr("width", (d) => Math.max(0, newX(d.value) - newX(0)));
    });

  zoomRect.call(zoom);

  // Reset on double-click
  svg.on("dblclick.zoom", () => {
    svg.transition().duration(400).call(zoom.transform, d3.zoomIdentity);
  });

  // Zoom hint
  svg
    .append("text")
    .attr("x", marginL + iW)
    .attr("y", marginT - 5)
    .attr("text-anchor", "end")
    .style("font-size", "9px")
    .attr("fill", "#aaa")
    .text("scroll/pinch to zoom x-axis · dbl-click to reset");
}

// Deduplicate rows by all non-internal fields — shared by table, panel, and chart.
// Builds the dedup key by direct string concat instead of constructing a new
// object + JSON.stringify per row. With thousands of protein/AMR rows the old
// approach showed up as a measurable chunk of VF/AMR tab load time.
/* Render any cell value as text. `taxids` and `pathogen` are nested objects on
     every hit row; they used to reach td.textContent verbatim and display as the
     literal "[object Object]", hiding the taxids that decide a hit's
     classification. _protExpandRow() flattens the two known ones into scalar
     columns; this is the backstop for anything else. */
function _protCellText(v) {
  if (v == null) return "";
  if (typeof v !== "object") return String(v);
  if (Array.isArray(v)) return v.map(_protCellText).filter(Boolean).join(", ");
  return Object.entries(v)
    .filter(([, val]) => val !== null && val !== undefined && val !== "")
    .map(([k, val]) => `${k}: ${_protCellText(val)}`)
    .join(" · ");
}

function _dedupRows(rows) {
  const seen = new Set();
  const out = [];
  for (let i = 0; i < rows.length; i++) {
    const r = rows[i];
    const keys = Object.keys(r).sort();
    let key = "";
    for (let j = 0; j < keys.length; j++) {
      const k = keys[j];
      if (k.charCodeAt(0) === 95 /* '_' */) continue;
      const v = r[k];
      key += k + "\u0001" + (v == null ? "" : v) + "\u0002";
    }
    if (seen.has(key)) continue;
    seen.add(key);
    out.push(r);
  }
  return out;
}

/* Protein annotation table with search */
(function () {
  let _protAllRows = [];
  let _protCols = [];
  let _protSortCol = null;
  let _protSortAsc = true;
  // ── Pagination state ─────────────────────────────────────────────────────
  // _protPageRows holds the full filtered+sorted result; only the current slice
  // is ever built into the DOM. Page size 0 means "All".
  let _protPageRows = [];
  let _protPage = 1;
  let _protPageSize = 100;

  // Category color map (Property field)
  const _catColors = {
    "Virulence Factor": "#e53935",
    "Antibiotic Resistance": "#fb8c00",
    AMR: "#fb8c00",
    "Drug Target": "#8e24aa",
    Transporter: "#00897b",
  };
  function _catColor(row) {
    // Go through _protRowProp so AMR-sheet rows (no Property column) resolve to
    // "Antibiotic Resistance" and get the same swatch as the CARD hits they
    // belong with, rather than falling through to the grey "Other".
    const prop = _protRowProp(row) || row["_source"] || "";
    return _catColors[prop] || _catColors[row["_source"]] || "#90a4ae";
  }

  // Keyword groups for the dropdown filter — matched against Property, Classification,
  // Gene, Product, Description, and Annotation fields (case-insensitive substring).
  const PROT_KEYWORD_GROUPS = {
    abx: [
      "antibiotic",
      "antimicrobial",
      "beta-lactam",
      "carbapenem",
      "tetracycline",
      "aminoglycoside",
      "macrolide",
      "quinolone",
      "fluoroquinolone",
      "vancomycin",
      "methicillin",
      "oxacillin",
      "rifampin",
      "colistin",
      "polymyxin",
      "resistance",
      "amr",
      "drug resistance",
    ],
    toxin: [
      "toxin",
      "exotoxin",
      "endotoxin",
      "enterotoxin",
      "cytotoxin",
      "leukotoxin",
      "hemolysin",
      "hemolysis",
      "cytolysin",
      "stx",
      "shiga",
      "pertussis toxin",
      "diphtheria",
      "botulinum",
      "clostridial",
    ],
    vf: ["virulence", "pathogenicity", "pathogenic", "virulence factor"],
    efflux: ["efflux", "efflux pump", "mdr", "multidrug", "ABC transporter", "RND", "MFS transporter", "drug efflux"],
    adhesin: [
      "adhesin",
      "invasion",
      "invasin",
      "fimbriae",
      "fimbrial",
      "pilus",
      "pili",
      "attachment",
      "colonization",
      "fim",
      "type iv pili",
    ],
    immune: [
      "immune evasion",
      "immune",
      "capsule",
      "complement",
      "serum resistance",
      "phagocytosis",
      "antiphagocytic",
      "opsonization",
      "IgA protease",
      "type iii secretion effector",
    ],
    iron: [
      "iron",
      "siderophore",
      "ferritin",
      "ferric",
      "fur",
      "iuc",
      "iut",
      "enterobactin",
      "aerobactin",
      "yersiniabactin",
      "catechol",
    ],
    secretion: [
      "secretion system",
      "type iii",
      "type-iii",
      "t3ss",
      "t4ss",
      "t6ss",
      "type vi",
      "type iv",
      "type ii secretion",
      "sec pathway",
      "tat pathway",
      "needle complex",
    ],
    biofilm: [
      "biofilm",
      "exopolysaccharide",
      "eps",
      "pellicle",
      "quorum sensing",
      "curli",
      "cellulose",
      "polysaccharide synthesis",
    ],
    mobile: [
      "transposon",
      "integron",
      "plasmid",
      "insertion sequence",
      "IS element",
      "mobile genetic",
      "conjugation",
      "horizontal gene transfer",
      "integrase",
      "recombinase",
      "phage",
    ],
  };

  // Fields to check when applying keyword group filter
  const _KW_FIELDS = [
    "Property",
    "property",
    "Classification",
    "classification",
    "Gene",
    "gene",
    "Gene Name",
    "gene_name",
    "Product",
    "product",
    "Antibiotics",
    "Description",
    "description",
    "Annotation",
    "annotation",
    "Function",
    "function",
    "Class",
    "class",
  ];

  /* ── Flatten the nested payloads into real columns ────────────────────────
         make_report.py attaches two objects to every hit: `taxids`
         ({species_taxid, taxon_id, genus_taxid}) and `pathogen`
         ({match, name, class, status, hc, genus, taxid}). _protCols is built from
         Object.keys(), so both became columns whose cells were rendered with
         td.textContent = <object> — i.e. the literal string "[object Object]",
         hiding the very taxids that decide how a hit is classified. Expand them
         into scalar columns and drop the raw objects.
         The originals are kept under underscore-prefixed keys, which _protCols
         filters out, so _protPathogenClass() and friends keep working. */
  function _protExpandRow(r) {
    const t = r.taxids && typeof r.taxids === "object" ? r.taxids : null;
    const p = r.pathogen && typeof r.pathogen === "object" ? r.pathogen : null;
    if (!t && !p) return r;
    const out = { ...r };
    delete out.taxids;
    delete out.pathogen;
    if (t) {
      out._taxids = t;
      out["Species Taxid"] = t.species_taxid || t.taxon_id || "";
      out["Genus Taxid"] = t.genus_taxid || "";
    }
    if (p) {
      out._pathogen = p;
      out["Pathogen"] = p.name || "";
      out["Pathogen Class"] = _PROT_CLASS_LABEL[String(p.class || "").toLowerCase()] || p.class || "";
      out["Pathogen Status"] = p.status || "";
      // How the pathogen-sheet match was made (taxid > name > genus). Worth
      // surfacing: a genus-level match is much weaker than a taxid match.
      out["Pathogen Match"] = p.match || "";
      if (!out["Species Taxid"] && p.taxid) out["Species Taxid"] = p.taxid;
    }
    return out;
  }

  function _buildProtTable() {
    // Combine per_gene_hits with amr_genes, removing exact duplicates
    // Apply column rename map so field names are unified regardless of source
    const sampleSet = new Set(DATA.map((r) => r["Specimen ID"] || "").filter(Boolean));
    const geneRows = _applyProtColRemap(
      (PROT.per_gene_hits || []).filter(
        (r) => !sampleHidden[r["Specimen ID"]] && (!r["Specimen ID"] || sampleSet.has(r["Specimen ID"])),
      ),
    );
    const amrRows = _applyProtColRemap(
      (PROT.amr_genes || [])
        .filter((r) => !sampleHidden[r["Specimen ID"]] && (!r["Specimen ID"] || sampleSet.has(r["Specimen ID"])))
        .map((r) => ({ ...r, _source: "AMR" })),
    );
    // When specimen merge is on, relabel each hit's Specimen ID onto its
    // specimen group before deduping so members of the same specimen collapse
    // into one row instead of listing every raw sample name separately.
    const _relabel = (r) => {
      if (typeof specimenKey !== "function" || !r["Specimen ID"]) return r;
      const grouped = specimenKey(r["Specimen ID"]);
      return grouped === r["Specimen ID"] ? r : { ...r, "Specimen ID": grouped };
    };
    // Expand BEFORE dedup: _dedupRows builds its key by concatenating raw values,
    // so while `taxids` / `pathogen` were still objects they stringified to the
    // constant "[object Object]" and two hits that differed only by taxid were
    // silently collapsed into one row. Flattening first puts the real taxids in
    // the key.
    _protAllRows = _dedupRows([...geneRows.map(_relabel), ...amrRows.map(_relabel)].map(_protExpandRow));
    if (!_protAllRows.length) {
      document.getElementById("prot-table-wrap").style.display = "none";
      return;
    }
    document.getElementById("prot-table-wrap").style.display = "block";

    // Gather columns (include Genus, Species/Organism if present)
    const colSet = new Set();
    _protAllRows.forEach((r) => Object.keys(r).forEach((k) => colSet.add(k)));
    _protCols = [...colSet].filter((c) => !c.startsWith("_"));

    // Build column selector for search
    const colSel = document.getElementById("prot-search-col");
    colSel.innerHTML = '<option value="">All columns</option>';
    _protCols.forEach((c) => {
      const opt = document.createElement("option");
      opt.value = c;
      opt.textContent = c;
      colSel.appendChild(opt);
    });

    _renderProtHeader();
    _renderProtCatLegend();
    _wireProtPager();
    // Go through _filterProt rather than rendering _protAllRows directly, so the
    // initial view already reflects the category / property / %id defaults.
    _filterProt();

    // Assign handlers rather than addEventListener: _buildProtTable re-runs on
    // every global filter change (redraw → _drawTab → drawProteins), and stacked
    // listeners would re-filter the table N times per keystroke.
    const searchEl = document.getElementById("prot-search");
    if (searchEl) searchEl.oninput = _filterProt;
    colSel.onchange = _filterProt;
    const tassFilter = document.getElementById("prot-tass-filter");
    // The two charts honour this control as well, so redraw them alongside the
    // table. Call the draw functions directly rather than drawProteins() — that
    // would re-enter _buildProtTable.
    if (tassFilter)
      tassFilter.onchange = () => {
        _filterProt();
        drawProtGenus();
        drawProtProperty();
      };
    const pidMinEl = document.getElementById("prot-pid-min");
    if (pidMinEl) pidMinEl.oninput = _filterProt;
    const kwGroupEl = document.getElementById("prot-keyword-group");
    if (kwGroupEl) kwGroupEl.onchange = _filterProt;
  }

  /* ── Category legend under the bar plot ──────────────────────────────────
         Was a hardcoded row of five colour swatches that did nothing. It now
         renders the categories actually present in the data and toggles the same
         PROT_HIDDEN_PROPS set the chart legend uses, so clicking a chip here
         hides that category in the table AND both charts (and vice versa). */
  function _renderProtCatLegend() {
    const host = document.getElementById("prot-cat-legend");
    if (!host) return;
    host.innerHTML = "";
    // Count per category across all rows, before the property filter is applied,
    // so a category you have hidden still shows how much it would bring back.
    const counts = new Map();
    const rowFor = new Map();
    _protAllRows.forEach((r) => {
      const p = _protRowProp(r);
      if (!p) return;
      counts.set(p, (counts.get(p) || 0) + 1);
      if (!rowFor.has(p)) rowFor.set(p, r);
    });
    if (!counts.size) return;

    const props = [...counts.keys()].sort();
    props.forEach((p) => {
      const off = PROT_HIDDEN_PROPS.has(p);
      const chip = document.createElement("span");
      chip.className = "prot-cat-chip" + (off ? " prot-cat-off" : "");
      chip.dataset.prop = p;
      chip.title = off ? `Show ${p} (${counts.get(p)} rows)` : `Hide ${p} (${counts.get(p)} rows)`;

      const sw = document.createElement("span");
      sw.className = "prot-cat-sw";
      sw.style.background = _catColor(rowFor.get(p));
      chip.appendChild(sw);

      const lbl = document.createElement("span");
      lbl.textContent = p;
      chip.appendChild(lbl);

      const ct = document.createElement("span");
      ct.className = "prot-cat-ct";
      ct.textContent = `(${counts.get(p).toLocaleString()})`;
      chip.appendChild(ct);

      chip.onclick = () => {
        if (PROT_HIDDEN_PROPS.has(p)) PROT_HIDDEN_PROPS.delete(p);
        else PROT_HIDDEN_PROPS.add(p);
        _syncProtCategoryViews();
      };
      host.appendChild(chip);
    });

    // Reset link — mirrors the one on the genus chart.
    const anyOff = props.some((p) => PROT_HIDDEN_PROPS.has(p));
    const reset = document.createElement("button");
    reset.type = "button";
    reset.className = "prot-cat-reset";
    reset.textContent = anyOff ? "show all" : "VF/AMR only";
    reset.title = anyOff ? "Show every annotation category" : "Show only Virulence Factor and AMR categories";
    reset.onclick = () => {
      if (anyOff) PROT_HIDDEN_PROPS.clear();
      else
        props.forEach((p) => {
          if (!protPropVisibleByDefault(p, rowFor.get(p))) PROT_HIDDEN_PROPS.add(p);
        });
      _syncProtCategoryViews();
    };
    host.appendChild(reset);
  }

  /* Single place that repaints everything keyed off PROT_HIDDEN_PROPS, so the
       table legend, the chart legend, both charts and the row list can never
       disagree about which categories are on. */
  function _syncProtCategoryViews() {
    _renderProtCatLegend();
    _filterProt();
    if (typeof drawProtGenus === "function") drawProtGenus();
    if (typeof drawProtProperty === "function") drawProtProperty();
  }
  // Let the chart legends refresh the chips after they toggle a category.
  window._renderProtCatLegend = _renderProtCatLegend;

  function _renderProtHeader() {
    const hr = document.getElementById("prot-header-row");
    hr.innerHTML = "";
    // Col 1: Category colour swatch
    const swTh = document.createElement("th");
    swTh.title = "Category colour";
    swTh.style.cssText =
      "background:#e8eaf6;padding:.28em .4em;white-space:nowrap;position:sticky;top:0;z-index:2;font-size:.8em;border:1px solid #e0e0e0;width:10px";
    hr.appendChild(swTh);
    // Col 2: TASS report indicator
    const tassTh = document.createElement("th");
    tassTh.title = "TASS report: whether this organism is detected in ANY sample in this run";
    tassTh.textContent = "TASS";
    tassTh.style.cssText =
      "background:#e8eaf6;padding:.28em .4em;white-space:nowrap;position:sticky;top:0;z-index:2;font-size:.75em;border:1px solid #e0e0e0;width:36px;text-align:center;cursor:help";
    hr.appendChild(tassTh);
    _protCols.forEach((c) => {
      const th = document.createElement("th");
      th.textContent = c;
      if (c === _protSortCol) th.classList.add(_protSortAsc ? "sort-asc" : "sort-desc");
      th.style.cssText =
        "background:#e8eaf6;padding:.28em .5em;white-space:nowrap;cursor:pointer;position:sticky;top:0;z-index:2;font-size:.8em;border:1px solid #e0e0e0";
      th.addEventListener("click", () => {
        if (_protSortCol === c) _protSortAsc = !_protSortAsc;
        else {
          _protSortCol = c;
          _protSortAsc = true;
        }
        _renderProtHeader();
        _filterProt();
      });
      hr.appendChild(th);
    });
  }

  function _filterProt() {
    const q = (document.getElementById("prot-search").value || "").toLowerCase();
    const col = document.getElementById("prot-search-col").value;
    const tassMode = (document.getElementById("prot-tass-filter") || {}).value || "all";
    const minPid = parseFloat((document.getElementById("prot-pid-min") || {}).value);
    // minPid is 0–100 scale; %id in data is also 0–100 → compare directly
    const pidThresh = !isNaN(minPid) ? minPid : 0;
    // Keyword group filter
    const kwGroup = (document.getElementById("prot-keyword-group") || {}).value || "";
    const kwTerms =
      kwGroup && PROT_KEYWORD_GROUPS[kwGroup] ? PROT_KEYWORD_GROUPS[kwGroup].map((t) => t.toLowerCase()) : null;

    // Use ALL data (not filtered) — organism is "in TASS" if found in ANY sample in this run.
    // Backed by a memoized prefix Set so this is O(1) per row instead of O(N).
    const _tassPrefixSet = _getTassPrefixSet();
    function _speciesInTass(species) {
      if (!species) return false;
      return _tassPrefixSet.has(species.trim().toLowerCase());
    }
    // Global Microbial Category selection — a hit inherits the category of the
    // detection it links to, so "Primary" alone hides Commensal/Unknown hits.
    const _catFilter = _protCategoryFilter();

    let rows = _protAllRows.filter((r) => {
      const rowProp = r["Property"] || r["Class"] || r["_source"] || "";
      if (rowProp && PROT_HIDDEN_PROPS.has(rowProp)) return false;
      if (!_protRowInCategories(r, _catFilter)) return false;
      // %id threshold (both data and pidThresh are 0–100)
      if (pidThresh > 0) {
        const pid = parseFloat(r["%id"] || r["pident"] || r["%ID"] || r["identity"] || 0);
        if (!isNaN(pid) && pid < pidThresh) return false;
      }
      // Keyword group filter — row must match at least one term in any checked field
      if (kwTerms) {
        const haystack = _KW_FIELDS.map((f) => String(r[f] || "").toLowerCase()).join(" ");
        if (!kwTerms.some((t) => haystack.includes(t))) return false;
      }
      // Text search
      if (q) {
        // _protCellText, not String(): object-valued fields stringify to
        // "[object object]" and would match the letter "j" but nothing useful.
        const match = col
          ? _protCellText(r[col]).toLowerCase().includes(q)
          : _protCols.some((c) => _protCellText(r[c]).toLowerCase().includes(q));
        if (!match) return false;
      }
      // TASS-presence filter (species-level)
      if (tassMode !== "all") {
        const rSpecies = (r["Species"] || "").trim();
        const inTass = _speciesInTass(rSpecies);
        if (tassMode === "in" && !inTass) return false;
        if (tassMode === "out" && inTass) return false;
      }
      // Bar chart click filter (sample × category), also used by "View VF/AMR"
      // jumps from Summary to pin the table to one sample.
      if (window._protBarFilter) {
        const bf = window._protBarFilter;
        const rSample = r["Specimen ID"] || r["Sample"] || r["sample"] || r["specimen_id"] || "";
        const rCat =
          r["Property"] ||
          r["property"] ||
          r["Category"] ||
          r["category"] ||
          r["Class"] ||
          r["class"] ||
          r["Type"] ||
          r["type"] ||
          r["Function"] ||
          r["function"] ||
          "";
        // `samples` carries the member sample ids of a merged specimen; fall
        // back to the single `sample` label when it is absent.
        const bfSamples = Array.isArray(bf.samples) && bf.samples.length ? bf.samples : bf.sample ? [bf.sample] : null;
        if (bfSamples && bfSamples.indexOf(rSample) === -1) return false;
        if (bf.cat && rCat !== bf.cat) return false;
      }
      return true;
    });
    if (_protSortCol) {
      rows = [...rows].sort((a, b) => {
        const va = a[_protSortCol] || "",
          vb = b[_protSortCol] || "";
        const na = parseFloat(va),
          nb = parseFloat(vb);
        const cmp = isNaN(na) || isNaN(nb) ? String(va).localeCompare(String(vb)) : na - nb;
        return _protSortAsc ? cmp : -cmp;
      });
    }
    // Hand the full result to the pager; it renders only the visible slice.
    _setProtPageRows(rows, true);
  }

  /* ── Pagination ──────────────────────────────────────────────────────────
         _setProtPageRows() is the single entry point: every filter/sort path calls
         it with the complete result set and it takes care of clamping the current
         page, rendering one slice and repainting the controls. */
  function _protPageCount() {
    if (!_protPageSize) return 1;
    return Math.max(1, Math.ceil(_protPageRows.length / _protPageSize));
  }
  function _setProtPageRows(rows, resetPage) {
    _protPageRows = rows || [];
    if (resetPage) _protPage = 1;
    _renderProtPage();
  }
  function _renderProtPage() {
    const total = _protPageRows.length;
    const pages = _protPageCount();
    // Clamp — the row count shrinks as filters tighten, so the active page can
    // fall off the end between renders.
    if (_protPage > pages) _protPage = pages;
    if (_protPage < 1) _protPage = 1;
    const start = _protPageSize ? (_protPage - 1) * _protPageSize : 0;
    const end = _protPageSize ? Math.min(start + _protPageSize, total) : total;
    _renderProtRows(_protPageRows.slice(start, end));

    const countEl = document.getElementById("prot-table-count");
    if (countEl) {
      countEl.textContent = total
        ? _protPageSize
          ? `${(start + 1).toLocaleString()}–${end.toLocaleString()} of ${total.toLocaleString()} rows`
          : `${total.toLocaleString()} rows`
        : "0 rows";
    }
    const status = document.getElementById("prot-page-status");
    if (status) status.textContent = total ? `Page ${_protPage.toLocaleString()} of ${pages.toLocaleString()}` : "—";
    const atFirst = _protPage <= 1;
    const atLast = _protPage >= pages;
    [
      ["prot-page-first", atFirst],
      ["prot-page-prev", atFirst],
      ["prot-page-next", atLast],
      ["prot-page-last", atLast],
    ].forEach(([id, disabled]) => {
      const b = document.getElementById(id);
      if (b) b.disabled = disabled || !total;
    });
    // Jumping pages leaves the scroll position deep in the previous page.
    const cont = document.getElementById("prot-table-container");
    if (cont) cont.scrollTop = 0;
  }
  function _gotoProtPage(p) {
    _protPage = p;
    _renderProtPage();
  }
  function _wireProtPager() {
    const on = (id, fn) => {
      const el = document.getElementById(id);
      // Replace rather than add: _buildProtTable can run more than once and
      // stacked listeners would advance the page several times per click.
      if (el) el.onclick = fn;
    };
    on("prot-page-first", () => _gotoProtPage(1));
    on("prot-page-prev", () => _gotoProtPage(_protPage - 1));
    on("prot-page-next", () => _gotoProtPage(_protPage + 1));
    on("prot-page-last", () => _gotoProtPage(_protPageCount()));
    const sizeEl = document.getElementById("prot-page-size");
    if (sizeEl) {
      sizeEl.onchange = () => {
        _protPageSize = parseInt(sizeEl.value, 10) || 0;
        _protPage = 1;
        _renderProtPage();
      };
    }
  }

  // Columns whose cells trigger gene-distribution mode on click
  const _GENE_CLICK_COLS = new Set([
    "Gene",
    "gene",
    "Gene Name",
    "gene_name",
    "Antibiotics", // renamed from Product via PROT_COL_REMAP
    "Product",
    "product",
    "annotation",
    "Annotation",
    "gene_id",
    "Description",
    "description",
  ]);

  function _renderProtRows(rows) {
    const tbody = document.getElementById("prot-table-body");
    tbody.innerHTML = "";
    // Build the TASS-prefix Set ONCE per render rather than per row.
    // Previously we ran DATA.map(...).filter(Boolean).some(...) inside
    // the per-row loop; with 14k DATA entries × thousands of protein
    // rows that was the dominant cost of opening the VF/AMR tab.
    const _tassPrefixSet = _getTassPrefixSet();
    // Use a DocumentFragment so we hit the live DOM exactly once at
    // the end instead of N times during the loop.
    const frag = document.createDocumentFragment();
    rows.forEach((r) => {
      const tr = document.createElement("tr");
      tr.style.cssText = "cursor:pointer;border-bottom:1px solid #f0f0f0";
      tr.addEventListener("mouseenter", () => (tr.style.background = "#e3f2fd"));
      tr.addEventListener("mouseleave", () => (tr.style.background = ""));
      // Row click (non-gene cells) → organism / category-breakdown mode
      tr.addEventListener("click", () => {
        const org = r["Organism"] || r["Species"] || r["Genus"] || r["gene"] || "";
        const genus = r["Genus"] || "";
        _showProtDetail(genus, org, null); // null → organism mode
      });
      // Col 1: Category colour swatch
      const swTd = document.createElement("td");
      const catLabel = r["Property"] || r["_source"] || "Other";
      swTd.title = catLabel;
      swTd.style.cssText = `background:${_catColor(r)};width:8px;min-width:8px;padding:0;border:1px solid #e0e0e0`;
      tr.appendChild(swTd);

      // Col 2: TASS report indicator — O(1) lookup via prefix set.
      const _rSpecies = (r["Species"] || "").trim().toLowerCase();
      const _rInTass = _rSpecies ? _tassPrefixSet.has(_rSpecies) : false;
      const tassTd = document.createElement("td");
      tassTd.style.cssText =
        "width:36px;text-align:center;padding:.2em .3em;border:1px solid #e0e0e0;white-space:nowrap;font-size:.72em;font-weight:600";
      if (_rInTass) {
        tassTd.textContent = "✓";
        tassTd.style.background = "#e8f5e9";
        tassTd.style.color = "#2e7d32";
        tassTd.title = "Organism detected in this run's TASS report";
      } else {
        tassTd.textContent = "Ext";
        tassTd.style.background = "#fff3e0";
        tassTd.style.color = "#e65100";
        tassTd.title = "Organism NOT found in any sample in this run's TASS report";
      }
      tr.appendChild(tassTd);

      _protCols.forEach((c) => {
        const td = document.createElement("td");
        const val = _protCellText(r[c]);
        td.textContent = val;

        if (_GENE_CLICK_COLS.has(c) && val) {
          // Gene / Product cell — clickable independently for gene-distribution popup
          td.style.cssText =
            "padding:.28em .5em;white-space:nowrap;font-size:.8em;border:1px solid #e0e0e0;max-width:220px;overflow:hidden;text-overflow:ellipsis;color:#1565c0;text-decoration:underline dotted;cursor:pointer";
          td.title = `See "${val}" distribution across all samples`;
          td.addEventListener("click", (e) => {
            e.stopPropagation(); // don't also trigger row (organism) click
            const org = r["Organism"] || r["Species"] || r["Genus"] || r["gene"] || "";
            const genus = r["Genus"] || "";
            _showProtDetail(genus, org, r); // pass row → gene-distribution mode
          });
        } else {
          td.style.cssText =
            "padding:.28em .5em;white-space:nowrap;font-size:.8em;border:1px solid #e0e0e0;max-width:220px;overflow:hidden;text-overflow:ellipsis";
        }
        tr.appendChild(td);
      });
      frag.appendChild(tr);
    });
    // Single DOM write — much faster than N append calls on the live tbody.
    tbody.appendChild(frag);
    // NOTE: the row count is written by _renderProtPage(), which knows the full
    // result size — `rows` here is only the current page.
  }

  // Expose so drawProteins can call it
  window._initProtTable = _buildProtTable;
  window._protBarFilter = null;
  window._filterProtExternal = function () {
    _filterProt();
  };
  window._clearProtBarFilter = function () {
    window._protBarFilter = null;
    window._protJumpSample = null;
    const badge = document.getElementById("prot-bar-filter-badge");
    if (badge) badge.style.display = "none";
    _filterProt();
  };
})();

function _showProtDetail(genus, organism, clickedRow) {
  const panel = document.getElementById("prot-detail-panel");
  const title = document.getElementById("prot-detail-title");
  const body = document.getElementById("prot-detail-body");
  if (!panel || !body) return;

  const label = organism || genus || "Unknown";
  // Title will be refined below once we know if a gene was clicked
  title.textContent = label;

  // Dynamically position panel below banner + tabbar
  const _bannerH = (document.getElementById("banner") || { offsetHeight: 60 }).offsetHeight;
  const _tabH = (document.getElementById("tabbar") || { offsetHeight: 44 }).offsetHeight;
  panel.style.top = _bannerH + _tabH + "px";

  panel.classList.add("open");
  // Store for resize-triggered redraw
  window._lastProtDetail = { g: genus, o: organism, r: clickedRow, gh: null };

  // Find matching DATA record(s) for TASS / coverage info
  const dataMatches = DATA.filter((r) => {
    const org = (r["Detected Organism"] || "").toLowerCase();
    const gen = (r["Genus"] || "").toLowerCase();
    if (organism && org.toLowerCase().includes(organism.toLowerCase())) return true;
    if (genus && gen.toLowerCase().includes(genus.toLowerCase())) return true;
    return false;
  });

  // Gene hits for this organism/genus (deduplicated)
  const allGeneRows = _dedupRows([...(PROT.per_gene_hits || []), ...(PROT.amr_genes || [])]);
  const _GENE_FIELDS_DET = [
    "Gene",
    "gene",
    "Gene Name",
    "gene_name",
    "Product",
    "product",
    "annotation",
    "Annotation",
    "gene_id",
    "Description",
    "description",
  ];
  let geneHits = allGeneRows.filter((r) => {
    const rOrg = (r["Organism"] || r["Species"] || "").toLowerCase();
    const rGen = (r["Genus"] || "").toLowerCase();
    if (organism && rOrg.includes(organism.toLowerCase())) return true;
    if (genus && rGen.includes(genus.toLowerCase())) return true;
    return false;
  });

  // When a specific gene/product cell was clicked, restrict the table to that gene only
  let _clickedGeneName = null;
  if (clickedRow) {
    _clickedGeneName = _GENE_FIELDS_DET.map((f) => clickedRow[f]).find((v) => v != null && v !== "") || null;
    if (_clickedGeneName) {
      const _cgl = _clickedGeneName.toLowerCase();
      geneHits = geneHits.filter((r) =>
        _GENE_FIELDS_DET.some((f) => r[f] != null && String(r[f]).toLowerCase() === _cgl),
      );
    }
  }

  // Refine panel title now that we know whether a gene was clicked
  if (_clickedGeneName) {
    title.textContent = `${_clickedGeneName}  ·  ${label}`;
  }

  let html = "";

  // 1) Gene categories bar chart placeholder — rendered first (top)
  html += `<div id="prot-compare-wrap">
    <div class="compare-title" id="prot-compare-title"></div>
    <div id="prot-compare-svg"></div>
  </div>`;

  // 2) Gene hits table — middle
  if (geneHits.length) {
    const geneCols = geneHits.length
      ? [...new Set(geneHits.flatMap((r) => Object.keys(r)).filter((k) => !k.startsWith("_")))]
      : [];
    html += `<div style="margin-top:.6em">
      <div style="font-size:.78em;font-weight:700;color:#1565c0;text-transform:uppercase;letter-spacing:.05em;margin-bottom:.35em">${
        _clickedGeneName
          ? `"${_clickedGeneName}" — all hits (${geneHits.length})`
          : `Associated Genes / Annotations (${geneHits.length})`
      }</div>
      <div style="max-height:280px;overflow:auto">
      <table style="border-collapse:collapse;width:100%;font-size:.76em">`;
    html += `<tr>${geneCols
      .map(
        (c) =>
          `<th style="background:#e8eaf6;padding:.22em .4em;text-align:left;border:1px solid #ddd;white-space:nowrap;position:sticky;top:0">${c}</th>`,
      )
      .join("")}</tr>`;
    geneHits.forEach((r) => {
      html += `<tr>${geneCols
        .map(
          (c) =>
            `<td style="padding:.22em .4em;border:1px solid #eee;white-space:nowrap;max-width:180px;overflow:hidden;text-overflow:ellipsis">${_protCellText(
              r[c],
            )}</td>`,
        )
        .join("")}</tr>`;
    });
    html += `</table></div></div>`;
  }

  // 3) TASS Detection Summary — bottom
  if (dataMatches.length) {
    html += `<div style="margin-top:.8em">
      <div style="font-size:.78em;font-weight:700;color:#1565c0;text-transform:uppercase;letter-spacing:.05em;margin-bottom:.35em">Detection Summary (${
        dataMatches.length
      } sample${dataMatches.length > 1 ? "s" : ""})</div>`;
    const statCols = [
      "Specimen ID",
      "TASS Score",
      "# Reads Aligned",
      "Breadth %",
      "Mean Depth",
      "Coverage",
      "Passes Threshold",
      "High Consequence",
      "Microbial Category",
      "Genus",
      "Phylum",
    ];
    html += `<table style="border-collapse:collapse;width:100%;font-size:.78em">`;
    html += `<tr>${statCols
      .map(
        (c) =>
          `<th style="background:#e8eaf6;padding:.22em .4em;text-align:left;border:1px solid #ddd;white-space:nowrap">${c}</th>`,
      )
      .join("")}</tr>`;
    dataMatches.forEach((r) => {
      html += `<tr>${statCols
        .map(
          (c) => `<td style="padding:.22em .4em;border:1px solid #eee;white-space:nowrap">${_protCellText(r[c])}</td>`,
        )
        .join("")}</tr>`;
    });
    html += `</table></div>`;
  }

  if (!dataMatches.length && !geneHits.length) {
    html = `<div id="prot-compare-wrap">
    <div class="compare-title" id="prot-compare-title"></div>
    <div id="prot-compare-svg"></div>
  </div><p style="color:#999;font-style:italic">No detailed data found for "${label}".</p>`;
  }

  body.innerHTML = html;

  // ── Search bar wiring ──────────────────────────────────────────────────
  const _searchIn = document.getElementById("prot-detail-search");
  const _searchClr = document.getElementById("prot-detail-search-clear");
  if (_searchIn) {
    _searchIn.value = ""; // reset on each open
    if (_searchClr) _searchClr.style.display = "none";

    const _filterDetailTable = () => {
      const q = _searchIn.value.trim().toLowerCase();
      // Collect all data <tr>s inside prot-detail-body tables (skip header rows)
      const rows = body.querySelectorAll("table tr:not(:first-child)");
      rows.forEach((row) => {
        if (!q) {
          row.style.display = "";
          return;
        }
        const text = Array.from(row.querySelectorAll("td"))
          .map((td) => td.textContent)
          .join(" ")
          .toLowerCase();
        row.style.display = text.includes(q) ? "" : "none";
      });
      if (_searchClr) _searchClr.style.display = q ? "" : "none";
    };

    _searchIn.oninput = _filterDetailTable;
    if (_searchClr) {
      _searchClr.onclick = () => {
        _searchIn.value = "";
        _filterDetailTable();
        _searchIn.focus();
      };
    }
  }

  // Store gene hits for resize redraw
  if (window._lastProtDetail) window._lastProtDetail.gh = geneHits;
  // Draw the comparison chart after DOM is updated
  _drawProtCompare(genus, organism, clickedRow, geneHits);
}

/* ─── Comparison chart for VF/AMR detail panel ───────────────────────────── */
function _drawProtCompare(genus, organism, clickedRow, geneHits) {
  const wrap = document.getElementById("prot-compare-svg");
  const title = document.getElementById("prot-compare-title");
  if (!wrap) return;
  wrap.innerHTML = "";

  const allGeneRows = _dedupRows([...(PROT.per_gene_hits || []), ...(PROT.amr_genes || [])]);

  // Detect gene name field
  const GENE_FIELDS = [
    "Gene",
    "gene",
    "Gene Name",
    "gene_name",
    "Product",
    "product",
    "annotation",
    "Annotation",
    "gene_id",
  ];
  const CAT_FIELDS = [
    "Property",
    "property",
    "Category",
    "category",
    "Class",
    "class",
    "Type",
    "type",
    "Function",
    "function",
  ];

  function _getField(row, candidates) {
    for (const f of candidates) if (row[f] !== undefined && row[f] !== null && row[f] !== "") return row[f];
    return null;
  }
  function _getSample(row) {
    return row["Specimen ID"] || row["Sample"] || row["sample"] || row["specimen_id"] || "Unknown";
  }
  function _getOrg(row) {
    return row["Organism"] || row["Species"] || row["organism"] || row["species"] || row["Detected Organism"] || "";
  }
  function _getGenus(row) {
    return row["Genus"] || row["genus"] || "";
  }

  const panelW = document.getElementById("prot-detail-panel").offsetWidth;
  const W = Math.max(280, panelW - 32);

  if (clickedRow) {
    // ── GENE MODE: compare this gene across all samples ───────────────────
    const clickedGene = _getField(clickedRow, GENE_FIELDS);
    if (!clickedGene) return;
    title.textContent = `"${clickedGene}" across all samples`;

    // Collect all rows matching this gene name
    const matches = allGeneRows.filter((r) => {
      const g = _getField(r, GENE_FIELDS);
      return g && g.toLowerCase() === clickedGene.toLowerCase();
    });

    if (!matches.length) {
      wrap.innerHTML = `<p style="color:#999;font-size:.78em;font-style:italic">No cross-sample data found for "${clickedGene}".</p>`;
      return;
    }

    // Group by sample
    const byS = {};
    matches.forEach((r) => {
      const s = _getSample(r);
      byS[s] = (byS[s] || 0) + 1;
    });
    const data = Object.entries(byS).sort((a, b) => b[1] - a[1]);

    const _maxLblGene = data.reduce((m, d) => Math.max(m, d[0].length), 0);
    const mL = Math.max(80, Math.min(200, _maxLblGene * 7 + 10));
    const mT = 8,
      mR = 16,
      mB = 40;
    const H = data.length * 24 + mT + mB;
    const iW = W - mL - mR;
    const xMax = d3.max(data, (d) => d[1]) || 1;

    const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H).style("overflow", "visible");
    const g = svg.append("g").attr("transform", `translate(${mL},${mT})`);
    const y = d3
      .scaleBand()
      .domain(data.map((d) => d[0]))
      .range([0, H - mT - mB])
      .padding(0.2);
    const x = d3.scaleLinear().domain([0, xMax]).range([0, iW]).nice();

    g.append("g").call(d3.axisLeft(y).tickSize(0)).select(".domain").remove();
    g.selectAll("text").style("font-size", "11px");

    g.append("g")
      .attr("transform", `translate(0,${H - mT - mB})`)
      .call(d3.axisBottom(x).ticks(Math.min(5, xMax)).tickFormat(d3.format("d")));

    g.selectAll("rect")
      .data(data)
      .join("rect")
      .attr("y", (d) => y(d[0]))
      .attr("height", y.bandwidth())
      .attr("x", 0)
      .attr("width", (d) => x(d[1]))
      .attr("fill", "#1565c0")
      .attr("rx", 2)
      .on("mouseover", (ev, d) => showTip(`<b>${d[0]}</b>: ${d[1]} hit${d[1] > 1 ? "s" : ""}`, ev))
      .on("mousemove", moveTip)
      .on("mouseout", hideTip);

    g.selectAll(".val-label")
      .data(data)
      .join("text")
      .attr("class", "val-label")
      .attr("x", (d) => x(d[1]) + 3)
      .attr("y", (d) => y(d[0]) + y.bandwidth() / 2 + 4)
      .style("font-size", "10px")
      .attr("fill", "#333")
      .text((d) => d[1]);

    svg
      .append("text")
      .attr("x", mL + iW / 2)
      .attr("y", H - 8)
      .attr("text-anchor", "middle")
      .style("font-size", "10px")
      .attr("fill", "#666")
      .text("Hit count");
  } else {
    // ── ORGANISM MODE: gene categories across organisms/samples ──────────
    const isGenus = !organism;
    title.textContent = isGenus
      ? `Gene categories for genus: ${genus}`
      : `Gene categories across samples: ${organism || genus}`;

    // Filter to this organism/genus
    const filtered = allGeneRows.filter((r) => {
      const rOrg = _getOrg(r).toLowerCase();
      const rGen = _getGenus(r).toLowerCase();
      if (organism && rOrg.includes(organism.toLowerCase())) return true;
      if (genus && rGen.includes(genus.toLowerCase())) return true;
      return false;
    });

    if (!filtered.length) {
      wrap.innerHTML = `<p style="color:#999;font-size:.78em;font-style:italic">No gene category data found.</p>`;
      return;
    }

    // Group by sample × category
    const samples = _orderedSamples([...new Set(filtered.map(_getSample))]);
    const cats = [...new Set(filtered.map((r) => _getField(r, CAT_FIELDS) || "Other"))].sort();
    const palette = d3.schemeTableau10;

    // count[sample][cat] = n
    const count = {};
    samples.forEach((s) => {
      count[s] = {};
      cats.forEach((c) => {
        count[s][c] = 0;
      });
    });
    filtered.forEach((r) => {
      const s = _getSample(r);
      const c = _getField(r, CAT_FIELDS) || "Other";
      if (count[s]) count[s][c] = (count[s][c] || 0) + 1;
    });

    // Build stacked data
    const stackData = samples.map((s) => {
      const obj = { sample: s };
      cats.forEach((c) => {
        obj[c] = count[s][c] || 0;
      });
      return obj;
    });
    const stack = d3.stack().keys(cats)(stackData);
    const totals = samples.map((s) => cats.reduce((t, c) => t + (count[s][c] || 0), 0));
    const xMax = d3.max(totals) || 1;

    // Dynamic left margin based on longest sample label
    const _maxLblOrg = samples.reduce((m, s) => Math.max(m, s.length), 0);
    const mL = Math.max(80, Math.min(200, _maxLblOrg * 7 + 10));
    // Legend in 2 columns below the chart
    const legendCols = 2;
    const legendRows = Math.ceil(cats.length / legendCols);
    const legendH = legendRows * 16 + 8;
    const mT = 8,
      mR = 16,
      mB = 30 + legendH;
    const chartH = samples.length * 28;
    const H = chartH + mT + mB;
    const iW = W - mL - mR;

    const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H).style("overflow", "visible");
    const g = svg.append("g").attr("transform", `translate(${mL},${mT})`);
    const y = d3.scaleBand().domain(samples).range([0, chartH]).padding(0.2);
    const x = d3.scaleLinear().domain([0, xMax]).range([0, iW]).nice();
    const col = d3.scaleOrdinal().domain(cats).range(palette);

    g.append("g").call(d3.axisLeft(y).tickSize(0)).select(".domain").remove();
    g.selectAll("text").style("font-size", "11px");
    g.append("g")
      .attr("transform", `translate(0,${chartH})`)
      .call(d3.axisBottom(x).ticks(Math.min(5, xMax)).tickFormat(d3.format("d")));

    stack.forEach((layer) => {
      g.selectAll(`.bar-${layer.key.replace(/\W/g, "_")}`)
        .data(layer)
        .join("rect")
        .attr("class", `bar-${layer.key.replace(/\W/g, "_")}`)
        .attr("y", (d) => y(d.data.sample))
        .attr("height", y.bandwidth())
        .attr("x", (d) => x(d[0]))
        .attr("width", (d) => x(d[1]) - x(d[0]))
        .attr("fill", col(layer.key))
        .style("cursor", "pointer")
        .on("mouseover", (ev, d) =>
          showTip(
            `<b>${d.data.sample}</b><br>${layer.key}: ${d[1] - d[0]}<br><small>Click to filter table</small>`,
            ev,
          ),
        )
        .on("mousemove", moveTip)
        .on("mouseout", hideTip)
        .on("click", (ev, d) => {
          ev.stopPropagation();
          window._protBarFilter = { sample: d.data.sample, cat: layer.key };
          const badge = document.getElementById("prot-bar-filter-badge");
          const badgeText = document.getElementById("prot-bar-filter-text");
          if (badge && badgeText) {
            badgeText.textContent = `${d.data.sample} \u00b7 ${layer.key}`;
            badge.style.display = "inline-flex";
          }
          if (window._filterProtExternal) window._filterProtExternal();
        });
    });

    // x-axis label
    svg
      .append("text")
      .attr("x", mL + iW / 2)
      .attr("y", mT + chartH + 28)
      .attr("text-anchor", "middle")
      .style("font-size", "10px")
      .attr("fill", "#666")
      .text("Gene hit count");

    // Legend — below chart in 2 columns
    const legY0 = mT + chartH + 38;
    const colW = iW / legendCols;
    cats.forEach((c, i) => {
      const col_ = i % legendCols,
        row_ = Math.floor(i / legendCols);
      const lx = mL + col_ * colW;
      const ly = legY0 + row_ * 16;
      svg
        .append("rect")
        .attr("x", lx)
        .attr("y", ly)
        .attr("width", 10)
        .attr("height", 10)
        .attr("fill", col(c))
        .attr("rx", 2);
      svg
        .append("text")
        .attr("x", lx + 13)
        .attr("y", ly + 9)
        .style("font-size", "9px")
        .attr("fill", "#333")
        .text(c.length > 18 ? c.slice(0, 17) + "…" : c);
    });
  }
}
