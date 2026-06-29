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
  const hasByRow = (PROT.per_gene_hits || []).length > 0 || (PROT.amr_genes || []).length > 0;
  const agg = {};
  const propSet = new Set();

  if (hasByRow) {
    const _pgFdPairs = new Set(filteredData().map((r) => `${r["Specimen ID"]}||${r["Genus"] || "Unknown"}`));
    const visibleHits = [
      ...(PROT.per_gene_hits || []).filter(
        (r) => !sampleHidden[r["Specimen ID"]] && _pgFdPairs.has(`${r["Specimen ID"]}||${r["Genus"] || "Unknown"}`),
      ),
      ...(PROT.amr_genes || [])
        .filter(
          (r) => !sampleHidden[r["Specimen ID"]] && _pgFdPairs.has(`${r["Specimen ID"]}||${r["Genus"] || "Unknown"}`),
        )
        .map((r) => ({ ...r, Property: r["Property"] || r["Class"] })),
    ];
    visibleHits.forEach((r) => {
      const genus = r["Genus"] || "Unknown";
      const prop = r["Property"] || r["Class"] || "Unknown";
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
  }

  if (!Object.keys(agg).length) {
    wrap.innerHTML = '<p style="color:#999">No genus annotation data.</p>';
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
        if (window._filterProtExternal) window._filterProtExternal();
      })
      .on("mouseover", (ev) => {
        const action = PROT_HIDDEN_PROPS.has(p) ? "Enable" : "Disable";
        showTip(`${action} ${p}<br><small>Click to ${action.toLowerCase()} in chart and table</small>`, ev);
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

  const _ppFdPairs = new Set(filteredData().map((r) => `${r["Specimen ID"]}||${r["Genus"] || "Unknown"}`));
  const _allHits = [
    ...(PROT.per_gene_hits || []).filter(
      (r) => !sampleHidden[r["Specimen ID"]] && _ppFdPairs.has(`${r["Specimen ID"]}||${r["Genus"] || "Unknown"}`),
    ),
    ...(PROT.amr_genes || [])
      .filter(
        (r) => !sampleHidden[r["Specimen ID"]] && _ppFdPairs.has(`${r["Specimen ID"]}||${r["Genus"] || "Unknown"}`),
      )
      .map((r) => ({ ...r, _source: "AMR" })),
  ];
  const _propGenus = {};
  _allHits.forEach((r) => {
    const prop = r["Property"] || r["Class"] || r["_source"] || "Other";
    const gen = r["Genus"] || "Unknown";
    if (!_propGenus[prop]) _propGenus[prop] = {};
    _propGenus[prop][gen] = (_propGenus[prop][gen] || 0) + 1;
  });

  // Always build from per-row hits for accurate sample-visibility filtering.
  // Only use pre-aggregated metadata_counts when no per-row data exists.
  const hasByRow = (PROT.per_gene_hits || []).length > 0 || (PROT.amr_genes || []).length > 0;
  let propRows;
  if (hasByRow) {
    const agg = {};
    _allHits.forEach((r) => {
      const prop = r["Property"] || r["Class"] || r["_source"] || "Other";
      agg[prop] = (agg[prop] || 0) + 1;
    });
    propRows = Object.entries(agg).map(([value, count]) => ({ field: "property", value, count }));
  } else {
    const counts = PROT.metadata_counts || [];
    propRows = counts.filter((r) => r["field"] === "property");
  }
  if (!propRows.length) {
    wrap.innerHTML = '<p style="color:#999">No property metadata.</p>';
    return;
  }

  const data = propRows
    .map((r) => ({ label: r["value"] || "", value: parseInt(r["count"]) || 0 }))
    .sort((a, b) => b.value - a.value)
    .slice(0, 20);

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

  // Category color map (Property field)
  const _catColors = {
    "Virulence Factor": "#e53935",
    "Antibiotic Resistance": "#fb8c00",
    AMR: "#fb8c00",
    "Drug Target": "#8e24aa",
    Transporter: "#00897b",
  };
  function _catColor(row) {
    const prop = row["Property"] || row["_source"] || "";
    return _catColors[prop] || "#90a4ae";
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
    _protAllRows = _dedupRows([...geneRows, ...amrRows]);
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
    _renderProtRows(_protAllRows);

    document.getElementById("prot-search").addEventListener("input", _filterProt);
    colSel.addEventListener("change", _filterProt);
    const tassFilter = document.getElementById("prot-tass-filter");
    if (tassFilter) tassFilter.addEventListener("change", _filterProt);
    const pidMinEl = document.getElementById("prot-pid-min");
    if (pidMinEl) pidMinEl.addEventListener("input", _filterProt);
    const kwGroupEl = document.getElementById("prot-keyword-group");
    if (kwGroupEl) kwGroupEl.addEventListener("change", _filterProt);
  }

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

    let rows = _protAllRows.filter((r) => {
      const rowProp = r["Property"] || r["Class"] || r["_source"] || "";
      if (rowProp && PROT_HIDDEN_PROPS.has(rowProp)) return false;
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
        const match = col
          ? String(r[col] || "")
              .toLowerCase()
              .includes(q)
          : _protCols.some((c) =>
              String(r[c] || "")
                .toLowerCase()
                .includes(q),
            );
        if (!match) return false;
      }
      // TASS-presence filter (species-level)
      if (tassMode !== "all") {
        const rSpecies = (r["Species"] || "").trim();
        const inTass = _speciesInTass(rSpecies);
        if (tassMode === "in" && !inTass) return false;
        if (tassMode === "out" && inTass) return false;
      }
      // Bar chart click filter (sample × category)
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
        if (bf.sample && rSample !== bf.sample) return false;
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
    document.getElementById("prot-table-count").textContent = `${rows.length} rows`;
    _renderProtRows(rows);
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
        const val = r[c] !== undefined ? r[c] : "";
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
    document.getElementById("prot-table-count").textContent = `${rows.length} rows`;
  }

  // Expose so drawProteins can call it
  window._initProtTable = _buildProtTable;
  window._protBarFilter = null;
  window._filterProtExternal = function () {
    _filterProt();
  };
  window._clearProtBarFilter = function () {
    window._protBarFilter = null;
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
            `<td style="padding:.22em .4em;border:1px solid #eee;white-space:nowrap;max-width:180px;overflow:hidden;text-overflow:ellipsis">${
              r[c] !== undefined ? r[c] : ""
            }</td>`,
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
          (c) =>
            `<td style="padding:.22em .4em;border:1px solid #eee;white-space:nowrap">${
              r[c] !== undefined ? r[c] : ""
            }</td>`,
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
