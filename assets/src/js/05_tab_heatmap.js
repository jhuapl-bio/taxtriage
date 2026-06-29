/* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB: HEATMAP            (data-tab="heatmap"  —  the default tab)
       -     drawHeatmap()  →  D3 heatmap grid. Rows are organisms grouped /
       -     labelled by the chosen taxonomy rank, columns are samples.
       -     Reads filter UI: rank, value column, color scale, search/min TASS.
       -     Performance: uses an orgRepRow Map to avoid O(N²) fd.find() scans
       -     during sort + label rendering.
═══════════════════════════════════════════════════════════════════════════ */
function drawHeatmap() {
  const wrap = document.getElementById("heatmap-svg-wrap");
  wrap.innerHTML = "";
  const fd = filteredData();
  const col = document.getElementById("hm-value-col").value || "# Reads Aligned";
  const rank =
    (document.getElementById("hm-rank") || document.getElementById("tass-rank") || { value: "Genus" }).value || "Genus";
  const showVals = document.getElementById("hm-show-vals").checked;
  const txtColor = document.getElementById("hm-cell-color").value;

  const samples = _orderedSamples(uniq(fd.map((r) => r["Specimen ID"])).filter(Boolean));

  // Row label = organism name; sort rows by taxonomy rank then organism name
  const rankField = rank === "Species" ? "Detected Organism" : rank;

  // Pre-build organism → representative-row map in a single pass over fd.
  // Previously the sort comparator and the label loops below each called
  // fd.find() per organism — with 14k organisms that's O(N²) ≈ 200M ops
  // and was the dominant cost of drawHeatmap. With the map this becomes
  // a single linear pass + O(1) lookups.
  const orgRepRow = new Map();
  for (let i = 0; i < fd.length; i++) {
    const o = fd[i]["Detected Organism"];
    if (o && !orgRepRow.has(o)) orgRepRow.set(o, fd[i]);
  }
  const organisms = Array.from(orgRepRow.keys()).sort((a, b) => {
    const ga = (orgRepRow.get(a)[rankField] || "") + "";
    const gb = (orgRepRow.get(b)[rankField] || "") + "";
    return ga.localeCompare(gb) || a.localeCompare(b);
  });

  if (!samples.length || !organisms.length) {
    wrap.innerHTML = '<p style="color:#999;padding:1em">No data to display.</p>';
    return;
  }

  // With horizontal scrolling enabled, allow a comfortable minimum cell width
  // so the chart doesn't get cramped — the container will scroll if needed.
  const cellW = Math.max(44, Math.min(80, Math.floor(((wrap.clientWidth || 900) * 0.85) / samples.length)));
  const cellH = 22;
  // Dynamic left margin — must fit organism name AND rank label side by side.
  // Rank label sits at x=-6 (right edge); organism name sits to its left with a
  // small gap. Compute the max combined width across all organisms.
  let _maxLabelWidth = 0;
  organisms.forEach((org) => {
    const r0 = orgRepRow.get(org);
    const rl = (r0 ? r0[rankField] : "") || "";
    const orgTrunc = Math.min(32, org.length);
    const rankW = rankField !== "Detected Organism" && rl ? rl.length * 6 + 14 : 6;
    _maxLabelWidth = Math.max(_maxLabelWidth, rankW + orgTrunc * 6.5);
  });
  const marginL = Math.max(140, _maxLabelWidth + 20);
  const marginT = 110,
    marginR = 90; // wider: 28px for mol-type badge column + original 60 for colorbar
  // Pre-compute rescued status so we can add bottom margin for the legend
  const _anyRescued = fd.some((r) => {
    const _pi = rowPassInfo(r);
    return _pi && _pi.rescued;
  });
  const marginB = 20;
  const W = marginL + samples.length * cellW + marginR;
  const H = marginT + organisms.length * cellH + marginB;

  const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H);

  const vals = fd.map((r) => num(r[col])).filter((v) => v > 0);
  const vmax = vals.length ? d3.max(vals) : 1;
  const scaleMode = (document.getElementById("hm-scale") || { value: "linear" }).value;
  let color;
  if (scaleMode === "log") {
    color = d3.scaleSequentialLog(d3.interpolateBlues).domain([Math.max(1, d3.min(vals) || 1), Math.max(1, vmax)]);
  } else if (scaleMode === "sqrt") {
    color = d3.scaleSequentialSqrt(d3.interpolateBlues).domain([0, vmax]);
  } else {
    color = d3.scaleSequential(d3.interpolateBlues).domain([0, vmax]);
  }

  const xScale = d3
    .scaleBand()
    .domain(samples)
    .range([marginL, marginL + samples.length * cellW])
    .padding(0.03);
  const yScale = d3
    .scaleBand()
    .domain(organisms)
    .range([marginT, marginT + organisms.length * cellH])
    .padding(0.03);

  // X axis (top)
  svg
    .append("g")
    .attr("class", "axis")
    .attr("transform", `translate(0,${marginT})`)
    .call(d3.axisTop(xScale).tickSize(0))
    .selectAll("text")
    .attr("transform", "rotate(-45)")
    .style("text-anchor", "start")
    .attr("dx", "0.4em")
    .attr("dy", "-0.2em");

  // Y axis — organism labels with rank prefix
  const yGroup = svg.append("g").attr("class", "axis").attr("transform", `translate(${marginL},0)`);
  // Build high-consequence organism set and mol_type map from filtered data
  const _hcOrgSet = new Set(fd.filter((r) => r["High Consequence"]).map((r) => r["Detected Organism"]));
  const _orgMolType = {};
  fd.forEach((r) => {
    if (r["Mol Type"]) _orgMolType[r["Detected Organism"]] = r["Mol Type"];
  });

  // Right-side mol-type badge group (circle badges after last cell)
  const mtBadgeGroup = svg.append("g").attr("class", "mt-badges");
  const _mtBadgeX = marginL + samples.length * cellW + 10; // 10px gap after last cell

  organisms.forEach((org) => {
    const r0 = orgRepRow.get(org);
    const rankLabel = (r0 ? r0[rankField] : "") || "";
    const isHC = _hcOrgSet.has(org);
    const mt = (_orgMolType[org] || "").toLowerCase();
    const y = yScale(org) + yScale.bandwidth() / 2 + 4;
    // rank label in muted colour
    if (rankField !== "Detected Organism" && rankLabel) {
      yGroup
        .append("text")
        .attr("x", -6)
        .attr("y", y)
        .attr("text-anchor", "end")
        .attr("font-size", 9)
        .attr("fill", "#aaa")
        .text(`[${rankLabel}]`);
    }
    const _nameX = rankField !== "Detected Organism" && rankLabel ? -(rankLabel.length * 6 + 14) : -6;
    const orgLabel = org.length > 32 ? org.slice(0, 31) + "…" : org;
    const isWatched = r0 ? _isWatched(r0) : false;
    const prefix = (isWatched ? "★ " : "") + (isHC ? "● " : "");
    const titleParts = [
      isWatched ? "On the follow-up list" : null,
      isHC ? "High Consequence Pathogen" : null,
      mt ? `Mol. type: ${mt.toUpperCase()}` : null,
      "Click to toggle follow-up",
    ].filter(Boolean);
    yGroup
      .append("text")
      .attr("x", _nameX)
      .attr("y", y)
      .attr("text-anchor", "end")
      .attr("font-size", 11)
      .attr("fill", isHC ? "#c62828" : isWatched ? "#b8860b" : "#333")
      .attr("font-weight", isHC || isWatched ? "600" : "normal")
      .style("cursor", "pointer")
      .attr("title", titleParts.length ? titleParts.join(" | ") : null)
      .text(prefix + orgLabel)
      .on("click", () => {
        if (r0) _toggleWatchKey(_watchKey(r0));
      });

    // Draw circle badge on the right if mol_type is set
    if (mt === "dna" || mt === "rna") {
      const cy = yScale(org) + yScale.bandwidth() / 2;
      const badgeColor = mt === "dna" ? "#1565c0" : "#6a1b9a";
      const badgeLetter = mt === "dna" ? "D" : "R";
      const badgeG = mtBadgeGroup
        .append("g")
        .attr("transform", `translate(${_mtBadgeX + 7}, ${cy})`)
        .attr("title", mt === "dna" ? "DNA pathogen" : "RNA pathogen");
      badgeG.append("circle").attr("r", 7).attr("fill", badgeColor).attr("opacity", 0.9);
      badgeG
        .append("text")
        .attr("text-anchor", "middle")
        .attr("dominant-baseline", "central")
        .attr("font-size", 7)
        .attr("font-weight", "bold")
        .attr("fill", "#fff")
        .attr("pointer-events", "none")
        .text(badgeLetter);
    }
  });

  // Build lookup
  const lookup = {};
  fd.forEach((r) => {
    const k = `${r["Specimen ID"]}|||${r["Detected Organism"]}`;
    if (!lookup[k] || num(r[col]) > num(lookup[k][col])) lookup[k] = r;
  });

  organisms.forEach((org) => {
    samples.forEach((sp) => {
      const r = lookup[`${sp}|||${org}`];
      const v = r ? num(r[col]) : null;
      const cell = svg
        .append("rect")
        .attr("x", xScale(sp))
        .attr("y", yScale(org))
        .attr("width", xScale.bandwidth())
        .attr("height", yScale.bandwidth())
        .attr("fill", v !== null && v > 0 ? color(v) : "#f0f0f0")
        .attr("rx", 2);

      cell
        .on("mouseover", (ev) => {
          if (!r) {
            showTip("Not detected", ev);
            return;
          }
          const tass = num(r["TASS Score"]);
          const hc = isTruthy(r["High Consequence"]) ? " 🔴 HC" : "";
          showTip(
            `<b>${org}</b>${hc}<br>Sample: ${sp}<br>${col}: <b>${
              v !== null ? v.toLocaleString() : "n/a"
            }</b><br>TASS: ${tass.toFixed(1)} · ${rankLabel(r, rank)}`,
            ev,
          );
        })
        .on("mousemove", moveTip)
        .on("mouseout", hideTip);

      if (showVals && v !== null && v > 0) {
        svg
          .append("text")
          .attr("x", xScale(sp) + xScale.bandwidth() / 2)
          .attr("y", yScale(org) + yScale.bandwidth() / 2 + 4)
          .attr("text-anchor", "middle")
          .attr("fill", txtColor)
          .attr("font-size", 9)
          .text(v >= 1000 ? `${(v / 1000).toFixed(1)}k` : v.toFixed(0));
      }

      // Rescued-strain marker: top-right corner triangle when the strain
      // is below cutoff but its species/genus rollup passes. Orange =
      // species rescue, purple = genus rescue.
      if (r) {
        const _pi = rowPassInfo(r);
        if (_pi.rescued) {
          const x0 = xScale(sp),
            y0 = yScale(org),
            bw = xScale.bandwidth(),
            s = Math.min(8, bw * 0.32);
          const _rl = _pi.rescueLevel === "species" ? "species" : "genus";
          const _tassNum = +r["TASS Score"];
          const _tassTxt = isFinite(_tassNum) ? _tassNum.toFixed(1) : "n/a";
          svg
            .append("path")
            .attr("d", `M${x0 + bw - s},${y0} L${x0 + bw},${y0} L${x0 + bw},${y0 + s} Z`)
            .attr("fill", _rl === "species" ? "#fb8c00" : "#8e24aa")
            .style("cursor", "help")
            .append("title")
            .text(
              `Rescued strain: this strain's TASS (${_tassTxt}) is below the cutoff, ` +
                `but its ${_rl}-level roll-up passes, so it is kept visible. ` +
                `${_rl === "species" ? "Orange" : "Purple"} corner = ${_rl} rescue. ` +
                `Turn off "Roll up threshold" to hide it.`,
            );
        }
      }
    });
  });

  // Colour bar
  const cbW = 12,
    cbH = Math.min(160, H - marginT - marginB);
  const cbX = W - marginR + 10 + 30, // +30 to clear the mol-type badge column (28px wide)
    cbY = marginT;
  const cbGrad = svg
    .append("defs")
    .append("linearGradient")
    .attr("id", "cb-grad")
    .attr("x1", "0%")
    .attr("y1", "100%")
    .attr("x2", "0%")
    .attr("y2", "0%");
  d3.range(0, 1.01, 0.1).forEach((t) => {
    cbGrad
      .append("stop")
      .attr("offset", `${t * 100}%`)
      .attr("stop-color", color(t * vmax));
  });
  svg
    .append("rect")
    .attr("x", cbX)
    .attr("y", cbY)
    .attr("width", cbW)
    .attr("height", cbH)
    .style("fill", "url(#cb-grad)")
    .attr("rx", 2);
  svg
    .append("text")
    .attr("x", cbX + cbW / 2)
    .attr("y", cbY - 4)
    .attr("text-anchor", "middle")
    .attr("font-size", 10)
    .text(vmax.toFixed(0));
  svg
    .append("text")
    .attr("x", cbX + cbW / 2)
    .attr("y", cbY + cbH + 12)
    .attr("text-anchor", "middle")
    .attr("font-size", 10)
    .text("0");

  svg
    .append("text")
    .attr("transform", `translate(${cbX + cbW + 12},${cbY + cbH / 2}) rotate(90)`)
    .attr("text-anchor", "middle")
    .attr("font-size", 9)
    .attr("fill", "#666")
    .text(`Scale: ${scaleMode}`);

  // Rescued-strain legend — rendered as HTML below the SVG so it can
  // never overlap the chart axes or cells regardless of chart dimensions.
  if (_anyRescued) {
    const lgDiv = document.createElement("div");
    lgDiv.style.cssText =
      "display:flex;align-items:center;gap:12px;padding:5px 8px 2px;" + "font-size:0.76em;color:#555;flex-wrap:wrap;";
    lgDiv.innerHTML =
      "<span style='font-weight:600;color:#444;'>Rescued strains:</span>" +
      "<span style='display:flex;align-items:center;gap:4px;'>" +
      "<svg width='10' height='10' style='flex-shrink:0'>" +
      "<path d='M0,0 L10,0 L10,10 Z' fill='#fb8c00'/>" +
      "</svg>" +
      "<span style='color:#666;'>species roll-up</span>" +
      "</span>" +
      "<span style='display:flex;align-items:center;gap:4px;'>" +
      "<svg width='10' height='10' style='flex-shrink:0'>" +
      "<path d='M0,0 L10,0 L10,10 Z' fill='#8e24aa'/>" +
      "</svg>" +
      "<span style='color:#666;'>genus roll-up</span>" +
      "</span>";
    wrap.appendChild(lgDiv);
  }
}
