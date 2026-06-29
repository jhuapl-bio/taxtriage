      /* ═══════════════════════════════════════════════════════════════════════
       -  §  COMPARE COVERAGES — overlapping depth histograms + profiles
       -     Opened from (a) the detections "Compare" modal (pinned rows) and
       -     (b) clicking an organism in Cross-Sample Feature Compare. Overlays
       -     per-sample depth-bin distributions and per-position coverage
       -     profiles; legend toggles hide/show each sample; when 2+ organisms
       -     are present a second "differences" row is added.
       ═══════════════════════════════════════════════════════════════════════ */
      let _COV = null;
      function _covFetch(sample, organism, taxid) {
        const map = _histMap();
        let cd = map.get(`${sample}||${organism}||${taxid}`);
        if (cd) return cd;
        for (const v of map.values()) {
          if (v.sample !== sample) continue;
          if ((taxid && String(v.taxon_id) === String(taxid)) || (organism && v.organism === organism)) return v;
        }
        return null;
      }
      function _covDepth(cd) {
        const hist = (cd && cd.depth_histogram) || {};
        const counts = _DEPTH_BINS.map((b) => num(hist[b]));
        const total = counts.reduce((s, v) => s + v, 0);
        const frac = counts.map((v) => (total ? (v / total) * 100 : 0));
        return { counts, frac, total };
      }
      function _covProfile(cd, N) {
        const pos = _posBins(cd);
        if (!pos || !pos.length) return null;
        N = N || 140;
        if (pos.length <= N) return pos.slice();
        const step = pos.length / N,
          out = [];
        for (let i = 0; i < N; i++) {
          const a = Math.floor(i * step),
            b = Math.max(a + 1, Math.floor((i + 1) * step));
          let s = 0,
            c = 0;
          for (let j = a; j < b && j < pos.length; j++) {
            s += pos[j];
            c++;
          }
          out.push(c ? s / c : 0);
        }
        return out;
      }
      // Overlapping bars: one semi-transparent bar per series at each category,
      // taller bars drawn first so smaller ones stay visible on top.
      function _covOverlapBars(series, cats, W, H, opts) {
        opts = opts || {};
        const mL = 44,
          mR = 12,
          mT = 12,
          mB = 34;
        const iw = W - mL - mR,
          ih = H - mT - mB;
        let maxV = 0;
        series.forEach((s) =>
          s.vals.forEach((v) => {
            if (v > maxV) maxV = v;
          }),
        );
        maxV = maxV || 1;
        const x0 = (i) => mL + (i + 0.5) * (iw / cats.length);
        const bw = Math.min(70, (iw / cats.length) * 0.7);
        const y = (v) => mT + ih - (v / maxV) * ih;
        let svg = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" style="max-width:100%;height:auto">`;
        for (let g = 0; g <= 4; g++) {
          const gv = (maxV * g) / 4,
            gy = y(gv);
          svg += `<line x1="${mL}" x2="${W - mR}" y1="${gy}" y2="${gy}" stroke="#eceff1"/>`;
          svg += `<text x="${mL - 6}" y="${gy + 3}" text-anchor="end" font-size="9" fill="#90a4ae">${gv.toFixed(0)}${
            opts.yunit || ""
          }</text>`;
        }
        cats.forEach((c, i) => {
          const cx = x0(i);
          const order = series.map((s, si) => ({ si, v: s.vals[i] || 0 })).sort((a, b) => b.v - a.v);
          order.forEach(({ si }) => {
            const s = series[si];
            const v = s.vals[i] || 0;
            const h = (v / maxV) * ih;
            svg += `<rect x="${(cx - bw / 2).toFixed(1)}" y="${y(v).toFixed(1)}" width="${bw.toFixed(
              1,
            )}" height="${Math.max(0, h).toFixed(1)}" fill="${s.color}" fill-opacity="0.4" stroke="${
              s.color
            }" stroke-width="1.2"/>`;
          });
          svg += `<text x="${cx}" y="${mT + ih + 14}" text-anchor="middle" font-size="9" fill="#546e7a">${c}</text>`;
        });
        svg += `</svg>`;
        return svg;
      }
      function _covDeltaBars(vals, cats, posColor, negColor, W, H) {
        const mL = 44,
          mR = 12,
          mT = 10,
          mB = 30;
        const iw = W - mL - mR,
          ih = H - mT - mB;
        let m = 0;
        vals.forEach((v) => {
          if (Math.abs(v) > m) m = Math.abs(v);
        });
        m = m || 1;
        const yc = mT + ih / 2;
        const y = (v) => yc - (v / m) * (ih / 2);
        const x0 = (i) => mL + (i + 0.5) * (iw / cats.length);
        const bw = Math.min(60, (iw / cats.length) * 0.6);
        let svg = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" style="max-width:100%;height:auto">`;
        svg += `<line x1="${mL}" x2="${W - mR}" y1="${yc}" y2="${yc}" stroke="#b0bec5"/>`;
        svg += `<text x="${mL - 6}" y="${mT + 8}" text-anchor="end" font-size="9" fill="#90a4ae">+${m.toFixed(
          0,
        )}</text>`;
        svg += `<text x="${mL - 6}" y="${mT + ih}" text-anchor="end" font-size="9" fill="#90a4ae">−${m.toFixed(
          0,
        )}</text>`;
        cats.forEach((c, i) => {
          const v = vals[i] || 0,
            cx = x0(i),
            yy = y(v),
            h = Math.abs(yy - yc);
          const col = v >= 0 ? posColor : negColor;
          svg += `<rect x="${(cx - bw / 2).toFixed(1)}" y="${Math.min(yy, yc).toFixed(1)}" width="${bw.toFixed(
            1,
          )}" height="${Math.max(1, h).toFixed(1)}" fill="${col}" fill-opacity="0.55" stroke="${col}"/>`;
          svg += `<text x="${cx}" y="${mT + ih + 14}" text-anchor="middle" font-size="9" fill="#546e7a">${c}</text>`;
        });
        svg += `</svg>`;
        return svg;
      }
      function _covProfileChart(series, W, H) {
        const mL = 36,
          mR = 12,
          mT = 10,
          mB = 22;
        const iw = W - mL - mR,
          ih = H - mT - mB;
        let svg = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" style="max-width:100%;height:auto">`;
        for (let g = 0; g <= 4; g++) {
          const gv = (100 * g) / 4,
            gy = mT + ih - (gv / 100) * ih;
          svg += `<line x1="${mL}" x2="${W - mR}" y1="${gy}" y2="${gy}" stroke="#eceff1"/>`;
          svg += `<text x="${mL - 5}" y="${gy + 3}" text-anchor="end" font-size="8" fill="#90a4ae">${gv}</text>`;
        }
        series.forEach((s) => {
          const n = s.vals.length;
          const px = (i) => mL + (n <= 1 ? 0 : (i / (n - 1)) * iw);
          const py = (v) => mT + ih - (Math.min(100, v) / 100) * ih;
          const d = s.vals.map((v, i) => `${i ? "L" : "M"}${px(i).toFixed(1)},${py(v).toFixed(1)}`).join(" ");
          const area = `${d} L${px(n - 1).toFixed(1)},${(mT + ih).toFixed(1)} L${px(0).toFixed(1)},${(mT + ih).toFixed(
            1,
          )} Z`;
          svg += `<path d="${area}" fill="${s.color}" fill-opacity="0.13"/>`;
          svg += `<path d="${d}" fill="none" stroke="${s.color}" stroke-width="1.4" stroke-opacity="0.9"/>`;
        });
        svg += `<text x="${mL}" y="${H - 5}" font-size="8" fill="#90a4ae">genome start</text>`;
        svg += `<text x="${W - mR}" y="${H - 5}" text-anchor="end" font-size="8" fill="#90a4ae">end</text>`;
        svg += `</svg>`;
        return svg;
      }
      // ── Aggregation helpers (mean / median across samples) ───────────────
      function _covMedian(arr) {
        const a = arr
          .filter((v) => v != null && !isNaN(v))
          .slice()
          .sort((x, y) => x - y);
        if (!a.length) return 0;
        const m = Math.floor(a.length / 2);
        return a.length % 2 ? a[m] : (a[m - 1] + a[m]) / 2;
      }
      // Aggregate a list of equal-length (or varied) numeric arrays element-wise.
      function _covAggArr(arrays, stat) {
        if (!arrays.length) return [];
        const L = Math.min(...arrays.map((a) => a.length));
        const out = [];
        for (let i = 0; i < L; i++) {
          const col = arrays.map((a) => a[i] || 0);
          out.push(stat === "median" ? _covMedian(col) : col.reduce((s, v) => s + v, 0) / col.length);
        }
        return out;
      }
      function _openCoverageModal(rawItems, opts) {
        opts = opts || {};
        const ov = document.getElementById("cov-overlay");
        if (!ov) return;
        const seen = new Set(),
          items = [];
        (rawItems || []).forEach((it) => {
          if (!it || !it.sample) return;
          const key = `${it.sample}||${it.organism}||${it.taxid}`;
          if (seen.has(key)) return;
          seen.add(key);
          items.push({ sample: it.sample, organism: it.organism || "", taxid: it.taxid || "" });
        });
        const orgNames = uniq(items.map((it) => it.organism));
        const orgColor = new Map();
        orgNames.forEach((nm, i) => orgColor.set(nm, _XS_SAMPLE_PALETTE[i % _XS_SAMPLE_PALETTE.length]));
        items.forEach((it, i) => {
          it.id = `cov${i}`;
          // Use the document-wide per-sample colour (same swatch shown in the
          // right-hand sidebar) so coverage-profile series line up with every
          // other chart. Fall back to the generic palette only for samples that
          // somehow have no assigned colour yet.
          it.color = sampleColors[it.sample] || _XS_SAMPLE_PALETTE[i % _XS_SAMPLE_PALETTE.length];
          it.cd = _covFetch(it.sample, it.organism, it.taxid);
          it.depth = it.cd ? _covDepth(it.cd) : null;
          it.profile = it.cd ? _covProfile(it.cd) : null;
        });
        _COV = {
          items,
          hidden: new Set(),
          orgNames,
          orgColor,
          multiOrg: orgNames.length > 1,
          mode: "aggregate",
          aggStat: "mean",
          selectedOrg: orgNames[0] || "",
        };
        const titleEl = document.getElementById("cov-title");
        if (titleEl) titleEl.innerHTML = opts.title || "Coverage comparison";
        _covRender();
        ov.style.display = "flex";
      }
      function _openCoverageModalForOrg(o, agg) {
        const samples = [...(o.samples || [])].sort();
        const items = samples.map((s) => ({ sample: s, organism: o.name, taxid: o.taxid }));
        _openCoverageModal(items, { title: `Coverage across samples — <i>${o.name}</i>` });
      }
      function _openCoverageModalForRows(rows) {
        const items = (rows || []).map((r) => ({
          sample: r["Specimen ID"],
          organism: r["Detected Organism"],
          taxid: r["Taxonomic ID #"],
        }));
        _openCoverageModal(items, { title: "Coverage comparison — pinned detections" });
      }
      // When several organisms are present we show ONE at a time (picked via the
      // searchable dropdown); this returns the items for the currently selected
      // organism. With a single organism it returns everything unchanged.
      function _covActiveItems() {
        const C = _COV;
        if (!C) return [];
        if (C.orgNames.length > 1 && C.selectedOrg) {
          return C.items.filter((it) => it.organism === C.selectedOrg);
        }
        return C.items;
      }
      function _covRender() {
        const body = document.getElementById("cov-body");
        if (!body || !_COV) return;
        const C = _COV;
        const act = _covActiveItems();
        const withData = act.filter((it) => it.cd);
        if (!withData.length) {
          body.innerHTML =
            `<div style="padding:1.4em;color:#777;font-size:.9em"><i class="fas fa-circle-info" style="color:#f59f00"></i> ` +
            `No per-contig depth/coverage data is attached for the selected organism(s)/sample(s). ` +
            `Depth histograms are only present when the report is generated by the pipeline (data loaded from XLSX/TSV uploads lacks them).</div>`;
          return;
        }
        // ── Organism picker (only when 2+ organisms) — one at a time ─────────
        const _esc = (s) => String(s).replace(/"/g, "&quot;");
        const orgPicker =
          C.orgNames.length > 1
            ? `<div class="cov-ctrl">` +
              `<span class="cov-ctrl-lbl">Organism:</span>` +
              `<select id="cov-org-input" class="cov-org-input">` +
              C.orgNames
                .map((n) => `<option value="${_esc(n)}"${n === C.selectedOrg ? " selected" : ""}>${_esc(n)}</option>`)
                .join("") +
              `</select>` +
              `<span class="cov-ctrl-lbl">${C.orgNames.length} organisms · showing one at a time</span>` +
              `</div>`
            : "";
        // ── Aggregate / Per-sample mode toggle ──────────────────────────────
        const ctrl =
          `<div class="cov-ctrl">` +
          `<span class="cov-ctrl-lbl">View:</span>` +
          `<span class="cov-seg">` +
          `<button class="cov-segbtn${
            C.mode === "aggregate" ? " on" : ""
          }" data-mode="aggregate" title="Combine all shown samples into one averaged coverage curve">Aggregate</button>` +
          `<button class="cov-segbtn${
            C.mode === "sample" ? " on" : ""
          }" data-mode="sample" title="Overlay every shown sample separately">Per-sample</button>` +
          `</span>` +
          (C.mode === "aggregate"
            ? `<span class="cov-ctrl-lbl">across samples:</span>` +
              `<span class="cov-seg">` +
              `<button class="cov-segbtn sm${C.aggStat === "mean" ? " on" : ""}" data-stat="mean">Mean</button>` +
              `<button class="cov-segbtn sm${C.aggStat === "median" ? " on" : ""}" data-stat="median">Median</button>` +
              `</span>`
            : "") +
          `</div>`;
        const legend =
          `<div class="cov-legend">` +
          act
            .map((it) => {
              const off = C.hidden.has(it.id) || !it.cd;
              const dis = !it.cd;
              return (
                `<span class="cov-leg${off ? " off" : ""}" data-id="${it.id}" title="${
                  dis ? "no coverage data" : "click to hide / show (also includes / excludes it from the aggregate)"
                }" style="${dis ? "opacity:.4;cursor:not-allowed" : "cursor:pointer"}">` +
                `<span class="cov-sw" style="background:${it.color}"></span>${it.sample}` +
                (dis ? ' <span style="color:#c0392b">(no data)</span>' : "") +
                `</span>`
              );
            })
            .join("") +
          `</div>`;
        body.innerHTML = orgPicker + ctrl + legend + `<div id="cov-charts"></div><div id="cov-stats"></div>`;
        // Organism dropdown wiring (plain select)
        const orgInput = document.getElementById("cov-org-input");
        if (orgInput) {
          orgInput.addEventListener("change", () => {
            const v = orgInput.value;
            if (v && v !== C.selectedOrg) {
              C.selectedOrg = v;
              _covRender();
            }
          });
        }
        // Mode + stat toggle wiring
        body.querySelectorAll(".cov-segbtn[data-mode]").forEach((el) => {
          el.addEventListener("click", () => {
            const m = el.getAttribute("data-mode");
            if (m && m !== C.mode) {
              C.mode = m;
              _covRender();
            }
          });
        });
        body.querySelectorAll(".cov-segbtn[data-stat]").forEach((el) => {
          el.addEventListener("click", () => {
            const s = el.getAttribute("data-stat");
            if (s && s !== C.aggStat) {
              C.aggStat = s;
              _covRender();
            }
          });
        });
        body.querySelectorAll(".cov-leg").forEach((el) => {
          const it = C.items.find((x) => x.id === el.getAttribute("data-id"));
          if (!it || !it.cd) return;
          el.addEventListener("click", () => {
            if (C.hidden.has(it.id)) C.hidden.delete(it.id);
            else C.hidden.add(it.id);
            _covRender();
          });
        });
        _covDraw();
        _covStats();
      }
      function _covDraw() {
        const host = document.getElementById("cov-charts");
        if (!host || !_COV) return;
        const C = _COV;
        const act = _covActiveItems();
        const vis = act.filter((it) => it.cd && !C.hidden.has(it.id));
        let html = "";

        // ── Aggregate mode: one averaged curve per organism across samples ──
        if (C.mode === "aggregate") {
          const byOrg = new Map();
          vis.forEach((it) => {
            if (!byOrg.has(it.organism)) byOrg.set(it.organism, []);
            byOrg.get(it.organism).push(it);
          });
          const statLbl = C.aggStat === "median" ? "Median" : "Mean";
          const depthSeries = [];
          const profSeries = [];
          byOrg.forEach((its, org) => {
            const color = byOrg.size > 1 ? C.orgColor.get(org) || its[0].color : "#4527a0";
            const tag = (byOrg.size > 1 ? org + " " : "") + `(${statLbl.toLowerCase()} of ${its.length})`;
            const dFrac = _covAggArr(
              its.map((it) => it.depth.frac),
              C.aggStat,
            );
            if (dFrac.length) depthSeries.push({ label: tag, color, vals: dFrac });
            const profs = its.filter((it) => it.profile && it.profile.length).map((it) => it.profile);
            if (profs.length) {
              const p = _covAggArr(profs, C.aggStat);
              if (p.length) profSeries.push({ label: tag, color, vals: p });
            }
          });
          const head = `${statLbl} across ${vis.length} sample(s)` + (byOrg.size > 1 ? " per organism" : "");
          html +=
            `<div class="cov-card"><div class="cov-card-h">Depth distribution — % of genome positions per depth bin · ${head}</div>` +
            (depthSeries.length
              ? _covOverlapBars(depthSeries, _DEPTH_BINS, 560, 230, { yunit: "%" })
              : `<div class="cov-empty">No visible samples — toggle one back on below.</div>`) +
            `</div>`;
          if (profSeries.length) {
            html +=
              `<div class="cov-card"><div class="cov-card-h">Coverage profile — % covered across genome positions · ${head}</div>` +
              _covProfileChart(profSeries, 560, 180) +
              `</div>`;
          }
          host.innerHTML = html;
          return;
        }

        // ── Per-sample mode: original overlaid-per-sample behaviour ──
        const depthSeries = vis.map((it) => ({ label: it.sample, color: it.color, vals: it.depth.frac }));
        html +=
          `<div class="cov-card"><div class="cov-card-h">Depth distribution — % of genome positions per depth bin (overlaid)</div>` +
          (depthSeries.length
            ? _covOverlapBars(depthSeries, _DEPTH_BINS, 560, 230, { yunit: "%" })
            : `<div class="cov-empty">No visible samples — toggle one back on below.</div>`) +
          `</div>`;
        const profSeries = vis
          .filter((it) => it.profile && it.profile.length)
          .map((it) => ({ label: it.sample, color: it.color, vals: it.profile }));
        if (profSeries.length) {
          html +=
            `<div class="cov-card"><div class="cov-card-h">Coverage profile — % covered across genome positions (overlaid)</div>` +
            _covProfileChart(profSeries, 560, 180) +
            `</div>`;
        }
        host.innerHTML = html;
        // With the one-organism-at-a-time picker, vis carries a single organism,
        // so the cross-organism difference card no longer applies. Guard on the
        // distinct organism count so it only ever shows if that changes.
        const _distinctOrgs = new Set(vis.map((it) => it.organism));
        if (_distinctOrgs.size > 1) _covDrawDiff(host, vis);
      }
      function _covDrawDiff(host, vis) {
        const C = _COV;
        const byOrg = new Map();
        vis.forEach((it) => {
          if (!byOrg.has(it.organism)) byOrg.set(it.organism, []);
          byOrg.get(it.organism).push(it);
        });
        const orgs = [...byOrg.keys()];
        if (orgs.length < 2) return;
        const meanFrac = (its) =>
          _DEPTH_BINS.map((_, bi) => {
            const a = its.map((it) => it.depth.frac[bi]);
            return a.reduce((s, v) => s + v, 0) / (a.length || 1);
          });
        const series = orgs.map((nm) => ({ label: nm, color: C.orgColor.get(nm), vals: meanFrac(byOrg.get(nm)) }));
        const a = series[0].vals,
          b = series[1].vals;
        const delta = a.map((v, i) => v - b[i]);
        let html = `<div class="cov-card cov-diff"><div class="cov-card-h"><i class="fas fa-code-compare"></i> Organism differences — mean depth distribution per organism (overlaid)</div>`;
        html += _covOverlapBars(series, _DEPTH_BINS, 560, 220, { yunit: "%" });
        html += `<div class="cov-card-h" style="margin-top:.5em">Δ <span style="color:${series[0].color}">${orgs[0]}</span> − <span style="color:${series[1].color}">${orgs[1]}</span> (percentage points per bin)</div>`;
        html += _covDeltaBars(delta, _DEPTH_BINS, series[0].color, series[1].color, 560, 150);
        const orgProf = orgs
          .map((nm) => {
            const its = byOrg.get(nm).filter((it) => it.profile && it.profile.length);
            if (!its.length) return null;
            const L = Math.min(...its.map((it) => it.profile.length));
            const vals = [];
            for (let i = 0; i < L; i++) vals.push(its.reduce((s, it) => s + it.profile[i], 0) / its.length);
            return { label: nm, color: C.orgColor.get(nm), vals };
          })
          .filter(Boolean);
        if (orgProf.length >= 2) {
          html +=
            `<div class="cov-card-h" style="margin-top:.5em">Coverage profile by organism (overlaid)</div>` +
            _covProfileChart(orgProf, 560, 160);
        }
        html += `</div>`;
        host.insertAdjacentHTML("beforeend", html);
      }
      function _covStats() {
        const host = document.getElementById("cov-stats");
        if (!host || !_COV) return;
        const C = _COV;
        const act = _covActiveItems();
        const mids = [0, 3, 7.5, 30, 75];
        const rows = act
          .filter((it) => it.cd)
          .map((it) => {
            const d = it.depth,
              tot = d.total;
            const breadth = 100 - d.frac[0];
            const deep = d.frac[4];
            const mean = tot ? d.counts.reduce((s, v, i) => s + v * mids[i], 0) / tot : 0;
            const dim = C.hidden.has(it.id) ? "opacity:.45" : "";
            return (
              `<tr style="${dim}">` +
              `<td style="padding:4px 8px;border-bottom:1px solid #eee"><span class="cov-sw" style="background:${it.color}"></span>${it.sample}</td>` +
              `<td style="padding:4px 8px;border-bottom:1px solid #eee"><i>${it.organism}</i></td>` +
              `<td style="padding:4px 8px;border-bottom:1px solid #eee;text-align:right">${_fmtInt(tot)}</td>` +
              `<td style="padding:4px 8px;border-bottom:1px solid #eee;text-align:right">${breadth.toFixed(1)}%</td>` +
              `<td style="padding:4px 8px;border-bottom:1px solid #eee;text-align:right">${deep.toFixed(1)}%</td>` +
              `<td style="padding:4px 8px;border-bottom:1px solid #eee;text-align:right">${mean.toFixed(1)}×</td>` +
              `</tr>`
            );
          })
          .join("");
        const vis = act.filter((it) => it.cd && !C.hidden.has(it.id));
        let summary = "";
        if (vis.length) {
          const br = vis.map((it) => 100 - it.depth.frac[0]);
          const avgBr = br.reduce((s, v) => s + v, 0) / br.length;
          summary = `<div style="font-size:.8em;color:#555;margin:.3em 0 .5em">${
            vis.length
          } sample(s) shown · mean breadth ${avgBr.toFixed(1)}% · spread ${Math.min(...br).toFixed(0)}–${Math.max(
            ...br,
          ).toFixed(0)}%</div>`;
        }

        // ── Aggregate summary (mean / median across the shown samples) ──
        let aggCard = "";
        if (C.mode === "aggregate" && vis.length) {
          const stat = C.aggStat;
          const statLbl = stat === "median" ? "Median" : "Mean";
          const aggOf = (arr) =>
            stat === "median" ? _covMedian(arr) : arr.reduce((s, v) => s + v, 0) / (arr.length || 1);
          const metric = (it) => {
            const d = it.depth,
              tot = d.total;
            return {
              breadth: 100 - d.frac[0],
              deep: d.frac[4],
              mean: tot ? d.counts.reduce((s, v, i) => s + v * mids[i], 0) / tot : 0,
            };
          };
          const byOrg = new Map();
          vis.forEach((it) => {
            if (!byOrg.has(it.organism)) byOrg.set(it.organism, []);
            byOrg.get(it.organism).push(it);
          });
          const aggRows = [...byOrg.entries()]
            .map(([org, its]) => {
              const ms = its.map(metric);
              const brs = ms.map((x) => x.breadth);
              const color = C.multiOrg ? C.orgColor.get(org) || its[0].color : "#4527a0";
              const label = C.multiOrg
                ? `<span class="cov-sw" style="background:${color}"></span><i>${org}</i>`
                : `<span class="cov-sw" style="background:${color}"></span>All samples`;
              return (
                `<tr>` +
                `<td style="padding:4px 8px;border-bottom:1px solid #eee">${label}</td>` +
                `<td style="padding:4px 8px;border-bottom:1px solid #eee;text-align:right">${its.length}</td>` +
                `<td style="padding:4px 8px;border-bottom:1px solid #eee;text-align:right">${aggOf(brs).toFixed(
                  1,
                )}%</td>` +
                `<td style="padding:4px 8px;border-bottom:1px solid #eee;text-align:right">${aggOf(
                  ms.map((x) => x.deep),
                ).toFixed(1)}%</td>` +
                `<td style="padding:4px 8px;border-bottom:1px solid #eee;text-align:right">${aggOf(
                  ms.map((x) => x.mean),
                ).toFixed(1)}×</td>` +
                `<td style="padding:4px 8px;border-bottom:1px solid #eee;text-align:right">${Math.min(...brs).toFixed(
                  0,
                )}–${Math.max(...brs).toFixed(0)}%</td>` +
                `</tr>`
              );
            })
            .join("");
          aggCard =
            `<div class="cov-card-h" style="margin-top:.6em">${statLbl} coverage across samples</div>` +
            `<div style="overflow:auto"><table style="border-collapse:collapse;width:100%;font-size:.82em">` +
            `<thead><tr style="background:#efeaf8">` +
            `<th style="padding:5px 8px;text-align:left;border-bottom:2px solid #4527a0">${
              C.multiOrg ? "Organism" : "Group"
            }</th>` +
            `<th style="padding:5px 8px;text-align:right;border-bottom:2px solid #4527a0"># Samples</th>` +
            `<th style="padding:5px 8px;text-align:right;border-bottom:2px solid #4527a0">${statLbl} Breadth</th>` +
            `<th style="padding:5px 8px;text-align:right;border-bottom:2px solid #4527a0">${statLbl} &gt;50×</th>` +
            `<th style="padding:5px 8px;text-align:right;border-bottom:2px solid #4527a0">${statLbl} depth≈</th>` +
            `<th style="padding:5px 8px;text-align:right;border-bottom:2px solid #4527a0">Breadth range</th>` +
            `</tr></thead><tbody>${aggRows}</tbody></table></div>`;
        }

        host.innerHTML =
          aggCard +
          `<div class="cov-card-h" style="margin-top:.6em">Per-sample coverage stats</div>` +
          summary +
          `<div style="overflow:auto"><table style="border-collapse:collapse;width:100%;font-size:.82em">` +
          `<thead><tr style="background:#f0f4fa">` +
          `<th style="padding:5px 8px;text-align:left;border-bottom:2px solid #1565c0">Sample</th>` +
          `<th style="padding:5px 8px;text-align:left;border-bottom:2px solid #1565c0">Organism</th>` +
          `<th style="padding:5px 8px;text-align:right;border-bottom:2px solid #1565c0">Genome positions</th>` +
          `<th style="padding:5px 8px;text-align:right;border-bottom:2px solid #1565c0">Breadth (&gt;0×)</th>` +
          `<th style="padding:5px 8px;text-align:right;border-bottom:2px solid #1565c0">&gt;50× positions</th>` +
          `<th style="padding:5px 8px;text-align:right;border-bottom:2px solid #1565c0">Mean depth≈</th>` +
          `</tr></thead><tbody>${rows}</tbody></table></div>`;
      }

      // ── View 1: frequency table (xlsx-style organism roll-up) ────────────
      const _XS_COLS = [
        { key: "name", label: "Detected Organism", left: true },
        { key: "cat", label: "Pathogen Type", left: true },
        { key: "taxid", label: "TaxID", left: true },
        { key: "sampleCount", label: "# Samples" },
        { key: "samplePct", label: "% Samples" },
        { key: "meanTass", label: "Mean TASS" },
        { key: "medianTass", label: "Median TASS" },
        { key: "maxTass", label: "Max TASS" },
        { key: "meanCov", label: "Mean Cov" },
        { key: "medianCov", label: "Median Cov" },
        { key: "maxCov", label: "Max Cov" },
        { key: "reads", label: "Total Reads" },
      ];
      function _xsRenderTable(rows, agg) {
        const k = _XS.sortKey,
          asc = _XS.sortAsc;
        const sorted = [...rows].sort((a, b) => {
          let av = a[k],
            bv = b[k];
          if (typeof av === "string") {
            av = av.toLowerCase();
            bv = (bv || "").toLowerCase();
            return asc ? av.localeCompare(bv) : bv.localeCompare(av);
          }
          return asc ? av - bv : bv - av;
        });
        const ind = (key) => (_XS.sortKey === key ? (_XS.sortAsc ? " ▲" : " ▼") : "");
        const fmt1 = (v) => (v == null ? "—" : v.toFixed(1));
        const head =
          "<tr>" +
          _XS_COLS
            .map(
              (c) =>
                `<th data-key="${c.key}" class="${c.left ? "xs-left" : ""}" title="click to sort">${c.label}${ind(
                  c.key,
                )}</th>`,
            )
            .join("") +
          "</tr>";
        const bodyRows = sorted
          .map((r, i) => {
            const pill = `<span class="xs-cat-pill" style="background:${_xsCatColor(r.cat)}">${r.cat}</span>`;
            const star = r.hc ? '<span title="high-consequence" style="color:#cc0000">● </span>' : "";
            const prev =
              `<span class="xs-prevbar"><span style="width:${Math.min(100, r.samplePct).toFixed(0)}%"></span></span>` +
              `${r.samplePct.toFixed(1)}%`;
            const samplesAttr = [...r.samples].sort().join("||");
            return (
              `<tr data-samples="${samplesAttr.replace(
                /"/g,
                "&quot;",
              )}" data-oi="${i}" title="Click to compare coverage across this organism's samples" style="cursor:pointer">` +
              `<td class="xs-left"><i>${star}${r.name}</i> <i class="fas fa-chart-column" style="color:#bcc6d0;font-size:.82em"></i></td>` +
              `<td class="xs-left">${pill}</td>` +
              `<td class="xs-left" style="color:#789">${r.taxid}</td>` +
              `<td class="xs-num"><b>${r.sampleCount}</b> / ${agg.totalSamples}</td>` +
              `<td class="xs-num">${prev}</td>` +
              `<td class="xs-num">${fmt1(r.meanTass)}</td>` +
              `<td class="xs-num">${fmt1(r.medianTass)}</td>` +
              `<td class="xs-num">${fmt1(r.maxTass)}</td>` +
              `<td class="xs-num">${fmt1(r.meanCov)}</td>` +
              `<td class="xs-num">${fmt1(r.medianCov)}</td>` +
              `<td class="xs-num">${fmt1(r.maxCov)}</td>` +
              `<td class="xs-num">${_fmtInt(r.reads)}</td>` +
              "</tr>"
            );
          })
          .join("");
        body_set(
          `<div class="xs-note">${rows.length} organism(s) across ${agg.totalSamples} sample(s). Click a column header to sort · <span style="color:#cc0000">●</span> = high-consequence · prevalence bar = % of samples detected · hover a row for sample list · click a row to compare coverage across samples.</div>` +
            `<div id="xs-table-wrap"><table><thead>${head}</thead><tbody>${bodyRows}</tbody></table></div>`,
        );
        // Delegated tooltip: show sample list on row hover
        const tbody = document.querySelector("#xs-table-wrap tbody");
        if (tbody) {
          tbody.addEventListener("mouseover", (ev) => {
            const tr = ev.target.closest("tr[data-samples]");
            if (!tr) return;
            const names = tr.dataset.samples.split("||").filter(Boolean);
            showTip(`<b>Detected in ${names.length} sample(s):</b><br>` + names.map((s) => `• ${s}`).join("<br>"), ev);
          });
          tbody.addEventListener("mousemove", moveTip);
          tbody.addEventListener("mouseleave", hideTip);
          // Click an organism row → overlay its depth coverage across samples.
          tbody.addEventListener("click", (ev) => {
            const tr = ev.target.closest("tr[data-oi]");
            if (!tr) return;
            const o = sorted[+tr.getAttribute("data-oi")];
            if (o) {
              hideTip();
              _openCoverageModalForOrg(o, agg);
            }
          });
        }
      }
      function body_set(html) {
        document.getElementById("xs-body").innerHTML = html;
      }

      // ── View 2: prevalence + TASS distribution (box & whisker) ───────────
      function _xsRenderDist(rows, agg) {
        const topN = parseInt((document.getElementById("xs-topn") || {}).value || "25", 10);
        const ordered = [...rows]
          .sort((a, b) => b.sampleCount - a.sampleCount || b.medianTass - a.medianTass)
          .slice(0, topN);
        body_set(
          '<div class="xs-note">Top organisms by prevalence. Bar = % of samples detected · box = TASS spread across samples (Q1–Q3, line = median, whiskers = min/max). Tighter boxes = more consistent scoring.</div><div id="xs-dist-svg"></div>',
        );
        const wrap = document.getElementById("xs-dist-svg");
        const W = wrap.clientWidth || 760;
        const rowH = 26,
          mT = 24,
          mB = 30,
          mL = 200,
          mGap = 24;
        const H = ordered.length * rowH + mT + mB;
        const prevW = 120; // width reserved for prevalence bar
        const boxX0 = mL + prevW + mGap;
        const boxW = Math.max(160, W - boxX0 - 60);
        const x = d3.scaleLinear().domain([0, 100]).range([0, boxW]);
        const xp = d3.scaleLinear().domain([0, 100]).range([0, prevW]);
        const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);

        // TASS axis ticks
        [0, 25, 50, 75, 100].forEach((t) => {
          svg
            .append("line")
            .attr("x1", boxX0 + x(t))
            .attr("x2", boxX0 + x(t))
            .attr("y1", mT - 6)
            .attr("y2", H - mB + 4)
            .attr("stroke", "#eceff1");
          svg
            .append("text")
            .attr("x", boxX0 + x(t))
            .attr("y", mT - 10)
            .attr("text-anchor", "middle")
            .style("font-size", "9px")
            .style("fill", "#90a4ae")
            .text(t);
        });
        svg
          .append("text")
          .attr("x", boxX0 + boxW / 2)
          .attr("y", H - 6)
          .attr("text-anchor", "middle")
          .style("font-size", "10px")
          .style("fill", "#607d8b")
          .text("TASS score across samples");
        svg
          .append("text")
          .attr("x", mL + prevW / 2)
          .attr("y", H - 6)
          .attr("text-anchor", "middle")
          .style("font-size", "10px")
          .style("fill", "#607d8b")
          .text("% of samples");

        ordered.forEach((r, i) => {
          const y = mT + i * rowH + rowH / 2;
          const col = _xsCatColor(r.cat);
          // organism label
          svg
            .append("text")
            .attr("x", mL - 8)
            .attr("y", y)
            .attr("text-anchor", "end")
            .attr("dominant-baseline", "middle")
            .style("font-size", "10px")
            .style("font-style", "italic")
            .style("fill", "#37474f")
            .text(r.name.length > 30 ? r.name.slice(0, 29) + "…" : r.name)
            .append("title")
            .text(`${r.name} (${r.cat})`);
          // prevalence bar
          svg
            .append("rect")
            .attr("x", mL)
            .attr("y", y - 7)
            .attr("width", prevW)
            .attr("height", 14)
            .attr("fill", "#eceff1")
            .attr("rx", 2);
          svg
            .append("rect")
            .attr("x", mL)
            .attr("y", y - 7)
            .attr("width", xp(r.samplePct))
            .attr("height", 14)
            .attr("fill", col)
            .attr("opacity", 0.85)
            .attr("rx", 2);
          // On-chart count label removed — the prevalence (count / total / %) and
          // full TASS spread are shown on hover over the row instead.
          // box & whisker for TASS
          const s = [...r.tassVals].sort((a, b) => a - b);
          if (!s.length) return;
          const q1 = _xsQuart(s, 0.25),
            med = _xsQuart(s, 0.5),
            q3 = _xsQuart(s, 0.75),
            lo = s[0],
            hi = s[s.length - 1];
          const yb = y;
          // whisker line
          svg
            .append("line")
            .attr("x1", boxX0 + x(lo))
            .attr("x2", boxX0 + x(hi))
            .attr("y1", yb)
            .attr("y2", yb)
            .attr("stroke", col)
            .attr("stroke-width", 1);
          [lo, hi].forEach((v) =>
            svg
              .append("line")
              .attr("x1", boxX0 + x(v))
              .attr("x2", boxX0 + x(v))
              .attr("y1", yb - 4)
              .attr("y2", yb + 4)
              .attr("stroke", col),
          );
          // box
          svg
            .append("rect")
            .attr("x", boxX0 + x(q1))
            .attr("y", yb - 7)
            .attr("width", Math.max(1, x(q3) - x(q1)))
            .attr("height", 14)
            .attr("fill", col)
            .attr("opacity", 0.28)
            .attr("stroke", col);
          // median
          svg
            .append("line")
            .attr("x1", boxX0 + x(med))
            .attr("x2", boxX0 + x(med))
            .attr("y1", yb - 7)
            .attr("y2", yb + 7)
            .attr("stroke", col)
            .attr("stroke-width", 2);
          // hover target — full row (label + bars) for easier hover
          const sampleListTip = [...r.samples]
            .sort()
            .map((s) => `• ${s}`)
            .join("<br>");
          const tipHtml =
            `<b>${r.name}</b><br>${r.cat} · ${r.sampleCount}/${agg.totalSamples} samples (${r.samplePct.toFixed(1)}%)` +
            `<br>TASS — min ${lo.toFixed(1)} · Q1 ${q1.toFixed(1)} · med ${med.toFixed(1)} · Q3 ${q3.toFixed(
              1,
            )} · max ${hi.toFixed(1)}` +
            `<br>mean cov ${r.meanCov.toFixed(1)}` +
            `<br><br><b>Samples:</b><br>${sampleListTip}`;
          svg
            .append("rect")
            .attr("x", 0)
            .attr("y", yb - rowH / 2)
            .attr("width", W)
            .attr("height", rowH)
            .attr("fill", "transparent")
            .style("cursor", "pointer")
            .on("mouseover", (ev) => showTip(tipHtml, ev))
            .on("mousemove", moveTip)
            .on("mouseout", hideTip);
        });
        _scheduleExportEnhance();
      }

      // ── View 3: sample clustering via in-browser PCA ─────────────────────
      // Build sample × organism matrix, mean-center columns, then power-iterate
      // the covariance to get the top-2 principal components. No deps.
      function _xsPCA(matrix) {
        const nS = matrix.length;
        if (!nS) return null;
        const nF = matrix[0].length;
        if (!nF) return null;
        // center columns
        const means = new Array(nF).fill(0);
        for (let j = 0; j < nF; j++) {
          let s = 0;
          for (let i = 0; i < nS; i++) s += matrix[i][j];
          means[j] = s / nS;
        }
        const X = matrix.map((row) => row.map((v, j) => v - means[j]));
        // covariance (nF×nF) — fine for the capped feature count we pass in
        const cov = Array.from({ length: nF }, () => new Array(nF).fill(0));
        for (let a = 0; a < nF; a++) {
          for (let b = a; b < nF; b++) {
            let s = 0;
            for (let i = 0; i < nS; i++) s += X[i][a] * X[i][b];
            s /= Math.max(1, nS - 1);
            cov[a][b] = s;
            cov[b][a] = s;
          }
        }
        function powerIter(C) {
          let v = new Array(C.length).fill(0).map(() => Math.random());
          let norm = Math.hypot(...v) || 1;
          v = v.map((x) => x / norm);
          for (let it = 0; it < 120; it++) {
            const w = new Array(C.length).fill(0);
            for (let a = 0; a < C.length; a++) {
              let s = 0;
              for (let b = 0; b < C.length; b++) s += C[a][b] * v[b];
              w[a] = s;
            }
            norm = Math.hypot(...w) || 1;
            v = w.map((x) => x / norm);
          }
          // eigenvalue (Rayleigh quotient)
          const Cv = new Array(C.length).fill(0);
          for (let a = 0; a < C.length; a++) {
            let s = 0;
            for (let b = 0; b < C.length; b++) s += C[a][b] * v[b];
            Cv[a] = s;
          }
          const lambda = v.reduce((s, x, i) => s + x * Cv[i], 0);
          return { vec: v, val: lambda };
        }
        const pc1 = powerIter(cov);
        // deflate then second component
        const cov2 = cov.map((row, a) => row.map((val, b) => val - pc1.val * pc1.vec[a] * pc1.vec[b]));
        const pc2 = powerIter(cov2);
        const totVar = cov.reduce((s, row, a) => s + row[a], 0) || 1;
        const proj = X.map((row) => [
          row.reduce((s, v, j) => s + v * pc1.vec[j], 0),
          row.reduce((s, v, j) => s + v * pc2.vec[j], 0),
        ]);
        return { proj, ev1: (pc1.val / totVar) * 100, ev2: (pc2.val / totVar) * 100 };
      }
      function _xsRenderPCA(rows, agg) {
        const metric = (document.getElementById("xs-metric") || {}).value || "TASS Score";
        const topN = parseInt((document.getElementById("xs-topn") || {}).value || "25", 10);
        const samples = agg.sampleList;
        if (samples.length < 3) {
          body_set(
            '<p style="color:#888;padding:1em">Need at least 3 samples in view to compute a PCA projection.</p>',
          );
          return;
        }
        // pick most informative organisms (by how many samples they vary across)
        const feats = [...rows]
          .filter((r) => r.sampleCount >= 1)
          .sort((a, b) => b.sampleCount - a.sampleCount)
          .slice(0, Math.max(2, Math.min(80, topN * 2)));
        if (feats.length < 2) {
          body_set('<p style="color:#888;padding:1em">Not enough organisms to build a clustering matrix.</p>');
          return;
        }
        const sIdx = Object.fromEntries(samples.map((s, i) => [s, i]));
        const matrix = samples.map(() => new Array(feats.length).fill(0));
        feats.forEach((f, j) => {
          f.tassMap.forEach((tval, sn) => {
            const i = sIdx[sn];
            if (i == null) return;
            if (metric === "presence") matrix[i][j] = 1;
            else if (metric === "Coverage") matrix[i][j] = f.covMap.get(sn) || 0;
            else matrix[i][j] = tval;
          });
        });
        const pca = _xsPCA(matrix);
        if (!pca) {
          body_set('<p style="color:#888;padding:1em">PCA could not be computed for this view.</p>');
          return;
        }
        body_set(
          `<div class="xs-note">Each point is a sample, positioned by a PCA of its organism ${metric} profile (top ${
            feats.length
          } organisms). Samples that sit close together have similar microbial composition; isolated points are outliers. Colour = sample type. PC1 explains ${pca.ev1.toFixed(
            1,
          )}% · PC2 ${pca.ev2.toFixed(1)}% of variance.</div><div id="xs-pca-svg"></div>`,
        );
        const wrap = document.getElementById("xs-pca-svg");
        const W = wrap.clientWidth || 760,
          H = 460,
          pad = 46;
        const xs = pca.proj.map((p) => p[0]),
          ys = pca.proj.map((p) => p[1]);
        const x = d3
          .scaleLinear()
          .domain([Math.min(...xs), Math.max(...xs)])
          .nice()
          .range([pad, W - pad]);
        const y = d3
          .scaleLinear()
          .domain([Math.min(...ys), Math.max(...ys)])
          .nice()
          .range([H - pad, pad]);
        const types = uniq(samples.map((s) => (SAMPLE_META[s] || {}).sample_type || "unknown"));
        const color = d3.scaleOrdinal().domain(types).range(PALETTE);
        const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
        // axes
        svg
          .append("line")
          .attr("x1", pad)
          .attr("x2", W - pad)
          .attr("y1", y(0))
          .attr("y2", y(0))
          .attr("stroke", "#e0e0e0");
        svg
          .append("line")
          .attr("x1", x(0))
          .attr("x2", x(0))
          .attr("y1", pad)
          .attr("y2", H - pad)
          .attr("stroke", "#e0e0e0");
        svg
          .append("text")
          .attr("x", W / 2)
          .attr("y", H - 12)
          .attr("text-anchor", "middle")
          .style("font-size", "11px")
          .style("fill", "#607d8b")
          .text(`PC1 (${pca.ev1.toFixed(1)}%)`);
        svg
          .append("text")
          .attr("transform", `translate(14,${H / 2}) rotate(-90)`)
          .attr("text-anchor", "middle")
          .style("font-size", "11px")
          .style("fill", "#607d8b")
          .text(`PC2 (${pca.ev2.toFixed(1)}%)`);
        samples.forEach((sn, i) => {
          const st = (SAMPLE_META[sn] || {}).sample_type || "unknown";
          svg
            .append("circle")
            .attr("cx", x(pca.proj[i][0]))
            .attr("cy", y(pca.proj[i][1]))
            .attr("r", 6)
            .attr("fill", color(st))
            .attr("opacity", 0.8)
            .attr("stroke", "#fff")
            .attr("stroke-width", 1)
            .style("cursor", "pointer")
            .on("mouseover", (ev) => showTip(`<b>${sn}</b><br>type: ${st}`, ev))
            .on("mousemove", moveTip)
            .on("mouseout", hideTip);
        });
        // legend
        const lg = svg.append("g").attr("transform", `translate(${pad + 6},${pad - 26})`);
        let lx = 0;
        types.forEach((t) => {
          const grp = lg.append("g").attr("transform", `translate(${lx},0)`);
          grp.append("circle").attr("r", 5).attr("fill", color(t));
          const txt = grp
            .append("text")
            .attr("x", 9)
            .attr("y", 4)
            .style("font-size", "10px")
            .style("fill", "#555")
            .text(t);
          lx += 26 + (txt.node().getComputedTextLength() || 30);
        });
        _scheduleExportEnhance();
      }

      // ── View 4: organism co-occurrence heatmap ───────────────────────────
      function _xsRenderCooc(rows, agg) {
        const topN = parseInt((document.getElementById("xs-topn") || {}).value || "25", 10);
        const metric = (document.getElementById("xs-metric") || {}).value || "TASS Score";
        const orgs = [...rows].sort((a, b) => b.sampleCount - a.sampleCount).slice(0, Math.min(40, topN));
        if (orgs.length < 2) {
          body_set('<p style="color:#888;padding:1em">Need at least 2 organisms to compute co-occurrence.</p>');
          return;
        }

        const isPresence = metric === "presence";
        const metricLabel = isPresence ? "Jaccard (presence)" : metric;
        const noteText = isPresence
          ? "How often each pair of organisms is detected in the same samples (Jaccard = shared ÷ union). Darker = more co-occurrence. Diagonal = organism's own prevalence."
          : `Cell colour = average <b>${metric}</b> across samples where both organisms were detected. Diagonal = mean ${metric} for that organism. Darker = higher value.`;

        body_set(
          `<div class="xs-note">${noteText}</div><div id="xs-cooc-svg" style="overflow:visible;padding-bottom:1em"></div>`,
        );

        const n = orgs.length;
        const J = Array.from({ length: n }, () => new Array(n).fill(0));
        const shared = Array.from({ length: n }, () => new Array(n).fill(0));

        // Build value matrix — Jaccard for presence, metric-weighted for TASS/Coverage
        let globalMax = 1;
        for (let i = 0; i < n; i++) {
          for (let j = i; j < n; j++) {
            const A = orgs[i].samples,
              B = orgs[j].samples;
            const mapI = metric === "Coverage" ? orgs[i].covMap : orgs[i].tassMap;
            const mapJ = metric === "Coverage" ? orgs[j].covMap : orgs[j].tassMap;
            let inter = 0;
            const sharedVals = [];
            A.forEach((s) => {
              if (B.has(s)) {
                inter++;
                if (!isPresence) sharedVals.push(((mapI.get(s) || 0) + (mapJ.get(s) || 0)) / 2);
              }
            });
            shared[i][j] = shared[j][i] = inter;

            let val;
            if (i === j) {
              if (isPresence) {
                val = A.size / agg.totalSamples;
              } else {
                const vals = [...mapI.values()];
                val = vals.length ? vals.reduce((s, v) => s + v, 0) / vals.length : 0;
              }
            } else {
              if (isPresence) {
                const uni = A.size + B.size - inter;
                val = uni ? inter / uni : 0;
              } else {
                val = sharedVals.length ? sharedVals.reduce((s, v) => s + v, 0) / sharedVals.length : 0;
              }
            }
            J[i][j] = val;
            J[j][i] = val;
            if (val > globalMax) globalMax = val;
          }
        }

        const colorDomain = isPresence ? [0, 1] : [0, globalMax];
        const color = d3.scaleSequential(d3.interpolateBlues).domain(colorDomain);
        const fmtVal = (v) => (isPresence ? v.toFixed(2) : v.toFixed(1));

        const wrap = document.getElementById("xs-cooc-svg");
        const W = wrap.clientWidth || 760;
        const mL = 210,
          mT = 210;
        const cell = Math.max(12, Math.min(26, Math.floor((W - mL - 20) / n)));
        const size = cell * n;
        // Extra right padding for -55° rotated column labels
        const labelPxEst = 26 * 5.5;
        const extraR = Math.ceil(labelPxEst * Math.cos((55 * Math.PI) / 180)) + 10;
        const svgW = mL + size + extraR;
        const Hsvg = size + mT + 20;

        const svg = d3
          .select(wrap)
          .append("svg")
          .attr("width", svgW)
          .attr("height", Hsvg)
          .attr("viewBox", `0 0 ${svgW} ${Hsvg}`)
          .style("overflow", "visible");

        const labels = orgs.map((o) => (o.name.length > 26 ? o.name.slice(0, 25) + "…" : o.name));

        // Column labels (rotated -55°)
        orgs.forEach((o, j) => {
          svg
            .append("text")
            .attr("transform", `translate(${mL + j * cell + cell / 2},${mT - 6}) rotate(-55)`)
            .attr("text-anchor", "start")
            .style("font-size", "9px")
            .style("font-style", "italic")
            .style("fill", "#455a64")
            .text(labels[j]);
        });

        orgs.forEach((o, i) => {
          // Row label
          svg
            .append("text")
            .attr("x", mL - 6)
            .attr("y", mT + i * cell + cell / 2)
            .attr("text-anchor", "end")
            .attr("dominant-baseline", "middle")
            .style("font-size", "9px")
            .style("font-style", "italic")
            .style("fill", "#455a64")
            .text(labels[i]);

          for (let j = 0; j < n; j++) {
            const v = J[i][j];
            const tipText =
              i === j
                ? isPresence
                  ? `<b>${orgs[i].name}</b><br>prevalence ${(v * 100).toFixed(0)}% (${orgs[i].sampleCount}/${
                      agg.totalSamples
                    } samples)`
                  : `<b>${orgs[i].name}</b><br>mean ${metricLabel}: ${fmtVal(v)} across ${
                      orgs[i].sampleCount
                    } sample(s)`
                : isPresence
                ? `<b>${orgs[i].name}</b> ✕ <b>${orgs[j].name}</b><br>${
                    shared[i][j]
                  } shared sample(s) · Jaccard ${fmtVal(v)}`
                : `<b>${orgs[i].name}</b> ✕ <b>${orgs[j].name}</b><br>${
                    shared[i][j]
                  } shared sample(s) · avg ${metricLabel}: ${fmtVal(v)}`;
            // Append the actual sample names behind the number so the detail lives
            // on hover. Diagonal = where this organism was detected; off-diagonal =
            // the samples the pair shares.
            const _capList = (arr, n = 14) => {
              const a = [...arr].sort();
              const s = a
                .slice(0, n)
                .map((x) => "• " + x)
                .join("<br>");
              return a.length > n ? s + `<br><span style="color:#999">+${a.length - n} more…</span>` : s;
            };
            let _coocExtra = "";
            if (i === j) {
              if (orgs[i].samples && orgs[i].samples.size)
                _coocExtra = `<br><span style="color:#9bb">Detected in:</span><br>${_capList(orgs[i].samples)}`;
            } else if (shared[i][j] > 0) {
              const _sh = [...orgs[i].samples].filter((s) => orgs[j].samples.has(s));
              _coocExtra = `<br><span style="color:#9bb">Shared samples:</span><br>${_capList(_sh)}`;
            }
            const _coocTip = tipText + _coocExtra;
            svg
              .append("rect")
              .attr("x", mL + j * cell)
              .attr("y", mT + i * cell)
              .attr("width", cell - 1)
              .attr("height", cell - 1)
              .attr("fill", v > 0 ? color(v) : "#f5f7fa")
              .attr("stroke", i === j ? "#37474f" : "none")
              .style("cursor", "pointer")
              .on("mouseover", (ev) => showTip(_coocTip, ev))
              .on("mousemove", moveTip)
              .on("mouseout", hideTip);
          }
        });

        // Colour-scale legend (right of grid)
        const cbW = 10,
          cbH = Math.min(140, size);
        const cbX = mL + size + 14,
          cbY = mT;
        const gradId = "xs-cooc-grad-" + Date.now();
        const cbGrad = svg
          .append("defs")
          .append("linearGradient")
          .attr("id", gradId)
          .attr("x1", "0%")
          .attr("y1", "100%")
          .attr("x2", "0%")
          .attr("y2", "0%");
        d3.range(0, 1.01, 0.1).forEach((t) =>
          cbGrad
            .append("stop")
            .attr("offset", `${t * 100}%`)
            .attr("stop-color", color(colorDomain[0] + t * (colorDomain[1] - colorDomain[0]))),
        );
        svg
          .append("rect")
          .attr("x", cbX)
          .attr("y", cbY)
          .attr("width", cbW)
          .attr("height", cbH)
          .style("fill", `url(#${gradId})`)
          .attr("rx", 2);
        svg
          .append("text")
          .attr("x", cbX + cbW / 2)
          .attr("y", cbY - 4)
          .attr("text-anchor", "middle")
          .style("font-size", "8px")
          .style("fill", "#607d8b")
          .text(isPresence ? "1" : colorDomain[1].toFixed(0));
        svg
          .append("text")
          .attr("x", cbX + cbW / 2)
          .attr("y", cbY + cbH + 11)
          .attr("text-anchor", "middle")
          .style("font-size", "8px")
          .style("fill", "#607d8b")
          .text("0");
        svg
          .append("text")
          .attr("transform", `translate(${cbX + cbW + 10},${cbY + cbH / 2}) rotate(90)`)
          .attr("text-anchor", "middle")
          .style("font-size", "8px")
          .style("fill", "#90a4ae")
          .text(metricLabel);
        _scheduleExportEnhance();
      }

      // ── CSV export of the frequency table (uses current filters/sort) ────
      function _xsExportCsv() {
        const agg = _XS.lastAgg;
        if (!agg) return;
        const rows = _xsFilterRows(agg);
        const k = _XS.sortKey,
          asc = _XS.sortAsc;
        rows.sort((a, b) => {
          let av = a[k],
            bv = b[k];
          if (typeof av === "string") return asc ? ("" + av).localeCompare("" + bv) : ("" + bv).localeCompare("" + av);
          return asc ? av - bv : bv - av;
        });
        const header = [
          "Detected Organism",
          "Pathogen Type",
          "TaxID",
          "sample_count",
          "total_samples",
          "sample_percent",
          "mean_tass",
          "median_tass",
          "max_tass",
          "mean_coverage",
          "median_coverage",
          "max_coverage",
          "total_reads",
          "high_consequence",
        ];
        const esc = (v) => {
          const s = "" + v;
          return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
        };
        const lines = [header.join(",")];
        rows.forEach((r) =>
          lines.push(
            [
              r.name,
              r.cat,
              r.taxid,
              r.sampleCount,
              agg.totalSamples,
              r.samplePct.toFixed(2),
              r.meanTass.toFixed(2),
              r.medianTass.toFixed(2),
              r.maxTass.toFixed(2),
              r.meanCov.toFixed(2),
              r.medianCov.toFixed(2),
              r.maxCov.toFixed(2),
              Math.round(r.reads),
              r.hc ? "yes" : "no",
            ]
              .map(esc)
              .join(","),
          ),
        );
        const blob = new Blob([lines.join("\n")], { type: "text/csv" });
        const a = document.createElement("a");
        a.href = URL.createObjectURL(blob);
        a.download = "organism_frequency_review.csv";
        a.click();
        setTimeout(() => URL.revokeObjectURL(a.href), 1000);
      }

      function _drawSummaryAnnotation(fd) {
        const wrap = document.getElementById("summary-annot-wrap");
        if (!wrap) return;
        // Show/hide the VF/AMR inner tab based on data availability
        const vfTab = document.querySelector('.sum-inner-tab[data-inner="vfamr"]');
        const vfPane = document.getElementById("sum-inner-vfamr");
        if (!HAS_PROT) {
          if (vfTab) vfTab.style.display = "none";
          if (vfPane) vfPane.classList.remove("active");
          return;
        }
        if (vfTab) vfTab.style.display = "";
        const genusRows = PROT.genus_summary || [];
        const geneHits = PROT.per_gene_hits || [];
        const amrRows = PROT.amr_genes || [];
        const annotOrgs = new Set(
          fd.filter((r) => isTruthy(r["IsAnnotated"])).map((r) => r["Taxonomic ID #"] || r["Detected Organism"]),
        ).size;

        const cards = [
          { label: "Annotated Organisms", value: _fmtInt(annotOrgs) },
          { label: "Gene Hits", value: _fmtInt(geneHits.length) },
          { label: "Genera w/ Hits", value: _fmtInt(new Set(genusRows.map((r) => r.Genus || r.genus)).size) },
          { label: "AMR Gene Hits", value: _fmtInt(amrRows.length) },
        ];
        document.getElementById("summary-annot-cards").innerHTML = cards
          .map(
            (c) =>
              `<div class="kpi-card" style="border-left-color:#2e7d32">` +
              `<div class="kpi-label">${c.label}</div>` +
              `<div class="kpi-value">${c.value}</div></div>`,
          )
          .join("");

        // Helper: best available gene/product label for an annotation row.
        const _geneLabel = (r) =>
          r["Gene"] ||
          r["Gene Name"] ||
          r.gene_name ||
          r["Product"] ||
          r.product ||
          r["Name"] ||
          r["Antibiotics"] ||
          "";

        // Helper: best available description/function for a gene row.
        const _geneDesc = (r) =>
          r["Description"] ||
          r.description ||
          r["Product"] ||
          r.product ||
          r["Antibiotics Class"] ||
          r["Antibiotics"] ||
          "";

        // Helper: category label for a gene row, with an optional fallback.
        const _geneCat = (r, fallback) =>
          r["Property"] || r["Class"] || r.property || r.class || r["Category"] || r.category || fallback || "";

        // Shared category → color map (mirrors _catColors in the VF/AMR tab).
        const _annotCatColors = {
          "Virulence Factor": "#e53935",
          "Antibiotic Resistance": "#fb8c00",
          AMR: "#fb8c00",
          "Drug Target": "#8e24aa",
          Transporter: "#00897b",
          Efflux: "#00897b",
        };

        // Top genera by # hits
        const byGenus = {};
        const _ensure = (g) => {
          if (!byGenus[g])
            // geneNames: Map<geneName, {samples: Set<sampleName>, desc: string}>  (unique genes, per-sample tracking + description)
            // amrNames:  Map<geneName, {samples: Set<sampleName>, desc: string}>
            byGenus[g] = { hits: 0, samples: new Set(), amr: 0, geneNames: new Map(), amrNames: new Map() };
          return byGenus[g];
        };
        genusRows.forEach((r) => {
          const g = r.Genus || r.genus || "Unknown";
          const n = num(r["# Hits"] || r["# hits"] || 0);
          const e = _ensure(g);
          e.hits += n;
          if (r.Sample || r.sample) e.samples.add(r.Sample || r.sample);
        });
        // Collect individual gene-hit names per genus — track which samples each gene was seen in.
        geneHits.forEach((r) => {
          const g = r.Genus || r.genus || "Unknown";
          const e = _ensure(g);
          const lbl = _geneLabel(r);
          const smp = r["Specimen ID"] || r.Sample || r.sample || null;
          if (lbl) {
            if (!e.geneNames.has(lbl))
              e.geneNames.set(lbl, { samples: new Set(), desc: _geneDesc(r), cat: _geneCat(r, "Virulence Factor") });
            if (smp) e.geneNames.get(lbl).samples.add(smp);
          }
          if (smp) e.samples.add(smp);
        });
        amrRows.forEach((r) => {
          const g = r.Genus || r.genus || "Unknown";
          const e = _ensure(g);
          e.amr += 1;
          const lbl = _geneLabel(r);
          const smp = r["Specimen ID"] || r.Sample || r.sample || null;
          if (lbl) {
            if (!e.amrNames.has(lbl))
              e.amrNames.set(lbl, { samples: new Set(), desc: _geneDesc(r), cat: _geneCat(r, "AMR") });
            if (smp) e.amrNames.get(lbl).samples.add(smp);
          }
          if (smp) e.samples.add(smp);
        });
        const ordered = Object.entries(byGenus)
          .sort((a, b) => b[1].geneNames.size - a[1].geneNames.size || b[1].hits - a[1].hits)
          .slice(0, 12);
        const cont = document.getElementById("summary-annot-genus");
        if (!ordered.length) {
          cont.innerHTML = '<p style="color:#888;padding:0.5em">No genus-level annotation hits.</p>';
          return;
        }
        let html =
          '<p style="font-size:0.74rem;color:#6c757d;margin:0 0 0.4em">Click a genus to open it in the VF / AMR tab.</p>' +
          "<table><thead><tr>" +
          ["Genus", "# Unique Genes", "AMR Hits", "# Samples", ""].map((h) => `<th>${h}</th>`).join("") +
          "</tr></thead><tbody>";
        // Tooltip HTML: geneNameMap is Map<name, {samples: Set, desc: string, cat: string}> — show gene names sorted by sample count with description and category badge.
        const _hitsTip = (title, geneNameMap) => {
          const items = Array.from(geneNameMap.entries())
            .map(([name, val]) => {
              const isSet = val instanceof Set;
              return [
                name,
                isSet ? val.size : val.samples.size,
                isSet ? "" : val.desc || "",
                isSet ? "" : val.cat || "",
              ];
            })
            .sort((a, b) => b[1] - a[1]);
          if (!items.length) return "";
          let t = `<b>${title}</b> <span style="color:#aaa;font-size:0.85em">(${items.length} unique gene${
            items.length === 1 ? "" : "s"
          })</span><br>`;
          t += `<table style="border-collapse:collapse;margin-top:4px;font-size:0.88em">`;
          t += `<tr style="color:#888"><td style="padding-right:12px">Gene</td><td style="text-align:right"># Samples</td></tr>`;
          items.slice(0, 15).forEach(([name, n, desc, cat]) => {
            const catColor = _annotCatColors[cat] || "#90a4ae";
            t += `<tr><td style="padding-right:12px;color:#e0e0e0;vertical-align:top"><i>${name}</i>`;
            if (desc && desc !== name) t += `<br><span style="font-size:0.8em;color:#aaa">${desc}</span>`;
            t += `</td><td style="vertical-align:top;white-space:nowrap">`;
            t += `<div style="display:flex;align-items:center;gap:6px;justify-content:space-between">`;
            if (cat)
              t += `<span style="font-size:0.72em;font-weight:600;color:${catColor};border:1px solid ${catColor};border-radius:3px;padding:0 4px;opacity:0.9;flex-shrink:0">${cat}</span>`;
            else t += `<span></span>`;
            t += `<span style="color:#90caf9;font-weight:600">${n}</span>`;
            t += `</div></td></tr>`;
          });
          if (items.length > 15)
            t += `<tr><td colspan="2" style="color:#aaa;padding-top:3px">…and ${items.length - 15} more</td></tr>`;
          t += `</table>`;
          return t;
        };
        // Sample list tooltip: show all sample names that had hits for this genus.
        const _samplesTip = (title, sampleSet) => {
          const smpArr = Array.from(sampleSet).sort();
          if (!smpArr.length) return "";
          let t = `<b>${title}</b> <span style="color:#aaa;font-size:0.85em">(${smpArr.length})</span><br>`;
          t += `<div style="margin-top:4px;max-height:220px;overflow:auto">`;
          smpArr.forEach((s) => {
            t += `<div style="color:#e0e0e0;padding:1px 0">${s}</div>`;
          });
          t += `</div>`;
          return t;
        };
        const _tipStore = {};
        ordered.forEach(([g, info], gi) => {
          const geneTip = _hitsTip(`${g} — unique gene hits`, info.geneNames);
          const amrTip = _hitsTip(`${g} — AMR genes`, info.amrNames);
          const samplesTip = _samplesTip(`${g} — samples with hits`, info.samples);
          _tipStore[gi] = { gene: geneTip, amr: amrTip, samples: samplesTip };
          const uniqueGenes = info.geneNames.size;
          html +=
            `<tr class="annot-genus-row" data-genus="${g}" data-gi="${gi}" style="cursor:pointer">` +
            `<td><i style="color:#1565c0;text-decoration:underline">${g}</i></td>` +
            `<td class="annot-gene-cell" style="text-align:right${
              geneTip ? ";cursor:help;text-decoration:underline dotted" : ""
            }">${uniqueGenes > 0 ? _fmtInt(uniqueGenes) : "—"}</td>` +
            `<td class="annot-amr-cell" style="text-align:right${
              amrTip ? ";cursor:help;text-decoration:underline dotted" : ""
            }">${info.amr ? _fmtInt(info.amr) : "—"}</td>` +
            `<td class="annot-samples-cell" style="text-align:right${
              samplesTip ? ";cursor:help;text-decoration:underline dotted" : ""
            }">${info.samples.size || "—"}</td>` +
            `<td style="text-align:right;color:#1565c0"><i class="fas fa-arrow-up-right-from-square" title="Open in VF / AMR tab"></i></td>` +
            "</tr>";
        });
        html += "</tbody></table>";
        cont.innerHTML = html;
        cont.querySelectorAll(".annot-genus-row").forEach((tr) => {
          const gi = tr.dataset.gi;
          const store = _tipStore[gi] || {};
          tr.addEventListener("click", () => _jumpToProteins(tr.dataset.genus));
          const wire = (cell, tip) => {
            if (!cell || !tip) return;
            cell.addEventListener("mouseover", (ev) => {
              ev.stopPropagation();
              showTip(tip, ev);
            });
            cell.addEventListener("mousemove", moveTip);
            cell.addEventListener("mouseout", hideTip);
          };
          wire(tr.querySelector(".annot-gene-cell"), store.gene);
          wire(tr.querySelector(".annot-amr-cell"), store.amr);
          wire(tr.querySelector(".annot-samples-cell"), store.samples);
        });
      }

      // Switch to the VF / AMR (proteins) tab and pre-filter it by genus.
      function _jumpToProteins(genus) {
        const btn = document.querySelector('.tab-btn[data-tab="proteins"]');
        if (!btn || btn.classList.contains("hidden")) return;
        btn.click();
        // After the pane is visible, set the protein search to the genus column.
        setTimeout(() => {
          const colSel = document.getElementById("prot-search-col");
          const search = document.getElementById("prot-search");
          if (colSel) {
            // pick a Genus column if present, else search all columns
            const opt = Array.from(colSel.options).find((o) => /genus/i.test(o.value));
            colSel.value = opt ? opt.value : "";
          }
          if (search) {
            search.value = genus;
            search.dispatchEvent(new Event("input", { bubbles: true }));
          }
        }, 60);
      }

      // ── NOVELTY DETECTION panel ───────────────────────────────────────────────
      // Reads BOOT.novelty ({samples: {<s>: {summary, candidates}}}) + BOOT.novelty_downloads.
      function _novEsc(v) {
        if (v == null) return "";
        return String(v).replace(
          /[&<>"']/g,
          (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c],
        );
      }
      function _novFmt(v) {
        if (v == null || v === "") return "";
        if (typeof v === "number") return Number.isInteger(v) ? String(v) : v.toFixed(4);
        return String(v);
      }
      const _NOV_SUMMARY_COLS = [
        ["sample", "Sample"],
        ["novelty_score", "Novelty score"],
        ["novelty_flag", "Flag"],
        ["dark_fraction", "Dark fraction"],
        ["highrank_only_fraction", "High-rank-only"],
        ["lowident_tail_mass", "Low-identity tail"],
        ["total_reads", "Residual reads"],
      ];
      const _NOV_SUMMARY_META = {
        novelty_score: {
          desc: "Composite score (higher = more divergent from known references). Typically 0–10+; samples above ~3 are flagged.",
          color: "#e8590c",
          range: [0, 10],
          pct: true,
        },
        dark_fraction: {
          desc: "Fraction of residual reads with no translated-search hit at any identity threshold. High values suggest truly unknown biology.",
          color: "#6741d9",
          range: [0, 1],
          pct: true,
        },
        highrank_only_fraction: {
          desc: "Fraction of reads where the best LCA could only be resolved to a high taxonomic rank (phylum/kingdom). Indicates distant but classifiable hits.",
          color: "#1c7ed6",
          range: [0, 1],
          pct: true,
        },
        lowident_tail_mass: {
          desc: "Mass of reads in the low-identity tail of the cosine-distance distribution. Elevated values point to divergent homologs.",
          color: "#2f9e44",
          range: [0, 0.5],
          pct: true,
        },
        total_reads: {
          desc: "Total residual reads entering novelty analysis (reads that did not align to any reference). For contig-based classifiers (kaiju, mmseqs2) these reads are assembled into contigs before classification.",
          color: "#1864ab",
          range: null,
          pct: false,
        },
      };
      const _NOV_CAND_COLS = [
        ["taxid", "Tax ID"],
        ["rank", "Rank"],
        ["name", "Name"],
        ["reads", _novUnitCap()],
        [
          "frac_of_sample",
          _novUnit() === "contigs"
            ? "Ctgs / residual reads"
            : _novUnit() === "genes"
            ? "Genes / residual reads"
            : "Frac of residual",
        ],
      ];
      let _novState = { sample: null, subtab: "overview" };
      let _novCandPage = 0;
      function _novCandPageSize() {
        const sel = document.getElementById("novelty-cand-page-size");
        return sel ? parseInt(sel.value) || 0 : 25;
      }
      // Wire candidate-table pager + page-size (static markup, so safe at parse time).
      (function _initNoveltyCandControls() {
        const go = (n) => {
          _novCandPage = n;
          _drawNoveltyCandidates();
        };
        const acts = {
          first: () => go(0),
          prev: () => go(_novCandPage - 1),
          next: () => go(_novCandPage + 1),
          last: () => {
            const rows = ((_novSamples()[_novState.sample] || {}).candidates || []).length;
            const ps = _novCandPageSize();
            go(ps > 0 ? Math.max(0, Math.ceil(rows / ps) - 1) : 0);
          },
        };
        Object.entries(acts).forEach(([k, fn]) => {
          const el = document.getElementById("novelty-cand-" + k);
          if (el) el.addEventListener("click", fn);
        });
        const ps = document.getElementById("novelty-cand-page-size");
        if (ps)
          ps.addEventListener("change", () => {
            _novCandPage = 0;
            _drawNoveltyCandidates();
          });
      })();

      // ── shared floating cell tooltip element ─────────────────────────────────
      function _novCellTip() {
        return document.getElementById("nov-cell-tip");
      }

      function _novShowCellTip(ev, label, value, meta, barFrac, ctx) {
        // ctx = optional { rank, nSamples, min, max, mean, runAvg } for sample-relative context
        const tip = _novCellTip();
        if (!tip) return;
        let html = `<div class="nct-label">${_novEsc(label)}</div>`;
        if (value !== null && value !== undefined) html += `<div class="nct-val">${_novEsc(String(value))}</div>`;
        if (barFrac !== null && barFrac !== undefined && isFinite(barFrac)) {
          const pct = Math.min(100, Math.max(0, barFrac * 100)).toFixed(1);
          const col = (meta && meta.color) || "#1c7ed6";
          html += `<div class="nct-bar-wrap"><div class="nct-bar" style="width:${pct}%;background:${col}"></div></div>`;
          html += `<div style="font-size:0.85em;color:#868e96;margin-bottom:2px">${pct}% of range</div>`;
        }
        if (ctx) {
          html += `<div style="margin-top:5px;padding-top:5px;border-top:1px solid #e9ecef;font-size:0.82em;color:#495057">`;
          if (ctx.rank != null && ctx.nSamples != null)
            html += `<div>Rank <b>${ctx.rank}</b> of ${ctx.nSamples} samples</div>`;
          if (ctx.mean != null)
            html += `<div>Run avg <b>${isFinite(ctx.mean) ? ctx.mean.toFixed(4).replace(/\.?0+$/, "") : "—"}</b></div>`;
          if (ctx.min != null && ctx.max != null)
            html += `<div style="color:#868e96">Range ${ctx.min.toFixed(4).replace(/\.?0+$/, "")} – ${ctx.max
              .toFixed(4)
              .replace(/\.?0+$/, "")}</div>`;
          if (ctx.runAvg != null) html += `<div>Run avg <b>${(ctx.runAvg * 100).toFixed(1)}%</b> for this bucket</div>`;
          html += `</div>`;
        }
        if (meta && meta.desc) html += `<div class="nct-desc">${_novEsc(meta.desc)}</div>`;
        tip.innerHTML = html;
        _novPositionTip(tip, ev);
        tip.style.display = "block";
      }
      function _novHideCellTip() {
        const tip = _novCellTip();
        if (tip) tip.style.display = "none";
      }
      function _novPositionTip(tip, ev) {
        const x = ev.clientX + 14,
          y = ev.clientY + 14;
        const vw = window.innerWidth,
          vh = window.innerHeight;
        tip.style.left = (x + 250 > vw ? x - 270 : x) + "px";
        tip.style.top = (y + 180 > vh ? y - 190 : y) + "px";
      }

      function _novSamples() {
        return (NOVELTY && NOVELTY.samples) || {};
      }

      function _drawNoveltySummary() {
        const samples = _novSamples();
        const names = Object.keys(samples).sort();
        const head = document.getElementById("novelty-summary-head");
        const body = document.getElementById("novelty-summary-body");
        if (!head || !body) return;
        // Compute per-column max for scaling bars
        const colMax = {};
        _NOV_SUMMARY_COLS.forEach(([k]) => {
          const vals = names.map((s) => +((samples[s] || {}).summary || {})[k]).filter(isFinite);
          colMax[k] = vals.length ? Math.max(...vals) : 1;
        });
        head.innerHTML = "<tr>" + _NOV_SUMMARY_COLS.map(([, lbl]) => `<th>${_novEsc(lbl)}</th>`).join("") + "</tr>";
        body.innerHTML = names
          .map((s) => {
            const sm = (samples[s] || {}).summary || { sample: s };
            const flagged = Number(sm.novelty_flag) === 1;
            return (
              `<tr${flagged ? ' style="background:#fff4e6"' : ""}>` +
              _NOV_SUMMARY_COLS
                .map(([k, lbl]) => {
                  if (k === "novelty_flag") {
                    return `<td style="text-align:center">${
                      flagged
                        ? '<span title="Flagged as novel" style="color:#e8590c;font-weight:600">● novel</span>'
                        : '<span style="color:#adb5bd">—</span>'
                    }</td>`;
                  }
                  const raw = sm[k];
                  const num = +raw;
                  const meta = _NOV_SUMMARY_META[k];
                  const displayed = _novFmt(raw);
                  // Inline mini-bar for numeric columns
                  let cell = `<td data-nov-col="${_novEsc(k)}" data-nov-label="${_novEsc(lbl)}" data-nov-val="${_novEsc(
                    displayed,
                  )}" style="cursor:default">`;
                  if (meta && isFinite(num) && num > 0) {
                    const range = meta.range || [0, colMax[k] || 1];
                    const frac = Math.min(1, (num - range[0]) / Math.max(1e-9, range[1] - range[0]));
                    const pct = (frac * 100).toFixed(0);
                    cell += `<div style="display:flex;align-items:center;gap:5px">`;
                    cell += `<span style="min-width:4.5em;text-align:right">${_novEsc(displayed)}</span>`;
                    cell += `<span style="flex:1;min-width:40px;max-width:90px;background:#f1f3f5;border-radius:3px;height:7px;overflow:hidden">`;
                    cell += `<span style="display:block;height:100%;width:${pct}%;background:${meta.color};border-radius:3px"></span>`;
                    cell += `</span></div>`;
                  } else {
                    cell += `<span>${_novEsc(displayed)}</span>`;
                  }
                  cell += `</td>`;
                  return cell;
                })
                .join("") +
              "</tr>"
            );
          })
          .join("");

        // Pre-compute run-wide stats per numeric column for tooltip context
        const colStats = {};
        _NOV_SUMMARY_COLS.forEach(([k]) => {
          const vals = names.map((s) => +((samples[s] || {}).summary || {})[k]).filter(isFinite);
          if (!vals.length) return;
          const sorted = [...vals].sort((a, b) => b - a); // descending for rank
          colStats[k] = {
            min: Math.min(...vals),
            max: Math.max(...vals),
            mean: vals.reduce((a, b) => a + b, 0) / vals.length,
            sorted,
            nSamples: vals.length,
          };
        });

        // Attach tooltip listeners to cells with data-nov-col
        body.querySelectorAll("td[data-nov-col]").forEach((td) => {
          td.addEventListener("mouseenter", (ev) => {
            const k = td.dataset.novCol,
              lbl = td.dataset.novLabel,
              valStr = td.dataset.novVal;
            const meta = _NOV_SUMMARY_META[k];
            const num = parseFloat(valStr);
            let frac = null;
            if (meta && isFinite(num)) {
              const range = meta.range || [0, 1];
              frac = Math.min(1, (num - range[0]) / Math.max(1e-9, range[1] - range[0]));
            }
            let ctx = null;
            const st = colStats[k];
            if (st && isFinite(num)) {
              const rank = st.sorted.findIndex((v) => v <= num) + 1 || st.nSamples;
              ctx = { rank, nSamples: st.nSamples, min: st.min, max: st.max, mean: st.mean };
            }
            _novShowCellTip(ev, lbl, valStr, meta, frac, ctx);
          });
          td.addEventListener("mousemove", (ev) => {
            const tip = _novCellTip();
            if (tip && tip.style.display !== "none") _novPositionTip(tip, ev);
          });
          td.addEventListener("mouseleave", _novHideCellTip);
        });
      }

      function _drawNoveltyCandidates() {
        const samples = _novSamples();
        const sel = document.getElementById("novelty-sample-sel");
        const head = document.getElementById("novelty-cand-head");
        const body = document.getElementById("novelty-cand-body");
        const countEl = document.getElementById("novelty-cand-count");
        if (!sel || !head || !body) return;
        const cur = _novState.sample || sel.value;
        const blocks = samples[cur] || {};
        const rows = blocks.candidates || [];
        // Pre-compute max reads and frac for bar scaling (over the full set so bars
        // stay comparable across pages, not just within the visible page).
        const maxReads = rows.length ? Math.max(1, ...rows.map((r) => +r.reads || 0)) : 1;
        const maxFrac = rows.length ? Math.max(1e-9, ...rows.map((r) => +r.frac_of_sample || 0)) : 1;
        head.innerHTML = "<tr>" + _NOV_CAND_COLS.map(([, lbl]) => `<th>${_novEsc(lbl)}</th>`).join("") + "</tr>";
        // ── pagination ──────────────────────────────────────────────────────
        const total = rows.length;
        const pageSize = _novCandPageSize();
        const pages = pageSize > 0 ? Math.max(1, Math.ceil(total / pageSize)) : 1;
        _novCandPage = Math.max(0, Math.min(_novCandPage, pages - 1));
        const start = pageSize > 0 ? _novCandPage * pageSize : 0;
        const end = pageSize > 0 ? Math.min(start + pageSize, total) : total;
        const pageRows = pageSize > 0 ? rows.slice(start, end) : rows;
        if (!rows.length) {
          body.innerHTML = `<tr><td colspan="${_NOV_CAND_COLS.length}" style="color:#888;text-align:center;padding:0.8em">No candidate taxa for this sample.</td></tr>`;
        } else {
          body.innerHTML = pageRows
            .map((r) => {
              const rankColor =
                {
                  species: "#2f9e44",
                  genus: "#1c7ed6",
                  family: "#f59f00",
                  order: "#9c36b5",
                  class: "#e03131",
                  phylum: "#1864ab",
                }[r.rank] || "#495057";
              return (
                "<tr>" +
                _NOV_CAND_COLS
                  .map(([k]) => {
                    if (k === "rank") {
                      return `<td><span style="color:${rankColor};font-weight:600">${_novEsc(
                        _novFmt(r[k]),
                      )}</span></td>`;
                    }
                    if (k === "name") {
                      return `<td style="font-style:italic">${_novEsc(_novFmt(r[k]))}</td>`;
                    }
                    if (k === "reads") {
                      const v = +r[k] || 0;
                      const pct = ((v / maxReads) * 100).toFixed(0);
                      return `<td data-cand-key="reads" data-cand-name="${_novEsc(
                        _novFmt(r.name),
                      )}" data-cand-rank="${_novEsc(
                        r.rank,
                      )}" data-cand-val="${v}" data-cand-pct="${pct}" style="cursor:default">
                  <div style="display:flex;align-items:center;gap:5px">
                    <span style="min-width:4em;text-align:right">${_novEsc(_novFmt(r[k]))}</span>
                    <span style="flex:1;min-width:40px;max-width:80px;background:#f1f3f5;border-radius:3px;height:7px;overflow:hidden">
                      <span style="display:block;height:100%;width:${pct}%;background:#1c7ed6;border-radius:3px"></span>
                    </span>
                  </div></td>`;
                    }
                    if (k === "frac_of_sample") {
                      const v = +r[k] || 0;
                      const pct = ((v / maxFrac) * 100).toFixed(0);
                      const dispPct = (v * 100).toFixed(3);
                      return `<td data-cand-key="frac" data-cand-name="${_novEsc(
                        _novFmt(r.name),
                      )}" data-cand-rank="${_novEsc(r.rank)}" data-cand-val="${_novEsc(
                        _novFmt(r[k]),
                      )}" data-cand-pct="${pct}" style="cursor:default">
                  <div style="display:flex;align-items:center;gap:5px">
                    <span style="min-width:5em;text-align:right">${dispPct}%</span>
                    <span style="flex:1;min-width:40px;max-width:80px;background:#f1f3f5;border-radius:3px;height:7px;overflow:hidden">
                      <span style="display:block;height:100%;width:${pct}%;background:#2f9e44;border-radius:3px"></span>
                    </span>
                  </div></td>`;
                    }
                    return `<td>${_novEsc(_novFmt(r[k]))}</td>`;
                  })
                  .join("") +
                "</tr>"
              );
            })
            .join("");
        }
        if (countEl) {
          countEl.textContent = total
            ? pageSize > 0 && pages > 1
              ? `${start + 1}–${end} of ${total} candidate${total === 1 ? "" : "s"}`
              : `${total} candidate${total === 1 ? "" : "s"}`
            : "";
        }

        // Pager controls (hidden when everything fits on one page).
        const pager = document.getElementById("novelty-cand-pager");
        if (pager) {
          const showPager = pageSize > 0 && pages > 1;
          pager.style.display = showPager ? "flex" : "none";
          if (showPager) {
            const lbl = document.getElementById("novelty-cand-page-label");
            if (lbl) lbl.textContent = `Page ${_novCandPage + 1} of ${pages}`;
            const first = document.getElementById("novelty-cand-first");
            const prev = document.getElementById("novelty-cand-prev");
            const next = document.getElementById("novelty-cand-next");
            const last = document.getElementById("novelty-cand-last");
            if (first) first.disabled = _novCandPage === 0;
            if (prev) prev.disabled = _novCandPage === 0;
            if (next) next.disabled = _novCandPage >= pages - 1;
            if (last) last.disabled = _novCandPage >= pages - 1;
          }
        }

        // Tooltip for candidate cells with inline bars
        body.querySelectorAll("td[data-cand-key]").forEach((td) => {
          td.addEventListener("mouseenter", (ev) => {
            const key = td.dataset.candKey,
              name = td.dataset.candName,
              rank = td.dataset.candRank;
            const val = td.dataset.candVal,
              pct = +td.dataset.candPct / 100;
            const isReads = key === "reads";
            const meta = {
              desc: `${isReads ? _novUnitCap() : "Fraction of sample reads"} assigned to ${rank} ${name}`,
              color: isReads ? "#1c7ed6" : "#2f9e44",
            };
            _novShowCellTip(ev, (isReads ? _novUnitCap() : "Frac of sample") + ` — ${name}`, val, meta, pct);
          });
          td.addEventListener("mousemove", (ev) => {
            const tip = _novCellTip();
            if (tip && tip.style.display !== "none") _novPositionTip(tip, ev);
          });
          td.addEventListener("mouseleave", _novHideCellTip);
        });
      }

      function _drawNoveltyDownloads() {
        const wrap = document.getElementById("novelty-downloads");
        if (!wrap) return;
        if (!NOVELTY_DL.length) {
          wrap.innerHTML = '<span style="color:#888">No downloadable files.</span>';
          return;
        }
        const icon = (k) => (k === "xlsx" ? "fa-file-excel" : "fa-file-code");
        wrap.innerHTML = NOVELTY_DL.map(
          (d) =>
            `<a href="${_novEsc(d.filename)}" download class="ncbi-link" ` +
            `style="display:inline-flex;align-items:center;gap:0.35em;padding:0.3em 0.6em;` +
            `border:1px solid #d0d7de;border-radius:6px;text-decoration:none;color:#0b5ed7" ` +
            `title="Download ${_novEsc(d.filename)}">` +
            `<i class="fas ${icon(d.kind)}"></i> ${_novEsc(d.label)} · ${_novEsc(d.kind.toUpperCase())}</a>`,
        ).join("");
      }

      // ── Method coverage + Genus comparison (TASS vs mmseqs rescue) ───────────
      let _novCmpMode = "coverage";
      let _novCmpToggleWired = false;
      let _novMcHidden = new Set(); // keys toggled off in method-coverage legend

      function _novWireCmpToggle() {
        if (_novCmpToggleWired) return;
        _novCmpToggleWired = true;
        document.querySelectorAll("#nov-cmp-toggle .nov-seg-btn").forEach((btn) => {
          btn.addEventListener("click", () => {
            _novCmpMode = btn.dataset.cmp;
            document
              .querySelectorAll("#nov-cmp-toggle .nov-seg-btn")
              .forEach((b) => b.classList.toggle("active", b === btn));
            const cov = document.getElementById("nov-cmp-coverage");
            const gen = document.getElementById("nov-cmp-genus");
            if (cov) cov.style.display = _novCmpMode === "coverage" ? "" : "none";
            if (gen) gen.style.display = _novCmpMode === "genus" ? "" : "none";
            if (_novCmpMode === "coverage") _drawNoveltyMethodCoverage();
            else _drawNoveltyGenusCompare();
          });
        });
      }

      // disjoint read-accounting buckets (sum ≈ 1). Older novelty JSON without the
      // extended count columns degrades to the fraction fields it does carry.
      // disjoint read-accounting buckets (sum ≈ 1). Built per-call so the classifier-specific
      // buckets ("<backend> → species", "<backend> → genus+") track the active --novelty mode.
      function _novMcSegs() {
        const cls = _novClsShort();
        const lc = _novClsLc();
        return [
          {
            key: "align",
            label: "Aligned to reference (TASS)",
            color: "#1c7ed6",
            desc: "Reads that mapped to a pulled reference — these are what the closed-set TASS scores are built from.",
          },
          {
            key: "k2only",
            label: "Kraken2-classified only",
            color: "#adb5bd",
            desc: `Classified by kraken2 but not aligned to a reference and not placed by ${lc} — known taxa with no reference in the run.`,
          },
          {
            key: "mmSp",
            label: `${cls} → species`,
            color: "#2f9e44",
            unit: _novUnit(),
            desc: `Residual (unaligned) ${_novUnit()} ${cls} placed at species-or-finer.`,
          },
          {
            key: "mmHr",
            label: `${cls} → genus+ (rescued)`,
            color: "#f59f00",
            unit: _novUnit(),
            desc: `Residual ${_novUnit()} ${lc} could only place above species (genus and higher) — the rescue bucket that the candidate table lists.`,
          },
          {
            key: "dark",
            label: "Dark matter",
            color: "#343a40",
            desc: `Explained by nothing: not aligned, not classified by kraken2, not placed by ${lc}.`,
          },
        ];
      }

      function _novMcFractions(sm) {
        const T = +sm.total_reads || 0;
        const R = +sm.ref_aligned,
          k2 = +sm.k2_classified;
        const Sp = +sm.mmseqs_assigned_species,
          Hr = +sm.mmseqs_assigned_highrank;
        const fAlign = isFinite(+sm.ref_aligned_frac) ? +sm.ref_aligned_frac : T ? (R || 0) / T : 0;
        const fHr = T && isFinite(Hr) ? Hr / T : +sm.highrank_only_fraction || 0;
        const fMm = isFinite(+sm.mmseqs_frac)
          ? +sm.mmseqs_frac
          : T && isFinite(Sp) && isFinite(Hr)
          ? (Sp + Hr) / T
          : fHr;
        const fSp = Math.max(0, fMm - fHr);
        const fK2 = isFinite(+sm.k2_frac) ? +sm.k2_frac : T ? (k2 || 0) / T : 0;
        const fK2only = Math.max(0, fK2 - fAlign);
        const fDark = Math.max(0, +sm.dark_fraction || 0);
        const clamp = (x) => Math.max(0, Math.min(1, x || 0));
        const reads = (f) => (T ? Math.round(f * T) : null);
        return {
          align: clamp(fAlign),
          k2only: clamp(fK2only),
          mmSp: clamp(fSp),
          mmHr: clamp(fHr),
          dark: clamp(fDark),
          _T: T,
          _reads: reads,
        };
      }

      function _drawNoveltyMethodCoverage() {
        const host = document.getElementById("nov-cmp-coverage");
        if (!host) return;
        const samples = _novSamples();
        const names = Object.keys(samples).sort();
        if (!names.length) {
          host.innerHTML = '<div style="color:#888;font-size:0.82em">No novelty samples.</div>';
          return;
        }
        const allSegs = _novMcSegs();
        const visibleSegs = allSegs.filter((s) => !_novMcHidden.has(s.key));
        const legend =
          '<div class="nov-mc-legend">' +
          allSegs
            .map(
              (s) =>
                `<span data-mc-key="${_novEsc(s.key)}" class="${_novMcHidden.has(s.key) ? "nov-mc-hidden" : ""}">` +
                `<span class="nov-mc-swatch" style="background:${s.color}"></span>${_novEsc(s.label)}</span>`,
            )
            .join("") +
          "</div>";
        const rows = names
          .map((s) => {
            const sm = (samples[s] || {}).summary || {};
            const fr = _novMcFractions(sm);
            const flagged = Number(sm.novelty_flag) === 1;
            // Sum visible fractions to rescale bars proportionally
            const visSum = visibleSegs.reduce((acc, seg) => acc + (fr[seg.key] || 0), 0) || 1;
            const segs = visibleSegs
              .map((seg) => {
                const f = fr[seg.key] || 0;
                if (f <= 0) return "";
                const rd = fr._reads(f);
                const rdTxt = rd != null ? `${rd.toLocaleString()} ${seg.unit || "reads"}` : "";
                const displayPct = ((f / visSum) * 100).toFixed(2);
                return `<div class="nov-mc-seg" style="width:${displayPct}%;background:${seg.color}"
              data-mc-label="${_novEsc(seg.label)}" data-mc-pct="${(f * 100).toFixed(1)}"
              data-mc-reads="${_novEsc(rdTxt)}" data-mc-color="${seg.color}"
              data-mc-desc="${_novEsc(seg.desc)}"></div>`;
              })
              .join("");
            const rescued = ((fr.mmSp + fr.mmHr) * 100).toFixed(1);
            return `<div class="nov-mc-row">
            <div class="nov-mc-name" title="${_novEsc(s)}">${
              flagged ? '<span style="color:#e8590c">●</span> ' : ""
            }${_novEsc(s)}</div>
            <div class="nov-mc-bar">${segs}</div>
            <div class="nov-mc-total" title="${_novEsc(_novClsShort())} rescued ${rescued}% of all reads">${
              fr._T ? fr._T.toLocaleString() : "—"
            } rds · +${rescued}%</div>
          </div>`;
          })
          .join("");
        host.innerHTML = legend + rows;

        // Wire up legend click-to-hide toggle
        host.querySelectorAll(".nov-mc-legend span[data-mc-key]").forEach((el) => {
          el.addEventListener("click", () => {
            const key = el.dataset.mcKey;
            if (_novMcHidden.has(key)) {
              _novMcHidden.delete(key);
            } else {
              // Don't allow hiding all segments
              const allKeys = allSegs.map((s) => s.key);
              if (_novMcHidden.size < allKeys.length - 1) _novMcHidden.add(key);
            }
            _drawNoveltyMethodCoverage();
          });
        });
        // Pre-compute run-wide average fraction per method bucket for tooltip context
        const mcRunAvg = {};
        _novMcSegs().forEach((seg) => {
          const fracs = names.map((s) => {
            const sm = (samples[s] || {}).summary || {};
            return _novMcFractions(sm)[seg.key] || 0;
          });
          mcRunAvg[seg.key] = fracs.length ? fracs.reduce((a, b) => a + b, 0) / fracs.length : 0;
        });

        host.querySelectorAll(".nov-mc-seg").forEach((el) => {
          el.addEventListener("mouseenter", (ev) => {
            const lbl = el.dataset.mcLabel,
              pct = +el.dataset.mcPct,
              reads = el.dataset.mcReads;
            // Find which bucket this segment belongs to for run-avg context
            const segKey = _novMcSegs().find((s) => s.label === lbl);
            const runAvg = segKey ? mcRunAvg[segKey.key] : null;
            _novShowCellTip(
              ev,
              lbl,
              `${pct}%${reads ? " · " + reads : ""}`,
              { color: el.dataset.mcColor, desc: el.dataset.mcDesc },
              pct / 100,
              runAvg != null ? { runAvg } : null,
            );
          });
          el.addEventListener("mousemove", (ev) => {
            const tip = _novCellTip();
            if (tip && tip.style.display !== "none") _novPositionTip(tip, ev);
          });
          el.addEventListener("mouseleave", _novHideCellTip);
        });
      }

      // closed-set genus rollup for one sample, from the main TASS records (BOOT.records)
      function _novClosedGenus(sample) {
        const out = {};
        (DATA || []).forEach((r) => {
          if ((r["Specimen ID"] || r.sample) !== sample) return;
          const g = (r["Genus Name"] || r["Genus"] || "").trim();
          if (!g) return;
          const key = g.toLowerCase();
          const o = out[key] || (out[key] = { name: g, tass: 0, cov: 0, reads: 0, pct: 0 });
          const gt = +r["Genus TASS"] || 0;
          if (gt > o.tass) o.tass = gt;
          // only sum strain-level rows so species/genus summary rows don't double count
          if ((r["Level"] || "Strain") === "Strain") {
            o.reads += +r["# Reads Aligned"] || 0;
            o.pct += +r["% Reads"] || 0;
            const cv = +r["Coverage"] || 0;
            if (cv > o.cov) o.cov = cv;
          }
        });
        return out;
      }

      function _drawNoveltyGenusCompare() {
        const host = document.getElementById("nov-cmp-genus");
        if (!host) return;
        const s = _novState.sample;
        const samples = _novSamples();
        const summary = (samples[s] || {}).summary || {};
        const cand = ((samples[s] || {}).candidates || []).filter((c) => (c.rank || "").toLowerCase() === "genus");
        const mm = {};
        cand.forEach((c) => {
          mm[(c.name || "").trim().toLowerCase()] = c;
        });
        const closed = _novClosedGenus(s);
        const keys = Array.from(new Set([...Object.keys(closed), ...Object.keys(mm)]));
        if (!keys.length) {
          host.innerHTML = `<div style="color:#888;font-size:0.82em;padding:0.4em 0">No genus-level closed-set rows or ${_novEsc(
            _novClsLc(),
          )} genus candidates for <b>${_novEsc(s || "—")}</b>.</div>`;
          return;
        }
        const maxReads = Math.max(1, ...cand.map((c) => +c.reads || 0));
        const rowsData = keys
          .map((k) => {
            const c = closed[k],
              m = mm[k];
            const tass = c ? c.tass : 0;
            const mFoh = m ? +m.frac_of_highrank || 0 : 0;
            const score = Math.max(tass / 100, mFoh);
            const bucket = c && m ? "both" : c ? "aligned_only" : "mmseqs_only";
            return { k, name: (c && c.name) || (m && m.name) || k, c, m, tass, mFoh, score, bucket };
          })
          .sort((a, b) => b.score - a.score);

        const _gcCls = _novClsShort();
        const _gcClsLc = _novClsLc();
        const _gcBucketMeta = {
          both: {
            label: "Both",
            short: "Both",
            color: "#1864ab",
            bg: "#e7f5ff",
            desc: `This genus is supported by both closed-set reference alignment (TASS) and ${_gcClsLc} assignment.`,
          },
          aligned_only: {
            label: "Aligned only",
            short: "Aligned only",
            color: "#d9480f",
            bg: "#fff4e6",
            desc: `This genus appears in reference-aligned closed-set detections but has no genus-level ${_gcClsLc} candidate in novelty output.`,
          },
          mmseqs_only: {
            label: `${_gcCls} only`,
            short: `${_gcCls} only`,
            color: "#2b8a3e",
            bg: "#ebfbee",
            desc: `This genus appears only in ${_gcClsLc} novelty candidates (reference set did not align it).`,
          },
          dark: {
            label: "Dark matter",
            short: "Dark matter",
            color: "#343a40",
            bg: "#f1f3f5",
            desc: `Residual reads not explained by reference alignment, kraken2 classification, or ${_gcClsLc} assignment; this bucket has no genus label.`,
          },
        };
        const bucketCounts = rowsData.reduce(
          (acc, rd) => {
            acc[rd.bucket] = (acc[rd.bucket] || 0) + 1;
            return acc;
          },
          { both: 0, aligned_only: 0, mmseqs_only: 0 },
        );
        const genusRows = Math.max(1, rowsData.length);
        const totalReads = +summary.total_reads || 0;
        const darkFrac = Math.max(0, +summary.dark_fraction || 0);
        const darkReads = totalReads ? Math.round(darkFrac * totalReads) : null;

        const _gcLegendChip = (key, label, tip, frac, extra = "") => {
          const m = _gcBucketMeta[key];
          return `<span class="nov-gc-chip" style="background:${m.bg};color:${m.color}" data-gc-badge="1"
            data-gc-label="${_novEsc(label)}" data-gc-tip="${_novEsc(tip)}" data-gc-frac="${(isFinite(frac)
              ? frac
              : 0
            ).toFixed(3)}" data-gc-color="${m.color}" data-gc-desc="${_novEsc(m.desc)}">${_novEsc(
              m.short,
            )}${extra}</span>`;
        };

        const legendHtml = `<div class="nov-gc-legend">
          ${_gcLegendChip(
            "dark",
            `Dark matter - ${_novEsc(s || "sample")}`,
            `${darkReads != null ? darkReads.toLocaleString() + " reads" : "No read count"} - ${(
              darkFrac * 100
            ).toFixed(1)}% of novelty input reads are unassigned to any method/genus.`,
            darkFrac,
          )}
          ${_gcLegendChip(
            "aligned_only",
            "Genus bucket - aligned only",
            `${bucketCounts.aligned_only} genus row(s) in this sample are closed-set only.`,
            bucketCounts.aligned_only / genusRows,
            `<span style="font-size:0.9em;opacity:0.9">(${bucketCounts.aligned_only})</span>`,
          )}
          ${_gcLegendChip(
            "mmseqs_only",
            `Genus bucket - ${_gcCls} only`,
            `${bucketCounts.mmseqs_only} genus row(s) in this sample are novelty-only (${_gcClsLc} rescue).`,
            bucketCounts.mmseqs_only / genusRows,
            `<span style="font-size:0.9em;opacity:0.9">(${bucketCounts.mmseqs_only})</span>`,
          )}
          ${_gcLegendChip(
            "both",
            "Genus bucket - both",
            `${bucketCounts.both} genus row(s) in this sample are supported by both aligned and ${_gcClsLc} evidence.`,
            bucketCounts.both / genusRows,
            `<span style="font-size:0.9em;opacity:0.9">(${bucketCounts.both})</span>`,
          )}
        </div>`;

        const srcTag = (rd) => {
          const m = _gcBucketMeta[rd.bucket] || _gcBucketMeta.mmseqs_only;
          const label = `Genus bucket - ${m.label} - ${rd.name}`;
          const tip =
            rd.bucket === "both"
              ? `Supported by both closed-set alignment and ${_gcClsLc} assignment.`
              : rd.bucket === "aligned_only"
              ? `Detected only by closed-set reference alignment (no ${_gcClsLc} genus candidate here).`
              : `Detected only by ${_gcClsLc} novelty assignment (no aligned genus row).`;
          return `<span class="nov-gc-chip" style="background:${m.bg};color:${m.color}" data-gc-badge="1"
            data-gc-label="${_novEsc(label)}" data-gc-tip="${_novEsc(tip)}"
            data-gc-frac="${(bucketCounts[rd.bucket] / genusRows).toFixed(3)}" data-gc-color="${m.color}"
            data-gc-desc="${_novEsc(m.desc)}">${_novEsc(m.short)}</span>`;
        };

        const tassCell = (rd) => {
          if (!rd.c) return '<td style="color:#adb5bd">—</td>';
          const tass = rd.tass,
            cov = rd.c.cov;
          const tcol = tass >= 70 ? "#2f9e44" : tass >= 40 ? "#f59f00" : "#e8590c";
          return `<td data-gc-tip="TASS ${tass.toFixed(1)} · best-member coverage ${cov.toFixed(1)}%"
            data-gc-frac="${(tass / 100).toFixed(3)}" data-gc-color="${tcol}" data-gc-label="Genus TASS — ${_novEsc(
              rd.name,
            )}" style="cursor:default;white-space:nowrap">
            <b style="color:${tcol}">${tass.toFixed(1)}</b>
            <span class="nov-gc-track"><span class="nov-gc-bar" style="width:${Math.min(100, tass).toFixed(
              0,
            )}%;background:${tcol}"></span></span>
            <span style="font-size:0.74em;color:#868e96">cov ${cov.toFixed(0)}%</span></td>`;
        };
        const pctCell = (rd) =>
          rd.c
            ? `<td style="white-space:nowrap;font-size:0.82em">${rd.c.pct.toFixed(3)}%</td>`
            : '<td style="color:#adb5bd">—</td>';
        const mmCell = (rd) => {
          if (!rd.m) return '<td style="color:#adb5bd">—</td>';
          const reads = +rd.m.reads || 0,
            w = ((reads / maxReads) * 100).toFixed(0);
          return `<td data-gc-tip="${reads} placed hit(s) · ${(rd.mFoh * 100).toFixed(1)}% of all genus+ placements"
            data-gc-frac="${rd.mFoh.toFixed(3)}" data-gc-color="#f59f00" data-gc-label="${_novEsc(
              _gcCls,
            )} hits — ${_novEsc(rd.name)}" style="cursor:default;white-space:nowrap">
            <b style="color:#e8590c">${reads}</b>
            <span class="nov-gc-track"><span class="nov-gc-bar" style="width:${w}%;background:#f59f00"></span></span>
            <span style="font-size:0.74em;color:#868e96">${(rd.mFoh * 100).toFixed(1)}%</span></td>`;
        };

        host.innerHTML = `
          <div style="font-size:0.74rem;color:#6c757d;margin:0 0 0.45em">
            Genus-level cross-check for <b>${_novEsc(
              s || "—",
            )}</b>: the closed-set <b>TASS</b> score / coverage / % reads (left)
            beside what the reference-free <b>${_novEsc(_gcCls)}</b> search places (right). Rows ${_novEsc(
              _gcClsLc,
            )} finds but the reference set
            never aligned are tagged <span style="color:#2b8a3e;font-weight:600">Rescued</span> — exactly the limited-reference / mock case.
          </div>
          ${legendHtml}
          <div style="overflow-x:auto">
          <table class="data-table" style="font-size:0.82em;width:100%">
            <thead><tr>
              <th>Genus</th><th>Bucket</th><th>Genus TASS - cov</th><th>% reads (aligned)</th><th>${_novEsc(
                _gcCls,
              )} hits - % placed</th>
            </tr></thead>
            <tbody>
            ${rowsData
              .map(
                (rd) => `<tr>
              <td style="font-style:italic;font-weight:600">${_novEsc(rd.name)}</td>
              <td>${srcTag(rd)}</td>
              ${tassCell(rd)}
              ${pctCell(rd)}
              ${mmCell(rd)}
            </tr>`,
              )
              .join("")}
            </tbody>
          </table></div>`;

        host.querySelectorAll("td[data-gc-tip]").forEach((td) => {
          td.addEventListener("mouseenter", (ev) => {
            _novShowCellTip(
              ev,
              td.dataset.gcLabel,
              td.dataset.gcTip,
              { color: td.dataset.gcColor },
              +td.dataset.gcFrac,
            );
          });
          td.addEventListener("mousemove", (ev) => {
            const tip = _novCellTip();
            if (tip && tip.style.display !== "none") _novPositionTip(tip, ev);
          });
          td.addEventListener("mouseleave", _novHideCellTip);
        });
        host.querySelectorAll("[data-gc-badge='1']").forEach((el) => {
          el.addEventListener("mouseenter", (ev) => {
            const frac = +el.dataset.gcFrac;
            _novShowCellTip(
              ev,
              el.dataset.gcLabel,
              el.dataset.gcTip,
              { color: el.dataset.gcColor, desc: el.dataset.gcDesc },
              isFinite(frac) ? frac : null,
            );
          });
          el.addEventListener("mousemove", (ev) => {
            const tip = _novCellTip();
            if (tip && tip.style.display !== "none") _novPositionTip(tip, ev);
          });
          el.addEventListener("mouseleave", _novHideCellTip);
        });
      }

      // ── Novelty score dot-chart (all samples) ───────────────────────────────
      function _drawNoveltyScoreChart() {
        const host = document.getElementById("nov-score-canvas");
        const resetBtn = document.getElementById("nov-score-reset-zoom");
        const zoomInBtn = document.getElementById("nov-score-zoom-in");
        const zoomOutBtn = document.getElementById("nov-score-zoom-out");
        if (!host) return;
        const samples = _novSamples();
        const names = Object.keys(samples).sort();
        if (!names.length) {
          host.innerHTML = "";
          host.style.display = "none";
          if (resetBtn) resetBtn.style.display = "none";
          if (zoomInBtn) zoomInBtn.style.display = "none";
          if (zoomOutBtn) zoomOutBtn.style.display = "none";
          return;
        }

        const data = names
          .map((s) => {
            const sm = (samples[s] || {}).summary || {};
            const score = +sm.novelty_score;
            const flagged = Number(sm.novelty_flag) === 1;
            const dark = +sm.dark_fraction || 0;
            const total = +sm.total_reads || 0;
            return { name: s, score: isFinite(score) ? score : null, flagged, dark, total };
          })
          .filter((d) => d.score !== null);

        if (!data.length) {
          host.innerHTML = "";
          host.style.display = "none";
          if (resetBtn) resetBtn.style.display = "none";
          if (zoomInBtn) zoomInBtn.style.display = "none";
          if (zoomOutBtn) zoomOutBtn.style.display = "none";
          return;
        }

        host.style.display = "block";
        if (zoomInBtn) zoomInBtn.style.display = "inline-block";
        if (zoomOutBtn) zoomOutBtn.style.display = "inline-block";
        if (host._novRo) {
          host._novRo.disconnect();
          host._novRo = null;
        }

        const THRESHOLD = 3;
        const fullXDomain = [-0.5, data.length - 0.5];
        let xDomain = [...fullXDomain];
        const globalMaxScore = Math.max(THRESHOLD * 1.5, d3.max(data, (d) => d.score) || THRESHOLD) * 1.05;
        const fullYDomain = [0, globalMaxScore];
        let yDomain = [...fullYDomain];
        const isFullXDomain = () =>
          Math.abs(xDomain[0] - fullXDomain[0]) < 1e-6 && Math.abs(xDomain[1] - fullXDomain[1]) < 1e-6;
        const isFullYDomain = () =>
          Math.abs(yDomain[0] - fullYDomain[0]) < 1e-6 && Math.abs(yDomain[1] - fullYDomain[1]) < 1e-6;
        const setZoomUi = () => {
          if (resetBtn) resetBtn.style.display = isFullXDomain() && isFullYDomain() ? "none" : "inline-block";
          if (zoomInBtn) zoomInBtn.disabled = data.length < 3;
          if (zoomOutBtn) zoomOutBtn.disabled = isFullXDomain() && isFullYDomain();
        };
        const applyXDomain = (next) => {
          if (!Array.isArray(next) || next.length !== 2) return false;
          const fullSpan = fullXDomain[1] - fullXDomain[0];
          const minSpan = Math.min(2, fullSpan);
          const span = Math.max(minSpan, Math.min(fullSpan, next[1] - next[0]));
          const center = (next[0] + next[1]) / 2;
          let lo = center - span / 2;
          let hi = center + span / 2;
          if (lo < fullXDomain[0]) {
            hi += fullXDomain[0] - lo;
            lo = fullXDomain[0];
          }
          if (hi > fullXDomain[1]) {
            lo -= hi - fullXDomain[1];
            hi = fullXDomain[1];
          }
          lo = Math.max(fullXDomain[0], lo);
          hi = Math.min(fullXDomain[1], hi);
          if (Math.abs(lo - xDomain[0]) < 1e-6 && Math.abs(hi - xDomain[1]) < 1e-6) return false;
          xDomain = [lo, hi];
          return true;
        };
        const applyYDomain = (next) => {
          if (!Array.isArray(next) || next.length !== 2) return false;
          const fullLo = fullYDomain[0];
          const fullHi = fullYDomain[1];
          const fullSpan = fullHi - fullLo;
          const minSpan = Math.max(fullSpan * 0.03, 0.1);
          const nextLo = Math.min(next[0], next[1]);
          const nextHi = Math.max(next[0], next[1]);
          const span = Math.max(minSpan, Math.min(fullSpan, nextHi - nextLo));
          const center = (nextLo + nextHi) / 2;
          let lo = center - span / 2;
          let hi = center + span / 2;
          if (lo < fullLo) {
            hi += fullLo - lo;
            lo = fullLo;
          }
          if (hi > fullHi) {
            lo -= hi - fullHi;
            hi = fullHi;
          }
          lo = Math.max(fullLo, lo);
          hi = Math.min(fullHi, hi);
          if (Math.abs(lo - yDomain[0]) < 1e-6 && Math.abs(hi - yDomain[1]) < 1e-6) return false;
          yDomain = [lo, hi];
          return true;
        };
        const trunc = (name, n = 24) => (name.length > n ? name.slice(0, n - 1) + "…" : name);
        const maxNameLen = data.reduce((m, d) => Math.max(m, d.name.length), 0);
        const labelPx = Math.min(maxNameLen, 38) * 6.2;
        const margin = {
          top: 14,
          right: 18,
          bottom: Math.max(56, Math.min(176, Math.ceil(labelPx * Math.SQRT1_2) + 24)),
          left: 50,
        };

        function _tooltipForDatum(ev, d) {
          const meta = _NOV_SUMMARY_META["novelty_score"];
          const frac = Math.min(1, d.score / 10);
          const ctx2 = {
            rank: [...data].sort((a, b) => b.score - a.score).findIndex((x) => x.name === d.name) + 1,
            nSamples: data.length,
            mean: data.reduce((a, b) => a + b.score, 0) / data.length,
            min: Math.min(...data.map((x) => x.score)),
            max: Math.max(...data.map((x) => x.score)),
          };
          _novShowCellTip(
            ev,
            d.name + (d.flagged ? " ● novel" : ""),
            `Score: ${d.score.toFixed(3)}  |  Dark: ${(d.dark * 100).toFixed(
              1,
            )}%  |  Reads: ${d.total.toLocaleString()}`,
            { color: d.flagged ? "#e8590c" : "#1c7ed6", desc: meta ? meta.desc : null },
            frac,
            ctx2,
          );
        }

        function paint() {
          const W = Math.max(host.clientWidth || 0, 380);
          const H = Math.max(250, 132 + margin.top + margin.bottom);
          const iW = Math.max(10, W - margin.left - margin.right);
          const iH = Math.max(10, H - margin.top - margin.bottom);
          host.innerHTML = "";

          const svg = d3
            .select(host)
            .append("svg")
            .attr("width", W)
            .attr("height", H)
            .style("display", "block")
            .style("overflow", "visible")
            .style("cursor", "crosshair");

          const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
          const x = d3.scaleLinear().domain(xDomain).range([0, iW]);
          const y = d3.scaleLinear().domain(yDomain).range([iH, 0]);

          g.append("rect").attr("width", iW).attr("height", iH).attr("fill", "#fff");
          g.append("g")
            .attr("class", "axis")
            .call(
              d3
                .axisLeft(y)
                .ticks(4)
                .tickSize(-iW)
                .tickFormat((v) => (+v).toFixed(1)),
            )
            .call((sel) => {
              sel.select(".domain").remove();
              sel.selectAll(".tick line").attr("stroke", "#e9ecef");
              sel.selectAll(".tick text").attr("fill", "#868e96").attr("font-size", 10);
            });

          if (THRESHOLD >= yDomain[0] && THRESHOLD <= yDomain[1]) {
            const thresholdY = y(THRESHOLD);
            g.append("line")
              .attr("x1", 0)
              .attr("x2", iW)
              .attr("y1", thresholdY)
              .attr("y2", thresholdY)
              .attr("stroke", "#f03e3e")
              .attr("stroke-width", 1.2)
              .attr("stroke-dasharray", "5,3");
            g.append("text")
              .attr("x", 4)
              .attr("y", thresholdY - 5)
              .attr("fill", "#f03e3e")
              .attr("font-size", 9)
              .text("flag threshold");
          }

          const clipId = `nov-score-clip-${Math.random().toString(36).slice(2, 10)}`;
          svg.append("defs").append("clipPath").attr("id", clipId).append("rect").attr("width", iW).attr("height", iH);

          const plot = g.append("g").attr("clip-path", `url(#${clipId})`);
          const pointData = data.map((d, i) => ({ ...d, i }));
          const barW = Math.max(3, Math.min(18, iW / Math.max(1, data.length) - 5));

          plot
            .selectAll("rect.nov-stem")
            .data(pointData)
            .join("rect")
            .attr("class", "nov-stem")
            .attr("x", (d) => x(d.i) - barW / 2)
            .attr("y", (d) => y(d.score))
            .attr("width", barW)
            .attr("height", (d) => Math.max(0, y(0) - y(d.score)))
            .attr("fill", (d) => (d.flagged ? "rgba(232,89,12,0.18)" : "rgba(28,126,214,0.15)"));

          plot
            .selectAll("circle.nov-dot")
            .data(pointData)
            .join("circle")
            .attr("class", "nov-dot")
            .attr("cx", (d) => x(d.i))
            .attr("cy", (d) => y(d.score))
            .attr("r", 4.8)
            .attr("fill", (d) => (d.flagged ? "#e8590c" : "#1c7ed6"))
            .attr("stroke", "#fff")
            .attr("stroke-width", 1.3);

          const lo = Math.max(0, Math.ceil(xDomain[0]));
          const hi = Math.min(data.length - 1, Math.floor(xDomain[1]));
          const visIdx = d3.range(lo, hi + 1);
          const maxTicks = Math.max(3, Math.floor(iW / 70));
          const step = Math.max(1, Math.ceil(visIdx.length / maxTicks));
          const tickValues = visIdx.filter((_, j) => j % step === 0);
          const xAxis = d3
            .axisBottom(x)
            .tickValues(tickValues)
            .tickFormat((i) => trunc(data[i]?.name || "", 24));

          g.append("g")
            .attr("class", "axis")
            .attr("transform", `translate(0,${iH})`)
            .call(xAxis)
            .call((sel) => {
              sel.select(".domain").attr("stroke", "#ced4da");
              sel.selectAll(".tick line").attr("stroke", "#ced4da");
              sel
                .selectAll(".tick text")
                .attr("fill", "#495057")
                .attr("font-size", 9)
                .attr("text-anchor", "end")
                .attr("dx", "-0.32em")
                .attr("dy", "0.42em")
                .attr("transform", "rotate(-45)");
            });

          svg
            .append("text")
            .attr("transform", `translate(11, ${margin.top + iH / 2}) rotate(-90)`)
            .attr("text-anchor", "middle")
            .attr("fill", "#6c757d")
            .attr("font-size", 10)
            .text("Novelty score");

          const brush = d3
            .brush()
            .extent([
              [0, 0],
              [iW, iH],
            ])
            .on("end", (ev) => {
              const s = ev.selection;
              if (!s) return;
              const x0 = x.invert(Math.min(s[0][0], s[1][0]));
              const x1 = x.invert(Math.max(s[0][0], s[1][0]));
              const y0 = y.invert(Math.max(s[0][1], s[1][1]));
              const y1 = y.invert(Math.min(s[0][1], s[1][1]));
              const xSpan = Math.abs(x1 - x0);
              const ySpan = Math.abs(y1 - y0);
              const nextX = [Math.max(-0.5, Math.floor(x0) - 0.5), Math.min(data.length - 0.5, Math.ceil(x1) + 0.5)];
              brushG.call(brush.move, null);
              const changedX = xSpan >= 0.45 ? applyXDomain(nextX) : false;
              const changedY = ySpan >= 0.02 ? applyYDomain([y0, y1]) : false;
              if (!changedX && !changedY) return;
              _novHideCellTip();
              setZoomUi();
              paint();
            });

          const brushG = g.append("g").attr("class", "nov-score-brush").call(brush);
          brushG.select(".selection").attr("fill", "#1971c233").attr("stroke", "#1971c2");
          brushG.select(".overlay").style("cursor", "crosshair");

          brushG
            .select(".overlay")
            .on("mousemove", (ev) => {
              const [mx, my] = d3.pointer(ev, g.node());
              if (mx < 0 || mx > iW || my < 0 || my > iH) {
                _novHideCellTip();
                return;
              }
              const idx = Math.round(x.invert(mx));
              if (idx < 0 || idx >= data.length) {
                _novHideCellTip();
                return;
              }
              _tooltipForDatum(ev, data[idx]);
            })
            .on("mouseleave", _novHideCellTip);
        }

        paint();
        setZoomUi();
        if (zoomInBtn) {
          zoomInBtn.onclick = () => {
            const xSpan = xDomain[1] - xDomain[0];
            const xCenter = (xDomain[0] + xDomain[1]) / 2;
            const ySpan = yDomain[1] - yDomain[0];
            const yCenter = (yDomain[0] + yDomain[1]) / 2;
            const changedX = applyXDomain([xCenter - xSpan * 0.35, xCenter + xSpan * 0.35]);
            const changedY = applyYDomain([yCenter - ySpan * 0.35, yCenter + ySpan * 0.35]);
            if (!changedX && !changedY) return;
            _novHideCellTip();
            setZoomUi();
            paint();
          };
        }
        if (zoomOutBtn) {
          zoomOutBtn.onclick = () => {
            const xSpan = xDomain[1] - xDomain[0];
            const xCenter = (xDomain[0] + xDomain[1]) / 2;
            const ySpan = yDomain[1] - yDomain[0];
            const yCenter = (yDomain[0] + yDomain[1]) / 2;
            const changedX = applyXDomain([xCenter - xSpan * 0.75, xCenter + xSpan * 0.75]);
            const changedY = applyYDomain([yCenter - ySpan * 0.75, yCenter + ySpan * 0.75]);
            if (!changedX && !changedY) return;
            _novHideCellTip();
            setZoomUi();
            paint();
          };
        }
        if (resetBtn) {
          resetBtn.onclick = () => {
            xDomain = [...fullXDomain];
            yDomain = [...fullYDomain];
            _novHideCellTip();
            setZoomUi();
            paint();
          };
        }

        const ro = new ResizeObserver(() => {
          requestAnimationFrame(paint);
        });
        ro.observe(host);
        host._novRo = ro;
      }

      // ── Novelty subtab switching ─────────────────────────────────────────────
      let _novSubtabsWired = false;
      function _novWireSubtabs() {
        if (_novSubtabsWired) return;
        _novSubtabsWired = true;
        document.querySelectorAll(".nov-subtab").forEach((btn) => {
          btn.addEventListener("click", () => {
            _novState.subtab = btn.dataset.nsub;
            document.querySelectorAll(".nov-subtab").forEach((b) => b.classList.toggle("active", b === btn));
            document.getElementById("novelty-sub-overview").style.display =
              _novState.subtab === "overview" ? "" : "none";
          });
        });
      }

      // Stamp the classifier badge + method-aware description onto the Novelty tab.
      function _drawNoveltyClassifier() {
        const info = _noveltyMethodInfo();
        const cls = _noveltyClassifier();
        const badge = document.getElementById("nov-classifier-badge");
        if (badge) {
          if (cls) {
            badge.dataset.cls = cls;
            badge.innerHTML = `<span class="nov-cls-dot"></span>Classifier: ${_novEsc(info.label)}`;
            badge.style.display = "inline-flex";
          } else {
            badge.style.display = "none";
          }
        }
        const desc = document.getElementById("novelty-overview-desc");
        if (desc && info.desc) {
          const queryDesc = _novGeneMode()
            ? "Pyrodigal-predicted genes from the de novo assembly"
            : "the de novo assembly";
          desc.innerHTML =
            `Reference-free novelty on ${queryDesc}. Backend: <b>${_novEsc(info.label)}</b> — ${_novEsc(info.desc)} ` +
            `Each sample gets a <b>novelty score</b> (higher = more divergent from known references) and a <b>flag</b> when ` +
            `that score crosses the run threshold. The candidate table lists the genus-or-higher taxa the classifier could ` +
            `still assign. Use the <b>downloads</b> for the raw per-sample / combined JSON and XLSX.`;
        }
        // Method-coverage / genus-compare section titles + the genus toggle label.
        const cmpTitle = document.getElementById("nov-cmp-title");
        if (cmpTitle && cls) {
          cmpTitle.textContent = `Where the reads landed — alignment / TASS vs ${info.short} rescue`;
        }
        const genusBtn = document.getElementById("nov-cmp-btn-genus");
        if (genusBtn && cls) {
          genusBtn.innerHTML = `<i class="fas fa-code-compare" style="margin-right: 0.3em"></i>Genus: TASS vs ${_novEsc(
            info.short,
          )}`;
        }
      }

      function drawNovelty() {
        const samples = _novSamples();
        const names = Object.keys(samples).sort();
        const sel = document.getElementById("novelty-sample-sel");
        if (sel) {
          // Prefer samples that actually have candidates; fall back to all.
          const withCand = names.filter((s) => (samples[s].candidates || []).length);
          const opts = withCand.length ? withCand : names;
          if (!_novState.sample || !samples[_novState.sample]) {
            _novState.sample = opts[0] || names[0] || null;
          }
          sel.innerHTML = opts.map((s) => `<option value="${_novEsc(s)}">${_novEsc(s)}</option>`).join("");
          if (_novState.sample) sel.value = _novState.sample;
          sel.onchange = () => {
            _novState.sample = sel.value;
            _novCandPage = 0;
            _drawNoveltyCandidates();
            _drawNoveltyGenusCompare();
          };
        }
        // Each sub-render is isolated so one failure can't blank the whole panel;
        // any error is surfaced in-place instead of silently leaving an empty tab.
        const _novSafe = (label, fn, targetId) => {
          try {
            fn();
          } catch (e) {
            console.error(`[novelty] ${label} failed:`, e);
            const el = targetId && document.getElementById(targetId);
            if (el)
              el.innerHTML =
                `<div style="color:#c92a2a;font-size:0.8em;padding:0.5em">` +
                `Could not render ${_novEsc(label)}: ${_novEsc((e && e.message) || String(e))}</div>`;
          }
        };
        _novSafe("sub-tabs", _novWireSubtabs);
        _novSafe("classifier badge", _drawNoveltyClassifier, "nov-classifier-badge");
        _novSafe("view toggle", _novWireCmpToggle);
        _novSafe("novelty score chart", _drawNoveltyScoreChart, "nov-score-canvas");
        _novSafe("per-sample summary", _drawNoveltySummary, "novelty-summary-body");
        _novSafe("candidate taxa", _drawNoveltyCandidates, "novelty-cand-body");
        _novSafe("method coverage", _drawNoveltyMethodCoverage, "nov-cmp-coverage");
        _novSafe("genus comparison", _drawNoveltyGenusCompare, "nov-cmp-genus");
        _novSafe("downloads", _drawNoveltyDownloads, "novelty-downloads");
      }

      function _drawTab(tab) {
        switch (tab) {
          case "summary":
            drawSummary();
            break;
          case "heatmap":
            drawHeatmap();
            break;
          case "tass":
            drawTassChart();
            break;
          case "sunburst":
            if (window.drawSunburst) window.drawSunburst();
            break;
          case "coverage":
            drawCoverage();
            break;
          case "proteins":
            if (HAS_PROT) drawProteins();
            break;
          case "histogram":
            if (window.drawHistogram) window.drawHistogram();
            break;
          case "explore":
            drawExplore();
            break;
          case "novelty":
            if (HAS_NOVELTY) drawNovelty();
            break;
          case "table":
            populateTable();
            break;
        }
        _TAB_DIRTY[tab] = false;
        _TAB_RENDERED[tab] = true;
        _scheduleExportEnhance();
      }
      function redraw() {
        _invalidateFilterCache();
        // The below-cutoff VF/AMR toggle + legends are only meaningful with
        // annotation data present.
        const _bvr = document.getElementById("below-vfamr-row");
        if (_bvr) _bvr.style.display = HAS_PROT ? "" : "none";
        // The sub-threshold / novelty toggle is meaningful whenever novelty data is
        // present (it also surfaces below-cutoff aligned rows regardless of novelty).
        const _nsr = document.getElementById("novelty-sub-row");
        if (_nsr) _nsr.style.display = HAS_NOVELTY ? "" : "none";
        ["sum-below-cutoff-legend", "tbl-below-cutoff-legend"].forEach((id) => {
          const el = document.getElementById(id);
          if (el) el.style.display = HAS_PROT ? "" : "none";
        });
        ["sum-novelty-only-legend", "tbl-novelty-only-legend"].forEach((id) => {
          const el = document.getElementById(id);
          if (el) el.style.display = HAS_NOVELTY ? "" : "none";
        });
        // Reset table to page 0 whenever filters change so the user
        // doesn't land on an empty page after narrowing results.
        _tblResetPage();
        // Mark every tab dirty so switches re-render lazily.
        for (const k in _TAB_DIRTY) _TAB_DIRTY[k] = true;
        // Render the currently visible tab now.
        _drawTab(activeTab);
        // Cheap shared updates that must reflect filter state immediately.
        _refreshMapMarkerColors();
        if (_selectedSample || _selectedGroup) _refreshMapPanelTable();
        if (_longiBuilt) _buildLongitudinalSection();
        // Keep the active metadata sub-tab in sync with the active filters.
        if (activeTab === "runmeta" && _activeMetaSub) {
          if (typeof _switchMetaSub === "function") _switchMetaSub(_activeMetaSub);
        }
        _updateCapabilityNotice();
      }

