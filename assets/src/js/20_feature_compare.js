      /* ═══════════════════════════════════════════════════════════════════════
       -  §  FEATURE COMPARE — per-cell sample tooltips + per-sample popup
       -     Hovering any matrix cell shows the relevant sample breakdown for
       -     that organism × metric. Clicking a cell opens a popup table listing
       -     every sample in view with that organism's per-sample values.
       ═══════════════════════════════════════════════════════════════════════ */
      const _XS_SAMPLE_PALETTE = [
        "#1565c0",
        "#2e7d32",
        "#c62828",
        "#6a1b9a",
        "#ef6c00",
        "#00838f",
        "#ad1457",
        "#558b2f",
        "#4527a0",
        "#d84315",
        "#0277bd",
        "#9e9d24",
        "#5d4037",
        "#00695c",
        "#7b1fa2",
        "#827717",
      ];
      function _cmpSampleTip(o, m, agg) {
        const all = (agg.sampleList || []).slice();
        const det = new Set(o.samples || []);
        // Lead with the exact value for this cell — it's no longer printed on the
        // matrix, so the tooltip is now the primary place to read it.
        const _cellVal = m.fmt(+o[m.key] || 0);
        const head =
          `<b><i>${o.name}</i></b><br>` +
          `<span style="color:#9bb">${m.label}:</span> <b>${_cellVal}</b>` +
          `<span style="color:#9bb;font-size:.85em"> — ${m.desc}</span><br>`;
        const cap = (arr, n = 14) => {
          const s = arr
            .slice(0, n)
            .map((x) => "• " + x)
            .join("<br>");
          return arr.length > n ? s + `<br><span style="color:#999">+${arr.length - n} more…</span>` : s;
        };
        let body = "";
        if (m.key === "samplePct" || m.key === "sampleCount") {
          const yes = [...det].sort();
          const no = all.filter((s) => !det.has(s)).sort();
          if (m.key === "samplePct") {
            body =
              `<span style="color:#7CFC9B"><b>Detected (${yes.length}/${all.length}, ${(o.samplePct || 0).toFixed(
                0,
              )}%)</b></span><br>${cap(yes)}` +
              (no.length
                ? `<br><span style="color:#ffc078"><b>Not detected (${no.length})</b></span><br>${cap(no)}`
                : "");
          } else {
            body = `<b>Detected in ${yes.length} of ${all.length} sample(s):</b><br>${cap(yes)}`;
          }
        } else if (m.key === "reads") {
          const yes = [...det].sort();
          body = `<b>Total aligned reads:</b> ${_fmtInt(o.reads || 0)}<br><span style="color:#9bb">across ${
            yes.length
          } sample(s):</span><br>${cap(yes)}`;
        } else if (m.key === "aniGroupSize") {
          const members = (agg.aniGroups && agg.aniGroups.get(o.aniGroup)) || [o];
          body =
            `<b>ANI group — ${members.length} similar reference(s):</b><br>` +
            cap(members.map((r) => r.name + (r.taxid ? ` (${r.taxid})` : "")));
        } else {
          const map = /tass/i.test(m.key) ? o.tassMap : /cov/i.test(m.key) ? o.covMap : null;
          if (map && map.size) {
            const rowsv = [...map.entries()].sort((a, b) => b[1] - a[1]);
            body = `<b>Per-sample ${m.label}:</b><br>` + cap(rowsv.map(([s, v]) => `${s}: <b>${(+v).toFixed(1)}</b>`));
          } else {
            body = `<b>${m.label}:</b> ${m.fmt(+o[m.key] || 0)}`;
          }
        }
        return head + body + `<br><span style="color:#90caf9;font-size:.85em">click for full per-sample table</span>`;
      }

      function _xsOpenSamplePopup(o, m, agg) {
        const ov = document.getElementById("xspop-overlay");
        const body = document.getElementById("xspop-body");
        const titleEl = document.getElementById("xspop-title");
        if (!ov || !body) return;
        const all = (agg.sampleList || []).slice().sort();
        const det = new Set(o.samples || []);
        const tass = o.tassMap || new Map();
        const cov = o.covMap || new Map();
        const hlTass = /tass/i.test(m.key);
        const hlCov = /cov/i.test(m.key);
        const hlDet = m.key === "samplePct" || m.key === "sampleCount";
        if (titleEl) titleEl.innerHTML = `<i>${o.name}</i> — per-sample detail`;
        const hl = (on) => (on ? "background:#fff3bf" : "");
        const rowsHtml = all
          .map((s) => {
            const isDet = det.has(s);
            const tv = tass.has(s) ? (+tass.get(s)).toFixed(1) : "—";
            const cv = cov.has(s) ? (+cov.get(s)).toFixed(1) : "—";
            return (
              `<tr style="${isDet ? "" : "color:#aaa"}">` +
              `<td style="padding:4px 10px;border-bottom:1px solid #eee;text-align:left">${s}</td>` +
              `<td style="padding:4px 10px;border-bottom:1px solid #eee;text-align:center;${hl(hlDet)}">${
                isDet ? '<span style="color:#2e7d32;font-weight:700">✓</span>' : '<span style="color:#c0392b">✗</span>'
              }</td>` +
              `<td style="padding:4px 10px;border-bottom:1px solid #eee;text-align:right;${hl(hlTass)}">${tv}</td>` +
              `<td style="padding:4px 10px;border-bottom:1px solid #eee;text-align:right;${hl(hlCov)}">${cv}</td>` +
              `</tr>`
            );
          })
          .join("");
        const nDet = det.size,
          nTot = all.length;
        const stats =
          `<div style="font-size:.82em;color:#555;margin-bottom:.5em">` +
          `Detected in <b>${nDet}</b> / ${nTot} samples (${((nDet / Math.max(1, nTot)) * 100).toFixed(0)}%) · ` +
          `Mean TASS ${(o.meanTass || 0).toFixed(1)} · Max TASS ${(o.maxTass || 0).toFixed(1)} · ` +
          `Mean Cov ${(o.meanCov || 0).toFixed(1)} · Total reads ${_fmtInt(o.reads || 0)}</div>`;
        body.innerHTML =
          stats +
          `<div style="overflow:auto;max-height:60vh"><table style="border-collapse:collapse;width:100%;font-size:.85em">` +
          `<thead><tr style="position:sticky;top:0;background:#f0f4fa">` +
          `<th style="padding:5px 10px;text-align:left;border-bottom:2px solid #1565c0">Sample</th>` +
          `<th style="padding:5px 10px;text-align:center;border-bottom:2px solid #1565c0;${hl(hlDet)}">Detected</th>` +
          `<th style="padding:5px 10px;text-align:right;border-bottom:2px solid #1565c0;${hl(hlTass)}">TASS</th>` +
          `<th style="padding:5px 10px;text-align:right;border-bottom:2px solid #1565c0;${hl(hlCov)}">Coverage %</th>` +
          `</tr></thead><tbody>${rowsHtml}</tbody></table></div>`;
        ov.style.display = "flex";
      }

