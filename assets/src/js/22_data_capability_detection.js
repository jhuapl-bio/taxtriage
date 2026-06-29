      /* ═══════════════════════════════════════════════════════════════════════
       -  §  DATA CAPABILITY DETECTION  (graceful fallback for old files)
       -     Some views need data that older TSV / JSON / XLSX exports never
       -     carried — ANI annotation (Feature Compare) and per-contig depth
       -     histograms (Histograms tab). We detect, per sample, what's present
       -     and surface any gaps as (a) a consolidated notice in this sidebar
       -     and (b) a warning icon + hover tooltip on the affected tab.
       ═══════════════════════════════════════════════════════════════════════ */
      function _capabilities() {
        const data = DATA || [];
        const samples = uniq(data.map((r) => r["Specimen ID"] || "").filter(Boolean));
        // samples that have any per-contig depth/breadth histogram attached
        const contigSamples = new Set(
          (CONTIG_DATA || [])
            .filter(
              (cd) =>
                (cd.depth_histogram && Object.keys(cd.depth_histogram).length) ||
                (cd.breadth_histogram && (cd.breadth_histogram.bins || cd.breadth_histogram.b64)),
            )
            .map((cd) => cd.sample),
        );
        const aniBySample = {};
        samples.forEach((s) => (aniBySample[s] = false));
        data.forEach((r) => {
          const s = r["Specimen ID"] || "";
          if (s && _aniAnnotated(r)) aniBySample[s] = true;
        });
        const aniMissing = samples.filter((s) => !aniBySample[s]);
        const histMissing = samples.filter((s) => !contigSamples.has(s));
        return {
          samples,
          aniMissing,
          histMissing,
          nSamples: samples.length,
        };
      }

      function _capSampleList(list, max = 8) {
        const shown = list.slice(0, max).join(", ");
        return list.length > max ? `${shown} +${list.length - max} more` : shown;
      }

      function _updateCapabilityNotice() {
        const cap = _capabilities();
        const box = document.getElementById("data-capability-notice");
        const items = [];
        if (cap.nSamples && cap.aniMissing.length) {
          const all = cap.aniMissing.length === cap.nSamples;
          items.push(
            `<div style="margin-bottom:.35em"><b>ANI grouping ${all ? "unavailable" : "partial"}</b> — ` +
              `${
                all ? "no samples carry" : cap.aniMissing.length + " of " + cap.nSamples + " sample(s) lack"
              } ANI data, ` +
              `so the <i>Feature Compare</i> ANI column is ${all ? "hidden" : "limited"}. ` +
              `ANI comparison requires running the pipeline with <code>--enable_matrix</code> (these files were run without it or predate ANI support).` +
              `<div style="color:#8a6d00;font-size:.92em;margin-top:.15em">Out of date: ${_capSampleList(
                cap.aniMissing,
              )}</div></div>`,
          );
        }
        if (cap.nSamples && cap.histMissing.length) {
          const all = cap.histMissing.length === cap.nSamples;
          items.push(
            `<div><b>Depth histograms ${all ? "unavailable" : "partial"}</b> — ` +
              `${
                all ? "no samples carry" : cap.histMissing.length + " of " + cap.nSamples + " sample(s) lack"
              } per-contig depth data ` +
              `(typical for data loaded from XLSX / TSV), so the <i>Histograms</i> tab is ${
                all ? "hidden" : "limited"
              } for them.` +
              `<div style="color:#8a6d00;font-size:.92em;margin-top:.15em">Out of date: ${_capSampleList(
                cap.histMissing,
              )}</div></div>`,
          );
        }
        if (box) {
          if (items.length) {
            box.innerHTML =
              `<div style="border:1px solid #ffe0a3;background:#fff8e6;border-radius:6px;padding:.5em .6em;margin:.4em 0;font-size:.78em;color:#5f4b00;line-height:1.45">` +
              `<div style="font-weight:700;margin-bottom:.3em"><i class="fas fa-triangle-exclamation" style="color:#f59f00"></i> Some views have limited data</div>` +
              items.join("") +
              `</div>`;
            box.style.display = "";
          } else {
            box.innerHTML = "";
            box.style.display = "none";
          }
        }

        // ── per-tab warning icon: Feature Compare (ANI) ──────────────────────
        const cmpWarn = document.getElementById("xs-compare-warn");
        const cmpBtn = cmpWarn ? cmpWarn.closest(".xs-subtab") : null;
        if (cmpWarn && cmpBtn) {
          const aniGap = cap.nSamples && cap.aniMissing.length;
          cmpWarn.style.display = aniGap ? "" : "none";
          if (cmpBtn._capTipOver) cmpBtn.removeEventListener("mouseover", cmpBtn._capTipOver);
          if (cmpBtn._capTipMove) cmpBtn.removeEventListener("mousemove", cmpBtn._capTipMove);
          if (cmpBtn._capTipOut) cmpBtn.removeEventListener("mouseout", cmpBtn._capTipOut);
          if (aniGap) {
            const all = cap.aniMissing.length === cap.nSamples;
            const _msg =
              `<b><i class="fas fa-triangle-exclamation" style="color:#f59f00"></i> ANI grouping ${
                all ? "unavailable" : "limited"
              }</b><br>` +
              `<span style="font-size:0.88em;color:#ccc">${
                all ? "No sample in view carries" : cap.aniMissing.length + " sample(s) lack"
              } ANI annotation, ` +
              `so ANI-based grouping is ${
                all ? "disabled" : "partial"
              } here. The other features still compare normally.<br>` +
              `<b>Re-run the pipeline with <code style="color:#ffd580">--enable_matrix</code></b> to enable ANI comparison.</span>`;
            cmpBtn._capTipOver = (e) => showTip(_msg, e);
            cmpBtn._capTipMove = (e) => moveTip(e);
            cmpBtn._capTipOut = () => hideTip();
            cmpBtn.addEventListener("mouseover", cmpBtn._capTipOver);
            cmpBtn.addEventListener("mousemove", cmpBtn._capTipMove);
            cmpBtn.addEventListener("mouseout", cmpBtn._capTipOut);
          }
        }
      }

