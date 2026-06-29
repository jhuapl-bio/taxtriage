/* ═══════════════════════════════════════════════════════════════════════════
       -  §  PINNED ROW COMPARISON
       -     _updateComparePinnedBtn() — refreshes the badge count and
       -     visibility of the Compare button.
       -     _openCompareModal()      — renders and shows the comparison popup.
═══════════════════════════════════════════════════════════════════════════ */
function _updateComparePinnedBtn() {
  const bar = document.getElementById("tbl-pin-bar");
  const countEl = document.getElementById("tbl-pin-bar-count");
  if (!bar) return;
  const n = _tblPinned.size;
  bar.style.display = n > 0 ? "flex" : "none";
  if (countEl) countEl.textContent = `${n} row${n === 1 ? "" : "s"} pinned`;
}

function _openCompareModal() {
  // Gather all DATA rows that match tbl-tab pinned keys
  const pinnedRows = DATA.filter((r) => _tblPinned.has(_tblRowKey(r)));
  if (!pinnedRows.length) return;
  _openCompareModalForRows(pinnedRows);
}

let _cmpLastPinnedRows = [];
function _openCompareModalForRows(pinnedRows) {
  _cmpLastPinnedRows = pinnedRows || [];
  const METRICS = [
    { key: "Specimen ID", label: "Sample" },
    { key: "Detected Organism", label: "Organism" },
    { key: "Taxonomic ID #", label: "Taxon ID" },
    { key: "TASS Score", label: "TASS Score" },
    { key: "# Reads Aligned", label: "Reads Aligned" },
    { key: "Genome Coverage %", label: "Genome Coverage %" },
    { key: "Microbial Category", label: "Category" },
    { key: "Genus", label: "Genus" },
    { key: "BSL Level", label: "BSL" },
    { key: "High Consequence", label: "High Consequence" },
    { key: "Mol Type", label: "Mol Type" },
  ].filter((m) => pinnedRows.some((r) => r[m.key] != null && r[m.key] !== ""));

  // Build transposed table: rows = metrics, cols = pinned organisms
  const bslColor = { "BSL-1": "#2e7d32", "BSL-2": "#b45309", "BSL-3": "#c2410c", "BSL-4": "#b91c1c" };
  const bslBg = { "BSL-1": "#e8f5e9", "BSL-2": "#fff9c4", "BSL-3": "#ffe0b2", "BSL-4": "#ffcdd2" };

  const colHeaders = pinnedRows
    .map((r, i) => {
      const org = r["Detected Organism"] || r["Taxonomic ID #"] || `Row ${i + 1}`;
      const sn = r["Specimen ID"] || "";
      const isHC = isTruthy(r["High Consequence"]);
      return `<th style="text-align:left;padding:6px 10px;background:#f0f4fa;white-space:normal;word-break:break-word;min-width:140px;max-width:220px;border-bottom:2px solid #1565c0">
            ${isHC ? '<span style="color:#c62828;margin-right:3px">⚠</span>' : ""}
            <span style="font-weight:700;font-size:0.9em">${org}</span><br>
            <span style="font-weight:400;font-size:0.78em;color:#666">${sn}</span>
          </th>`;
    })
    .join("");

  const bodyRows = METRICS.map((m) => {
    const cells = pinnedRows
      .map((r) => {
        let raw = r[m.key] != null ? r[m.key] : "—";
        let style = "padding:6px 10px;white-space:nowrap;border-bottom:1px solid #eee";
        if (m.key === "TASS Score") {
          const sc = parseFloat(raw);
          if (!isNaN(sc) && sc > 0) {
            const t = Math.min(1, sc / 100);
            const hue = Math.round(210 - t * 65);
            const alpha = (0.1 + t * 0.3).toFixed(2);
            style += `;background:hsla(${hue},70%,45%,${alpha});font-weight:700;text-align:right`;
            raw = sc.toFixed(2);
          } else {
            style += ";text-align:right";
          }
        } else if (m.key === "# Reads Aligned") {
          const n = parseFloat(raw);
          if (!isNaN(n)) {
            style += ";text-align:right";
            raw = _fmtBig(n).short;
          }
        } else if (m.key === "BSL Level") {
          const bg = bslBg[raw] || "";
          const fg = bslColor[raw] || "";
          if (bg) style += `;background:${bg};color:${fg};font-weight:700;text-align:center`;
        } else if (m.key === "High Consequence") {
          const hc = isTruthy(raw);
          style += ";text-align:center";
          raw = hc ? '<span style="color:#c62828;font-weight:700">⚠ Yes</span>' : '<span style="color:#888">No</span>';
        } else if (m.key === "Genome Coverage %") {
          const n = parseFloat(raw);
          if (!isNaN(n)) {
            style += ";text-align:right";
            raw = n.toFixed(1) + "%";
          }
        }
        return `<td style="${style}">${raw}</td>`;
      })
      .join("");

    return `<tr>
            <td style="padding:6px 10px;font-size:0.8em;color:#6c757d;font-weight:600;white-space:nowrap;background:#fafafa;border-bottom:1px solid #eee;border-right:2px solid #e0e0e0;text-transform:uppercase;letter-spacing:0.02em">${m.label}</td>
            ${cells}
          </tr>`;
  }).join("");

  // Summary row: highest TASS per organism
  const maxTass = pinnedRows.map((r) => parseFloat(r["TASS Score"]) || 0);
  const maxIdx = maxTass.indexOf(Math.max(...maxTass));

  const body = document.getElementById("compare-body");
  body.innerHTML =
    `<p style="font-size:0.82em;color:#888;margin-bottom:0.6em">${pinnedRows.length} pinned row${
      pinnedRows.length !== 1 ? "s" : ""
    } — click any row in the table to toggle pinning</p>` +
    `<div style="overflow-x:auto"><table style="border-collapse:collapse;width:100%;font-size:0.85em">` +
    `<thead><tr>` +
    `<th style="text-align:left;padding:6px 10px;background:#f0f4fa;border-bottom:2px solid #1565c0;width:130px;font-size:0.78em;color:#555;text-transform:uppercase;letter-spacing:0.02em">Target</th>` +
    colHeaders +
    `</tr></thead><tbody>${bodyRows}</tbody></table></div>`;

  const overlay = document.getElementById("compare-overlay");
  overlay.style.display = "flex";
}

// Wire up compare button and modal close
(function _initCompare() {
  const btn = document.getElementById("compare-pinned-btn");
  if (btn) btn.addEventListener("click", _openCompareModal);

  // Toolbar "Unpin All" button
  const unpinBarBtn = document.getElementById("tbl-unpin-all-btn");
  if (unpinBarBtn)
    unpinBarBtn.addEventListener("click", () => {
      _tblPinned.clear();
      populateTable();
      _updateComparePinnedBtn();
    });

  const overlay = document.getElementById("compare-overlay");
  const closeModal = () => {
    overlay.style.display = "none";
  };
  document.getElementById("compare-close-btn").addEventListener("click", closeModal);
  document.getElementById("compare-close-btn2").addEventListener("click", closeModal);
  overlay.addEventListener("click", (e) => {
    if (e.target === overlay) closeModal();
  });
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape" && overlay.style.display !== "none") closeModal();
  });

  // Modal "Unpin All" button
  document.getElementById("compare-unpin-all-btn").addEventListener("click", () => {
    _tblPinned.clear();
    populateTable();
    _updateComparePinnedBtn();
    closeModal();
  });

  // "Compare Coverages" — overlay the pinned detections' depth histograms.
  const covBtn = document.getElementById("compare-cov-btn");
  if (covBtn)
    covBtn.addEventListener("click", () => {
      if (_cmpLastPinnedRows && _cmpLastPinnedRows.length) _openCoverageModalForRows(_cmpLastPinnedRows);
    });

  // Per-sample detail popup (Feature Compare cell click) close wiring.
  const xspop = document.getElementById("xspop-overlay");
  if (xspop) {
    const closeXs = () => (xspop.style.display = "none");
    const xb = document.getElementById("xspop-close");
    if (xb) xb.addEventListener("click", closeXs);
    xspop.addEventListener("click", (e) => {
      if (e.target === xspop) closeXs();
    });
    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape" && xspop.style.display !== "none") closeXs();
    });
  }

  // Compare Coverages overlay close wiring.
  const covOv = document.getElementById("cov-overlay");
  if (covOv) {
    const closeCov = () => (covOv.style.display = "none");
    ["cov-close", "cov-close2"].forEach((id) => {
      const b = document.getElementById(id);
      if (b) b.addEventListener("click", closeCov);
    });
    covOv.addEventListener("click", (e) => {
      if (e.target === covOv) closeCov();
    });
    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape" && covOv.style.display !== "none") closeCov();
    });
  }
})();
