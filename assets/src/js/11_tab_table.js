      /* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB: TABLE              (data-tab="table")
       -     buildTable() / populateTable() / renderTableHeaders().
       -     Flat tabular view of the filtered DATA — sortable headers, BSL
       -     level derivation (_computeBslLevels), supports CSV/XLSX export
       -     buttons defined in the tab UI.
═══════════════════════════════════════════════════════════════════════════ */
      /* ─── BSL level derivation ───────────────────────────────────────────────── */
      function _computeBslLevels() {
        if (!HAS_PROT) return;

        // Parse "['BSL-2']" / "['BSL-2','BSL-3']" / "BSL-3" → max numeric level
        const _parseBslMax = (str) => {
          if (!str) return 0;
          const nums = String(str).match(/BSL-(\d+)/gi);
          return nums ? Math.max(...nums.map((s) => parseInt(s.replace(/\D/g, "")))) : 0;
        };

        // Build species → maxBSL and genus → maxBSL from protein annotation rows
        const bySpecies = {},
          byGenus = {};
        [...(PROT.per_gene_hits || []), ...(PROT.amr_genes || [])].forEach((r) => {
          const lvl = _parseBslMax(r["Level"] || r["level"]);
          if (!lvl) return;
          const sp = (r["Species"] || r["species"] || "").trim().toLowerCase();
          const gn = (r["Genus"] || r["genus"] || "").trim().toLowerCase();
          if (sp) bySpecies[sp] = Math.max(bySpecies[sp] || 0, lvl);
          if (gn) byGenus[gn] = Math.max(byGenus[gn] || 0, lvl);
        });

        // Stamp each DATA record with its max BSL level
        DATA.forEach((rec) => {
          const org = (rec["Detected Organism"] || "").toLowerCase();
          const gn = (rec["Genus"] || "").toLowerCase();
          let max = 0;
          for (const [sp, lvl] of Object.entries(bySpecies)) {
            if (org.startsWith(sp) || org.includes(sp)) max = Math.max(max, lvl);
          }
          if (!max && gn && byGenus[gn]) max = byGenus[gn];
          rec["BSL Level"] = max ? `BSL-${max}` : "";
        });

        // Insert 'BSL Level' into ALL_COLS right after 'Detected Organism'
        if (!ALL_COLS.includes("BSL Level")) {
          const idx = ALL_COLS.indexOf("Detected Organism");
          idx >= 0 ? ALL_COLS.splice(idx + 1, 0, "BSL Level") : ALL_COLS.push("BSL Level");
        }
      }

      function buildTable() {
        visibleCols = [...ALL_COLS];
        // Append any user-added annotation columns (Notes etc.) so they render
        // as ordinary trailing columns. _ttApplyRowAnnotations mirrors their
        // values onto DATA rows so sorting and DOM export pick them up.
        if (TT_ANNOT.rowCols.length) {
          TT_ANNOT.rowCols.forEach((c) => {
            if (!visibleCols.includes(c)) visibleCols.push(c);
          });
          _ttApplyRowAnnotations();
        }
        // Hide "Specimen ID" column when "Group by Sample" is active — it's
        // redundant because rows are already separated by sample group headers.
        const _grpChecked = document.getElementById("tbl-group-sample")?.checked;
        if (_grpChecked) {
          visibleCols = visibleCols.filter((c) => c !== "Specimen ID");
        }
      }

      function renderTableHeaders() {
        const hr = document.getElementById("header-row");
        hr.innerHTML = "";
        visibleCols.forEach((c, i) => {
          const th = document.createElement("th");
          th.textContent = c;
          if (c === sortCol) th.classList.add(sortAsc ? "sort-asc" : "sort-desc");
          th.addEventListener("click", () => {
            if (sortCol === c) sortAsc = !sortAsc;
            else {
              sortCol = c;
              sortAsc = true;
            }
            _tblResetPage();
            renderTableHeaders();
            populateTable();
          });
          hr.appendChild(th);
        });
      }

