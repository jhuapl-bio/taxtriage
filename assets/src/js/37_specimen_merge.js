/* ═══════════════════════════════════════════════════════════════════════════
       -  §  SPECIMEN MERGE  (right-panel control + drag-to-group modal)
       -     Lets an analyst treat several samples (e.g. a DNA and an RNA
       -     library from one swab) as a single biological specimen. Grouping
       -     defaults to the samplesheet `specimen` column (via SAMPLE_META →
       -     specimenOf) and can be edited live by dragging samples together.
       -
       -     • Right-panel bar  — a Merge On/Off toggle + a "Group…" button.
       -     • Group modal      — drag sample chips into specimen bins, name a
       -                          new specimen when you drop onto "+ New".
       -     • Aggregation      — reversible & view-level: filteredData()
       -                          collapses rows per specimen (reads summed,
       -                          TASS/coverage taken as the max). Nothing in
       -                          the underlying DATA array is mutated, so the
       -                          toggle restores per-sample rows instantly.
       -     • Metadata guard   — before applying, conflicting per-sample
       -                          metadata (lat/long, run, collection time, …)
       -                          is surfaced and the user picks which value
       -                          the merged specimen should keep.
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  "use strict";

  // Pipeline bookkeeping is not biological/run metadata and should never
  // trigger a merge warning. Every other populated RUN_META field is checked
  // dynamically, so schema-added fields get the same protection automatically.
  const _META_SKIP = new Set([
    "sample_name",
    "sample_type",
    "workflow_revision",
    "commit_id",
    "total_reads",
    "aligned_reads",
    "total_organism_reads",
    "num_species_groups",
    "num_keys",
    "num_subkeys",
    "num_toplevelkeys",
    "weights",
    "best_cutoffs",
    "best_cutoffs_source",
    "best_cutoffs_by_domain",
    "preferred_granularity",
  ]);
  const _META_FIRST = [
    "latitude",
    "longitude",
    "sample_origin_country",
    "sample_origin_state_province_territory",
    "location",
    "host_scientific_name",
    "host_disease",
    "environmental_site",
    "collection_time",
    "run_id",
    "platform",
    "depth",
    "salinity",
  ];

  /* ── data helpers ──────────────────────────────────────────────────── */

  // Every known sample id (union of SAMPLE_META, DATA and RUN_META).
  function _samples() {
    const s = new Set();
    if (typeof SAMPLE_META !== "undefined") Object.keys(SAMPLE_META || {}).forEach((k) => s.add(k));
    if (typeof DATA !== "undefined") DATA.forEach((r) => r["Specimen ID"] && s.add(r["Specimen ID"]));
    if (typeof RUN_META !== "undefined") RUN_META.forEach((r) => r.sample_name && s.add(r.sample_name));
    return Array.from(s).filter(Boolean);
  }

  // Merge whatever metadata we have for a sample from SAMPLE_META + RUN_META.
  function _metaFor(sample) {
    const sm = (typeof SAMPLE_META !== "undefined" && SAMPLE_META[sample]) || {};
    const rm = (typeof RUN_META !== "undefined" && RUN_META.find((r) => r.sample_name === sample)) || {};
    return Object.assign({}, sm, rm);
  }

  const _norm = (v) => (v == null ? "" : String(v).trim());
  const _metaLabel = (field) =>
    typeof _metaKeyLabel === "function"
      ? _metaKeyLabel(field)
      : field.replace(/_/g, " ").replace(/\b\w/g, (c) => c.toUpperCase());
  const _isMultiMeta = (field) =>
    typeof _MULTI_VALUE_META_FIELDS !== "undefined" && _MULTI_VALUE_META_FIELDS.includes(field);
  const _metaValues = (v, multi) => {
    if (v == null || v === "") return [];
    const vals = Array.isArray(v)
      ? v
      : multi && typeof _splitMultiMetaValue === "function"
      ? _splitMultiMetaValue(v)
      : [v];
    return (Array.isArray(vals) ? vals : [vals]).map((x) => String(x).trim()).filter(Boolean);
  };

  // A sample's molecule type, derived from its DATA rows: "dna" or "rna" when
  // unambiguous, "" when the sample mixes both, carries neither, or is unknown
  // (rendered as a white/empty indicator).
  function _sampleMolType(sample) {
    if (typeof DATA === "undefined") return "";
    const seen = new Set();
    for (let i = 0; i < DATA.length; i++) {
      if (DATA[i]["Specimen ID"] !== sample) continue;
      const mt = (DATA[i]["Mol Type"] || "").toLowerCase();
      if (mt === "dna" || mt === "rna") seen.add(mt);
    }
    return seen.size === 1 ? Array.from(seen)[0] : "";
  }
  // Colours for the mol-type dot on a chip. Both / none / unknown → white.
  const _MOLTYPE_DOT = { dna: "#1565c0", rna: "#e8590c" };
  function _molDotStyle(sample) {
    const mt = _sampleMolType(sample);
    const bg = _MOLTYPE_DOT[mt] || "#fff";
    const title = mt ? mt.toUpperCase() : "DNA + RNA / unspecified";
    return { bg, title, mt };
  }

  // Ensure specimen colours exist, then group the right-panel sample rows
  // under a draggable background-box header per specimen — instead of
  // appending a per-row name chip that had no width limit and could shove
  // the sample name off to the right. The header + its member rows read as
  // one enveloped box; dragging the header (its grip) moves every member of
  // that specimen together to a new position, resorting the whole group in
  // one drag instead of dragging each sample individually.
  function _annotateSidebar() {
    const cont = document.getElementById("sample-list");
    if (!cont) return;
    if (typeof _ensureSpecimenColors === "function") _ensureSpecimenColors();

    // buildSampleList() rebuilds every row on each call and always re-runs
    // this afterward, so start from a clean slate every time.
    cont.querySelectorAll(".specimen-group-header").forEach((h) => h.remove());
    cont.querySelectorAll(".sample-entry[data-sid]").forEach((row) => {
      row.style.background = "";
      row.style.borderLeft = "";
      row.style.borderRight = "";
      row.style.borderBottom = "";
      row.style.borderRadius = "";
      row.style.paddingLeft = "";
    });

    const on = typeof specimenMergeEnabled !== "undefined" && specimenMergeEnabled;
    if (!on) return;
    const groups = typeof specimenGroups === "function" ? specimenGroups() : new Map();

    // A merged specimen is one sortable unit. Keep all of its member rows
    // contiguous in both the shared order and the DOM, preserving the order in
    // which specimen groups first appear and the member order within each one.
    _coalesceSpecimenOrder();
    const rowsById = new Map(
      Array.from(cont.querySelectorAll(".sample-entry[data-sid]")).map((row) => [row.getAttribute("data-sid"), row]),
    );
    _sampleOrder.forEach((s) => {
      const row = rowsById.get(s);
      if (row) cont.appendChild(row);
    });

    // Walk rows in current DOM order; wherever a run of consecutive rows
    // shares a multi-member specimen, insert one header above the run and
    // tint the run so it reads as an enveloping box beneath it.
    const rows = Array.from(cont.querySelectorAll(".sample-entry[data-sid]"));
    let i = 0;
    while (i < rows.length) {
      const sid = rows[i].getAttribute("data-sid");
      const spec = typeof specimenOf === "function" ? specimenOf(sid) : sid;
      const grp = groups.get(spec);
      if (!spec || spec === sid || !grp || grp.length <= 1) {
        i++;
        continue;
      }
      const run = [rows[i]];
      let j = i + 1;
      while (
        j < rows.length &&
        (typeof specimenOf === "function" ? specimenOf(rows[j].getAttribute("data-sid")) : "") === spec
      ) {
        run.push(rows[j]);
        j++;
      }
      const col = (typeof sampleColors !== "undefined" && sampleColors[spec]) || "#1565c0";

      const header = document.createElement("div");
      header.className = "specimen-group-header";
      header.dataset.spec = spec;
      header.draggable = true;
      header.title = "Drag to reorder this specimen group";
      header.style.cssText =
        "display:flex;align-items:center;gap:5px;font-size:0.68em;font-weight:700;color:#fff;" +
        "padding:2px 6px;margin-top:6px;border-radius:4px 4px 0 0;background:" +
        col +
        ";border:2px solid " +
        col +
        ";border-bottom:none;box-shadow:0 2px 5px rgba(0,0,0,0.35);position:relative;z-index:1;cursor:grab;";
      header.innerHTML =
        `<i class="fas fa-grip-vertical" style="opacity:0.75;flex-shrink:0;"></i>` +
        `<span style="overflow:hidden;text-overflow:ellipsis;white-space:nowrap;flex:1;min-width:0;" title="Merged specimen: ${_esc(
          spec,
        )}">${_esc(spec)}</span>` +
        `<span style="opacity:0.85;flex-shrink:0;" title="${run.length} shown on this page">${grp.length} sample${
          grp.length > 1 ? "s" : ""
        }</span>`;
      cont.insertBefore(header, run[0]);

      // Color swatch — same mechanism as an individual sample row's color
      // picker (an <input type="color"> bound to sampleColors[id]), just
      // keyed by the specimen name so the whole envelope box can be
      // recolored right where it's shown, not only from the Group… modal.
      const colorIn = document.createElement("input");
      colorIn.type = "color";
      colorIn.className = "specimen-group-color";
      colorIn.value = col.length === 7 ? col : "#1565c0"; // input[type=color] needs #rrggbb
      colorIn.title = `Specimen color: ${spec}`;
      colorIn.style.cssText =
        "width:14px;height:14px;padding:0;border:1px solid rgba(255,255,255,0.7);border-radius:3px;flex-shrink:0;cursor:pointer;background:transparent;";
      colorIn.draggable = false;
      colorIn.addEventListener("mousedown", (e) => e.stopPropagation());
      colorIn.addEventListener("click", (e) => e.stopPropagation());
      colorIn.addEventListener("dragstart", (e) => {
        e.preventDefault();
        e.stopPropagation();
      });
      colorIn.addEventListener("input", (e) => {
        e.stopPropagation();
        const newCol = e.target.value;
        if (typeof sampleColors !== "undefined") sampleColors[spec] = newCol;
        header.style.background = newCol;
        header.style.border = "2px solid " + newCol;
        header.style.borderBottom = "none";
        run.forEach((row, idx) => {
          const isLast = idx === run.length - 1;
          row.style.background = newCol + "26";
          row.style.borderLeft = "4px solid " + newCol;
          row.style.borderRight = "2px solid " + newCol;
          row.style.borderBottom = isLast ? "2px solid " + newCol : "none";
        });
        if (typeof _refreshMapMarkerColors === "function") _refreshMapMarkerColors();
        if (typeof redraw === "function") redraw();
      });
      header.insertBefore(colorIn, header.firstChild.nextSibling); // after the grip icon

      // The envelope: a stronger tint + full left/right border across the
      // whole run (not just a left accent), so the box reads as a distinct
      // enclosure rather than blending into whatever sits above it. The last
      // row also gets the bottom edge closed off to complete the box.
      run.forEach((row, idx) => {
        const isLast = idx === run.length - 1;
        row.style.background = col + "26"; // ~15% tint — up from ~8%
        row.style.borderLeft = "4px solid " + col;
        row.style.borderRight = "2px solid " + col;
        row.style.borderBottom = isLast ? "2px solid " + col : "none";
        row.style.paddingLeft = "2px";
        row.style.borderRadius = isLast ? "0 0 4px 4px" : "0";
      });

      header.addEventListener("dragstart", (e) => {
        e.dataTransfer.effectAllowed = "move";
        e.dataTransfer.setData("application/x-specimen-group", spec);
        header.style.opacity = "0.5";
      });
      header.addEventListener("dragend", () => (header.style.opacity = ""));
      header.addEventListener("dragover", (e) => {
        e.preventDefault();
        e.dataTransfer.dropEffect = "move";
        // Give visual + positional feedback on which side of this header the
        // dragged group will land — dropping on the lower half moves it
        // *below* this specimen (needed to drag a group all the way down,
        // e.g. past the last group), the upper half moves it above.
        const after = _dropIsAfter(e, header);
        header.style.outline = "2px solid #fff";
        header.style.borderTop = after ? "" : "3px solid #fff";
        header.style.borderBottom = after ? "3px solid #fff" : "none";
      });
      header.addEventListener("dragleave", () => {
        header.style.outline = "";
        header.style.borderTop = "";
        header.style.borderBottom = "none";
      });
      header.addEventListener("drop", (e) => {
        e.preventDefault();
        header.style.outline = "";
        header.style.borderTop = "";
        header.style.borderBottom = "none";
        const draggedSpec = e.dataTransfer.getData("application/x-specimen-group");
        if (!draggedSpec || draggedSpec === spec) return;
        _reorderSpecimenGroup(draggedSpec, spec, groups, _dropIsAfter(e, header));
      });

      i = j;
    }

    // A merged header must also be droppable onto a singleton row; previously
    // groups could only be sorted relative to other merged headers, which made
    // sorting appear broken in runs containing singleton specimens.
    cont.querySelectorAll(".sample-entry[data-sid]").forEach((row) => {
      row.addEventListener("dragover", (e) => {
        if (!e.dataTransfer.types.includes("application/x-specimen-group")) return;
        e.preventDefault();
        e.dataTransfer.dropEffect = "move";
        const after = _dropIsAfter(e, row);
        row.style.borderTop = after ? "" : "3px solid #1565c0";
        row.style.borderBottom = after ? "3px solid #1565c0" : "";
      });
      row.addEventListener("dragleave", () => {
        row.style.borderTop = "";
        row.style.borderBottom = "";
      });
      row.addEventListener("drop", (e) => {
        const draggedSpec = e.dataTransfer.getData("application/x-specimen-group");
        if (!draggedSpec) return;
        e.preventDefault();
        e.stopImmediatePropagation();
        row.style.borderTop = "";
        row.style.borderBottom = "";
        const after = _dropIsAfter(e, row);
        const targetSample = row.getAttribute("data-sid");
        const targetSpec = typeof specimenOf === "function" ? specimenOf(targetSample) : targetSample;
        if (draggedSpec !== targetSpec) _reorderSpecimenGroup(draggedSpec, targetSpec, groups, after);
      });
    });
  }

  // True when the drag event's pointer sits in the lower half of `el` — used
  // to decide whether a dropped specimen group should land above or below
  // the row/header it was dropped on, so a group can be dragged to *any*
  // position (including all the way to the bottom, past the last group).
  function _dropIsAfter(e, el) {
    const rect = el.getBoundingClientRect();
    return e.clientY - rect.top > rect.height / 2;
  }

  function _coalesceSpecimenOrder() {
    if (typeof _sampleOrder === "undefined" || typeof specimenOf !== "function") return;
    const bySpec = new Map();
    _sampleOrder.forEach((s) => {
      const spec = specimenOf(s);
      if (!bySpec.has(spec)) bySpec.set(spec, []);
      bySpec.get(spec).push(s);
    });
    _sampleOrder = Array.from(bySpec.values()).flat();
  }

  // Move every sample belonging to `draggedSpec` to sit immediately before
  // (or, when `after` is true, immediately after) `targetSpec`'s samples in
  // the global sample display order, then rebuild the sidebar so the DOM
  // (and every downstream view that reads _sampleOrder) reflects the new
  // grouping. Honoring `after` matters: without it, a drop always inserted
  // *before* the target, so a group could never be moved past the last
  // specimen in the list — dragging "down" beyond the last group looked like
  // a no-op.
  function _reorderSpecimenGroup(draggedSpec, targetSpec, groups, after) {
    if (typeof _sampleOrder === "undefined") return;
    const draggedMembers = (groups.get(draggedSpec) || []).filter((s) => _sampleOrder.includes(s));
    if (!draggedMembers.length) return;
    const rest = _sampleOrder.filter((s) => !draggedMembers.includes(s));
    const isTarget = (s) => (typeof specimenOf === "function" ? specimenOf(s) : s) === targetSpec;
    let insertAt;
    if (after) {
      // Last index belonging to targetSpec, insert right after it.
      let lastIdx = -1;
      rest.forEach((s, i) => {
        if (isTarget(s)) lastIdx = i;
      });
      insertAt = lastIdx === -1 ? rest.length : lastIdx + 1;
    } else {
      const anchorIdx = rest.findIndex(isTarget);
      insertAt = anchorIdx === -1 ? rest.length : anchorIdx;
    }
    rest.splice(insertAt, 0, ...draggedMembers);
    _sampleOrder = rest;
    if (typeof buildSampleList === "function") buildSampleList();
    if (typeof redraw === "function") redraw();
  }

  /* ── right-panel bar ───────────────────────────────────────────────── */

  function refreshBar() {
    const bar = document.getElementById("specimen-merge-bar");
    if (!bar) return;
    const samples = _samples();
    // Show the bar whenever there is more than one sample to combine.
    bar.style.display = samples.length > 1 ? "flex" : "none";

    const grouped = typeof hasSpecimenGrouping === "function" && hasSpecimenGrouping();
    const tog = document.getElementById("specimen-merge-toggle");
    if (tog) {
      const on = typeof specimenMergeEnabled !== "undefined" && specimenMergeEnabled;
      tog.textContent = "Merge: " + (on ? "On" : "Off");
      tog.style.background = on ? "#1565c0" : "#fff";
      tog.style.color = on ? "#fff" : "#1565c0";
      tog.style.opacity = grouped ? "1" : "0.75";
      // Hover explanation is provided by the richer showTip tooltip wired in
      // wire() (_mergeToggleTip); no native title attribute needed here.
    }

    const cnt = document.getElementById("specimen-merge-count");
    if (cnt) {
      const groups = typeof specimenGroups === "function" ? specimenGroups() : new Map();
      const nSpec = groups.size;
      const nMerged = Array.from(groups.values()).filter((m) => m.length > 1).length;
      cnt.textContent = nMerged ? `${nSpec} specimen(s), ${nMerged} merged` : `${samples.length} sample(s)`;
    }
    if (typeof _ensureSpecimenColors === "function") _ensureSpecimenColors();
    _annotateSidebar();
  }

  function setEnabled(on) {
    if (typeof specimenMergeEnabled === "undefined") return;
    specimenMergeEnabled = !!on;
    if (typeof _ensureSpecimenColors === "function") _ensureSpecimenColors();
    if (typeof _invalidateFilterCache === "function") _invalidateFilterCache();
    refreshBar();
    const bannerSub = document.getElementById("banner-sub");
    if (bannerSub && typeof _buildBannerSub === "function") bannerSub.textContent = _buildBannerSub();
    if (typeof _rebuildMapMarkers === "function") _rebuildMapMarkers();
    if (typeof _buildRunMetaTable === "function") _buildRunMetaTable();
    if (typeof redraw === "function") redraw();
  }

  /* ── metadata conflict detection ───────────────────────────────────── */

  // For a proposed assignment (sample → specimen), find fields whose value
  // differs across the members of any multi-sample specimen. Single-valued
  // fields require one selected value; declared multiselect fields are safely
  // union-merged and shown as an informational notice.
  function findConflicts(assign) {
    const bySpec = new Map();
    Object.keys(assign).forEach((s) => {
      const g = assign[s];
      if (!bySpec.has(g)) bySpec.set(g, []);
      bySpec.get(g).push(s);
    });
    const conflicts = [];
    bySpec.forEach((members, spec) => {
      if (members.length < 2) return;
      const fieldSet = new Set(_META_FIRST);
      members.forEach((s) => {
        const rm = (typeof RUN_META !== "undefined" && RUN_META.find((r) => r.sample_name === s)) || {};
        Object.keys(rm).forEach((field) => {
          if (!_META_SKIP.has(field)) fieldSet.add(field);
        });
      });
      const fields = Array.from(fieldSet).sort((a, b) => {
        const ai = _META_FIRST.indexOf(a),
          bi = _META_FIRST.indexOf(b);
        if (ai !== -1 || bi !== -1) return (ai === -1 ? 999 : ai) - (bi === -1 ? 999 : bi);
        return a.localeCompare(b);
      });
      fields.forEach((field) => {
        const multi = _isMultiMeta(field);
        const seen = new Map(); // value → samples carrying it
        members.forEach((s) => {
          _metaValues(_metaFor(s)[field], multi).forEach((v) => {
            if (!seen.has(v)) seen.set(v, []);
            seen.get(v).push(s);
          });
        });
        if (multi && seen.size) {
          conflicts.push({
            spec,
            field,
            label: _metaLabel(field),
            kind: "multi",
            merged: Array.from(seen.keys()).sort((a, b) => a.localeCompare(b)),
            options: Array.from(seen.entries()).map(([value, samples]) => ({ value, samples })),
          });
        } else if (seen.size > 1) {
          conflicts.push({
            spec,
            field,
            label: _metaLabel(field),
            kind: "single",
            options: Array.from(seen.entries()).map(([value, samples]) => ({ value, samples })),
          });
        }
      });
    });
    return conflicts;
  }

  /* ── modal plumbing ────────────────────────────────────────────────── */

  function _closeModal() {
    const m = document.getElementById("specimen-merge-modal");
    if (m) m.remove();
  }

  function _esc(s) {
    return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
  }

  // Build the initial staged assignment from the live grouping. Samples that
  // are their own singleton specimen start "unassigned"; anything already in a
  // multi-sample (or renamed) specimen starts inside that bin.
  function _initialAssign() {
    const assign = {};
    _samples().forEach((s) => {
      const g = typeof specimenOf === "function" ? specimenOf(s) : s;
      if (g && g !== s) assign[s] = g;
    });
    // Also honour metadata groups where >1 sample share a name even == sample.
    if (typeof specimenGroups === "function") {
      specimenGroups().forEach((members, g) => {
        if (members.length > 1) members.forEach((s) => (assign[s] = g));
      });
    }
    return assign;
  }

  function openModal(opts) {
    opts = opts || {};
    _closeModal();
    const assign = opts.assign ? Object.assign({}, opts.assign) : _initialAssign(); // staged; committed on Apply
    // Display order of specimen boxes — user-draggable, preserved across
    // re-renders; new groups append at the end, removed groups drop out.
    let _binOrder = [];
    let _selectedSpec = null;

    const back = document.createElement("div");
    back.id = "specimen-merge-modal";
    back.style.cssText =
      "position:fixed;inset:0;background:rgba(0,0,0,0.45);z-index:10001;display:flex;align-items:center;justify-content:center;";

    back.innerHTML =
      `<div style="background:#fff;border-radius:8px;max-width:760px;width:94%;height:94vh;max-height:820px;display:flex;flex-direction:column;overflow:hidden;box-shadow:0 8px 30px rgba(0,0,0,0.3);">` +
      `<div style="padding:14px 18px;border-bottom:1px solid #eee;display:flex;align-items:center;gap:8px;">` +
      `<b style="font-size:1.02em;color:#1565c0;"><i class="fas fa-object-group"></i> ${_esc(
        opts.title || "Group samples into specimens",
      )}</b>` +
      `<span style="flex:1"></span>` +
      `<span style="cursor:pointer;font-size:1.3em;color:#888;" id="sm-close">&times;</span>` +
      `</div>` +
      `<div style="padding:6px 18px;font-size:0.82em;color:#666;">${
        opts.message ||
        "Drag samples into a specimen to combine them (e.g. a DNA and RNA library from one swab). Reads are summed and the TASS score becomes the strongest across the group. Drop onto <b>+ New specimen</b> to create and name one."
      }` +
      `<span style="display:inline-flex;align-items:center;gap:10px;margin-left:6px;">` +
      `<span style="display:inline-flex;align-items:center;gap:3px;"><span style="width:9px;height:9px;border-radius:50%;background:#1565c0;display:inline-block;"></span>DNA</span>` +
      `<span style="display:inline-flex;align-items:center;gap:3px;"><span style="width:9px;height:9px;border-radius:50%;background:#e8590c;display:inline-block;"></span>RNA</span>` +
      `<span style="display:inline-flex;align-items:center;gap:3px;"><span style="width:9px;height:9px;border-radius:50%;background:#fff;border:1px solid #b8c8e0;display:inline-block;"></span>both / unspecified</span>` +
      `</span></div>` +
      `<div style="display:flex;align-items:center;gap:7px;padding:0 18px 4px;">` +
      `<div style="position:relative;flex:1;min-width:0;"><i class="fas fa-magnifying-glass" aria-hidden="true" style="position:absolute;left:8px;top:50%;transform:translateY(-50%);font-size:0.72em;color:#90a4ae;"></i>` +
      `<input id="sm-sample-search" type="search" placeholder="Search samples or specimen groups…" aria-label="Search samples or specimen groups" style="box-sizing:border-box;width:100%;padding:5px 8px 5px 25px;border:1px solid #c5d3e3;border-radius:5px;font-size:0.78em;"></div>` +
      `<span id="sm-search-count" style="font-size:0.7em;color:#78909c;white-space:nowrap;"></span></div>` +
      `<div style="display:flex;gap:12px;padding:8px 18px 4px;overflow:auto;flex:1 1 auto;min-height:150px;">` +
      `<div style="flex:1;min-width:210px;display:flex;flex-direction:column;">` +
      `<div style="font-size:0.75em;color:#888;font-weight:600;margin-bottom:4px;">Unassigned samples</div>` +
      `<div id="sm-pool" class="sm-drop" data-bin="" style="flex:1;min-height:120px;border:1px dashed #cbd5e6;border-radius:6px;padding:6px;overflow:auto;background:#fafcff;"></div>` +
      `</div>` +
      `<div style="flex:1.2;min-width:240px;display:flex;flex-direction:column;">` +
      `<div style="font-size:0.75em;color:#888;font-weight:600;margin-bottom:4px;">Specimens <span style="font-weight:400;color:#aab;">(drag a box by its <i class="fas fa-grip-vertical"></i> handle to reorder)</span></div>` +
      `<div id="sm-bins" style="flex:1;overflow:auto;"></div>` +
      `<div id="sm-newbin" class="sm-drop" data-bin="__new__" style="margin-top:6px;border:1px dashed #1565c0;border-radius:6px;padding:8px;text-align:center;color:#1565c0;font-size:0.8em;font-weight:600;background:#f0f7ff;cursor:copy;">+ New specimen <span style="color:#88a;">(drop here to create)</span></div>` +
      `</div>` +
      `</div>` +
      `<div id="sm-conflicts" style="padding:0 18px;max-height:28vh;overflow:auto;flex:0 1 auto;"></div>` +
      `<div style="padding:12px 18px;border-top:1px solid #eee;display:flex;gap:8px;align-items:center;">` +
      `<button id="sm-reset" style="padding:5px 12px;border:1px solid #ccc;border-radius:5px;background:#fff;color:#666;cursor:pointer;font-size:0.85em;">Reset to samplesheet</button>` +
      `<span style="flex:1"></span>` +
      `<button id="sm-cancel" style="padding:5px 12px;border:1px solid #ccc;border-radius:5px;background:#fff;color:#666;cursor:pointer;font-size:0.85em;">Cancel</button>` +
      `<button id="sm-apply" style="padding:5px 14px;border:1px solid #1565c0;border-radius:5px;background:#1565c0;color:#fff;cursor:pointer;font-size:0.85em;font-weight:600;">Apply &amp; merge</button>` +
      `</div></div>`;

    document.body.appendChild(back);

    const chip = (s) => {
      const d = _molDotStyle(s);
      return (
        `<span class="sm-chip" draggable="true" data-sample="${_esc(s)}" ` +
        `style="display:inline-flex;align-items:center;gap:5px;margin:3px;padding:3px 8px;border:1px solid #b8c8e0;` +
        `border-radius:12px;background:#fff;font-size:0.78em;cursor:grab;">` +
        `<span title="Mol type: ${d.title}" style="display:inline-block;width:9px;height:9px;border-radius:50%;` +
        `flex-shrink:0;background:${d.bg};border:1px solid ${d.mt ? d.bg : "#b8c8e0"};"></span>` +
        `${_esc(s)}</span>`
      );
    };

    function render() {
      const query = String((back.querySelector("#sm-sample-search") || {}).value || "")
        .trim()
        .toLowerCase();
      // Pool = samples with no bin assignment.
      const allSamples = _samples();
      const pool = allSamples.filter((s) => !assign[s] && (!query || s.toLowerCase().includes(query)));
      const poolEl = back.querySelector("#sm-pool");
      poolEl.innerHTML = pool.length
        ? pool.map(chip).join("")
        : `<div style="color:#aab;font-size:0.78em;padding:8px;">${
            query ? "No matching unassigned samples." : "All samples assigned."
          }</div>`;

      // Bins = distinct specimen names in assign.
      const bins = new Map();
      Object.keys(assign).forEach((s) => {
        if (!bins.has(assign[s])) bins.set(assign[s], []);
        bins.get(assign[s]).push(s);
      });

      // Keep a stable, user-orderable sequence of boxes: preserve whatever
      // order the user last dragged into, append newly-created groups
      // (alphabetically) at the end, and drop groups that no longer exist.
      const liveKeys = Array.from(bins.keys());
      _binOrder = _binOrder.filter((g) => liveKeys.includes(g));
      liveKeys
        .filter((g) => !_binOrder.includes(g))
        .sort()
        .forEach((g) => _binOrder.push(g));

      const binsEl = back.querySelector("#sm-bins");
      if (!bins.size) {
        binsEl.innerHTML = `<div style="color:#aab;font-size:0.78em;padding:8px;">No specimens yet — drag samples onto “+ New specimen”.</div>`;
      } else {
        const renderedBins = _binOrder
          .map((g) => {
            const members = bins.get(g).sort();
            const shownMembers =
              !query || g.toLowerCase().includes(query)
                ? members
                : members.filter((s) => s.toLowerCase().includes(query));
            if (!shownMembers.length) return "";
            const multi = members.length > 1;
            // The whole box is the draggable/droppable unit: a background
            // envelope around the group name + all its member sample chips,
            // rather than a standalone name chip that can overflow. Drag the
            // grip handle to reorder boxes; drop a sample chip anywhere in
            // the box to assign it to this specimen.
            // Give the box a color the same way an individual sample row does —
            // a native <input type="color"> bound to sampleColors[g], so the
            // merged specimen can be recolored the same way as any sample.
            if (typeof sampleColors !== "undefined" && !sampleColors[g]) {
              sampleColors[g] =
                (sampleColors[members[0]] && sampleColors[members[0]]) ||
                (typeof PALETTE !== "undefined" && PALETTE.length
                  ? PALETTE[_binOrder.indexOf(g) % PALETTE.length]
                  : "#607d8b");
            }
            const binColor = (typeof sampleColors !== "undefined" && sampleColors[g]) || "#607d8b";
            return (
              `<div class="sm-drop sm-bin" draggable="true" data-bin="${_esc(g)}" role="button" tabindex="0" ` +
              `title="Select this specimen to review its metadata warnings" ` +
              `style="border:1px solid ${
                multi ? "#1565c0" : "#d6e2f5"
              };border-radius:6px;padding:6px 8px;margin-bottom:6px;background:${
                multi ? "#f0f7ff" : "#fff"
              };cursor:grab;">` +
              `<div style="display:flex;align-items:center;gap:6px;margin-bottom:2px;">` +
              `<i class="fas fa-grip-vertical" style="color:#b0bec5;font-size:0.78em;flex-shrink:0;" title="Drag box to reorder"></i>` +
              `<input type="color" class="sm-bin-color" data-bin="${_esc(g)}" value="${_esc(binColor)}" ` +
              `title="Specimen color" style="width:16px;height:16px;padding:0;border:1px solid #b8c8e0;border-radius:3px;flex-shrink:0;cursor:pointer;">` +
              `<b style="font-size:0.8em;color:#0d47a1;max-width:150px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;" title="${_esc(
                g,
              )}">${_esc(g)}</b>` +
              `<span style="font-size:0.7em;color:#88a;flex-shrink:0;">${members.length} sample(s)</span>` +
              `<span style="flex:1"></span>` +
              `<span class="sm-rename" data-bin="${_esc(
                g,
              )}" title="Rename specimen" style="cursor:pointer;color:#1565c0;font-size:0.78em;flex-shrink:0;"><i class="fas fa-pen"></i></span>` +
              `</div>` +
              `<div>${shownMembers.map(chip).join("")}</div>` +
              `</div>`
            );
          })
          .join("");
        binsEl.innerHTML =
          renderedBins ||
          `<div style="color:#aab;font-size:0.78em;padding:8px;">No matching samples or specimen groups.</div>`;
      }
      const matched = query
        ? allSamples.filter(
            (s) =>
              s.toLowerCase().includes(query) ||
              String(assign[s] || "")
                .toLowerCase()
                .includes(query),
          ).length
        : allSamples.length;
      const searchCount = back.querySelector("#sm-search-count");
      if (searchCount)
        searchCount.textContent = query ? `${matched} of ${allSamples.length}` : `${allSamples.length} samples`;
      _wireDnd();
      // Live conflict preview.
      _renderConflicts(findConflicts(_multiOnly(assign)));
    }

    // Assignment restricted to specimens that actually have >1 member.
    function _multiOnly(a) {
      const counts = {};
      Object.values(a).forEach((g) => (counts[g] = (counts[g] || 0) + 1));
      const out = {};
      Object.keys(a).forEach((s) => {
        if (counts[a[s]] > 1) out[s] = a[s];
      });
      return out;
    }

    let _resolved = {}; // spec → { field: value } chosen in the conflict panel

    function _paintSelectedBin() {
      back.querySelectorAll(".sm-bin").forEach((bin) => {
        const selected = bin.getAttribute("data-bin") === _selectedSpec;
        bin.style.boxShadow = selected ? "0 0 0 2px rgba(21,101,192,0.28)" : "";
        bin.style.borderColor = selected ? "#0d47a1" : "";
        bin.setAttribute("aria-pressed", selected ? "true" : "false");
      });
    }

    function _renderConflicts(conflicts) {
      const el = back.querySelector("#sm-conflicts");
      const multiSpecs = Array.from(new Set(Object.values(_multiOnly(assign))));
      if (!conflicts.length) {
        el.innerHTML = `<div style="border:1px solid #b7dfc5;background:#f3fbf6;border-radius:6px;padding:7px 10px;margin:4px 0 2px;font-size:0.78em;color:#2b6f44;"><i class="fas fa-circle-check"></i> No differing metadata${
          _selectedSpec ? ` for <b>${_esc(_selectedSpec)}</b>` : " in the current specimen groups"
        }.</div>`;
        _resolved = {};
        if (_selectedSpec && !multiSpecs.includes(_selectedSpec)) _selectedSpec = null;
        _paintSelectedBin();
        return;
      }
      const conflictSpecs = Array.from(new Set(conflicts.map((c) => c.spec)));
      if (!_selectedSpec || !multiSpecs.includes(_selectedSpec)) _selectedSpec = conflictSpecs[0] || multiSpecs[0];

      // Reconcile defaults without discarding choices already made for other
      // groups. This lets users click through every specimen and retain each
      // dropdown selection until Apply is pressed.
      Object.keys(_resolved).forEach((spec) => {
        if (!conflictSpecs.includes(spec)) delete _resolved[spec];
      });
      conflictSpecs.forEach((spec) => {
        const activeFields = new Set(conflicts.filter((c) => c.spec === spec).map((c) => c.field));
        Object.keys(_resolved[spec] || {}).forEach((field) => {
          if (!activeFields.has(field)) delete _resolved[spec][field];
        });
      });
      conflicts.forEach((c) => {
        _resolved[c.spec] = _resolved[c.spec] || {};
        const valid = c.kind === "multi" ? c.merged : c.options.map((o) => o.value);
        const old = _resolved[c.spec][c.field];
        if (old == null || (c.kind !== "multi" && !valid.includes(old))) {
          _resolved[c.spec][c.field] = c.kind === "multi" ? c.merged.slice() : c.options[0].value;
        } else if (c.kind === "multi") {
          _resolved[c.spec][c.field] = c.merged.slice();
        }
      });

      const shown = conflicts.filter((c) => c.spec === _selectedSpec);
      if (!shown.length) {
        el.innerHTML =
          `<div style="border:1px solid #b7dfc5;background:#f3fbf6;border-radius:6px;padding:7px 10px;margin:4px 0 2px;font-size:0.78em;color:#2b6f44;">` +
          `<i class="fas fa-circle-check"></i> No differing metadata for <b>${_esc(_selectedSpec)}</b>.</div>`;
        _paintSelectedBin();
        return;
      }
      const singles = shown.filter((c) => c.kind !== "multi");
      const multis = shown.filter((c) => c.kind === "multi");
      el.innerHTML =
        `<div style="border:1px solid #f0c36d;background:#fff9e8;border-radius:6px;padding:8px 10px;margin:4px 0 2px;">` +
        `<div style="font-size:0.8em;color:#8a6d1a;font-weight:600;margin-bottom:4px;">` +
        `<i class="fas fa-triangle-exclamation"></i> Metadata warnings for <b>${_esc(_selectedSpec)}</b> ` +
        `<span style="font-weight:400;color:#9b843f">(${conflictSpecs.indexOf(_selectedSpec) + 1} of ${
          conflictSpecs.length
        } group${conflictSpecs.length === 1 ? "" : "s"} with differences)</span></div>` +
        (singles.length
          ? `<div style="font-size:0.76em;color:#735c16;margin-bottom:6px;">Choose the one value the merged specimen should keep for each single-value field.</div>`
          : "") +
        singles
          .map((c, i) => {
            const opts = c.options
              .map(
                (o, j) =>
                  `<option value="${_esc(o.value)}" ${
                    _resolved[c.spec] && _resolved[c.spec][c.field] === o.value ? "selected" : ""
                  }>${_esc(o.value)} (${_esc(o.samples.join(", "))})</option>`,
              )
              .join("");
            return (
              `<label style="display:grid;grid-template-columns:minmax(130px,0.65fr) minmax(180px,1.35fr);gap:8px;align-items:center;margin:4px 0;font-size:0.8em;">` +
              `<span><b style="color:#8a6d1a;">${_esc(c.spec)}</b> · ${_esc(c.label)}</span>` +
              `<select class="sm-cf-select" data-spec="${_esc(c.spec)}" data-field="${_esc(
                c.field,
              )}" style="min-width:0;width:100%;padding:3px 5px;border:1px solid #d3b25c;border-radius:4px;background:#fff;">${opts}</select></label>`
            );
          })
          .join("") +
        (multis.length
          ? `<div style="border-top:${singles.length ? "1px solid #ead9a4" : "0"};margin-top:7px;padding-top:6px;">` +
            `<div style="font-size:0.78em;color:#2b6f44;font-weight:700;margin-bottom:3px;"><i class="fas fa-code-merge"></i> Multiselect metadata will be merged</div>` +
            multis
              .map(
                (c) =>
                  `<div style="font-size:0.76em;color:#47634f;margin:3px 0;"><b>${_esc(c.spec)} · ${_esc(
                    c.label,
                  )}:</b> ${c.merged.map(_esc).join(", ")}</div>`,
              )
              .join("") +
            `</div>`
          : "") +
        `</div>`;
      _paintSelectedBin();
      el.querySelectorAll(".sm-cf-select").forEach((select) => {
        select.addEventListener("change", () => {
          const spec = select.getAttribute("data-spec");
          const field = select.getAttribute("data-field");
          _resolved[spec] = _resolved[spec] || {};
          _resolved[spec][field] = select.value;
        });
      });
    }

    function _wireDnd() {
      back.querySelectorAll(".sm-chip").forEach((el) => {
        el.addEventListener("dragstart", (e) => {
          e.stopPropagation(); // don't let the enclosing box's dragstart also fire
          e.dataTransfer.effectAllowed = "move";
          e.dataTransfer.setData("text/plain", el.getAttribute("data-sample"));
          el.style.opacity = "0.5";
        });
        el.addEventListener("dragend", (e) => {
          e.stopPropagation();
          el.style.opacity = "";
        });
      });
      // Whole specimen boxes: dragging the box itself (grip handle or
      // background, not a sample chip inside it) reorders the boxes.
      back.querySelectorAll(".sm-bin").forEach((el) => {
        el.addEventListener("dragstart", (e) => {
          e.dataTransfer.effectAllowed = "move";
          e.dataTransfer.setData("application/x-sm-bin", el.getAttribute("data-bin"));
          el.style.opacity = "0.5";
        });
        el.addEventListener("dragend", () => (el.style.opacity = ""));
        const selectGroup = (e) => {
          if (e && e.target && e.target.closest(".sm-rename, .sm-bin-color")) return;
          _selectedSpec = el.getAttribute("data-bin");
          _renderConflicts(findConflicts(_multiOnly(assign)));
        };
        el.addEventListener("click", selectGroup);
        el.addEventListener("keydown", (e) => {
          if (e.key === "Enter" || e.key === " ") {
            e.preventDefault();
            selectGroup(e);
          }
        });
      });
      back.querySelectorAll(".sm-drop").forEach((el) => {
        el.addEventListener("dragover", (e) => {
          e.preventDefault();
          e.dataTransfer.dropEffect = "move";
          el.style.outline = "2px solid #1565c0";
        });
        el.addEventListener("dragleave", () => (el.style.outline = ""));
        el.addEventListener("drop", (e) => {
          e.preventDefault();
          el.style.outline = "";
          const targetBin = el.getAttribute("data-bin");

          // Reordering a specimen box (dropped onto another box).
          const draggedBin = e.dataTransfer.getData("application/x-sm-bin");
          if (draggedBin) {
            if (!targetBin || targetBin === "__new__" || targetBin === draggedBin) return;
            const from = _binOrder.indexOf(draggedBin);
            const to = _binOrder.indexOf(targetBin);
            if (from === -1 || to === -1) return;
            _binOrder.splice(from, 1);
            _binOrder.splice(to, 0, draggedBin);
            render();
            return;
          }

          // Assigning a sample chip to a bin (or back to the pool).
          const sample = e.dataTransfer.getData("text/plain");
          if (!sample) return;
          let bin = targetBin;
          if (bin === "__new__") {
            const name = (window.prompt("Name for the new specimen:", _suggestName(sample)) || "").trim();
            if (!name) return;
            bin = name;
          }
          if (bin === "") delete assign[sample]; // back to the pool
          else assign[sample] = bin;
          render();
        });
      });
      back.querySelectorAll(".sm-rename").forEach((el) => {
        el.addEventListener("click", (e) => {
          e.stopPropagation();
          const old = el.getAttribute("data-bin");
          const name = (window.prompt("Rename specimen:", old) || "").trim();
          if (!name || name === old) return;
          Object.keys(assign).forEach((s) => {
            if (assign[s] === old) assign[s] = name;
          });
          if (_selectedSpec === old) _selectedSpec = name;
          render();
        });
      });
      // Specimen color swatch — same mechanism as the individual sample color
      // picker in the right-panel sidebar (an <input type="color"> bound to
      // sampleColors[id]), just keyed by the specimen name instead of a raw
      // sample id. Takes effect immediately across the report, same as
      // recoloring a sample.
      back.querySelectorAll(".sm-bin-color").forEach((el) => {
        // Prevent the box's own draggable/click handlers from hijacking the
        // native color-picker interaction.
        el.addEventListener("mousedown", (e) => e.stopPropagation());
        el.addEventListener("click", (e) => e.stopPropagation());
        el.addEventListener("dragstart", (e) => e.preventDefault());
        el.addEventListener("input", (e) => {
          e.stopPropagation();
          const g = el.getAttribute("data-bin");
          if (typeof sampleColors !== "undefined") sampleColors[g] = e.target.value;
          if (typeof _refreshMapMarkerColors === "function") _refreshMapMarkerColors();
          if (typeof redraw === "function") redraw();
        });
      });
    }

    // Suggest a specimen name from a sample id by trimming a trailing
    // library/mol-type suffix (…-DNA, _RNA, .r1 …) — purely a convenience.
    function _suggestName(sample) {
      return String(sample).replace(/([._-])(dna|rna|r[12]|lib\d*|s\d+)$/i, "") || sample;
    }

    back.addEventListener("click", (e) => {
      if (e.target === back) _closeModal();
    });
    back.querySelector("#sm-close").addEventListener("click", _closeModal);
    back.querySelector("#sm-cancel").addEventListener("click", _closeModal);
    back.querySelector("#sm-sample-search").addEventListener("input", render);

    back.querySelector("#sm-reset").addEventListener("click", () => {
      if (typeof SPECIMEN_OVERRIDE !== "undefined") {
        Object.keys(SPECIMEN_OVERRIDE).forEach((k) => delete SPECIMEN_OVERRIDE[k]);
      }
      if (typeof SPECIMEN_META_RESOLVED !== "undefined") {
        Object.keys(SPECIMEN_META_RESOLVED).forEach((k) => delete SPECIMEN_META_RESOLVED[k]);
      }
      _closeModal();
      setEnabled(specimenMergeEnabled); // refresh views with defaults restored
    });

    back.querySelector("#sm-apply").addEventListener("click", () => {
      // Commit staged assignment → SPECIMEN_OVERRIDE, keeping only genuine
      // overrides (anything that differs from the samplesheet default).
      if (typeof SPECIMEN_OVERRIDE !== "undefined") {
        _samples().forEach((s) => {
          const meta = _metaFor(s);
          const metaSp =
            meta.specimen != null ? meta.specimen : meta.specimen_id != null ? meta.specimen_id : meta.specimen_group;
          const dflt = _norm(metaSp) || s;
          const val = assign[s] || s; // unassigned ⇒ its own specimen
          if (!val || val === dflt) delete SPECIMEN_OVERRIDE[s];
          else SPECIMEN_OVERRIDE[s] = val;
        });
      }
      // Persist chosen metadata resolutions for merged specimens.
      if (typeof SPECIMEN_META_RESOLVED !== "undefined") {
        Object.keys(SPECIMEN_META_RESOLVED).forEach((k) => delete SPECIMEN_META_RESOLVED[k]);
        Object.keys(_resolved).forEach((spec) => {
          SPECIMEN_META_RESOLVED[spec] = Object.assign({}, _resolved[spec]);
        });
      }
      _closeModal();
      setEnabled(true); // show the effect immediately
      // If the grouping produced no multi-sample specimen, nothing visibly
      // changes — say so rather than leave the user wondering.
      const grouped = typeof hasSpecimenGrouping === "function" && hasSpecimenGrouping();
      if (!grouped) {
        window.alert("No specimens were merged — drag two or more samples into the same specimen to combine them.");
      }
    });

    render();
  }

  /* ── one-click: combine every sample into a single specimen ─────────── */

  // Assign all known samples to one specimen group and turn merge on — the
  // fast path for "just treat this whole run as one specimen" without dragging
  // chips in the modal. Reversible: clears cleanly via Group… → Reset.
  function combineAll() {
    const samples = _samples();
    if (samples.length < 2) return;
    const dflt = "All samples";
    const name = (window.prompt("Combine ALL samples into one specimen named:", dflt) || "").trim();
    if (!name) return;
    const assign = {};
    samples.forEach((s) => (assign[s] = name));
    // Review before committing: single-valued metadata conflicts require an
    // explicit dropdown choice; multiselect fields show the union that will be
    // retained. Cancel leaves the current grouping untouched.
    openModal({
      assign,
      title: "Review metadata before combining all samples",
      message:
        `<b>Warning:</b> these samples may have different metadata. Review the choices below before merging. ` +
        `Single-value fields (for example latitude, longitude, country and state/territory) keep one selected value; ` +
        `multiselect fields (for example host disease) are union-merged.`,
    });
  }

  /* ── wiring ────────────────────────────────────────────────────────── */

  // Rich hover explanation of exactly what changes when the toggle flips,
  // mirroring the aggregation rules applied in _collapseSpecimens().
  function _mergeToggleTip() {
    const grouped = typeof hasSpecimenGrouping === "function" && hasSpecimenGrouping();
    const on = typeof specimenMergeEnabled !== "undefined" && specimenMergeEnabled;
    if (!grouped) {
      return (
        `<b style="color:#1565c0">Specimen merge</b><br>` +
        `<span style="color:#ccc">No specimens grouped yet — click to open <b>Group…</b> and drag samples together first.</span>`
      );
    }
    const action = on
      ? "Toggling <b>off</b> restores per-sample rows — every sample is shown and scored individually again."
      : "Toggling <b>on</b> collapses each specimen's samples (e.g. a DNA + RNA library from one swab) into one row.";
    return (
      `<b style="color:#1565c0">Specimen merge — ${on ? "On" : "Off"}</b><br>` +
      `<span style="color:#ccc">${action}</span>` +
      `<table style="border-collapse:collapse;font-size:0.85em;margin-top:5px">` +
      `<tr><td style="padding:1px 8px 1px 0;color:#90caf9"># Reads Aligned / Unique Reads / # Reads</td><td style="color:#fff">summed across samples</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#90caf9">TASS, Coverage, Breadth %, evidence scores</td><td style="color:#fff">max across samples</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#90caf9">Mean Depth</td><td style="color:#fff">summed across samples</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#90caf9">Mean MapQ / Mean BaseQ</td><td style="color:#fff">read-weighted mean</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#90caf9">High Consequence / Passes Threshold</td><td style="color:#fff">true if any sample is true</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#90caf9">Mol Type</td><td style="color:#fff">"both" if DNA + RNA mixed</td></tr>` +
      `<tr><td style="padding:1px 8px 1px 0;color:#90caf9">Everything else</td><td style="color:#fff">taken from the highest-TASS sample</td></tr>` +
      `</table>` +
      `<span style="color:#aaa;font-size:0.85em">Nothing in the underlying data is mutated — toggle off any time to see individual samples.</span>`
    );
  }

  function wire() {
    const tog = document.getElementById("specimen-merge-toggle");
    if (tog && !tog._wired) {
      tog._wired = true;
      tog.removeAttribute("title"); // superseded by the richer showTip hover below
      tog.addEventListener("mouseover", (ev) => {
        if (typeof showTip === "function") showTip(_mergeToggleTip(), ev);
      });
      tog.addEventListener("mousemove", (ev) => {
        if (typeof moveTip === "function") moveTip(ev);
      });
      tog.addEventListener("mouseout", () => {
        if (typeof hideTip === "function") hideTip();
      });
      tog.addEventListener("click", () => {
        if (typeof hideTip === "function") hideTip();
        const on = typeof specimenMergeEnabled !== "undefined" && specimenMergeEnabled;
        // Turning merge ON with nothing grouped yet would be a no-op — send the
        // user straight to the grouping modal instead of silently doing nothing.
        if (!on && !(typeof hasSpecimenGrouping === "function" && hasSpecimenGrouping())) {
          openModal();
          return;
        }
        setEnabled(!on);
      });
    }
    const grp = document.getElementById("specimen-merge-group-btn");
    if (grp && !grp._wired) {
      grp._wired = true;
      grp.addEventListener("click", openModal);
    }
    const combineBtn = document.getElementById("specimen-merge-combine-all-btn");
    if (combineBtn && !combineBtn._wired) {
      combineBtn._wired = true;
      combineBtn.addEventListener("click", combineAll);
    }

    // Default the merge ON when the samplesheet already groups samples into
    // multi-sample specimens (requirement: merge by the specimen column if
    // present) — but only once, and never fight a later user toggle.
    if (!wire._defaulted) {
      wire._defaulted = true;
      if (
        typeof hasSpecimenGrouping === "function" &&
        hasSpecimenGrouping() &&
        typeof specimenMergeEnabled !== "undefined" &&
        !specimenMergeEnabled
      ) {
        specimenMergeEnabled = true;
        if (typeof _invalidateFilterCache === "function") _invalidateFilterCache();
        if (typeof _rebuildMapMarkers === "function") _rebuildMapMarkers();
      }
    }

    // Keep the bar (toggle label + counts) in sync whenever the sample list is
    // rebuilt (uploads, deletions) by wrapping buildSampleList once.
    if (typeof buildSampleList === "function" && !buildSampleList._smWrapped) {
      const orig = buildSampleList;
      window.buildSampleList = function () {
        const r = orig.apply(this, arguments);
        refreshBar();
        return r;
      };
      window.buildSampleList._smWrapped = true;
    }

    refreshBar();
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", wire);
  else wire();

  // Expose a few helpers for other modules / debugging.
  window.specimenMerge = { refreshBar, setEnabled, openModal, combineAll, findConflicts };
})();
