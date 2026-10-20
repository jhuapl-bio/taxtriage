/* ═══════════════════════════════════════════════════════════════════════════
       -  §  RESIZABLE RIGHT SIDEBAR
       -     Drag the splitter between #content and #sidebar to widen /
       -     narrow the right panel. Width is clamped (200 px ↔ 70 vw) and
       -     persisted across reloads via localStorage key
       -     "taxtriage:sidebarWidth".
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  const SIDEBAR_W_KEY = "taxtriage:sidebarWidth";
  const sidebar = document.getElementById("sidebar");
  const handle = document.getElementById("sidebar-resizer");
  if (!sidebar || !handle) return;

  // Restore last width if saved.
  try {
    const saved = parseInt(localStorage.getItem(SIDEBAR_W_KEY), 10);
    if (!isNaN(saved) && saved >= 200) {
      sidebar.style.width = Math.min(saved, Math.floor(window.innerWidth * 0.7)) + "px";
    }
  } catch (e) {
    /* localStorage may be blocked — silently ignore */
  }

  let _sbDrag = null;
  let _saveTimer = null;

  handle.addEventListener("mousedown", (e) => {
    _sbDrag = { startX: e.clientX, startW: sidebar.offsetWidth };
    handle.classList.add("dragging");
    document.body.classList.add("resizing-sidebar");
    e.preventDefault();
  });

  document.addEventListener("mousemove", (e) => {
    if (!_sbDrag) return;
    // Sidebar is on the RIGHT, so dragging left grows it.
    const dx = _sbDrag.startX - e.clientX;
    const minW = 200;
    const maxW = Math.floor(window.innerWidth * 0.7);
    const newW = Math.max(minW, Math.min(maxW, _sbDrag.startW + dx));
    sidebar.style.width = newW + "px";
  });

  document.addEventListener("mouseup", () => {
    if (!_sbDrag) return;
    _sbDrag = null;
    handle.classList.remove("dragging");
    document.body.classList.remove("resizing-sidebar");
    // Persist (debounced so a flurry of mouseups doesn't hammer storage).
    clearTimeout(_saveTimer);
    _saveTimer = setTimeout(() => {
      try {
        localStorage.setItem(SIDEBAR_W_KEY, String(sidebar.offsetWidth));
      } catch (e) {}
    }, 100);
    // Charts in the active tab read clientWidth — re-render so the
    // visualisations expand / contract to match the new content area.
    if (typeof redraw === "function") redraw();
  });

  // Re-clamp on window resize so the sidebar never exceeds 70vw.
  window.addEventListener("resize", () => {
    const maxW = Math.floor(window.innerWidth * 0.7);
    if (sidebar.offsetWidth > maxW) sidebar.style.width = maxW + "px";
  });
})();

/* ── Banner subtitle helper ─────────────────────────────────────────────────
         Builds the "N samples × M organisms × X% reads aligned" string from
         the current DATA and SAMPLE_META globals.  Called from init() and from
         the upload / clear handlers so all three stay in sync.
      ────────────────────────────────────────────────────────────────────────── */
function _buildBannerSub(opts) {
  const rawSamples = uniq([
    ...DATA.map((r) => r["Specimen ID"] || ""),
    ...Object.keys((typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {}),
  ]).filter(Boolean);
  const grouped =
    typeof specimenMergeEnabled !== "undefined" &&
    specimenMergeEnabled &&
    typeof hasSpecimenGrouping === "function" &&
    hasSpecimenGrouping();
  const entities = grouped ? uniq(rawSamples.map((s) => specimenOf(s))).filter(Boolean) : rawSamples;

  // Unique organisms: deduplicate on Taxonomic ID # (key), not subkey
  const uniqueOrgs = new Set(DATA.map((r) => r["Taxonomic ID #"] || "").filter(Boolean)).size;

  // % reads aligned (to an organism reference) = Σ(total_organism_reads from SAMPLE_META) / Σ(total_reads from SAMPLE_META)
  // Uses per-sample metadata from the JSON. Falls back to summing "# Reads Aligned"
  // from DATA rows if SAMPLE_META doesn't carry total_reads (e.g. plain TSV upload).
  const totalInputReads = rawSamples.reduce((s, sn) => s + (parseFloat((SAMPLE_META[sn] || {}).total_reads) || 0), 0);
  const totalOrgReads =
    rawSamples.reduce((s, sn) => s + (parseFloat((SAMPLE_META[sn] || {}).total_organism_reads) || 0), 0) ||
    DATA.reduce((s, r) => s + (parseFloat(r["# Reads Aligned"]) || 0), 0);
  let pctClass = totalInputReads > 0 ? ((totalOrgReads / totalInputReads) * 100).toFixed(1) + "%" : "N/A";
  if (pctClass === "0.0%" && totalOrgReads > 0) pctClass = "<0.1%";
  const totalReadsFmt = _fmtBig(totalInputReads).short;

  const nEnt = entities.length;
  const unit = grouped ? (nEnt === 1 ? "specimen" : "specimens") : nEnt === 1 ? "sample" : "samples";
  const stats = [
    `${nEnt} ${unit}`,
    `${uniqueOrgs} organism${uniqueOrgs === 1 ? "" : "s"}`,
    `${pctClass} aligned`,
    `${totalReadsFmt} reads`,
  ].join(" \u00B7 ");
  const mask = _maskBannerNote(rawSamples);

  // Banner: compact stats; the full wording (and masking detail) lives in the tooltip.
  const el = document.getElementById("banner-sub");
  if (el) {
    el.title =
      `${nEnt} ${grouped ? "specimen(s)" : "sample(s)"} \u2022 ${uniqueOrgs} unique organism(s) \u2022 ` +
      `${pctClass} reads aligned \u2022 ${totalReadsFmt} total reads` +
      (mask.text ? `\n\u26A0 Region masking: ${mask.text}` : "");
  }
  // opts.long (PDF cover): keep the masking note inline, since there is no pill there.
  return opts && opts.long && mask.text ? `${stats} \u00B7 \u26A0 Region masking: ${mask.text}` : stats;
}

/* ── Region-masking notice ──────────────────────────────────────────────────
   When match_paths ran with --mask_bed / --fcs, each sample's metadata carries
   a `masking` block (applied, mode, masked_bp, reads_masked, ...). Scores in
   this report were then computed only over unmasked / targeted bases, so say
   so: a short amber pill in the banner (#banner-mask, per-sample detail in its
   tooltip) and the full sentence on the PDF cover (_buildBannerSub({long})). */
function _maskBannerNote(samples) {
  const meta = (typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {};
  const masked = samples.filter((s) => ((meta[s] || {}).masking || {}).applied);
  const pill = document.getElementById("banner-mask");
  const label = document.getElementById("banner-mask-label");
  if (!masked.length) {
    if (pill) pill.hidden = true;
    return { text: "" };
  }
  const modes = uniq(masked.map((s) => meta[s].masking.mode || "exclude"));
  const fcs = masked.some((s) => {
    const m = meta[s].masking;
    return (m.bed || []).concat(m.exclude_bed || []).some((b) => /fcs|mask\.exclude/i.test(b));
  });
  const kind =
    modes.length > 1 ? "mixed" : modes[0] === "target" ? "target" : fcs ? "fcs" : "exclude";
  const what = {
    mixed: "mixed modes",
    target: "TASS from targeted regions only",
    fcs: "FCS contaminant regions excluded",
    exclude: "masked regions excluded",
  }[kind];
  const short = { mixed: "Masked", target: "Targeted", fcs: "FCS-masked", exclude: "Masked" }[kind];
  const count = `${masked.length}/${samples.length}`;
  const esc = typeof _esc === "function" ? _esc : (v) => String(v);
  const fmt = (v) => (Number(v) || 0).toLocaleString();
  const rows = masked
    .map((s) => {
      const m = meta[s].masking;
      const pct = m.reads_considered ? ((100 * (m.reads_masked || 0)) / m.reads_considered).toFixed(2) : "0";
      return (
        `<tr><td>${esc(s)}</td><td>${esc(m.mode || "exclude")}</td><td class="num">${fmt(m.references_masked)}</td>` +
        `<td class="num">${fmt(m.masked_bp)}</td><td class="num">${fmt(m.reads_masked)} (${pct}%)</td>` +
        `<td>${m.penalty_power ? esc(m.penalty_power) : "\u2014"}</td></tr>`
      );
    })
    .join("");
  const bams = masked.map((s) => meta[s].masking.masked_bam).filter(Boolean);
  const text = `${what} (${count} sample${samples.length === 1 ? "" : "s"})`;
  if (pill) {
    pill.hidden = false;
    pill.removeAttribute("title"); // native tooltips don't show in every viewer; use the popover
    if (label) label.textContent = `${short} ${count}`;
    let pop = document.getElementById("banner-mask-pop");
    if (!pop) {
      pop = document.createElement("div");
      pop.id = "banner-mask-pop";
      pop.className = "banner-mask-pop";
      pop.setAttribute("role", "tooltip");
      pill.appendChild(pop);
      pill.setAttribute("aria-describedby", "banner-mask-pop");
      // Click / tap pins it open (touch screens have no hover); Esc or a click elsewhere closes.
      pill.addEventListener("click", (ev) => {
        ev.stopPropagation();
        if (ev.target.closest && ev.target.closest(".banner-mask-pop")) return; // clicks inside keep it open
        pill.classList.toggle("open");
      });
      document.addEventListener("click", () => pill.classList.remove("open"));
      pill.addEventListener("keydown", (ev) => {
        if (ev.key === "Escape") pill.classList.remove("open");
        if (ev.key === "Enter" || ev.key === " ") {
          ev.preventDefault();
          pill.classList.toggle("open");
        }
      });
    }
    pop.innerHTML =
      `<div class="bmp-title">Region masking: ${esc(text)}</div>` +
      `<div class="bmp-note">Breadth, depth, Gini and TASS were computed only over unmasked bases.</div>` +
      `<table><thead><tr><th>Sample</th><th>Mode</th><th>Refs</th><th>bp masked</th><th>Reads removed</th><th>Penalty</th></tr></thead>` +
      `<tbody>${rows}</tbody></table>` +
      (bams.length ? `<div class="bmp-note">Masked BAMs: ${bams.map(esc).join(", ")}</div>` : "");
    pill.setAttribute("aria-label", `Region masking: ${text}`);
  }
  return { text };
}
