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
         Builds the "N samples × M organisms × X% reads classified" string from
         the current DATA and SAMPLE_META globals.  Called from init() and from
         the upload / clear handlers so all three stay in sync.
      ────────────────────────────────────────────────────────────────────────── */
function _buildBannerSub() {
  const samples = uniq(DATA.map((r) => r["Specimen ID"] || "")).filter(Boolean);

  // Unique organisms: deduplicate on Taxonomic ID # (key), not subkey
  const uniqueOrgs = new Set(DATA.map((r) => r["Taxonomic ID #"] || "").filter(Boolean)).size;

  // % reads classified = Σ(total_organism_reads from SAMPLE_META) / Σ(total_reads from SAMPLE_META)
  // Uses per-sample metadata from the JSON. Falls back to summing "# Reads Aligned"
  // from DATA rows if SAMPLE_META doesn't carry total_reads (e.g. plain TSV upload).
  const totalInputReads = samples.reduce((s, sn) => s + (parseFloat((SAMPLE_META[sn] || {}).total_reads) || 0), 0);
  const totalOrgReads =
    samples.reduce((s, sn) => s + (parseFloat((SAMPLE_META[sn] || {}).total_organism_reads) || 0), 0) ||
    DATA.reduce((s, r) => s + (parseFloat(r["# Reads Aligned"]) || 0), 0);
  let pctClass = totalInputReads > 0 ? ((totalOrgReads / totalInputReads) * 100).toFixed(1) + "%" : "N/A";
  if (pctClass === "0.0%" && totalOrgReads > 0) pctClass = "<0.1%";

  return `${samples.length} sample(s) \u2022 ${uniqueOrgs} unique organism(s) \u2022 ${pctClass} reads classified`;
}
