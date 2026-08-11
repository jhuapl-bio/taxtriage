/* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB: RUN METADATA       (data-tab="runmeta"  —  hidden if no meta)
       -     _buildRunMetaTable()  — renders the metadata grid (one row per
       -     sample, columns from the union of metadata keys present in DATA).
       -     _metaKeyLabel()       — pretty-prints raw keys for column headers.
       -     _switchMetaSub(id)    — activates a run-metadata analysis sub-tab.
       -     _updateMetaSubTabStates() — enables/disables sub-tabs based on data.
═══════════════════════════════════════════════════════════════════════════ */

// Activate a sub-tab pane and trigger its chart build.
function _switchMetaSub(id) {
  _activeMetaSub = id;
  document
    .querySelectorAll(".meta-subtab")
    .forEach((b) => b.classList.toggle("active", b.getAttribute("data-metasub") === id));
  document
    .querySelectorAll(".meta-subpane")
    .forEach((p) => p.classList.toggle("active", p.id === `meta-subpane-${id}`));
  // Trigger the relevant build now that the pane is visible
  if (id === "longi" && typeof _buildLongitudinalSection === "function") _buildLongitudinalSection();
  if (id === "geo" && typeof _buildGeoComparison === "function") _buildGeoComparison();
  if (id === "host" && typeof _buildHostDisease === "function") _buildHostDisease();
  if (id === "ghm" && typeof _mgBuildGroupHeatmap === "function") _mgBuildGroupHeatmap();
  if (id === "net" && typeof _mgBuildNetwork === "function") _mgBuildNetwork();
  if (id === "cmp" && typeof _buildComparison === "function") _buildComparison();
}

// Inspect RUN_META and enable/disable each sub-tab; auto-switch away from
// a disabled tab.  Also shows/hides the per-pane warning banner.
function _updateMetaSubTabStates() {
  // ── Longitudinal: needs at least one record with a parseable collection_time
  const hasLongi = (RUN_META || []).some((r) => {
    try {
      return !!_parseLongiDate(r.collection_time);
    } catch (e) {
      return false;
    }
  });

  // Inline field lists to avoid TDZ if this is called before the const declarations
  // lower in the script (_GEO_META_FIELDS / _HOST_META_FIELDS at ~line 19970).
  const _geoFields = ["sample_origin_country", "sample_origin_state_province_territory"];
  const _hostFields = ["host_scientific_name", "host_disease", "environmental_site"];

  // ── Mapping & Geography: ALWAYS available so users can add latitude /
  //    longitude or country / state directly in the table above, even when
  //    none was supplied in the samplesheet. The pane shows its own per-level
  //    "no data yet" hint when nothing is populated.
  const hasGeo = true;

  // ── Host & Disease: needs host / disease / site field
  const hasHost = typeof _anyMetaValue === "function" && _anyMetaValue(_hostFields);

  // ── Cross-Entry: always available (shows internal no-data message when < 2 entries)
  const hasCmp = true;

  // ── Group views: available as soon as there is at least one column worth
  //    grouping on. With none, both views would fall back to per-sample rows,
  //    which the heatmap / cross-sample tabs already cover better — so the
  //    warning points the user at the "Group by" bar instead.
  let _groupCands = 0;
  try {
    if (window.metaGrouping) _groupCands = window.metaGrouping.candidates().length;
  } catch (e) {}
  const hasGroups = _groupCands > 0;
  const _groupWarn =
    "This tab needs at least one categorical metadata column to group on. Add a column in the table " +
    "above (e.g. <b>sample_type</b>, <b>environmental_site</b>) or upload a metadata CSV, then pick it " +
    'in the <b>"Group by"</b> bar.';

  // Keep the bar's option list in step with whatever metadata is loaded now.
  try {
    if (typeof _mgSyncGroupBar === "function") _mgSyncGroupBar();
    if (typeof _mgWireGroupBar === "function") _mgWireGroupBar();
  } catch (e) {}
  const _mgBar = document.getElementById("mg-bar");
  if (_mgBar) _mgBar.style.display = hasGroups ? "" : "none";

  const configs = [
    {
      id: "longi",
      ok: hasLongi,
      warn: "This tab requires a <b>collection_time</b> column in your metadata samplesheet (ISO date or M/D/YYYY).",
    },
    {
      id: "geo",
      ok: hasGeo,
      warn: "This tab requires at least one of: <b>sample_origin_country</b> or <b>sample_origin_state_province_territory</b> in your metadata samplesheet.",
    },
    {
      id: "host",
      ok: hasHost,
      warn: "This tab requires at least one of: <b>host_scientific_name</b>, <b>host_disease</b>, or <b>environmental_site</b> in your metadata samplesheet.",
    },
    { id: "ghm", ok: hasGroups, warn: _groupWarn },
    { id: "net", ok: hasGroups, warn: _groupWarn },
    { id: "cmp", ok: hasCmp, warn: "" },
  ];

  let firstEnabled = null;
  configs.forEach(({ id, ok, warn }) => {
    const btn = document.querySelector(`.meta-subtab[data-metasub="${id}"]`);
    const warnEl = document.getElementById(`meta-subtab-warn-${id}`);
    const warnMsg = document.getElementById(`meta-subtab-warn-${id}-msg`);
    if (btn) btn.disabled = !ok;
    if (warnEl) warnEl.style.display = !ok && warn ? "flex" : "none";
    if (warnMsg) warnMsg.innerHTML = warn;
    if (ok && !firstEnabled) firstEnabled = id;
  });

  // If no active sub-tab yet, or the active one just became disabled → switch.
  // Default to Mapping & Geography (the map) when it's available, otherwise
  // fall back to the first enabled sub-tab.
  const _defaultSub = (configs.find((c) => c.id === "geo" && c.ok) || {}).id || firstEnabled;
  const activeBtn = _activeMetaSub && document.querySelector(`.meta-subtab[data-metasub="${_activeMetaSub}"]`);
  if (!_activeMetaSub || (activeBtn && activeBtn.disabled)) {
    if (_defaultSub) _switchMetaSub(_defaultSub);
  }
}

// Pretty-print a raw key name → display label
function _metaKeyLabel(k) {
  const KNOWN = {
    sample_name: "Sample",
    run_id: "Run ID",
    sample_id: "Sample ID",
    organism: "Organism",
    submitter_organization_name: "Submitter Organization",
    sample_origin_country: "Country",
    sample_origin_state_province_territory: "State / Province / Territory",
    host_scientific_name: "Host (Scientific Name)",
    environmental_site: "Environmental Site",
    host_disease: "Host Disease",
    library_preparation_kit: "Library Prep Kit",
    sequencing_instrument: "Sequencing Instrument",
    sequencing_protocol_primer_set: "Primer Set",
    sequencing_platform: "Sequencing Platform",
    latitude: "Latitude",
    longitude: "Longitude",
    depth: "Depth (m)",
    salinity: "Salinity (PSU)",
    collection_time: "Collection Time",
    location: "Location",
  };
  return KNOWN[k] || k.replace(/_/g, " ").replace(/\b\w/g, (c) => c.toUpperCase());
}
