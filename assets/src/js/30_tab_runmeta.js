/* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB: RUN METADATA       (data-tab="runmeta"  —  hidden if no meta)
       -     _buildRunMetaTable()  — renders the metadata grid (one row per
       -     sample, columns from the union of metadata keys present in DATA).
       -     _metaKeyLabel()       — pretty-prints raw keys for column headers.
       -     _switchMetaSub(id)    — activates a run-metadata analysis sub-tab.
       -     _updateMetaSubTabStates() — enables/disables sub-tabs based on data.
═══════════════════════════════════════════════════════════════════════════ */

/* Activate a Trends sub-tab pane and trigger its chart build.

   "geo" is no longer one of the sub-tabs — Mapping is its own top-level tab
   now — but it is still addressed through this function (by the PDF renderer,
   the "plot these coordinates" prompts, and the filter redraw) so callers
   don't have to know where a view lives. Asking for "geo" therefore builds the
   map and, when `opts.activate` is set, brings the Mapping tab forward.

   Only the panes inside #pane-trends take part in the active/inactive toggle;
   the geo pane is permanently active inside its own tab.                     */
function _switchMetaSub(id, opts) {
  opts = opts || {};
  if (id === "geo") {
    if (opts.activate) {
      const mb = document.getElementById("map-tab-btn");
      if (mb) mb.click();
    }
    if (typeof _buildGeoComparison === "function") _buildGeoComparison();
    return;
  }
  _activeMetaSub = id;
  document
    .querySelectorAll(".meta-subtab")
    .forEach((b) => b.classList.toggle("active", b.getAttribute("data-metasub") === id));
  document
    .querySelectorAll("#pane-trends .meta-subpane")
    .forEach((p) => p.classList.toggle("active", p.id === `meta-subpane-${id}`));
  // Trigger the relevant build now that the pane is visible
  if (id === "longi" && typeof _buildLongitudinalSection === "function") _buildLongitudinalSection();
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
  // Keep in step with _SITE_META_FIELDS / _HOST_META_FIELDS in
  // 33_categorical_metadata.js — inlined here for the same TDZ reason.
  const _hostFields = [
    "host_scientific_name",
    "host_disease",
    "environmental_site",
    "site",
    "site_name",
    "sampling_site",
    "collection_site",
    "site_type",
    "site_id",
  ];

  // ── Mapping is ALWAYS available (its own top-level tab now) so users can
  //    add latitude / longitude or country / state in the Metadata tab even
  //    when none was supplied in the samplesheet. The pane shows its own
  //    per-level "no data yet" hint when nothing is populated.

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
    "This tab needs at least one categorical metadata column to group on. Add a column in the " +
    '<a href="#" class="mg-goto-meta-link">Metadata tab</a> (e.g. <b>sample_type</b>, <b>site</b>, ' +
    "<b>environmental_site</b>) or upload a metadata CSV, then pick it " +
    'in the <b>"Group by"</b> bar.';

  // Keep the bar's option list in step with whatever metadata is loaded now.
  try {
    if (typeof _mgSyncGroupBar === "function") _mgSyncGroupBar();
    if (typeof _mgWireGroupBar === "function") _mgWireGroupBar();
  } catch (e) {}
  // Every mounted copy of the bar (Metadata / Mapping / Trends) hides together.
  document.querySelectorAll(".mg-bar-mount").forEach((el) => {
    el.style.display = hasGroups ? "" : "none";
  });

  // Mapping lives in its own tab now, so its "no data yet" banner is updated
  // here rather than through the sub-tab config list below.
  const _geoWarn = document.getElementById("meta-subtab-warn-geo");
  if (_geoWarn) _geoWarn.style.display = "none";

  const configs = [
    {
      id: "longi",
      ok: hasLongi,
      warn: "This tab requires a <b>collection_time</b> column in your metadata samplesheet (ISO date or M/D/YYYY).",
    },
    {
      id: "host",
      ok: hasHost,
      warn: "This tab requires at least one of: <b>host_scientific_name</b>, <b>host_disease</b>, or a site column (<b>environmental_site</b> / <b>site</b>) in your metadata samplesheet.",
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

  // If no active sub-tab yet, or the active one just became disabled → switch
  // to the first enabled Trends sub-tab.
  const _defaultSub = firstEnabled;
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
    // Generic site spellings — same concept as environmental_site, labelled
    // distinctly so a run carrying both doesn't show two identical headers.
    site: "Site",
    site_name: "Site Name",
    site_id: "Site ID",
    site_type: "Site Type",
    sampling_site: "Sampling Site",
    collection_site: "Collection Site",
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
    map_region: "Map Region",
  };
  return KNOWN[k] || k.replace(/_/g, " ").replace(/\b\w/g, (c) => c.toUpperCase());
}
