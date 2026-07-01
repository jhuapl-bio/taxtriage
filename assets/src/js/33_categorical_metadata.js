/* ═══════════════════════════════════════════════════════════════════════════
       -  §  CATEGORICAL METADATA COMPARISONS  (geo + host/disease bar charts)
       -     Driven by the AMD-P/Talos standard metadata fields. Each chart groups
       -     the currently-visible detections by a categorical metadata field
       -     (country / state, or host / disease / site) and aggregates a
       -     chosen metric. Sections self-hide when no relevant metadata is present.
       -     Functions: _aggregateMetaByField, _renderMetaBarChart,
       -     _buildGeoComparison, _buildHostDisease.
═══════════════════════════════════════════════════════════════════════════ */
const _META_METRIC_LABELS = {
  samples: "Samples",
  detections: "Detections",
  reads: "# Reads Aligned",
  tass: "Mean TASS Score",
};
const _GEO_META_FIELDS = ["sample_origin_country", "sample_origin_state_province_territory"];
const _HOST_META_FIELDS = ["host_scientific_name", "host_disease", "environmental_site"];

// True if any RUN_META record has a non-empty value for one of `fields`.
function _anyMetaValue(fields) {
  return (RUN_META || []).some((r) =>
    fields.some((f) => {
      const _v = r[f];
      if (_v == null) return false;
      if (Array.isArray(_v)) return _v.some((s) => String(s).trim() !== "");
      return String(_v).trim() !== "";
    }),
  );
}

// Group visible detections by a categorical metadata field and aggregate a metric.
// Returns [{label, value}] sorted descending by value.
// Array field values (e.g. host_disease split on commas) expand into multiple groups.
function _aggregateMetaByField(fieldKey, metric, sampleSet) {
  // sample_name → array of categorical values (always an array for uniform handling)
  // sampleSet (optional): restrict to these sample_names (used when drilling
  // into a single country's states).
  const valBy = {};
  (RUN_META || []).forEach((r) => {
    if (sampleSet && !sampleSet.has(r.sample_name)) return;
    const v = r[fieldKey];
    if (v == null) return;
    const vals = Array.isArray(v) ? v.map((s) => String(s).trim()).filter(Boolean) : [String(v).trim()].filter(Boolean);
    if (vals.length) valBy[r.sample_name] = vals;
  });
  if (!Object.keys(valBy).length) return [];

  const groups = {};
  const ensure = (val) =>
    (groups[val] = groups[val] || { samples: new Set(), detections: 0, reads: 0, tassSum: 0, tassN: 0 });

  if (metric === "samples") {
    // Count visible (non-hidden) samples per category value
    Object.entries(valBy).forEach(([s, vals]) => {
      if (sampleHidden[s]) return;
      vals.forEach((val) => ensure(val).samples.add(s));
    });
  } else {
    filteredData().forEach((r) => {
      const s = r["Specimen ID"] || "";
      const vals = valBy[s];
      if (!vals) return;
      vals.forEach((val) => {
        const g = ensure(val);
        g.samples.add(s);
        g.detections += 1;
        const reads = parseFloat(r["# Reads Aligned"]);
        if (!isNaN(reads)) g.reads += reads;
        const tass = parseFloat(r["TASS Score"]);
        if (!isNaN(tass)) {
          g.tassSum += tass;
          g.tassN += 1;
        }
      });
    });
  }

  const out = Object.entries(groups).map(([label, g]) => {
    let value;
    if (metric === "samples") value = g.samples.size;
    else if (metric === "detections") value = g.detections;
    else if (metric === "reads") value = g.reads;
    else value = g.tassN ? g.tassSum / g.tassN : 0;
    return { label, value };
  });
  out.sort((a, b) => b.value - a.value);
  return out;
}

// Horizontal bar chart for [{label, value}] data into the element `wrapId`.
// opts (optional): { width, fitHeight } — when given the chart sizes itself
// to the supplied box (used by the resizable map overlay) instead of using
// the container's natural width and a fixed row height.
function _renderMetaBarChart(wrapId, data, metric, opts) {
  opts = opts || {};
  const wrap = document.getElementById(wrapId);
  if (!wrap) return;
  wrap.innerHTML = "";
  if (!data.length) return;

  const isFloat = metric === "tass";
  const fmt = (v) => (isFloat ? v.toFixed(1) : Math.round(v).toLocaleString());

  const top = data.slice(0, 25);
  const compact = !!(opts.width || opts.fitHeight);
  const marginT = 12,
    marginR = compact ? 46 : 70,
    marginB = 30;
  // Width: explicit (overlay) → container → fallback. No hard 420 floor so
  // the chart can shrink to fit a narrow, resizable panel.
  const W = Math.max(170, Math.floor(opts.width || wrap.clientWidth || 600));
  // Row height: when a fit height is supplied, distribute rows to fill it
  // (clamped so labels stay legible); otherwise use the original sizing.
  const avail = opts.fitHeight ? Math.max(60, opts.fitHeight - marginT - marginB) : 460;
  const rowH = Math.max(opts.fitHeight ? 15 : 22, Math.min(38, Math.floor(avail / Math.max(top.length, 1))));
  // Left margin tracks label length but never eats more than ~55% of width.
  const marginL = Math.min(
    Math.max(70, top.reduce((m, d) => Math.max(m, String(d.label).length * 6.0), 0) + 12),
    Math.floor(W * 0.55),
  );
  const chartH = top.length * rowH;
  const H = chartH + marginT + marginB;
  const iW = Math.max(20, W - marginL - marginR);
  const xMax = d3.max(top, (d) => d.value) || 1;

  const svg = d3.select(wrap).append("svg").attr("width", W).attr("height", H).style("overflow", "visible");
  const g = svg.append("g").attr("transform", `translate(${marginL},${marginT})`);

  const x = d3.scaleLinear().domain([0, xMax]).range([0, iW]).nice();
  const y = d3
    .scaleBand()
    .domain(top.map((d) => d.label))
    .range([0, chartH])
    .padding(0.2);

  g.append("g").attr("class", "axis").attr("transform", `translate(0,${chartH})`).call(d3.axisBottom(x).ticks(5));
  g.append("g").attr("class", "axis").call(d3.axisLeft(y).tickSize(0)).select(".domain").remove();

  const color = d3.scaleOrdinal(d3.schemeTableau10);
  const metricLabel = _META_METRIC_LABELS[metric] || metric;

  top.forEach((d, i) => {
    g.append("rect")
      .attr("x", 0)
      .attr("y", y(d.label))
      .attr("height", y.bandwidth())
      .attr("width", x(d.value))
      .attr("fill", color(i))
      .attr("rx", 2)
      .on("mouseover", (ev) => showTip(`<b>${d.label}</b><br>${metricLabel}: ${fmt(d.value)}`, ev))
      .on("mousemove", moveTip)
      .on("mouseout", hideTip);
    g.append("text")
      .attr("x", x(d.value) + 5)
      .attr("y", y(d.label) + y.bandwidth() / 2)
      .attr("dy", "0.35em")
      .style("font-size", "10px")
      .attr("fill", "#444")
      .text(fmt(d.value));
  });
}

/* ── Choropleth support ─────────────────────────────────────────────────
       -   Renders a d3-geo Natural-Earth choropleth of the visible detections
       -   shaded by the selected metric. Country level uses lightweight
       -   world-atlas TopoJSON; State / Province / Territory level lazily loads
       -   Natural Earth admin-1 GeoJSON. Boundaries are fetched from public CDNs
       -   at view time (consistent with the Map tab's online tiles) and the view
       -   degrades to the ranked bar list if the network is unavailable.
      ──────────────────────────────────────────────────────────────────────── */
const _GEO_CACHE = {};
function _geoFetchJSON(key, url) {
  if (_GEO_CACHE[key] !== undefined) return _GEO_CACHE[key];
  // Try each URL in turn so a single CDN hiccup doesn't disable the map.
  const tryNext = (i) => {
    if (i >= url.length) return Promise.reject(new Error("all geo sources failed"));
    return d3.json(url[i]).catch(() => tryNext(i + 1));
  };
  const p = tryNext(0).catch((e) => {
    _GEO_CACHE[key] = undefined;
    throw e;
  });
  _GEO_CACHE[key] = p;
  return p;
}
// Plain GeoJSON sources (depend only on d3; no topojson needed).
function _geoCountries() {
  return _geoFetchJSON("countries", [
    "https://cdn.jsdelivr.net/gh/nvkelso/natural-earth-vector@master/geojson/ne_110m_admin_0_countries.geojson",
    "https://cdn.jsdelivr.net/gh/johan/world.geo.json@master/countries.geo.json",
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_110m_admin_0_countries.geojson",
  ]).then((fc) => (fc && fc.features) || []);
}
function _geoAdmin1() {
  return _geoFetchJSON("admin1", [
    "https://cdn.jsdelivr.net/gh/nvkelso/natural-earth-vector@master/geojson/ne_50m_admin_1_states_provinces.geojson",
    "https://cdn.jsdelivr.net/gh/nvkelso/natural-earth-vector@v5.1.2/geojson/ne_50m_admin_1_states_provinces.geojson",
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_1_states_provinces.geojson",
  ]).then((fc) => (fc && fc.features) || []);
}
function _geoNorm(s) {
  return String(s == null ? "" : s)
    .normalize("NFD")
    .replace(/[\u0300-\u036f]/g, "")
    .toLowerCase()
    .replace(/[.,'’`]/g, "")
    .replace(/&/g, "and")
    .replace(/\s+/g, " ")
    .trim()
    .replace(/^the /, "");
}
const _COUNTRY_ALIASES = {
  "united state of america": "united states of america",
  "united states": "united states of america",
  "united states of america": "united states of america",
  usa: "united states of america",
  us: "united states of america",
  "u s a": "united states of america",
  "u s": "united states of america",
  uk: "united kingdom",
  "great britain": "united kingdom",
  england: "united kingdom",
  "czech republic": "czechia",
  "ivory coast": "cote divoire",
  burma: "myanmar",
  "south korea": "south korea",
  "republic of korea": "south korea",
  "north korea": "north korea",
  "russian federation": "russia",
  "democratic republic of the congo": "dem rep congo",
};
const _US_STATE_ABBR = {
  AL: "Alabama",
  AK: "Alaska",
  AZ: "Arizona",
  AR: "Arkansas",
  CA: "California",
  CO: "Colorado",
  CT: "Connecticut",
  DE: "Delaware",
  FL: "Florida",
  GA: "Georgia",
  HI: "Hawaii",
  ID: "Idaho",
  IL: "Illinois",
  IN: "Indiana",
  IA: "Iowa",
  KS: "Kansas",
  KY: "Kentucky",
  LA: "Louisiana",
  ME: "Maine",
  MD: "Maryland",
  MA: "Massachusetts",
  MI: "Michigan",
  MN: "Minnesota",
  MS: "Mississippi",
  MO: "Missouri",
  MT: "Montana",
  NE: "Nebraska",
  NV: "Nevada",
  NH: "New Hampshire",
  NJ: "New Jersey",
  NM: "New Mexico",
  NY: "New York",
  NC: "North Carolina",
  ND: "North Dakota",
  OH: "Ohio",
  OK: "Oklahoma",
  OR: "Oregon",
  PA: "Pennsylvania",
  RI: "Rhode Island",
  SC: "South Carolina",
  SD: "South Dakota",
  TN: "Tennessee",
  TX: "Texas",
  UT: "Utah",
  VT: "Vermont",
  VA: "Virginia",
  WA: "Washington",
  WV: "West Virginia",
  WI: "Wisconsin",
  WY: "Wyoming",
  DC: "District of Columbia",
  PR: "Puerto Rico",
  GU: "Guam",
  VI: "United States Virgin Islands",
  AS: "American Samoa",
  MP: "Northern Mariana Islands",
};
// Normalized lookup keys for a metadata label, given the geo level.
function _geoMatchKeys(field, label) {
  const out = [];
  const raw = String(label == null ? "" : label).trim();
  let n = _geoNorm(raw);
  if (field === "sample_origin_country") {
    n = _COUNTRY_ALIASES[n] || n;
    out.push(n);
  } else if (field === "sample_origin_state_province_territory") {
    out.push(n);
    const up = raw.toUpperCase();
    if (_US_STATE_ABBR[up]) out.push(_geoNorm(_US_STATE_ABBR[up]));
  } else {
    out.push(n);
  }
  return out.filter(Boolean);
}
// Candidate normalized names carried by a boundary feature.
function _geoFeatureNames(field, f) {
  const p = f.properties || {};
  const cands = [];
  if (field === "sample_origin_state_province_territory") {
    [p.name, p.name_en, p.gn_name, p.woe_name, p.name_alt, p.abbrev, p.postal].forEach(
      (v) => v && cands.push(_geoNorm(v)),
    );
  } else {
    [p.name, p.NAME, p.name_long, p.NAME_LONG, p.admin, p.ADMIN, p.geounit, p.sovereignt, p.SOVEREIGNT].forEach((v) => {
      if (!v) return;
      const n = _geoNorm(v);
      cands.push(n);
      if (_COUNTRY_ALIASES[n]) cands.push(_COUNTRY_ALIASES[n]); // normalize map-side names too
    });
  }
  return cands;
}
function _geoFeatureMatch(field, f, metaByNorm, stateCountries) {
  const cands = _geoFeatureNames(field, f);
  for (const c of cands) {
    if (!metaByNorm.has(c)) continue;
    // State / territory disambiguation: an abbreviation like "GA" matches
    // Georgia (US) AND Gorno-Altay (RU). When the metadata carries a country
    // for this state label, only match boundary features in that country.
    if (field === "sample_origin_state_province_territory" && stateCountries) {
      const allowed = stateCountries.get(c);
      if (allowed && allowed.size && !allowed.has(_geoFeatureCountryNorm(f))) continue;
    }
    return metaByNorm.get(c);
  }
  return null;
}
// normalized state key → Set of normalized country names (from the metadata),
// used to constrain admin-1 matches to the sample's own country.
function _geoStateCountryMap() {
  const m = new Map();
  (RUN_META || []).forEach((r) => {
    let cnorm = _geoNorm(r.sample_origin_country);
    cnorm = _COUNTRY_ALIASES[cnorm] || cnorm;
    if (!cnorm) return;
    const states = r.sample_origin_state_province_territory;
    const arr = states == null ? [] : Array.isArray(states) ? states : [states];
    arr.forEach((sv) => {
      _geoMatchKeys("sample_origin_state_province_territory", sv).forEach((k) => {
        if (!m.has(k)) m.set(k, new Set());
        m.get(k).add(cnorm);
      });
    });
  });
  return m;
}
// Normalized country name a boundary feature belongs to (for state→country
// filtering when drilling). Admin-1 features carry an `admin` country field.
function _geoFeatureCountryNorm(f) {
  const p = f.properties || {};
  const c = p.admin || p.ADMIN || p.sovereignt || p.SOVEREIGNT || p.geonunit || p.GEOUNIT || p.name || p.NAME || "";
  const n = _geoNorm(c);
  return _COUNTRY_ALIASES[n] || n;
}
// One-time dialog when a state/territory label matches admin-1 regions in
// more than one country (e.g. "GA" → Georgia US + Gorno-Altay RU). Explains
// whether the metadata country resolved it.
const _geoAmbigWarned = new Set();
function _geoWarnAmbiguousStates(flat, allFeatures, stateCountries) {
  if (!Array.isArray(flat) || !Array.isArray(allFeatures)) return;
  const msgs = [];
  flat.forEach((d) => {
    const keys = _geoMatchKeys("sample_origin_state_province_territory", d.label);
    const countriesHit = new Set();
    allFeatures.forEach((f) => {
      const names = _geoFeatureNames("sample_origin_state_province_territory", f);
      if (names.some((nm) => keys.includes(nm))) countriesHit.add(_geoFeatureCountryNorm(f));
    });
    if (countriesHit.size <= 1) return; // unambiguous
    const sig = _geoNorm(d.label) + "::" + [...countriesHit].sort().join("|");
    if (_geoAmbigWarned.has(sig)) return;
    _geoAmbigWarned.add(sig);
    const allowed = new Set();
    keys.forEach((k) => (stateCountries.get(k) || new Set()).forEach((c) => allowed.add(c)));
    const pretty = (s) => String(s).replace(/\b\w/g, (m2) => m2.toUpperCase());
    const hitList = [...countriesHit].map(pretty).join(", ");
    if (allowed.size)
      msgs.push(
        `“${d.label}” matches regions in multiple countries (${hitList}). ` +
          `Restricting to the sample's country: ${[...allowed].map(pretty).join(", ")}.`,
      );
    else
      msgs.push(
        `“${d.label}” matches regions in multiple countries (${hitList}) and no sample_origin_country is set, ` +
          `so all are shown. Set the country in the Run Metadata tab to disambiguate.`,
      );
  });
  if (msgs.length) setTimeout(() => alert("Ambiguous state / territory names:\n\n• " + msgs.join("\n• ")), 0);
}
// sample_names whose origin country matches a normalized country name.
function _geoCountrySampleSet(countryNorm) {
  const set = new Set();
  (RUN_META || []).forEach((r) => {
    let n = _geoNorm(r.sample_origin_country);
    n = _COUNTRY_ALIASES[n] || n;
    if (n && n === countryNorm) set.add(r.sample_name);
  });
  return set;
}
// ── Reverse-geocode: fill country / state / location from lat-long ───────
// Uses the same Natural-Earth GeoJSON the choropleth loads, doing a
// point-in-polygon test (d3.geoContains) per record. Only fills fields that
// are currently blank, so user-entered values are never overwritten.
let _geoAutofillRunning = false;
function _ttRecordsNeedingGeoFill() {
  return (RUN_META || []).filter((r) => {
    const lat = parseFloat(r.latitude),
      lon = parseFloat(r.longitude);
    if (isNaN(lat) || isNaN(lon)) return false;
    const noCountry = !String(r.sample_origin_country == null ? "" : r.sample_origin_country).trim();
    const noState = !String(
      r.sample_origin_state_province_territory == null ? "" : r.sample_origin_state_province_territory,
    ).trim();
    return noCountry || noState;
  });
}
function _ttAutofillGeoFromCoords(opts) {
  opts = opts || {};
  if (_geoAutofillRunning) return;
  const need = _ttRecordsNeedingGeoFill();
  if (!need.length) {
    if (opts.notify)
      alert("Nothing to fill — every sample with coordinates already has a country and state/territory.");
    return;
  }
  if (typeof d3 === "undefined" || typeof d3.geoContains !== "function") {
    if (opts.notify) alert("Geographic boundaries library (d3-geo) is unavailable.");
    return;
  }
  _geoAutofillRunning = true;
  const needState = need.some(
    (r) =>
      !String(r.sample_origin_state_province_territory == null ? "" : r.sample_origin_state_province_territory).trim(),
  );
  Promise.all([_geoCountries().catch(() => []), needState ? _geoAdmin1().catch(() => []) : Promise.resolve([])])
    .then(([countries, states]) => {
      let filled = 0;
      need.forEach((r) => {
        const lat = parseFloat(r.latitude),
          lon = parseFloat(r.longitude);
        const pt = [lon, lat];
        if (!String(r.sample_origin_country == null ? "" : r.sample_origin_country).trim() && countries.length) {
          const f = countries.find((ft) => {
            try {
              return d3.geoContains(ft, pt);
            } catch (e) {
              return false;
            }
          });
          if (f) {
            const p = f.properties || {};
            const nm = p.name || p.NAME || p.admin || p.ADMIN || p.name_long || p.NAME_LONG || p.sovereignt;
            if (nm) {
              r.sample_origin_country = String(nm);
              filled++;
            }
          }
        }
        if (
          !String(
            r.sample_origin_state_province_territory == null ? "" : r.sample_origin_state_province_territory,
          ).trim() &&
          states.length
        ) {
          const f = states.find((ft) => {
            try {
              return d3.geoContains(ft, pt);
            } catch (e) {
              return false;
            }
          });
          if (f) {
            const p = f.properties || {};
            const nm = p.name || p.name_en || p.gn_name || p.woe_name || p.name_alt;
            if (nm) {
              r.sample_origin_state_province_territory = String(nm);
              filled++;
            }
          }
        }
        // Compose a human-readable `location` string when still blank.
        if (!String(r.location == null ? "" : r.location).trim()) {
          const parts = [r.sample_origin_state_province_territory, r.sample_origin_country].filter(
            (x) => x && String(x).trim(),
          );
          if (parts.length) r.location = parts.join(", ");
        }
        if (typeof _normalizeMetaRecord === "function") _normalizeMetaRecord(r);
      });
      _geoAutofillRunning = false;
      if (filled) {
        const hasGeo = (RUN_META || []).some((r) => r.latitude != null && r.longitude != null);
        const mapBtn = document.getElementById("map-tab-btn");
        if (mapBtn) mapBtn.classList.toggle("hidden", !hasGeo);
        if (typeof _buildRunMetaTable === "function") _buildRunMetaTable();
        if (typeof _updateMetaSubTabStates === "function") _updateMetaSubTabStates();
        // Switch to the geo sub-tab whenever autofill just populated geographic
        // fields — the most useful thing to show after a lat/lon edit is the
        // updated choropleth.  Fall back to the previously-active sub-tab if
        // the geo sub-tab is still disabled (e.g. no boundary data available).
        const _geoSubBtn = document.querySelector('.meta-subtab[data-metasub="geo"]');
        const _targetSub = _geoSubBtn && !_geoSubBtn.disabled ? "geo" : _activeMetaSub;
        if (typeof _targetSub !== "undefined" && _targetSub && typeof _switchMetaSub === "function") {
          try {
            _switchMetaSub(_targetSub);
          } catch (e) {}
        }
        if (typeof _rebuildMapMarkers === "function") {
          try {
            _rebuildMapMarkers();
          } catch (e) {}
        }
        if (opts.notify) alert("Filled location fields for " + filled + " value(s) from coordinates.");
      } else if (opts.notify) {
        alert("Could not match any coordinates to a region (boundaries unavailable or points off-map).");
      }
    })
    .catch(() => {
      _geoAutofillRunning = false;
      if (opts.notify) alert("Could not load geographic boundaries (offline?).");
    });
}
window._ttAutofillGeoFromCoords = _ttAutofillGeoFromCoords;

// ── Per-region metadata breakdown for the choropleth tooltip ─────────────
// Sample-name set whose value for `field` matches the aggregated datum label.
function _geoSamplesForDatum(field, datumLabel) {
  const keys = new Set(_geoMatchKeys(field, datumLabel));
  const set = new Set();
  (RUN_META || []).forEach((r) => {
    const v = r[field];
    const vals = v == null ? [] : Array.isArray(v) ? v : [v];
    for (const x of vals) {
      for (const k of _geoMatchKeys(field, x)) {
        if (keys.has(k)) {
          set.add(r.sample_name);
          break;
        }
      }
    }
  });
  return set;
}
// Tally a metadata field's values across a set of samples → sorted [label,count].
function _geoMetaTally(sampleSet, field) {
  const m = {};
  (RUN_META || []).forEach((r) => {
    if (!sampleSet.has(r.sample_name)) return;
    const v = r[field];
    if (v == null) return;
    (Array.isArray(v) ? v : [v]).forEach((x) => {
      const s = String(x).trim();
      if (s) m[s] = (m[s] || 0) + 1;
    });
  });
  return Object.entries(m).sort((a, b) => b[1] - a[1]);
}
// Top detected organisms (visible detections) for a set of samples.
function _geoTopOrganisms(sampleSet, n) {
  const m = {};
  filteredData().forEach((r) => {
    const s = r["Specimen ID"] || "";
    if (!sampleSet.has(s)) return;
    const o = (r["Detected Organism"] || "").trim();
    if (o) m[o] = (m[o] || 0) + 1;
  });
  return Object.entries(m)
    .sort((a, b) => b[1] - a[1])
    .slice(0, n || 5);
}
// Per-organism detection stats (TASS / coverage / % reads / reads) across a
// set of samples — the data behind the pinned-region heatmap table.
function _geoRegionOrganismRows(sampleSet) {
  const by = {};
  filteredData().forEach((r) => {
    const s = r["Specimen ID"] || "";
    if (!sampleSet.has(s)) return;
    const name = (r["Detected Organism"] || "").trim();
    if (!name) return;
    const o = (by[name] = by[name] || { name, tass: 0, cov: 0, reads: 0, pct: 0, samples: new Set() });
    const t = parseFloat(r["TASS Score"]);
    if (!isNaN(t)) o.tass = Math.max(o.tass, t);
    const cands = [parseFloat(r["Breadth %"]), parseFloat(r["Genome Coverage %"]), parseFloat(r["Coverage"])].filter(
      (x) => !isNaN(x),
    );
    if (cands.length) o.cov = Math.max(o.cov, Math.max.apply(null, cands));
    const rd = parseFloat(r["# Reads Aligned"]);
    if (!isNaN(rd)) o.reads += rd;
    const pc = parseFloat(r["% Reads"]);
    if (!isNaN(pc)) o.pct = Math.max(o.pct, pc);
    o.samples.add(s);
  });
  return Object.values(by).sort((a, b) => b.tass - a.tass || b.reads - a.reads);
}

// ── Choropleth region pinning ────────────────────────────────────────────
// Clicking a country drills into its states, so pinning uses Shift/Ctrl-click
// or right-click. A pinned region shows a full info table in the ranked panel.
// Pins are independent floating windows (tooltip frozen in place), one per
// region, so several can be open at once for side-by-side comparison. Keyed
// by `${field}::${normalized feature name}`.
const _geoPins = new Map();
let _geoPinSeq = 0;
let _geoPinZ = 600;
function _geoPinMapKey(field, key) {
  return field + "::" + key;
}
function _geoIsPinnedFeature(field, f) {
  const k = _geoFeatureNames(field, f)[0] || "";
  return _geoPins.has(_geoPinMapKey(field, k));
}
// Aggregate metric stats for the samples in a region.
function _geoRegionStats(field, label) {
  const sset = _geoSamplesForDatum(field, label);
  let det = 0,
    reads = 0,
    tassSum = 0,
    tassN = 0;
  filteredData().forEach((r) => {
    if (!sset.has(r["Specimen ID"] || "")) return;
    det++;
    const rd = parseFloat(r["# Reads Aligned"]);
    if (!isNaN(rd)) reads += rd;
    const t = parseFloat(r["TASS Score"]);
    if (!isNaN(t)) {
      tassSum += t;
      tassN++;
    }
  });
  return { samples: sset, detections: det, reads, meanTass: tassN ? tassSum / tassN : 0 };
}
// Build the rich region-info HTML (same content as the hover tooltip, plus
// the detected-organism heatmap table). Used as the body of a pin window.
function _geoRegionInfoHTML(field, label, opts) {
  opts = opts || {};
  const maxOrg = opts.maxOrg || 0; // 0 = all
  const esc = (x) =>
    String(x == null ? "" : x)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  const stats = _geoRegionStats(field, label);
  const sset = stats.samples;
  const fmtN = (v) => Math.round(v).toLocaleString();
  const isCountry = field === "sample_origin_country";
  let countryLine = "";
  if (!isCountry) {
    const cs = new Set();
    (RUN_META || []).forEach((r) => {
      if (!sset.has(r.sample_name)) return;
      const c = String(r.sample_origin_country == null ? "" : r.sample_origin_country).trim();
      if (c) cs.add(c);
    });
    if (cs.size) countryLine = `<tr><td>Country</td><td>${[...cs].map(esc).join(", ")}</td></tr>`;
  }
  const rowsHtml = [];
  rowsHtml.push(`<tr><td>Samples</td><td>${sset.size}</td></tr>`);
  if (sset.size)
    rowsHtml.push(
      `<tr><td>Sample IDs</td><td>${[...sset].slice(0, 8).map(esc).join(", ")}${
        sset.size > 8 ? ` +${sset.size - 8}` : ""
      }</td></tr>`,
    );
  rowsHtml.push(`<tr><td>Detections</td><td>${fmtN(stats.detections)}</td></tr>`);
  rowsHtml.push(`<tr><td>Reads aligned</td><td>${fmtN(stats.reads)}</td></tr>`);
  rowsHtml.push(`<tr><td>Mean TASS</td><td>${stats.meanTass.toFixed(1)}</td></tr>`);
  const tally = (f2, lbl) => {
    const t = _geoMetaTally(sset, f2);
    if (t.length)
      rowsHtml.push(
        `<tr><td>${lbl}</td><td>${t
          .slice(0, 5)
          .map(([k, v]) => `${esc(k)} (${v})`)
          .join(", ")}${t.length > 5 ? `, +${t.length - 5}` : ""}</td></tr>`,
      );
  };
  tally("host_scientific_name", "Hosts");
  tally("host_disease", "Host disease");
  tally("environmental_site", "Env. site");
  if (isCountry) {
    const nested = _aggregateMetaNested(
      "sample_origin_country",
      "sample_origin_state_province_territory",
      _geoBarMetric || "detections",
    )[_geoNorm(label)];
    if (nested && nested.children && nested.children.length)
      rowsHtml.push(
        `<tr><td>States / territories</td><td>${nested.children
          .slice(0, 8)
          .map((c) => `${esc(c.label)}: ${(+c.value).toLocaleString()}`)
          .join("; ")}</td></tr>`,
      );
  }
  const orgRows = _geoRegionOrganismRows(sset);
  const tassBg = (v) => {
    if (isNaN(v) || v <= 0) return "";
    const t = Math.min(1, v / 100);
    return `background:hsla(${Math.round(210 - t * 65)},70%,45%,${(0.08 + t * 0.27).toFixed(3)})`;
  };
  const covBg = (v) => {
    if (isNaN(v) || v <= 0) return "";
    const t = Math.min(1, v / 100);
    return `background:hsla(145,60%,40%,${(0.06 + t * 0.3).toFixed(3)})`;
  };
  const pctBg = (v) => {
    if (isNaN(v) || v <= 0) return "";
    const t = Math.min(1, v / 100);
    return `background:hsla(28,90%,52%,${(0.06 + t * 0.34).toFixed(3)})`;
  };
  let orgTableHtml = "";
  if (orgRows.length) {
    const shown = maxOrg > 0 ? orgRows.slice(0, maxOrg) : orgRows;
    const body = shown
      .map(
        (o) =>
          `<tr><td class="org" title="${esc(o.name)}"><i>${esc(o.name)}</i></td>` +
          `<td style="${tassBg(o.tass)};text-align:center;font-weight:600">${o.tass ? o.tass.toFixed(1) : "—"}</td>` +
          `<td style="${covBg(o.cov)};text-align:center">${o.cov ? o.cov.toFixed(1) : "—"}</td>` +
          `<td style="${pctBg(o.pct)};text-align:right">${o.pct ? o.pct.toFixed(2) : "—"}</td>` +
          `<td style="text-align:right">${Math.round(o.reads).toLocaleString()}</td></tr>`,
      )
      .join("");
    const moreNote =
      maxOrg > 0 && orgRows.length > maxOrg
        ? `<div style="color:#9aa6b4;font-size:.72em;margin-top:.2em">…and ${
            orgRows.length - maxOrg
          } more — pin to see all</div>`
        : "";
    orgTableHtml =
      `<div style="font-weight:700;color:#607089;margin:.5em 0 .25em">Detected organisms (${orgRows.length})</div>` +
      `<div class="geo-org-scroll"><table class="geo-org-table">` +
      `<thead><tr><th style="text-align:left">Organism</th><th>TASS</th><th>Cov&nbsp;%</th><th>%&nbsp;Reads</th><th>Reads</th></tr></thead>` +
      `<tbody>${body}</tbody></table></div>${moreNote}`;
  } else {
    orgTableHtml = `<div style="color:#9aa6b4;font-size:.78em;margin-top:.4em">No detected organisms for this region in the current view.</div>`;
  }
  return (
    `<table class="geo-pin-table" style="width:100%;border-collapse:collapse;font-size:.78em">${countryLine}${rowsHtml.join(
      "",
    )}</table>` + orgTableHtml
  );
}
// Bring a pin window to the front.
function _geoRaisePin(win) {
  win.style.zIndex = String(++_geoPinZ);
}
// Make a pin window draggable by its header, clamped inside the map wrap.
function _geoMakePinDraggable(win, handle, wrap) {
  let sx = 0,
    sy = 0,
    ox = 0,
    oy = 0,
    dragging = false;
  const onMove = (e) => {
    if (!dragging) return;
    const r = wrap.getBoundingClientRect();
    let nx = ox + (e.clientX - sx);
    let ny = oy + (e.clientY - sy);
    nx = Math.max(0, Math.min(nx, r.width - win.offsetWidth));
    ny = Math.max(0, Math.min(ny, r.height - win.offsetHeight));
    win.style.left = nx + "px";
    win.style.top = ny + "px";
    e.preventDefault();
  };
  const onUp = () => {
    dragging = false;
    document.removeEventListener("pointermove", onMove);
    document.removeEventListener("pointerup", onUp);
  };
  handle.addEventListener("pointerdown", (e) => {
    if (e.target.closest(".geo-pin-win-close")) return;
    dragging = true;
    sx = e.clientX;
    sy = e.clientY;
    ox = parseFloat(win.style.left) || 0;
    oy = parseFloat(win.style.top) || 0;
    _geoRaisePin(win);
    document.addEventListener("pointermove", onMove);
    document.addEventListener("pointerup", onUp);
    e.preventDefault();
  });
}
function _geoRemovePin(mk) {
  const rec = _geoPins.get(mk);
  if (rec && rec.el && rec.el.parentNode) rec.el.parentNode.removeChild(rec.el);
  _geoPins.delete(mk);
  if (typeof _geoRedraw === "function") _geoRedraw();
}
function _geoCreatePinWindow(field, key, label, ev) {
  const wrap = document.getElementById("geo-map-wrap");
  if (!wrap) return;
  const mk = _geoPinMapKey(field, key);
  const esc = (x) =>
    String(x == null ? "" : x)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  const win = document.createElement("div");
  win.className = "geo-pin-win";
  const rect = wrap.getBoundingClientRect();
  const W = 250;
  let x = ev && ev.clientX != null ? ev.clientX - rect.left + 10 : 18 + (_geoPinSeq % 6) * 22;
  let y = ev && ev.clientY != null ? ev.clientY - rect.top + 10 : 18 + (_geoPinSeq % 6) * 22;
  _geoPinSeq++;
  win.style.left = Math.max(0, Math.min(x, rect.width - W - 6)) + "px";
  win.style.top = Math.max(0, Math.min(y, rect.height - 80)) + "px";
  win.innerHTML =
    `<div class="geo-pin-win-head"><span class="geo-pin-win-title"><i class="fas fa-thumbtack"></i> ${esc(
      label,
    )}</span><button type="button" class="geo-pin-win-close" title="Close">✕</button></div>` +
    `<div class="geo-pin-win-body">${_geoRegionInfoHTML(field, label)}</div>` +
    `<div class="geo-pin-win-grip" title="Drag to resize"></div>`;
  wrap.appendChild(win);
  _geoPins.set(mk, { field, key, label, el: win });
  win.querySelector(".geo-pin-win-close").addEventListener("click", () => _geoRemovePin(mk));
  win.addEventListener("pointerdown", () => _geoRaisePin(win));
  _geoMakePinDraggable(win, win.querySelector(".geo-pin-win-head"), wrap);
  _geoMakePinResizable(win, win.querySelector(".geo-pin-win-grip"), wrap);
  _geoRaisePin(win);
}
// Resize a pin window by dragging its bottom-right grip, clamped to the map.
function _geoMakePinResizable(win, grip, wrap) {
  if (!grip) return;
  let sx = 0,
    sy = 0,
    ow = 0,
    oh = 0,
    resizing = false;
  const onMove = (e) => {
    if (!resizing) return;
    const r = wrap.getBoundingClientRect();
    const left = parseFloat(win.style.left) || 0;
    const top = parseFloat(win.style.top) || 0;
    let nw = Math.max(180, Math.min(ow + (e.clientX - sx), r.width - left - 4));
    let nh = Math.max(110, Math.min(oh + (e.clientY - sy), r.height - top - 4));
    win.style.width = nw + "px";
    win.style.height = nh + "px";
    e.preventDefault();
  };
  const onUp = () => {
    resizing = false;
    document.removeEventListener("pointermove", onMove);
    document.removeEventListener("pointerup", onUp);
  };
  grip.addEventListener("pointerdown", (e) => {
    resizing = true;
    sx = e.clientX;
    sy = e.clientY;
    ow = win.offsetWidth;
    oh = win.offsetHeight;
    _geoRaisePin(win);
    document.addEventListener("pointermove", onMove);
    document.addEventListener("pointerup", onUp);
    e.preventDefault();
    e.stopPropagation();
  });
}
function _geoTogglePin(field, f, ev) {
  const key = _geoFeatureNames(field, f)[0] || "";
  if (!key) return;
  const mk = _geoPinMapKey(field, key);
  if (_geoPins.has(mk)) {
    _geoRemovePin(mk);
    return;
  }
  const p = f.properties || {};
  const label = p.name || p.NAME || p.name_en || p.admin || key;
  _geoCreatePinWindow(field, key, label, ev);
  if (typeof _geoRedraw === "function") _geoRedraw();
}
// Re-render the body of every open pin window (after a metric / data change).
function _geoRefreshAllPins() {
  _geoPins.forEach((rec) => {
    const body = rec.el && rec.el.querySelector(".geo-pin-win-body");
    if (body) body.innerHTML = _geoRegionInfoHTML(rec.field, rec.label);
  });
}

// parentField → childField nested aggregation, keyed by normalized parent.
function _aggregateMetaNested(parentField, childField, metric) {
  const toArr = (v) =>
    v == null
      ? []
      : Array.isArray(v)
      ? v.map((x) => String(x).trim()).filter(Boolean)
      : [String(v).trim()].filter(Boolean);
  const recBy = {};
  (RUN_META || []).forEach((r) => {
    const p = toArr(r[parentField]);
    if (!p.length) return;
    recBy[r.sample_name] = { parents: p, children: toArr(r[childField]) };
  });
  const parents = {};
  const eP = (p) =>
    (parents[p] = parents[p] || { samples: new Set(), detections: 0, reads: 0, tassSum: 0, tassN: 0, kids: {} });
  const eK = (P, c) => (P.kids[c] = P.kids[c] || { samples: new Set(), detections: 0, reads: 0, tassSum: 0, tassN: 0 });
  if (metric === "samples") {
    Object.entries(recBy).forEach(([s, o]) => {
      if (sampleHidden[s]) return;
      o.parents.forEach((p) => {
        const P = eP(p);
        P.samples.add(s);
        (o.children.length ? o.children : ["—"]).forEach((c) => eK(P, c).samples.add(s));
      });
    });
  } else {
    filteredData().forEach((row) => {
      const s = row["Specimen ID"] || "";
      const o = recBy[s];
      if (!o) return;
      const reads = parseFloat(row["# Reads Aligned"]);
      const tass = parseFloat(row["TASS Score"]);
      o.parents.forEach((p) => {
        const P = eP(p);
        P.samples.add(s);
        P.detections += 1;
        if (!isNaN(reads)) P.reads += reads;
        if (!isNaN(tass)) {
          P.tassSum += tass;
          P.tassN += 1;
        }
        (o.children.length ? o.children : ["—"]).forEach((c) => {
          const K = eK(P, c);
          K.samples.add(s);
          K.detections += 1;
          if (!isNaN(reads)) K.reads += reads;
          if (!isNaN(tass)) {
            K.tassSum += tass;
            K.tassN += 1;
          }
        });
      });
    });
  }
  const val = (g) =>
    metric === "samples"
      ? g.samples.size
      : metric === "detections"
      ? g.detections
      : metric === "reads"
      ? g.reads
      : g.tassN
      ? g.tassSum / g.tassN
      : 0;
  const out = {};
  Object.entries(parents).forEach(([pLabel, P]) => {
    const children = Object.entries(P.kids)
      .map(([c, g]) => ({ label: c, value: val(g) }))
      .filter((d) => d.label !== "—")
      .sort((a, b) => b.value - a.value);
    out[_geoNorm(pLabel)] = { label: pLabel, value: val(P), children };
  });
  return out;
}

let _geoNested = {};
function _geoRenderLegend(el, color, maxV, metric) {
  if (!el) return;
  const fmt = (v) => (metric === "tass" ? (+v).toFixed(1) : Math.round(v).toLocaleString());
  const stops = [0, 0.25, 0.5, 0.75, 1].map((t) => color(t * maxV));
  el.style.display = "block";
  el.innerHTML =
    `<div style="font-weight:600;margin-bottom:1px">${_META_METRIC_LABELS[metric] || metric}</div>` +
    `<div class="geo-legend-bar" style="background:linear-gradient(to right, ${stops.join(",")})"></div>` +
    `<div class="geo-legend-row"><span>0</span><span>${fmt(maxV)}</span></div>`;
}
function _geoHover(ev, f, field, metric, metaByNorm, stateCountries) {
  const p = f.properties || {};
  const nm = p.name || p.NAME || p.name_en || p.admin || "—";
  const lbl = _META_METRIC_LABELS[metric] || metric;
  const fmt = (v) => (metric === "tass" ? (+v).toFixed(1) : Math.round(v).toLocaleString());
  const m = _geoFeatureMatch(field, f, metaByNorm, stateCountries);
  const _esc = (s) =>
    String(s == null ? "" : s)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  const _lblSpan = (t) => `<br><span style="color:#9fb0c3;font-size:0.9em">${t}:</span> `;
  let html = `<b>${_esc(nm)}</b>`;
  if (m) {
    html += `<br>${lbl}: <b>${fmt(m.value)}</b>`;
    // Which samples sit in this region, and their key metadata.
    const sset = _geoSamplesForDatum(field, m.label);
    if (sset.size) {
      const slist = [...sset];
      html += _lblSpan(`Samples (${sset.size})`) + _esc(slist.slice(0, 5).join(", "));
      if (slist.length > 5) html += ` +${slist.length - 5} more`;

      const hosts = _geoMetaTally(sset, "host_scientific_name");
      if (hosts.length)
        html +=
          _lblSpan("Hosts") +
          hosts
            .slice(0, 3)
            .map(([k, v]) => `${_esc(k)} (${v})`)
            .join(", ");

      const diseases = _geoMetaTally(sset, "host_disease");
      if (diseases.length)
        html +=
          _lblSpan("Host disease") +
          diseases
            .slice(0, 4)
            .map(([k, v]) => `${_esc(k)} (${v})`)
            .join(", ");

      const sites = _geoMetaTally(sset, "environmental_site");
      if (sites.length)
        html +=
          _lblSpan("Environmental site") +
          sites
            .slice(0, 3)
            .map(([k, v]) => `${_esc(k)} (${v})`)
            .join(", ");

      const orgs = _geoTopOrganisms(sset, 4);
      if (orgs.length) html += _lblSpan("Top organisms") + orgs.map(([k, v]) => `<i>${_esc(k)}</i> (${v})`).join(", ");
    }
    // Nested child regions (states within a country, etc.).
    const nested = _geoNested[_geoNorm(m.label)];
    if (nested && nested.children && nested.children.length) {
      const childName = field === "sample_origin_country" ? "States / territories" : "Sub-regions";
      html += _lblSpan(childName);
      html += nested.children
        .slice(0, 6)
        .map((c) => `${_esc(c.label)}: ${fmt(c.value)}`)
        .join("&nbsp;&nbsp;•&nbsp;&nbsp;");
      if (nested.children.length > 6) html += `&nbsp;&nbsp;…and ${nested.children.length - 6} more`;
    }
  } else {
    html += `<br><span style="color:#9fb0c3">no data</span>`;
  }
  // Clear, always-visible hint on how to pin this region into a standalone,
  // draggable / resizable info window.
  const pinHint =
    field === "sample_origin_country"
      ? "Shift-click (or right-click) to pin &nbsp;·&nbsp; click to drill into states"
      : "Click (or right-click) to pin this region";
  html +=
    `<div style="margin-top:.5em;padding-top:.4em;border-top:1px solid rgba(255,255,255,.18);` +
    `color:#ffd54f;font-size:.84em;font-weight:600">` +
    `<i class="fas fa-thumbtack" style="margin-right:.35em"></i>${pinHint}</div>`;
  showTip(html, ev);
}
function _geoDrawChoropleth(field, metric, ctx) {
  ctx = ctx || {};
  const wrap = document.getElementById("geo-map-wrap");
  const svgEl = document.getElementById("geo-map-svg");
  const statusEl = document.getElementById("geo-map-status");
  const legendEl = document.getElementById("geo-legend");
  if (!wrap || !svgEl || typeof d3.geoNaturalEarth1 !== "function") return;
  const W = wrap.clientWidth || 760;
  const H = wrap.clientHeight || 560;
  const setStatus = (msg) => {
    if (!statusEl) return;
    statusEl.style.display = msg ? "flex" : "none";
    statusEl.textContent = msg || "";
  };
  const flat = _aggregateMetaByField(field, metric, ctx.sampleSet);
  const metaByNorm = new Map();
  flat.forEach((d) => {
    _geoMatchKeys(field, d.label).forEach((k) => {
      if (!metaByNorm.has(k)) metaByNorm.set(k, d);
    });
  });
  const maxV = d3.max(flat, (d) => d.value) || 1;
  const color = d3.scaleSequential(d3.interpolateYlOrRd).domain([0, maxV]);
  const childField = !ctx.drill && field === "sample_origin_country" ? "sample_origin_state_province_territory" : null;
  _geoNested = childField ? _aggregateMetaNested(field, childField, metric) : {};

  const svg = d3.select(svgEl).attr("viewBox", `0 0 ${W} ${H}`).attr("preserveAspectRatio", "xMidYMid meet");
  svg.selectAll("*").remove();
  if (legendEl) legendEl.style.display = "none";
  const isState = field === "sample_origin_state_province_territory";
  setStatus("Loading boundaries…");
  const loader = isState ? _geoAdmin1() : _geoCountries();
  loader
    .then((all) => {
      if (!all || !all.length) {
        setStatus("Map boundaries unavailable — use the ranked list →");
        return;
      }
      setStatus("");
      // When drilling into one country, restrict admin-1 polygons to it and
      // zoom the projection to those regions.
      let features = all;
      if (ctx.drill) {
        const sub = all.filter((f) => _geoFeatureCountryNorm(f) === ctx.drill.norm);
        if (sub.length) features = sub;
      }
      // State/territory → country constraint (so "GA" doesn't light up both
      // Georgia-US and Gorno-Altay-RU). Built from the metadata; when drilled
      // into a country the features are already restricted so it's redundant.
      const stateCountries = isState ? _geoStateCountryMap() : null;
      if (isState && !ctx.drill) _geoWarnAmbiguousStates(flat, all, stateCountries);
      const fc = { type: "FeatureCollection", features };
      const proj = d3.geoNaturalEarth1().fitSize([W, H], fc);
      const path = d3.geoPath(proj);
      const sel = svg
        .append("g")
        .selectAll("path")
        .data(features)
        .join("path")
        .attr("class", "geo-region")
        .attr("d", path)
        .attr("stroke", (f) => (_geoIsPinnedFeature(field, f) ? "#1565c0" : "#ffffff"))
        .attr("stroke-width", (f) => (_geoIsPinnedFeature(field, f) ? 2.2 : isState ? 0.3 : 0.4))
        .attr("fill", (f) => {
          const m = _geoFeatureMatch(field, f, metaByNorm, stateCountries);
          return m ? color(m.value) : "#dfe6ef";
        })
        .on("mouseover", (ev, f) => _geoHover(ev, f, field, metric, metaByNorm, stateCountries))
        .on("mousemove", moveTip)
        .on("mouseout", hideTip);
      // Right-click pins any region (drill-free way to inspect a country).
      sel.on("contextmenu", (ev, f) => {
        ev.preventDefault();
        hideTip();
        _geoTogglePin(field, f, ev);
      });
      // Click behaviour:
      //  • Country (top level): plain click drills into its states; Shift/
      //    Ctrl/⌘-click pins (so you can still inspect a country without zooming).
      //  • State / territory (or drilled view): plain click OR Shift-click pins,
      //    since there's nothing deeper to drill into.
      const _canDrill = field === "sample_origin_country" && !ctx.drill;
      sel.style("cursor", "pointer").on("click", (ev, f) => {
        if (!_canDrill || ev.shiftKey || ev.ctrlKey || ev.metaKey) {
          hideTip();
          _geoTogglePin(field, f, ev);
          return;
        }
        // Country, plain click → drill.
        hideTip();
        let n = _geoFeatureNames(field, f)[0] || "";
        n = _COUNTRY_ALIASES[n] || n;
        const p = f.properties || {};
        const label = p.name || p.NAME || p.admin || n;
        _geoDrill = { norm: n, label };
        if (typeof _geoRedraw === "function") _geoRedraw();
      });
      _geoRenderLegend(legendEl, color, maxV, metric);
      // Refresh any open pin windows against the current data/metric.
      _geoRefreshAllPins();
    })
    .catch(() => {
      setStatus("Map boundaries unavailable (offline?) — use the ranked list →");
      if (legendEl) legendEl.style.display = "none";
    });
}

// Current ranked-bar data, kept so the overlay can re-render to fit on resize
// without recomputing the aggregation.
let _geoBarData = [];
let _geoBarMetric = "detections";
// Drill state: when set, the map is showing one country's states/territories.
let _geoDrill = null;
let _geoRedraw = null; // the active draw() closure, callable from map clicks
function _geoFitOverlayBar() {
  const body = document.getElementById("geo-overlay-body");
  if (!body) return;
  _renderMetaBarChart("geo-cmp-chart", _geoBarData, _geoBarMetric, {
    width: body.clientWidth - 2,
    fitHeight: body.clientHeight - 2,
  });
}
// Make the overlay panel draggable (by its header) and resizable (corner
// grip), constrained to stay within the map div. Re-fits the bar on resize.
function _geoMakeOverlayInteractive() {
  const wrap = document.getElementById("geo-map-wrap");
  const ov = document.getElementById("geo-overlay");
  const head = document.getElementById("geo-overlay-toggle");
  const grip = document.getElementById("geo-overlay-grip");
  if (!wrap || !ov || !head || ov._interactive) return;
  ov._interactive = true;
  let mode = null,
    sx = 0,
    sy = 0,
    ox = 0,
    oy = 0,
    ow = 0,
    oh = 0,
    moved = false;
  const onMove = (e) => {
    if (!mode) return;
    const wr = wrap.getBoundingClientRect();
    const dx = e.clientX - sx;
    const dy = e.clientY - sy;
    if (Math.abs(dx) + Math.abs(dy) > 3) moved = true;
    if (mode === "move") {
      let nl = Math.max(0, Math.min(ox + dx, wr.width - ow));
      let nt = Math.max(0, Math.min(oy + dy, wr.height - oh));
      ov.style.left = nl + "px";
      ov.style.top = nt + "px";
      ov.style.right = "auto";
    } else {
      let nw = Math.max(190, Math.min(ow + dx, wr.width - ox));
      let nh = Math.max(150, Math.min(oh + dy, wr.height - oy));
      ov.style.width = nw + "px";
      ov.style.height = nh + "px";
      _geoFitOverlayBar();
    }
    e.preventDefault();
  };
  const onUp = () => {
    if (mode === "resize") _geoFitOverlayBar();
    mode = null;
    document.removeEventListener("pointermove", onMove);
    document.removeEventListener("pointerup", onUp);
  };
  const start = (e, m) => {
    mode = m;
    moved = false;
    const r = ov.getBoundingClientRect();
    const wr = wrap.getBoundingClientRect();
    sx = e.clientX;
    sy = e.clientY;
    ox = r.left - wr.left;
    oy = r.top - wr.top;
    ow = r.width;
    oh = r.height;
    // pin to left/top so movement is absolute within the map div
    ov.style.left = ox + "px";
    ov.style.top = oy + "px";
    ov.style.right = "auto";
    document.addEventListener("pointermove", onMove);
    document.addEventListener("pointerup", onUp);
    e.preventDefault();
  };
  head.addEventListener("pointerdown", (e) => {
    if (e.target.closest(".geo-overlay-grip")) return;
    start(e, "move");
  });
  if (grip) grip.addEventListener("pointerdown", (e) => start(e, "resize"));
  // Collapse on a genuine click (not the end of a drag).
  head.addEventListener("click", () => {
    if (moved) {
      moved = false;
      return;
    }
    ov.classList.toggle("collapsed");
    if (!ov.classList.contains("collapsed")) _geoFitOverlayBar();
  });
}

function _buildGeoComparison() {
  const section = document.getElementById("geo-cmp-section");
  if (!section) return;
  if (typeof _updateMetaSubTabStates === "function") _updateMetaSubTabStates();
  if (!_anyMetaValue(_GEO_META_FIELDS)) return;
  const levelSel = document.getElementById("geo-cmp-level");
  const metricSel = document.getElementById("geo-cmp-metric");
  const noData = document.getElementById("geo-cmp-no-data");
  const lvlLbl = document.getElementById("geo-overlay-level-lbl");
  const LBL = {
    sample_origin_country: "Countries",
    sample_origin_state_province_territory: "States / Territories",
  };
  const backBtn = document.getElementById("geo-back");
  const draw = () => {
    let field = (levelSel && levelSel.value) || "sample_origin_country";
    const metric = (metricSel && metricSel.value) || "detections";
    // A drill overrides to state level, restricted to the picked country.
    let sampleSet = null;
    let drill = null;
    if (_geoDrill) {
      field = "sample_origin_state_province_territory";
      if (levelSel) levelSel.value = field;
      sampleSet = _geoCountrySampleSet(_geoDrill.norm);
      drill = _geoDrill;
    }
    const data = _aggregateMetaByField(field, metric, sampleSet);
    if (noData) noData.style.display = data.length ? "none" : "block";
    if (lvlLbl) lvlLbl.textContent = drill ? `${drill.label} — states` : LBL[field] || "Regions";
    if (backBtn) backBtn.style.display = drill ? "inline-flex" : "none";
    _geoBarData = data;
    _geoBarMetric = metric;
    _geoFitOverlayBar();
    _geoDrawChoropleth(field, metric, { drill, sampleSet });
  };
  _geoRedraw = draw;
  if (!section._wired) {
    section._wired = true;
    // Manually changing the level clears any active country drill.
    if (levelSel)
      levelSel.addEventListener("change", () => {
        _geoDrill = null;
        draw();
      });
    if (metricSel) metricSel.addEventListener("change", draw);
    if (backBtn)
      backBtn.addEventListener("click", () => {
        _geoDrill = null;
        if (levelSel) levelSel.value = "sample_origin_country";
        draw();
      });
    _geoMakeOverlayInteractive();
  }
  draw();
}

function _buildHostDisease() {
  const section = document.getElementById("host-cmp-section");
  if (!section) return;
  if (typeof _updateMetaSubTabStates === "function") _updateMetaSubTabStates();
  if (!_anyMetaValue(_HOST_META_FIELDS)) return;
  const fieldSel = document.getElementById("host-cmp-field");
  const metricSel = document.getElementById("host-cmp-metric");
  const noData = document.getElementById("host-cmp-no-data");
  const draw = () => {
    const field = (fieldSel && fieldSel.value) || "host_disease";
    const metric = (metricSel && metricSel.value) || "detections";
    const data = _aggregateMetaByField(field, metric);
    if (noData) noData.style.display = data.length ? "none" : "block";
    _renderMetaBarChart("host-cmp-chart", data, metric);
  };
  if (!section._wired) {
    section._wired = true;
    if (fieldSel) fieldSel.addEventListener("change", draw);
    if (metricSel) metricSel.addEventListener("change", draw);
  }
  draw();
  if (typeof _buildHostMatrix === "function") _buildHostMatrix();
}

/* ── Symptom × Organism matrix + drill-down ──────────────────────────────
       -   For each symptom (host_disease value) counts the number of distinct
       -   visible samples in which each organism was detected — e.g. how many
       -   "cold-like symptoms" samples contain RSV, or how many "diarrhea"
       -   samples contain E. coli. Rows = symptoms, columns = organisms, cell =
       -   sample count. Clicking a symptom row drills into a per-organism bar.
      ──────────────────────────────────────────────────────────────────────── */
function _hostMatrixData(symField, orgField) {
  const toArr = (v) =>
    v == null
      ? []
      : Array.isArray(v)
      ? v.map((x) => String(x).trim()).filter(Boolean)
      : [String(v).trim()].filter(Boolean);
  const symBy = {};
  (RUN_META || []).forEach((r) => {
    const vals = toArr(r[symField]);
    if (vals.length) symBy[r.sample_name] = vals;
  });
  if (!Object.keys(symBy).length) return null;
  // sample → set of detected organisms (visible detections only)
  const orgBy = {};
  filteredData().forEach((row) => {
    const s = row["Specimen ID"] || "";
    if (!symBy[s]) return;
    const o = String(row[orgField] == null ? "" : row[orgField]).trim();
    if (!o) return;
    (orgBy[s] = orgBy[s] || new Set()).add(o);
  });
  // symptom → organism → Set(samples)
  const cells = {};
  Object.entries(symBy).forEach(([s, syms]) => {
    if (sampleHidden[s]) return;
    const orgs = orgBy[s];
    if (!orgs || !orgs.size) return;
    syms.forEach((sym) => {
      cells[sym] = cells[sym] || {};
      orgs.forEach((o) => {
        (cells[sym][o] = cells[sym][o] || new Set()).add(s);
      });
    });
  });
  const symList = [];
  const orgCount = {};
  Object.entries(cells).forEach(([sym, om]) => {
    let t = 0;
    Object.entries(om).forEach(([o, set]) => {
      const n = set.size;
      om[o] = n; // collapse Set → count
      t += n;
      orgCount[o] = (orgCount[o] || 0) + n;
    });
    symList.push({ sym, total: t });
  });
  if (!symList.length || !Object.keys(orgCount).length) return null;
  symList.sort((a, b) => b.total - a.total);
  const orgList = Object.entries(orgCount)
    .map(([o, n]) => ({ o, n }))
    .sort((a, b) => b.n - a.n);
  return { cells, symList, orgList };
}

function _renderHostMatrix(model, topN, onRowClick) {
  const wrap = document.getElementById("host-matrix-chart");
  if (!wrap) return;
  wrap.innerHTML = "";
  const orgs = topN > 0 ? model.orgList.slice(0, topN) : model.orgList;
  let maxCell = 1;
  model.symList.forEach((s) =>
    orgs.forEach((o) => {
      const v = model.cells[s.sym][o.o] || 0;
      if (v > maxCell) maxCell = v;
    }),
  );
  const color = d3.scaleSequential(d3.interpolateBlues).domain([0, maxCell * 1.15]);
  const tbl = document.createElement("table");
  tbl.className = "hx-matrix";
  const thead = document.createElement("thead");
  const hr = document.createElement("tr");
  const corner = document.createElement("th");
  corner.className = "hx-corner";
  corner.textContent = "Symptom \\ Organism";
  hr.appendChild(corner);
  orgs.forEach((o) => {
    const th = document.createElement("th");
    th.className = "hx-col";
    const d = document.createElement("div");
    d.textContent = o.o;
    d.title = o.o + " — " + o.n + " sample-hits total";
    th.appendChild(d);
    hr.appendChild(th);
  });
  thead.appendChild(hr);
  tbl.appendChild(thead);
  const tb = document.createElement("tbody");
  model.symList.forEach((s) => {
    const tr = document.createElement("tr");
    const rh = document.createElement("th");
    rh.className = "hx-row";
    rh.textContent = s.sym;
    rh.title = s.sym + " — click to drill in";
    rh.addEventListener("click", () => {
      tb.querySelectorAll("th.hx-row").forEach((x) => x.classList.remove("hx-active"));
      rh.classList.add("hx-active");
      onRowClick(s.sym);
    });
    tr.appendChild(rh);
    orgs.forEach((o) => {
      const td = document.createElement("td");
      td.className = "hx-cell";
      const v = model.cells[s.sym][o.o] || 0;
      if (v > 0) {
        td.style.background = color(v);
        const sp = document.createElement("span");
        if (v > maxCell * 0.6) sp.style.color = "#fff";
        sp.textContent = v;
        td.appendChild(sp);
        td.title = `${s.sym} × ${o.o}: ${v} sample${v === 1 ? "" : "s"}`;
      }
      tr.appendChild(td);
    });
    tb.appendChild(tr);
  });
  tbl.appendChild(tb);
  wrap.appendChild(tbl);
}

function _hostDrill(model, sym) {
  const box = document.getElementById("host-drill");
  const lbl = document.getElementById("host-drill-label");
  if (!box) return;
  const om = model.cells[sym] || {};
  const data = Object.entries(om)
    .map(([o, n]) => ({ label: o, value: n }))
    .sort((a, b) => b.value - a.value);
  if (lbl) lbl.textContent = sym;
  box.style.display = data.length ? "block" : "none";
  _renderMetaBarChart("host-drill-chart", data, "samples");
}

function _buildHostMatrix() {
  const section = document.getElementById("host-matrix-section");
  if (!section) return;
  const fieldSel = document.getElementById("host-mx-field");
  const levelSel = document.getElementById("host-mx-level");
  const topnSel = document.getElementById("host-mx-topn");
  const noData = document.getElementById("host-matrix-no-data");
  const drill = document.getElementById("host-drill");
  const draw = () => {
    const symField = (fieldSel && fieldSel.value) || "host_disease";
    const orgField = (levelSel && levelSel.value) || "Detected Organism";
    const topN = parseInt((topnSel && topnSel.value) || "15", 10) || 0;
    if (drill) drill.style.display = "none";
    const model = _hostMatrixData(symField, orgField);
    const ok = !!(model && model.symList.length && model.orgList.length);
    if (noData) noData.style.display = ok ? "none" : "block";
    const wrap = document.getElementById("host-matrix-chart");
    if (wrap && !ok) wrap.innerHTML = "";
    if (!ok) return;
    _renderHostMatrix(model, topN, (sym) => _hostDrill(model, sym));
  };
  if (!section._wired) {
    section._wired = true;
    [fieldSel, levelSel, topnSel].forEach((sel) => sel && sel.addEventListener("change", draw));
  }
  draw();
}
