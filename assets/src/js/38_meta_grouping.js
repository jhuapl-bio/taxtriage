/* ═══════════════════════════════════════════════════════════════════════════
       -  §  ARBITRARY METADATA GROUPING  (shared "Group by" engine)
       -     One grouping selection, shared by every metadata view. The user
       -     picks N metadata columns; samples are bucketed into COMPOSITE
       -     groups keyed by the union of those columns' values, e.g.
       -         ["environmental_site", "sequencing_platform"]
       -           → "Agricultural │ MiSeq", "Water │ MinION", "Soil │ MiSeq"
       -     Any column present on RUN_META is eligible — the standard
       -     AMD-P/Talos fields and whatever arbitrary columns arrived from a
       -     samplesheet or an uploaded meta.csv. Groupability is decided
       -     heuristically (see _mgColumnStats) so free-text notes and
       -     coordinates don't clutter the picker.
       -
       -     Public surface (all on window, all null-safe before init):
       -       metaGrouping.fields()          → active field list (array)
       -       metaGrouping.setFields(a)      → set + notify + redraw
       -       metaGrouping.active()          → true when ≥1 field chosen
       -       metaGrouping.candidates()      → [{key,label,distinct,coverage}]
       -       metaGrouping.groupOf(sample)   → composite key | null
       -       metaGrouping.groupsOf(sample)  → [key] (multi-valued fields)
       -       metaGrouping.model()           → full {groups, order, bySample}
       -       metaGrouping.color(key)        → stable hex
       -       metaGrouping.shape(key)        → stable shape id
       -       metaGrouping.onChange(fn)      → subscribe
       -
       -     Multi-valued cells (host_disease is split to an array upstream by
       -     _normalizeMetaRecord) fan a sample into EVERY combination, which
       -     is what makes "compare Agricultural vs Water vs Soil all at once"
       -     work when a sample legitimately belongs to more than one bucket.
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  "use strict";

  // Separator between the parts of a composite key. U+2502 is used rather than
  // a plain "|" so a value that itself contains a pipe can't forge a key.
  const SEP = " │ ";

  // Columns that are never offered as a grouping dimension: identifiers (one
  // group per sample is meaningless), coordinates (continuous), and the
  // internal bookkeeping keys the report attaches to metadata records.
  const EXCLUDE = new Set([
    "sample_name",
    "sample",
    "sample_id",
    "specimen",
    "latitude",
    "longitude",
    "__mergedmembers",
    "__index",
  ]);

  // Fields that are *known* to be good grouping dimensions. These sort to the
  // top of the picker and keep their curated labels; everything else is
  // discovered at runtime.
  //
  // Being listed here is not just about ordering — curated fields are exempt
  // from the "one group per sample" and prose guards in _mgColumnStats. That
  // matters for site columns in particular: a survey with one site per sample
  // is completely normal, and the generic `site` column would otherwise be
  // rejected as an identifier on exactly the runs where it is most useful.
  //
  // `environmental_site` is the AMD-P/Talos standard name; the rest are the
  // generic spellings submitters actually use for the same concept.
  const PREFERRED = [
    // Published by the map's region drawing tool (44_map_regions.js) when the
    // user turns on "Group by region" — curated so it is never rejected as an
    // identifier on a run where every sample sits in its own region.
    "map_region",
    "sample_type",
    "specimen_type",
    "environmental_site",
    "site",
    "site_name",
    "sampling_site",
    "collection_site",
    "site_id",
    "site_type",
    "host_scientific_name",
    "host_disease",
    "sample_origin_country",
    "sample_origin_state_province_territory",
    "sequencing_platform",
    "sequencing_instrument",
    "library_preparation_kit",
    "sequencing_protocol_primer_set",
    "submitter_organization_name",
    "organism",
    "run_id",
    "location",
    "collection_time",
  ];

  // Marker shapes cycled alongside the palette. Combining colour AND shape
  // keeps groups distinguishable past the 15-colour palette and survives
  // greyscale printing / colour-vision deficiency.
  const SHAPES = ["circle", "square", "triangle", "diamond", "star", "hexagon", "cross", "wedge"];

  // ── internal state ──────────────────────────────────────────────────────
  let _fields = []; // active grouping columns, in user-chosen order
  let _colorBy = {}; // composite key → hex
  let _shapeBy = {}; // composite key → shape id
  /* ── Per-group display state ────────────────────────────────────────────
     Each group cycles through three states as the user clicks its chip:

       "normal"     (default) — drawn like everything else
         ↓ click
       "highlight"            — this group is emphasised, every OTHER group
                                fades back. The map additionally draws the
                                similarity network between its members.
         ↓ click
       "hidden"               — removed from the views and the map entirely
         ↓ click
       "normal"               — back to the start

     Only the non-default states are stored, so the default (empty map) shows
     everything, and a group that appears later — from an upload or an edited
     cell — starts normal without touching this state.

     Deliberately separate from the sidebar's `sampleHidden`: that removes a
     sample from the whole analysis; this only affects the grouped views.

     Several groups can be highlighted at once (they are compared against each
     other, and everything else fades) — the cycle is per-chip, not exclusive. */
  const STATE_NORMAL = "normal";
  const STATE_HIGHLIGHT = "highlight";
  const STATE_HIDDEN = "hidden";
  const CYCLE = [STATE_NORMAL, STATE_HIGHLIGHT, STATE_HIDDEN];

  let _groupState = {}; // key → "highlight" | "hidden"  (normal is absent)
  // Bumped on every display-state mutation. A size-based cache token would
  // treat "hide A, show A, hide B" as unchanged and serve a stale list.
  let _hiddenSeq = 0;

  function _stateOf(key) {
    return _groupState[key] || STATE_NORMAL;
  }
  function _hiddenKeys() {
    return Object.keys(_groupState).filter((k) => _groupState[k] === STATE_HIDDEN);
  }
  function _highlightKeys() {
    return Object.keys(_groupState).filter((k) => _groupState[k] === STATE_HIGHLIGHT);
  }
  const _listeners = [];
  let _modelCache = null;
  let _modelKey = ""; // cache token: fields + merge state + hidden samples

  function _palette() {
    return typeof PALETTE !== "undefined" && Array.isArray(PALETTE) && PALETTE.length
      ? PALETTE
      : ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#117733", "#882255"];
  }

  function _label(key) {
    return typeof _metaKeyLabel === "function" ? _metaKeyLabel(key) : String(key).replace(/_/g, " ");
  }

  function _records() {
    return typeof RUN_META !== "undefined" && Array.isArray(RUN_META) ? RUN_META : [];
  }

  // Values of one field on one record, always as an array of trimmed strings.
  function _cellValues(rec, field) {
    const v = rec ? rec[field] : null;
    if (v == null) return [];
    const arr = Array.isArray(v) ? v : [v];
    const out = [];
    arr.forEach((x) => {
      const s = String(x).trim();
      if (s && s.toLowerCase() !== "null" && s.toLowerCase() !== "na" && s.toLowerCase() !== "nan") out.push(s);
    });
    return out;
  }

  /* ── Heuristic column profiling ────────────────────────────────────────
     Decides which columns are worth offering as a grouping dimension. The
     failure mode to avoid in both directions:
       • too permissive → the picker fills with free-text notes columns whose
         "groups" are one sample each, which is just an unglamorous no-op;
       • too strict → the arbitrary columns the user actually cares about
         (their own site codes, plate numbers, project names) get hidden, and
         the whole point of arbitrary grouping is lost.

     A column qualifies when it:
       • is not an identifier / coordinate (EXCLUDE),
       • has ≥ 2 distinct non-empty values (one group is not a comparison),
       • does not give a group per sample — a column where (nearly) every row
         is unique produces no aggregation at all. Curated fields are exempt:
         `sample_origin_country` on a 9-country run is all-unique but still a
         meaningful way to colour a map.
       • does not look like prose (long values / many words per value), which
         is the signature of a notes or description column,
       • is not continuous-numeric (a measurement or a numeric ID). Low-
         cardinality numerics — plate, replicate, year — are kept, because
         those are real grouping dimensions.
     Returns [{key,label,distinct,coverage,values,preferred}] sorted with the
     curated fields first, then by coverage, then alphabetically.               */
  const MAX_DISTINCT = 50;

  function _mgColumnStats() {
    const recs = _records();
    if (!recs.length) return [];
    const keys = new Set();
    recs.forEach((r) => Object.keys(r || {}).forEach((k) => keys.add(k)));

    const out = [];
    keys.forEach((k) => {
      const kl = String(k).toLowerCase();
      if (EXCLUDE.has(kl)) return;
      if (kl.startsWith("__")) return;

      const values = new Set();
      let filled = 0;
      let numeric = 0;
      let cells = 0;
      let charSum = 0;
      let wordSum = 0;
      recs.forEach((r) => {
        const vals = _cellValues(r, k);
        if (!vals.length) return;
        filled += 1;
        vals.forEach((v) => {
          values.add(v);
          cells += 1;
          charSum += v.length;
          wordSum += v.split(/\s+/).length;
          if (!isNaN(Number(v))) numeric += 1;
        });
      });

      const distinct = values.size;
      if (distinct < 2) return; // nothing to compare
      if (distinct > MAX_DISTINCT && distinct > filled / 2) return;

      const pi = PREFERRED.indexOf(kl);
      const preferred = pi >= 0;

      // Continuous-numeric guard applies to every column: a lat-like or
      // measurement column is never a category regardless of its name.
      if (numeric >= cells * 0.9 && distinct > 12) return;

      if (!preferred) {
        // One group per sample → no aggregation happens. Needs a few rows
        // before it can be judged; a 2-sample run is legitimately all-unique.
        if (filled >= 4 && distinct >= filled * 0.9) return;
        // Prose guard: notes / descriptions / comments read as sentences.
        const avgChars = cells ? charSum / cells : 0;
        const avgWords = cells ? wordSum / cells : 0;
        if (avgWords > 2.5 || avgChars > 32) return;
      }

      out.push({
        key: k,
        label: _label(k),
        distinct,
        coverage: filled / recs.length,
        values: Array.from(values).sort((a, b) => a.localeCompare(b)),
        preferred,
        _rank: preferred ? pi : 999,
      });
    });

    out.sort((a, b) => a._rank - b._rank || b.coverage - a.coverage || a.label.localeCompare(b.label));
    return out;
  }

  /* ── Identity resolution ───────────────────────────────────────────────
     Grouping keys are built against the identity the report currently shows:
     the library sample_name normally, the specimen when merge is on. Values
     are unioned across a specimen's libraries unless the merge dialog stored
     an explicit resolution for that field, matching _metaValuesByEntity.      */
  function _viewId(sample) {
    return typeof specimenMergeEnabled !== "undefined" && specimenMergeEnabled && typeof specimenOf === "function"
      ? specimenOf(sample)
      : sample;
  }

  function _resolvedOverride(id, field) {
    if (typeof SPECIMEN_META_RESOLVED === "undefined" || !SPECIMEN_META_RESOLVED[id]) return null;
    const v = SPECIMEN_META_RESOLVED[id][field];
    if (v == null || v === "") return null;
    return (Array.isArray(v) ? v : [v]).map((x) => String(x).trim()).filter(Boolean);
  }

  // entity id → { field → [values] }, honouring specimen merge + resolutions.
  function _valuesByEntity(fields) {
    const by = {};
    _records().forEach((r) => {
      const sn = r && r.sample_name;
      if (!sn) return;
      if (typeof sampleHidden !== "undefined" && sampleHidden[sn]) return;
      const id = _viewId(sn);
      if (!by[id]) by[id] = {};
      fields.forEach((f) => {
        if (!by[id][f]) by[id][f] = new Set();
        _cellValues(r, f).forEach((v) => by[id][f].add(v));
      });
    });
    Object.keys(by).forEach((id) => {
      fields.forEach((f) => {
        const ov = _resolvedOverride(id, f);
        if (ov && ov.length) by[id][f] = new Set(ov);
      });
    });
    return by;
  }

  // Cartesian product of the per-field value sets → composite keys.
  // A sample missing a value for one of the fields contributes "(unset)" for
  // that dimension rather than dropping out entirely, so a partially-annotated
  // run still shows every sample somewhere.
  const UNSET = "(unset)";

  function _combosFor(perField, fields) {
    let combos = [[]];
    for (let i = 0; i < fields.length; i++) {
      const vals = perField[fields[i]] && perField[fields[i]].size ? Array.from(perField[fields[i]]) : [UNSET];
      const next = [];
      combos.forEach((c) => vals.forEach((v) => next.push(c.concat([v]))));
      combos = next;
      // Runaway guard: a pathological multi-valued combination could explode.
      if (combos.length > 400) break;
    }
    return combos;
  }

  /* ── The grouping model ────────────────────────────────────────────────
     {
       fields:    [field, …],
       order:     [key, …]    every group, display order, largest first
       visible:   [key, …]    `order` minus the hidden groups
       highlight: [key, …]    the emphasised groups ([] when none)
       groups:    { key: {key, parts, label, samples:[id], color, shape,
                          state, hidden, highlighted, dimmed} }
       bySample:  { id: [key, …] }
     }
     `dimmed` is the derived flag the renderers actually want: true when
     something else is highlighted and this group is not. Computing it once
     here keeps every view from re-deriving (and disagreeing about) it.
     `order` deliberately keeps the hidden groups: the legend has to keep
     listing them (greyed out) or there would be no way to switch one back on.
     Views iterate `visible`; legends iterate `order`.
     Cached until the fields, the merge state, or the hidden-sample set change. */
  function _cacheToken(fields) {
    const merge = typeof specimenMergeEnabled !== "undefined" && specimenMergeEnabled ? "1" : "0";
    let hidden = 0;
    try {
      hidden = Object.keys(sampleHidden || {}).filter((k) => sampleHidden[k]).length;
    } catch (e) {}
    // The toggled-off groups are part of the token: without them a toggle
    // would leave `visible` stale until something else invalidated the cache.
    return fields.join(",") + "|" + merge + "|" + hidden + "|" + _records().length + "|" + _hiddenSeq;
  }

  function _buildModel() {
    const fields = _fields.slice();
    const model = { fields, order: [], visible: [], highlight: [], groups: {}, bySample: {} };
    if (!fields.length) return model;

    const byEntity = _valuesByEntity(fields);
    Object.keys(byEntity).forEach((id) => {
      const combos = _combosFor(byEntity[id], fields);
      combos.forEach((parts) => {
        const key = parts.join(SEP);
        if (!model.groups[key]) model.groups[key] = { key, parts, label: key, samples: [] };
        if (model.groups[key].samples.indexOf(id) === -1) model.groups[key].samples.push(id);
        if (!model.bySample[id]) model.bySample[id] = [];
        if (model.bySample[id].indexOf(key) === -1) model.bySample[id].push(key);
      });
    });

    model.order = Object.keys(model.groups).sort(
      (a, b) => model.groups[b].samples.length - model.groups[a].samples.length || a.localeCompare(b),
    );

    // Stable colour/shape: assigned on first sight and never reshuffled, so a
    // group keeps its colour when the user narrows a filter or toggles merge.
    const pal = _palette();
    model.order.forEach((key) => {
      if (!_colorBy[key]) {
        const used = Object.keys(_colorBy).length;
        _colorBy[key] = pal[used % pal.length];
        _shapeBy[key] = SHAPES[Math.floor(used / pal.length) % SHAPES.length];
      }
      model.groups[key].color = _colorBy[key];
      model.groups[key].shape = _shapeBy[key] || SHAPES[0];
    });

    model.visible = model.order.filter((k) => _stateOf(k) !== STATE_HIDDEN);
    // A highlight on a group that no longer exists must not dim the whole map,
    // so only live groups count towards "is anything highlighted".
    model.highlight = model.visible.filter((k) => _stateOf(k) === STATE_HIGHLIGHT);
    const anyHi = model.highlight.length > 0;
    model.order.forEach((key) => {
      const st = _stateOf(key);
      const g = model.groups[key];
      g.state = st;
      g.hidden = st === STATE_HIDDEN;
      g.highlighted = st === STATE_HIGHLIGHT;
      g.dimmed = anyHi && !g.hidden && !g.highlighted;
    });
    return model;
  }

  function _model() {
    const token = _cacheToken(_fields);
    if (!_modelCache || _modelKey !== token) {
      _modelCache = _buildModel();
      _modelKey = token;
    }
    return _modelCache;
  }

  function _invalidate() {
    _modelCache = null;
    _modelKey = "";
  }

  // Common tail for every visibility mutation: invalidate, then tell the views.
  function _bumpVisibility(opts) {
    _hiddenSeq += 1;
    _invalidate();
    if (!(opts && opts.silent)) _notify();
  }

  function _notify() {
    _listeners.forEach((fn) => {
      try {
        fn(_fields.slice());
      } catch (e) {
        /* a broken subscriber must not break the others */
      }
    });
  }

  /* ── Public API ────────────────────────────────────────────────────────── */
  const api = {
    SEP,
    UNSET,
    SHAPES,

    fields: function () {
      return _fields.slice();
    },

    // Set the active grouping. Unknown / duplicate fields are dropped so a
    // restored session that predates a metadata change can't break the views.
    setFields: function (fields, opts) {
      opts = opts || {};
      const valid = new Set(_mgColumnStats().map((c) => c.key));
      const next = [];
      (Array.isArray(fields) ? fields : []).forEach((f) => {
        if (valid.has(f) && next.indexOf(f) === -1) next.push(f);
      });
      const changed = next.join(",") !== _fields.join(",");
      _fields = next;
      // Changing the grouping changes the key space entirely ("Water" becomes
      // "Water │ MiSeq"), so carrying the old states over would either do
      // nothing or silently re-hide something the user forgot about. A new
      // grouping starts with every group normal.
      if (changed && Object.keys(_groupState).length) {
        _groupState = {};
        _hiddenSeq += 1;
      }
      _invalidate();
      if (changed && !opts.silent) _notify();
      return changed;
    },

    active: function () {
      return _fields.length > 0;
    },

    /* ── Group display state ────────────────────────────────────────────
       normal → highlight → hidden → normal, cycled by clicking a legend chip.
       Highlighting emphasises a group and fades the rest (and draws the map's
       similarity network for it); hiding drops it from the views and the map
       entirely. Non-normal groups stay in `model().order` so the legend can
       always list them and cycle them back. */
    STATE_NORMAL,
    STATE_HIGHLIGHT,
    STATE_HIDDEN,

    stateOf: _stateOf,

    isHidden: function (key) {
      return _stateOf(key) === STATE_HIDDEN;
    },

    isHighlighted: function (key) {
      return _stateOf(key) === STATE_HIGHLIGHT;
    },

    // True when this group should be drawn faded: something else is
    // highlighted and this is not it.
    isDimmed: function (key) {
      const m = _model();
      return !!(m.groups[key] && m.groups[key].dimmed);
    },

    anyHighlighted: function () {
      return _model().highlight.length > 0;
    },

    highlighted: function () {
      return _model().highlight.slice();
    },

    hidden: _hiddenKeys,

    // Only count groups that still exist, or a stale key from a previous
    // grouping would show a "3 hidden" badge with nothing to restore.
    hiddenCount: function () {
      const live = _model().groups;
      return _hiddenKeys().filter((k) => live[k]).length;
    },

    // Advance one step around the cycle. `back` steps the other way, so a
    // shift-click can undo an over-shoot instead of going all the way round.
    cycle: function (key, opts) {
      if (!key) return _stateOf(key);
      const i = CYCLE.indexOf(_stateOf(key));
      const step = opts && opts.back ? -1 : 1;
      const next = CYCLE[(i + step + CYCLE.length) % CYCLE.length];
      return api.setState(key, next, opts);
    },

    setState: function (key, state, opts) {
      if (!key) return STATE_NORMAL;
      const st = CYCLE.indexOf(state) === -1 ? STATE_NORMAL : state;
      if (st === STATE_NORMAL) delete _groupState[key];
      else _groupState[key] = st;
      _bumpVisibility(opts);
      return st;
    },

    setGroupStates: function (map, opts) {
      _groupState = {};
      Object.keys(map || {}).forEach((k) => {
        if (CYCLE.indexOf(map[k]) > 0) _groupState[k] = map[k];
      });
      _bumpVisibility(opts);
    },

    // Back-compat shim for the plain 2-state call sites: hide ⇄ show.
    setHiddenGroups: function (keys, opts) {
      _groupState = {};
      (Array.isArray(keys) ? keys : []).forEach((k) => (_groupState[k] = STATE_HIDDEN));
      _bumpVisibility(opts);
    },

    // Show this group and nothing else. Re-soloing the group that is already
    // alone restores everything, so the same call both enters and leaves the
    // isolated view.
    solo: function (key, opts) {
      const m = _model();
      const alreadySolo = m.visible.length === 1 && m.visible[0] === key;
      _groupState = {};
      if (!alreadySolo) m.order.forEach((k) => k !== key && (_groupState[k] = STATE_HIDDEN));
      _bumpVisibility(opts);
      return !alreadySolo;
    },

    // Reset every group to normal — clears highlights AND un-hides.
    showAllGroups: function (opts) {
      if (!Object.keys(_groupState).length) return false;
      _groupState = {};
      _bumpVisibility(opts);
      return true;
    },

    // Clear only the highlights, leaving hidden groups hidden.
    clearHighlight: function (opts) {
      const hi = _highlightKeys();
      if (!hi.length) return false;
      hi.forEach((k) => delete _groupState[k]);
      _bumpVisibility(opts);
      return true;
    },

    // True when a sample still belongs to at least one visible group. A
    // multi-valued sample stays visible while ANY of its groups is on —
    // hiding "Soil" must not remove a sample that is also "Water".
    sampleVisible: function (sample) {
      if (!_fields.length) return true;
      const keys = _model().bySample[_viewId(sample)];
      if (!keys || !keys.length) return true; // ungrouped samples are never hidden by a group toggle
      return keys.some((k) => _stateOf(k) !== STATE_HIDDEN);
    },

    // True when a sample belongs to a highlighted group. Used by the map to
    // decide which markers stay vivid.
    sampleHighlighted: function (sample) {
      if (!_fields.length) return false;
      const keys = _model().bySample[_viewId(sample)];
      if (!keys || !keys.length) return false;
      return keys.some((k) => _stateOf(k) === STATE_HIGHLIGHT);
    },

    // True when a sample should be faded: a highlight is active and this
    // sample is in none of the highlighted groups.
    sampleDimmed: function (sample) {
      if (!api.anyHighlighted()) return false;
      return !api.sampleHighlighted(sample);
    },

    candidates: _mgColumnStats,

    // Re-profile the columns and drop any active field that no longer exists
    // (e.g. after a meta.csv upload replaced the metadata).
    refresh: function () {
      _invalidate();
      const valid = new Set(_mgColumnStats().map((c) => c.key));
      const kept = _fields.filter((f) => valid.has(f));
      if (kept.length !== _fields.length) {
        _fields = kept;
        _invalidate();
      }
      return _fields.slice();
    },

    model: _model,

    // Primary (first) group for a sample — what a single-colour marker uses.
    groupOf: function (sample) {
      const g = _model().bySample[_viewId(sample)];
      return g && g.length ? g[0] : null;
    },

    // Every group a sample belongs to (>1 only for multi-valued cells).
    groupsOf: function (sample) {
      const g = _model().bySample[_viewId(sample)];
      return g ? g.slice() : [];
    },

    samplesOf: function (key) {
      const g = _model().groups[key];
      return g ? g.samples.slice() : [];
    },

    color: function (key) {
      if (_colorBy[key]) return _colorBy[key];
      const m = _model();
      return (m.groups[key] && m.groups[key].color) || "#8899aa";
    },

    shape: function (key) {
      if (_shapeBy[key]) return _shapeBy[key];
      const m = _model();
      return (m.groups[key] && m.groups[key].shape) || "circle";
    },

    // Human-readable description of the current grouping, for chart subtitles.
    describe: function () {
      if (!_fields.length) return "";
      return _fields.map(_label).join(" × ");
    },

    onChange: function (fn) {
      if (typeof fn === "function") _listeners.push(fn);
    },

    invalidate: _invalidate,

    // Serialize / restore with the rest of the session state.
    capture: function () {
      return {
        fields: _fields.slice(),
        colors: Object.assign({}, _colorBy),
        shapes: Object.assign({}, _shapeBy),
        groupState: Object.assign({}, _groupState),
      };
    },

    restore: function (state) {
      state = state && typeof state === "object" ? state : {};
      _colorBy = Object.assign({}, state.colors || {});
      _shapeBy = Object.assign({}, state.shapes || {});
      // Fields FIRST: setFields() clears the toggles whenever the grouping
      // changes (the key space changes with it), so restoring the hidden set
      // before this would immediately wipe it.
      api.setFields(state.fields || [], { silent: true });
      _groupState = {};
      // Accept both the current {key: state} map and the older hiddenGroups
      // array, so a session saved before the 3-state cycle still restores.
      const gs = state.groupState;
      if (gs && typeof gs === "object") {
        Object.keys(gs).forEach((k) => {
          if (CYCLE.indexOf(gs[k]) > 0) _groupState[k] = gs[k];
        });
      } else if (Array.isArray(state.hiddenGroups)) {
        state.hiddenGroups.forEach((k) => (_groupState[k] = STATE_HIDDEN));
      }
      _hiddenSeq += 1;
      _invalidate();
    },
  };

  window.metaGrouping = api;

  /* ── SVG shape path helpers ────────────────────────────────────────────
     Shared by the map markers, the network nodes, and the legends so a group
     looks identical everywhere it appears. Returns an SVG element string
     centred on (cx, cy) with a bounding radius r.                            */
  function _shapePath(shape, cx, cy, r) {
    const pts = (n, rot) => {
      const a = [];
      for (let i = 0; i < n; i++) {
        const t = rot + (i * 2 * Math.PI) / n;
        a.push(`${(cx + r * Math.cos(t)).toFixed(2)},${(cy + r * Math.sin(t)).toFixed(2)}`);
      }
      return a.join(" ");
    };
    switch (shape) {
      case "square":
        return `<rect x="${cx - r * 0.88}" y="${cy - r * 0.88}" width="${r * 1.76}" height="${r * 1.76}"`;
      case "triangle":
        return `<polygon points="${pts(3, -Math.PI / 2)}"`;
      case "diamond":
        return `<polygon points="${pts(4, -Math.PI / 2)}"`;
      case "hexagon":
        return `<polygon points="${pts(6, -Math.PI / 2)}"`;
      case "star": {
        const a = [];
        for (let i = 0; i < 10; i++) {
          const rr = i % 2 ? r * 0.45 : r;
          const t = -Math.PI / 2 + (i * Math.PI) / 5;
          a.push(`${(cx + rr * Math.cos(t)).toFixed(2)},${(cy + rr * Math.sin(t)).toFixed(2)}`);
        }
        return `<polygon points="${a.join(" ")}"`;
      }
      case "cross": {
        const w = r * 0.4;
        return (
          `<polygon points="${cx - w},${cy - r} ${cx + w},${cy - r} ${cx + w},${cy - w} ${cx + r},${cy - w} ` +
          `${cx + r},${cy + w} ${cx + w},${cy + w} ${cx + w},${cy + r} ${cx - w},${cy + r} ` +
          `${cx - w},${cy + w} ${cx - r},${cy + w} ${cx - r},${cy - w} ${cx - w},${cy - w}"`
        );
      }
      case "wedge":
        return `<polygon points="${cx},${cy + r} ${cx + r},${cy - r * 0.6} ${cx - r},${cy - r * 0.6}"`;
      default:
        return `<circle cx="${cx}" cy="${cy}" r="${r}"`;
    }
  }

  // Full <svg> glyph for a group — used in legends and picker chips.
  window._mgShapeSvg = function (shape, color, size, opts) {
    opts = opts || {};
    const s = size || 14;
    const r = (s / 2) * 0.86;
    const stroke = opts.stroke || "#37474f";
    const sw = opts.strokeWidth == null ? 1.1 : opts.strokeWidth;
    return (
      `<svg width="${s}" height="${s}" viewBox="0 0 ${s} ${s}" style="vertical-align:middle;flex:none">` +
      `${_shapePath(shape, s / 2, s / 2, r)} fill="${color}" stroke="${stroke}" stroke-width="${sw}"` +
      (opts.opacity != null ? ` opacity="${opts.opacity}"` : "") +
      ` /></svg>`
    );
  };
  window._mgShapeEl = _shapePath;

  /* ── Group picker ──────────────────────────────────────────────────────
     Rendered into any container by id, as a COLLAPSED SEARCHABLE DROPDOWN
     (see 43_multiselect.js). It used to be a row of chips, one per group —
     fine for a handful of groups, unusable for the fifty-plus a real
     multi-site run produces: the row wrapped into a wall that pushed the map
     and the heatmap off screen, and finding one group meant reading all of
     them. The dropdown is one line high whatever the run size and has a
     search box.

     Every group starts CHECKED (visible). Three states survive the move:

       checkbox off  → HIDDEN     removed from the views and the map
       row click     → HIGHLIGHT  this group is emphasised, everything else
                                  fades, and the map draws its similarity
                                  network; clicking again steps on through the
                                  normal → highlight → hidden cycle
       shift-click   → step BACKWARDS through the cycle (undo an over-shoot)
       alt-click     → solo (show only this group)

     A "Reset (n)" button appears beside the dropdown once anything is
     off-default. Registered pickers re-render together, so a change in the
     Group Heatmap is immediately reflected in the bar, the map and the
     network.

     opts.onPick — optional extra callback fired after a change (the map uses
     it to refit its bounds and redraw the overlay).                          */
  const _STATE_ICON = { normal: "fa-circle", highlight: "fa-bullseye", hidden: "fa-eye-slash" };
  const _NEXT_LABEL = { normal: "highlight it", highlight: "hide it", hidden: "restore it" };
  const _legendTargets = new Map(); // id → opts, for cross-view re-render

  window._mgRenderLegend = function (elOrId, opts) {
    opts = opts || {};
    const el = typeof elOrId === "string" ? document.getElementById(elOrId) : elOrId;
    if (!el) return;
    if (typeof elOrId === "string") _legendTargets.set(elOrId, opts);

    const m = api.model();
    if (!m.fields.length || !m.order.length) {
      el.innerHTML = "";
      el.style.display = "none";
      return;
    }
    el.style.display = "";
    const nOff = m.order.filter((k) => m.groups[k].state !== STATE_NORMAL).length;
    const nShown = m.order.filter((k) => !m.groups[k].hidden).length;

    const items = m.order.map((k) => {
      const g = m.groups[k];
      const n = g.samples.length;
      return {
        key: k,
        label: k,
        count: n,
        //  Unchecked == hidden, so the box reads as "is this group in the
        //  views?" rather than as a separate selection concept.
        checked: !g.hidden,
        stateIcon: _STATE_ICON[g.state],
        swatchHTML: window._mgShapeSvg(
          g.shape,
          g.hidden ? "#cfd8dc" : g.color,
          13,
          g.hidden ? { stroke: "#b0bec5" } : null,
        ),
        rowClass: (g.hidden ? "mg-off " : "") + (g.highlighted ? "mg-hi " : "") + (g.dimmed ? "mg-dim" : ""),
        title:
          `${k} — ${n} sample${n === 1 ? "" : "s"}\n` +
          `Currently: ${g.state}.\n` +
          `Untick to hide it · click the row to ${_NEXT_LABEL[g.state]}\n` +
          `Shift-click steps back · Alt-click shows only this group`,
      };
    });

    window.ttMultiSelect.render(el, {
      icon: "fa-layer-group",
      title: api.describe() + " — choose which groups the views show",
      placeholder: "Search groups…",
      emptyText: "No groups.",
      summary: (n, total) => `${_esc(api.describe())}: ${nShown} of ${total} group${total === 1 ? "" : "s"} shown`,
      items,
      onToggle: (key, checked) => {
        api.setState(key, checked ? STATE_NORMAL : STATE_HIDDEN);
        if (typeof opts.onPick === "function") opts.onPick(key);
      },
      onRowClick: (key, ev) => {
        if (ev.altKey) api.solo(key);
        else api.cycle(key, { back: ev.shiftKey });
        if (typeof opts.onPick === "function") opts.onPick(key);
      },
      onAll: () => {
        api.showAllGroups();
        if (typeof opts.onPick === "function") opts.onPick(null);
      },
      onNone: () => {
        //  One call, one broadcast — looping setState would re-render every
        //  registered picker once per group.
        api.setHiddenGroups(m.order);
        if (typeof opts.onPick === "function") opts.onPick(null);
      },
      footNote:
        '<i class="fas fa-circle-info"></i> Untick to hide a group · click a row to cycle ' +
        "normal \u2192 highlight \u2192 hidden",
      actions: nOff
        ? [
            {
              label: '<i class="fas fa-rotate-left"></i> Reset (' + nOff + ")",
              title: "Reset every group to normal",
              className: "mg-legend-showall",
              onClick: () => api.showAllGroups(),
            },
          ]
        : [],
    });
  };

  // Re-render every legend that has been rendered at least once. Called from
  // the change broadcast so one toggle updates all of them.
  function _mgRefreshLegends() {
    _legendTargets.forEach((o, id) => window._mgRenderLegend(id, o));
  }
  window._mgRefreshLegends = _mgRefreshLegends;

  function _esc(s) {
    return String(s == null ? "" : s).replace(
      /[&<>"']/g,
      (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c],
    );
  }
  window._mgEsc = _esc;

  /* ── Shared "Group by" bar (multi-instance) ────────────────────────────
     A multi-select popover over the detected candidate columns. Selecting
     more than one column means UNION (composite key), which is what lets a
     user combine, say, environmental_site with sequencing_platform.

     The bar is MOUNTED, not authored: every `<div class="mg-bar-mount"
     data-mg-mount="…">` in the document gets its own copy of the same
     control, and all copies read and write the ONE selection held in
     `_fields`. That is what lets the Metadata tab own the grouping while the
     Mapping and Trends tabs carry a live mirror of it — change it in any of
     them and all three update together.

     Instance ids: "meta" (the primary, in the Metadata tab), "map", "trends".
     Mirrors also get a "Manage columns in Metadata →" link so the table that
     defines the columns is one click away from wherever you're looking.      */

  function _mgBarHtml(inst) {
    const isPrimary = inst === "meta";
    const link = isPrimary
      ? '<span class="mg-shared-note" title="This selection also drives the Mapping and Trends tabs">' +
        '<i class="fas fa-link"></i> shared with Mapping &amp; Trends</span>'
      : '<a href="#" class="mg-goto-meta mg-goto-meta-link" ' +
        'title="Open the Metadata tab to add, edit or import the columns you group on">' +
        '<i class="fas fa-table-list"></i> Manage columns in Metadata</a>';
    return (
      '<div class="mg-bar" id="mg-bar-' +
      inst +
      '" data-mg-inst="' +
      inst +
      '">' +
      '<div class="mg-bar-row">' +
      '<span class="mg-bar-label"><i class="fas fa-layer-group"></i> Group by</span>' +
      '<div class="mg-field-wrap">' +
      '<button class="mg-field-btn" type="button" title="Choose one or more metadata columns to group on">' +
      '<span>Select columns</span><i class="fas fa-caret-down"></i></button>' +
      '<div class="mg-field-pop">' +
      '<div class="mg-field-pop-head">Pick one or more columns. Multiple columns combine as a ' +
      "<strong>union</strong> &mdash; one group per observed combination.</div>" +
      '<div class="mg-field-list"></div>' +
      "</div></div>" +
      '<div class="mg-summary"></div>' +
      '<button class="mg-clear" type="button" style="display:none" title="Clear the grouping">' +
      '<i class="fas fa-xmark"></i> Clear</button>' +
      '<span class="mg-bar-spacer"></span>' +
      link +
      "</div>" +
      '<div class="mg-legend" id="mg-legend-' +
      inst +
      '"></div>' +
      "</div>"
    );
  }

  /* Fill every .mg-bar-mount that has not been filled yet, then wire it.
     Safe to call repeatedly (mount points remember they are done).           */
  window._mgMountGroupBars = function () {
    document.querySelectorAll(".mg-bar-mount").forEach((mount) => {
      if (mount.getAttribute("data-mg-mounted")) return;
      const inst = mount.getAttribute("data-mg-mount") || "meta";
      mount.innerHTML = _mgBarHtml(inst);
      mount.setAttribute("data-mg-mounted", "1");
    });
    window._mgWireGroupBar();
  };

  function _mgBars() {
    return Array.from(document.querySelectorAll(".mg-bar"));
  }

  function _renderGroupBar() {
    const lists = document.querySelectorAll(".mg-bar .mg-field-list");
    if (!lists.length) return;
    const cands = _mgColumnStats();
    const active = new Set(_fields);

    lists.forEach((list) => {
      if (!cands.length) {
        list.innerHTML =
          '<div class="mg-empty">No groupable metadata columns found. Add columns in the Metadata tab, ' +
          "or upload a metadata CSV with categorical fields.</div>";
        return;
      }
      list.innerHTML = cands
        .map((c) => {
          const pct = Math.round(c.coverage * 100);
          return (
            `<label class="mg-field-opt${c.preferred ? " mg-preferred" : ""}">` +
            `<input type="checkbox" data-field="${_esc(c.key)}"${active.has(c.key) ? " checked" : ""} />` +
            `<span class="mg-field-name">${_esc(c.label)}</span>` +
            `<span class="mg-field-meta">${c.distinct} value${c.distinct === 1 ? "" : "s"} · ${pct}%</span>` +
            `</label>`
          );
        })
        .join("");
      list.querySelectorAll("input[data-field]").forEach((cb) => {
        cb.addEventListener("change", () => {
          const picked = Array.from(list.querySelectorAll("input[data-field]:checked")).map((x) =>
            x.getAttribute("data-field"),
          );
          // Preserve the user's click order so composite labels read the way
          // they built them ("Site │ Platform", not alphabetical).
          const ordered = _fields.filter((f) => picked.indexOf(f) !== -1);
          picked.forEach((f) => {
            if (ordered.indexOf(f) === -1) ordered.push(f);
          });
          api.setFields(ordered);
        });
      });
    });
    _renderGroupSummary();
  }

  function _renderGroupSummary() {
    const html = !_fields.length
      ? '<span class="mg-none">No grouping — views show individual samples</span>'
      : (() => {
          const m = api.model();
          return (
            `<span class="mg-chips">` +
            _fields
              .map((f) => `<span class="mg-chip">${_esc(_label(f))}</span>`)
              .join('<span class="mg-plus">+</span>') +
            `</span><span class="mg-count">${m.order.length} group${m.order.length === 1 ? "" : "s"}</span>`
          );
        })();
    document.querySelectorAll(".mg-bar .mg-summary").forEach((sum) => {
      sum.innerHTML = html;
    });
    document.querySelectorAll(".mg-bar .mg-clear").forEach((clr) => {
      clr.style.display = _fields.length ? "inline-flex" : "none";
    });
    // One legend per mounted bar, each registered under its own id so the
    // shared refresh loop keeps every copy in step.
    _mgBars().forEach((bar) => {
      const lg = bar.querySelector(".mg-legend[id]");
      if (lg) _mgRenderLegend(lg.id);
    });
  }

  window._mgSyncGroupBar = function () {
    api.refresh();
    _renderGroupBar();
  };

  /* Jump to the Metadata tab and flash the table — used by the mirrors' link
     and by the "add a column" prompts scattered through the metadata views. */
  window._mgGotoMetaTable = function () {
    const btn = document.getElementById("runmeta-tab-btn");
    if (btn) btn.click();
    const wrap = document.getElementById("runmeta-table-wrap");
    const bar = document.getElementById("runmeta-toolbar");
    setTimeout(() => {
      if (wrap && wrap.scrollIntoView) wrap.scrollIntoView({ behavior: "smooth", block: "center" });
      [bar, wrap].forEach((el) => {
        if (!el) return;
        el.classList.remove("mg-flash");
        // reflow so the animation restarts on a repeat click
        void el.offsetWidth;
        el.classList.add("mg-flash");
        setTimeout(() => el.classList.remove("mg-flash"), 1600);
      });
    }, 80);
  };

  // Any link anywhere in the report can carry .mg-goto-meta-link to reach the
  // metadata table — delegated so markup added later still works.
  document.addEventListener("click", (e) => {
    const a = e.target.closest && e.target.closest(".mg-goto-meta-link, .mg-goto-meta");
    if (!a) return;
    e.preventDefault();
    window._mgGotoMetaTable();
  });

  window._mgWireGroupBar = function () {
    _mgBars().forEach((bar) => {
      if (bar.getAttribute("data-mg-wired")) return;
      const btn = bar.querySelector(".mg-field-btn");
      const pop = bar.querySelector(".mg-field-pop");
      if (!btn || !pop) return;
      bar.setAttribute("data-mg-wired", "1");

      btn.addEventListener("click", (e) => {
        e.stopPropagation();
        const open = pop.style.display === "block";
        // Only one popover open at a time, whichever bar it belongs to.
        document.querySelectorAll(".mg-bar .mg-field-pop").forEach((p) => (p.style.display = "none"));
        pop.style.display = open ? "none" : "block";
        if (!open) _renderGroupBar();
      });
      document.addEventListener("click", (e) => {
        if (pop.style.display === "block" && !pop.contains(e.target) && e.target !== btn) pop.style.display = "none";
      });
      const clr = bar.querySelector(".mg-clear");
      if (clr) clr.addEventListener("click", () => api.setFields([]));
    });

    if (!_barWired && _mgBars().length) {
      _barWired = true;
      api.onChange(() => {
        _renderGroupSummary();
        _mgBroadcastRedraw();
      });
    }
    if (_mgBars().length) _renderGroupBar();
  };

  let _barWired = false;

  // Mount as soon as the DOM is there — the panes are static markup, so this
  // runs once and every later _mgSyncGroupBar() call just refreshes content.
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", () => window._mgMountGroupBars());
  } else {
    window._mgMountGroupBars();
  }

  // Ask every grouping-aware view to re-render. Each guard is defensive: the
  // views live in later files and may not exist in a stripped build.
  function _mgBroadcastRedraw() {
    try {
      _mgRefreshLegends();
    } catch (e) {}
    try {
      if (typeof _refreshMapMarkerColors === "function") _refreshMapMarkerColors();
    } catch (e) {}
    // The map's similarity overlay + its table only exist while something is
    // highlighted, so they are rebuilt (or torn down) on every state change.
    try {
      if (typeof _mgRefreshMapNetwork === "function") _mgRefreshMapNetwork();
    } catch (e) {}
    try {
      if (typeof _mgBuildGroupHeatmap === "function" && _activeMetaSub === "ghm") _mgBuildGroupHeatmap();
    } catch (e) {}
    try {
      if (typeof _mgBuildNetwork === "function" && _activeMetaSub === "net") _mgBuildNetwork();
    } catch (e) {}
    try {
      if (typeof _geoRedraw === "function") _geoRedraw();
    } catch (e) {}
  }
  window._mgBroadcastRedraw = _mgBroadcastRedraw;
})();
