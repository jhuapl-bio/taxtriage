/* ═══════════════════════════════════════════════════════════════════════════
       -  §  CONTRAST TAB
       -     Group samples into 2–4 sets and compare what was detected in each.
       -     A set is built from any categorical metadata column (sample type,
       -     site, run, host, country, … — the same columns the Group by picker
       -     offers), the shared Group by selection, the collection date
       -     (date, month, season, year or a date range) or hand-picked samples.
       -
       -     Views (sub-tabs): Mirrored bars · Change (Δ, difference or log₂
       -     fold) · Radar · A-vs-B Scatter · Table (with Mann–Whitney U p and
       -     Benjamini–Hochberg q). Clicking an organism opens a per-sample
       -     strip plot of it under the chart.
       -
       -     Reads filteredData(), so every sidebar filter applies: hidden
       -     samples, TASS cutoff, search, categories, the Strain / Species /
       -     Genus view level and specimen merge. A sample with no detections
       -     after filtering still counts in its set's denominator.
       -
       -     Public: window.ttContrast = { render, state, compute, setSets }.
       -     Drawn by _drawTab("contrast") (21_compare_coverages.js).
       -     Settings persist per report file in localStorage.
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  "use strict";
  var COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"];
  var LETTERS = ["A", "B", "C", "D"];
  var TIME_FIELDS = { day: "Collection date", month: "Month", season: "Season", year: "Year", range: "Date range" };
  var METRICS = {
    rel: { label: "Mean % of sample reads", short: "Mean % of reads", kind: "pct" },
    prev: { label: "Prevalence (% of samples detected)", short: "Prevalence", kind: "pct" },
    mreads: { label: "Mean reads aligned per sample", short: "Mean reads", kind: "count" },
    rpm: { label: "Mean RPM (reads per million)", short: "Mean RPM", kind: "count" },
    tass: { label: "Mean TASS score (0 where not detected)", short: "Mean TASS", kind: "score" },
  };
  var MK = ["rel", "mreads", "rpm", "tass"]; // per-sample value arrays kept for each of these
  var GROUP_BY = [
    ["Microbial Category", "Microbial category"],
    ["Kingdom", "Kingdom"],
    ["Phylum", "Phylum"],
    ["Class", "Class"],
    ["Order", "Order"],
    ["Family", "Family"],
    ["Genus", "Genus"],
  ];
  var VIEWS = [
    ["mirror", "Mirrored bars"],
    ["delta", "Change (Δ)"],
    ["radar", "Radar"],
    ["scatter", "Scatter"],
    ["table", "Table"],
  ];
  var SEASONS = ["Winter", "Spring", "Summer", "Autumn"];
  // Metadata columns read as the collection date (first one with a parseable value wins).
  var DATE_FIELDS = [
    "collection_time",
    "collection_date",
    "collection date",
    "date",
    "sample_date",
    "sampling_date",
    "eventdate",
  ];
  var S = {
    sets: null,
    level: "org",
    gby: "Microbial Category",
    metric: "rel",
    view: "mirror",
    topN: 20,
    tgroup: "",
    exclCtl: true,
    left: 0,
    right: 1,
    scale: "lin",
    sort: "delta",
    dmode: "diff",
    raxes: "grp",
    rk: 8,
    tsort: "absd",
    tdir: -1,
    tall: false,
    pfield: "",
    sscale: "log",
    ssort: "sum",
    sel: "",
    q: {},
  };
  var HOST = null,
    D = null,
    C = null;
  var LS = "tt.contrast.v1:" + ((location.pathname || "").split("/").pop() || "report");

  function d3() {
    return window.d3;
  }
  function esc(s) {
    return String(s == null ? "" : s).replace(/[&<>"]/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c];
    });
  }
  function natural(a, b) {
    return String(a).localeCompare(String(b), undefined, { numeric: true, sensitivity: "base" });
  }
  function fmtN(n) {
    return Number(n || 0).toLocaleString();
  }
  function n0(v) {
    var x = Number(v);
    return isFinite(x) ? x : 0;
  }
  function tip(e, h) {
    if (typeof showTip === "function") showTip(h, e);
  }
  function untip() {
    if (typeof hideTip === "function") hideTip();
  }
  function pad(n) {
    return (n < 10 ? "0" : "") + n;
  }
  function sum(a) {
    var t = 0;
    for (var i = 0; i < a.length; i++) t += a[i];
    return t;
  }
  function isPlaceholder(k) {
    return k.charAt(0) === "(";
  }

  function load() {
    try {
      var o = JSON.parse(localStorage.getItem(LS) || "null");
      if (o)
        Object.keys(S).forEach(function (k) {
          if (k !== "sel" && k !== "q" && o[k] !== undefined) S[k] = o[k];
        });
    } catch (e) {
      /* storage blocked */
    }
  }
  function save() {
    try {
      var o = {};
      Object.keys(S).forEach(function (k) {
        if (k !== "sel" && k !== "q") o[k] = S[k];
      });
      localStorage.setItem(LS, JSON.stringify(o));
    } catch (e) {
      /* storage blocked */
    }
  }
  load();

  /* ── data ─────────────────────────────────────────────────────────────── */
  function cellValues(rec, field) {
    var v = rec ? rec[field] : null;
    if (v == null) return [];
    var arr = Array.isArray(v) ? v : [v],
      out = [];
    arr.forEach(function (x) {
      if (x == null || typeof x === "object") return;
      var s = String(x).trim(),
        l = s.toLowerCase();
      if (s && l !== "null" && l !== "na" && l !== "nan" && l !== "none" && l !== "n/a") out.push(s);
    });
    return out;
  }
  function parseDate(v) {
    if (v == null || v === "") return null;
    var d = typeof _parseLongiDate === "function" ? _parseLongiDate(v) : new Date(String(v));
    return d && !isNaN(d) ? d.getTime() : null;
  }
  function dateOf(rec) {
    for (var i = 0; i < DATE_FIELDS.length; i++) {
      var vs = cellValues(rec, DATE_FIELDS[i]);
      for (var j = 0; j < vs.length; j++) {
        var t = parseDate(vs[j]);
        if (t != null) return t;
      }
    }
    return null;
  }
  function grouping() {
    var G = window.metaGrouping;
    try {
      return G && G.active && G.active() ? G : null;
    } catch (e) {
      return null;
    }
  }
  // Metadata columns offered as set fields: the Group by picker's candidates
  // (curated first), minus the date columns (they get the time fields instead).
  function metaCandidates() {
    var c = [];
    try {
      if (window.metaGrouping && window.metaGrouping.candidates) c = window.metaGrouping.candidates() || [];
    } catch (e) {
      c = [];
    }
    var out = c
      .filter(function (x) {
        return DATE_FIELDS.indexOf(String(x.key).toLowerCase()) < 0;
      })
      .map(function (x) {
        return { key: x.key, label: x.label || String(x.key).replace(/_/g, " ") };
      });
    if (
      !out.some(function (x) {
        return x.key === "sample_type";
      })
    )
      out.unshift({ key: "sample_type", label: "Sample type" });
    return out;
  }
  function attrOf(r, k) {
    var v;
    if (k === "Microbial Category") {
      v = String(r["Microbial Category"] || "")
        .split(";")[0]
        .trim();
      return v || (typeof UNCATEGORIZED !== "undefined" ? UNCATEGORIZED : "Uncategorized");
    }
    if (k === "Kingdom") v = r.Kingdom || r.Domain || r.Superkingdom;
    else if (k === "Genus") v = r.Genus || r["Genus Name"];
    else v = r[k];
    v = v == null ? "" : String(v).trim();
    return v;
  }
  function prepare() {
    var ids = typeof _allSampleIds === "function" ? _allSampleIds() : Object.keys(SAMPLE_META || {});
    var SM = (typeof SAMPLE_META !== "undefined" && SAMPLE_META) || {};
    var RM = {};
    ((typeof RUN_META !== "undefined" && RUN_META) || []).forEach(function (r) {
      if (r && r.sample_name != null) RM[r.sample_name] = r;
    });
    var cands = metaCandidates(),
      G = grouping(),
      ents = {},
      order = [];
    var keyOfId = function (id) {
      return typeof specimenKey === "function" ? specimenKey(id) : String(id);
    };
    var ent = function (eid) {
      var e = ents[eid];
      if (!e) {
        e = ents[eid] = { sid: eid, name: eid, members: [], t: null, total: 0, control: false, meta: {}, sgroup: [] };
        order.push(eid);
      }
      return e;
    };
    ids.forEach(function (id) {
      if (!id || (typeof sampleHidden !== "undefined" && sampleHidden[id])) return;
      var sm = SM[id] || {};
      if (sm.insilico === true) return;
      var e = ent(keyOfId(id));
      e.members.push(id);
      var rec = Object.assign({}, sm, RM[id] || {});
      var tr = Number(sm.total_reads != null ? sm.total_reads : sm.input_data && sm.input_data.reads);
      if (isFinite(tr) && tr > 0) e.total += tr;
      cands.forEach(function (c) {
        var vs = cellValues(rec, c.key);
        if (!vs.length) return;
        var arr = e.meta[c.key] || (e.meta[c.key] = []);
        vs.forEach(function (v) {
          if (arr.indexOf(v) < 0) arr.push(v);
        });
      });
      var t = dateOf(rec);
      if (t != null && (e.t == null || t < e.t)) e.t = t;
      if (typeof ttControlType === "function" && ttControlType(id)) e.control = true;
      if (G) {
        try {
          G.groupsOf(id).forEach(function (g) {
            if (e.sgroup.indexOf(g) < 0) e.sgroup.push(g);
          });
        } catch (err) {
          /* ignore */
        }
      }
    });
    // Detections in the current filtered view, summed per (organism, sample/specimen).
    var fd = typeof filteredData === "function" ? filteredData() || [] : [];
    var orgs = {},
      hasAttr = {};
    fd.forEach(function (r) {
      var sid = r["Specimen ID"];
      if (sid == null || sid === "") return;
      var nm = r["Detected Organism"] || r["Species Name"] || "";
      if (!nm) return;
      var e = ents[keyOfId(sid)] || ents[sid];
      if (!e) {
        e = ent(keyOfId(sid));
        e.members.push(sid);
      }
      var o = orgs[nm] || (orgs[nm] = { species: nm, attrs: {}, hc: false, ps: {} });
      GROUP_BY.forEach(function (g) {
        if (o.attrs[g[0]]) return;
        var v = attrOf(r, g[0]);
        if (v) {
          o.attrs[g[0]] = v;
          hasAttr[g[0]] = 1;
        }
      });
      if (r["High Consequence"] === true || /^(true|yes|1)$/i.test(String(r["High Consequence"] || ""))) o.hc = true;
      var p = o.ps[e.sid] || (o.ps[e.sid] = { reads: 0, rel: 0, rpm: 0, tass: 0 });
      p.reads += n0(r["# Reads Aligned"]);
      p.rel += n0(r["% Reads"]);
      p.rpm += n0(r.RPM);
      p.tass = Math.max(p.tass, n0(r["TASS Score"]));
    });
    var samples = order.map(function (eid) {
      var e = ents[eid],
        day = "";
      if (e.t != null && isFinite(e.t)) {
        var d = new Date(e.t);
        day = d.getFullYear() + "-" + pad(d.getMonth() + 1) + "-" + pad(d.getDate());
      }
      var mo = day ? +day.slice(5, 7) : 0;
      e.day = day;
      e.month = day ? day.slice(0, 7) : "";
      e.year = day ? day.slice(0, 4) : "";
      e.season = day ? SEASONS[Math.floor((mo % 12) / 3)] : "";
      e.stype = (e.meta.sample_type || [])[0] || "";
      return e;
    });
    var gbys = GROUP_BY.filter(function (g) {
      return hasAttr[g[0]];
    });
    if (!gbys.length) gbys = [GROUP_BY[0]];
    if (
      !gbys.some(function (g) {
        return g[0] === S.gby;
      })
    )
      S.gby = gbys[0][0];
    var taxa = Object.keys(orgs).map(function (k) {
      var o = orgs[k];
      return {
        species: o.species,
        group: o.attrs[S.gby] || "Unknown",
        cat: o.attrs["Microbial Category"] || "",
        hc: o.hc,
        ps: o.ps,
      };
    });
    var lvl = (document.getElementById("view-level") || {}).value || "Strain";
    return { samples: samples, taxa: taxa, cands: cands, hasGroups: !!G, gbys: gbys, viewLevel: lvl };
  }
  function pool() {
    return D.samples.filter(function (s) {
      return !(S.exclCtl && s.control);
    });
  }
  function fieldLabel(f) {
    if (TIME_FIELDS[f]) return TIME_FIELDS[f];
    if (f === "sgroup") return "Group by selection";
    if (f === "samples") return "Hand-picked samples";
    if (f.indexOf("m:") === 0) {
      var k = f.slice(2),
        c = D.cands.filter(function (x) {
          return x.key === k;
        })[0];
      return c ? c.label : k.replace(/_/g, " ");
    }
    return f;
  }
  function keysOf(s, f) {
    var v;
    if (f === "sgroup") v = s.sgroup;
    else if (f.indexOf("m:") === 0) v = s.meta[f.slice(2)] || [];
    else v = s[f] ? [s[f]] : [];
    if (v.length) return v;
    return TIME_FIELDS[f] ? ["(undated)"] : f === "sgroup" ? ["(ungrouped)"] : ["(none)"];
  }
  // values of a field with counts, in a sensible order (time-like fields chronologically)
  function fieldValues(f) {
    var cnt = {};
    pool().forEach(function (s) {
      keysOf(s, f).forEach(function (k) {
        cnt[k] = (cnt[k] || 0) + 1;
      });
    });
    var keys = Object.keys(cnt);
    keys.sort(function (a, b) {
      var na = isPlaceholder(a),
        nb = isPlaceholder(b);
      if (na !== nb) return na ? 1 : -1;
      if (f === "season") return SEASONS.indexOf(a) - SEASONS.indexOf(b);
      return natural(a, b);
    });
    return keys.map(function (k) {
      return { v: k, n: cnt[k] };
    });
  }
  function inSet(s, st) {
    if (st.field === "range") return !!s.day && (!st.from || s.day >= st.from) && (!st.to || s.day <= st.to);
    if (st.field === "samples") return (st.vals || []).indexOf(s.sid) >= 0;
    var ks = keysOf(s, st.field),
      vals = st.vals || [];
    for (var i = 0; i < ks.length; i++) if (vals.indexOf(ks[i]) >= 0) return true;
    return false;
  }
  function setSamples(st) {
    return pool().filter(function (s) {
      return inSet(s, st);
    });
  }
  function realValues(f) {
    return fieldValues(f).filter(function (x) {
      return !isPlaceholder(x.v);
    });
  }
  function availFields() {
    var f = [];
    if (D.hasGroups) f.push("sgroup");
    D.cands.forEach(function (c) {
      f.push("m:" + c.key);
    });
    f = f.concat(["day", "month", "season", "year", "range"]);
    f.push("samples");
    var dated = pool().some(function (s) {
      return s.day;
    });
    return f.filter(function (k) {
      if (k === "range") return dated;
      if (k === "samples") return true;
      return realValues(k).length > 0;
    });
  }
  function bestField() {
    var fields = availFields();
    for (var i = 0; i < fields.length; i++) {
      var f = fields[i];
      if (f === "samples" || f === "range") continue;
      if (realValues(f).length >= 2) return f;
    }
    return "samples";
  }
  function preset(kind, f) {
    f = f || bestField();
    var vals =
      f === "samples" || f === "range"
        ? pool()
            .slice()
            .sort(function (a, b) {
              return (a.t || 0) - (b.t || 0) || natural(a.name, b.name);
            })
            .map(function (s) {
              return { v: s.sid, n: 1, day: s.day };
            })
        : realValues(f);
    var mk = function (vs) {
      if (f === "range") {
        var days = vs
          .map(function (x) {
            return x.day;
          })
          .filter(Boolean)
          .sort();
        return { field: "range", from: days[0] || "", to: days[days.length - 1] || "", label: "" };
      }
      return {
        field: f,
        vals: vs.map(function (x) {
          return x.v;
        }),
        label: "",
      };
    };
    if (vals.length < 2) return [mk(vals), mk([])];
    if (kind === "each")
      return vals.slice(0, 4).map(function (x) {
        return mk([x]);
      });
    if (kind === "half") {
      var h = Math.ceil(vals.length / 2);
      return [mk(vals.slice(0, h)), mk(vals.slice(h))];
    }
    return [mk([vals[0]]), mk([vals[vals.length - 1]])];
  }
  function validate() {
    var fields = availFields();
    if (!Array.isArray(S.sets) || S.sets.length < 2) {
      S.sets = preset("fl");
      return;
    }
    S.sets = S.sets.slice(0, 4).map(function (st) {
      st = Object.assign({ field: "samples", vals: [], label: "" }, st);
      // A field with no values in this view falls back; chosen values are kept even when the current
      // filters don't have them, so switching back restores the set.
      if (fields.indexOf(st.field) < 0 && !(st.vals || []).length && st.field !== "range")
        st = { field: bestField(), vals: [], label: st.label || "" };
      return st;
    });
    if (
      S.tgroup &&
      !D.taxa.some(function (t) {
        return t.group === S.tgroup;
      })
    )
      S.tgroup = "";
    if (S.left >= S.sets.length) S.left = 0;
    if (S.right >= S.sets.length || S.right === S.left) S.right = S.left === 1 ? 0 : 1;
  }
  function setLabel(st, i) {
    if (st.label) return st.label;
    var L = "Set " + LETTERS[i];
    if (st.field === "range") return st.from || st.to ? (st.from || "…") + " → " + (st.to || "…") : L;
    var v = st.vals || [];
    if (!v.length) return L;
    if (st.field === "samples") return v.length + " sample" + (v.length === 1 ? "" : "s");
    return v.length === 1 ? v[0] : v.length === 2 ? v[0] + " + " + v[1] : v[0] + " + " + (v.length - 1) + " more";
  }

  /* ── statistics ───────────────────────────────────────────────────────── */
  function erf(x) {
    // Abramowitz & Stegun 7.1.26
    var s = x < 0 ? -1 : 1;
    x = Math.abs(x);
    var t = 1 / (1 + 0.3275911 * x);
    return (
      s *
      (1 -
        ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) *
          t *
          Math.exp(-x * x))
    );
  }
  // two-sided Mann–Whitney U (normal approximation with tie and continuity correction)
  function mwu(x, y) {
    var n1 = x.length,
      n2 = y.length,
      n = n1 + n2;
    if (n1 < 3 || n2 < 3) return NaN;
    var all = x
      .map(function (v) {
        return [v, 0];
      })
      .concat(
        y.map(function (v) {
          return [v, 1];
        }),
      )
      .sort(function (a, b) {
        return a[0] - b[0];
      });
    var r1 = 0,
      ties = 0,
      i = 0;
    while (i < n) {
      var j = i;
      while (j + 1 < n && all[j + 1][0] === all[i][0]) j++;
      var rank = (i + j) / 2 + 1,
        tc = j - i + 1;
      if (tc > 1) ties += tc * tc * tc - tc;
      for (var k = i; k <= j; k++) if (all[k][1] === 0) r1 += rank;
      i = j + 1;
    }
    var u = r1 - (n1 * (n1 + 1)) / 2,
      mu = (n1 * n2) / 2;
    var sd = Math.sqrt(((n1 * n2) / 12) * (n + 1 - ties / (n * (n - 1))));
    if (!(sd > 0)) return 1;
    var z = Math.max(0, Math.abs(u - mu) - 0.5) / sd;
    return Math.min(1, 2 * (1 - 0.5 * (1 + erf(z / Math.SQRT2))));
  }
  function bh(rows) {
    // Benjamini–Hochberg q-values
    var t = rows
        .filter(function (r) {
          return isFinite(r.p);
        })
        .sort(function (a, b) {
          return a.p - b.p;
        }),
      m = t.length,
      prev = 1;
    for (var i = m - 1; i >= 0; i--) {
      prev = Math.min(prev, (t[i].p * m) / (i + 1));
      t[i].q = prev;
    }
  }
  function median(a) {
    if (!a.length) return 0;
    var s = a.slice().sort(function (x, y) {
        return x - y;
      }),
      h = s.length >> 1;
    return s.length % 2 ? s[h] : (s[h - 1] + s[h]) / 2;
  }

  function units(level, all) {
    var u = {};
    D.taxa.forEach(function (t) {
      if (!all && S.tgroup && t.group !== S.tgroup) return;
      if (level === "org") {
        u[t.species] = { key: t.species, name: t.species, group: t.group, cat: t.cat, hc: t.hc, ps: t.ps };
        return;
      }
      var g =
        u[t.group] ||
        (u[t.group] = { key: t.group, name: t.group, group: t.group, cat: "", hc: false, ps: {}, members: 0 });
      g.members++;
      if (t.hc) g.hc = true;
      Object.keys(t.ps).forEach(function (sid) {
        var p = t.ps[sid],
          q = g.ps[sid] || (g.ps[sid] = { reads: 0, rel: 0, rpm: 0, tass: 0 });
        q.reads += p.reads;
        q.rel += p.rel;
        q.rpm += p.rpm;
        q.tass = Math.max(q.tass, p.tass);
      });
    });
    return Object.keys(u).map(function (k) {
      return u[k];
    });
  }
  function pv(p, m) {
    return !p ? 0 : m === "mreads" ? p.reads : p[m];
  }
  function rowsFor(level, sets, all) {
    return units(level, all)
      .map(function (u) {
        var v = sets.map(function (se) {
          var n = se.samples.length,
            det = 0,
            A = { det: [] };
          MK.forEach(function (m) {
            A[m] = [];
          });
          se.samples.forEach(function (s) {
            var p = u.ps[s.sid];
            if (p) det++;
            A.det.push(p ? 1 : 0);
            MK.forEach(function (m) {
              A[m].push(pv(p, m));
            });
          });
          var o = { n: n, det: det, A: A, prev: n ? (det / n) * 100 : 0 };
          MK.forEach(function (m) {
            o[m] = n ? sum(A[m]) / n : 0;
          });
          return o;
        });
        return { key: u.key, name: u.name, group: u.group, cat: u.cat, hc: u.hc, members: u.members, v: v };
      })
      .filter(function (r) {
        return r.v.some(function (x) {
          return x.det > 0;
        });
      });
  }
  function val(r, i) {
    return r.v[i] ? r.v[i][S.metric] : 0;
  }
  function eps() {
    var m = S.metric;
    if (m === "prev") {
      var n = Math.max.apply(
        null,
        C.sets.map(function (s) {
          return s.samples.length || 1;
        }),
      );
      return 50 / n;
    }
    return m === "rel" ? 0.001 : 0.5;
  }
  function arrOf(r, i) {
    var x = r.v[i];
    return S.metric === "prev" ? x.A.det : x.A[S.metric];
  }
  function compute() {
    var sets = S.sets.map(function (st, i) {
      return { i: i, st: st, color: COLORS[i], label: setLabel(st, i), samples: setSamples(st) };
    });
    C = { sets: sets };
    var L = S.left,
      R = S.right;
    var rows = rowsFor(S.level, sets);
    var e = eps();
    rows.forEach(function (r) {
      var a = val(r, L),
        b = val(r, R);
      r.a = a;
      r.b = b;
      r.d = b - a;
      r.lfc = Math.log2((b + e) / (a + e));
      r.onlyL = r.v[L].det > 0 && r.v[R].det === 0;
      r.onlyR = r.v[R].det > 0 && r.v[L].det === 0;
      r.p = mwu(arrOf(r, L), arrOf(r, R));
      r.q = NaN;
    });
    bh(rows);
    C.rows = rows;
    var seen = {},
      dup = 0;
    sets.forEach(function (se) {
      se.samples.forEach(function (s) {
        if (seen[s.sid]) dup++;
        seen[s.sid] = 1;
      });
    });
    C.dup = dup;
    var org = S.level === "org" && !S.tgroup ? rows : rowsFor("org", sets, true); // set KPIs use every organism
    sets.forEach(function (se, i) {
      se.rich = org.filter(function (r) {
        return r.v[i].det > 0;
      }).length;
      var p = org.map(function (r) {
          return r.v[i].rel;
        }),
        t = sum(p),
        h = 0;
      p.forEach(function (x) {
        if (x > 0 && t > 0) {
          var q = x / t;
          h -= q * Math.log(q);
        }
      });
      se.shannon = h;
      se.reads = se.samples.reduce(function (x, s) {
        return x + s.total;
      }, 0);
      se.aligned = org.reduce(function (x, r) {
        return x + sum(r.v[i].A.mreads);
      }, 0);
    });
    var sh = 0,
      oL = 0,
      oR = 0,
      mn = 0,
      sm = 0;
    org.forEach(function (r) {
      var x = r.v[L],
        y = r.v[R];
      if (x.det && y.det) sh++;
      else if (x.det) oL++;
      else if (y.det) oR++;
      mn += Math.min(x.rel, y.rel);
      sm += x.rel + y.rel;
    });
    C.pair = {
      shared: sh,
      onlyL: oL,
      onlyR: oR,
      jac: sh + oL + oR ? sh / (sh + oL + oR) : 0,
      bc: sm ? 1 - (2 * mn) / sm : 0,
    };
  }

  /* ── formatting ───────────────────────────────────────────────────────── */
  function kind(m) {
    return METRICS[m || S.metric].kind;
  }
  function fv(v, m) {
    if (!isFinite(v)) return "—";
    var k = kind(m),
      f = d3().format;
    if (k === "count") return v >= 100 ? fmtN(Math.round(v)) : f(".3~g")(v);
    if (k === "score") return f(".1f")(v);
    if (v === 0) return "0%";
    return (v < 0.01 ? f(".1~e")(v) : f(".3~g")(v)) + "%";
  }
  function fd(v) {
    var s = v > 0 ? "+" : v < 0 ? "−" : "±",
      a = Math.abs(v),
      k = kind(),
      f = d3().format;
    if (k === "count") return s + (a >= 100 ? fmtN(Math.round(a)) : f(".3~g")(a));
    if (k === "score") return s + f(".1f")(a);
    return s + (a < 0.01 && a > 0 ? f(".1~e")(a) : f(".3~g")(a)) + " pp";
  }
  function fold(l) {
    if (!isFinite(l)) return "—";
    var x = Math.pow(2, Math.abs(l));
    return (l >= 0 ? "×" : "÷") + (x >= 100 ? fmtN(Math.round(x)) : d3().format(".3~g")(x));
  }
  function fp(p) {
    return !isFinite(p) ? "—" : p < 0.001 ? "<0.001" : d3().format(".2~g")(p);
  }
  function tickFmt() {
    var k = kind(),
      f = d3().format;
    return k === "count"
      ? f("~s")
      : k === "score"
      ? f("~g")
      : function (v) {
          return f("~g")(v) + "%";
        };
  }
  function opt(v, l, cur) {
    return (
      '<option value="' + esc(v) + '"' + (String(v) === String(cur) ? " selected" : "") + ">" + esc(l) + "</option>"
    );
  }
  function trunc(s, n) {
    s = String(s);
    return s.length > n ? s.slice(0, n - 1) + "…" : s;
  }
  function setName(i) {
    return C.sets[i] ? C.sets[i].label : "";
  }
  function sw(i) {
    return '<i class="ab-sw" style="background:' + COLORS[i] + '"></i>';
  }
  function gbyLabel() {
    var g = GROUP_BY.filter(function (x) {
      return x[0] === S.gby;
    })[0];
    return g ? g[1] : S.gby;
  }
  var _gcol = {};
  function gcol(g) {
    if (S.gby === "Microbial Category" && typeof _xsCatColor === "function") return _xsCatColor(g);
    if (_gcol[g]) return _gcol[g];
    var pal = typeof PALETTE !== "undefined" && PALETTE.length ? PALETTE : COLORS;
    _gcol[g] = pal[Object.keys(_gcol).length % pal.length];
    return _gcol[g];
  }
  function sampleSub(s) {
    return [s.stype, s.day, s.members.length > 1 ? s.members.length + " libraries" : ""].filter(Boolean).join(" · ");
  }
  function hcBadge(r) {
    return r.hc ? ' <span class="ab-hc" title="High consequence">HC</span>' : "";
  }
  function rowTip(r) {
    var h =
      "<b>" +
      esc(r.name) +
      "</b>" +
      (r.hc ? ' <span style="color:#ffb4a8">· high consequence</span>' : "") +
      '<div style="color:#bbb;font-size:11px">' +
      esc(S.level === "org" ? gbyLabel() + ": " + r.group : (r.members || 0) + " organisms") +
      '</div><table style="margin-top:4px;font-size:11px">';
    C.sets.forEach(function (se, i) {
      var x = r.v[i];
      h +=
        "<tr><td>" +
        sw(i) +
        " " +
        esc(trunc(se.label, 28)) +
        '</td><td style="text-align:right;padding-left:8px"><b>' +
        fv(x[S.metric]) +
        '</b></td><td style="color:#bbb;padding-left:6px">' +
        x.det +
        "/" +
        x.n +
        " samples</td></tr>";
    });
    h +=
      '</table><div style="margin-top:3px">' +
      esc(trunc(setName(S.right), 20)) +
      " vs " +
      esc(trunc(setName(S.left), 20)) +
      ": <b>" +
      fd(r.d) +
      "</b> · " +
      fold(r.lfc) +
      (isFinite(r.p) ? " · p " + fp(r.p) : "") +
      '</div><div style="color:#bbb;font-size:10px;margin-top:2px">Click for the per-sample view</div>';
    return h;
  }

  /* ── render ───────────────────────────────────────────────────────────── */
  var BUSY = false;
  // Replacing the pane's HTML blurs a focused input, which can fire its change
  // event and ask for another render while the first is still running.
  function render(host) {
    if (BUSY || !host) return;
    BUSY = true;
    try {
      renderNow(host);
    } finally {
      BUSY = false;
    }
  }
  function renderNow(host) {
    HOST = host;
    untip();
    if (!d3()) {
      host.innerHTML =
        '<div class="ab-empty">The charting library (d3) did not load, so the Contrast view is unavailable.</div>';
      return;
    }
    D = prepare();
    _gcol = {};
    wire(host);
    if (D.samples.length < 2) {
      host.innerHTML =
        '<div class="ab-wrap"><div class="ab-h">Contrast</div><div class="ab-empty">Contrast compares sets of samples — it needs at least two samples in view (check hidden samples in the right panel).</div></div>';
      return;
    }
    validate();
    compute();
    save();
    var keepScroll = {};
    Array.prototype.forEach.call(host.querySelectorAll(".ab-vals"), function (el) {
      keepScroll[el.getAttribute("data-i")] = el.scrollTop;
    });
    host.innerHTML =
      '<div class="ab-wrap">' +
      '<div class="ab-h">Contrast sets</div>' +
      '<p class="ab-sub">Group samples into sets — by any metadata column (sample type, site, run, host, …), the shared Group by selection, collection date or period, or hand-picked samples — and compare what was detected in each. Uses the detections and samples in the current view (right-panel filters, TASS cutoff, <b>' +
      esc(D.viewLevel) +
      "</b> level" +
      (typeof specimenMergeEnabled !== "undefined" && specimenMergeEnabled ? ", merged specimens" : "") +
      ").</p>" +
      setsHtml() +
      barHtml() +
      kpiHtml() +
      '<div class="ab-card ab-main"><div class="ab-vtabs">' +
      VIEWS.map(function (v) {
        return (
          '<button type="button" data-view="' +
          v[0] +
          '" class="' +
          (S.view === v[0] ? "on" : "") +
          '">' +
          v[1] +
          "</button>"
        );
      }).join("") +
      '<span class="ab-vctl">' +
      viewCtl() +
      '</span></div><div class="ab-chart"></div><div class="ab-note"></div></div>' +
      '<div class="ab-detail-host"></div></div>';
    Array.prototype.forEach.call(host.querySelectorAll(".ab-vals"), function (el) {
      var k = el.getAttribute("data-i");
      if (keepScroll[k]) el.scrollTop = keepScroll[k];
    });
    drawView();
    drawDetail();
  }
  function setsHtml() {
    var fields = availFields(),
      pf =
        S.pfield && fields.indexOf(S.pfield) >= 0 && S.pfield !== "samples" && S.pfield !== "range"
          ? S.pfield
          : bestField();
    var h =
      '<div class="ab-card ab-sets"><div class="ab-preset"><b>Quick sets</b> by <select data-k="pfield">' +
      fields
        .filter(function (f) {
          return f !== "samples" && f !== "range";
        })
        .map(function (f) {
          return opt(f, fieldLabel(f), pf);
        })
        .join("") +
      "</select>" +
      '<button type="button" class="ab-btn" data-preset="fl" title="Set A = the first value, Set B = the last">First vs last</button>' +
      '<button type="button" class="ab-btn" data-preset="half" title="Earlier half of the values vs the later half">Earlier vs later half</button>' +
      '<button type="button" class="ab-btn" data-preset="each" title="One set per value (up to 4)">One set each (≤4)</button>' +
      (C.dup
        ? '<span class="ab-warn" title="A sample in two sets counts in both">⚠ ' +
          C.dup +
          " sample" +
          (C.dup === 1 ? " is" : "s are") +
          " in more than one set</span>"
        : "") +
      '</div><div class="ab-setgrid">';
    S.sets.forEach(function (st, i) {
      var n = C.sets[i].samples.length;
      h +=
        '<div class="ab-set" style="--c:' +
        COLORS[i] +
        '"><div class="ab-sh">' +
        sw(i) +
        '<input type="text" class="ab-lab" data-i="' +
        i +
        '" value="' +
        esc(st.label) +
        '" placeholder="' +
        esc(setLabel(Object.assign({}, st, { label: "" }), i)) +
        '" title="Name this set">' +
        '<span class="ab-n' +
        (n ? "" : " zero") +
        '">' +
        n +
        " sample" +
        (n === 1 ? "" : "s") +
        "</span>" +
        (S.sets.length > 2
          ? '<button type="button" class="ab-x" data-rm="' + i + '" title="Remove this set">×</button>'
          : "") +
        "</div>" +
        '<select data-i="' +
        i +
        '" data-k="field" class="ab-field">' +
        (fields.indexOf(st.field) < 0 ? opt(st.field, fieldLabel(st.field) + " (not in view)", st.field) : "") +
        fields
          .map(function (f) {
            return opt(f, fieldLabel(f), st.field);
          })
          .join("") +
        "</select>" +
        setBody(st, i) +
        "</div>";
    });
    if (S.sets.length < 4)
      h += '<button type="button" class="ab-add" data-add="1">+ Add set ' + LETTERS[S.sets.length] + "</button>";
    return h + "</div></div>";
  }
  function setBody(st, i) {
    if (st.field === "range") {
      var days = pool()
        .map(function (s) {
          return s.day;
        })
        .filter(Boolean)
        .sort();
      var mn = days[0] || "",
        mx = days[days.length - 1] || "";
      return (
        '<div class="ab-range"><input type="date" data-i="' +
        i +
        '" data-k="from" value="' +
        esc(st.from || "") +
        '" min="' +
        mn +
        '" max="' +
        mx +
        '"><span>to</span>' +
        '<input type="date" data-i="' +
        i +
        '" data-k="to" value="' +
        esc(st.to || "") +
        '" min="' +
        mn +
        '" max="' +
        mx +
        '"></div>' +
        '<div class="ab-hint">Data spans ' +
        esc(mn || "—") +
        " → " +
        esc(mx || "—") +
        "</div>"
      );
    }
    var sel = {};
    (st.vals || []).forEach(function (v) {
      sel[v] = 1;
    });
    var items,
      q = (S.q[i] || "").toLowerCase();
    if (st.field === "samples") {
      items = pool()
        .slice()
        .sort(function (a, b) {
          return natural(a.stype, b.stype) || (a.t || 0) - (b.t || 0) || natural(a.name, b.name);
        })
        .filter(function (s) {
          if (!q) return true;
          var hay =
            s.name +
            " " +
            s.members.join(" ") +
            " " +
            s.day +
            " " +
            Object.keys(s.meta)
              .map(function (k) {
                return s.meta[k].join(" ");
              })
              .join(" ");
          return hay.toLowerCase().indexOf(q) >= 0;
        })
        .map(function (s) {
          return { v: s.sid, l: s.name, sub: sampleSub(s) };
        });
    } else {
      items = fieldValues(st.field)
        .map(function (x) {
          return { v: x.v, l: x.v, sub: x.n + " sample" + (x.n === 1 ? "" : "s") };
        })
        .filter(function (x) {
          return !q || x.l.toLowerCase().indexOf(q) >= 0;
        });
    }
    var many = st.field === "samples" || items.length > 8;
    var present = {};
    if (st.field === "samples")
      pool().forEach(function (s) {
        present[s.sid] = 1;
      });
    else
      fieldValues(st.field).forEach(function (x) {
        present[x.v] = 1;
      });
    var chosen = (st.vals || []).filter(function (v) {
        return present[v];
      }).length,
      away = (st.vals || []).length - chosen;
    return (
      (many
        ? '<input type="search" class="ab-q" data-i="' +
          i +
          '" value="' +
          esc(S.q[i] || "") +
          '" placeholder="Filter…">'
        : "") +
      '<div class="ab-vacts"><a data-all="' +
      i +
      '">All' +
      (q ? " listed" : "") +
      '</a> · <a data-none="' +
      i +
      '">None</a> · <span>' +
      chosen +
      " chosen</span>" +
      (away
        ? ' · <span class="ab-away" title="Chosen, but not in the current filters / hidden samples — they come back when you switch">' +
          away +
          " not in this view</span>"
        : "") +
      "</div>" +
      '<div class="ab-vals" data-i="' +
      i +
      '">' +
      (items.length
        ? items
            .map(function (x) {
              return (
                '<label class="ab-v"><input type="checkbox" data-i="' +
                i +
                '" data-v="' +
                esc(x.v) +
                '"' +
                (sel[x.v] ? " checked" : "") +
                '><span class="ab-vl" title="' +
                esc(x.l) +
                '">' +
                esc(x.l) +
                '</span><span class="ab-vs">' +
                esc(x.sub) +
                "</span></label>"
              );
            })
            .join("")
        : '<div class="ab-hint">Nothing matches.</div>') +
      "</div>"
    );
  }
  function barHtml() {
    var groups = {};
    D.taxa.forEach(function (t) {
      groups[t.group] = 1;
    });
    var pair =
      S.sets.length > 2
        ? '<label>Compare <select data-k="left">' +
          S.sets
            .map(function (st, i) {
              return opt(i, LETTERS[i] + " · " + trunc(C.sets[i].label, 22), S.left);
            })
            .join("") +
          "</select> with " +
          '<select data-k="right">' +
          S.sets
            .map(function (st, i) {
              return opt(i, LETTERS[i] + " · " + trunc(C.sets[i].label, 22), S.right);
            })
            .join("") +
          "</select></label>"
        : "";
    return (
      '<div class="ab-bar">' +
      '<span class="ab-seg" role="group" aria-label="Level"><button type="button" data-level="org" class="' +
      (S.level === "org" ? "on" : "") +
      '">Organisms</button><button type="button" data-level="grp" class="' +
      (S.level === "grp" ? "on" : "") +
      '">Groups</button></span>' +
      '<label>Group by <select data-k="gby">' +
      D.gbys
        .map(function (g) {
          return opt(g[0], g[1], S.gby);
        })
        .join("") +
      "</select></label>" +
      '<label>Measure <select data-k="metric">' +
      Object.keys(METRICS)
        .map(function (k) {
          return opt(k, METRICS[k].label, S.metric);
        })
        .join("") +
      "</select></label>" +
      '<label>Taxa <select data-k="tgroup">' +
      opt("", "All " + gbyLabel().toLowerCase() + " groups", S.tgroup) +
      Object.keys(groups)
        .sort(natural)
        .map(function (g) {
          return opt(g, g, S.tgroup);
        })
        .join("") +
      "</select></label>" +
      '<label>Top <input type="number" data-k="topN" min="3" max="200" value="' +
      S.topN +
      '"></label>' +
      '<label>Scale <select data-k="scale">' +
      opt("lin", "Linear", S.scale) +
      opt("sqrt", "Square root", S.scale) +
      "</select></label>" +
      '<label class="ab-chk" title="Samples whose metadata marks them as a negative / positive control"><input type="checkbox" data-k="exclCtl"' +
      (S.exclCtl ? " checked" : "") +
      "> Leave out controls</label>" +
      pair +
      "</div>"
    );
  }
  function kpiHtml() {
    var h = '<div class="ab-kpis">',
      f = d3().format;
    C.sets.forEach(function (se, i) {
      h +=
        '<div class="ab-kpi" style="--c:' +
        COLORS[i] +
        '"><div class="ab-kt">' +
        sw(i) +
        esc(trunc(se.label, 34)) +
        "</div>" +
        '<div class="ab-kv"><b>' +
        se.samples.length +
        "</b> sample" +
        (se.samples.length === 1 ? "" : "s") +
        " · <b>" +
        fmtN(se.rich) +
        "</b> organisms</div>" +
        '<div class="ab-kv">' +
        (se.reads ? fmtN(se.reads) + " reads · " : "") +
        fmtN(se.aligned) +
        " aligned · Shannon " +
        f(".2f")(se.shannon || 0) +
        "</div></div>";
    });
    var P = C.pair,
      L = S.left,
      R = S.right;
    h +=
      '<div class="ab-kpi ab-kpair"><div class="ab-kt">' +
      sw(L) +
      LETTERS[L] +
      " vs " +
      sw(R) +
      LETTERS[R] +
      "</div>" +
      '<div class="ab-kv"><b>' +
      P.shared +
      '</b> shared · <b style="color:' +
      COLORS[L] +
      '">' +
      P.onlyL +
      "</b> only " +
      LETTERS[L] +
      ' · <b style="color:' +
      COLORS[R] +
      '">' +
      P.onlyR +
      "</b> only " +
      LETTERS[R] +
      "</div>" +
      '<div class="ab-kv" title="Jaccard = shared ÷ all organisms detected in either set. Bray–Curtis dissimilarity on the mean % of reads profiles: 0 = identical, 1 = nothing in common.">Jaccard ' +
      f(".2f")(P.jac) +
      " · Bray–Curtis " +
      f(".2f")(P.bc) +
      "</div></div>";
    return h + "</div>";
  }
  function viewCtl() {
    if (S.view === "mirror")
      return (
        '<label>Order <select data-k="sort">' +
        opt("delta", "Biggest change", S.sort) +
        opt("sum", "Largest overall", S.sort) +
        opt("left", "Largest in " + LETTERS[S.left], S.sort) +
        opt("right", "Largest in " + LETTERS[S.right], S.sort) +
        opt("name", "A → Z", S.sort) +
        "</select></label>"
      );
    if (S.view === "delta")
      return (
        '<label>Show <select data-k="dmode">' +
        opt("diff", "Difference (" + LETTERS[S.right] + " − " + LETTERS[S.left] + ")", S.dmode) +
        opt("fold", "Fold change (log₂ " + LETTERS[S.right] + " ÷ " + LETTERS[S.left] + ")", S.dmode) +
        "</select></label>"
      );
    if (S.view === "radar")
      return (
        '<label>Spokes <select data-k="raxes">' +
        opt("grp", gbyLabel() + " groups", S.raxes) +
        opt("org", "Top organisms", S.raxes) +
        '</select></label><label>How many <input type="number" data-k="rk" min="3" max="16" value="' +
        S.rk +
        '"></label>'
      );
    if (S.view === "scatter")
      return (
        "<label>Top " +
        S.topN +
        ' by <select data-k="ssort">' +
        opt("sum", "Largest overall", S.ssort) +
        opt("delta", "Biggest change", S.ssort) +
        '</select></label><label>Axes <select data-k="sscale">' +
        opt("log", "Log-like (spreads small values)", S.sscale) +
        opt("sqrt", "Square root", S.sscale) +
        opt("lin", "Linear", S.sscale) +
        "</select></label>"
      );
    if (S.view === "table")
      return (
        '<label class="ab-chk"><input type="checkbox" data-k="tall"' +
        (S.tall ? " checked" : "") +
        '> Show every row</label><button type="button" class="ab-btn" data-csv="1">⬇ CSV</button>'
      );
    return "";
  }

  /* ── views ────────────────────────────────────────────────────────────── */
  function width(el) {
    return Math.max(520, (el.clientWidth || (HOST && HOST.clientWidth) || 900) - 4);
  }
  function ranked(by) {
    var rows = C.rows.slice();
    var f = {
      delta: function (r) {
        return -Math.abs(r.d);
      },
      sum: function (r) {
        return -(r.a + r.b);
      },
      left: function (r) {
        return -r.a;
      },
      right: function (r) {
        return -r.b;
      },
      absl: function (r) {
        return -Math.abs(r.lfc);
      },
    }[by];
    if (by === "name" || !f)
      rows.sort(function (a, b) {
        return natural(a.name, b.name);
      });
    else
      rows.sort(function (a, b) {
        return f(a) - f(b) || natural(a.name, b.name);
      });
    return rows;
  }
  function scaleFn() {
    return S.scale === "sqrt" ? d3().scaleSqrt() : d3().scaleLinear();
  }
  function note(t) {
    var n = HOST.querySelector(".ab-note");
    if (n) n.innerHTML = t || "";
  }
  function itCls() {
    return S.level === "org" ? " it" : "";
  }
  function needTwo(el) {
    var L = C.sets[S.left],
      R = C.sets[S.right];
    if (!L.samples.length || !R.samples.length) {
      el.innerHTML =
        '<div class="ab-empty">Choose samples for ' +
        (!L.samples.length ? "Set " + LETTERS[S.left] : "") +
        (!L.samples.length && !R.samples.length ? " and " : "") +
        (!R.samples.length ? "Set " + LETTERS[S.right] : "") +
        " above to compare.</div>";
      return true;
    }
    if (!C.rows.length) {
      el.innerHTML = '<div class="ab-empty">No organisms detected in these sets with the current filters.</div>';
      return true;
    }
    return false;
  }
  function drawView() {
    var el = HOST.querySelector(".ab-chart");
    if (!el) return;
    el.innerHTML = "";
    note("");
    if (S.view === "radar") return viewRadar(el);
    if (needTwo(el)) return;
    if (S.view === "mirror") return viewMirror(el);
    if (S.view === "delta") return viewDelta(el);
    if (S.view === "scatter") return viewScatter(el);
    if (S.view === "table") return viewTable(el);
  }
  function bindRow(sel) {
    sel
      .style("cursor", "pointer")
      .on("mousemove", function (ev, r) {
        tip(ev, rowTip(r));
      })
      .on("mouseleave", untip)
      .on("click", function (ev, r) {
        untip();
        S.sel = r.key;
        drawDetail(true);
      });
  }
  function svgFor(el, W, H) {
    return d3()
      .select(el)
      .append("svg")
      .attr("class", "ab-svg")
      .attr("width", W)
      .attr("height", H)
      .attr("viewBox", "0 0 " + W + " " + H);
  }
  function viewMirror(el) {
    var d = d3(),
      L = S.left,
      R = S.right,
      rows = ranked(S.sort).slice(0, S.topN);
    var W = width(el),
      lab = Math.min(250, Math.max(150, W * 0.24)),
      half = (W - lab) / 2 - 8,
      valW = 58;
    var rowH = 22,
      top = 34,
      H = top + rows.length * rowH + 10,
      cx = W / 2;
    var mx =
      d.max(rows, function (r) {
        return Math.max(r.a, r.b);
      }) || 1;
    var x = scaleFn()
      .domain([0, mx])
      .range([0, half - valW]);
    var svg = svgFor(el, W, H);
    svg
      .append("text")
      .attr("x", cx - lab / 2 - 4)
      .attr("y", 16)
      .attr("text-anchor", "end")
      .attr("class", "ab-hd")
      .attr("fill", COLORS[L])
      .text("◀ " + trunc(setName(L), 40));
    svg
      .append("text")
      .attr("x", cx + lab / 2 + 4)
      .attr("y", 16)
      .attr("class", "ab-hd")
      .attr("fill", COLORS[R])
      .text(trunc(setName(R), 40) + " ▶");
    svg
      .append("text")
      .attr("x", cx)
      .attr("y", 16)
      .attr("text-anchor", "middle")
      .attr("class", "ab-sub")
      .text(METRICS[S.metric].short);
    var g = svg
      .selectAll("g.ab-r")
      .data(rows)
      .enter()
      .append("g")
      .attr("class", "ab-r")
      .attr("transform", function (r, i) {
        return "translate(0," + (top + i * rowH) + ")";
      });
    g.append("rect")
      .attr("class", "ab-band")
      .attr("x", 0)
      .attr("width", W)
      .attr("height", rowH - 2)
      .attr("fill", function (r, i) {
        return i % 2 ? "#fafbfc" : "#fff";
      });
    g.append("rect")
      .attr("x", function (r) {
        return cx - lab / 2 - x(r.a);
      })
      .attr("y", 3)
      .attr("height", rowH - 8)
      .attr("width", function (r) {
        return Math.max(r.a > 0 ? 1 : 0, x(r.a));
      })
      .attr("fill", COLORS[L])
      .attr("rx", 2);
    g.append("rect")
      .attr("x", cx + lab / 2)
      .attr("y", 3)
      .attr("height", rowH - 8)
      .attr("width", function (r) {
        return Math.max(r.b > 0 ? 1 : 0, x(r.b));
      })
      .attr("fill", COLORS[R])
      .attr("rx", 2);
    g.append("text")
      .attr("class", "ab-val")
      .attr("x", function (r) {
        return cx - lab / 2 - x(r.a) - 4;
      })
      .attr("y", rowH / 2 + 2)
      .attr("text-anchor", "end")
      .text(function (r) {
        return r.v[L].det ? fv(r.a) : "—";
      });
    g.append("text")
      .attr("class", "ab-val")
      .attr("x", function (r) {
        return cx + lab / 2 + x(r.b) + 4;
      })
      .attr("y", rowH / 2 + 2)
      .text(function (r) {
        return r.v[R].det ? fv(r.b) : "—";
      });
    g.append("text")
      .attr("class", "ab-nm" + itCls())
      .attr("x", cx)
      .attr("y", rowH / 2 + 2)
      .attr("text-anchor", "middle")
      .style("fill", function (r) {
        return r.onlyL ? COLORS[L] : r.onlyR ? COLORS[R] : null;
      })
      .style("font-weight", function (r) {
        return r.onlyL || r.onlyR ? 700 : null;
      })
      .text(function (r) {
        return trunc(r.name, Math.floor(lab / 7));
      })
      .append("title")
      .text(function (r) {
        return r.name;
      });
    svg
      .append("line")
      .attr("x1", cx - lab / 2)
      .attr("x2", cx - lab / 2)
      .attr("y1", top - 4)
      .attr("y2", H - 8)
      .attr("stroke", "#cfd8de");
    svg
      .append("line")
      .attr("x1", cx + lab / 2)
      .attr("x2", cx + lab / 2)
      .attr("y1", top - 4)
      .attr("y2", H - 8)
      .attr("stroke", "#cfd8de");
    bindRow(g);
    note(
      "Top " +
        rows.length +
        " of " +
        C.rows.length +
        " " +
        (S.level === "org" ? "organisms" : "groups") +
        " · names in a set colour were found only in that set · " +
        (S.scale === "sqrt" ? "square-root" : "linear") +
        " bar scale, shared by both sides.",
    );
  }
  function viewDelta(el) {
    var d = d3(),
      L = S.left,
      R = S.right,
      isFold = S.dmode === "fold";
    var rows = ranked(isFold ? "absl" : "delta")
      .slice(0, S.topN)
      .sort(function (a, b) {
        return (isFold ? b.lfc - a.lfc : b.d - a.d) || natural(a.name, b.name);
      });
    var W = width(el),
      lab = Math.min(260, Math.max(150, W * 0.26)),
      right = 120,
      rowH = 22,
      top = 40,
      H = top + rows.length * rowH + 12;
    var m =
      d.max(rows, function (r) {
        return Math.abs(isFold ? r.lfc : r.d);
      }) || 1;
    var x = d
      .scaleLinear()
      .domain([-m, m])
      .range([lab, W - right])
      .nice();
    var svg = svgFor(el, W, H);
    var ax = d
      .axisTop(x)
      .ticks(Math.max(3, Math.floor((W - lab - right) / 90)))
      .tickFormat(function (v) {
        return isFold ? (v === 0 ? "×1" : fold(v)) : kind() === "count" ? d.format("~s")(v) : d.format("~g")(v);
      });
    svg
      .append("g")
      .attr("transform", "translate(0," + (top - 6) + ")")
      .attr("class", "ab-axis")
      .call(ax);
    svg
      .append("text")
      .attr("x", x(0) - 6)
      .attr("y", 12)
      .attr("text-anchor", "end")
      .attr("class", "ab-hd")
      .attr("fill", COLORS[L])
      .text("◀ higher in " + trunc(setName(L), 30));
    svg
      .append("text")
      .attr("x", x(0) + 6)
      .attr("y", 12)
      .attr("class", "ab-hd")
      .attr("fill", COLORS[R])
      .text("higher in " + trunc(setName(R), 30) + " ▶");
    var g = svg
      .selectAll("g.ab-r")
      .data(rows)
      .enter()
      .append("g")
      .attr("class", "ab-r")
      .attr("transform", function (r, i) {
        return "translate(0," + (top + i * rowH) + ")";
      });
    var vOf = function (r) {
      return isFold ? r.lfc : r.d;
    };
    g.append("rect")
      .attr("x", 0)
      .attr("width", W)
      .attr("height", rowH - 2)
      .attr("fill", function (r, i) {
        return i % 2 ? "#fafbfc" : "#fff";
      });
    g.append("rect")
      .attr("y", 4)
      .attr("height", rowH - 10)
      .attr("rx", 2)
      .attr("x", function (r) {
        return Math.min(x(0), x(vOf(r)));
      })
      .attr("width", function (r) {
        return Math.max(1, Math.abs(x(vOf(r)) - x(0)));
      })
      .attr("fill", function (r) {
        return vOf(r) >= 0 ? COLORS[R] : COLORS[L];
      });
    g.append("text")
      .attr("class", "ab-nm" + itCls())
      .attr("x", lab - 8)
      .attr("y", rowH / 2 + 2)
      .attr("text-anchor", "end")
      .text(function (r) {
        return trunc(r.name, Math.floor(lab / 7));
      })
      .append("title")
      .text(function (r) {
        return r.name;
      });
    g.append("text")
      .attr("class", "ab-val")
      .attr("y", rowH / 2 + 2)
      .attr("x", function (r) {
        return vOf(r) >= 0 ? x(vOf(r)) + 4 : x(vOf(r)) - 4;
      })
      .attr("text-anchor", function (r) {
        return vOf(r) >= 0 ? "start" : "end";
      })
      .text(function (r) {
        return (
          (isFold ? fold(r.lfc) : fd(r.d)) + (r.onlyR ? "  only " + LETTERS[R] : r.onlyL ? "  only " + LETTERS[L] : "")
        );
      });
    svg
      .append("line")
      .attr("x1", x(0))
      .attr("x2", x(0))
      .attr("y1", top - 6)
      .attr("y2", H - 8)
      .attr("stroke", "#56707f");
    bindRow(g);
    var unit = kind() === "pct" ? " (percentage points)" : "";
    note(
      isFold
        ? "Fold change of " +
            METRICS[S.metric].short.toLowerCase() +
            " (log₂ scale; a small offset keeps organisms found in only one set finite). Top " +
            rows.length +
            " by size of change."
        : "Difference in " +
            METRICS[S.metric].short.toLowerCase() +
            unit +
            ", " +
            LETTERS[R] +
            " minus " +
            LETTERS[L] +
            ". Top " +
            rows.length +
            " by size of change.",
    );
  }
  function logTicks(mx) {
    var t = [0],
      m = S.metric;
    var base =
      kind() === "count"
        ? [1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7]
        : m === "rel"
        ? [0.01, 0.1, 1, 10, 100]
        : [1, 2, 5, 10, 20, 50, 100];
    base.forEach(function (v) {
      if (v <= mx) t.push(v);
    });
    return t;
  }
  function viewRadar(el) {
    var d = d3();
    if (
      !C.sets.some(function (se) {
        return se.samples.length;
      })
    ) {
      el.innerHTML = '<div class="ab-empty">Choose samples for the sets above.</div>';
      return;
    }
    var grpSpokes = S.raxes === "grp" && S.level === "org",
      fellBack = false;
    var spokes = function (src) {
      src.forEach(function (r) {
        r.tot = C.sets.reduce(function (a, se, i) {
          return a + val(r, i);
        }, 0);
      });
      return src
        .filter(function (r) {
          return r.tot > 0;
        })
        .sort(function (a, b) {
          return b.tot - a.tot;
        })
        .slice(0, Math.max(3, Math.min(16, S.rk)));
    };
    var rows = spokes(grpSpokes ? rowsFor("grp", C.sets) : C.rows.slice());
    // Few taxtriage runs span 3+ microbial categories: fall back to organism spokes rather than a dead end.
    if (grpSpokes && rows.length < 3) {
      var alt = spokes(C.rows.slice());
      if (alt.length >= 3) {
        rows = alt;
        grpSpokes = false;
        fellBack = true;
      }
    }
    if (rows.length < 3) {
      el.innerHTML =
        '<div class="ab-empty">A radar needs at least 3 ' +
        (S.raxes === "grp" ? gbyLabel().toLowerCase() + " groups" : "organisms") +
        " — try another Group by, fewer filters, or Top organisms.</div>";
      return;
    }
    var W = width(el),
      size = Math.min(W, 640),
      H = size,
      cx = W / 2,
      cy = H / 2 + 6,
      R0 = size / 2 - 92;
    var mx =
      d.max(rows, function (r) {
        return d.max(C.sets, function (se, i) {
          return val(r, i);
        });
      }) || 1;
    var rs = scaleFn().domain([0, mx]).range([0, R0]);
    var ang = function (i) {
      return (i / rows.length) * 2 * Math.PI;
    };
    var svg = svgFor(el, W, H);
    var g = svg.append("g").attr("transform", "translate(" + cx + "," + cy + ")");
    [0.25, 0.5, 0.75, 1].forEach(function (f) {
      var v = S.scale === "sqrt" ? mx * f * f : mx * f;
      g.append("circle").attr("r", rs(v)).attr("fill", "none").attr("stroke", "#e3e9ee");
      g.append("text")
        .attr("x", 3)
        .attr("y", -rs(v) - 2)
        .attr("class", "ab-tick")
        .text(fv(v));
    });
    rows.forEach(function (r, i) {
      var a = ang(i);
      g.append("line")
        .attr("x2", Math.sin(a) * R0)
        .attr("y2", -Math.cos(a) * R0)
        .attr("stroke", "#e3e9ee");
      var t = g
        .append("text")
        .attr("x", Math.sin(a) * (R0 + 14))
        .attr("y", -Math.cos(a) * (R0 + 14) + 4)
        .attr("class", "ab-nm" + (!grpSpokes && S.level === "org" ? " it" : ""))
        .attr("text-anchor", Math.abs(Math.sin(a)) < 0.2 ? "middle" : Math.sin(a) > 0 ? "start" : "end")
        .text(trunc(r.name, 24));
      t.append("title").text(r.name);
      if (!grpSpokes)
        t.style("cursor", "pointer").on("click", function () {
          S.sel = r.key;
          drawDetail(true);
        });
    });
    C.sets.forEach(function (se, k) {
      if (!se.samples.length) return;
      var pts = rows.map(function (r, i) {
        return [ang(i), rs(val(r, k))];
      });
      g.append("path")
        .attr("d", d.lineRadial().curve(d.curveLinearClosed)(pts))
        .attr("fill", COLORS[k])
        .attr("fill-opacity", 0.12)
        .attr("stroke", COLORS[k])
        .attr("stroke-width", 2);
      g.selectAll(null)
        .data(rows)
        .enter()
        .append("circle")
        .attr("r", 3.5)
        .attr("fill", COLORS[k])
        .attr("cx", function (r, i) {
          return Math.sin(ang(i)) * rs(val(r, k));
        })
        .attr("cy", function (r, i) {
          return -Math.cos(ang(i)) * rs(val(r, k));
        })
        .on("mousemove", function (ev, r) {
          tip(
            ev,
            "<b>" +
              esc(r.name) +
              "</b><br>" +
              sw(k) +
              " " +
              esc(se.label) +
              ": <b>" +
              fv(val(r, k)) +
              "</b> (" +
              r.v[k].det +
              "/" +
              r.v[k].n +
              " samples)",
          );
        })
        .on("mouseleave", untip);
    });
    el.insertAdjacentHTML(
      "beforeend",
      '<div class="ab-legend">' +
        C.sets
          .map(function (se, i) {
            return se.samples.length ? "<span>" + sw(i) + esc(se.label) + "</span>" : "";
          })
          .join("") +
        "</div>",
    );
    note(
      (fellBack
        ? "Fewer than 3 " + gbyLabel().toLowerCase() + " groups in these sets, so the spokes are organisms instead. "
        : "") +
        (grpSpokes
          ? "One spoke per " + gbyLabel().toLowerCase() + " group"
          : "One spoke per " + (S.level === "org" ? "organism" : "group")) +
        " (top " +
        rows.length +
        " by " +
        METRICS[S.metric].short.toLowerCase() +
        " over all sets) · every set drawn" +
        (S.scale === "sqrt" ? " · square-root radius" : "") +
        ".",
    );
  }
  function viewScatter(el) {
    var d = d3(),
      L = S.left,
      R = S.right,
      rows = ranked(S.ssort === "delta" ? "delta" : "sum").slice(0, S.topN);
    var W = width(el),
      H = Math.min(560, Math.max(380, W * 0.55)),
      m = { l: 64, r: 24, t: 16, b: 46 };
    var mx =
      d.max(rows, function (r) {
        return Math.max(r.a, r.b);
      }) || 1;
    var sfn = function () {
      return S.sscale === "log"
        ? d.scaleSymlog().constant(S.metric === "rel" ? 0.01 : 1)
        : S.sscale === "sqrt"
        ? d.scaleSqrt()
        : d.scaleLinear();
    };
    var x = sfn()
        .domain([0, mx * 1.05])
        .range([m.l, W - m.r]),
      y = sfn()
        .domain([0, mx * 1.05])
        .range([H - m.b, m.t]);
    if (S.sscale !== "log") {
      x.nice();
      y.nice();
    }
    var svg = svgFor(el, W, H),
      tf = tickFmt();
    var tv = S.sscale === "log" ? logTicks(mx * 1.05) : null;
    svg
      .append("g")
      .attr("class", "ab-axis")
      .attr("transform", "translate(0," + (H - m.b) + ")")
      .call(tv ? d.axisBottom(x).tickValues(tv).tickFormat(tf) : d.axisBottom(x).ticks(6).tickFormat(tf));
    svg
      .append("g")
      .attr("class", "ab-axis")
      .attr("transform", "translate(" + m.l + ",0)")
      .call(tv ? d.axisLeft(y).tickValues(tv).tickFormat(tf) : d.axisLeft(y).ticks(6).tickFormat(tf));
    var top = Math.min(x.domain()[1], y.domain()[1]);
    svg
      .append("line")
      .attr("x1", x(0))
      .attr("y1", y(0))
      .attr("x2", x(top))
      .attr("y2", y(top))
      .attr("stroke", "#9aa9b5")
      .attr("stroke-dasharray", "4 4");
    svg
      .append("text")
      .attr("x", (m.l + W - m.r) / 2)
      .attr("y", H - 8)
      .attr("text-anchor", "middle")
      .attr("class", "ab-hd")
      .attr("fill", COLORS[L])
      .text(trunc(setName(L), 50) + " — " + METRICS[S.metric].short);
    svg
      .append("text")
      .attr("transform", "translate(14," + (m.t + H - m.b) / 2 + ") rotate(-90)")
      .attr("text-anchor", "middle")
      .attr("class", "ab-hd")
      .attr("fill", COLORS[R])
      .text(trunc(setName(R), 50) + " — " + METRICS[S.metric].short);
    svg
      .append("text")
      .attr("x", x(top) - 4)
      .attr("y", y(top) + 14)
      .attr("text-anchor", "end")
      .attr("class", "ab-sub")
      .text("equal");
    var pr = function (r) {
      var a = r.v[L],
        b = r.v[R];
      return (a.det + b.det) / Math.max(1, a.n + b.n);
    };
    var pts = svg
      .append("g")
      .selectAll("circle")
      .data(
        rows.slice().sort(function (a, b) {
          return pr(b) - pr(a);
        }),
      )
      .enter()
      .append("circle")
      .attr("cx", function (r) {
        return x(r.a);
      })
      .attr("cy", function (r) {
        return y(r.b);
      })
      .attr("r", function (r) {
        return 3 + 7 * Math.sqrt(pr(r));
      })
      .attr("fill", function (r) {
        return gcol(r.group);
      })
      .attr("fill-opacity", 0.7)
      .attr("stroke", function (r) {
        return r.hc ? "#b42318" : "#fff";
      })
      .attr("stroke-width", function (r) {
        return r.hc ? 2 : 1;
      });
    bindRow(pts);
    var lab = rows
      .slice()
      .sort(function (a, b) {
        return Math.abs(b.d) - Math.abs(a.d);
      })
      .slice(0, Math.min(8, rows.length));
    // Labels near the right edge flip to the left of their dot so they are never clipped.
    var flip = function (r) {
      return x(r.a) + 8 + Math.min(26, r.name.length) * 6 > W - m.r;
    };
    svg
      .append("g")
      .selectAll("text")
      .data(lab)
      .enter()
      .append("text")
      .attr("class", "ab-pl" + itCls())
      .attr("x", function (r) {
        return flip(r) ? x(r.a) - 8 : x(r.a) + 8;
      })
      .attr("text-anchor", function (r) {
        return flip(r) ? "end" : "start";
      })
      .attr("y", function (r) {
        return y(r.b) + 4;
      })
      .text(function (r) {
        return trunc(r.name, 26);
      });
    var used = {};
    rows.forEach(function (r) {
      used[r.group] = 1;
    });
    el.insertAdjacentHTML(
      "beforeend",
      '<div class="ab-legend">' +
        Object.keys(used)
          .sort(natural)
          .map(function (g) {
            return '<span><i class="ab-sw" style="background:' + gcol(g) + '"></i>' + esc(g) + "</span>";
          })
          .join("") +
        "</div>",
    );
    note(
      (S.sscale === "log" ? "Log-like axes (symmetric log, so zeros still show). " : "") +
        "Each dot is " +
        (S.level === "org" ? "an organism" : "a group") +
        " (colour = " +
        gbyLabel().toLowerCase() +
        ", size = share of samples where found, red ring = high consequence). Above the dashed line = higher in " +
        esc(setName(R)) +
        "; below = higher in " +
        esc(setName(L)) +
        ". Showing the top " +
        rows.length +
        " of " +
        C.rows.length +
        " by " +
        (S.ssort === "delta" ? "size of change" : "combined value") +
        " (Top box above); labels: the " +
        Math.min(8, rows.length) +
        " biggest changes among them.",
    );
  }
  function viewTable(el) {
    var L = S.left,
      R = S.right,
      rows = C.rows.slice(),
      k = S.tsort,
      dir = S.tdir;
    var get = {
      name: function (r) {
        return r.name;
      },
      group: function (r) {
        return r.group;
      },
      d: function (r) {
        return r.d;
      },
      absd: function (r) {
        return Math.abs(r.d);
      },
      lfc: function (r) {
        return r.lfc;
      },
      p: function (r) {
        return isFinite(r.p) ? r.p : 9;
      },
      q: function (r) {
        return isFinite(r.q) ? r.q : 9;
      },
    };
    var fn =
      get[k] ||
      (k.charAt(0) === "v"
        ? function (r) {
            return val(r, +k.slice(1));
          }
        : get.absd);
    rows.sort(function (a, b) {
      var x = fn(a),
        y = fn(b);
      return (typeof x === "string" ? natural(x, y) : x - y) * dir || natural(a.name, b.name);
    });
    var shown = S.tall ? rows : rows.slice(0, 300);
    var th = function (key, label, title) {
      return (
        '<th data-sort="' +
        key +
        '" title="' +
        esc(title || "Sort") +
        '" class="' +
        (k === key ? "on" : "") +
        '">' +
        label +
        (k === key ? (dir < 0 ? " ▾" : " ▴") : "") +
        "</th>"
      );
    };
    var h =
      '<div class="ab-twrap"><table class="ab-table"><thead><tr>' +
      th("name", S.level === "org" ? "Organism" : gbyLabel()) +
      (S.level === "org" ? th("group", esc(gbyLabel())) : "") +
      C.sets
        .map(function (se, i) {
          return (
            th("v" + i, sw(i) + LETTERS[i], se.label + " — " + METRICS[S.metric].label) + '<th class="ab-sm">found</th>'
          );
        })
        .join("") +
      th("d", "Δ " + LETTERS[R] + "−" + LETTERS[L], "Difference") +
      th("lfc", "Fold", "Fold change " + LETTERS[R] + " ÷ " + LETTERS[L]) +
      th(
        "p",
        "p",
        "Mann–Whitney U on the per-sample values, " + LETTERS[L] + " vs " + LETTERS[R] + " (needs ≥3 samples per set)",
      ) +
      th("q", "q", "Benjamini–Hochberg adjusted p") +
      "</tr></thead><tbody>";
    shown.forEach(function (r) {
      h +=
        '<tr data-key="' +
        esc(r.key) +
        '" class="' +
        (S.sel === r.key ? "sel" : "") +
        '"><td class="nm"><span class="' +
        (S.level === "org" ? "it" : "") +
        '">' +
        esc(r.name) +
        "</span>" +
        hcBadge(r) +
        "</td>" +
        (S.level === "org" ? '<td class="gp">' + esc(r.group) + "</td>" : "") +
        C.sets
          .map(function (se, i) {
            var x = r.v[i];
            return (
              '<td class="num">' +
              (x.det ? fv(x[S.metric]) : '<span class="z">0</span>') +
              '</td><td class="ab-sm">' +
              x.det +
              "/" +
              x.n +
              "</td>"
            );
          })
          .join("") +
        '<td class="num" style="color:' +
        (r.d > 0 ? COLORS[R] : r.d < 0 ? COLORS[L] : "#56707f") +
        '">' +
        fd(r.d) +
        '</td><td class="num">' +
        (r.onlyL
          ? '<span style="color:' + COLORS[L] + '">only ' + LETTERS[L] + "</span>"
          : r.onlyR
          ? '<span style="color:' + COLORS[R] + '">only ' + LETTERS[R] + "</span>"
          : fold(r.lfc)) +
        "</td>" +
        '<td class="num">' +
        fp(r.p) +
        '</td><td class="num' +
        (r.q < 0.05 ? " sig" : "") +
        '">' +
        fp(r.q) +
        "</td></tr>";
    });
    h += "</tbody></table></div>";
    if (rows.length > shown.length)
      h += '<div class="ab-hint">Showing ' + shown.length + " of " + rows.length + " — tick “Show every row”.</div>";
    el.innerHTML = h;
    note(
      "p / q are indicative only: samples from one site, run or batch are not independent replicates, and a set with fewer than 3 samples has no test. q &lt; 0.05 is shown in bold.",
    );
  }
  function csv() {
    var L = S.left,
      R = S.right;
    var q = function (v) {
      v = String(v == null ? "" : v);
      return /[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v;
    };
    var head = [S.level === "org" ? "organism" : "group", S.gby, "high_consequence"];
    C.sets.forEach(function (se, i) {
      head.push(
        LETTERS[i] + " " + se.label + " — " + METRICS[S.metric].short,
        LETTERS[i] + " samples detected",
        LETTERS[i] + " samples",
      );
    });
    head.push(
      "difference " + LETTERS[R] + "-" + LETTERS[L],
      "log2 fold " + LETTERS[R] + "/" + LETTERS[L],
      "p (Mann-Whitney)",
      "q (BH)",
    );
    var lines = [head.map(q).join(",")];
    C.rows.forEach(function (r) {
      var a = [r.name, r.group, r.hc ? "yes" : ""];
      C.sets.forEach(function (se, i) {
        a.push(val(r, i), r.v[i].det, r.v[i].n);
      });
      a.push(r.d, r.lfc, isFinite(r.p) ? r.p : "", isFinite(r.q) ? r.q : "");
      lines.push(a.map(q).join(","));
    });
    var name = "contrast_" + LETTERS[L] + "_vs_" + LETTERS[R] + ".csv",
      text = lines.join("\n");
    if (typeof _downloadText === "function") {
      _downloadText(text, name, "text/csv");
      return;
    }
    var blob = new Blob([text], { type: "text/csv" }),
      url = URL.createObjectURL(blob),
      aEl = document.createElement("a");
    aEl.href = url;
    aEl.download = name;
    document.body.appendChild(aEl);
    aEl.click();
    setTimeout(function () {
      URL.revokeObjectURL(url);
      aEl.remove();
    }, 500);
  }

  /* ── detail: per-sample strip plot ────────────────────────────────────── */
  function drawDetail(scroll) {
    var host = HOST && HOST.querySelector(".ab-detail-host");
    if (!host) return;
    var r =
      S.sel &&
      C &&
      C.rows.filter(function (x) {
        return x.key === S.sel;
      })[0];
    if (!r) {
      host.innerHTML = "";
      return;
    }
    Array.prototype.forEach.call(HOST.querySelectorAll(".ab-table tr.sel"), function (t) {
      t.classList.remove("sel");
    });
    var trs = HOST.querySelector(
      '.ab-table tr[data-key="' + (window.CSS && CSS.escape ? CSS.escape(r.key) : r.key) + '"]',
    );
    if (trs) trs.classList.add("sel");
    var isOrg = S.level === "org",
      vk = S.metric === "prev" ? "rel" : S.metric;
    var h =
      '<div class="ab-card ab-detail"><div class="ab-dh"><div><div class="ab-dt' +
      (isOrg ? " it" : "") +
      '"><span>' +
      esc(r.name) +
      "</span>" +
      hcBadge(r) +
      '</div><div class="ab-ds">' +
      esc(isOrg ? gbyLabel() + ": " + r.group : (r.members || 0) + " organisms") +
      " · " +
      esc(METRICS[vk].label) +
      '</div></div><div class="ab-db">' +
      '<button type="button" class="ab-btn" data-dact="close" title="Close">×</button></div></div>' +
      '<div class="ab-dstats">' +
      C.sets
        .map(function (se, i) {
          var x = r.v[i],
            arr = x.A[vk];
          return (
            "<div>" +
            sw(i) +
            "<b>" +
            esc(trunc(se.label, 30)) +
            "</b><span>" +
            (S.metric === "prev"
              ? fv(x.prev) + " of samples"
              : "mean " + fv(x[vk], vk) + " · median " + fv(median(arr), vk)) +
            " · found in " +
            x.det +
            "/" +
            x.n +
            "</span></div>"
          );
        })
        .join("") +
      '</div><div class="ab-strip"></div>' +
      '<div class="ab-hint">Each dot is one sample (' +
      (vk === "rel" ? "% of that sample’s reads" : METRICS[vk].short.replace(/^Mean /, "")) +
      "); faded dots = not detected; the bar marks the mean. Hover for the sample.</div></div>";
    host.innerHTML = h;
    var el = host.querySelector(".ab-strip"),
      d = d3(),
      W = width(el),
      lab = 190,
      rowH = 44,
      H = 26 + C.sets.length * rowH;
    var mx =
      d.max(C.sets, function (se, i) {
        return d.max(r.v[i].A[vk]);
      }) || 1;
    var x = scaleFn()
      .domain([0, mx])
      .range([lab, W - 20])
      .nice();
    var svg = svgFor(el, W, H);
    var tf =
      kind(vk) === "count"
        ? d.format("~s")
        : kind(vk) === "score"
        ? d.format("~g")
        : function (v) {
            return d.format("~g")(v) + "%";
          };
    svg
      .append("g")
      .attr("class", "ab-axis")
      .attr("transform", "translate(0," + (H - 20) + ")")
      .call(d.axisBottom(x).ticks(6).tickFormat(tf));
    C.sets.forEach(function (se, i) {
      var y0 = 6 + i * rowH,
        xv = r.v[i];
      svg
        .append("rect")
        .attr("x", lab)
        .attr("y", y0)
        .attr("width", W - 20 - lab)
        .attr("height", rowH - 8)
        .attr("fill", i % 2 ? "#fafbfc" : "#f5f8fa");
      svg
        .append("text")
        .attr("x", lab - 8)
        .attr("y", y0 + rowH / 2)
        .attr("text-anchor", "end")
        .attr("class", "ab-hd")
        .attr("fill", COLORS[i])
        .text(LETTERS[i] + " · " + trunc(se.label, 22))
        .append("title")
        .text(se.label);
      var data = se.samples.map(function (s, j) {
        return { s: s, v: xv.A[vk][j], det: xv.A.det[j], reads: xv.A.mreads[j], pct: xv.A.rel[j], tass: xv.A.tass[j] };
      });
      svg
        .append("g")
        .selectAll("circle")
        .data(data)
        .enter()
        .append("circle")
        .attr("cx", function (p) {
          return x(p.v);
        })
        .attr("cy", function (p) {
          var hsh = 0,
            k = String(p.s.sid);
          for (var c = 0; c < k.length; c++) hsh = (hsh * 31 + k.charCodeAt(c)) % 997;
          return y0 + 6 + (hsh / 997) * (rowH - 20);
        })
        .attr("r", 4)
        .attr("fill", COLORS[i])
        .attr("fill-opacity", function (p) {
          return p.det ? 0.75 : 0.18;
        })
        .attr("stroke", "#fff")
        .on("mousemove", function (ev, p) {
          tip(
            ev,
            "<b>" +
              esc(p.s.name) +
              "</b>" +
              (sampleSub(p.s) ? "<br>" + esc(sampleSub(p.s)) : "") +
              "<br>" +
              (p.det
                ? fmtN(p.reads) +
                  " reads aligned · " +
                  d.format(".3~g")(p.pct) +
                  "% of the sample · TASS " +
                  d.format(".1f")(p.tass)
                : "not detected"),
          );
        })
        .on("mouseleave", untip);
      var mean = xv[vk];
      svg
        .append("line")
        .attr("x1", x(mean))
        .attr("x2", x(mean))
        .attr("y1", y0 + 2)
        .attr("y2", y0 + rowH - 10)
        .attr("stroke", COLORS[i])
        .attr("stroke-width", 3);
    });
    if (scroll) host.scrollIntoView({ behavior: "smooth", block: "nearest" });
  }

  /* ── events ───────────────────────────────────────────────────────────── */
  function rerender() {
    untip();
    render(HOST);
  }
  function wire(host) {
    if (host._ab) return;
    host._ab = 1;
    host.addEventListener("click", function (e) {
      var t = e.target,
        b;
      if ((b = t.closest("[data-view]"))) {
        S.view = b.getAttribute("data-view");
        save();
        return rerender();
      }
      if ((b = t.closest("[data-level]"))) {
        S.level = b.getAttribute("data-level");
        S.sel = "";
        save();
        return rerender();
      }
      if ((b = t.closest("[data-preset]"))) {
        var f = (host.querySelector('select[data-k="pfield"]') || {}).value;
        S.sets = preset(b.getAttribute("data-preset"), f);
        S.left = 0;
        S.right = 1;
        S.q = {};
        save();
        return rerender();
      }
      if ((b = t.closest("[data-add]"))) {
        var used = {};
        S.sets.forEach(function (st) {
          (st.vals || []).forEach(function (v) {
            used[st.field + ":" + v] = 1;
          });
        });
        var f2 = S.sets[S.sets.length - 1].field;
        var next =
          f2 === "range" || f2 === "samples"
            ? []
            : realValues(f2)
                .filter(function (x) {
                  return !used[f2 + ":" + x.v];
                })
                .slice(0, 1)
                .map(function (x) {
                  return x.v;
                });
        S.sets.push({ field: f2, vals: next, label: "" });
        save();
        return rerender();
      }
      if ((b = t.closest("[data-rm]"))) {
        S.sets.splice(+b.getAttribute("data-rm"), 1);
        if (S.left >= S.sets.length) S.left = 0;
        if (S.right >= S.sets.length || S.right === S.left) S.right = S.left === 0 ? 1 : 0;
        S.q = {};
        save();
        return rerender();
      }
      if ((b = t.closest("[data-all]")) || (b = t.closest("[data-none]"))) {
        var k = +(b.getAttribute("data-all") || b.getAttribute("data-none")),
          on = b.hasAttribute("data-all");
        var st = S.sets[k],
          cur = {};
        var vis = Array.prototype.map.call(
          host.querySelectorAll('.ab-vals[data-i="' + k + '"] input[data-v]'),
          function (x) {
            return x.getAttribute("data-v");
          },
        );
        (st.vals || []).forEach(function (v) {
          cur[v] = 1;
        });
        vis.forEach(function (v) {
          if (on) cur[v] = 1;
          else delete cur[v];
        });
        st.vals = Object.keys(cur);
        save();
        return rerender();
      }
      if ((b = t.closest("[data-sort]"))) {
        var sk = b.getAttribute("data-sort");
        if (S.tsort === sk) S.tdir = -S.tdir;
        else {
          S.tsort = sk;
          S.tdir = sk === "name" || sk === "group" || sk === "p" || sk === "q" ? 1 : -1;
        }
        save();
        drawView();
        return;
      }
      if (t.closest("[data-csv]")) {
        csv();
        return;
      }
      if ((b = t.closest("[data-dact]"))) {
        if (b.getAttribute("data-dact") === "close") {
          S.sel = "";
          drawDetail();
        }
        return;
      }
      var tr = t.closest(".ab-table tbody tr");
      if (tr) {
        S.sel = tr.getAttribute("data-key");
        drawDetail(true);
      }
    });
    host.addEventListener("change", function (e) {
      var t = e.target,
        k = t.getAttribute("data-k"),
        i = t.getAttribute("data-i");
      if (t.matches("input[data-v]")) {
        var st = S.sets[+i],
          v = t.getAttribute("data-v"),
          vals = (st.vals || []).filter(function (x) {
            return x !== v;
          });
        if (t.checked) vals.push(v);
        st.vals = vals;
        save();
        return rerender();
      }
      if (t.classList.contains("ab-lab")) {
        S.sets[+i].label = t.value.trim();
        save();
        return rerender();
      }
      if (k === "field") {
        var nf = t.value,
          st2 = S.sets[+i];
        if (nf === "range") {
          var days = pool()
            .map(function (s) {
              return s.day;
            })
            .filter(Boolean)
            .sort();
          S.sets[+i] = { field: "range", from: days[0] || "", to: days[days.length - 1] || "", label: st2.label };
        } else S.sets[+i] = { field: nf, vals: [], label: st2.label };
        S.q[i] = "";
        save();
        return rerender();
      }
      if (k === "from" || k === "to") {
        S.sets[+i][k] = t.value;
        save();
        return rerender();
      }
      if (k === "pfield") {
        S.pfield = t.value;
        save();
        return;
      }
      if (k === "exclCtl" || k === "tall") {
        S[k] = t.checked;
        save();
        return k === "tall" ? drawView() : rerender();
      }
      if (k === "topN" || k === "rk") {
        var n = parseInt(t.value, 10);
        if (isFinite(n)) S[k] = Math.max(3, Math.min(k === "rk" ? 16 : 200, n));
        save();
        return rerender();
      }
      if (k === "left" || k === "right") {
        var n2 = +t.value,
          other = k === "left" ? "right" : "left";
        if (S[other] === n2) S[other] = S[k];
        S[k] = n2;
        save();
        return rerender();
      }
      if (k === "gby") {
        S.gby = t.value;
        S.tgroup = "";
        S.sel = S.level === "grp" ? "" : S.sel;
        save();
        return rerender();
      }
      if (k && Object.prototype.hasOwnProperty.call(S, k)) {
        S[k] = t.value;
        save();
        return rerender();
      }
    });
    host.addEventListener("input", function (e) {
      var t = e.target;
      if (!t.classList.contains("ab-q")) return;
      var i = t.getAttribute("data-i");
      S.q[i] = t.value;
      var box = host.querySelector('.ab-vals[data-i="' + i + '"]');
      var tmp = document.createElement("div");
      tmp.innerHTML = setBody(S.sets[+i], +i);
      var nb = tmp.querySelector(".ab-vals"),
        na = tmp.querySelector(".ab-vacts");
      if (box && nb) box.innerHTML = nb.innerHTML;
      var oa = box && box.previousElementSibling;
      if (oa && na && oa.classList.contains("ab-vacts")) oa.innerHTML = na.innerHTML;
    });
    var ro = window.ResizeObserver
      ? new ResizeObserver(function () {
          var w = host.clientWidth;
          if (!w || Math.abs(w - (host._abW || 0)) < 24) return;
          host._abW = w;
          clearTimeout(host._abT);
          host._abT = setTimeout(function () {
            if (host.offsetParent && C) {
              drawView();
              drawDetail();
            }
          }, 150);
        })
      : null;
    if (ro) ro.observe(host);
  }

  // The shared Group by selection feeds the "Group by selection" set field.
  try {
    if (window.metaGrouping && window.metaGrouping.onChange) {
      window.metaGrouping.onChange(function () {
        if (typeof _TAB_DIRTY !== "undefined") _TAB_DIRTY.contrast = true;
        if (typeof activeTab !== "undefined" && activeTab === "contrast" && HOST) rerender();
      });
    }
  } catch (e) {
    /* grouping engine not loaded */
  }

  window.ttContrast = {
    render: render,
    state: S,
    compute: function () {
      return C;
    },
    setSets: function (sets) {
      S.sets = sets;
      save();
      if (HOST) rerender();
    },
  };
})();
