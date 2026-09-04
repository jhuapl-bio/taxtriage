/* ═══════════════════════════════════════════════════════════════════════════
       -  §  IN-SILICO CROSS-REFERENCE  (Detections table ⚗ badge)
       -     Puts the in-silico dilution series next to a REAL detection.
       -
       -     A detection row says "organism X had N reads in sample S". The
       -     subsampling suite says "at depth D this organism yields R reads and
       -     TASS T". This module places the real sample ON that dilution axis at
       -     its own sequencing depth and answers the question the suite exists to
       -     answer: given how deep this sample was sequenced, is the organism
       -     showing up as strongly as the dilution series says it should?
       -
       -     Everything is depth-relative. If the series runs 5k / 8k / 10k /
       -     15k / 20k / 30k and the sample carries 20k reads, the sample is drawn
       -     at 20k and compared against the series value interpolated at 20k —
       -     never against the deepest dataset.
       -
       -     Public API (all no-ops when the run had no --sim_subsample):
       -       insilicoCompareFor(row)  -> comparison object | null
       -       insilicoBadgeHTML(row)   -> ⚗ chip HTML for the organism cell
       -       insilicoTipHTML(row)     -> hover card (mini plot + key numbers)
       -       openInsilicoCompare(row) -> full modal (reads + TASS panels)
       -       wireInsilicoBadges(root, keyToRow)  -> hover/click wiring
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  "use strict";

  var ACCENT = "#4527a0";
  var GOOD = "#2e7d32";
  var WARN = "#ef6c00";
  var BAD = "#c62828";
  var MUTED = "#777";
  var GRID = "#e9e5f3";
  var AXIS = "#bdb6d4";
  var REAL = "#0277bd";   // the real sample always reads as blue

  function suite() {
    var s = typeof INSILICO_SUITE !== "undefined" ? INSILICO_SUITE : (window.HEATMAP_BOOT || {}).insilico_suite;
    return s && s.enabled && s.groups && s.groups.length ? s : null;
  }

  function esc(s) {
    return String(s == null ? "" : s).replace(/[&<>"]/g, function (m) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[m];
    });
  }

  function kfmt(n) {
    if (n == null || isNaN(n)) return "—";
    n = +n;
    var a = Math.abs(n);
    if (a >= 1e6) return (n / 1e6).toFixed(a % 1e6 ? 1 : 0) + "M";
    if (a >= 1000) return (n / 1000).toFixed(a % 1000 ? 1 : 0) + "k";
    return String(Math.round(n * 10) / 10);
  }

  function icomma(n) {
    if (n == null || isNaN(n)) return "—";
    return Math.round(n).toLocaleString();
  }

  function pct(x, d) {
    if (x == null || isNaN(x) || !isFinite(x)) return "—";
    return (x * 100).toFixed(d == null ? 0 : d) + "%";
  }

  // ── locate the row's organism in the suite ────────────────────────────────
  // Prefer the group seeded from this very sample; otherwise any group that
  // simulated the same organism. Matching is by taxid, falling back to an exact
  // (case-insensitive) name match for rows whose taxid was not carried through.
  // A detection can sit BELOW the level the series was built at: the suite
  // collapses to Species (or Strain), so a Strain-level detection of a monkeypox
  // isolate has no series of its own while its species does. Matching therefore
  // walks the row's lineage — the row's own taxon first, then its species, then
  // its genus (rolling the series up to genus, since a genus is usually several
  // simulated organisms). Every hit records which rung matched so the UI can say
  // so plainly instead of implying a like-for-like comparison.
  function lc(v) {
    return String(v == null ? "" : v).trim().toLowerCase();
  }

  function genusRollup(group) {
    if (group.__genusRollup) return group.__genusRollup;
    if (typeof window.insilicoRollupOrganisms !== "function") return [];  // 45 not loaded; don't cache
    var out = window.insilicoRollupOrganisms(group, "Genus");
    group.__genusRollup = out;
    return out;
  }

  function matchOrganism(sample, row) {
    var s = suite();
    if (!s || !row) return [];
    var tid = String(row["Taxonomic ID #"] == null ? "" : row["Taxonomic ID #"]).trim();
    var nm = lc(row["Detected Organism"]);
    var spName = lc(row["Species Name"]);
    var genName = lc(row["Genus Name"]);
    var rowLevel = row["Level"] || "Strain";
    var hits = [];
    // Only climb: a Strain may borrow its species' or genus' series, a Species
    // may borrow its genus', a Genus borrows nothing. Without this gate a
    // Genus-level row would match on its own Species Name field (which the
    // flattener fills in for every level) and silently show a species series.
    var canSpecies = rowLevel === "Strain";
    var canGenus = rowLevel === "Strain" || rowLevel === "Species";

    s.groups.forEach(function (g) {
      var seriesLevel = g.level || "Strain";
      (g.organisms || []).forEach(function (o) {
        var how = null;
        if (tid && String(o.taxid) === tid) how = "exact";
        else if (nm && lc(o.name) === nm) how = "name";
        // The row's species matched a series organism (or that organism's own
        // species label) — the usual Strain -> Species case.
        else if (canSpecies && spName && (lc(o.species) === spName || lc(o.name) === spName)) how = "species";
        if (how) {
          hits.push({
            group: g, org: o, own: g.parent === sample,
            how: how, seriesLevel: seriesLevel, rowLevel: rowLevel,
            matchLabel: how === "species" ? "Species" : seriesLevel,
          });
        }
      });
      // Genus rung: roll the whole group up so the series covers every member of
      // the genus, not just whichever one happens to share the name. A
      // Genus-level detection lands here as its natural match, not as a climb.
      if (genName && (canGenus || rowLevel === "Genus")) {
        genusRollup(g).forEach(function (o) {
          if (lc(o.name) !== genName) return;
          hits.push({
            group: g, org: o, own: g.parent === sample,
            // For a Genus-level detection this IS its own level — the series is
            // merely assembled from the members, not climbed to.
            how: rowLevel === "Genus" ? "genus_native" : "genus",
            seriesLevel: seriesLevel, rowLevel: rowLevel, matchLabel: "Genus",
          });
        });
      }
    });

    var rank = { exact: 0, genus_native: 0, name: 1, species: 2, genus: 3 };
    hits.sort(function (a, b) {
      if (a.own !== b.own) return a.own ? -1 : 1;                 // this sample's own series first
      if (rank[a.how] !== rank[b.how]) return rank[a.how] - rank[b.how];  // closest rung wins
      return (b.org.series || []).length - (a.org.series || []).length;
    });
    // Keep at most one hit per (group, rung) so the selector is not flooded.
    var seen = {};
    return hits.filter(function (h) {
      var k = h.group.parent + "|" + h.group.platform + "|" + h.how + "|" + h.org.taxid;
      if (seen[k]) return false;
      seen[k] = 1;
      return true;
    });
  }

  // ── interpolation along the series ────────────────────────────────────────
  // Piecewise-linear in read count. Below the shallowest dataset reads are
  // scaled through the origin (reads grow with depth); above the deepest the
  // last segment's slope is continued. `kind:"clamp"` (used for TASS, which is
  // not linear in depth) holds the endpoint value instead and flags it.
  function interpAt(series, x, key, kind) {
    if (!series || !series.length) return null;
    var pts = series
      .map(function (s) { return { x: +s.count, y: +s[key] }; })
      .filter(function (p) { return isFinite(p.x) && isFinite(p.y); })
      .sort(function (a, b) { return a.x - b.x; });
    if (!pts.length) return null;
    var lo = pts[0], hi = pts[pts.length - 1];
    if (pts.length === 1) {
      if (kind === "clamp") return { v: lo.y, extrap: x !== lo.x };
      return { v: lo.x > 0 ? lo.y * (x / lo.x) : lo.y, extrap: x !== lo.x };
    }
    if (x <= lo.x) {
      if (kind === "clamp") return { v: lo.y, extrap: x < lo.x };
      return { v: lo.x > 0 ? lo.y * (x / lo.x) : lo.y, extrap: x < lo.x };
    }
    if (x >= hi.x) {
      if (kind === "clamp") return { v: hi.y, extrap: x > hi.x };
      var p = pts[pts.length - 2];
      var m = (hi.y - p.y) / ((hi.x - p.x) || 1);
      return { v: Math.max(0, hi.y + m * (x - hi.x)), extrap: x > hi.x };
    }
    for (var i = 0; i < pts.length - 1; i++) {
      var a = pts[i], b = pts[i + 1];
      if (x >= a.x && x <= b.x) {
        var t = (x - a.x) / ((b.x - a.x) || 1);
        return { v: a.y + t * (b.y - a.y), extrap: false };
      }
    }
    return { v: hi.y, extrap: true };
  }

  // ── build the comparison ──────────────────────────────────────────────────
  // When the series is a rung above the detection, the row's own read count is
  // only part of the picture: a strain's reads are a subset of its species'. If
  // the report also holds the matching Species/Genus row for this sample (the
  // rollup levels are all present in DATA), use that row's reads and TASS for
  // the comparison and keep the original row's numbers alongside it.
  function rollupRowFor(row, level) {
    if (typeof DATA === "undefined" || !DATA || !DATA.length) return null;
    var sample = row["Specimen ID"];
    var want = level === "Genus" ? lc(row["Genus Name"]) : lc(row["Species Name"]);
    if (!want) return null;
    var best = null;
    for (var i = 0; i < DATA.length; i++) {
      var r = DATA[i];
      if (r["Specimen ID"] !== sample) continue;
      if ((r["Level"] || "Strain") !== level) continue;
      var key = level === "Genus" ? lc(r["Genus Name"]) : lc(r["Species Name"]);
      if (key !== want) continue;
      if (!best || (+r["# Reads Aligned"] || 0) > (+best["# Reads Aligned"] || 0)) best = r;
    }
    return best;
  }

  function insilicoCompareFor(row, groupIdx) {
    if (!row || !suite()) return null;
    var sample = row["Specimen ID"];
    var hits = matchOrganism(sample, row);
    if (!hits.length) return null;
    var hit = hits[Math.min(groupIdx || 0, hits.length - 1)];
    var g = hit.group, o = hit.org;
    var rolledUp = hit.how === "species" || hit.how === "genus";
    var cmpLevel = hit.how === "genus" || hit.how === "genus_native" ? "Genus" : "Species";
    var cmpRow = rolledUp ? rollupRowFor(row, cmpLevel) : null;
    var srcRow = cmpRow || row;

    var meta = (typeof SAMPLE_META !== "undefined" && SAMPLE_META[sample]) || {};
    var platform = String(meta.platform || "").toUpperCase();
    // The series counts read PAIRS for paired-end groups, while the report's
    // read columns count each mate. Halve the real sample so both sides of the
    // comparison are in the same unit.
    var paired = (g.read_unit || "reads") === "read pairs";
    var rpr = paired && (platform === "ILLUMINA" || platform === "") ? 2 : 1;
    var unitMismatch = paired && platform && platform !== "ILLUMINA";

    var totalReads = +meta.total_reads || 0;
    var depth = totalReads > 0 ? totalReads / rpr : 0;
    var orgReads = (+srcRow["# Reads Aligned"] || 0) / rpr;
    var tass = +srcRow["TASS Score"] || 0;
    var rowReads = (+row["# Reads Aligned"] || 0) / rpr;
    var rowTass = +row["TASS Score"] || 0;

    var thr = +(suite().params || {}).detection_threshold || 0;
    var predObs = depth > 0 ? interpAt(o.series, depth, "observed_reads") : null;
    var predTass = depth > 0 ? interpAt(o.series, depth, "tass", "clamp") : null;
    // The "expected" model is compositional and therefore exact at any depth:
    // this organism's share of the simulated pool times the sample's depth.
    var expReads = depth > 0 ? o.expected_fraction * depth : null;

    var vsSeries = predObs && predObs.v > 0 ? orgReads / predObs.v : null;
    var vsExpected = expReads > 0 ? orgReads / expReads : null;

    return {
      row: row, sample: sample, group: g, org: o,
      // How the row was matched to the series, and at what level the comparison
      // is therefore being made.
      how: hit.how, rolledUp: rolledUp, genusSeries: hit.how === "genus" || hit.how === "genus_native",
      rowLevel: row["Level"] || "Strain",
      seriesLevel: hit.seriesLevel, compareLevel: rolledUp ? cmpLevel : (row["Level"] || "Strain"),
      compareName: rolledUp ? (cmpLevel === "Genus" ? row["Genus Name"] : row["Species Name"]) : row["Detected Organism"],
      usedRollupRow: !!cmpRow, rowReads: rowReads, rowTass: rowTass,
      rowShare: orgReads > 0 ? rowReads / orgReads : null,
      groupIdx: hits.indexOf(hit), nGroups: hits.length, matchedByName: hit.how === "name", ownSeries: hit.own,
      unit: g.read_unit || "reads", paired: paired, rpr: rpr, unitMismatch: unitMismatch,
      platform: platform, totalReads: totalReads, depth: depth,
      orgReads: orgReads, tass: tass, threshold: thr,
      expectedReads: expReads,
      predictedReads: predObs ? predObs.v : null, predictedExtrap: !!(predObs && predObs.extrap),
      predictedTass: predTass ? predTass.v : null, predictedTassExtrap: !!(predTass && predTass.extrap),
      vsSeries: vsSeries, vsExpected: vsExpected,
      lod: o.lod_count == null ? null : o.lod_count,
      aboveLod: o.lod_count != null && depth >= o.lod_count,
      lodMargin: o.lod_count ? depth / o.lod_count : null,
    };
  }

  // Verdict drives the badge colour and the hover headline.
  function verdict(c) {
    if (!c || !c.depth) return { key: "nodepth", color: MUTED, label: "no sample depth recorded" };
    if (c.lod != null && !c.aboveLod) {
      return { key: "belowlod", color: BAD,
               label: "sequenced below this organism's limit of detection (" + kfmt(c.lod) + " " + c.unit + ")" };
    }
    if (c.lod == null) {
      return { key: "nolod", color: WARN, label: "never detected anywhere in the dilution series" };
    }
    if (c.vsSeries == null) return { key: "unknown", color: MUTED, label: "series has no value at this depth" };
    if (c.rolledUp && !c.usedRollupRow && c.vsSeries < 0.5) {
      // Comparing a subset against its parent level — being below the series is
      // the expected arithmetic, not a finding.
      return { key: "subset", color: MUTED,
               label: "below the " + c.compareLevel.toLowerCase() + " series, as expected for a single " +
                      c.rowLevel.toLowerCase() };
    }
    if (c.vsSeries >= 0.5 && c.vsSeries <= 2) {
      return { key: "match", color: GOOD, label: "in line with the dilution series at this depth" };
    }
    if (c.vsSeries < 0.5) {
      return { key: "under", color: WARN, label: "under-recovered vs the dilution series at this depth" };
    }
    return { key: "over", color: ACCENT, label: "over-recovered vs the dilution series at this depth" };
  }

  // ── scales ────────────────────────────────────────────────────────────────
  function makeScale(min, max, log, p0, p1) {
    var t = log
      ? function (v) { return Math.log10(Math.max(0, v) + 1); }
      : function (v) { return v; };
    var a = t(min), b = t(max);
    if (!(b > a)) b = a + 1;
    return {
      log: log, min: min, max: max,
      f: function (v) {
        var p = (t(v) - a) / (b - a);
        return p0 + Math.max(-0.05, Math.min(1.05, p)) * (p1 - p0);
      },
    };
  }

  function ticksFor(min, max, log, n) {
    var out = [];
    if (log) {
      var lo = Math.floor(Math.log10(Math.max(min, 1)));
      var hi = Math.ceil(Math.log10(Math.max(max, 10)));
      for (var e = lo; e <= hi; e++) {
        var v = Math.pow(10, e);
        if (v >= min * 0.5 && v <= max * 1.5) out.push(v);
      }
      if (out.length < 2) out = [min, max];
      return out;
    }
    n = n || 4;
    for (var i = 0; i <= n; i++) out.push(min + ((max - min) * i) / n);
    return out;
  }

  // Log axes once the span is wide enough that a linear one would flatten the
  // shallow end into the baseline.
  function wantLog(min, max) {
    return min > 0 ? max / min >= 50 : max >= 5000;
  }

  // ── the plot ──────────────────────────────────────────────────────────────
  // One panel. `key` is the series field ("observed_reads" or "tass"); the real
  // sample is drawn as a diamond at its own depth with a drop line, so the
  // comparison is read off the x position, not off a bar order.
  function panel(c, key, opts) {
    opts = opts || {};
    var W = opts.W || 380, plotH = opts.plotH || 132;
    var padL = opts.padL || 52, padR = 14, padT = 12, labelH = 16, titleH = 13;
    var H = padT + plotH + labelH + titleH + 6;
    var yb = padT + plotH;
    var series = (c.org.series || []).slice().sort(function (a, b) { return a.count - b.count; });
    var isReads = key === "observed_reads";

    var xs = series.map(function (s) { return s.count; }).concat([c.depth]);
    var xmin = Math.min.apply(null, xs.filter(function (v) { return v > 0; }));
    var xmax = Math.max.apply(null, xs);
    var xlog = wantLog(xmin, xmax);
    var X = makeScale(xlog ? xmin / 1.6 : 0, xmax * 1.12, xlog, padL, W - padR);

    // Reads scale to the data; TASS is a fixed 0-100 score, so its axis is
    // pinned to that range — a per-plot autoscale made two TASS panels
    // incomparable and exaggerated small differences.
    var ymax, ylog, Y, yTicks;
    if (isReads) {
      var yvals = series.map(function (s) { return +s[key] || 0; })
        .concat(series.map(function (s) { return +s.expected_reads || 0; }));
      yvals.push(c.orgReads);
      if (c.expectedReads != null) yvals.push(c.expectedReads);
      ymax = Math.max.apply(null, yvals.concat([1]));
      ylog = wantLog(1, ymax);
      Y = makeScale(0, ymax * 1.1, ylog, yb, padT);
      yTicks = ticksFor(0, ymax * 1.1, ylog, 4);
    } else {
      ymax = 100;
      ylog = false;
      Y = makeScale(0, 100, false, yb, padT);
      yTicks = [0, 25, 50, 75, 100];
    }

    var s = '<svg viewBox="0 0 ' + W + " " + H + '" width="' + W + '" height="' + H +
            '" role="img" style="display:block">';

    // grid + y labels
    yTicks.forEach(function (v) {
      var y = Y.f(v);
      if (y < padT - 2 || y > yb + 2) return;
      s += '<line x1="' + padL + '" y1="' + y.toFixed(1) + '" x2="' + (W - padR) + '" y2="' + y.toFixed(1) +
           '" stroke="' + GRID + '" stroke-width="1"/>' +
           '<text x="' + (padL - 6) + '" y="' + (y + 3).toFixed(1) + '" text-anchor="end" font-size="9" fill="' +
           MUTED + '">' + esc(isReads ? kfmt(v) : String(Math.round(v))) + "</text>";
    });
    s += '<line x1="' + padL + '" y1="' + yb + '" x2="' + (W - padR) + '" y2="' + yb +
         '" stroke="' + AXIS + '" stroke-width="1.2"/>' +
         '<line x1="' + padL + '" y1="' + padT + '" x2="' + padL + '" y2="' + yb +
         '" stroke="' + AXIS + '" stroke-width="1.2"/>' +
         '<text transform="translate(11,' + (padT + plotH / 2) + ') rotate(-90)" text-anchor="middle" font-size="9" fill="' +
         MUTED + '">' + esc(isReads ? c.unit + (ylog ? " (log)" : "") : "TASS") + "</text>";

    // x ticks
    ticksFor(X.min, X.max, xlog, 4).forEach(function (v) {
      var x = X.f(v);
      if (x < padL - 2 || x > W - padR + 2) return;
      s += '<line x1="' + x.toFixed(1) + '" y1="' + yb + '" x2="' + x.toFixed(1) + '" y2="' + (yb + 3) +
           '" stroke="' + AXIS + '" stroke-width="1"/>' +
           '<text x="' + x.toFixed(1) + '" y="' + (yb + 13) + '" text-anchor="middle" font-size="9" fill="#666">' +
           esc(kfmt(v)) + "</text>";
    });
    s += '<text x="' + ((padL + W - padR) / 2) + '" y="' + (yb + labelH + 10) +
         '" text-anchor="middle" font-size="8.5" fill="#999">sequencing depth (' + esc(c.unit) + ")</text>";

    // LoD rule
    if (c.lod != null) {
      var xl = X.f(c.lod);
      s += '<line x1="' + xl.toFixed(1) + '" y1="' + padT + '" x2="' + xl.toFixed(1) + '" y2="' + yb +
           '" stroke="' + ACCENT + '" stroke-width="1" stroke-dasharray="3 3" opacity=".5"/>' +
           // flip the label to the left of the rule when the rule sits near the
           // right edge, so it never lands on the last series marker
           '<text x="' + (xl + (xl > W - padR - 40 ? -3 : 3)).toFixed(1) + '" y="' + (padT + 8) +
           '" text-anchor="' + (xl > W - padR - 40 ? "end" : "start") + '" font-size="8" fill="' + ACCENT +
           '" opacity=".85">LoD</text>';
    }

    // expected (compositional) line — reads panel only
    if (isReads) {
      var ep = series.map(function (p) {
        return X.f(p.count).toFixed(1) + "," + Y.f(c.org.expected_fraction * p.count).toFixed(1);
      });
      s += '<polyline fill="none" stroke="' + ACCENT + '" stroke-width="1.2" stroke-dasharray="4 3" opacity=".75" points="' +
           ep.join(" ") + '"/>';
    }
    // detection threshold — TASS panel only
    if (!isReads && c.threshold) {
      var yt = Y.f(c.threshold);
      s += '<line x1="' + padL + '" y1="' + yt.toFixed(1) + '" x2="' + (W - padR) + '" y2="' + yt.toFixed(1) +
           '" stroke="' + BAD + '" stroke-width="1" stroke-dasharray="4 3" opacity=".7"/>' +
           '<text x="' + (W - padR - 2) + '" y="' + (yt - 3).toFixed(1) + '" text-anchor="end" font-size="8" fill="' +
           BAD + '">cutoff ' + esc(String(c.threshold)) + "</text>";
    }

    // observed series
    var op = series.map(function (p) {
      return X.f(p.count).toFixed(1) + "," + Y.f(+p[key] || 0).toFixed(1);
    });
    s += '<polyline fill="none" stroke="' + GOOD + '" stroke-width="1.8" stroke-linejoin="round" points="' +
         op.join(" ") + '"/>';
    series.forEach(function (p) {
      s += '<circle cx="' + X.f(p.count).toFixed(1) + '" cy="' + Y.f(+p[key] || 0).toFixed(1) +
           '" r="3" fill="' + (p.detected ? GOOD : "#bbb") + '" stroke="' + GOOD + '" stroke-width="1"/>';
    });

    // the real sample
    if (c.depth > 0) {
      var rx = X.f(c.depth);
      var ry = Y.f(isReads ? c.orgReads : c.tass);
      s += '<line x1="' + rx.toFixed(1) + '" y1="' + padT + '" x2="' + rx.toFixed(1) + '" y2="' + yb +
           '" stroke="' + REAL + '" stroke-width="1" stroke-dasharray="2 2" opacity=".55"/>';
      // predicted value at the sample's depth, as a hollow marker on the series
      var pv = isReads ? c.predictedReads : c.predictedTass;
      if (pv != null) {
        // Sized to clear the diamond so the two stay distinguishable when the
        // sample lands right on the series value (the "in line" case).
        s += '<circle cx="' + rx.toFixed(1) + '" cy="' + Y.f(pv).toFixed(1) +
             '" r="7" fill="none" stroke="' + GOOD + '" stroke-width="1.8"/>';
      }
      var d = 5;
      s += '<polygon points="' + [
             rx.toFixed(1) + "," + (ry - d).toFixed(1),
             (rx + d).toFixed(1) + "," + ry.toFixed(1),
             rx.toFixed(1) + "," + (ry + d).toFixed(1),
             (rx - d).toFixed(1) + "," + ry.toFixed(1),
           ].join(" ") + '" fill="' + REAL + '" stroke="#fff" stroke-width="1"/>';
    }
    s += "</svg>";
    return s;
  }

  // Swatch kinds: solid = filled square, hdash = horizontal dashed rule,
  // ring = hollow circle, diamond = the real sample's marker, vdash = the
  // vertical dashed rule the diamond and ring sit on.
  function swatch(kind, color) {
    if (kind === "hdash") {
      return '<span style="display:inline-block;width:13px;border-top:1.6px dashed ' + color +
             ';vertical-align:3px;margin-right:4px"></span>';
    }
    if (kind === "ring") {
      return '<span style="display:inline-block;width:9px;height:9px;border-radius:50%;border:1.8px solid ' +
             color + ';background:#fff;margin-right:4px;vertical-align:-1px"></span>';
    }
    if (kind === "diamond") {
      return '<span style="display:inline-block;width:9px;height:9px;background:' + color +
             ';transform:rotate(45deg);margin-right:5px;vertical-align:-1px"></span>';
    }
    if (kind === "vdash") {
      return '<span style="display:inline-block;width:0;height:11px;border-left:1.6px dashed ' + color +
             ';margin:0 8px 0 4px;vertical-align:-2px"></span>';
    }
    return '<span style="display:inline-block;width:9px;height:9px;border-radius:2px;background:' + color +
           ';margin-right:4px;vertical-align:-1px"></span>';
  }

  function legendHTML(isReads, c) {
    var unit = (c && c.unit) || "reads";
    var depth = c && c.depth ? icomma(c.depth) + " " + unit : "its own depth";
    var it = [[isReads ? "series observed" : "series TASS", GOOD, "solid"]];
    if (isReads) it.push(["expected from pool share", ACCENT, "hdash"]);
    else it.push(["detection cutoff", BAD, "hdash"]);
    it.push([
      "this sample, plotted at its sequencing depth (" + depth + ")", REAL, "diamond",
    ]);
    it.push([
      "what the series gives at that same depth — the value the diamond is measured against",
      GOOD, "ring",
    ]);
    it.push(["the sample's depth (both markers sit on this rule)", REAL, "vdash"]);
    if (c && c.lod != null) it.push(["limit of detection", ACCENT, "vdash"]);
    return (
      '<div style="font-size:.72em;color:' + MUTED + ';margin-top:4px;display:flex;flex-direction:column;gap:.18em">' +
      it.map(function (o) {
        return '<span style="padding-left:17px;text-indent:-17px">' + swatch(o[2], o[1]) + esc(o[0]) + "</span>";
      }).join("") + "</div>"
    );
  }

  // ── numbers table shared by the hover card and the modal ──────────────────
  function statRows(c, dark) {
    var lab = dark ? "opacity:.72" : "color:" + MUTED;
    var v = verdict(c);
    var rows = [
      ["Sample depth", icomma(c.depth) + " " + c.unit + (c.rpr === 2 ? " (" + icomma(c.totalReads) + " reads)" : "")],
      [(c.rolledUp ? c.compareLevel + " reads here" : "Organism reads here"), icomma(c.orgReads) + " " + c.unit],
      ["Series at this depth", c.predictedReads == null ? "—" : icomma(c.predictedReads) + " " + c.unit +
        (c.predictedExtrap ? " (extrapolated)" : "")],
      ["Expected from pool share", c.expectedReads == null ? "—" : icomma(c.expectedReads) + " " + c.unit +
        " (" + pct(c.org.expected_fraction, 1) + " of pool)"],
      ["Sample vs series", c.vsSeries == null ? "—" :
        '<b style="color:' + v.color + '">' + (c.vsSeries).toFixed(2) + "×</b>"],
      ["Sample vs expected", c.vsExpected == null ? "—" : c.vsExpected.toFixed(2) + "×"],
      ["TASS here", c.tass.toFixed(1) + (c.predictedTass == null ? "" :
        "  vs series " + c.predictedTass.toFixed(1) + (c.predictedTassExtrap ? " (clamped)" : ""))],
      ["Limit of detection", c.lod == null ? "never detected in the series" :
        kfmt(c.lod) + " " + c.unit + (c.lodMargin ? " · this sample is " + c.lodMargin.toFixed(1) + "× that depth" : "")],
    ];
    if (c.rolledUp && c.usedRollupRow) {
      // Keep the row the user actually clicked visible, and say how much of the
      // rolled-up total it accounts for.
      rows.splice(2, 0, [
        "This " + c.rowLevel.toLowerCase() + "'s reads",
        icomma(c.rowReads) + " " + c.unit +
          (c.rowShare != null ? " (" + pct(c.rowShare) + " of the " + c.compareLevel.toLowerCase() + ")" : "") +
          " · TASS " + c.rowTass.toFixed(1),
      ]);
    }
    if (c.rolledUp && c.org.n_members > 1 && c.org.members) {
      rows.push(["Genus members simulated", c.org.n_members + " · " + c.org.members.join(", ")]);
    }
    return (
      '<table style="border-collapse:collapse;font-size:.9em">' +
      rows.map(function (r) {
        return '<tr><td style="padding:1px 9px 1px 0;' + lab + ';white-space:nowrap">' + esc(r[0]) +
               '</td><td style="padding:1px 0;font-weight:600;white-space:nowrap">' + r[1] + "</td></tr>";
      }).join("") + "</table>"
    );
  }

  // Spelled out wherever the comparison is shown: which level the detection is,
  // which level the series is, and what that means for reading the plot.
  function levelNotice(c) {
    if (!c.rolledUp) {
      // Same-level comparison, but a genus series is still an assembly of its
      // simulated members — say so rather than implying one simulated taxon.
      if (c.genusSeries && c.org.n_members > 1) {
        return "This is a genus detection. No genus was simulated directly, so the series is assembled from the " +
               c.org.n_members + " simulated members of this genus (" + (c.org.members || []).join(", ") +
               "): their reads and expected share add, and the genus counts as detected wherever any member was.";
      }
      return null;
    }
    var lvl = c.rowLevel.toLowerCase();
    var txt =
      "You are looking at a " + lvl + " detection. " +
      "The dilution series exists at " + c.compareLevel.toLowerCase() + " level" +
      (c.compareName ? " (" + c.compareName + ")" : "") + ", so this plot is the " +
      c.compareLevel.toLowerCase() + "'s series" +
      (c.how === "genus" ? ", rolled up from every simulated member of the genus" : "") + ".";
    txt += c.usedRollupRow
      ? " The sample side uses this sample's " + c.compareLevel.toLowerCase() +
        "-level row, so both sides are at the same level."
      : " No " + c.compareLevel.toLowerCase() + "-level row exists for this sample, so the " + lvl +
        "'s own reads are plotted — a " + lvl + " is a subset of its " + c.compareLevel.toLowerCase() +
        ", so expect it to sit below the series.";
    return txt;
  }

  function levelBanner(c, dark) {
    var t = levelNotice(c);
    if (!t) return "";
    var head = c.rolledUp
      ? c.rowLevel + " detection · " + c.compareLevel + " series"
      : c.compareLevel + " series · assembled from " + c.org.n_members + " members";
    return (
      '<div style="margin:4px 0;padding:5px 8px;border-radius:5px;font-size:.85em;line-height:1.35;' +
      (dark ? "background:rgba(255,255,255,.09);color:#e8e2ff" : "background:#f3efff;color:#3b2a72") +
      ';border-left:3px solid ' + ACCENT + '">' +
      "<b>" + esc(head) + "</b><br>" + esc(t) + "</div>"
    );
  }

  function provenanceNote(c) {
    var bits = [];
    if (!c.ownSeries) bits.push("series simulated from " + c.group.parent);
    if (c.matchedByName) bits.push("matched by organism name, not taxid");
    if (c.unitMismatch) bits.push("sample is " + c.platform + " but the series is paired — depths are not directly comparable");
    if (c.nGroups > 1) bits.push(c.nGroups + " series available for this organism");
    return bits.length ? bits.join(" · ") : "";
  }

  // ── hover card ────────────────────────────────────────────────────────────
  function insilicoTipHTML(row) {
    var c = insilicoCompareFor(row);
    if (!c) return "";
    var v = verdict(c);
    return (
      '<div style="font-weight:700;margin-bottom:1px">' + esc(c.org.name) +
      " <span style='opacity:.7;font-weight:400'>vs dilution series" +
      (c.rolledUp ? " · " + esc(c.compareLevel) + " level" : "") + "</span></div>" +
      '<div style="font-size:.85em;color:' + v.color + ';font-weight:600;margin-bottom:4px">' + esc(v.label) + "</div>" +
      levelBanner(c, true) +
      '<div style="background:#fff;border-radius:5px;padding:2px 2px 0;margin-bottom:5px">' +
      panel(c, "observed_reads", { W: 360, plotH: 118 }) +
      '<div style="font-size:.72em;color:#555;padding:0 4px 4px">' +
      swatch("diamond", REAL) + "this sample &nbsp; " + swatch("ring", GOOD) + "series at this depth</div>" +
      "</div>" +
      statRows(c, true) +
      (provenanceNote(c) ? '<div style="font-size:.78em;opacity:.6;margin-top:4px">' + esc(provenanceNote(c)) + "</div>" : "") +
      '<div style="font-size:.78em;opacity:.6;margin-top:3px">click for the full comparison</div>'
    );
  }

  // ── modal ─────────────────────────────────────────────────────────────────
  function ensureModal() {
    var o = document.getElementById("insilico-cmp-overlay");
    if (o) return o;
    o = document.createElement("div");
    o.id = "insilico-cmp-overlay";
    o.setAttribute(
      "style",
      "display:none;position:fixed;inset:0;background:rgba(20,16,40,.45);z-index:9500;align-items:center;justify-content:center;padding:2vh 2vw"
    );
    o.innerHTML =
      '<div id="insilico-cmp-box" role="dialog" aria-modal="true" style="background:#fff;border-radius:10px;' +
      'max-width:980px;width:100%;max-height:96vh;overflow:auto;box-shadow:0 10px 40px rgba(0,0,0,.3);padding:1.1em 1.3em">' +
      '<div id="insilico-cmp-body"></div></div>';
    o.addEventListener("click", function (e) {
      if (e.target === o) closeModal();
    });
    document.addEventListener("keydown", function (e) {
      if (e.key === "Escape" && o.style.display !== "none") closeModal();
    });
    document.body.appendChild(o);
    return o;
  }

  function closeModal() {
    var o = document.getElementById("insilico-cmp-overlay");
    if (o) o.style.display = "none";
  }

  function seriesTable(c) {
    var s = (c.org.series || []).slice().sort(function (a, b) { return a.count - b.count; });
    var h =
      '<table style="border-collapse:collapse;width:100%;font-size:.84em;margin-top:.3em">' +
      "<thead><tr>" +
      ["Depth (" + c.unit + ")", "Expected", "Observed", "Recovery", "TASS", "Detected", "vs this sample"]
        .map(function (t) {
          return '<th style="text-align:left;padding:.35em .6em;border-bottom:2px solid ' + ACCENT +
                 ';white-space:nowrap;color:#333">' + esc(t) + "</th>";
        })
        .join("") +
      "</tr></thead><tbody>";
    s.forEach(function (p, i) {
      var rec = p.expected_reads ? p.observed_reads / p.expected_reads : null;
      var rel = c.depth > 0 ? p.count / c.depth : null;
      h +=
        '<tr style="' + (i % 2 ? "background:#faf9ff" : "") + '">' +
        '<td style="padding:.3em .6em;border-bottom:1px solid #eee">' + icomma(p.count) + "</td>" +
        '<td style="padding:.3em .6em;border-bottom:1px solid #eee">' + icomma(p.expected_reads) + "</td>" +
        '<td style="padding:.3em .6em;border-bottom:1px solid #eee">' + icomma(p.observed_reads) + "</td>" +
        '<td style="padding:.3em .6em;border-bottom:1px solid #eee">' + pct(rec) + "</td>" +
        '<td style="padding:.3em .6em;border-bottom:1px solid #eee">' + p.tass + "</td>" +
        '<td style="padding:.3em .6em;border-bottom:1px solid #eee;color:' + (p.detected ? GOOD : BAD) + '">' +
          (p.detected ? "yes" : "no") + (p.detection_rate > 0 && p.detection_rate < 1 ? " (" + pct(p.detection_rate) + " of reps)" : "") + "</td>" +
        '<td style="padding:.3em .6em;border-bottom:1px solid #eee;color:' + MUTED + '">' +
          (rel == null ? "—" : rel.toFixed(2) + "× this sample's depth") + "</td>" +
        "</tr>";
    });
    return h + "</tbody></table>";
  }

  function openInsilicoCompare(row, groupIdx) {
    var c = insilicoCompareFor(row, groupIdx);
    if (!c) return;
    var o = ensureModal();
    var v = verdict(c);
    var body = document.getElementById("insilico-cmp-body");
    var hits = matchOrganism(c.sample, row);
    var sel = "";
    if (hits.length > 1) {
      sel =
        '<select id="insilico-cmp-group" style="font-size:.85em;padding:.2em .4em;border:1px solid #ddd;border-radius:5px">' +
        hits.map(function (h, i) {
          return '<option value="' + i + '"' + (i === c.groupIdx ? " selected" : "") + ">" +
                 esc(h.matchLabel + " · " + h.group.parent + " · " + h.group.platform) + "</option>";
        }).join("") + "</select>";
    }
    body.innerHTML =
      '<div style="position:relative;z-index:8;display:flex;flex-wrap:wrap;align-items:baseline;' +
      'justify-content:space-between;gap:.6em;background:#fff">' +
      '<div><span style="font-size:1.1em;font-weight:700;color:' + ACCENT + '">' + esc(c.org.name) + "</span>" +
      '<span style="color:' + MUTED + ';font-size:.85em"> in ' + esc(c.sample) + " vs the in-silico dilution series</span></div>" +
      '<div style="display:flex;gap:.5em;align-items:center">' + sel +
      '<button type="button" id="insilico-cmp-close" style="border:1px solid #ddd;background:#fff;border-radius:6px;padding:.25em .6em;cursor:pointer">Close</button></div></div>' +
      '<div style="color:' + v.color + ';font-weight:600;margin:.3em 0 .1em">' + esc(v.label) + "</div>" +
      levelBanner(c, false) +
      (provenanceNote(c) ? '<div style="font-size:.8em;color:' + MUTED + '">' + esc(provenanceNote(c)) + "</div>" : "") +
      '<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:.9em;margin-top:.8em">' +
      '<div style="border:1px solid #eee;border-radius:8px;padding:.6em .7em">' +
      '<div style="font-size:.85em;font-weight:600;color:#333;margin-bottom:.2em">Reads vs depth</div>' +
      panel(c, "observed_reads", { W: 420, plotH: 150 }) + legendHTML(true, c) + "</div>" +
      '<div style="border:1px solid #eee;border-radius:8px;padding:.6em .7em">' +
      '<div style="font-size:.85em;font-weight:600;color:#333;margin-bottom:.2em">TASS vs depth</div>' +
      panel(c, "tass", { W: 420, plotH: 150 }) + legendHTML(false, c) + "</div>" +
      '<div style="border:1px solid #eee;border-radius:8px;padding:.6em .7em">' +
      '<div style="font-size:.85em;font-weight:600;color:#333;margin-bottom:.35em">This sample on the series</div>' +
      '<div style="position:relative">' + statRows(c, false) + "</div></div></div>" +
      '<div style="margin-top:1em;font-size:.9em;font-weight:600;color:#333">Dilution series datapoints</div>' +
      '<div style="position:relative;overflow-x:auto">' + seriesTable(c) + "</div>";
    o.style.display = "flex";
    var btn = document.getElementById("insilico-cmp-close");
    if (btn) btn.addEventListener("click", closeModal);
    var gs = document.getElementById("insilico-cmp-group");
    if (gs) gs.addEventListener("change", function () { openInsilicoCompare(row, +gs.value); });
  }

  // ── badge ─────────────────────────────────────────────────────────────────
  function insilicoBadgeHTML(row) {
    var c = insilicoCompareFor(row);
    if (!c) return "";
    var v = verdict(c);
    var arrow = c.vsSeries == null ? "" : c.vsSeries >= 2 ? "↑" : c.vsSeries <= 0.5 ? "↓" : "";
    // A rung marker when the series is above the detection's own level, so the
    // level shift is visible before the user even hovers.
    var rung = c.how === "genus" ? " gen" : c.how === "species" ? " sp" : "";
    var tip = c.rolledUp
      ? "In-silico dilution series available at " + c.compareLevel + " level (this row is " +
        c.rowLevel + ") — hover to compare, click for detail"
      : "In-silico dilution series available — hover to compare, click for detail";
    return (
      '<span class="insilico-badge" style="display:inline-block;margin-left:5px;font-size:9px;font-weight:700;' +
      "padding:0 4px;border-radius:3px;vertical-align:middle;line-height:1.5;cursor:pointer;" +
      "background:" + v.color + "1a;color:" + v.color + ";border:1px solid " + v.color + '55" ' +
      'title="' + esc(tip) + '">&#x2697;' + esc(rung) + " " + arrow + "</span>"
    );
  }

  // Hover/click wiring for a rendered table. `keyToRow` maps a row element's
  // data-key back to the data row (the summary table already keeps such a map).
  function wireInsilicoBadges(root, keyToRow) {
    if (!root || !suite()) return;
    root.querySelectorAll(".insilico-badge").forEach(function (b) {
      var tr = b.closest("tr[data-key]");
      if (!tr) return;
      var row = keyToRow(tr.dataset.key);
      if (!row) return;
      b.addEventListener("mouseover", function (ev) {
        if (typeof showTip === "function") showTip(insilicoTipHTML(row), ev);
      });
      b.addEventListener("mousemove", function (ev) {
        if (typeof moveTip === "function") moveTip(ev);
      });
      b.addEventListener("mouseout", function () {
        if (typeof hideTip === "function") hideTip();
      });
      b.addEventListener("click", function (ev) {
        ev.stopPropagation();
        if (typeof hideTip === "function") hideTip();
        openInsilicoCompare(row);
      });
    });
  }

  window.insilicoCompareFor = insilicoCompareFor;
  window.insilicoBadgeHTML = insilicoBadgeHTML;
  window.insilicoTipHTML = insilicoTipHTML;
  window.openInsilicoCompare = openInsilicoCompare;
  window.wireInsilicoBadges = wireInsilicoBadges;
})();
