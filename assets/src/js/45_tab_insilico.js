/* ═══════════════════════════════════════════════════════════════════════════
       -  §  IN-SILICO SUITE TAB
       -     Renders the spike-in / dilution-series subsampling suite:
       -       (1) a run-parameter provenance panel + suite export controls,
       -       (2) per (parent × platform) group, a per-dataset "expected vs
       -           reality" table (target/observed reads, TP/FP/FN, P/R/F1),
       -       (3) group-level metric charts (performance vs depth, detection
       -           composition vs depth, read recovery vs depth), and
       -       (4) a per-organism dilution-series limit-of-detection view with a
       -           compact expected-vs-observed chart per organism.
       -     Every plot is drawn inside a padded axis frame (own y-axis gutter,
       -     baseline rule, detection strip and x-label band are separate rows,
       -     band width is sized from the widest x label) so nothing overlaps the
       -     marks. Hovering any band/point/row raises the shared report tooltip
       -     with the full per-point statistics.
       -     Data source: BOOT.insilico_suite (built server-side in make_report.py
       -     from the subsample datasets, which each flow through as their own
       -     sample). The tab is hidden unless --sim_subsample was enabled.
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  "use strict";

  var ACCENT = "#4527a0";
  var ACCENT_L = "#b39ddb";
  var GOOD = "#2e7d32";
  var WARN = "#ef6c00";
  var BAD = "#c62828";
  var MUTED = "#777";
  var GRID = "#e9e5f3";
  var AXIS = "#bdb6d4";

  // ── DOM helper ────────────────────────────────────────────────────────────
  function el(tag, attrs, children) {
    var e = document.createElement(tag);
    if (attrs) {
      Object.keys(attrs).forEach(function (k) {
        if (k === "style") e.setAttribute("style", attrs[k]);
        else if (k === "class") e.className = attrs[k];
        else if (k === "html") e.innerHTML = attrs[k];
        else e.setAttribute(k, attrs[k]);
      });
    }
    (children || []).forEach(function (c) {
      if (c == null) return;
      e.appendChild(typeof c === "string" ? document.createTextNode(c) : c);
    });
    return e;
  }

  // ── formatting ────────────────────────────────────────────────────────────
  function fmt(n) {
    if (n == null || n === "") return "—";
    if (typeof n === "number") return n.toLocaleString();
    return String(n);
  }

  // Compact count label used on axes ("12.5k", "3M").
  function kfmt(n) {
    if (n == null || n === "" || isNaN(n)) return "—";
    n = +n;
    var a = Math.abs(n);
    if (a >= 1e6) return (n / 1e6).toFixed(a % 1e6 ? 1 : 0) + "M";
    if (a >= 1000) return (n / 1000).toFixed(a % 1000 ? 1 : 0) + "k";
    return String(Math.round(n * 100) / 100);
  }

  function pct(x, d) {
    if (x == null || isNaN(x)) return "—";
    return (x * 100).toFixed(d == null ? 1 : d) + "%";
  }

  function esc(s) {
    return String(s == null ? "" : s).replace(/[&<>"]/g, function (m) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[m];
    });
  }

  function ratio(obs, exp) {
    if (!exp) return null;
    return obs / exp;
  }

  function log2fc(obs, exp) {
    if (!exp || obs == null) return null;
    return Math.log((obs + 0.5) / (exp + 0.5)) / Math.LN2;
  }

  // Human labels for known param keys (unknown keys are shown verbatim).
  var PARAM_LABELS = {
    mode: "Subsample mode",
    series_counts: "Read-count series",
    replicates: "Replicates (per count)",
    detection_threshold: "Detection TASS cutoff",
    seed: "Random seed",
    sim_nreads: "Master reads (per sample)",
    iss_model: "ISS error model",
    iss_mode: "ISS mode",
    sim_ont_divisor: "ONT read divisor",
    abundance_source: "Abundance source",
    keep_fastq: "Keep subsampled FASTQs",
    sim_subsample_mode: "Subsample mode",
  };

  // ── tooltip plumbing ──────────────────────────────────────────────────────
  // Reuses the report's shared #tooltip (showTip/moveTip/hideTip in
  // 02_utilities.js). Marks carry their rendered HTML in a URI-encoded data-tt
  // attribute; one delegated listener per card raises it.
  function tipAttr(html) {
    return ' data-tt="' + encodeURIComponent(html) + '"';
  }

  function tipBody(title, sub, rows) {
    var h = '<div style="font-weight:700;margin-bottom:2px">' + esc(title) + "</div>";
    if (sub) h += '<div style="font-size:.85em;opacity:.75;margin-bottom:5px">' + esc(sub) + "</div>";
    h += '<table style="border-collapse:collapse;font-size:.9em">';
    rows.forEach(function (r) {
      if (!r) return;
      h +=
        '<tr><td style="padding:1px 9px 1px 0;opacity:.72;white-space:nowrap">' +
        esc(r[0]) +
        '</td><td style="padding:1px 0;font-weight:600;white-space:nowrap">' +
        (r[2] ? r[1] : esc(r[1])) +
        "</td></tr>";
    });
    return h + "</table>";
  }

  function attachTips(root) {
    if (!root || root._insilTips) return;
    root._insilTips = true;
    var find = function (e) {
      var t = e.target;
      if (!t || !t.closest) return null;
      return t.closest("[data-tt]");
    };
    root.addEventListener("mouseover", function (e) {
      var t = find(e);
      if (!t || typeof showTip !== "function") return;
      showTip(decodeURIComponent(t.getAttribute("data-tt")), e);
    });
    root.addEventListener("mousemove", function (e) {
      var t = find(e);
      if (t && typeof moveTip === "function") moveTip(e);
    });
    root.addEventListener("mouseout", function (e) {
      var t = find(e);
      if (t && typeof hideTip === "function") hideTip();
    });
  }

  // Hover affordance for the invisible band hit-areas (inline styles can't do
  // :hover, so one stylesheet is injected once).
  function ensureStyles() {
    if (document.getElementById("insilico-tab-style")) return;
    var s = document.createElement("style");
    s.id = "insilico-tab-style";
    s.textContent =
      ".insil-band{fill:" + ACCENT + ";opacity:0;cursor:crosshair}" +
      ".insil-band:hover{opacity:.07}" +
      ".insil-pt{cursor:crosshair}" +
      ".insil-pt:hover{stroke:#000;stroke-width:1.4}" +
      ".insil-row{cursor:crosshair}" +
      ".insil-row:hover{background:#f1ecff !important}" +
      ".insil-scroll{overflow-x:auto;overflow-y:hidden}" +
      ".insil-exp-btn{border:1px solid " + ACCENT_L + ";background:#fff;color:" + ACCENT +
      ";border-radius:6px;padding:.3em .7em;font-size:.8em;font-weight:600;cursor:pointer}" +
      ".insil-exp-btn:hover{background:#f3efff}";
    document.head.appendChild(s);
  }

  // ── axis frame ────────────────────────────────────────────────────────────
  // One shared geometry helper so every chart reserves real estate for its axes
  // instead of drawing labels on top of the marks. `extraB` is an extra band
  // between the baseline and the x labels (used for the detection strip).
  function axisFrame(n, labels, opts) {
    opts = opts || {};
    var yLabels = opts.yLabels || [];
    var maxY = yLabels.reduce(function (a, b) {
      return Math.max(a, String(b == null ? "" : b).length);
    }, 1);
    // The left gutter is sized from the widest y tick label (plus room for the
    // rotated axis title when there is one), so the title, the tick labels and
    // the plot area can never land on top of each other.
    var padL = opts.padL == null ? Math.ceil(8 + (opts.yTitle ? 13 : 0) + maxY * 6.2 + 6) : opts.padL;
    var padR = opts.padR == null ? 12 : opts.padR;
    var padT = opts.padT == null ? 14 : opts.padT;
    var plotH = opts.plotH == null ? 150 : opts.plotH;
    var labelH = opts.labelH == null ? 16 : opts.labelH;
    var xTitleH = opts.xTitle ? 13 : 0;
    var extraB = opts.extraB || 0;
    var maxLab = (labels || []).reduce(function (a, b) {
      return Math.max(a, String(b == null ? "" : b).length);
    }, 1);
    // Band is never narrower than the widest x label -> x labels never collide.
    var minBand = Math.max(opts.minBand || 44, maxLab * 6.6 + 12);
    var innerW = Math.max(opts.minInner || 200, n * minBand);
    return {
      padL: padL,
      padR: padR,
      padT: padT,
      plotH: plotH,
      labelH: labelH,
      xTitleH: xTitleH,
      extraB: extraB,
      innerW: innerW,
      W: padL + padR + innerW,
      // Height reserves every row explicitly: plot, detection strip (extraB),
      // x tick labels, x axis title.
      H: padT + plotH + extraB + labelH + xTitleH + 6,
      band: innerW / Math.max(1, n),
      yb: padT + plotH,
    };
  }

  function cx(f, i) {
    return f.padL + f.band * (i + 0.5);
  }

  // Y axis: gridlines + right-aligned labels inside the left gutter.
  function yAxis(f, ticks, labelFn, title) {
    var s = "";
    ticks.forEach(function (t) {
      var y = f.yb - t.p * f.plotH;
      s +=
        '<line x1="' + f.padL + '" y1="' + y.toFixed(1) + '" x2="' + (f.W - f.padR) +
        '" y2="' + y.toFixed(1) + '" stroke="' + GRID + '" stroke-width="1"/>' +
        '<text x="' + (f.padL - 6) + '" y="' + (y + 3).toFixed(1) +
        '" text-anchor="end" font-size="9" fill="' + MUTED + '">' + esc(labelFn(t.v)) + "</text>";
    });
    // baseline
    s +=
      '<line x1="' + f.padL + '" y1="' + f.yb + '" x2="' + (f.W - f.padR) +
      '" y2="' + f.yb + '" stroke="' + AXIS + '" stroke-width="1.2"/>' +
      '<line x1="' + f.padL + '" y1="' + f.padT + '" x2="' + f.padL +
      '" y2="' + f.yb + '" stroke="' + AXIS + '" stroke-width="1.2"/>';
    if (title) {
      s +=
        '<text transform="translate(11,' + (f.padT + f.plotH / 2) + ') rotate(-90)" text-anchor="middle" ' +
        'font-size="9" fill="' + MUTED + '">' + esc(title) + "</text>";
    }
    return s;
  }

  // X labels sit in their own band below the baseline (and below `extraB`).
  function xLabels(f, labels, title) {
    var y = f.yb + f.extraB + 13;
    var s = "";
    labels.forEach(function (lab, i) {
      s +=
        '<text x="' + cx(f, i).toFixed(1) + '" y="' + y +
        '" text-anchor="middle" font-size="9" fill="#666">' + esc(lab) + "</text>";
    });
    if (title) {
      s +=
        '<text x="' + (f.padL + f.innerW / 2) + '" y="' + (y + 12) +
        '" text-anchor="middle" font-size="8.5" fill="#999">' + esc(title) + "</text>";
    }
    return s;
  }

  function svgOpen(f) {
    return (
      '<svg viewBox="0 0 ' + f.W + " " + f.H + '" width="' + f.W + '" height="' + f.H +
      '" role="img" style="display:block">'
    );
  }

  function chartCard(title, hint, svgHtml) {
    var wrap = el("div", {
      class: "chart-wrap",
      style: "position:relative;border:1px solid #eee;border-radius:8px;padding:.55em .65em;background:#fff",
    });
    wrap.appendChild(
      el("div", { class: "chart-title", style: "font-size:.85em;font-weight:600;color:#333;margin-bottom:.15em" }, [title])
    );
    if (hint) wrap.appendChild(el("div", { style: "font-size:.72em;color:" + MUTED + ";margin-bottom:.3em" }, [hint]));
    wrap.appendChild(el("div", { class: "insil-scroll", html: svgHtml }));
    return wrap;
  }

  function legend(items) {
    return el(
      "div",
      { style: "font-size:.74em;color:" + MUTED + ";margin-top:.3em;display:flex;flex-wrap:wrap;gap:.8em" },
      items.map(function (it) {
        return el("span", {
          html:
            '<span style="display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:4px;vertical-align:-1px;' +
            (it.outline ? "border:1.4px solid " + it.color : "background:" + it.color) +
            '"></span>' + esc(it.label),
        });
      })
    );
  }

  // ── params panel + export toolbar ─────────────────────────────────────────
  function renderParams(suite) {
    var params = suite.params || {};
    var host = document.getElementById("insilico-params-panel");
    if (!host) return;
    host.innerHTML = "";
    host.setAttribute(
      "style",
      "background:#faf9ff;border:1px solid #e6e1f5;border-radius:8px;padding:1em 1.2em;margin-bottom:1.1em"
    );
    var head = el("div", {
      style: "display:flex;flex-wrap:wrap;align-items:center;justify-content:space-between;gap:.6em",
    });
    head.appendChild(
      el("div", {
        class: "chart-title",
        html: '<i class="fas fa-flask"></i> In-Silico Subsampling — parameters used',
      })
    );
    head.appendChild(exportBar(suite));
    host.appendChild(head);

    var grid = el("div", {
      style:
        "display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:.5em .9em;margin-top:.6em",
    });
    var keys = Object.keys(params);
    var ordered = Object.keys(PARAM_LABELS).filter(function (k) {
      return keys.indexOf(k) !== -1;
    });
    keys.forEach(function (k) {
      if (ordered.indexOf(k) === -1) ordered.push(k);
    });
    if (!ordered.length) {
      grid.appendChild(el("div", { style: "color:" + MUTED }, ["No parameters recorded."]));
    }
    ordered.forEach(function (k) {
      var v = params[k];
      if (Array.isArray(v)) v = v.join(", ");
      if (typeof v === "boolean") v = v ? "yes" : "no";
      var cell = el("div", {
        style:
          "display:flex;flex-direction:column;padding:.35em .5em;background:#fff;border:1px solid #ece8f7;border-radius:6px",
      });
      cell.appendChild(
        el("span", { style: "font-size:.72em;color:" + MUTED + ";text-transform:uppercase;letter-spacing:.03em" }, [
          PARAM_LABELS[k] || k,
        ])
      );
      cell.appendChild(
        el("span", { style: "font-size:.95em;font-weight:600;color:#222;margin-top:2px" }, [
          v == null || v === "" ? "—" : String(v),
        ])
      );
      grid.appendChild(cell);
    });
    host.appendChild(grid);
  }

  // ── export ────────────────────────────────────────────────────────────────
  function dl(text, filename, type) {
    if (typeof _downloadText === "function") {
      _downloadText(text, filename, type);
      return;
    }
    var blob = new Blob([text], { type: type || "text/plain;charset=utf-8" });
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    setTimeout(function () {
      URL.revokeObjectURL(a.href);
      a.remove();
    }, 250);
  }

  function delimEscape(v, d) {
    var s = v == null ? "" : String(v);
    return /["\r\n]/.test(s) || s.indexOf(d) !== -1 ? '"' + s.replace(/"/g, '""') + '"' : s;
  }

  function toDelim(rows, d) {
    return rows
      .map(function (r) {
        return r
          .map(function (v) {
            return delimEscape(v, d);
          })
          .join(d);
      })
      .join("\r\n");
  }

  function paramRows(suite) {
    var p = suite.params || {};
    var rows = [["Parameter", "Value"]];
    Object.keys(p).forEach(function (k) {
      var v = p[k];
      if (Array.isArray(v)) v = v.join(", ");
      rows.push([PARAM_LABELS[k] || k, v == null ? "" : String(v)]);
    });
    return rows;
  }

  function datasetRows(suite) {
    var rows = [
      ["Parent", "Platform", "Mode", "Read unit", "Dataset", "Replicate", "Target", "Actual",
       "Master total", "Seed", "Observed aligned", "Recovery", "Detected", "TP", "FP", "FN",
       "Precision", "Recall", "F1"],
    ];
    (suite.groups || []).forEach(function (g) {
      (g.datasets || []).forEach(function (d) {
        var rec = ratio(d.observed_total_reads, d.actual_count);
        rows.push([
          g.parent, g.platform, g.mode, g.read_unit || "reads", d.id, d.replicate,
          d.target_count, d.actual_count, d.total_master_reads == null ? "" : d.total_master_reads,
          d.seed == null ? "" : d.seed, d.observed_total_reads,
          rec == null ? "" : rec.toFixed(4), d.n_detected, d.tp, d.fp, d.fn,
          d.precision, d.recall, d.f1,
        ]);
      });
    });
    return rows;
  }

  function organismRows(suite) {
    var rows = [
      ["Parent", "Platform", "Level", "Taxid", "Organism", "Members", "Category", "Expected fraction",
       "LoD count", "Target count", "Expected reads", "Observed reads", "Recovery", "log2 FC", "TASS",
       "Detection rate", "Detected", "Replicates"],
    ];
    (suite.groups || []).forEach(function (g) {
      // Export what is on screen: the rollup the user is looking at.
      rollupOrganisms(g).forEach(function (o) {
        (o.series || []).forEach(function (s) {
          var rec = ratio(s.observed_reads, s.expected_reads);
          var lfc = log2fc(s.observed_reads, s.expected_reads);
          rows.push([
            g.parent, g.platform, o.rolled ? "Genus" : (g.level || "Strain"), o.taxid, o.name,
            o.members ? o.members.join("; ") : o.name, o.category, o.expected_fraction,
            o.lod_count == null ? "" : o.lod_count, s.count, s.expected_reads, s.observed_reads,
            rec == null ? "" : rec.toFixed(4), lfc == null ? "" : lfc.toFixed(3),
            s.tass, s.detection_rate, s.detected ? "yes" : "no", s.n_reps,
          ]);
        });
      });
    });
    return rows;
  }

  function exportSuite(format, suite) {
    var stamp = new Date().toISOString().slice(0, 10);
    var base = "taxtriage-insilico-suite-" + stamp;
    if (format === "json") {
      dl(JSON.stringify(suite, null, 2), base + ".json", "application/json;charset=utf-8");
      return;
    }
    if (format === "xlsx") {
      if (typeof XLSX === "undefined") {
        alert("XLSX export needs the SheetJS library, which this report did not load. Use CSV, TSV or JSON.");
        return;
      }
      var wb = XLSX.utils.book_new();
      XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(paramRows(suite)), "Parameters");
      XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(datasetRows(suite)), "Datasets");
      XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(organismRows(suite)), "Organism series");
      XLSX.writeFile(wb, base + ".xlsx");
      return;
    }
    var d = format === "tsv" ? "\t" : ",";
    var ext = format === "tsv" ? "tsv" : "csv";
    dl(toDelim(datasetRows(suite), d), base + "-datasets." + ext, "text/plain;charset=utf-8");
    // Second file follows a beat later so browsers don't collapse the two
    // programmatic downloads into one.
    setTimeout(function () {
      dl(toDelim(organismRows(suite), d), base + "-organism-series." + ext, "text/plain;charset=utf-8");
    }, 400);
  }

  function exportBar(suite) {
    var bar = el("div", { style: "display:flex;align-items:center;gap:.4em" });
    bar.appendChild(el("span", { style: "font-size:.75em;color:" + MUTED }, ["Export suite:"]));
    [
      ["CSV", "csv", "Datasets + organism series as two CSV files"],
      ["TSV", "tsv", "Datasets + organism series as two TSV files"],
      ["XLSX", "xlsx", "One workbook: Parameters, Datasets, Organism series"],
      ["JSON", "json", "Raw suite payload exactly as the report holds it"],
    ].forEach(function (o) {
      var b = el("button", { type: "button", class: "insil-exp-btn", title: o[2] }, [o[0]]);
      b.addEventListener("click", function () {
        exportSuite(o[1], suite);
      });
      bar.appendChild(b);
    });
    return bar;
  }

  // ── misc UI bits ──────────────────────────────────────────────────────────
  function pill(txt, color) {
    return el(
      "span",
      {
        style:
          "display:inline-block;padding:1px 7px;border-radius:10px;font-size:.78em;font-weight:600;color:#fff;background:" +
          color,
      },
      [txt]
    );
  }

  function f1Color(f1) {
    if (f1 >= 0.9) return GOOD;
    if (f1 >= 0.6) return WARN;
    return BAD;
  }

  // ── per-dataset table ─────────────────────────────────────────────────────
  function renderDatasetTable(group) {
    var unit = group.read_unit || "reads";
    var wrap = el("div", { style: "position:relative;overflow-x:auto;margin:.4em 0 1em" });
    var t = el("table", {
      style: "border-collapse:collapse;width:100%;font-size:.85em;min-width:720px",
    });
    var heads = [
      "Dataset (rep)",
      "Target " + unit,
      "Actual " + unit,
      "Observed aligned (" + unit + ")",
      "Recovery",
      "Detected",
      "TP",
      "FP",
      "FN",
      "Precision",
      "Recall",
      "F1",
    ];
    var thead = el("thead");
    var hr = el("tr");
    heads.forEach(function (h) {
      hr.appendChild(
        el(
          "th",
          {
            style:
              "text-align:left;padding:.4em .6em;border-bottom:2px solid " +
              ACCENT +
              ";white-space:nowrap;font-size:.92em;color:#333",
          },
          [h]
        )
      );
    });
    thead.appendChild(hr);
    t.appendChild(thead);
    var tb = el("tbody");
    group.datasets.forEach(function (d, i) {
      var rec = ratio(d.observed_total_reads, d.actual_count);
      var tr = el("tr", {
        class: "insil-row",
        style: i % 2 ? "background:#faf9ff" : "",
      });
      tr.setAttribute(
        "data-tt",
        encodeURIComponent(
          tipBody(d.id, group.parent + " · " + group.platform + " · " + group.mode + " · rep " + d.replicate, [
            ["Target depth", fmt(d.target_count) + " " + unit],
            ["Actual depth", fmt(d.actual_count) + " " + unit],
            d.total_master_reads != null ? ["Master pool", fmt(d.total_master_reads) + " " + unit] : null,
            d.seed != null ? ["Seed", String(d.seed)] : null,
            ["Observed aligned", fmt(d.observed_total_reads) + " " + unit],
            ["Alignment recovery", rec == null ? "—" : pct(rec)],
            ["Organisms detected", fmt(d.n_detected)],
            ["True positives", String(d.tp)],
            ["False positives", String(d.fp)],
            ["False negatives", String(d.fn)],
            ["Precision", d.precision.toFixed(3)],
            ["Recall", d.recall.toFixed(3)],
            [
              "F1",
              '<span style="color:' + f1Color(d.f1) + '">' + d.f1.toFixed(3) + "</span>",
              true,
            ],
          ])
        )
      );
      var shortId = "c" + d.target_count + " · r" + d.replicate;
      var cells = [
        el("span", { style: "font-family:monospace;font-size:.92em" }, [shortId]),
        fmt(d.target_count),
        fmt(d.actual_count),
        fmt(d.observed_total_reads),
        rec == null ? "—" : pct(rec, 0),
        fmt(d.n_detected),
        String(d.tp),
        String(d.fp),
        String(d.fn),
        d.precision.toFixed(2),
        d.recall.toFixed(2),
        el("span", { style: "font-weight:700;color:" + f1Color(d.f1) }, [d.f1.toFixed(2)]),
      ];
      cells.forEach(function (c) {
        tr.appendChild(
          el("td", { style: "padding:.35em .6em;border-bottom:1px solid #eee;white-space:nowrap" }, [
            typeof c === "string" ? document.createTextNode(c) : c,
          ])
        );
      });
      tb.appendChild(tr);
    });
    t.appendChild(tb);
    wrap.appendChild(t);
    return wrap;
  }

  // ── group-level aggregation across replicates ─────────────────────────────
  function aggByCount(group) {
    var m = {};
    (group.datasets || []).forEach(function (d) {
      var c = d.target_count;
      if (!m[c]) {
        m[c] = { count: c, n: 0, tp: 0, fp: 0, fn: 0, precision: 0, recall: 0, f1: 0,
                 observed: 0, actual: 0, detected: 0, f1s: [] };
      }
      var a = m[c];
      a.n++;
      a.tp += d.tp; a.fp += d.fp; a.fn += d.fn;
      a.precision += d.precision; a.recall += d.recall; a.f1 += d.f1;
      a.f1s.push(d.f1);
      a.observed += d.observed_total_reads;
      a.actual += d.actual_count;
      a.detected += d.n_detected;
    });
    return Object.keys(m)
      .map(function (k) { return m[k]; })
      .sort(function (a, b) { return a.count - b.count; })
      .map(function (a) {
        var sd = 0;
        if (a.f1s.length > 1) {
          var mu = a.f1 / a.n;
          sd = Math.sqrt(a.f1s.reduce(function (s, v) { return s + (v - mu) * (v - mu); }, 0) / (a.f1s.length - 1));
        }
        return {
          count: a.count, n: a.n,
          tp: a.tp / a.n, fp: a.fp / a.n, fn: a.fn / a.n,
          precision: a.precision / a.n, recall: a.recall / a.n, f1: a.f1 / a.n, f1_sd: sd,
          observed: a.observed / a.n, actual: a.actual / a.n, detected: a.detected / a.n,
        };
      });
  }

  // Chart 1 — precision / recall / F1 vs sequencing depth.
  function metricChart(agg, group) {
    var labels = agg.map(function (a) { return kfmt(a.count); });
    var ticks = [0, 0.25, 0.5, 0.75, 1].map(function (v) { return { v: v, p: v }; });
    var yLab = function (v) { return (v * 100).toFixed(0) + "%"; };
    var f = axisFrame(agg.length, labels, {
      plotH: 132,
      yLabels: ticks.map(function (t) { return yLab(t.v); }),
      yTitle: true,
      xTitle: true,
    });
    var s = svgOpen(f);
    s += yAxis(f, ticks, yLab, "score");
    var seriesDefs = [
      { key: "precision", color: "#1565c0", label: "Precision" },
      { key: "recall", color: "#ef6c00", label: "Recall" },
      { key: "f1", color: ACCENT, label: "F1" },
    ];
    seriesDefs.forEach(function (sd) {
      var pts = agg.map(function (a, i) {
        return [cx(f, i), f.yb - Math.max(0, Math.min(1, a[sd.key])) * f.plotH];
      });
      s +=
        '<polyline fill="none" stroke="' + sd.color + '" stroke-width="1.8" stroke-linejoin="round" points="' +
        pts.map(function (p) { return p[0].toFixed(1) + "," + p[1].toFixed(1); }).join(" ") + '"/>';
      pts.forEach(function (p) {
        s += '<circle class="insil-pt" cx="' + p[0].toFixed(1) + '" cy="' + p[1].toFixed(1) +
             '" r="3.2" fill="' + sd.color + '"/>';
      });
    });
    // band hit-areas last so they sit above the marks
    agg.forEach(function (a, i) {
      var tip = tipBody(
        kfmt(a.count) + " " + (group.read_unit || "reads"),
        group.parent + " · " + group.platform + " · " + a.n + " replicate" + (a.n === 1 ? "" : "s"),
        [
          ["Precision", a.precision.toFixed(3)],
          ["Recall", a.recall.toFixed(3)],
          ["F1 (mean)", a.f1.toFixed(3)],
          a.n > 1 ? ["F1 SD", a.f1_sd.toFixed(3)] : null,
          ["TP / FP / FN", a.tp.toFixed(1) + " / " + a.fp.toFixed(1) + " / " + a.fn.toFixed(1)],
          ["Organisms detected", a.detected.toFixed(1)],
        ]
      );
      s += '<rect class="insil-band" x="' + (f.padL + f.band * i).toFixed(1) + '" y="' + f.padT +
           '" width="' + f.band.toFixed(1) + '" height="' + (f.plotH + f.extraB) + '"' + tipAttr(tip) + "/>";
    });
    s += xLabels(f, labels, "target depth (" + (group.read_unit || "reads") + ")");
    s += "</svg>";
    var card = chartCard("Performance vs depth", "Mean across replicates at each target depth.", s);
    card.appendChild(
      legend([
        { label: "Precision", color: "#1565c0" },
        { label: "Recall", color: "#ef6c00" },
        { label: "F1", color: ACCENT },
      ])
    );
    return card;
  }

  // Chart 2 — detection composition (TP / FP / FN) vs depth, stacked.
  function compositionChart(agg, group) {
    var labels = agg.map(function (a) { return kfmt(a.count); });
    var maxV = 1;
    agg.forEach(function (a) { maxV = Math.max(maxV, a.tp + a.fp + a.fn); });
    var ticks = [0, 0.5, 1].map(function (p) { return { v: maxV * p, p: p }; });
    var f = axisFrame(agg.length, labels, {
      plotH: 132,
      yLabels: ticks.map(function (t) { return kfmt(t.v); }),
      yTitle: true,
      xTitle: true,
    });
    var sc = function (v) { return (v / maxV) * f.plotH; };
    var s = svgOpen(f);
    s += yAxis(f, ticks, function (v) { return kfmt(v); }, "organisms");
    var bw = Math.min(26, f.band * 0.5);
    agg.forEach(function (a, i) {
      var x = cx(f, i) - bw / 2;
      var y = f.yb;
      [["tp", GOOD], ["fp", BAD], ["fn", "#9e9e9e"]].forEach(function (p) {
        var h = sc(a[p[0]]);
        if (h <= 0) return;
        y -= h;
        s += '<rect x="' + x.toFixed(1) + '" y="' + y.toFixed(1) + '" width="' + bw.toFixed(1) +
             '" height="' + h.toFixed(1) + '" fill="' + p[1] + '" opacity=".85"/>';
      });
      var tip = tipBody(
        kfmt(a.count) + " " + (group.read_unit || "reads"),
        group.parent + " · " + group.platform,
        [
          ["True positives", a.tp.toFixed(1)],
          ["False positives", a.fp.toFixed(1)],
          ["False negatives", a.fn.toFixed(1)],
          ["Total detected", a.detected.toFixed(1)],
          ["Precision", a.precision.toFixed(3)],
          ["Recall", a.recall.toFixed(3)],
          ["Replicates", String(a.n)],
        ]
      );
      s += '<rect class="insil-band" x="' + (f.padL + f.band * i).toFixed(1) + '" y="' + f.padT +
           '" width="' + f.band.toFixed(1) + '" height="' + f.plotH + '"' + tipAttr(tip) + "/>";
    });
    s += xLabels(f, labels, "target depth (" + (group.read_unit || "reads") + ")");
    s += "</svg>";
    var card = chartCard("Detection composition vs depth", "Mean organism counts per depth.", s);
    card.appendChild(
      legend([
        { label: "True positive", color: GOOD },
        { label: "False positive", color: BAD },
        { label: "False negative (missed)", color: "#9e9e9e" },
      ])
    );
    return card;
  }

  // Chart 3 — how much of each subsample actually aligned.
  function recoveryChart(agg, group) {
    var unit = group.read_unit || "reads";
    var labels = agg.map(function (a) { return kfmt(a.count); });
    var maxV = 1;
    agg.forEach(function (a) { maxV = Math.max(maxV, a.actual, a.observed); });
    var ticks = [0, 0.25, 1].map(function (p) { return { v: maxV * p, p: Math.sqrt(p) }; });
    var f = axisFrame(agg.length, labels, {
      plotH: 132,
      yLabels: ticks.map(function (t) { return kfmt(t.v); }),
      yTitle: true,
      xTitle: true,
    });
    var sc = function (v) { return (Math.sqrt(Math.max(0, v)) / Math.sqrt(maxV)) * f.plotH; };
    var s = svgOpen(f);
    s += yAxis(f, ticks, function (v) { return kfmt(v); }, unit + " (√ scale)");
    var bw = Math.min(16, f.band * 0.3);
    agg.forEach(function (a, i) {
      var c = cx(f, i);
      var ah = sc(a.actual), oh = sc(a.observed);
      s += '<rect x="' + (c - bw - 1).toFixed(1) + '" y="' + (f.yb - ah).toFixed(1) + '" width="' + bw.toFixed(1) +
           '" height="' + ah.toFixed(1) + '" fill="none" stroke="' + ACCENT + '" stroke-width="1.2" rx="1"/>';
      s += '<rect x="' + (c + 1).toFixed(1) + '" y="' + (f.yb - oh).toFixed(1) + '" width="' + bw.toFixed(1) +
           '" height="' + oh.toFixed(1) + '" fill="' + ACCENT + '" opacity=".8" rx="1"/>';
      var rec = ratio(a.observed, a.actual);
      var tip = tipBody(
        kfmt(a.count) + " " + unit,
        group.parent + " · " + group.platform,
        [
          ["Actual in dataset", fmt(Math.round(a.actual)) + " " + unit],
          ["Observed aligned", fmt(Math.round(a.observed)) + " " + unit],
          ["Recovery", rec == null ? "—" : pct(rec)],
          ["Unaligned", fmt(Math.max(0, Math.round(a.actual - a.observed))) + " " + unit],
          ["Replicates", String(a.n)],
        ]
      );
      s += '<rect class="insil-band" x="' + (f.padL + f.band * i).toFixed(1) + '" y="' + f.padT +
           '" width="' + f.band.toFixed(1) + '" height="' + f.plotH + '"' + tipAttr(tip) + "/>";
    });
    s += xLabels(f, labels, "target depth (" + unit + ")");
    s += "</svg>";
    var card = chartCard("Read recovery vs depth", "How much of each subsample aligned to a reference.", s);
    card.appendChild(
      legend([
        { label: "In dataset", color: ACCENT, outline: true },
        { label: "Aligned", color: ACCENT },
      ])
    );
    return card;
  }

  function renderGroupCharts(group) {
    var agg = aggByCount(group);
    if (!agg.length) return null;
    var wrap = el("div", { style: "margin:.2em 0 1.1em" });
    wrap.appendChild(
      el("div", { style: "font-weight:600;color:#333;margin:.2em 0 .5em;font-size:.95em" }, [
        "Run metrics across the dilution series",
      ])
    );
    var grid = el("div", {
      style: "display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:.8em",
    });
    grid.appendChild(metricChart(agg, group));
    grid.appendChild(compositionChart(agg, group));
    grid.appendChild(recoveryChart(agg, group));
    wrap.appendChild(grid);
    return wrap;
  }

  // ── taxonomic rollup ──────────────────────────────────────────────────────
  // The suite is built at ONE level (group.level — Species when available, else
  // Strain). Rolling up to Genus merges every organism sharing a genus into a
  // single series: reads add, TASS takes the strongest member (a genus is as
  // detectable as its best-supported member, not the sum of their scores), and
  // the genus counts as detected at a depth if any member was.
  var INSIL_ROLLUP = "native";   // "native" | "Genus"

  function rollupOrganisms(group) {
    var orgs = group.organisms || [];
    if (INSIL_ROLLUP !== "Genus") return orgs;
    var buckets = {};
    orgs.forEach(function (o) {
      var g = String(o.genus || "").trim();
      // An organism with no genus label stays on its own rather than being
      // silently merged into an "Unknown" bucket with unrelated taxa.
      var key = g ? "g:" + g.toLowerCase() : "t:" + o.taxid;
      if (!buckets[key]) buckets[key] = { genus: g, members: [] };
      buckets[key].members.push(o);
    });
    var out = Object.keys(buckets).map(function (k) {
      var b = buckets[k];
      var ms = b.members;
      if (ms.length === 1 && b.genus) {
        // Single member: keep its real numbers, just relabel to the genus.
        var only = ms[0];
        return {
          taxid: only.taxid, name: b.genus, category: only.category,
          expected_fraction: only.expected_fraction, lod_count: only.lod_count,
          series: only.series, rolled: true, n_members: 1,
          members: [only.name], member_taxids: [only.taxid],
        };
      }
      if (!b.genus) return ms[0];
      // Union the count axis (all members share it, but be defensive).
      var counts = {};
      ms.forEach(function (o) {
        (o.series || []).forEach(function (p) { counts[p.count] = 1; });
      });
      var cs = Object.keys(counts).map(Number).sort(function (a, b2) { return a - b2; });
      var lod = null;
      var series = cs.map(function (c) {
        var exp = 0, obsr = 0, tass = 0, dr = 0, det = false, nreps = 1;
        ms.forEach(function (o) {
          var p = (o.series || []).find(function (q) { return q.count === c; });
          if (!p) return;
          exp += +p.expected_reads || 0;
          obsr += +p.observed_reads || 0;
          tass = Math.max(tass, +p.tass || 0);
          dr = Math.max(dr, +p.detection_rate || 0);
          det = det || !!p.detected;
          nreps = Math.max(nreps, +p.n_reps || 1);
        });
        if (det && lod === null) lod = c;
        return {
          count: c,
          expected_reads: Math.round(exp * 10) / 10,
          observed_reads: Math.round(obsr * 10) / 10,
          tass: Math.round(tass * 100) / 100,
          detection_rate: Math.round(dr * 1000) / 1000,
          detected: det,
          n_reps: nreps,
        };
      });
      var cats = {};
      ms.forEach(function (o) { cats[o.category || "Unknown"] = 1; });
      var catKeys = Object.keys(cats);
      return {
        taxid: ms.map(function (o) { return o.taxid; }).join("+"),
        name: b.genus,
        category: catKeys.length === 1 ? catKeys[0] : catKeys.join(" / "),
        expected_fraction: ms.reduce(function (a, o) { return a + (+o.expected_fraction || 0); }, 0),
        lod_count: lod,
        series: series,
        rolled: true,
        n_members: ms.length,
        members: ms.map(function (o) { return o.name; }),
        member_taxids: ms.map(function (o) { return o.taxid; }),
      };
    });
    out.sort(function (a, b) { return (b.expected_fraction || 0) - (a.expected_fraction || 0); });
    return out;
  }

  // ── per-organism dilution chart ───────────────────────────────────────────
  // Expected (outline) vs observed (filled) bars over the count series, with a
  // separate detection strip below the baseline and x labels below that. Band
  // width is derived from the widest label, so marks and axis never collide.
  function organismChart(o, group) {
    var series = o.series || [];
    var unit = group.read_unit || "reads";
    var labels = series.map(function (s) { return kfmt(s.count); });
    var maxV = 1;
    series.forEach(function (s) {
      maxV = Math.max(maxV, s.expected_reads, s.observed_reads);
    });
    var ticks = [0, 0.25, 1].map(function (p) { return { v: maxV * p, p: Math.sqrt(p) }; });
    var f = axisFrame(series.length, labels, {
      plotH: 92,
      padR: 10,
      padT: 10,
      extraB: 18,   // detection strip
      minBand: 46,
      minInner: 170,
      yLabels: ticks.map(function (t) { return kfmt(t.v); }),
      yTitle: true,
      xTitle: true,
    });
    // sqrt scale keeps small dilutions visible next to large ones
    var sc = function (v) { return (Math.sqrt(Math.max(0, v)) / Math.sqrt(maxV)) * f.plotH; };
    var s = svgOpen(f);
    s += yAxis(f, ticks, function (v) { return kfmt(v); }, unit + " (√)");
    var dotY = f.yb + f.extraB / 2 + 1;
    var bw = Math.min(16, f.band * 0.3);

    series.forEach(function (pt, i) {
      var c = cx(f, i);
      var eh = sc(pt.expected_reads);
      var oh = sc(pt.observed_reads);
      // expected (outline)
      s += '<rect x="' + (c - bw - 1).toFixed(1) + '" y="' + (f.yb - eh).toFixed(1) + '" width="' + bw.toFixed(1) +
           '" height="' + eh.toFixed(1) + '" fill="none" stroke="' + ACCENT + '" stroke-width="1.2" rx="1"/>';
      // observed (filled, coloured by detection)
      s += '<rect x="' + (c + 1).toFixed(1) + '" y="' + (f.yb - oh).toFixed(1) + '" width="' + bw.toFixed(1) +
           '" height="' + oh.toFixed(1) + '" fill="' + (pt.detected ? GOOD : BAD) + '" opacity=".85" rx="1"/>';
      // detection strip (partial detection across replicates -> amber)
      var dr = pt.detection_rate;
      var dotFill = dr >= 1 ? GOOD : dr > 0 ? WARN : "#e4e4e4";
      var dotStroke = dr > 0 ? (dr >= 1 ? GOOD : WARN) : "#c4c4c4";
      s += '<circle cx="' + c.toFixed(1) + '" cy="' + dotY.toFixed(1) + '" r="4" fill="' + dotFill +
           '" stroke="' + dotStroke + '" stroke-width="1"/>';
      // LoD marker
      if (o.lod_count != null && pt.count === o.lod_count) {
        s += '<line x1="' + c.toFixed(1) + '" y1="' + f.padT + '" x2="' + c.toFixed(1) + '" y2="' + f.yb +
             '" stroke="' + ACCENT + '" stroke-width="1" stroke-dasharray="3 3" opacity=".55"/>';
      }
    });

    // hit-areas on top
    series.forEach(function (pt, i) {
      var rec = ratio(pt.observed_reads, pt.expected_reads);
      var lfc = log2fc(pt.observed_reads, pt.expected_reads);
      var tip = tipBody(o.name, "taxid " + o.taxid + " · " + kfmt(pt.count) + " " + unit + " target", [
        ["Expected", fmt(pt.expected_reads) + " " + unit],
        ["Observed", fmt(pt.observed_reads) + " " + unit],
        ["Recovery", rec == null ? "—" : pct(rec)],
        ["log2 fold change", lfc == null ? "—" : (lfc > 0 ? "+" : "") + lfc.toFixed(2)],
        ["TASS (mean)", String(pt.tass)],
        ["Detection rate", pct(pt.detection_rate, 0) + " of " + pt.n_reps + " rep" + (pt.n_reps === 1 ? "" : "s")],
        ["Called detected", pt.detected ? "yes" : "no"],
        ["Expected share of pool", pct(o.expected_fraction)],
        ["Limit of detection", o.lod_count == null ? "not reached" : kfmt(o.lod_count) + " " + unit],
      ]);
      s += '<rect class="insil-band" x="' + (f.padL + f.band * i).toFixed(1) + '" y="' + f.padT +
           '" width="' + f.band.toFixed(1) + '" height="' + (f.plotH + f.extraB) + '"' + tipAttr(tip) + "/>";
    });

    s += xLabels(f, labels, "target depth (" + unit + ")");
    s += "</svg>";
    return s;
  }

  function renderOrganisms(group) {
    var wrap = el("div", { style: "margin-top:.3em" });
    var nativeLevel = group.level || "Strain";
    var head = el("div", {
      style: "display:flex;flex-wrap:wrap;align-items:baseline;justify-content:space-between;gap:.6em;margin:.2em 0 .5em",
    });
    head.appendChild(
      el("div", { style: "font-weight:600;color:#333;font-size:.95em" }, [
        "Per-organism dilution series — expected vs observed (limit of detection)",
      ])
    );
    // Rollup control: the series is built at one level; Genus merges members.
    var ctl = el("div", { style: "display:flex;align-items:center;gap:.4em;font-size:.8em;color:" + MUTED });
    ctl.appendChild(el("span", {}, ["Level:"]));
    [[nativeLevel, "native"], ["Genus", "Genus"]].forEach(function (o) {
      var on = INSIL_ROLLUP === o[1];
      var b = el("button", {
        type: "button",
        style:
          "border:1px solid " + (on ? ACCENT : "#ddd") + ";background:" + (on ? ACCENT : "#fff") +
          ";color:" + (on ? "#fff" : "#555") + ";border-radius:6px;padding:.2em .6em;font-size:.95em;" +
          "font-weight:600;cursor:pointer",
        title: o[1] === "Genus"
          ? "Merge members of the same genus into one series (reads add, TASS takes the strongest member)"
          : "The level the series was simulated and scored at",
      }, [o[0]]);
      b.addEventListener("click", function () {
        INSIL_ROLLUP = o[1];
        window.drawInsilico();
      });
      ctl.appendChild(b);
    });
    head.appendChild(ctl);
    wrap.appendChild(head);
    wrap.appendChild(
      el("div", {
        style: "font-size:.76em;color:" + MUTED + ";margin-bottom:.6em",
        html:
          '<span style="border:1.2px solid ' + ACCENT + ';padding:0 6px;border-radius:2px;margin-right:4px">&nbsp;</span>expected &nbsp;&nbsp;' +
          '<span style="background:' + GOOD + ';color:#fff;padding:0 6px;border-radius:2px;margin-right:4px">obs</span>detected &nbsp;&nbsp;' +
          '<span style="background:' + BAD + ';color:#fff;padding:0 6px;border-radius:2px;margin-right:4px">obs</span>below threshold &nbsp;&nbsp;' +
          '<span style="display:inline-block;width:9px;height:9px;border-radius:50%;background:' + WARN + ';margin-right:4px"></span>detected in some replicates &nbsp;&nbsp;' +
          '<span style="color:' + ACCENT + '">┆</span> limit of detection &nbsp;&nbsp;' +
          '<span style="opacity:.8">hover any depth for full statistics</span>',
      })
    );

    var orgs = rollupOrganisms(group);
    if (!orgs.length) {
      wrap.appendChild(el("div", { style: "color:" + MUTED }, ["No simulated organisms recovered."]));
      return wrap;
    }
    if (INSIL_ROLLUP === "Genus") {
      wrap.appendChild(
        el("div", { style: "font-size:.76em;color:" + MUTED + ";margin:-.3em 0 .6em" }, [
          "Rolled up from " + nativeLevel + " to Genus: reads and expected share add across members; " +
          "TASS is the strongest member; the genus counts as detected wherever any member was.",
        ])
      );
    }
    var grid = el("div", {
      style: "display:grid;grid-template-columns:repeat(auto-fill,minmax(330px,1fr));gap:.9em",
    });
    orgs.forEach(function (o) {
      var card = el("div", {
        class: "chart-wrap",
        style: "position:relative;border:1px solid #eee;border-radius:8px;padding:.6em .7em;background:#fff",
      });
      var head = el("div", {
        class: "chart-title",
        style: "display:flex;justify-content:space-between;align-items:baseline;gap:.5em",
      });
      head.appendChild(el("span", { style: "font-weight:600;color:#222" }, [o.name]));
      var lodTxt =
        o.lod_count == null
          ? "not detected"
          : "LoD " + kfmt(o.lod_count) + " " + (group.read_unit === "read pairs" ? "pr" : "rd");
      head.appendChild(pill(lodTxt, o.lod_count == null ? BAD : ACCENT));
      card.appendChild(head);
      var sub =
        (o.rolled && o.n_members > 1
          ? "Genus · " + o.n_members + " members"
          : (o.rolled ? "Genus" : nativeLevel) + " · taxid " + o.taxid) +
        " · " + o.category + " · expected " + (o.expected_fraction * 100).toFixed(1) + "% of pool";
      var subEl = el("div", { style: "font-size:.75em;color:" + MUTED + ";margin:1px 0 4px" }, [sub]);
      if (o.rolled && o.n_members > 1) {
        subEl.setAttribute("style", subEl.getAttribute("style") + ";cursor:help");
        subEl.setAttribute(
          "data-tt",
          encodeURIComponent(
            tipBody(o.name, "rolled up from " + o.n_members + " " + nativeLevel.toLowerCase() + "-level series",
              o.members.map(function (m, i) { return [m, "taxid " + o.member_taxids[i]]; }))
          )
        );
      }
      card.appendChild(subEl);
      card.appendChild(el("div", { class: "insil-scroll", html: organismChart(o, group) }));
      grid.appendChild(card);
    });
    wrap.appendChild(grid);
    return wrap;
  }

  // ── group ─────────────────────────────────────────────────────────────────
  function renderGroup(group) {
    var box = el("div", {
      style: "border:1px solid #e6e1f5;border-radius:10px;padding:1em 1.1em;margin-bottom:1.2em;background:#fff",
    });
    var plat =
      group.platform === "iss"
        ? "InSilicoSeq (Illumina)"
        : group.platform === "nanosim"
        ? "NanoSim (ONT)"
        : group.platform === "background"
        ? "Natural background (real reads)"
        : group.platform;
    var isBg = group.platform === "background" || group.source === "background";
    var head = el("div", { style: "display:flex;flex-wrap:wrap;align-items:center;gap:.6em;margin-bottom:.4em" });
    head.appendChild(el("span", { style: "font-size:1.05em;font-weight:700;color:" + ACCENT }, [group.parent]));
    head.appendChild(pill(plat, isBg ? "#00695c" : "#5c6bc0"));
    head.appendChild(pill(group.mode, group.mode === "consistent" ? "#00897b" : "#8e24aa"));
    head.appendChild(
      el("span", { style: "font-size:.8em;color:" + MUTED }, [
        group.n_datasets +
          " datasets · counts " +
          group.counts
            .map(function (c) { return kfmt(c); })
            .join(", ") +
          " " +
          (group.read_unit || "reads"),
      ])
    );
    box.appendChild(head);
    box.appendChild(renderDatasetTable(group));
    var charts = renderGroupCharts(group);
    if (charts) box.appendChild(charts);
    box.appendChild(renderOrganisms(group));
    attachTips(box);
    return box;
  }

  // ── entry point ───────────────────────────────────────────────────────────
  // Shared with 46_insilico_compare.js so a Detections row can be matched
  // against a genus-rolled series with exactly the same maths.
  window.insilicoRollupOrganisms = function (group, mode) {
    var prev = INSIL_ROLLUP;
    INSIL_ROLLUP = mode || "native";
    try { return rollupOrganisms(group); } finally { INSIL_ROLLUP = prev; }
  };

  window.drawInsilico = function drawInsilico() {
    var suite =
      typeof INSILICO_SUITE !== "undefined" ? INSILICO_SUITE : (window.HEATMAP_BOOT || {}).insilico_suite;
    var groupsHost = document.getElementById("insilico-groups");
    var emptyHost = document.getElementById("insilico-empty");
    if (!groupsHost) return;
    if (!suite || !suite.enabled || !(suite.groups && suite.groups.length)) {
      groupsHost.innerHTML = "";
      var pp = document.getElementById("insilico-params-panel");
      if (pp) pp.innerHTML = "";
      if (emptyHost) emptyHost.style.display = "";
      return;
    }
    ensureStyles();
    if (emptyHost) emptyHost.style.display = "none";
    renderParams(suite);
    groupsHost.innerHTML = "";
    suite.groups.forEach(function (g) {
      groupsHost.appendChild(renderGroup(g));
    });
  };

  // Reveal the tab when subsample data is present (BOOT is already parsed by the
  // time this end-of-body script runs).
  function _unhideInsilicoTab() {
    var has = typeof HAS_INSILICO !== "undefined" ? HAS_INSILICO : !!((window.HEATMAP_BOOT || {}).has_insilico_suite);
    var btn = document.getElementById("insilico-tab-btn");
    if (btn) btn.classList.toggle("hidden", !has);
  }
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", _unhideInsilicoTab);
  } else {
    _unhideInsilicoTab();
  }
})();
