/* ═══════════════════════════════════════════════════════════════════════════
       -  §  IN-SILICO SUITE TAB
       -     Renders the spike-in / dilution-series subsampling suite:
       -       (1) a run-parameter provenance panel,
       -       (2) per (parent × platform) group, a per-dataset "expected vs
       -           reality" table (target/observed reads, TP/FP/FN, P/R/F1), and
       -       (3) a per-organism dilution-series limit-of-detection view with a
       -           compact expected-vs-observed chart per organism.
       -     Data source: BOOT.insilico_suite (built server-side in make_report.py
       -     from the subsample datasets, which each flow through as their own
       -     sample). The tab is hidden unless --sim_subsample was enabled.
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  "use strict";

  var ACCENT = "#4527a0";
  var GOOD = "#2e7d32";
  var WARN = "#ef6c00";
  var BAD = "#c62828";
  var MUTED = "#777";

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

  function fmt(n) {
    if (n == null || n === "") return "—";
    if (typeof n === "number") return n.toLocaleString();
    return String(n);
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

  function renderParams(params) {
    var host = document.getElementById("insilico-params-panel");
    if (!host) return;
    host.innerHTML = "";
    host.setAttribute(
      "style",
      "background:#faf9ff;border:1px solid #e6e1f5;border-radius:8px;padding:1em 1.2em;margin-bottom:1.1em"
    );
    host.appendChild(
      el("div", {
        class: "chart-title",
        html: '<i class="fas fa-flask"></i> In-Silico Subsampling — parameters used',
      })
    );
    var grid = el("div", {
      style:
        "display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:.5em .9em;margin-top:.6em",
    });
    var keys = Object.keys(params || {});
    // Show known keys first (in label order), then any extras.
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
        style: "display:flex;flex-direction:column;padding:.35em .5em;background:#fff;border:1px solid #ece8f7;border-radius:6px",
      });
      cell.appendChild(el("span", { style: "font-size:.72em;color:" + MUTED + ";text-transform:uppercase;letter-spacing:.03em" }, [PARAM_LABELS[k] || k]));
      cell.appendChild(el("span", { style: "font-size:.95em;font-weight:600;color:#222;margin-top:2px" }, [v == null || v === "" ? "—" : String(v)]));
      grid.appendChild(cell);
    });
    host.appendChild(grid);
  }

  function pill(txt, color) {
    return el("span", {
      style:
        "display:inline-block;padding:1px 7px;border-radius:10px;font-size:.78em;font-weight:600;color:#fff;background:" +
        color,
    }, [txt]);
  }

  function f1Color(f1) {
    if (f1 >= 0.9) return GOOD;
    if (f1 >= 0.6) return WARN;
    return BAD;
  }

  function renderDatasetTable(group) {
    var unit = group.read_unit || "reads";
    var wrap = el("div", { style: "overflow-x:auto;margin:.4em 0 1em" });
    var t = el("table", {
      style: "border-collapse:collapse;width:100%;font-size:.85em;min-width:640px",
    });
    var heads = [
      "Dataset (rep)",
      "Target " + unit,
      "Actual " + unit,
      "Observed aligned (" + unit + ")",
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
        el("th", {
          style:
            "text-align:left;padding:.4em .6em;border-bottom:2px solid " +
            ACCENT +
            ";white-space:nowrap;font-size:.92em;color:#333",
        }, [h])
      );
    });
    thead.appendChild(hr);
    t.appendChild(thead);
    var tb = el("tbody");
    group.datasets.forEach(function (d, i) {
      var tr = el("tr", { style: i % 2 ? "background:#faf9ff" : "" });
      var shortId = "c" + (d.target_count) + " · r" + d.replicate;
      var cells = [
        el("span", { style: "font-family:monospace;font-size:.92em" }, [shortId]),
        fmt(d.target_count),
        fmt(d.actual_count),
        fmt(d.observed_total_reads),
        fmt(d.n_detected),
        String(d.tp),
        String(d.fp),
        String(d.fn),
        d.precision.toFixed(2),
        d.recall.toFixed(2),
        el("span", { style: "font-weight:700;color:" + f1Color(d.f1) }, [d.f1.toFixed(2)]),
      ];
      cells.forEach(function (c) {
        tr.appendChild(el("td", { style: "padding:.35em .6em;border-bottom:1px solid #eee;white-space:nowrap" }, [
          typeof c === "string" ? document.createTextNode(c) : c,
        ]));
      });
      tb.appendChild(tr);
    });
    t.appendChild(tb);
    wrap.appendChild(t);
    return wrap;
  }

  // Compact expected-vs-observed chart for one organism across the count series.
  // Inline SVG (no external chart lib). X = ordinal count positions; two bars
  // per count (expected outline, observed filled) on a shared log-ish scale, plus
  // a detection dot row beneath.
  function organismChart(series) {
    var W = Math.max(220, series.length * 64);
    var H = 96;
    var padL = 4, padR = 4, padT = 8, padB = 22;
    var innerW = W - padL - padR;
    var innerH = H - padT - padB;
    var maxV = 1;
    series.forEach(function (s) {
      maxV = Math.max(maxV, s.expected_reads, s.observed_reads);
    });
    var scale = function (v) {
      // sqrt scale keeps small dilutions visible next to large ones
      return (Math.sqrt(Math.max(0, v)) / Math.sqrt(maxV)) * innerH;
    };
    var band = innerW / series.length;
    var bw = Math.min(18, band * 0.32);
    var svg =
      '<svg viewBox="0 0 ' + W + " " + H + '" width="' + W + '" height="' + H + '" role="img">';
    series.forEach(function (s, i) {
      var cx = padL + band * (i + 0.5);
      var eh = scale(s.expected_reads);
      var oh = scale(s.observed_reads);
      var yb = padT + innerH;
      // expected (outline)
      svg +=
        '<rect x="' + (cx - bw - 1) + '" y="' + (yb - eh) + '" width="' + bw + '" height="' + eh +
        '" fill="none" stroke="' + ACCENT + '" stroke-width="1.2" rx="1"/>';
      // observed (filled)
      var oc = s.detected ? GOOD : BAD;
      svg +=
        '<rect x="' + (cx + 1) + '" y="' + (yb - oh) + '" width="' + bw + '" height="' + oh +
        '" fill="' + oc + '" opacity="0.85" rx="1"/>';
      // detection dot
      svg +=
        '<circle cx="' + cx + '" cy="' + (H - 12) + '" r="4" fill="' +
        (s.detected ? GOOD : "#ddd") + '" stroke="' + (s.detected ? GOOD : "#bbb") + '"/>';
      // count label
      svg +=
        '<text x="' + cx + '" y="' + (H - 2) + '" text-anchor="middle" font-size="9" fill="#666">' +
        (s.count >= 1000 ? s.count / 1000 + "k" : s.count) + "</text>";
    });
    svg += "</svg>";
    return svg;
  }

  function renderOrganisms(group) {
    var wrap = el("div", { style: "margin-top:.3em" });
    wrap.appendChild(
      el("div", { style: "font-weight:600;color:#333;margin:.2em 0 .5em;font-size:.95em" }, [
        "Per-organism dilution series — expected vs observed (limit of detection)",
      ])
    );
    // legend
    wrap.appendChild(
      el("div", { style: "font-size:.76em;color:" + MUTED + ";margin-bottom:.6em", html:
        '<span style="border:1.2px solid ' + ACCENT + ';padding:0 6px;border-radius:2px;margin-right:4px">&nbsp;</span>expected &nbsp;&nbsp;' +
        '<span style="background:' + GOOD + ';color:#fff;padding:0 6px;border-radius:2px;margin-right:4px">obs</span>detected &nbsp;&nbsp;' +
        '<span style="background:' + BAD + ';color:#fff;padding:0 6px;border-radius:2px;margin-right:4px">obs</span>below threshold' }));

    if (!group.organisms.length) {
      wrap.appendChild(el("div", { style: "color:" + MUTED }, ["No simulated organisms recovered."]));
      return wrap;
    }
    var grid = el("div", {
      style: "display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:.9em",
    });
    group.organisms.forEach(function (o) {
      var card = el("div", {
        style: "border:1px solid #eee;border-radius:8px;padding:.6em .7em;background:#fff",
      });
      var head = el("div", { style: "display:flex;justify-content:space-between;align-items:baseline;gap:.5em" });
      head.appendChild(el("span", { style: "font-weight:600;color:#222" }, [o.name]));
      var lodTxt = o.lod_count == null ? "not detected" : "LoD " + (o.lod_count >= 1000 ? o.lod_count / 1000 + "k" : o.lod_count) + " rd";
      head.appendChild(pill(lodTxt, o.lod_count == null ? BAD : ACCENT));
      card.appendChild(head);
      card.appendChild(
        el("div", { style: "font-size:.75em;color:" + MUTED + ";margin:1px 0 4px" }, [
          "taxid " + o.taxid + " · " + o.category + " · expected " + (o.expected_fraction * 100).toFixed(1) + "% of pool",
        ])
      );
      card.appendChild(el("div", { html: organismChart(o.series) }));
      card.appendChild(el("div", { style: "font-size:.7em;color:#999;margin-top:2px" }, ["read count →"]));
      grid.appendChild(card);
    });
    wrap.appendChild(grid);
    return wrap;
  }

  function renderGroup(group) {
    var box = el("div", {
      style: "border:1px solid #e6e1f5;border-radius:10px;padding:1em 1.1em;margin-bottom:1.2em;background:#fff",
    });
    var plat = group.platform === "iss" ? "InSilicoSeq (Illumina)"
      : group.platform === "nanosim" ? "NanoSim (ONT)"
      : group.platform === "background" ? "Natural background (real reads)"
      : group.platform;
    var isBg = group.platform === "background" || group.source === "background";
    var head = el("div", { style: "display:flex;flex-wrap:wrap;align-items:center;gap:.6em;margin-bottom:.4em" });
    head.appendChild(el("span", { style: "font-size:1.05em;font-weight:700;color:" + ACCENT }, [group.parent]));
    head.appendChild(pill(plat, isBg ? "#00695c" : "#5c6bc0"));
    head.appendChild(pill(group.mode, group.mode === "consistent" ? "#00897b" : "#8e24aa"));
    head.appendChild(el("span", { style: "font-size:.8em;color:" + MUTED }, [
      group.n_datasets + " datasets · counts " + group.counts.map(function (c) { return c >= 1000 ? c / 1000 + "k" : c; }).join(", ") +
      " " + (group.read_unit || "reads"),
    ]));
    box.appendChild(head);
    box.appendChild(renderDatasetTable(group));
    box.appendChild(renderOrganisms(group));
    return box;
  }

  window.drawInsilico = function drawInsilico() {
    var suite = typeof INSILICO_SUITE !== "undefined" ? INSILICO_SUITE : (window.HEATMAP_BOOT || {}).insilico_suite;
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
    if (emptyHost) emptyHost.style.display = "none";
    renderParams(suite.params || {});
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
