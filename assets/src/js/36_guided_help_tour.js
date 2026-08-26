/* ═══════════════════════════════════════════════════════════════════════════
       -  §  GUIDED HELP TOUR
       -     A right-aligned "Help & Tour" button in the banner opens a dimmed
       -     overlay with a single glowing card that walks through every tab —
       -     description, how-to-operate tips and a small schematic figure. The
       -     matching tab button is spotlit and the tab is activated so the real
       -     content shows behind the dim. Users step Back / Next, jump via the
       -     progress dots, or Skip out at any time (Esc / click backdrop also
       -     close). The tab that was active before the tour is restored on exit.
      ═══════════════════════════════════════════════════════════════════════════ */
(function initHelpTour() {
  // ── small schematic SVG figures (theme blues) ──
  const FIG = {
    overview: `<svg viewBox="0 0 200 120"><rect x="14" y="18" width="172" height="84" rx="8" fill="#fff" stroke="#bcd6f2"/><rect x="26" y="30" width="60" height="10" rx="3" fill="#1565c0"/><rect x="26" y="48" width="148" height="7" rx="3" fill="#cfe2f7"/><rect x="26" y="60" width="120" height="7" rx="3" fill="#cfe2f7"/><rect x="110" y="74" width="64" height="18" rx="4" fill="#e3f0fd"/><circle cx="158" cy="40" r="12" fill="#2196f3"/><path d="M153 40l4 4 7-8" stroke="#fff" stroke-width="2.4" fill="none" stroke-linecap="round" stroke-linejoin="round"/></svg>`,
    summary: `<svg viewBox="0 0 200 120"><rect x="46" y="16" width="108" height="92" rx="8" fill="#fff" stroke="#bcd6f2"/><rect x="78" y="10" width="44" height="14" rx="4" fill="#90b8e0"/>${[
      0, 1, 2, 3,
    ]
      .map(
        (r) =>
          `<circle cx="64" cy="${40 + r * 16}" r="5" fill="#2196f3"/><path d="M61 ${
            40 + r * 16
          }l2 2 4-4" stroke="#fff" stroke-width="1.6" fill="none"/><rect x="76" y="${36 + r * 16}" width="${
            64 - r * 8
          }" height="7" rx="3" fill="#cfe2f7"/>`,
      )
      .join("")}</svg>`,
    heatmap: `<svg viewBox="0 0 200 120">${Array.from({ length: 4 })
      .map((_, r) =>
        Array.from({ length: 6 })
          .map((__, c) => {
            const cols = ["#e3f2fd", "#90caf9", "#42a5f5", "#1565c0", "#0d47a1"];
            const v = (r * 6 + c) % 5;
            return `<rect x="${26 + c * 26}" y="${22 + r * 22}" width="22" height="18" rx="3" fill="${cols[v]}"/>`;
          })
          .join(""),
      )
      .join("")}</svg>`,
    tass: `<svg viewBox="0 0 200 120">${[60, 84, 48, 96, 72]
      .map(
        (h, i) =>
          `<rect x="${28 + i * 32}" y="${104 - h}" width="20" height="${h}" rx="3" fill="#1565c0" opacity="${
            0.55 + i * 0.09
          }"/>`,
      )
      .join("")}<line x1="18" y1="104" x2="186" y2="104" stroke="#9bb8d6" stroke-width="2"/></svg>`,
    sunburst: `<svg viewBox="0 0 200 120"><g transform="translate(100,62)"><circle r="16" fill="#1565c0"/><path d="M0 0 L0 -34 A34 34 0 0 1 29 -17 Z" fill="#42a5f5"/><path d="M0 0 L29 -17 A34 34 0 0 1 29 17 Z" fill="#90caf9"/><path d="M0 0 L29 17 A34 34 0 0 1 -29 17 Z" fill="#64b5f6"/><path d="M0 0 L-29 17 A34 34 0 0 1 0 -34 Z" fill="#2196f3"/><path d="M0 -34 L0 -52 A52 52 0 0 1 45 -26 Z" fill="#bbdefb"/><path d="M0 -52 A52 52 0 0 1 45 -26" fill="none"/></g></svg>`,
    coverage: `<svg viewBox="0 0 200 120"><line x1="22" y1="100" x2="186" y2="100" stroke="#9bb8d6" stroke-width="2"/><line x1="22" y1="100" x2="22" y2="16" stroke="#9bb8d6" stroke-width="2"/>${[
      [50, 40],
      [70, 55],
      [90, 38],
      [110, 70],
      [130, 50],
      [150, 78],
      [95, 60],
      [60, 72],
    ]
      .map(
        ([x, y], i) =>
          `<circle cx="${x}" cy="${y}" r="6" fill="${
            ["#1565c0", "#ff7f0e", "#2ca02c", "#d62728"][i % 4]
          }" opacity="0.85"/>`,
      )
      .join("")}</svg>`,
    proteins: `<svg viewBox="0 0 200 120"><path d="M70 20 C110 40 90 80 130 100 M130 20 C90 40 110 80 70 100" stroke="#1565c0" stroke-width="3" fill="none"/>${[
      28, 44, 60, 76,
    ]
      .map(
        (y, i) =>
          `<line x1="${78 + (i % 2) * 6}" y1="${y}" x2="${122 - (i % 2) * 6}" y2="${
            y + 6
          }" stroke="#42a5f5" stroke-width="3"/>`,
      )
      .join(
        "",
      )}<rect x="150" y="46" width="34" height="9" rx="3" fill="#d62728"/><rect x="150" y="60" width="24" height="9" rx="3" fill="#ff9800"/></svg>`,
    histogram: `<svg viewBox="0 0 200 120">${[20, 38, 62, 80, 70, 52, 34, 18]
      .map(
        (h, i) =>
          `<rect x="${24 + i * 20}" y="${104 - h}" width="16" height="${h}" fill="#2196f3" opacity="${
            0.5 + i * 0.05
          }"/>`,
      )
      .join("")}<line x1="18" y1="104" x2="186" y2="104" stroke="#9bb8d6" stroke-width="2"/></svg>`,
    explore: `<svg viewBox="0 0 200 120"><line x1="22" y1="100" x2="186" y2="100" stroke="#9bb8d6" stroke-width="2"/><polyline points="30,80 60,50 90,64 120,32 150,46 178,24" fill="none" stroke="#1565c0" stroke-width="3"/>${[
      [30, 80],
      [60, 50],
      [90, 64],
      [120, 32],
      [150, 46],
      [178, 24],
    ]
      .map(([x, y]) => `<circle cx="${x}" cy="${y}" r="4.5" fill="#ff7f0e"/>`)
      .join("")}</svg>`,
    table: `<svg viewBox="0 0 200 120"><rect x="20" y="22" width="160" height="78" rx="6" fill="#fff" stroke="#bcd6f2"/><rect x="20" y="22" width="160" height="18" fill="#1565c0"/>${[
      0, 1, 2, 3,
    ]
      .map((r) => `<line x1="20" y1="${40 + r * 15}" x2="180" y2="${40 + r * 15}" stroke="#e1ebf6"/>`)
      .join(
        "",
      )}<line x1="74" y1="22" x2="74" y2="100" stroke="#e1ebf6"/><line x1="128" y1="22" x2="128" y2="100" stroke="#e1ebf6"/></svg>`,
    map: `<svg viewBox="0 0 200 120"><rect x="24" y="20" width="152" height="80" rx="8" fill="#e8f3ff" stroke="#bcd6f2"/><path d="M40 70 q20 -18 44 -6 t52 -4" stroke="#a8c8ea" stroke-width="2" fill="none"/><path d="M100 36c-13 0-23 10-23 23 0 17 23 33 23 33s23-16 23-33c0-13-10-23-23-23z" fill="#d62728"/><circle cx="100" cy="59" r="8" fill="#fff"/></svg>`,
    runmeta: `<svg viewBox="0 0 200 120">${[
      [34, 30, "#1565c0"],
      [104, 30, "#2196f3"],
      [34, 58, "#43a047"],
      [104, 58, "#fb8c00"],
    ]
      .map(
        ([x, y, c]) =>
          `<rect x="${x}" y="${y}" width="62" height="20" rx="6" fill="${c}" opacity="0.85"/><circle cx="${
            x + 11
          }" cy="${y + 10}" r="3.5" fill="#fff"/>`,
      )
      .join("")}</svg>`,
    sidebar: `<svg viewBox="0 0 200 120"><rect x="118" y="14" width="70" height="92" rx="8" fill="#fff" stroke="#bcd6f2"/><rect x="128" y="24" width="50" height="9" rx="3" fill="#1565c0"/><rect x="128" y="40" width="50" height="8" rx="4" fill="#e3f0fd" stroke="#bcd6f2"/>${[
      0, 1, 2, 3,
    ]
      .map(
        (r) =>
          `<circle cx="134" cy="${64 + r * 12}" r="5" fill="${
            ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"][r]
          }"/><rect x="144" y="${60 + r * 12}" width="34" height="7" rx="3" fill="#cfe2f7"/>`,
      )
      .join("")}<rect x="14" y="14" width="92" height="92" rx="8" fill="#eef5fc" stroke="#dbe8f5"/></svg>`,
  };

  // ── tour content. tab:'x' targets a .tab-btn[data-tab=x]; sel targets an element. ──
  const STEPS = [
    {
      badge: "Welcome",
      icon: "fa-wand-magic-sparkles",
      title: "Welcome to TaxTriage",
      fig: FIG.overview,
      desc: "This report combines every sample from your run into one interactive view. The tabs along the top each focus on a different question — from a quick overview to deep per-genome coverage. This quick tour shows what each tab does and how to drive it.",
      tips: [
        "Use Back / Next or the dots below to move through the tour.",
        "Hit Skip tour (top-right) any time — you'll return to where you started.",
      ],
    },
    {
      tab: "summary",
      icon: "fa-clipboard-list",
      title: "Summary",
      fig: FIG.summary,
      desc: "The landing view: headline KPIs plus sub-tabs for Detections, Genera, VF/AMR and Cross-Sample patterns. Start here for a run-at-a-glance read on what was found and how confident the calls are.",
      tips: [
        "Switch the inner sub-tabs to pivot between detections, genera and cross-sample views.",
        "KPI cards summarise totals and the recommended TASS cutoff.",
      ],
    },
    {
      tab: "heatmap",
      icon: "fa-th",
      title: "Heatmap",
      fig: FIG.heatmap,
      desc: "A samples × organisms grid where cell colour encodes the selected metric (TASS, reads, coverage…). It's the fastest way to spot which organisms are shared across samples and which are unique.",
      tips: [
        "Pick the metric and scaling from the controls above the grid.",
        "Hover any cell for the exact value; darker = higher.",
      ],
    },
    {
      tab: "tass",
      icon: "fa-chart-bar",
      title: "TASS Comparison",
      fig: FIG.tass,
      desc: "Compares TASS confidence scores side-by-side so you can rank detections and see how each sits against the per-sample-type cutoff.",
      tips: [
        "Bars are coloured by sample using the same colours as the right panel.",
        "Use it to justify which calls clear the confidence threshold.",
      ],
    },
    {
      tab: "sunburst",
      icon: "fa-circle-nodes",
      title: "Sunburst",
      fig: FIG.sunburst,
      desc: "A radial taxonomy view — inner rings are higher ranks (kingdom→genus), outer rings are species/strains, sized by the metric you choose. Great for seeing taxonomic composition at a glance.",
      tips: [
        "Click a wedge to zoom into that branch; click the centre to zoom out.",
        "Switch the sizing metric from the controls above.",
      ],
    },
    {
      tab: "coverage",
      icon: "fa-layer-group",
      title: "Coverage",
      fig: FIG.coverage,
      desc: "Plots coverage against TASS or read count so outliers jump out — e.g. high reads but low breadth. Each point is a detection, coloured by its sample.",
      tips: [
        "Point colours match the sample colours in the right panel.",
        "Open a detection's coverage comparison to overlay genome-position profiles.",
      ],
    },
    {
      tab: "proteins",
      icon: "fa-dna",
      title: "VF / AMR",
      fig: FIG.proteins,
      desc: "Virulence-factor, AMR and transporter gene hits, with charts plus a searchable table. Appears only when protein annotations are present in the run.",
      tips: ["Click a bar to filter the table to those gene hits.", "Click a table row for the per-gene detail panel."],
    },
    {
      tab: "histogram",
      icon: "fa-chart-column",
      title: "Histograms",
      fig: FIG.histogram,
      desc: "Per-contig / per-assembly read-distribution histograms, showing how depth is spread across each genome. Appears when per-contig data is available.",
      tips: [
        "Pick the sample from the selector to redraw its distributions.",
        "A warning icon flags assemblies with uneven or sparse coverage.",
      ],
    },
    {
      tab: "explore",
      icon: "fa-magnifying-glass-chart",
      title: "Explore",
      fig: FIG.explore,
      desc: "Cross-sample exploratory charts for spotting multi-metric patterns and trends across the whole run — the place to slice the data freely.",
      tips: [
        "Choose axes/metrics from the controls to build your own comparison.",
        "Use it to hunt for batch effects or sample-type trends.",
      ],
    },
    {
      tab: "table",
      icon: "fa-table",
      title: "Table",
      fig: FIG.table,
      desc: "The full detections table behind every chart. Sort, search and paginate the raw rows, and pin detections to compare their coverage profiles.",
      tips: [
        "Click a column header to sort; use the sidebar search to filter.",
        "Pin rows, then open Compare to overlay their coverage profiles.",
      ],
    },
    {
      tab: "runmeta",
      icon: "fa-table-list",
      title: "Metadata",
      fig: FIG.runmeta,
      desc: "Per-sample run details — type, platform and any uploaded metadata fields — in one editable table. The <b>Group by</b> bar underneath is shared: whatever you pick here also drives the Mapping and Trends tabs.",
      tips: [
        "The View button on a sample jumps straight here and highlights its row.",
        "Upload a metadata CSV to enrich these fields.",
        "Add latitude / longitude to light up Mapping; add collection_time for Trends.",
      ],
    },
    {
      tab: "map",
      icon: "fa-map-location-dot",
      title: "Mapping",
      fig: FIG.runmeta,
      desc: "Samples plotted at their precise coordinates, or aggregated by country / state. Draw a lasso or a rectangle to cut out an area and get a summary for it — draw several and they are compared side by side.",
      tips: [
        "Regions you draw stay put: the Show regions switch hides them without deleting them.",
        "Group by region turns your drawn areas into a metadata column the other views can use.",
      ],
    },
    {
      tab: "trends",
      icon: "fa-chart-line",
      title: "Trends & group analysis",
      fig: FIG.runmeta,
      desc: "Longitudinal change over time, host & disease breakdowns, and the group heatmap / network / cross-entry comparison — all reading the shared Group by selection.",
      tips: ["Group by columns are defined in the Metadata tab; the bar here mirrors that selection."],
    },
    {
      sel: "#sidebar",
      icon: "fa-sliders",
      title: "The right panel — filters & sample colours",
      fig: FIG.sidebar,
      desc: "This panel drives the entire report. Search by sample or organism (regex), set a minimum TASS cutoff, and show/hide samples. Crucially, the colour swatch next to each sample is that sample's colour everywhere — heatmap, charts and coverage profiles all follow it.",
      tips: [
        "Click a sample's swatch to recolour it consistently across every tab.",
        "Toggle a sample off to drop it from all views at once.",
        "Export Report PDF (top of the panel) captures the current filtered state.",
      ],
    },
    {
      badge: "You're set",
      icon: "fa-circle-check",
      title: "That's the tour!",
      fig: FIG.overview,
      desc: "You've seen every tab and how the right panel ties them together. Re-open this tour any time from the Help & Tour button in the top-right corner.",
      tips: [
        "Hover icons and chart controls for inline tooltips with more detail.",
        "Filters and sample colours apply across all tabs simultaneously.",
      ],
    },
  ];

  const overlay = document.getElementById("help-overlay");
  if (!overlay) return;
  const elBadge = document.getElementById("help-step-badge");
  const elFig = document.getElementById("help-figure");
  const elIcon = document.getElementById("help-card-icon");
  const elTitle = document.getElementById("help-card-title");
  const elDesc = document.getElementById("help-card-desc");
  const elTips = document.getElementById("help-card-tips");
  const elDots = document.getElementById("help-dots");
  const elProg = document.getElementById("help-progress-text");
  const btnPrev = document.getElementById("help-prev");
  const btnNext = document.getElementById("help-next");
  const btnSkip = document.getElementById("help-skip");
  const card = document.getElementById("help-card");

  let steps = []; // active steps after filtering hidden tabs
  let idx = 0;
  let returnTab = null;

  function tabBtn(tab) {
    return document.querySelector(`.tab-btn[data-tab="${tab}"]`);
  }
  function tabAvailable(tab) {
    const b = tabBtn(tab);
    return b && !b.classList.contains("hidden") && b.offsetParent !== null;
  }
  function buildSteps() {
    return STEPS.filter((s) => {
      if (s.tab) return tabAvailable(s.tab);
      if (s.sel) return !!document.querySelector(s.sel);
      return true;
    });
  }
  function clearSpot() {
    document.querySelectorAll(".tab-btn.help-spotlight").forEach((b) => b.classList.remove("help-spotlight"));
    document.querySelectorAll(".help-spot-el").forEach((e) => e.classList.remove("help-spot-el"));
  }
  function gotoTab(tab) {
    const b = tabBtn(tab);
    if (!b || b.classList.contains("hidden")) return null;
    if (!b.classList.contains("active")) {
      try {
        b.click();
      } catch (e) {}
    }
    return b;
  }
  function render() {
    const s = steps[idx];
    if (!s) return;
    clearSpot();
    // spotlight + activate target
    if (s.tab) {
      const b = gotoTab(s.tab);
      if (b) b.classList.add("help-spotlight");
    } else if (s.sel) {
      const t = document.querySelector(s.sel);
      if (t) t.classList.add("help-spot-el");
    }
    elBadge.textContent = s.badge || (s.tab ? "Tab guide" : "Guide");
    elFig.innerHTML = s.fig || "";
    elIcon.innerHTML = `<i class="fas ${s.icon || "fa-circle-question"}"></i>`;
    elTitle.textContent = s.title || "";
    elDesc.textContent = s.desc || "";
    elTips.innerHTML = (s.tips || []).map((t) => `<li>${t}</li>`).join("");
    // dots
    elDots.innerHTML = steps
      .map((_, i) => `<span class="help-dot${i === idx ? " active" : ""}" data-i="${i}"></span>`)
      .join("");
    elProg.textContent = `${idx + 1} / ${steps.length}`;
    btnPrev.disabled = idx === 0;
    const last = idx === steps.length - 1;
    btnNext.innerHTML = last ? `Done <i class="fas fa-check"></i>` : `Next <i class="fas fa-arrow-right"></i>`;
    card.scrollTop = 0;
  }
  function show(i) {
    idx = Math.max(0, Math.min(steps.length - 1, i));
    render();
  }
  function open() {
    steps = buildSteps();
    if (!steps.length) return;
    const cur = document.querySelector(".tab-btn.active");
    returnTab = cur ? cur.dataset.tab : null;
    document.body.classList.add("help-active");
    overlay.classList.add("open");
    overlay.setAttribute("aria-hidden", "false");
    idx = 0;
    render();
    card.focus && card.focus();
  }
  function close() {
    clearSpot();
    overlay.classList.remove("open");
    overlay.setAttribute("aria-hidden", "true");
    document.body.classList.remove("help-active");
    if (returnTab) gotoTab(returnTab);
  }

  // ── Pipeline info button: hover tooltip ───────────────────────────────
  (function initPipelineInfo() {
    const btn = document.getElementById("pipeline-info-btn");
    if (!btn) return;

    // Prefer global BOOT fields (set by make_report.py); fall back to
    // scanning per-sample SAMPLE_META for older reports that lack them.
    const _naVal = (v) => !v || v === "NA" || v === "null" || v === "NULL" || v === "none";

    let revisionDisplay,
      commits = [];
    const _bootRev = (typeof BOOT !== "undefined" && BOOT.pipeline_revision) || null;
    const _bootCommit = (typeof BOOT !== "undefined" && BOOT.pipeline_commit) || null;

    // Helper: is a git hash (7–40 hex chars)?
    const _isHash = (v) => v && /^[0-9a-f]{7,40}$/i.test(v);
    // Format commit: clickable link for real hashes, plain text for "local"/short ids
    const _fmtCommit = (v) =>
      _isHash(v)
        ? `<a href="https://github.com/jhuapl-bio/taxtriage/commit/${v}" ` +
          `target="_blank" style="color:#90caf9">${v}</a>`
        : v || "";

    if (_bootRev !== null && _bootRev !== undefined) {
      // Global value from BOOT (preferred path — always present in new reports)
      revisionDisplay = _bootRev && !_naVal(_bootRev) ? _bootRev : "Not Specified or Local Build";
      if (_bootCommit && !_naVal(_bootCommit)) commits = [_bootCommit];
    } else {
      // Fallback: scan SAMPLE_META (older reports / uploaded JSON)
      const allRevisions = [
        ...new Set(
          Object.values(SAMPLE_META)
            .map((m) => m.workflow_revision)
            .filter((v) => v && !_naVal(v)),
        ),
      ];
      const allCommits = [
        ...new Set(
          Object.values(SAMPLE_META)
            .map((m) => m.commit_id)
            .filter((v) => v && !_naVal(v)),
        ),
      ];
      revisionDisplay = allRevisions.length ? allRevisions.join(", ") : "Not Specified or Local Build";
      commits = allCommits;
    }
    const platforms = [
      ...new Set(
        Object.values(SAMPLE_META)
          .map((m) => m.platform)
          .filter((v) => v && v !== "unknown"),
      ),
    ];
    const generatedAt = (typeof BOOT !== "undefined" && BOOT.report_generated_at) || null;

    const _fmtDate = (iso) => {
      try {
        return new Date(iso).toLocaleString(undefined, {
          year: "numeric",
          month: "short",
          day: "numeric",
          hour: "2-digit",
          minute: "2-digit",
          timeZoneName: "short",
        });
      } catch (e) {
        return iso;
      }
    };

    const row = (label, val, color) =>
      `<tr><td style="padding:2px 10px 2px 0;color:${color || "#90caf9"};white-space:nowrap">${label}</td>` +
      `<td style="padding:2px 0;color:#e0e0e0;font-family:monospace;font-size:0.95em">${val}</td></tr>`;

    let tip =
      `<b style="color:#90caf9">TaxTriage Run Info</b><br>` +
      `<table style="margin-top:6px;border-collapse:collapse;min-width:220px">`;
    tip += row("Branch / revision", revisionDisplay || "—");
    if (commits.length) tip += row("Commit", commits.map(_fmtCommit).join(", "));
    tip += row("Samples", Object.keys(SAMPLE_META).length || "—", "#ffd580");
    if (platforms.length) tip += row("Platform(s)", platforms.join(", "), "#a5d6a7");
    tip += row("Report built", generatedAt ? _fmtDate(generatedAt) : "—", "#ce93d8");
    tip += `</table>`;

    btn.addEventListener("mouseover", (ev) => showTip(tip, ev));
    btn.addEventListener("mousemove", moveTip);
    btn.addEventListener("mouseout", hideTip);
  })();

  document.getElementById("help-btn")?.addEventListener("click", open);
  btnSkip?.addEventListener("click", close);
  btnPrev?.addEventListener("click", () => show(idx - 1));
  btnNext?.addEventListener("click", () => (idx === steps.length - 1 ? close() : show(idx + 1)));
  elDots?.addEventListener("click", (e) => {
    const d = e.target.closest(".help-dot");
    if (d) show(parseInt(d.dataset.i, 10));
  });
  document.getElementById("help-dim")?.addEventListener("click", close);
  document.addEventListener("keydown", (e) => {
    if (!overlay.classList.contains("open")) return;
    if (e.key === "Escape") close();
    else if (e.key === "ArrowRight") idx === steps.length - 1 ? close() : show(idx + 1);
    else if (e.key === "ArrowLeft") show(idx - 1);
  });

  /* ══════════════════════════════════════════════════════════════════════
     SPOTLIGHT WALKTHROUGH ENGINE
        A reusable stepped highlighter modelled on the attached report: a dim
        backdrop, a glowing ring that moves to each target (measured live from
        getBoundingClientRect) and a compact popover pinned next to it with
        Back / Next / dots. Used by per-tab Help mode and the right-panel "?"
        button. Leaves the little bottom-left #tab-help-panel visible on top.
     ══════════════════════════════════════════════════════════════════════ */
  const Spotlight = (function () {
    let dim, ring, pop, popIco, popTitle, popStep, popBody, popTips, popDots, btnBack, btnNext, btnClose;
    let steps = [],
      idx = 0,
      target = null,
      onExit = null,
      built = false,
      reposition = null,
      settleTimer = null;

    function build() {
      if (built) return;
      built = true;
      dim = document.createElement("div");
      dim.id = "hs-dim";
      ring = document.createElement("div");
      ring.id = "hs-ring";
      pop = document.createElement("div");
      pop.id = "hs-pop";
      pop.setAttribute("role", "dialog");
      pop.innerHTML =
        '<div id="hs-pop-head"><span id="hs-pop-ico"><i class="fas fa-circle-question"></i></span>' +
        '<h4 id="hs-pop-title"></h4><span id="hs-pop-step"></span>' +
        '<button id="hs-pop-x" type="button" title="Close" aria-label="Close help">&times;</button></div>' +
        '<div id="hs-pop-body"></div><ul id="hs-pop-tips"></ul>' +
        '<div id="hs-pop-foot"><div id="hs-pop-dots"></div><span class="hs-spacer"></span>' +
        '<button id="hs-back" type="button"><i class="fas fa-arrow-left"></i> Back</button>' +
        '<button id="hs-next" type="button">Next <i class="fas fa-arrow-right"></i></button></div>';
      document.body.appendChild(dim);
      document.body.appendChild(ring);
      document.body.appendChild(pop);
      popIco = pop.querySelector("#hs-pop-ico");
      popTitle = pop.querySelector("#hs-pop-title");
      popStep = pop.querySelector("#hs-pop-step");
      popBody = pop.querySelector("#hs-pop-body");
      popTips = pop.querySelector("#hs-pop-tips");
      popDots = pop.querySelector("#hs-pop-dots");
      btnBack = pop.querySelector("#hs-back");
      btnNext = pop.querySelector("#hs-next");
      btnClose = pop.querySelector("#hs-pop-x");
      btnBack.addEventListener("click", () => go(idx - 1));
      btnNext.addEventListener("click", () => (idx >= steps.length - 1 ? stop() : go(idx + 1)));
      btnClose.addEventListener("click", stop);
      dim.addEventListener("click", stop);
      popDots.addEventListener("click", (e) => {
        const d = e.target.closest(".hs-dot");
        if (d) go(parseInt(d.dataset.i, 10));
      });
      document.addEventListener("keydown", (e) => {
        if (!isOpen()) return;
        if (e.key === "Escape") stop();
        else if (e.key === "ArrowRight") idx >= steps.length - 1 ? stop() : go(idx + 1);
        else if (e.key === "ArrowLeft") go(idx - 1);
      });
    }

    function isOpen() {
      return built && dim.classList.contains("show");
    }
    // first selector that exists (visible preferred); sel may be string or array
    function resolve(sel) {
      if (!sel) return null;
      const list = Array.isArray(sel) ? sel : [sel];
      for (const s of list) {
        const el = document.querySelector(s);
        if (el && el.offsetParent !== null) return el;
      }
      for (const s of list) {
        const el = document.querySelector(s);
        if (el) return el;
      }
      return null;
    }
    function clearTarget() {
      if (target) {
        target.classList.remove("hs-target");
        target = null;
      }
    }
    function place() {
      if (!target) {
        ring.style.display = "none";
        return;
      }
      const r = target.getBoundingClientRect();
      const pad = 6;
      ring.style.display = "block";
      ring.style.top = r.top - pad + "px";
      ring.style.left = r.left - pad + "px";
      ring.style.width = r.width + pad * 2 + "px";
      ring.style.height = r.height + pad * 2 + "px";
      const pw = pop.offsetWidth || 320,
        ph = pop.offsetHeight || 210,
        gap = 16,
        vw = window.innerWidth,
        vh = window.innerHeight;
      let left,
        top,
        beside = true;
      if (r.right + gap + pw <= vw) left = r.right + gap; // right of target
      else if (r.left - gap - pw >= 0) left = r.left - gap - pw; // left of target
      else {
        beside = false;
        left = Math.max(8, Math.min(vw - pw - 8, r.left));
      }
      if (beside) top = r.top + r.height / 2 - ph / 2;
      else top = r.bottom + gap + ph <= vh ? r.bottom + gap : r.top - gap - ph;
      top = Math.max(8, Math.min(vh - ph - 8, top));
      pop.style.left = left + "px";
      pop.style.top = top + "px";
    }
    function render() {
      const s = steps[idx];
      if (!s) return;
      clearTarget();
      target = resolve(s.sel);
      if (target) {
        target.classList.add("hs-target");
        try {
          target.scrollIntoView({ behavior: "smooth", block: "center", inline: "nearest" });
        } catch (e) {}
      }
      popIco.innerHTML = `<i class="fas ${s.icon || "fa-circle-dot"}"></i>`;
      popTitle.textContent = s.title || "";
      popStep.textContent = `${idx + 1} / ${steps.length}`;
      popBody.innerHTML = s.desc || "";
      popTips.innerHTML = (s.tips || []).map((t) => `<li>${t}</li>`).join("");
      popDots.innerHTML = steps
        .map((_, i) => `<span class="hs-dot${i === idx ? " on" : ""}" data-i="${i}"></span>`)
        .join("");
      btnBack.disabled = idx === 0;
      btnNext.innerHTML =
        idx >= steps.length - 1 ? 'Done <i class="fas fa-check"></i>' : 'Next <i class="fas fa-arrow-right"></i>';
      requestAnimationFrame(place);
      clearTimeout(settleTimer);
      settleTimer = setTimeout(place, 360); // re-measure after the smooth scroll settles
    }
    function go(i) {
      idx = Math.max(0, Math.min(steps.length - 1, i));
      render();
    }
    function start(stepList, opts) {
      opts = opts || {};
      build();
      const resolved = (stepList || []).filter((s) => resolve(s.sel));
      if (!resolved.length) return false;
      steps = resolved;
      idx = 0;
      onExit = opts.onExit || null;
      document.body.classList.add("hs-open");
      dim.classList.add("show");
      ring.classList.add("show");
      pop.classList.add("show");
      if (!reposition) {
        reposition = () => {
          if (isOpen()) place();
        };
        window.addEventListener("resize", reposition);
        window.addEventListener("scroll", reposition, true);
      }
      render();
      return true;
    }
    function stop() {
      if (!built) return;
      clearTarget();
      dim.classList.remove("show");
      ring.classList.remove("show");
      pop.classList.remove("show");
      document.body.classList.remove("hs-open");
      const cb = onExit;
      onExit = null;
      if (cb) cb();
    }
    return { start, stop, isOpen };
  })();

  // ── Major-area step lists for each tab (the per-tab Help-mode walkthrough) ──
  const TAB_AREAS = {
    summary: [
      {
        sel: "#summary-kpi-row",
        icon: "fa-gauge-high",
        title: "Headline KPIs",
        desc: "At-a-glance totals for the run — samples, detections and the recommended TASS cutoff.",
        tips: ["These recompute live as you filter samples in the right panel."],
      },
      {
        sel: "#sum-inner-tabs",
        icon: "fa-folder-tree",
        title: "Summary sub-views",
        desc: "Switch between <b>Detections</b>, <b>Genera</b>, <b>VF/AMR</b> and <b>Cross-Sample</b> without leaving this tab.",
        tips: ["Each sub-view answers a different first-look question."],
      },
      {
        sel: "#watch-panel",
        icon: "fa-star",
        title: "Follow-up watchlist",
        desc: "Star organisms here to flag them; the star markers then appear across every other tab.",
        tips: ["Use it to build a shortlist while you triage."],
      },
    ],
    heatmap: [
      {
        sel: "#pane-heatmap .tab-controls",
        icon: "fa-sliders",
        title: "Metric & scaling",
        desc: "Choose the value shown (TASS, reads, coverage…), the taxonomic rank, colour scaling and whether cells print their value.",
        tips: ["Log / sqrt scaling helps when a few cells dominate."],
      },
      {
        sel: "#heatmap-svg-wrap",
        icon: "fa-th",
        title: "Samples × organisms grid",
        desc: "Each cell is one organism in one sample; darker = higher. The fastest way to spot shared vs unique organisms.",
        tips: ["Hover any cell for the exact value."],
      },
    ],
    tass: [
      {
        sel: "#pane-tass .tab-controls",
        icon: "fa-sliders",
        title: "Comparison controls",
        desc: "Set the rank, colour scaling, taxonomic level and whether the per-sample-type cutoff line is drawn.",
        tips: ["Toggle the cutoff to see which calls clear the threshold."],
      },
      {
        sel: "#tass-svg-wrap",
        icon: "fa-chart-bar",
        title: "TASS chart",
        desc: "Confidence scores side-by-side, coloured by sample using the same colours as the right panel.",
        tips: [],
      },
    ],
    sunburst: [
      {
        sel: "#sun-metric",
        icon: "fa-sliders",
        title: "Sizing metric",
        desc: "Pick which metric sizes the wedges — reads, TASS, coverage and so on.",
        tips: [],
      },
      {
        sel: "#sun-panels-container",
        icon: "fa-circle-nodes",
        title: "Radial taxonomy",
        desc: "Inner rings are higher ranks (kingdom→genus), outer rings species/strains. Click a wedge to zoom in; click the centre to zoom out.",
        tips: ["Add a panel to compare two samples side-by-side."],
      },
    ],
    coverage: [
      {
        sel: "#pane-coverage .tab-controls",
        icon: "fa-sliders",
        title: "Axes & scaling",
        desc: "Choose what the X, Y and bubble-size axes represent, plus the colour scaling.",
        tips: ["High reads + low breadth outliers jump out here."],
      },
      {
        sel: "#coverage-svg-wrap",
        icon: "fa-layer-group",
        title: "Coverage scatter",
        desc: "Each point is a detection, coloured by its sample. Open a point's comparison to overlay genome-position profiles.",
        tips: [],
      },
    ],
    proteins: [
      {
        sel: "#prot-genus-wrap",
        icon: "fa-dna",
        title: "VF / AMR by genus",
        desc: "Virulence-factor, AMR and transporter gene hits summarised as charts.",
        tips: ["Click a bar to filter the table below to those hits."],
      },
      {
        sel: "#prot-table-wrap",
        icon: "fa-table",
        title: "Gene table",
        desc: "The searchable per-gene table. Click a row for the detailed panel.",
        tips: [],
      },
    ],
    histogram: [
      {
        sel: "#hist-controls",
        icon: "fa-sliders",
        title: "Sample selector",
        desc: "Pick the sample to redraw its per-contig / per-assembly read-distribution histograms.",
        tips: ["A warning icon flags uneven or sparse coverage."],
      },
      {
        sel: "#pane-histogram",
        icon: "fa-chart-column",
        title: "Read distributions",
        desc: "Shows how sequencing depth is spread across each genome.",
        tips: [],
      },
    ],
    explore: [
      {
        sel: "#pane-explore .tab-controls",
        icon: "fa-sliders",
        title: "Build a comparison",
        desc: "Choose the colour-by, size-by and scaling to slice the run across any metrics you like.",
        tips: ["Good for hunting batch effects or sample-type trends."],
      },
      {
        sel: "#explore-bubble-wrap",
        icon: "fa-magnifying-glass-chart",
        title: "Exploratory charts",
        desc: "Cross-sample bubble and radar views for spotting multi-metric patterns.",
        tips: [],
      },
    ],
    table: [
      {
        sel: "#tbl-toolbar",
        icon: "fa-sliders",
        title: "Table toolbar",
        desc: "Search, choose columns and export the raw detections that sit behind every chart.",
        tips: ["Sort by clicking a column header."],
      },
      {
        sel: "#pane-table",
        icon: "fa-table",
        title: "Detections table",
        desc: "Every detection row. Pin rows, then open Compare to overlay their coverage profiles.",
        tips: [],
      },
    ],
    runmeta: [
      {
        sel: "#runmeta-toolbar",
        icon: "fa-sliders",
        title: "Metadata toolbar",
        desc: "Edit per-sample details and upload a metadata CSV to enrich the fields.",
        tips: [],
      },
      {
        sel: '.mg-bar-mount[data-mg-mount="meta"]',
        icon: "fa-layer-group",
        title: "Group by (shared)",
        desc: "Pick one or more metadata columns to bucket samples into groups. The Mapping and Trends tabs mirror this exact selection.",
        tips: ["Two or more columns combine as a union — one group per observed combination."],
      },
    ],
    map: [
      {
        sel: "#mapdraw-bar",
        icon: "fa-draw-polygon",
        title: "Draw regions",
        desc: "Arm the lasso or rectangle, trace an area on the map, and a summary card for that region appears underneath.",
        tips: [
          "Show regions hides them without deleting them.",
          "Group by region publishes the regions as a metadata column.",
        ],
      },
      {
        sel: "#geo-cmp-level",
        icon: "fa-earth-americas",
        title: "Precise vs aggregated",
        desc: "Switch between per-sample pins and a country / state choropleth of the selected metric.",
        tips: [],
      },
    ],
    trends: [
      {
        sel: "#meta-subtabs",
        icon: "fa-chart-line",
        title: "Trends sub-views",
        desc: "Longitudinal, host & disease, group heatmap, group network and cross-entry comparison.",
        tips: ["Sub-tabs grey out until the metadata they need is present."],
      },
    ],
  };

  // ── Right-panel walkthrough (the inline "?" button next to Export Report) ──
  const SIDEBAR_STEPS = [
    {
      sel: ["#report-pdf-btn", "#sidebar h3"],
      icon: "fa-file-pdf",
      title: "Filters & Export",
      desc: "This right panel drives the entire report. <b>Export Report PDF</b> captures the current filtered state as a printable layout.",
      tips: ["Every filter you set here applies across all tabs at once."],
    },
    {
      sel: ["#filter-search-box", "#filter-text"],
      icon: "fa-magnifying-glass",
      title: "Search",
      desc: "Filter by <b>Sample</b> or <b>Organism</b>. Plain words do a substring match; power users can use regex.",
      tips: ["Use the options toggle to set the scope and case-sensitivity."],
    },
    {
      sel: ["#per-type-tass-wrap", "#filter-min-wrap"],
      icon: "fa-crosshairs",
      title: "Confidence cutoff",
      desc: "Set the minimum TASS score. Detections below the cutoff drop out of every view.",
      tips: ["Per-sample-type cutoffs can be set independently."],
    },
    {
      sel: ["#filter-kingdom-wrap", "#filter-mc", "#filter-mt-dna"],
      icon: "fa-filter",
      title: "More filters",
      desc: "Narrow by kingdom, molecular type (DNA / RNA), high-confidence only and other toggles.",
      tips: [],
    },
    {
      sel: ["#toggle-all-samples-btn", "#sidebar"],
      icon: "fa-palette",
      title: "Samples & colours",
      desc: "Show or hide individual samples, and click a sample's colour swatch to recolour it — that colour then follows the sample across the heatmap, charts and coverage profiles.",
      tips: ["Toggle a sample off to drop it from every view at once."],
    },
    {
      sel: ["#specimen-merge-bar", "#specimen-merge-toggle"],
      icon: "fa-object-group",
      title: "Group samples into a specimen",
      desc: "Turn on <b>Specimen</b> mode to merge related samples (e.g. replicates or timepoints) into one specimen. <b>Group selected</b> combines the ones you tick; <b>Combine all</b> rolls everything into a single specimen.",
      tips: [
        "Grouped samples share a row and colour across every tab.",
        "The count shows how many specimens are currently defined.",
      ],
    },
    {
      sel: ["#sample-list", "#sample-list-tools"],
      icon: "fa-pen",
      title: "Edit & rename samples",
      desc: 'Each sample row has a <i class="fas fa-pen"></i> pencil to rename it and controls to attach files. Renames and edits flow through to every chart, table and export.',
      tips: [
        "Use the search box above the list to find a sample fast.",
        "Page through long sample lists with the ‹ › arrows.",
      ],
    },
    {
      sel: ["#upload-zone", "#file-upload-input"],
      icon: "fa-file-arrow-up",
      title: "Load run data (JSON)",
      desc: "Drop or browse for <code>all.samples.json</code>, <code>.paths.json</code>, <code>.xlsx</code>, <code>.tsv</code> or <code>.txt</code> files to load a run into this report. A metadata CSV can be added just below to enrich sample fields.",
      tips: [
        "Loaded data merges into the current view — filters still apply.",
        "Clear uploaded data to return to the bundled run.",
      ],
    },
    {
      sel: ["#state-io-row", "#state-load-input"],
      icon: "fa-floppy-disk",
      title: "Save / load state",
      desc: "<b>Export State</b> saves a self-contained snapshot of the current data <em>and</em> every filter, grouping and colour choice as a JSON file. <b>Load State</b> reopens one to return to exactly this view.",
      tips: [
        "Great for sharing an annotated view with a colleague.",
        "State files reload offline — no pipeline rerun needed.",
      ],
    },
  ];

  /* ── Per-tab Help mode ───────────────────────────────────────────────
           A toggle that pins a small contextual panel showing help for whatever
           tab the user is on, and runs the stepped Spotlight walkthrough over
           that tab's major areas. The right-panel "?" button reuses the same
           machinery, scoped to the sidebar. */
  (function initHelpMode() {
    const btn = document.getElementById("helpmode-btn");
    const panel = document.getElementById("tab-help-panel");
    if (!btn || !panel) return;
    const sideBtn = document.getElementById("sidebar-help-btn");
    const elFig = document.getElementById("tab-help-fig");
    const elGroup = document.getElementById("tab-help-group");
    const elEyebrow = document.getElementById("tab-help-eyebrow");
    const elTitle = document.getElementById("tab-help-title");
    const elDesc = document.getElementById("tab-help-desc");
    const elTips = document.getElementById("tab-help-tips");
    const elClose = document.getElementById("tab-help-close");
    const elTour = document.getElementById("tab-help-tour");
    let on = false;
    let scope = "tab"; // "tab" = current tab · "sidebar" = right panel ("?" button)

    // One-time: a "replay walkthrough" link at the top of the panel foot.
    const foot = document.getElementById("tab-help-foot");
    let elWalk = document.getElementById("tab-help-walk");
    if (foot && !elWalk) {
      elWalk = document.createElement("a");
      elWalk.href = "#";
      elWalk.id = "tab-help-walk";
      elWalk.innerHTML = '<i class="fas fa-play"></i> Walk me through this view';
      foot.insertBefore(elWalk, foot.firstChild);
    }

    function currentTab() {
      const active = document.querySelector(".tab-btn.active");
      return active ? active.dataset.tab : typeof activeTab !== "undefined" ? activeTab : "summary";
    }
    function stepFor(tab) {
      return STEPS.find((s) => s.tab === tab) || null;
    }
    // Steps for the stepped Spotlight walkthrough (current scope).
    function walkStepsFor() {
      if (scope === "sidebar") return SIDEBAR_STEPS;
      const tab = currentTab();
      if (TAB_AREAS[tab]) return TAB_AREAS[tab];
      const s = stepFor(tab);
      return s ? [{ sel: "#pane-" + tab, icon: s.icon, title: s.title, desc: s.desc, tips: s.tips }] : [];
    }
    function renderPanel() {
      if (scope === "sidebar") {
        if (elEyebrow) elEyebrow.innerHTML = '<i class="fas fa-circle-info"></i> Right panel';
        elGroup.textContent = "Filters & sample colours";
        elFig.innerHTML = FIG.sidebar || "";
        elTitle.textContent = "The right panel";
        elDesc.textContent =
          "Filters, search, the confidence cutoff and sample colours all live here and apply across every tab at once.";
        elTips.innerHTML = [
          "Search by sample or organism (regex supported).",
          "A sample's colour swatch sets its colour everywhere.",
          "Export Report PDF captures the current filtered state.",
        ]
          .map((t) => `<li>${t}</li>`)
          .join("");
        return;
      }
      if (elEyebrow) elEyebrow.innerHTML = '<i class="fas fa-circle-info"></i> Help mode';
      const tab = currentTab();
      const step = stepFor(tab);
      const info = (typeof PDF_SECTION_INFO !== "undefined" && PDF_SECTION_INFO[tab]) || {};
      elGroup.textContent = info.group || "";
      if (!step) {
        elFig.innerHTML = FIG[tab] || "";
        elTitle.textContent = ((tabBtn(tab) && tabBtn(tab).textContent) || tab).trim();
        elDesc.textContent = info.what || "Help for this view isn't available yet.";
        elTips.innerHTML = info.how ? `<li>${info.how}</li>` : "";
        return;
      }
      elFig.innerHTML = step.fig || FIG[tab] || "";
      elTitle.textContent = step.title || "";
      elDesc.textContent = step.desc || "";
      elTips.innerHTML = (step.tips || []).map((t) => `<li>${t}</li>`).join("");
    }
    function startWalk() {
      const steps = walkStepsFor();
      if (steps && steps.length) Spotlight.start(steps, {});
    }
    function syncButtons() {
      btn.classList.toggle("on", on && scope === "tab");
      btn.setAttribute("aria-pressed", on && scope === "tab" ? "true" : "false");
      if (sideBtn) {
        sideBtn.classList.toggle("on", on && scope === "sidebar");
        sideBtn.setAttribute("aria-pressed", on && scope === "sidebar" ? "true" : "false");
      }
    }
    function setOn(v, nextScope) {
      on = v;
      scope = v && nextScope ? nextScope : v ? scope : "tab";
      syncButtons();
      panel.classList.toggle("open", on);
      panel.setAttribute("aria-hidden", on ? "false" : "true");
      if (on) {
        renderPanel();
        startWalk();
      } else {
        Spotlight.stop();
      }
    }
    btn.addEventListener("click", () => setOn(!(on && scope === "tab"), "tab"));
    if (sideBtn) sideBtn.addEventListener("click", () => setOn(!(on && scope === "sidebar"), "sidebar"));
    elClose.addEventListener("click", () => setOn(false));
    if (elWalk)
      elWalk.addEventListener("click", (e) => {
        e.preventDefault();
        startWalk();
      });
    elTour.addEventListener("click", (e) => {
      e.preventDefault();
      setOn(false);
      open();
    });
    // While help mode is on, follow tab switches: snap back to tab scope,
    // re-render the panel and restart the walkthrough for the new tab.
    const tabbar = document.getElementById("tabbar");
    if (tabbar)
      tabbar.addEventListener("click", () => {
        if (!on) return;
        scope = "tab";
        syncButtons();
        setTimeout(() => {
          renderPanel();
          startWalk();
        }, 90);
      });
  })();
})();
