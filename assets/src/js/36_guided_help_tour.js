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
      tab: "map",
      icon: "fa-map-location-dot",
      title: "Map",
      fig: FIG.map,
      desc: "Plots samples geographically when location metadata is supplied — useful for surveillance and tracking spread across sites.",
      tips: ["Appears once samples carry latitude/longitude metadata.", "Click a pin to focus that sample."],
    },
    {
      tab: "runmeta",
      icon: "fa-tags",
      title: "Run Metadata",
      fig: FIG.runmeta,
      desc: "Per-sample run details — type, platform and any uploaded metadata fields — collected in one table for provenance and QC.",
      tips: [
        "The View button on a sample jumps straight here and highlights its row.",
        "Upload a metadata CSV to enrich these fields.",
      ],
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

  /* ── Per-tab Help mode ───────────────────────────────────────────────
           A toggle that pins a small contextual panel showing help for whatever
           tab the user is on. Reuses the guided-tour STEP content + the PDF
           section groups, and follows tab switches live. */
  (function initHelpMode() {
    const btn = document.getElementById("helpmode-btn");
    const panel = document.getElementById("tab-help-panel");
    if (!btn || !panel) return;
    const elFig = document.getElementById("tab-help-fig");
    const elGroup = document.getElementById("tab-help-group");
    const elTitle = document.getElementById("tab-help-title");
    const elDesc = document.getElementById("tab-help-desc");
    const elTips = document.getElementById("tab-help-tips");
    const elClose = document.getElementById("tab-help-close");
    const elTour = document.getElementById("tab-help-tour");
    let on = false;

    function currentTab() {
      const active = document.querySelector(".tab-btn.active");
      return active ? active.dataset.tab : typeof activeTab !== "undefined" ? activeTab : "summary";
    }
    function stepFor(tab) {
      return STEPS.find((s) => s.tab === tab) || null;
    }
    function renderPanel() {
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
    function setOn(v) {
      on = v;
      btn.classList.toggle("on", on);
      btn.setAttribute("aria-pressed", on ? "true" : "false");
      panel.classList.toggle("open", on);
      panel.setAttribute("aria-hidden", on ? "false" : "true");
      if (on) renderPanel();
    }
    btn.addEventListener("click", () => setOn(!on));
    elClose.addEventListener("click", () => setOn(false));
    elTour.addEventListener("click", (e) => {
      e.preventDefault();
      setOn(false);
      open();
    });
    // Re-render when the user switches tabs while help mode is on.
    const tabbar = document.getElementById("tabbar");
    if (tabbar) tabbar.addEventListener("click", () => on && setTimeout(renderPanel, 60));
  })();
})();
