/* ═══════════════════════════════════════════════════════════════════════════
       -  §  ADMIN DIALOG + REPORT CONFIG
       -     BOOT.report_config  — normalised --report_config (bin/report_admin.py):
       -                           hidden / disabled / ordered / renamed tabs,
       -                           the opening tab, Admin dialog options.
       -     BOOT.run_info       — Nextflow run metadata, every param grouped by
       -                           schema section, software versions. null when
       -                           the dialog is switched off (--report_admin false
       -                           or admin.enabled: false).
       -
       -     Tab state uses its OWN classes (tab-admin-hidden / tab-admin-disabled)
       -     so the data-driven .hidden / .tab-disabled toggles elsewhere can never
       -     re-show a tab the config removed. 15_tab_switching.js refuses clicks
       -     on either class, which covers every programmatic btn.click() too.
       -
       -     Entry points: the "Admin" button in the right panel, a click on the
       -     banner's "Run info" button (opens the Run info section), and
       -     TTAdmin.open(section).
═══════════════════════════════════════════════════════════════════════════ */
const TTAdmin = (function () {
  const CFG = (typeof BOOT !== "undefined" && BOOT.report_config) || {};
  const RUN = (typeof BOOT !== "undefined" && BOOT.run_info) || null;
  const ALL_SECTIONS = ["tabs", "params", "run", "versions", "config"];
  const cfgTabs = Object.assign({ hidden: [], disabled: [], default: null, order: [], labels: {} }, CFG.tabs || {});
  const cfgAdmin = Object.assign({ enabled: true, allow_tab_toggle: true, sections: {}, redact: [] }, CFG.admin || {});
  ALL_SECTIONS.forEach((s) => {
    if (cfgAdmin.sections[s] === undefined) cfgAdmin.sections[s] = true;
  });
  // Older / hand-built reports carry no run_info at all; the dialog still
  // offers tab controls + config, and says why the rest is empty.
  const ENABLED = cfgAdmin.enabled !== false && !(BOOT && BOOT.report_config && BOOT.run_info === null);

  const STORE_KEY = "taxtriage:adminTabs:" + ((typeof BOOT !== "undefined" && BOOT.report_generated_at) || "report");

  const SECTION_META = {
    tabs: { icon: "fa-table-columns", label: "Tabs", sub: "Show, hide, disable" },
    params: { icon: "fa-sliders", label: "Parameters", sub: "Every pipeline param" },
    run: { icon: "fa-circle-info", label: "Run info", sub: "Nextflow + report build" },
    versions: { icon: "fa-cubes", label: "Software", sub: "Tool versions" },
    config: { icon: "fa-file-code", label: "Report config", sub: "--report_config file" },
  };

  const esc = (s) =>
    String(s == null ? "" : s)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");

  // ── Tab state ────────────────────────────────────────────────────────────
  function fromConfig() {
    return {
      hidden: (cfgTabs.hidden || []).slice(),
      disabled: (cfgTabs.disabled || []).slice(),
      default: cfgTabs.default || null,
    };
  }
  function loadSaved() {
    if (!cfgAdmin.allow_tab_toggle) return null;
    try {
      const raw = localStorage.getItem(STORE_KEY);
      if (!raw) return null;
      const s = JSON.parse(raw);
      if (!s || !Array.isArray(s.hidden) || !Array.isArray(s.disabled)) return null;
      return { hidden: s.hidden, disabled: s.disabled, default: s.default || null };
    } catch (e) {
      return null;
    }
  }
  function save() {
    try {
      localStorage.setItem(STORE_KEY, JSON.stringify(state));
    } catch (e) {}
  }
  function clearSaved() {
    try {
      localStorage.removeItem(STORE_KEY);
    } catch (e) {}
  }
  let state = loadSaved() || fromConfig();
  const viewerChanged = () => JSON.stringify(state) !== JSON.stringify(fromConfig());

  const tabButtons = () => Array.from(document.querySelectorAll("#tabbar .tab-btn"));
  const tabBtn = (id) => document.querySelector(`#tabbar .tab-btn[data-tab="${id}"]`);

  // Label text = the button's first non-blank text node (icons are <i> siblings).
  function labelNode(btn) {
    for (const n of btn.childNodes) if (n.nodeType === 3 && n.textContent.trim()) return n;
    return null;
  }
  function originalLabel(btn) {
    if (btn.dataset.ttLabel === undefined) {
      const n = labelNode(btn);
      btn.dataset.ttLabel = n ? n.textContent.trim() : btn.dataset.tab;
    }
    return btn.dataset.ttLabel;
  }
  function tabLabel(btn) {
    return (cfgTabs.labels || {})[btn.dataset.tab] || originalLabel(btn);
  }
  function tabIcon(btn) {
    const i = btn.querySelector("i");
    return i ? i.className : "fas fa-square";
  }

  function applyLabels() {
    tabButtons().forEach((btn) => {
      const want = (cfgTabs.labels || {})[btn.dataset.tab];
      originalLabel(btn);
      if (!want) return;
      const n = labelNode(btn);
      if (n) n.textContent = " " + want + " ";
    });
  }
  function applyOrder() {
    const bar = document.getElementById("tabbar");
    if (!bar || !(cfgTabs.order || []).length) return;
    const btns = tabButtons();
    const ordered = [];
    cfgTabs.order.forEach((id) => {
      const b = btns.find((x) => x.dataset.tab === id);
      if (b && !ordered.includes(b)) ordered.push(b);
    });
    btns.forEach((b) => ordered.includes(b) || ordered.push(b));
    ordered.forEach((b) => bar.appendChild(b));
  }
  function applyTitle() {
    if (!CFG.title) return;
    const banner = document.getElementById("banner");
    if (banner) {
      for (const n of banner.childNodes) {
        if (n.nodeType === 3 && n.textContent.trim()) {
          n.textContent = " " + CFG.title + " ";
          break;
        }
      }
    }
    const h1 = document.querySelector("#pdf-report-cover h1");
    if (h1) h1.textContent = CFG.title;
    document.title = CFG.title;
  }
  function applyStates() {
    const hid = new Set(state.hidden),
      dis = new Set(state.disabled);
    tabButtons().forEach((btn) => {
      const id = btn.dataset.tab;
      btn.classList.toggle("tab-admin-hidden", hid.has(id));
      btn.classList.toggle("tab-admin-disabled", !hid.has(id) && dis.has(id));
      if (dis.has(id) && !hid.has(id)) {
        if (btn.dataset.ttTitle === undefined) btn.dataset.ttTitle = btn.getAttribute("title") || "";
        btn.setAttribute("title", "Disabled by the report configuration (Admin → Tabs)");
        btn.setAttribute("aria-disabled", "true");
      } else if (btn.dataset.ttTitle !== undefined) {
        btn.setAttribute("title", btn.dataset.ttTitle);
        btn.removeAttribute("aria-disabled");
        delete btn.dataset.ttTitle;
      }
    });
  }
  const usable = (btn) =>
    btn &&
    !btn.classList.contains("hidden") &&
    !btn.classList.contains("tab-disabled") &&
    !btn.classList.contains("tab-admin-hidden") &&
    !btn.classList.contains("tab-admin-disabled");

  // If the active tab is now hidden/disabled, move to the default (or first
  // usable) tab. `preferDefault` also jumps to the default when it is usable.
  function ensureActiveAllowed(preferDefault) {
    const active = document.querySelector("#tabbar .tab-btn.active");
    const dflt = state.default && tabBtn(state.default);
    if (preferDefault && usable(dflt) && dflt !== active) return dflt.click();
    if (usable(active)) return;
    const target = usable(dflt) ? dflt : tabButtons().find(usable);
    if (target) target.click();
  }
  function applyAll() {
    applyStates();
  }

  // Early pass (parse time): hide/disable + labels + order before first paint
  // so a removed tab never flashes. No clicks yet — the report isn't built.
  try {
    applyTitle();
    applyLabels();
    applyOrder();
    applyStates();
  } catch (e) {
    console.warn("[taxtriage] report config:", e);
  }

  // ── Dialog ───────────────────────────────────────────────────────────────
  let overlay = null,
    current = null;
  const visibleSections = () => ALL_SECTIONS.filter((s) => cfgAdmin.sections[s] !== false);

  function ensureDialog() {
    if (overlay) return overlay;
    overlay = document.createElement("div");
    overlay.id = "tt-admin-overlay";
    overlay.className = "export-modal-overlay";
    const nav = visibleSections()
      .map((s) => {
        const m = SECTION_META[s];
        return (
          `<button type="button" class="tt-admin-navbtn" data-sec="${s}" role="tab" aria-selected="false">` +
          `<i class="fas ${m.icon}"></i><span class="tt-admin-nav-main">${m.label}</span>` +
          `<span class="tt-admin-nav-sub">${m.sub}</span></button>`
        );
      })
      .join("");
    overlay.innerHTML = `
      <div class="export-modal tt-admin-modal" role="dialog" aria-modal="true" aria-labelledby="tt-admin-title">
        <header>
          <i class="fas fa-user-gear"></i>
          <span id="tt-admin-title">Admin</span>
          <span class="tt-admin-head-sub">${esc((RUN && RUN.workflow && RUN.workflow.runName) || "")}</span>
          <button type="button" id="tt-admin-close" title="Close" aria-label="Close">&times;</button>
        </header>
        <nav class="tt-admin-nav" role="tablist" aria-label="Admin sections">${nav}</nav>
        <div class="tt-admin-body" id="tt-admin-body"></div>
      </div>`;
    document.body.appendChild(overlay);
    const close = () => (overlay.style.display = "none");
    overlay.querySelector("#tt-admin-close").addEventListener("click", close);
    overlay.addEventListener("mousedown", (e) => {
      if (e.target === overlay) close();
    });
    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape" && overlay.style.display === "flex") close();
    });
    overlay.querySelector(".tt-admin-nav").addEventListener("click", (e) => {
      const b = e.target.closest(".tt-admin-navbtn");
      if (b) show(b.dataset.sec);
    });
    overlay.querySelector("#tt-admin-body").addEventListener("click", onBodyClick);
    overlay.querySelector("#tt-admin-body").addEventListener("change", onBodyChange);
    overlay.querySelector("#tt-admin-body").addEventListener("input", onBodyInput);
    return overlay;
  }

  function show(sec) {
    const secs = visibleSections();
    if (!secs.includes(sec)) sec = secs[0];
    current = sec;
    overlay.querySelectorAll(".tt-admin-navbtn").forEach((b) => {
      const on = b.dataset.sec === sec;
      b.classList.toggle("active", on);
      b.setAttribute("aria-selected", on ? "true" : "false");
    });
    const body = overlay.querySelector("#tt-admin-body");
    body.innerHTML = RENDER[sec] ? RENDER[sec]() : "";
    body.scrollTop = 0;
    if (sec === "params") filterParams();
  }

  function open(sec) {
    if (!ENABLED) return;
    ensureDialog();
    overlay.style.display = "flex";
    show(sec || current || visibleSections()[0]);
  }

  const notCaptured = (what) =>
    `<div class="tt-admin-empty"><i class="fas fa-circle-exclamation"></i> ${what} ` +
    `was not captured for this report. It is recorded when the report is built by the pipeline ` +
    `(<code>--report_admin true</code>, the default).</div>`;

  function fmtVal(v) {
    if (v === null || v === undefined || v === "") return `<span class="tt-admin-na">—</span>`;
    if (v === true || v === false) return `<span class="tt-admin-bool tt-admin-bool-${v}">${v}</span>`;
    if (typeof v === "object") return `<code>${esc(JSON.stringify(v))}</code>`;
    return `<code>${esc(v)}</code>`;
  }

  // ── Section: Tabs ─────────────────────────────────────────────────────────
  function renderTabs() {
    const locked = !cfgAdmin.allow_tab_toggle;
    const rows = tabButtons()
      .map((btn) => {
        const id = btn.dataset.tab;
        const hidden = state.hidden.includes(id);
        const disabled = state.disabled.includes(id);
        const noData = btn.classList.contains("hidden");
        const dis = locked ? "disabled" : "";
        return `<tr data-tab="${id}" class="${hidden ? "is-hidden" : ""}">
          <td><i class="${esc(tabIcon(btn))}"></i> ${esc(tabLabel(btn))}
            <span class="tt-admin-tabid">${id}</span></td>
          <td>${
            noData
              ? '<span class="tt-admin-pill muted" title="This run has no data for this tab, so it stays hidden whatever you set here">no data</span>'
              : '<span class="tt-admin-pill ok">available</span>'
          }</td>
          <td class="c"><input type="checkbox" data-act="visible" ${
            hidden ? "" : "checked"
          } ${dis} aria-label="Show ${esc(tabLabel(btn))}"></td>
          <td class="c"><input type="checkbox" data-act="enabled" ${disabled ? "" : "checked"} ${
            dis || (hidden ? "disabled" : "")
          } aria-label="Enable ${esc(tabLabel(btn))}"></td>
          <td class="c"><input type="radio" name="tt-admin-default" data-act="default" ${
            state.default === id ? "checked" : ""
          } ${dis || (hidden || disabled ? "disabled" : "")} aria-label="Open on ${esc(tabLabel(btn))}"></td>
        </tr>`;
      })
      .join("");
    const note = locked
      ? `<div class="tt-admin-note warn"><i class="fas fa-lock"></i> Tab visibility is locked by this report's configuration (<code>admin.allow_tab_toggle: false</code>).</div>`
      : viewerChanged()
      ? `<div class="tt-admin-note"><i class="fas fa-user-pen"></i> You have changed the tab layout from the report's defaults. Your changes are kept in this browser only.</div>`
      : "";
    return `
      <p class="tt-admin-intro">Choose which tabs this report shows. <b>Visible</b> removes a tab from the bar,
        <b>Enabled</b> greys it out, and <b>Opens on</b> sets the tab the report starts on. To make a layout the
        default for future runs, download it as a config and pass it to the pipeline with
        <code>--report_config</code>.</p>
      ${note}
      <table class="tt-admin-table tt-admin-tabs">
        <thead><tr><th>Tab</th><th>Data</th><th class="c">Visible</th><th class="c">Enabled</th><th class="c">Opens on</th></tr></thead>
        <tbody>${rows}</tbody>
      </table>
      <div class="tt-admin-actions">
        <button type="button" data-act="tabs-reset" ${
          locked ? "disabled" : ""
        }><i class="fas fa-rotate-left"></i> Reset to report defaults</button>
        <span class="tt-admin-spacer"></span>
        <button type="button" data-act="dl-config-yml"><i class="fas fa-download"></i> Config (YAML)</button>
        <button type="button" data-act="dl-config-json"><i class="fas fa-download"></i> Config (JSON)</button>
      </div>`;
  }

  // ── Section: Parameters ───────────────────────────────────────────────────
  function renderParams() {
    if (!RUN || !RUN.param_groups || !RUN.param_groups.length) return notCaptured("The pipeline parameter list");
    const nChanged = RUN.param_groups.reduce((n, g) => n + g.params.filter((p) => p.changed).length, 0);
    const nAll = RUN.param_groups.reduce((n, g) => n + g.params.length, 0);
    const groups = RUN.param_groups
      .map((g, gi) => {
        const rows = g.params
          .map(
            (p) => `<tr class="tt-admin-prow${p.changed ? " changed" : ""}${p.hidden ? " schema-hidden" : ""}"
                data-name="${esc(p.name)}" data-search="${esc(
                  (p.name + " " + JSON.stringify(p.value) + " " + (p.description || "")).toLowerCase(),
                )}">
              <td class="pname"><code>${esc(p.name)}</code>${
                p.changed
                  ? ' <span class="tt-admin-pill chg" title="Differs from the schema default">changed</span>'
                  : ""
              }
                ${p.description ? `<div class="pdesc">${esc(p.description)}</div>` : ""}</td>
              <td>${fmtVal(p.value)}</td>
              <td class="pdef">${fmtVal(p.default)}</td>
            </tr>`,
          )
          .join("");
        return `<details class="tt-admin-group" open data-gi="${gi}">
            <summary><span>${esc(g.title)}</span> <span class="tt-admin-count"></span></summary>
            ${g.description ? `<div class="tt-admin-gdesc">${esc(g.description)}</div>` : ""}
            <table class="tt-admin-table"><thead><tr><th>Parameter</th><th>Value</th><th>Default</th></tr></thead>
            <tbody>${rows}</tbody></table>
          </details>`;
      })
      .join("");
    return `
      <p class="tt-admin-intro">Every pipeline parameter for this run, grouped as in <code>nextflow_schema.json</code>
        — the complete version of MultiQC's "Workflow Summary". <b>${nChanged}</b> of ${nAll} differ from their defaults.
        ${RUN.captured_at ? `Captured ${esc(fmtDate(RUN.captured_at))}.` : ""}</p>
      <div class="tt-admin-toolbar">
        <input type="search" id="tt-admin-psearch" placeholder="Filter parameters…" autocomplete="off" spellcheck="false">
        <label><input type="checkbox" id="tt-admin-pchanged"> Changed only</label>
        <label title="Params the schema marks as hidden (generic nf-core / institutional options)"><input type="checkbox" id="tt-admin-phidden"> Include hidden</label>
        <span class="tt-admin-spacer"></span>
        <button type="button" data-act="dl-params-file" title="The changed params as a Nextflow -params-file, to reproduce this run"><i class="fas fa-download"></i> Params file</button>
        <button type="button" data-act="dl-params-all" title="Every param with its default, description and group"><i class="fas fa-download"></i> All (JSON)</button>
      </div>
      <div id="tt-admin-pgroups">${groups}</div>
      <div class="tt-admin-empty" id="tt-admin-pnone" hidden>No parameters match.</div>`;
  }
  function filterParams() {
    const body = overlay && overlay.querySelector("#tt-admin-body");
    if (!body || current !== "params") return;
    const q = ((body.querySelector("#tt-admin-psearch") || {}).value || "").trim().toLowerCase();
    const onlyChanged = !!(body.querySelector("#tt-admin-pchanged") || {}).checked;
    const withHidden = !!(body.querySelector("#tt-admin-phidden") || {}).checked;
    let any = 0;
    body.querySelectorAll(".tt-admin-group").forEach((g) => {
      let n = 0;
      g.querySelectorAll(".tt-admin-prow").forEach((r) => {
        const ok =
          (!q || r.dataset.search.includes(q)) &&
          (!onlyChanged || r.classList.contains("changed")) &&
          (withHidden || !r.classList.contains("schema-hidden") || r.classList.contains("changed"));
        r.hidden = !ok;
        if (ok) n++;
      });
      g.hidden = n === 0;
      const c = g.querySelector(".tt-admin-count");
      if (c) c.textContent = `(${n})`;
      any += n;
    });
    const none = body.querySelector("#tt-admin-pnone");
    if (none) none.hidden = any > 0;
  }

  // ── Section: Run info ─────────────────────────────────────────────────────
  const WF_LABELS = [
    ["runName", "Run name"],
    ["sessionId", "Session ID"],
    ["start", "Started"],
    ["resume", "Resumed"],
    ["pipeline", "Pipeline"],
    ["pipelineVersion", "Pipeline version"],
    ["repository", "Repository"],
    ["revision", "Revision"],
    ["commitId", "Commit"],
    ["nextflowVersion", "Nextflow version"],
    ["nextflowBuild", "Nextflow build"],
    ["profile", "Profile(s)"],
    ["containerEngine", "Container engine"],
    ["container", "Container(s)"],
    ["userName", "User"],
    ["launchDir", "Launch dir"],
    ["workDir", "Work dir"],
    ["projectDir", "Project dir"],
    ["outdir", "Output dir"],
    ["configFiles", "Config files"],
  ];
  function fmtDate(iso) {
    try {
      const d = new Date(iso);
      if (isNaN(d)) return iso;
      return d.toLocaleString(undefined, {
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
  }
  const kv = (label, html) => `<tr><th>${esc(label)}</th><td>${html}</td></tr>`;

  function renderRun() {
    const wf = (RUN && RUN.workflow) || {};
    const B = typeof BOOT !== "undefined" ? BOOT : {};
    let wfHtml;
    if (!RUN || !RUN.captured) {
      wfHtml = notCaptured("Nextflow run metadata");
    } else {
      const rows = WF_LABELS.filter(([k]) => wf[k] !== undefined && wf[k] !== null && wf[k] !== "")
        .map(([k, label]) => {
          let v = wf[k];
          if (k === "start") return kv(label, `<code>${esc(fmtDate(v))}</code>`);
          if (k === "commitId" && /^[0-9a-f]{7,40}$/i.test(v))
            return kv(
              label,
              `<a href="https://github.com/jhuapl-bio/taxtriage/commit/${esc(
                v,
              )}" target="_blank" rel="noopener"><code>${esc(v)}</code></a>`,
            );
          if (Array.isArray(v))
            return kv(
              label,
              v.length ? v.map((x) => `<code>${esc(x)}</code>`).join("<br>") : '<span class="tt-admin-na">—</span>',
            );
          return kv(label, fmtVal(v));
        })
        .join("");
      wfHtml = `<table class="tt-admin-kv">${rows}</table>`;
      if (wf.commandLine) {
        wfHtml += `<div class="tt-admin-cmd"><div class="tt-admin-cmd-head">Command line
          <button type="button" data-act="copy" data-copy="cmd"><i class="fas fa-copy"></i> Copy</button></div>
          <pre id="tt-admin-cmdline">${esc(wf.commandLine)}</pre></div>`;
      }
    }

    // Report contents — everything here comes from BOOT, so it works even
    // when run info was not captured.
    const samples = Object.keys(B.sample_meta || {}).length;
    const recs = (B.records || []).length;
    const feat = [
      ["VF / AMR annotations", !!B.has_prot],
      ["Novelty detection", !!B.has_novelty],
      ["Pathogen sheet", !!B.has_pathogens],
      ["HMP healthy-range reference", !!B.hmp],
      ["In-silico suite", !!B.has_insilico_suite],
      ["Geographic coordinates", !!B.has_geo],
      ["Sample QC default rules", !!(B.sample_flags && (B.sample_flags.rules || []).length)],
      ["Organism QC default rules", !!(B.organism_flags && (B.organism_flags.rules || []).length)],
    ]
      .map(
        ([l, on]) =>
          `<span class="tt-admin-pill ${on ? "ok" : "muted"}"><i class="fas ${
            on ? "fa-check" : "fa-minus"
          }"></i> ${l}</span>`,
      )
      .join(" ");
    const tassSrc = (() => {
      const m = Object.values(B.sample_meta || {}).find((x) => x && x.best_cutoffs_source);
      return m ? m.best_cutoffs_source : null;
    })();
    const bc = ((B.best_cutoffs || {}).subkey || {}).best_threshold;
    const build = (RUN && RUN.report_build) || {};
    const rep = [
      kv("Report built", B.report_generated_at ? `<code>${esc(fmtDate(B.report_generated_at))}</code>` : fmtVal(null)),
      kv("Branch / revision", fmtVal(B.pipeline_revision || "Not specified or local build")),
      kv("Commit", fmtVal(B.pipeline_commit)),
      kv("Samples", fmtVal(samples)),
      kv("Detection rows", fmtVal(recs)),
      kv(
        "Default TASS cutoff",
        bc != null
          ? `<code>${esc(bc)}</code>${tassSrc ? ` <span class="tt-admin-na">(${esc(tassSrc)})</span>` : ""}`
          : fmtVal(null),
      ),
      build.python ? kv("Build Python", fmtVal(build.python)) : "",
      build.argv ? kv("make_report.py", `<pre class="tt-admin-pre">${esc(build.argv)}</pre>`) : "",
    ].join("");
    const red =
      RUN && RUN.redacted && RUN.redacted.length
        ? `<div class="tt-admin-note"><i class="fas fa-user-secret"></i> Redacted by the report config: ${RUN.redacted
            .map((x) => `<code>${esc(x)}</code>`)
            .join(", ")}</div>`
        : "";
    return `
      <h4 class="tt-admin-h">Nextflow run</h4>${red}${wfHtml}
      <h4 class="tt-admin-h">This report</h4>
      <table class="tt-admin-kv">${rep}</table>
      <div class="tt-admin-feats">${feat}</div>
      <div class="tt-admin-actions"><span class="tt-admin-spacer"></span>
        <button type="button" data-act="dl-runinfo"><i class="fas fa-download"></i> Run info (JSON)</button></div>`;
  }

  // ── Section: Software versions ────────────────────────────────────────────
  function renderVersions() {
    const v = (RUN && RUN.software_versions) || [];
    if (!v.length) return notCaptured("The software version list");
    const tools = {};
    v.forEach((r) => (tools[r.tool] = tools[r.tool] || new Set()).add(r.version));
    const multi = Object.entries(tools).filter(([, s]) => s.size > 1);
    const rows = v
      .map(
        (r) => `<tr class="tt-admin-vrow" data-search="${esc(
          (r.process + " " + r.tool + " " + r.version).toLowerCase(),
        )}">
          <td><code title="${esc(r.process_full || r.process)}">${esc(r.process)}</code></td>
          <td>${esc(r.tool)}</td><td><code>${esc(r.version)}</code></td></tr>`,
      )
      .join("");
    return `
      <p class="tt-admin-intro">Tool versions reported by the processes that fed this report (nf-core <code>versions.yml</code>).
        ${Object.keys(tools).length} tools across ${new Set(v.map((r) => r.process)).size} processes.
        ${
          multi.length
            ? `<br><i class="fas fa-triangle-exclamation" style="color:#f59f00"></i> Several versions of: ${multi
                .map(([t]) => `<code>${esc(t)}</code>`)
                .join(", ")}.`
            : ""
        }</p>
      <div class="tt-admin-toolbar">
        <input type="search" id="tt-admin-vsearch" placeholder="Filter tools / processes…" autocomplete="off" spellcheck="false">
        <span class="tt-admin-spacer"></span>
        <button type="button" data-act="dl-versions"><i class="fas fa-download"></i> CSV</button>
      </div>
      <table class="tt-admin-table"><thead><tr><th>Process</th><th>Tool</th><th>Version</th></tr></thead><tbody>${rows}</tbody></table>`;
  }

  // ── Section: Report config ────────────────────────────────────────────────
  function currentConfig() {
    const out = {};
    if (CFG.title) out.title = CFG.title;
    const tabs = {};
    if (state.hidden.length) tabs.hidden = state.hidden.slice();
    if (state.disabled.length) tabs.disabled = state.disabled.slice();
    if (state.default) tabs.default = state.default;
    if ((cfgTabs.order || []).length) tabs.order = cfgTabs.order.slice();
    if (Object.keys(cfgTabs.labels || {}).length) tabs.labels = Object.assign({}, cfgTabs.labels);
    out.tabs = tabs;
    const admin = { enabled: cfgAdmin.enabled !== false, allow_tab_toggle: !!cfgAdmin.allow_tab_toggle };
    const secs = ALL_SECTIONS.filter((s) => cfgAdmin.sections[s] !== false);
    if (secs.length !== ALL_SECTIONS.length) admin.sections = secs;
    if ((cfgAdmin.redact || []).length) admin.redact = cfgAdmin.redact.slice();
    out.admin = admin;
    return out;
  }
  function toYaml(obj, ind) {
    ind = ind || "";
    const scalar = (v) => {
      if (v === null || v === undefined) return "null";
      if (typeof v === "boolean" || typeof v === "number") return String(v);
      const s = String(v);
      return /^[A-Za-z0-9_][A-Za-z0-9_ .\/-]*$/.test(s) &&
        !/^(true|false|null|yes|no|on|off)$/i.test(s) &&
        !/^[\d.]+$/.test(s)
        ? s
        : JSON.stringify(s);
    };
    return Object.entries(obj)
      .map(([k, v]) => {
        if (Array.isArray(v)) return `${ind}${k}: [${v.map(scalar).join(", ")}]`;
        if (v && typeof v === "object") {
          return Object.keys(v).length ? `${ind}${k}:\n${toYaml(v, ind + "  ")}` : `${ind}${k}: {}`;
        }
        return `${ind}${k}: ${scalar(v)}`;
      })
      .join("\n");
  }
  function configYaml() {
    return (
      "# TaxTriage interactive report config — pass with:  --report_config report_config.yml\n" +
      "# Tab ids: " +
      tabButtons()
        .map((b) => b.dataset.tab)
        .join(", ") +
      "\n" +
      toYaml(currentConfig()) +
      "\n"
    );
  }
  function renderConfig() {
    const src = CFG.source
      ? `This report was built with <code>${esc(CFG.source)}</code>.`
      : "This report was built without a <code>--report_config</code> file (all defaults).";
    return `
      <p class="tt-admin-intro">${src} Below is the configuration as it currently applies${
        viewerChanged() ? ", <b>including your tab changes</b>" : ""
      }.
        Save it and pass it to the pipeline with <code>--report_config report_config.yml</code>
        (JSON works too), or check one locally with <code>bin/report_admin.py --validate FILE</code>.</p>
      <pre class="tt-admin-pre tt-admin-yaml" id="tt-admin-yaml">${esc(configYaml())}</pre>
      <div class="tt-admin-actions"><span class="tt-admin-spacer"></span>
        <button type="button" data-act="copy" data-copy="yaml"><i class="fas fa-copy"></i> Copy</button>
        <button type="button" data-act="dl-config-yml"><i class="fas fa-download"></i> YAML</button>
        <button type="button" data-act="dl-config-json"><i class="fas fa-download"></i> JSON</button>
      </div>
      <details class="tt-admin-group"><summary>Reference</summary>
<pre class="tt-admin-pre">title: "My run"                 # banner + page title
tabs:
  hidden:   [novelty, insilico]   # removed from the tab bar
  disabled: [explore]             # visible but greyed out
  show:     [summary, heatmap]    # whitelist: hide every other tab
  default:  heatmap               # tab the report opens on
  order:    [summary, table]      # these first, the rest keep their order
  labels:   {table: Detections}   # rename tabs
admin:
  enabled: true                   # false: no Admin button, no run info in the HTML
  allow_tab_toggle: true          # false: viewers cannot change tab visibility
  sections: [tabs, params, run, versions, config]
  redact: [paths, userName]       # "paths" = dirs, config files, command line, user</pre></details>`;
  }

  const RENDER = {
    tabs: renderTabs,
    params: renderParams,
    run: renderRun,
    versions: renderVersions,
    config: renderConfig,
  };

  // ── Downloads / copy ──────────────────────────────────────────────────────
  function download(name, text, type) {
    const blob = new Blob([text], { type: type || "text/plain" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1500);
  }
  function copy(text, btn) {
    const done = () => {
      if (!btn) return;
      const old = btn.innerHTML;
      btn.innerHTML = '<i class="fas fa-check"></i> Copied';
      setTimeout(() => (btn.innerHTML = old), 1200);
    };
    try {
      navigator.clipboard.writeText(text).then(done, () => fallbackCopy(text, done));
    } catch (e) {
      fallbackCopy(text, done);
    }
  }
  function fallbackCopy(text, done) {
    const ta = document.createElement("textarea");
    ta.value = text;
    document.body.appendChild(ta);
    ta.select();
    try {
      document.execCommand("copy");
      done();
    } catch (e) {}
    ta.remove();
  }
  const csvCell = (v) => (/[",\n]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : String(v));
  const stem = () => {
    const rn = RUN && RUN.workflow && RUN.workflow.runName;
    return "taxtriage" + (rn && rn !== "(redacted)" ? "_" + String(rn).replace(/[^\w.-]+/g, "_") : "");
  };

  function onBodyClick(e) {
    const b = e.target.closest("button[data-act]");
    if (!b) return;
    const act = b.dataset.act;
    if (act === "tabs-reset") {
      state = fromConfig();
      clearSaved();
      applyAll();
      ensureActiveAllowed(false);
      show("tabs");
    } else if (act === "dl-config-yml") {
      download("report_config.yml", configYaml(), "text/yaml");
    } else if (act === "dl-config-json") {
      download("report_config.json", JSON.stringify(currentConfig(), null, 2) + "\n", "application/json");
    } else if (act === "dl-params-file") {
      const p = {};
      (RUN.param_groups || []).forEach((g) =>
        g.params.forEach((x) => {
          if (x.changed && x.value !== "(redacted)") p[x.name] = x.value;
        }),
      );
      download(stem() + "_params.json", JSON.stringify(p, null, 2) + "\n", "application/json");
    } else if (act === "dl-params-all") {
      download(stem() + "_all_params.json", JSON.stringify(RUN.param_groups, null, 2) + "\n", "application/json");
    } else if (act === "dl-runinfo") {
      download(stem() + "_run_info.json", JSON.stringify(RUN || {}, null, 2) + "\n", "application/json");
    } else if (act === "dl-versions") {
      const rows = [["process", "tool", "version"]].concat(
        (RUN.software_versions || []).map((r) => [r.process_full || r.process, r.tool, r.version]),
      );
      download(
        stem() + "_software_versions.csv",
        rows.map((r) => r.map(csvCell).join(",")).join("\n") + "\n",
        "text/csv",
      );
    } else if (act === "copy") {
      const src = b.dataset.copy === "cmd" ? (RUN.workflow || {}).commandLine : configYaml();
      copy(src || "", b);
    }
  }
  function onBodyChange(e) {
    const el = e.target;
    if (current === "params") return filterParams();
    if (current !== "tabs" || !el.dataset.act) return;
    const id = el.closest("tr").dataset.tab;
    const drop = (arr) => arr.filter((x) => x !== id);
    if (el.dataset.act === "visible") {
      state.hidden = el.checked ? drop(state.hidden) : state.hidden.concat(id);
    } else if (el.dataset.act === "enabled") {
      state.disabled = el.checked ? drop(state.disabled) : state.disabled.concat(id);
    } else if (el.dataset.act === "default") {
      state.default = id;
    }
    if (state.default && (state.hidden.includes(state.default) || state.disabled.includes(state.default)))
      state.default = null;
    // Never let a viewer hide every tab.
    if (!tabButtons().some((b) => !state.hidden.includes(b.dataset.tab) && !state.disabled.includes(b.dataset.tab))) {
      state.hidden = drop(state.hidden);
      state.disabled = drop(state.disabled);
    }
    save();
    applyAll();
    ensureActiveAllowed(false);
    show("tabs");
  }
  function onBodyInput(e) {
    if (e.target.id === "tt-admin-psearch") return filterParams();
    if (e.target.id === "tt-admin-vsearch") {
      const q = e.target.value.trim().toLowerCase();
      overlay.querySelectorAll(".tt-admin-vrow").forEach((r) => (r.hidden = !!q && !r.dataset.search.includes(q)));
    }
  }

  // ── Wiring ────────────────────────────────────────────────────────────────
  function wire() {
    const btn = document.getElementById("admin-btn");
    if (btn) {
      if (!ENABLED) btn.remove();
      else btn.addEventListener("click", () => open());
    }
    // Banner "Run info" keeps its hover tooltip; a click opens the full view.
    const info = document.getElementById("pipeline-info-btn");
    if (info && ENABLED && cfgAdmin.sections.run !== false) {
      info.style.cursor = "pointer";
      info.addEventListener("click", () => {
        if (typeof hideTip === "function") hideTip();
        open("run");
      });
    }
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", wire, { once: true });
  else wire();

  return {
    open,
    enabled: ENABLED,
    // Called by the deferred init scheduler (28_meta_csv.js) once the report
    // is built: re-apply (init code may have toggled tabs) and move to the
    // configured opening tab.
    afterInit() {
      applyAll();
      ensureActiveAllowed(true);
    },
    usable,
    config: () => currentConfig(),
  };
})();
