/* ═══════════════════════════════════════════════════════════════════════════
       -  §  SAMPLE FLAGS — UI
       -     The rule builder dialog (#flag-modal-overlay) and the compact
       -     sidebar summary that opens it. The rule MODEL and evaluation live
       -     in 41_sample_flags.js; this file only reads/writes TT_FLAG_RULES
       -     and calls ttFlagRefresh() after every edit.
       -
       -     Edits apply immediately (the report behind the dialog updates as
       -     you type) — "Done" just closes. That way the flagged-sample list in
       -     the preview and the charts never disagree.
═══════════════════════════════════════════════════════════════════════════ */

/* ── Sidebar summary ─────────────────────────────────────────────────────── */
function ttFlagRenderSummary() {
  const line = document.getElementById("flag-summary-line");
  const toggle = document.getElementById("flag-hide-toggle");
  const btn = document.getElementById("flag-open-btn");
  if (!line) return;

  const c = ttFlagCounts();
  if (toggle) toggle.checked = !!TT_FLAG_HIDE_ALL;
  if (toggle) toggle.disabled = c.flagged === 0;
  if (btn) btn.classList.toggle("has-flags", c.flagged > 0);

  if (!c.rules) {
    line.className = "flag-line quiet";
    line.innerHTML = `No sample rules active. <span class="flag-line-sub">Add one to flag samples by reads, organism count or metadata.</span>`;
  } else if (!c.total) {
    //  Rules configured but no data yet (a freshly emptied report, or the
    //  GitHub Pages "Start empty" path). "All 0 samples pass" would be true
    //  but useless.
    line.className = "flag-line quiet";
    line.innerHTML =
      `${c.rules} sample rule${c.rules === 1 ? "" : "s"} ready. ` +
      `<span class="flag-line-sub">They apply as soon as samples are loaded.</span>`;
  } else if (c.flagged === 0) {
    line.className = "flag-line ok";
    line.innerHTML =
      `<i class="fas fa-circle-check"></i> All ${c.total} sample${c.total === 1 ? "" : "s"} pass ` +
      `${c.rules} rule${c.rules === 1 ? "" : "s"}.`;
  } else {
    line.className = "flag-line warn";
    line.innerHTML =
      `<i class="fas fa-flag"></i> <b>${c.flagged}</b> of ${c.total} sample${c.total === 1 ? "" : "s"} flagged` +
      (c.hidden ? ` · <b>${c.hidden}</b> hidden` : "") +
      `<span class="flag-line-sub">by ${c.rules} rule${c.rules === 1 ? "" : "s"}</span>`;
  }
}

/* ── Dialog: open / close ────────────────────────────────────────────────── */
function ttFlagOpenModal() {
  const ov = document.getElementById("flag-modal-overlay");
  if (!ov) return;
  ov.style.display = "flex";
  ov.setAttribute("aria-hidden", "false");
  _ttFlagSyncTopControls();
  ttFlagRenderRules();
  ttFlagRenderPreview();
  const first = ov.querySelector("select, input, button");
  if (first) first.focus();
}

function ttFlagCloseModal() {
  const ov = document.getElementById("flag-modal-overlay");
  if (!ov) return;
  ov.style.display = "none";
  ov.setAttribute("aria-hidden", "true");
}

function _ttFlagSyncTopControls() {
  const en = document.getElementById("flag-enabled");
  const lg = document.getElementById("flag-logic");
  const mf = document.getElementById("flag-missing-fails");
  const ha = document.getElementById("flag-hide-all");
  if (en) en.checked = !!TT_FLAG_ENABLED;
  if (lg) lg.value = TT_FLAG_LOGIC;
  if (mf) mf.checked = !!TT_FLAG_MISSING_FAILS;
  if (ha) ha.checked = !!TT_FLAG_HIDE_ALL;
}

/* ── Dialog: the rule rows ───────────────────────────────────────────────── */
const _TT_FLAG_SOURCES = [
  { key: "meta", label: "Sample metrics" },
  { key: "derived", label: "Detection counts" },
  { key: "runmeta", label: "Metadata column" },
  { key: "data", label: "Detection column (aggregated)" },
];

function _ttFlagOptions(list, selected, valueKey, labelKey) {
  return list
    .map((o) => {
      const v = o[valueKey || "key"];
      const l = o[labelKey || "label"];
      return `<option value="${_ttFlagEsc(v)}"${String(v) === String(selected) ? " selected" : ""}>${_ttFlagEsc(
        l,
      )}</option>`;
    })
    .join("");
}

/** Options for a rule's field select.
 *  A rule can legitimately name a field the current catalogue does not list —
 *  a pipeline default targeting a metadata column that this run never
 *  populated, for instance. Dropping it would show the wrong field as
 *  selected (or an empty select) while the rule quietly kept evaluating the
 *  real one, so the rule's own field is always injected. */
function _ttFlagFieldOptions(list, rule) {
  const opts = list.slice();
  if (rule.field && !opts.some((f) => f.key === rule.field)) {
    opts.unshift({ key: rule.field, label: rule.field + " (not in this run)", num: true });
  }
  if (!opts.length) {
    return `<option value="">(no ${rule.source === "runmeta" ? "metadata columns" : "fields"} available yet)</option>`;
  }
  return _ttFlagOptions(opts, rule.field);
}

/** Operators that make sense for a field: numeric fields lose "contains",
 *  text fields lose the inequalities. Both keep equality and empty checks. */
function _ttFlagOpsFor(def) {
  const isNum = !def || def.num !== false;
  return TT_FLAG_OPS.filter((o) => o.num === null || o.num === isNum);
}

function ttFlagRenderRules() {
  const host = document.getElementById("flag-rule-list");
  if (!host) return;
  if (!TT_FLAG_RULES.length) {
    host.innerHTML =
      `<div class="flag-empty"><i class="fas fa-flag"></i> No rules yet. ` +
      `Add one below, or start from a preset.</div>`;
    return;
  }
  const fields = ttFlagFields();
  host.innerHTML = TT_FLAG_RULES.map((r, i) => {
    const list = fields[r.source] || [];
    const def = list.find((f) => f.key === r.field) || list[0] || { key: r.field, label: r.field, num: true };
    const ops = _ttFlagOpsFor(def);
    const opDef = TT_FLAG_OPS.find((o) => o.op === r.op) || {};
    const joiner =
      i === 0
        ? "<span class='flag-joiner lead'>Flag&nbsp;if</span>"
        : `<span class="flag-joiner">${TT_FLAG_LOGIC === "all" ? "and" : "or"}</span>`;
    return (
      `<div class="flag-rule${r.on ? "" : " off"}" data-rule="${_ttFlagEsc(r.id)}">` +
      joiner +
      `<input type="checkbox" class="flag-rule-on" ${r.on ? "checked" : ""} title="Enable / disable this rule" />` +
      `<select class="flag-rule-source" title="Where the value comes from">${_ttFlagOptions(
        _TT_FLAG_SOURCES,
        r.source,
      )}</select>` +
      `<select class="flag-rule-field" title="Which value to test">${_ttFlagFieldOptions(list, r)}</select>` +
      (r.source === "data"
        ? `<select class="flag-rule-agg" title="How the column is rolled up across the sample's detections">${_ttFlagOptions(
            TT_FLAG_AGGS,
            r.agg || "max",
          )}</select>`
        : "") +
      (def.needsTass
        ? `<input type="number" class="flag-rule-tass" value="${_ttFlagNum(
            r.tass,
          )}" min="0" max="100" step="1" title="TASS cutoff used when counting organisms" />`
        : "") +
      `<select class="flag-rule-op">${_ttFlagOptions(ops, r.op, "op", "label")}</select>` +
      (opDef.novalue
        ? `<span class="flag-rule-novalue">—</span>`
        : `<input type="${def.num === false ? "text" : "number"}" class="flag-rule-value" value="${_ttFlagEsc(
            r.value,
          )}" placeholder="value" />`) +
      `<select class="flag-rule-action" title="What happens to a sample that matches">` +
      `<option value="flag"${r.action !== "hide" ? " selected" : ""}>flag it</option>` +
      `<option value="hide"${r.action === "hide" ? " selected" : ""}>flag &amp; hide it</option></select>` +
      `<button type="button" class="flag-rule-del" title="Remove this rule"><i class="fas fa-trash"></i></button>` +
      `<span class="flag-rule-hits"></span>` +
      `</div>`
    );
  }).join("");
  _ttFlagUpdateRuleHits();
}

/** Per-rule "n samples" counter shown at the end of each row, so you can see
 *  which rule is doing the work without reading the preview list. */
function _ttFlagUpdateRuleHits() {
  const st = ttFlagEvaluate();
  const per = {};
  st.forEach((v) => {
    if (!v.flagged) return;
    v.hits.forEach((h) => {
      per[h.rule.id] = (per[h.rule.id] || 0) + 1;
    });
  });
  document.querySelectorAll("#flag-rule-list .flag-rule").forEach((row) => {
    const el = row.querySelector(".flag-rule-hits");
    if (!el) return;
    const n = per[row.dataset.rule] || 0;
    el.textContent = n ? `${n} sample${n === 1 ? "" : "s"}` : "";
    el.className = "flag-rule-hits" + (n ? " on" : "");
  });
}

/* ── Dialog: live preview of what the rules catch ────────────────────────── */
function ttFlagRenderPreview() {
  const host = document.getElementById("flag-preview");
  if (!host) return;
  const c = ttFlagCounts();
  const flagged = ttFlagFlaggedSamples();
  if (!c.rules) {
    host.innerHTML = `<div class="flag-preview-head quiet">No active rules — nothing is flagged.</div>`;
    return;
  }
  if (!flagged.length) {
    host.innerHTML =
      `<div class="flag-preview-head ok"><i class="fas fa-circle-check"></i> ` +
      `All ${c.total} sample${c.total === 1 ? "" : "s"} pass.</div>`;
    return;
  }
  host.innerHTML =
    `<div class="flag-preview-head warn"><i class="fas fa-flag"></i> ` +
    `<b>${flagged.length}</b> of ${c.total} sample${c.total === 1 ? "" : "s"} flagged` +
    (c.hidden ? ` · <b>${c.hidden}</b> will be hidden from every view` : "") +
    `</div><div class="flag-preview-list">` +
    flagged
      .map((s) => {
        const st = ttFlagStateFor(s) || { hits: [] };
        return (
          `<div class="flag-preview-row${st.hide ? " hidden" : ""}">` +
          `<div class="flag-preview-name"><i class="fas ${st.hide ? "fa-eye-slash" : "fa-flag"}"></i> ${_ttFlagEsc(
            s,
          )}</div>` +
          `<ul class="flag-preview-why">${st.hits.map((h) => `<li>${_ttFlagEsc(h.text)}</li>`).join("")}</ul></div>`
        );
      })
      .join("") +
    `</div>`;
}

/* ── Presets ─────────────────────────────────────────────────────────────── */
/*  Deliberately the same shapes the pipeline params produce, so a preset and
    a --report_flag_* default are indistinguishable once loaded.            */
const TT_FLAG_PRESETS = {
  lowreads: { source: "meta", field: "total_reads", op: "<", value: "100000", action: "flag" },
  noalign: { source: "meta", field: "aligned_reads", op: "<", value: "1", action: "flag" },
  fewtaxa: {
    source: "derived",
    field: "unique_taxids_above_tass",
    op: "<",
    value: "1",
    tass: 75,
    action: "flag",
  },
  control: { source: "meta", field: "control_type", op: "!empty", value: "", action: "flag" },
  nodetect: { source: "derived", field: "passing_detections", op: "<", value: "1", action: "flag" },
};

/* ── Wiring ──────────────────────────────────────────────────────────────── */
(function _ttFlagWireUI() {
  function _apply(opts) {
    if (typeof ttFlagRefresh === "function") ttFlagRefresh(opts);
    ttFlagRenderPreview();
    _ttFlagUpdateRuleHits();
  }

  let _t = null;
  function _applyDebounced() {
    clearTimeout(_t);
    _t = setTimeout(() => _apply(), 220);
  }

  function ready() {
    const openBtn = document.getElementById("flag-open-btn");
    const closeBtn = document.getElementById("flag-modal-close");
    const doneBtn = document.getElementById("flag-apply-btn");
    const overlay = document.getElementById("flag-modal-overlay");
    const list = document.getElementById("flag-rule-list");

    if (openBtn && !openBtn._wired) {
      openBtn._wired = true;
      openBtn.addEventListener("click", ttFlagOpenModal);
    }
    [closeBtn, doneBtn].forEach((b) => {
      if (b && !b._wired) {
        b._wired = true;
        b.addEventListener("click", ttFlagCloseModal);
      }
    });
    if (overlay && !overlay._wired) {
      overlay._wired = true;
      //  Click on the backdrop (not the card) closes, matching the other modals.
      overlay.addEventListener("click", (e) => {
        if (e.target === overlay) ttFlagCloseModal();
      });
      document.addEventListener("keydown", (e) => {
        if (e.key === "Escape" && overlay.style.display === "flex") ttFlagCloseModal();
      });
    }

    // ── Top-level controls ────────────────────────────────────────────
    const en = document.getElementById("flag-enabled");
    if (en && !en._wired) {
      en._wired = true;
      en.addEventListener("change", () => {
        TT_FLAG_ENABLED = en.checked;
        _apply();
      });
    }
    const lg = document.getElementById("flag-logic");
    if (lg && !lg._wired) {
      lg._wired = true;
      lg.addEventListener("change", () => {
        TT_FLAG_LOGIC = lg.value === "all" ? "all" : "any";
        ttFlagRenderRules(); // the and/or joiners between rows change
        _apply();
      });
    }
    const mf = document.getElementById("flag-missing-fails");
    if (mf && !mf._wired) {
      mf._wired = true;
      mf.addEventListener("change", () => {
        TT_FLAG_MISSING_FAILS = mf.checked;
        _apply();
      });
    }
    const ha = document.getElementById("flag-hide-all");
    if (ha && !ha._wired) {
      ha._wired = true;
      ha.addEventListener("change", () => {
        TT_FLAG_HIDE_ALL = ha.checked;
        _apply();
      });
    }
    const sideToggle = document.getElementById("flag-hide-toggle");
    if (sideToggle && !sideToggle._wired) {
      sideToggle._wired = true;
      sideToggle.addEventListener("change", () => {
        TT_FLAG_HIDE_ALL = sideToggle.checked;
        _ttFlagSyncTopControls();
        _apply();
      });
    }

    // ── Add / preset / reset / clear ──────────────────────────────────
    const addBtn = document.getElementById("flag-add-rule");
    if (addBtn && !addBtn._wired) {
      addBtn._wired = true;
      addBtn.addEventListener("click", () => {
        TT_FLAG_RULES.push(ttFlagNewRule());
        ttFlagRenderRules();
        _apply();
      });
    }
    document.querySelectorAll(".flag-preset").forEach((b) => {
      if (b._wired) return;
      b._wired = true;
      b.addEventListener("click", () => {
        const seed = TT_FLAG_PRESETS[b.dataset.preset];
        if (!seed) return;
        TT_FLAG_RULES.push(ttFlagNewRule(seed));
        ttFlagRenderRules();
        _apply();
      });
    });
    const resetBtn = document.getElementById("flag-reset-btn");
    if (resetBtn && !resetBtn._wired) {
      resetBtn._wired = true;
      resetBtn.addEventListener("click", () => {
        ttFlagLoadConfig((typeof BOOT !== "undefined" && BOOT && BOOT.sample_flags) || null);
        _ttFlagSyncTopControls();
        ttFlagRenderRules();
        _apply();
      });
    }
    const clearBtn = document.getElementById("flag-clear-btn");
    if (clearBtn && !clearBtn._wired) {
      clearBtn._wired = true;
      clearBtn.addEventListener("click", () => {
        TT_FLAG_RULES = [];
        ttFlagRenderRules();
        _apply();
      });
    }

    // ── Rule rows (delegated — rows are re-rendered constantly) ───────
    if (list && !list._wired) {
      list._wired = true;
      const ruleOf = (el) => {
        const row = el.closest(".flag-rule");
        if (!row) return null;
        return TT_FLAG_RULES.find((r) => r.id === row.dataset.rule) || null;
      };
      list.addEventListener("change", (e) => {
        const r = ruleOf(e.target);
        if (!r) return;
        const t = e.target;
        if (t.classList.contains("flag-rule-on")) r.on = t.checked;
        else if (t.classList.contains("flag-rule-source")) {
          r.source = t.value;
          //  Field lists differ per source — land on the first valid field and
          //  an operator that field actually supports.
          const fields = ttFlagFields()[r.source] || [];
          const def = fields[0] || { key: "", num: true };
          r.field = def.key;
          const ops = _ttFlagOpsFor(def);
          if (!ops.some((o) => o.op === r.op)) r.op = ops[0] ? ops[0].op : "==";
        } else if (t.classList.contains("flag-rule-field")) {
          r.field = t.value;
          const def = ttFlagFieldDef(r);
          const ops = _ttFlagOpsFor(def);
          if (!ops.some((o) => o.op === r.op)) r.op = ops[0] ? ops[0].op : "==";
        } else if (t.classList.contains("flag-rule-op")) r.op = t.value;
        else if (t.classList.contains("flag-rule-agg")) r.agg = t.value;
        else if (t.classList.contains("flag-rule-action")) r.action = t.value === "hide" ? "hide" : "flag";
        else if (t.classList.contains("flag-rule-tass")) r.tass = Number(t.value);
        else if (t.classList.contains("flag-rule-value")) r.value = t.value;
        ttFlagRenderRules(); // shape of the row can change (agg / tass / value)
        _apply();
      });
      //  Typing in a value or TASS box updates live, debounced.
      list.addEventListener("input", (e) => {
        const r = ruleOf(e.target);
        if (!r) return;
        if (e.target.classList.contains("flag-rule-value")) r.value = e.target.value;
        else if (e.target.classList.contains("flag-rule-tass")) r.tass = Number(e.target.value);
        else return;
        _applyDebounced();
      });
      list.addEventListener("click", (e) => {
        const del = e.target.closest(".flag-rule-del");
        if (!del) return;
        const row = del.closest(".flag-rule");
        if (!row) return;
        TT_FLAG_RULES = TT_FLAG_RULES.filter((r) => r.id !== row.dataset.rule);
        ttFlagRenderRules();
        _apply();
      });
    }

    // ── Clicking a flagged-sample chip anywhere scrolls the sidebar to it ──
    document.addEventListener("click", (e) => {
      const chip = e.target.closest(".tt-flag-badge[data-flag-sample]");
      if (!chip) return;
      const sid = chip.getAttribute("data-flag-sample");
      const row = document.querySelector(
        `.sample-entry[data-sid="${window.CSS && CSS.escape ? CSS.escape(sid) : sid}"]`,
      );
      if (row) {
        row.scrollIntoView({ block: "center", behavior: "smooth" });
        row.classList.add("flag-pulse");
        setTimeout(() => row.classList.remove("flag-pulse"), 1400);
      } else {
        ttFlagOpenModal();
      }
    });

    ttFlagRenderSummary();
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", ready);
  else ready();
})();
