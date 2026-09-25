/* ═══════════════════════════════════════════════════════════════════════════
       -  §  ORGANISM FLAGS — UI
       -     The sidebar summary + Highlight / Hide / Only control, and the
       -     rule builder dialog (#oflag-modal-overlay). The model lives in
       -     48_organism_flags.js; this file edits TT_OFLAG_RULES and calls
       -     ttOFlagRefresh() after every change. Edits apply live.
       -
       -     Layout of one rule card:
       -        [on] Rule 1 ··············· [flag it ▾] [n rows] [🗑]
       -        IF   [source] [field] [peer] [op] [value]      [×]
       -        AND  …
       -        + and condition
═══════════════════════════════════════════════════════════════════════════ */

/* ── Sidebar summary ─────────────────────────────────────────────────────── */
function ttOFlagRenderSummary() {
  const line = document.getElementById("oflag-summary-line");
  const btn = document.getElementById("oflag-open-btn");
  _ttOFlagSyncView();
  if (!line) return;
  const c = ttOFlagCounts();
  if (btn) btn.classList.toggle("has-flags", c.flagged > 0);
  const lvl = c.level.toLowerCase();
  if (!c.rules) {
    line.className = "flag-line quiet";
    line.innerHTML =
      `No organism rules active. <span class="flag-line-sub">Add one to flag or hide hits by reads, genus, ` +
      `shared ANI and more.</span>`;
  } else if (!c.total) {
    line.className = "flag-line quiet";
    line.innerHTML = `${c.rules} organism rule${c.rules === 1 ? "" : "s"} ready.`;
  } else if (!c.flagged) {
    line.className = "flag-line ok";
    line.innerHTML = `<i class="fas fa-circle-check"></i> No ${lvl} rows match ${c.rules} rule${
      c.rules === 1 ? "" : "s"
    }.`;
  } else {
    line.className = "flag-line warn";
    line.innerHTML =
      `<i class="fas fa-flag"></i> <b>${c.flagged}</b> of ${c.total.toLocaleString()} ${lvl} row${
        c.total === 1 ? "" : "s"
      } flagged` +
      (c.hidden ? ` · <b>${c.hidden}</b> ${c.onlyActive ? "hidden (unflagged)" : "hidden"}` : "") +
      `<span class="flag-line-sub">${c.organisms} distinct organism${c.organisms === 1 ? "" : "s"} · ${c.rules} rule${
        c.rules === 1 ? "" : "s"
      }</span>`;
  }
}

function _ttOFlagSyncView() {
  const st = typeof ttOFlagEvaluate === "function" ? ttOFlagEvaluate() : null;
  ["oflag-view-mode", "oflag-view-mode-modal"].forEach((id) => {
    const el = document.getElementById(id);
    if (!el) return;
    el.value = TT_OFLAG_VIEW;
    el.classList.toggle("filtering", TT_OFLAG_VIEW !== "all");
    const only = el.querySelector('option[value="only"]');
    if (only) only.disabled = !(st && st.anyFlagged);
  });
  const en = document.getElementById("oflag-enabled");
  if (en) en.checked = !!TT_OFLAG_ENABLED;
  const sw = document.getElementById("oflag-enabled-side");
  if (sw) sw.checked = !!TT_OFLAG_ENABLED;
}

/* ── Dialog open / close ─────────────────────────────────────────────────── */
function ttOFlagOpenModal() {
  const ov = document.getElementById("oflag-modal-overlay");
  if (!ov) return;
  ov.style.display = "flex";
  ov.setAttribute("aria-hidden", "false");
  _ttOFlagSyncView();
  ttOFlagRenderRules();
  ttOFlagRenderPreview();
}
function ttOFlagCloseModal() {
  const ov = document.getElementById("oflag-modal-overlay");
  if (!ov) return;
  ov.style.display = "none";
  ov.setAttribute("aria-hidden", "true");
}

/* ── Rule cards ──────────────────────────────────────────────────────────── */
const _TT_OFLAG_SOURCE_OPTS = [
  { key: "col", label: "Detection column" },
  { key: "derived", label: "In-sample context" },
];

function _ttOFlagOpts(list, selected, vk, lk) {
  return list
    .map((o) => {
      const v = o[vk || "key"];
      return `<option value="${_ttOFlagEsc(v)}"${String(v) === String(selected) ? " selected" : ""}${
        o.disabled ? " disabled" : ""
      }>${_ttOFlagEsc(o[lk || "label"])}</option>`;
    })
    .join("");
}

function _ttOFlagOpsFor(def) {
  if (!def || def.num == null) return TT_OFLAG_OPS.slice();
  return TT_OFLAG_OPS.filter((o) => o.num === null || o.num === def.num);
}

/** Does any row in this report carry ANI? Drives the "(needs --enable_matrix)"
 *  hint so an ANI rule that can never match is obvious. */
function _ttOFlagHasAni() {
  const src = (typeof DATA !== "undefined" && DATA) || [];
  for (let i = 0; i < src.length; i++) if (_ttOFlagAniAnnotated(src[i])) return true;
  return false;
}

function _ttOFlagCondHTML(c, i, fields, hasAni) {
  let list = fields[c.source] || [];
  if (c.source === "derived" && !hasAni) {
    list = list.map((f) => (f.ani ? Object.assign({}, f, { label: f.label + " — needs --enable_matrix" }) : f));
  }
  if (c.field && !list.some((f) => f.key === c.field)) {
    list = [{ key: c.field, label: c.field + " (not in this run)" }].concat(list);
  }
  const def = ttOFlagFieldDef(c) || {};
  const ops = _ttOFlagOpsFor(def);
  const opDef = TT_OFLAG_OPS.find((o) => o.op === c.op) || {};
  const placeholder = c.op === "in" || c.op === "!in" ? "a, b, c" : "value";
  return (
    `<div class="flag-rule oflag-cond" data-cond="${_ttOFlagEsc(c.id)}">` +
    `<span class="flag-joiner${i === 0 ? " lead" : ""}">${i === 0 ? "if" : "and"}</span>` +
    `<select class="oflag-source" title="Where the value comes from">${_ttOFlagOpts(
      _TT_OFLAG_SOURCE_OPTS,
      c.source,
    )}</select>` +
    `<select class="oflag-field flag-rule-field">${_ttOFlagOpts(list, c.field)}</select>` +
    (def.needsPeer
      ? `<select class="oflag-peer" title="'stronger hit' only counts a partner with more reads (then higher TASS), so of two near-identical references only the weaker one matches. ANI partners are the ones the pipeline recorded at or above --ani_threshold (default 95%).">${_ttOFlagOpts(
          TT_OFLAG_PEERS,
          c.peer || "stronger",
        )}</select>`
      : "") +
    `<select class="oflag-op">${_ttOFlagOpts(ops, c.op, "op", "label")}</select>` +
    (opDef.novalue
      ? `<span class="flag-rule-novalue">—</span>`
      : `<input type="${
          def.num === true && c.op !== "in" && c.op !== "!in" ? "number" : "text"
        }" class="oflag-value flag-rule-value" value="${_ttOFlagEsc(c.value)}" placeholder="${placeholder}" />`) +
    `<button type="button" class="flag-rule-del oflag-cond-del" title="Remove this condition"><i class="fas fa-xmark"></i></button>` +
    `</div>`
  );
}

function ttOFlagRenderRules() {
  const host = document.getElementById("oflag-rule-list");
  if (!host) return;
  if (!TT_OFLAG_RULES.length) {
    host.innerHTML =
      `<div class="flag-empty"><i class="fas fa-flag"></i> No organism rules yet. ` +
      `Add one below, or start from a preset.</div>`;
    return;
  }
  const fields = ttOFlagFields();
  const hasAni = _ttOFlagHasAni();
  host.innerHTML = TT_OFLAG_RULES.map(
    (r, i) =>
      `<div class="oflag-card${r.on ? "" : " off"}" data-rule="${_ttOFlagEsc(r.id)}">` +
      `<div class="oflag-card-head">` +
      `<label><input type="checkbox" class="oflag-rule-on" ${r.on ? "checked" : ""}/> Rule ${i + 1}</label>` +
      `<span class="oflag-card-spacer"></span>` +
      `<select class="oflag-action" title="What happens to a matching organism row">` +
      `<option value="flag"${r.action !== "hide" ? " selected" : ""}>flag it</option>` +
      `<option value="hide"${r.action === "hide" ? " selected" : ""}>flag &amp; hide it</option></select>` +
      `<span class="flag-rule-hits"></span>` +
      `<button type="button" class="flag-rule-del oflag-rule-del" title="Remove this rule"><i class="fas fa-trash"></i></button>` +
      `</div>` +
      r.conds.map((c, j) => _ttOFlagCondHTML(c, j, fields, hasAni)).join("") +
      `<button type="button" class="oflag-add-cond"><i class="fas fa-plus"></i> and condition</button>` +
      `</div>`,
  ).join("");
  _ttOFlagUpdateRuleHits();
}

function _ttOFlagUpdateRuleHits() {
  const st = ttOFlagEvaluate();
  const lvl = _ttOFlagViewLevel();
  const per = {};
  st.forEach((v) => {
    if ((v.row["Level"] || "Strain") !== lvl) return;
    v.hits.forEach((h) => (per[h.rule.id] = (per[h.rule.id] || 0) + 1));
  });
  document.querySelectorAll("#oflag-rule-list .oflag-card").forEach((card) => {
    const el = card.querySelector(".flag-rule-hits");
    if (!el) return;
    const n = per[card.dataset.rule] || 0;
    el.textContent = n ? `${n} row${n === 1 ? "" : "s"}` : "0 rows";
    el.className = "flag-rule-hits" + (n ? " on" : "");
  });
}

/* ── Live preview ────────────────────────────────────────────────────────── */
const _TT_OFLAG_PREVIEW_MAX = 150;

function ttOFlagRenderPreview() {
  const host = document.getElementById("oflag-preview");
  if (!host) return;
  const c = ttOFlagCounts();
  if (!c.rules) {
    host.innerHTML = `<div class="flag-preview-head quiet">No active rules — nothing is flagged.</div>`;
    return;
  }
  const st = ttOFlagEvaluate();
  const rows = [];
  st.forEach((v) => {
    if ((v.row["Level"] || "Strain") === c.level) rows.push(v);
  });
  if (!rows.length) {
    host.innerHTML = `<div class="flag-preview-head ok"><i class="fas fa-circle-check"></i> No ${c.level.toLowerCase()} rows match.</div>`;
    return;
  }
  rows.sort(
    (a, b) =>
      String(a.row["Specimen ID"]).localeCompare(String(b.row["Specimen ID"])) ||
      _ttOFlagNum(b.row["# Reads Aligned"]) - _ttOFlagNum(a.row["# Reads Aligned"]),
  );
  host.innerHTML =
    `<div class="flag-preview-head warn"><i class="fas fa-flag"></i> <b>${
      rows.length
    }</b> ${c.level.toLowerCase()} row${rows.length === 1 ? "" : "s"} flagged (${c.organisms} organism${
      c.organisms === 1 ? "" : "s"
    })` +
    (c.hidden ? ` · <b>${c.hidden}</b> hidden from every view` : "") +
    `</div><div class="flag-preview-list">` +
    rows
      .slice(0, _TT_OFLAG_PREVIEW_MAX)
      .map(
        (v) =>
          `<div class="flag-preview-row${v.hide ? " hidden" : ""}">` +
          `<div class="flag-preview-name"><i class="fas ${v.hide ? "fa-eye-slash" : "fa-flag"}"></i> <i>${_ttOFlagEsc(
            v.row["Detected Organism"],
          )}</i><div class="oflag-preview-sample">${_ttOFlagEsc(v.row["Specimen ID"])} · ${_ttOFlagFmt(
            _ttOFlagNum(v.row["# Reads Aligned"]),
          )} reads</div></div>` +
          `<ul class="flag-preview-why">${v.hits.map((h) => `<li>${_ttOFlagEsc(h.text)}</li>`).join("")}</ul></div>`,
      )
      .join("") +
    (rows.length > _TT_OFLAG_PREVIEW_MAX
      ? `<div class="flag-preview-row quiet">… and ${rows.length - _TT_OFLAG_PREVIEW_MAX} more</div>`
      : "") +
    `</div>`;
}

/* ── Presets ─────────────────────────────────────────────────────────────── */
/*  Same shapes as the pipeline's --report_org_flag_* params produce.        */
const TT_OFLAG_PRESETS = {
  lowreads: { conds: [{ source: "col", field: "# Reads Aligned", op: "<", value: "10" }] },
  lowtass: { conds: [{ source: "col", field: "TASS Score", op: "<", value: "50" }] },
  ani: { conds: [{ source: "derived", field: "ani_max", op: ">=", value: "99", peer: "stronger" }] },
  genuslow: {
    conds: [
      { source: "col", field: "Genus", op: "in", value: "" },
      { source: "col", field: "# Reads Aligned", op: "<", value: "50" },
    ],
  },
  genusminor: {
    conds: [
      { source: "derived", field: "genus_rank", op: ">", value: "1" },
      { source: "derived", field: "genus_share", op: "<", value: "5" },
    ],
  },
  k2only: {
    conds: [
      { source: "derived", field: "k2_ratio", op: "<", value: "0.05" },
      { source: "col", field: "K2 Reads", op: ">=", value: "50" },
    ],
  },
};

/* ── Wiring ──────────────────────────────────────────────────────────────── */
(function _ttOFlagWireUI() {
  function apply() {
    ttOFlagRefresh();
    ttOFlagRenderPreview();
    _ttOFlagUpdateRuleHits();
  }
  let _t = null;
  function applyDebounced() {
    clearTimeout(_t);
    _t = setTimeout(apply, 250);
  }
  const ruleOf = (el) => {
    const card = el.closest(".oflag-card");
    return card ? TT_OFLAG_RULES.find((r) => r.id === card.dataset.rule) || null : null;
  };
  const condOf = (rule, el) => {
    const row = el.closest(".oflag-cond");
    return rule && row ? rule.conds.find((c) => c.id === row.dataset.cond) || null : null;
  };
  /*  After a field / source change, land on an operator the field supports
      and clear a value that no longer makes sense (text → number).         */
  function fixOp(c) {
    const def = ttOFlagFieldDef(c);
    const ops = _ttOFlagOpsFor(def);
    if (!ops.some((o) => o.op === c.op)) {
      c.op = ops[0] ? ops[0].op : "==";
      c.value = "";
    }
  }

  function ready() {
    const on = (id, ev, fn) => {
      const el = document.getElementById(id);
      if (el && !el["_w_" + ev]) {
        el["_w_" + ev] = true;
        el.addEventListener(ev, fn);
      }
      return el;
    };
    on("oflag-open-btn", "click", ttOFlagOpenModal);
    on("oflag-modal-close", "click", ttOFlagCloseModal);
    on("oflag-apply-btn", "click", ttOFlagCloseModal);
    const ov = on("oflag-modal-overlay", "click", (e) => {
      if (e.target === e.currentTarget) ttOFlagCloseModal();
    });
    if (ov && !ov._escWired) {
      ov._escWired = true;
      document.addEventListener("keydown", (e) => {
        if (e.key === "Escape" && ov.style.display === "flex") ttOFlagCloseModal();
      });
    }
    ["oflag-enabled", "oflag-enabled-side"].forEach((id) =>
      on(id, "change", (e) => {
        TT_OFLAG_ENABLED = e.target.checked;
        apply();
      }),
    );
    ["oflag-view-mode", "oflag-view-mode-modal"].forEach((id) =>
      on(id, "change", (e) => {
        TT_OFLAG_VIEW = ttOFlagNormalizeView(e.target.value);
        apply();
      }),
    );
    on("oflag-add-rule", "click", () => {
      TT_OFLAG_RULES.push(ttOFlagNewRule());
      ttOFlagRenderRules();
      apply();
    });
    document.querySelectorAll(".oflag-preset").forEach((b) => {
      if (b._wired) return;
      b._wired = true;
      b.addEventListener("click", () => {
        const seed = TT_OFLAG_PRESETS[b.dataset.preset];
        if (!seed) return;
        TT_OFLAG_RULES.push(ttOFlagNewRule(JSON.parse(JSON.stringify(seed))));
        ttOFlagRenderRules();
        apply();
        //  Presets with a blank value (the genus list) want typing next.
        const cards = document.querySelectorAll("#oflag-rule-list .oflag-card");
        const last = cards[cards.length - 1];
        const blank = last && Array.from(last.querySelectorAll(".oflag-value")).find((i) => !i.value);
        if (blank) blank.focus();
      });
    });
    on("oflag-reset-btn", "click", () => {
      ttOFlagLoadConfig((typeof BOOT !== "undefined" && BOOT && BOOT.organism_flags) || null);
      ttOFlagRenderRules();
      apply();
    });
    on("oflag-clear-btn", "click", () => {
      TT_OFLAG_RULES = [];
      ttOFlagRenderRules();
      apply();
    });

    const list = document.getElementById("oflag-rule-list");
    if (list && !list._wired) {
      list._wired = true;
      list.addEventListener("change", (e) => {
        const t = e.target;
        const r = ruleOf(t);
        if (!r) return;
        if (t.classList.contains("oflag-rule-on")) r.on = t.checked;
        else if (t.classList.contains("oflag-action")) r.action = t.value === "hide" ? "hide" : "flag";
        else {
          const c = condOf(r, t);
          if (!c) return;
          if (t.classList.contains("oflag-source")) {
            c.source = t.value;
            const f = (ttOFlagFields()[c.source] || [])[0];
            c.field = f ? f.key : "";
            fixOp(c);
          } else if (t.classList.contains("oflag-field")) {
            c.field = t.value;
            fixOp(c);
          } else if (t.classList.contains("oflag-peer")) c.peer = t.value === "any" ? "any" : "stronger";
          else if (t.classList.contains("oflag-op")) c.op = t.value;
          else if (t.classList.contains("oflag-value")) c.value = t.value;
        }
        ttOFlagRenderRules();
        apply();
      });
      list.addEventListener("input", (e) => {
        if (!e.target.classList.contains("oflag-value")) return;
        const c = condOf(ruleOf(e.target), e.target);
        if (!c) return;
        c.value = e.target.value;
        applyDebounced();
      });
      list.addEventListener("click", (e) => {
        const r = ruleOf(e.target);
        if (!r) return;
        if (e.target.closest(".oflag-rule-del")) {
          TT_OFLAG_RULES = TT_OFLAG_RULES.filter((x) => x !== r);
        } else if (e.target.closest(".oflag-cond-del")) {
          const c = condOf(r, e.target);
          r.conds = r.conds.filter((x) => x !== c);
          //  A rule with no conditions left is gone, not "matches everything".
          if (!r.conds.length) TT_OFLAG_RULES = TT_OFLAG_RULES.filter((x) => x !== r);
        } else if (e.target.closest(".oflag-add-cond")) {
          r.conds.push(ttOFlagNewCond({ source: "col", field: "TASS Score", op: "<", value: "50" }));
        } else return;
        ttOFlagRenderRules();
        apply();
      });
    }
    ttOFlagRenderSummary();
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", ready);
  else ready();
})();
