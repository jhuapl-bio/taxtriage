/* ═══════════════════════════════════════════════════════════════════════════
       -  §  MULTI-SELECT DROPDOWN  (shared control)
       -     A collapsed button that opens a searchable checkbox list.
       -
       -     Replaces the long chip rows that used to sit above several views.
       -     Chips are fine for five things and unusable for fifty: they wrap
       -     into a wall that pushes the actual chart off-screen, and there is
       -     no way to find one item without reading all of them. A dropdown
       -     costs one click and stays one line high whatever the run size.
       -
       -     Used by:
       -       • the group legends (shared Group-by bar, map, group heatmap,
       -         group network) — see _mgRenderLegend in 38_meta_grouping.js
       -       • the Longitudinal organism picker — see 32_longitudinal_second.js
       -
       -     The host element is re-rendered wholesale by its owner on every
       -     state change, so THIS module owns the open/query state (keyed by
       -     host id) and re-applies it after each render. Without that, ticking
       -     a box would close the panel you are working in.
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  const _state = new Map(); // hostId → { open, query }

  function _st(id) {
    let s = _state.get(id);
    if (!s) _state.set(id, (s = { open: false, query: "" }));
    return s;
  }

  function _esc(s) {
    return String(s == null ? "" : s).replace(
      /[&<>"']/g,
      (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c],
    );
  }

  /** Render a multi-select into `host`.
   *
   *  opts.items      [{ key, label, checked, count, swatchHTML, title,
   *                     rowClass, stateIcon }]
   *  opts.summary    (nChecked, nTotal) → button label text
   *  opts.title      tooltip for the button
   *  opts.placeholder search-box placeholder
   *  opts.emptyText  shown when there are no items at all
   *  opts.onToggle   (key, checked, ev)  — a checkbox changed
   *  opts.onAll / opts.onNone            — the All / None links
   *  opts.onRowClick (key, ev)           — clicking the row body, not the box
   *                                        (the group legends use this for the
   *                                        highlight cycle)
   *  opts.actions    [{ label, title, className, onClick }] rendered beside the
   *                  button (e.g. "Reset (3)")
   */
  function render(host, opts) {
    const el = typeof host === "string" ? document.getElementById(host) : host;
    if (!el) return;
    const id = el.id || (el.id = "tt-ms-" + Math.abs(hashString(String(opts.key || Math.random()))));
    const o = opts || {};
    const items = o.items || [];
    const st = _st(id);

    if (!items.length && o.hideWhenEmpty !== false && !o.emptyText) {
      el.innerHTML = "";
      el.style.display = "none";
      return;
    }
    el.style.display = "";
    el.classList.add("tt-ms");

    const nChecked = items.filter((i) => i.checked).length;
    const label = o.summary ? o.summary(nChecked, items.length) : `${nChecked} of ${items.length} selected`;
    const q = st.query.trim().toLowerCase();
    const shown = q ? items.filter((i) => String(i.label).toLowerCase().includes(q)) : items;

    el.innerHTML =
      `<div class="tt-ms-wrap">` +
      `<button type="button" class="tt-ms-btn${st.open ? " open" : ""}" aria-expanded="${st.open}" ` +
      `title="${_esc(o.title || "")}">` +
      (o.icon ? `<i class="fas ${_esc(o.icon)}" aria-hidden="true"></i>` : "") +
      `<span class="tt-ms-label">${_esc(label)}</span>` +
      `<svg class="tt-ms-caret" width="10" height="10" viewBox="0 0 10 10" fill="none" aria-hidden="true">` +
      `<path d="M2 3.5l3 3 3-3" stroke="currentColor" stroke-width="1.5" stroke-linecap="round"/></svg>` +
      `</button>` +
      `<div class="tt-ms-pop${st.open ? " open" : ""}" role="group">` +
      `<div class="tt-ms-head">` +
      `<input type="search" class="tt-ms-search" placeholder="${_esc(o.placeholder || "Search…")}" ` +
      `value="${_esc(st.query)}" aria-label="Filter the list" />` +
      `<div class="tt-ms-quick">` +
      (o.onAll ? `<button type="button" class="tt-ms-all">All</button>` : "") +
      (o.onNone ? `<button type="button" class="tt-ms-none">None</button>` : "") +
      `</div></div>` +
      `<div class="tt-ms-list">` +
      (shown.length
        ? shown
            .map(
              (i) =>
                `<label class="tt-ms-opt${i.rowClass ? " " + i.rowClass : ""}" title="${_esc(i.title || i.label)}" ` +
                `data-key="${_esc(i.key)}">` +
                `<input type="checkbox" class="tt-ms-cb" data-key="${_esc(i.key)}"${i.checked ? " checked" : ""} />` +
                (i.stateIcon ? `<i class="fas ${_esc(i.stateIcon)} tt-ms-state"></i>` : "") +
                (i.swatchHTML || "") +
                `<span class="tt-ms-opt-label">${_esc(i.label)}</span>` +
                (i.count != null ? `<span class="tt-ms-opt-count">${_esc(i.count)}</span>` : "") +
                `</label>`,
            )
            .join("")
        : `<div class="tt-ms-empty">${_esc(q ? "No matches." : o.emptyText || "Nothing to show.")}</div>`) +
      `</div>` +
      (o.footNote ? `<div class="tt-ms-foot">${o.footNote}</div>` : "") +
      `</div>` +
      `</div>` +
      (o.actions || [])
        .map(
          (a, n) =>
            `<button type="button" class="tt-ms-action${a.className ? " " + a.className : ""}" data-act="${n}" ` +
            `title="${_esc(a.title || "")}">${a.label}</button>`,
        )
        .join("");

    // ── wiring ──────────────────────────────────────────────────────────
    const btn = el.querySelector(".tt-ms-btn");
    const pop = el.querySelector(".tt-ms-pop");
    btn.addEventListener("click", (ev) => {
      ev.stopPropagation();
      st.open = !st.open;
      pop.classList.toggle("open", st.open);
      btn.classList.toggle("open", st.open);
      btn.setAttribute("aria-expanded", String(st.open));
      if (st.open) {
        const s = el.querySelector(".tt-ms-search");
        if (s) s.focus();
      }
    });

    const search = el.querySelector(".tt-ms-search");
    if (search) {
      search.addEventListener("input", () => {
        st.query = search.value;
        render(el, opts); // re-render the list; open state survives via `st`
      });
      // Clicks inside the panel must not bubble to the document closer below.
      search.addEventListener("click", (ev) => ev.stopPropagation());
    }
    if (pop) pop.addEventListener("click", (ev) => ev.stopPropagation());

    el.querySelectorAll(".tt-ms-cb").forEach((cb) => {
      cb.addEventListener("change", (ev) => {
        ev.stopPropagation();
        if (typeof o.onToggle === "function") o.onToggle(cb.getAttribute("data-key"), cb.checked, ev);
      });
      cb.addEventListener("click", (ev) => ev.stopPropagation());
    });

    //  Clicking the row (but not its checkbox) is a secondary action. The
    //  group legends use it for the normal → highlight → hidden cycle, so the
    //  three-state behaviour survives the move off chips.
    if (typeof o.onRowClick === "function") {
      el.querySelectorAll(".tt-ms-opt").forEach((row) => {
        row.addEventListener("click", (ev) => {
          if (ev.target.classList.contains("tt-ms-cb")) return;
          ev.preventDefault();
          ev.stopPropagation();
          o.onRowClick(row.getAttribute("data-key"), ev);
        });
      });
    }

    (o.actions || []).forEach((a, n) => {
      const b = el.querySelector(`.tt-ms-action[data-act="${n}"]`);
      if (b && typeof a.onClick === "function") b.addEventListener("click", a.onClick);
    });
  }

  function hashString(s) {
    let h = 0;
    for (let i = 0; i < s.length; i++) h = (h * 31 + s.charCodeAt(i)) | 0;
    return h;
  }

  /** Close every open panel — one document listener for all instances. */
  document.addEventListener("click", () => {
    let any = false;
    _state.forEach((s) => {
      if (s.open) {
        s.open = false;
        any = true;
      }
    });
    if (!any) return;
    document.querySelectorAll(".tt-ms-pop.open").forEach((p) => p.classList.remove("open"));
    document.querySelectorAll(".tt-ms-btn.open").forEach((b) => {
      b.classList.remove("open");
      b.setAttribute("aria-expanded", "false");
    });
  });
  document.addEventListener("keydown", (ev) => {
    if (ev.key !== "Escape") return;
    _state.forEach((s) => (s.open = false));
    document.querySelectorAll(".tt-ms-pop.open").forEach((p) => p.classList.remove("open"));
    document.querySelectorAll(".tt-ms-btn.open").forEach((b) => b.classList.remove("open"));
  });

  window.ttMultiSelect = { render, isOpen: (id) => !!_st(id).open };
})();
