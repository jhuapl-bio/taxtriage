/* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB REORDERING  (drag + drop, keyboard)
       -     Every tab in the top bar can be dragged left or right to change its
       -     place. Keyboard: focus a tab and press Alt+Shift+← / → to move it.
       -
       -     The order is a viewer preference saved in localStorage (one order
       -     for every report opened in this browser), applied on top of the
       -     report config's tabs.order (51_admin_panel.js, which runs first).
       -     Tabs the saved order does not know (added in a newer report) keep
       -     their default slot. A "↺ Reset order" button at the end of the bar appears once
       -     the order differs from the default and puts it back.
       -
       -     Pointer events (not HTML5 drag-and-drop) so it works on <button>s
       -     in every browser. A real drag swallows the click that follows, so
       -     dropping a tab never switches to it. Touch is left alone: on touch
       -     screens a horizontal swipe keeps scrolling the bar.
       -
       -     Public: window.ttTabOrder = { get, set, reset, custom }.
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  "use strict";
  const KEY = "taxtriage:tabOrder:v1";
  const DRAG_PX = 6; // movement before a press becomes a drag
  const EDGE_PX = 36; // auto-scroll zone at each end of the bar
  const bar = document.getElementById("tabbar");
  if (!bar) return;

  const tabs = () => Array.from(bar.querySelectorAll(".tab-btn"));
  const ids = () => tabs().map((b) => b.dataset.tab);
  const canon = (id) => (typeof ttCanonicalTabId === "function" ? ttCanonicalTabId(id) : id);
  const DEFAULT = ids(); // after the report config's order (51_admin_panel.js) has been applied
  const same = (a, b) => a.length === b.length && a.every((x, i) => x === b[i]);

  // ── reset button (not a .tab-btn, so tab switching / a11y code ignore it) ──
  const resetBtn = document.createElement("button");
  resetBtn.type = "button";
  resetBtn.className = "tab-order-reset";
  resetBtn.title = "Reset the tab order";
  resetBtn.setAttribute("aria-label", "Reset the tab order");
  resetBtn.textContent = "↺ Reset order"; // plain text: no icon-font dependency
  resetBtn.hidden = true;
  bar.appendChild(resetBtn);
  resetBtn.addEventListener("click", () => reset());

  function read() {
    try {
      const a = JSON.parse(localStorage.getItem(KEY) || "null");
      return Array.isArray(a) ? Array.from(new Set(a.map(canon))) : null;
    } catch (e) {
      return null;
    }
  }
  function write(order) {
    try {
      if (!order || same(order, DEFAULT)) localStorage.removeItem(KEY);
      else localStorage.setItem(KEY, JSON.stringify(order));
    } catch (e) {
      /* storage blocked: the order still holds for this page view */
    }
  }
  function syncReset() {
    resetBtn.hidden = same(ids(), DEFAULT);
  }

  // Saved ids that still exist, with any tab the saved order doesn't know
  // slotted back in at its default index.
  function merge(saved) {
    const have = new Set(DEFAULT);
    const out = saved.filter((id) => have.has(id));
    DEFAULT.forEach((id, i) => {
      if (!out.includes(id)) out.splice(Math.min(i, out.length), 0, id);
    });
    return out;
  }
  function apply(order) {
    const byId = {};
    tabs().forEach((b) => (byId[b.dataset.tab] = b));
    merge(order).forEach((id) => byId[id] && bar.insertBefore(byId[id], resetBtn));
    bar.appendChild(resetBtn);
    syncReset();
  }
  function set(order) {
    apply(order);
    write(ids());
  }
  function reset() {
    apply(DEFAULT);
    write(null);
    const active = bar.querySelector(".tab-btn.active");
    if (active) active.scrollIntoView({ block: "nearest", inline: "nearest" });
  }

  const saved = read();
  if (saved) apply(saved);

  // ── pointer drag ─────────────────────────────────────────────────────────
  let drag = null;
  let swallowClick = false;
  let scrollTimer = null;

  const visibleOthers = (b) => tabs().filter((x) => x !== b && x.offsetParent !== null);

  function placeAt(clientX) {
    const b = drag.btn;
    const others = visibleOthers(b);
    const next = others.find((x) => {
      const r = x.getBoundingClientRect();
      return clientX < r.left + r.width / 2;
    });
    const ref = next || resetBtn;
    if (b.nextElementSibling !== ref) bar.insertBefore(b, ref);
  }
  function autoScroll() {
    clearInterval(scrollTimer);
    scrollTimer = null;
    if (!drag || !drag.on) return;
    const r = bar.getBoundingClientRect();
    const dir = drag.x < r.left + EDGE_PX ? -1 : drag.x > r.right - EDGE_PX ? 1 : 0;
    if (!dir || bar.scrollWidth <= bar.clientWidth) return;
    scrollTimer = setInterval(() => {
      if (!drag) return clearInterval(scrollTimer);
      bar.scrollLeft += dir * 14;
      placeAt(drag.x);
    }, 30);
  }
  function startDrag() {
    drag.on = true;
    drag.btn.classList.add("tt-tab-dragging");
    document.body.classList.add("tt-tab-drag-active");
    try {
      bar.setPointerCapture(drag.pid);
    } catch (e) {}
    if (typeof hideTip === "function") hideTip();
  }
  function endDrag(commit) {
    clearInterval(scrollTimer);
    scrollTimer = null;
    if (!drag) return;
    const d = drag;
    drag = null;
    if (!d.on) return;
    d.btn.classList.remove("tt-tab-dragging");
    document.body.classList.remove("tt-tab-drag-active");
    try {
      bar.releasePointerCapture(d.pid);
    } catch (e) {}
    if (!commit) apply(d.before);
    else if (!same(ids(), d.before)) write(ids());
    syncReset();
    // The click that follows pointerup must not switch tabs.
    swallowClick = true;
    setTimeout(() => (swallowClick = false), 60);
  }

  bar.addEventListener("pointerdown", (e) => {
    if (e.button !== 0 || e.pointerType === "touch") return;
    const btn = e.target.closest(".tab-btn");
    if (!btn || !bar.contains(btn)) return;
    drag = { btn, pid: e.pointerId, x0: e.clientX, x: e.clientX, on: false, before: ids() };
  });
  bar.addEventListener("pointermove", (e) => {
    if (!drag || e.pointerId !== drag.pid) return;
    drag.x = e.clientX;
    if (!drag.on) {
      if (Math.abs(e.clientX - drag.x0) < DRAG_PX) return;
      startDrag();
    }
    e.preventDefault();
    placeAt(e.clientX);
    autoScroll();
  });
  bar.addEventListener("pointerup", (e) => drag && e.pointerId === drag.pid && endDrag(true));
  bar.addEventListener("pointercancel", () => endDrag(true));
  bar.addEventListener("lostpointercapture", () => drag && drag.on && endDrag(true));
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape" && drag && drag.on) endDrag(false);
  });
  // Capture phase on the bar runs before each tab's own click handler.
  bar.addEventListener(
    "click",
    (e) => {
      if (!swallowClick) return;
      swallowClick = false;
      e.stopImmediatePropagation();
      e.preventDefault();
    },
    true,
  );
  // Buttons don't start a native drag, but their icons / text selection can.
  bar.addEventListener("dragstart", (e) => e.preventDefault());

  // ── keyboard: Alt+Shift+←/→ moves the focused tab ────────────────────────
  // Capture phase so the tablist's own arrow-key focus handling
  // (15_tab_switching.js) doesn't also run.
  bar.addEventListener(
    "keydown",
    (e) => {
      if (!e.altKey || !e.shiftKey || (e.key !== "ArrowLeft" && e.key !== "ArrowRight")) return;
      const btn =
        document.activeElement && document.activeElement.closest && document.activeElement.closest(".tab-btn");
      if (!btn || !bar.contains(btn)) return;
      e.preventDefault();
      e.stopPropagation();
      const vis = tabs().filter((x) => x.offsetParent !== null);
      const i = vis.indexOf(btn);
      if (e.key === "ArrowLeft" && i > 0) bar.insertBefore(btn, vis[i - 1]);
      else if (e.key === "ArrowRight" && i < vis.length - 1) bar.insertBefore(btn, vis[i + 1].nextElementSibling);
      else return;
      write(ids());
      syncReset();
      btn.focus();
      btn.scrollIntoView({ block: "nearest", inline: "nearest" });
    },
    true,
  );

  // Hint on every tab's tooltip, once.
  tabs().forEach((b) => {
    const t = b.getAttribute("title");
    b.setAttribute("title", (t ? t + "\n\n" : "") + "Drag to reorder tabs (or Alt+Shift+←/→).");
  });

  window.ttTabOrder = {
    get: ids,
    set: set,
    reset: reset,
    // The viewer's saved order, or null when it matches the default.
    custom: () => (same(ids(), DEFAULT) ? null : ids()),
  };
})();
