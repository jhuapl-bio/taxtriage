/* ═══════════════════════════════════════════════════════════════════════════
       -  §  RESIZABLE DETAIL PANEL
       -     Drag handle for the inline VF/AMR per-organism detail panel.
       -     Persists user-chosen height to localStorage.
═══════════════════════════════════════════════════════════════════════════ */
(function () {
  // Use proper object (not boolean) so we can store panel ref without losing it
  let _drag = null,
    _rTimer = null;
  document.addEventListener("mousedown", (e) => {
    if (e.target && e.target.classList.contains("panel-resize-handle")) {
      const panel = e.target.closest('[id$="-detail-panel"]') || e.target.parentElement;
      if (!panel) return;
      _drag = { panel, startX: e.clientX, startW: panel.offsetWidth };
      document.body.style.cursor = "ew-resize";
      document.body.style.userSelect = "none";
      e.preventDefault();
    }
  });
  document.addEventListener("mousemove", (e) => {
    if (!_drag) return;
    const dx = _drag.startX - e.clientX;
    const newW = Math.max(280, Math.min(window.innerWidth * 0.8, _drag.startW + dx));
    _drag.panel.style.width = newW + "px";
    // Redraw comparison chart at new width (debounced)
    clearTimeout(_rTimer);
    _rTimer = setTimeout(() => {
      if (window._lastProtDetail) {
        const { g, o, r, gh } = window._lastProtDetail;
        _drawProtCompare(g, o, r, gh);
      }
    }, 200);
  });
  document.addEventListener("mouseup", () => {
    if (_drag) {
      _drag = null;
      document.body.style.cursor = "";
      document.body.style.userSelect = "";
    }
  });
})();
