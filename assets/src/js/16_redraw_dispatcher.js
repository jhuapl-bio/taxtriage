      /* ═══════════════════════════════════════════════════════════════════════════
       -  §  REDRAW DISPATCHER
       -     _TAB_DIRTY  — { tab → bool } flags tracking which tabs need a
       -                   re-render after state changes.
       -     _drawTab()  — calls the right draw* function for a single tab.
       -     redraw()    — invalidates filteredData() cache, marks ALL tabs
       -                   dirty, then renders only the currently active tab.
       -     This is what makes filter changes feel instant — hidden tabs are
       -     not redrawn until the user switches to them.
═══════════════════════════════════════════════════════════════════════════ */

