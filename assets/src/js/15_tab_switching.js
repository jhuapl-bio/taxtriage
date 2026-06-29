      /* ═══════════════════════════════════════════════════════════════════════════
       -  §  TAB SWITCHING
       -     Wires every .tab-btn click → toggles the .active class on the
       -     button + matching .pane, then calls _drawTab(activeTab) only if
       -     that tab is currently dirty (lazy render — see REDRAW DISPATCHER).
═══════════════════════════════════════════════════════════════════════════ */
      document.querySelectorAll(".tab-btn").forEach((btn) => {
        btn.addEventListener("click", () => {
          if (btn.classList.contains("tab-disabled")) return;
          document.querySelectorAll(".tab-btn").forEach((b) => b.classList.remove("active"));
          document.querySelectorAll(".pane").forEach((p) => p.classList.remove("active"));
          btn.classList.add("active");
          activeTab = btn.dataset.tab;
          document.getElementById(`pane-${activeTab}`).classList.add("active");
          // Render the new tab if its data is stale OR it has never been drawn.
          // The "never drawn" guard is critical: if report init threw before
          // redraw() marked tabs dirty, the dirty flag is still its initial false
          // and the pane would otherwise render nothing on click.
          if (_TAB_DIRTY[activeTab] || !_TAB_RENDERED[activeTab]) _drawTab(activeTab);

          // ── Tab-specific init (runs after pane is visible) ──────────────
          if (activeTab === "map") {
            // Delay so the pane layout is computed before Leaflet measures the container
            setTimeout(() => {
              if (typeof L === "undefined") {
                const c = document.getElementById("map-container");
                if (c)
                  c.innerHTML =
                    '<p style="padding:2em;text-align:center;color:#888">' +
                    '<i class="fas fa-triangle-exclamation" style="font-size:2em;display:block;margin-bottom:.5em;opacity:.4"></i>' +
                    "Map requires an internet connection to load the Leaflet library.</p>";
                return;
              }
              if (typeof _initMap === "function") _initMap();
            }, 80);
          } else if (activeTab === "runmeta") {
            // Rebuild the metadata table and update sub-tab enabled states
            if (typeof _buildRunMetaTable === "function") _buildRunMetaTable();
            if (typeof _updateMetaSubTabStates === "function") _updateMetaSubTabStates();
          }
        });
      });

      /* ══════════════════════════════════════════════════════════════════════
         Tab bar accessibility + keyboard navigation.
         Additive only — the click handler above still owns activation. This
         layers on the ARIA tablist pattern (roles, aria-selected) and the
         expected keyboard behaviour (←/→ move between tabs, Home/End jump to
         first/last, Enter/Space activate). Hidden + disabled tabs are skipped.
         A MutationObserver keeps roving tabindex / aria-selected in sync no
         matter how a tab is activated or shown/hidden elsewhere in the app.
         ══════════════════════════════════════════════════════════════════════ */
      (function () {
        const bar = document.getElementById("tabbar");
        if (!bar) return;
        bar.setAttribute("role", "tablist");
        bar.setAttribute("aria-label", "Analysis views");

        const allTabs = () => Array.from(bar.querySelectorAll(".tab-btn"));
        // Tabs a keyboard user can actually land on: visible + not disabled.
        const navTabs = () =>
          allTabs().filter(
            (b) => !b.classList.contains("hidden") && !b.classList.contains("tab-disabled") && b.offsetParent !== null,
          );

        function syncA11y() {
          allTabs().forEach((b) => {
            b.setAttribute("role", "tab");
            const sel = b.classList.contains("active");
            b.setAttribute("aria-selected", sel ? "true" : "false");
            // Roving tabindex: only the active tab is in the tab order.
            b.tabIndex = sel ? 0 : -1;
            const pane = document.getElementById(`pane-${b.dataset.tab}`);
            if (pane) {
              b.setAttribute("aria-controls", `pane-${b.dataset.tab}`);
              pane.setAttribute("role", "tabpanel");
            }
          });
          // Guarantee at least one focusable tab even if none is active yet.
          if (!allTabs().some((b) => b.tabIndex === 0)) {
            const first = navTabs()[0];
            if (first) first.tabIndex = 0;
          }
        }

        function focusTab(tab) {
          if (!tab) return;
          tab.tabIndex = 0;
          tab.focus();
        }

        bar.addEventListener("keydown", (e) => {
          const tabs = navTabs();
          if (!tabs.length) return;
          const current = document.activeElement;
          let idx = tabs.indexOf(current);
          if (idx === -1) idx = tabs.findIndex((b) => b.classList.contains("active"));
          let next = null;
          switch (e.key) {
            case "ArrowRight":
            case "ArrowDown":
              next = tabs[(idx + 1) % tabs.length];
              break;
            case "ArrowLeft":
            case "ArrowUp":
              next = tabs[(idx - 1 + tabs.length) % tabs.length];
              break;
            case "Home":
              next = tabs[0];
              break;
            case "End":
              next = tabs[tabs.length - 1];
              break;
            case "Enter":
            case " ":
              if (current && current.classList.contains("tab-btn")) {
                e.preventDefault();
                current.click();
              }
              return;
            default:
              return;
          }
          if (next) {
            e.preventDefault();
            focusTab(next);
            next.click(); // activate-on-focus, matching the existing click flow
          }
        });

        // Re-sync whenever a tab's class (active / hidden / disabled) changes.
        new MutationObserver(syncA11y).observe(bar, {
          subtree: true,
          attributes: true,
          attributeFilter: ["class", "style"],
        });
        syncA11y();
      })();

      /* ── Summary inner sub-tab switching (Detections / Cross-Sample) ── */
      (function () {
        document.querySelectorAll(".sum-inner-tab").forEach((btn) => {
          btn.addEventListener("click", () => {
            const target = btn.dataset.inner;
            document.querySelectorAll(".sum-inner-tab").forEach((b) => b.classList.remove("active"));
            document.querySelectorAll(".sum-inner-pane").forEach((p) => p.classList.remove("active"));
            btn.classList.add("active");
            const pane = document.getElementById(`sum-inner-${target}`);
            if (pane) pane.classList.add("active");
            // Trigger xs redraw if switching to cross-sample for the first time
            if (target === "xs" && typeof _drawXS === "function") _drawXS();
          });
        });
      })();

      /* ── Run Metadata analysis sub-tab click handlers ── */
      (function () {
        document.querySelectorAll(".meta-subtab").forEach((btn) => {
          btn.addEventListener("click", () => {
            if (btn.disabled) return;
            if (typeof _switchMetaSub === "function") _switchMetaSub(btn.getAttribute("data-metasub"));
          });
        });
      })();

