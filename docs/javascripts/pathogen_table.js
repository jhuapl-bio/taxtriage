/* Interactive browser for assets/pathogen_sheet.csv.
 *
 * Data is fetched once as docs/data/pathogen_sheet.json (columnar arrays) and
 * everything after that is client side: full-text search, multi-select facets,
 * a dependent taxonomy drilldown, sorting, pagination and CSV export of the
 * current selection.
 *
 * Facet counts are computed "excluding self" — the counts shown next to the
 * options of facet X reflect every active filter *except* X. That is what makes
 * multi-select feel right: options in the group you are working in don't
 * vanish as you tick them.
 */

(function () {
  "use strict";

  // Pages are served at directory URLs (/pathogen-sheet/), so a bare relative
  // path would resolve inside the page directory rather than the site root.
  // MkDocs rewrites this script's own src per page, so derive the data URL from
  // it — correct at any nesting depth and under any base path.
  var self =
    document.currentScript || document.querySelector('script[src*="pathogen_table.js"]');
  var DATA_URL = "data/pathogen_sheet.json";
  if (self && self.src) {
    var derived = self.src.replace(
      /javascripts\/pathogen_table\.js(\?.*)?$/,
      "data/pathogen_sheet.json"
    );
    if (derived !== self.src) DATA_URL = derived;
  }

  // Facets rendered as checkbox lists: [key, label, isMultiValue, searchable]
  var FACETS = [
    ["general_classification", "Classification", false, false],
    ["status", "Evidence status", false, false],
    ["high_consequence", "High consequence", false, false],
    ["kingdom", "Kingdom", false, false],
    ["mol_type", "Molecule type", false, false],
    ["pathogenic_sites", "Pathogenic sites", true, false],
    ["commensal_sites", "Commensal sites", true, false],
    ["host_organism", "Host organism", false, true],
  ];

  // Dependent dropdowns, ordered broad -> narrow.
  var TAXONOMY = ["phylum", "class", "order", "family", "genus"];

  var TABLE_COLS = [
    ["name", "Organism"],
    ["taxid", "Tax ID"],
    ["general_classification", "Class."],
    ["status", "Status"],
    ["high_consequence", "HC"],
    ["pathogenic_sites", "Pathogenic sites"],
    ["commensal_sites", "Commensal sites"],
    ["genus", "Genus"],
    ["family", "Family"],
    ["mol_type", "Mol"],
    ["host_organism", "Host"],
  ];

  var DETAIL_COLS = [
    ["taxid", "Tax ID"],
    ["general_classification", "Classification"],
    ["status", "Evidence status"],
    ["high_consequence", "High consequence"],
    ["alternative_names", "Also known as"],
    ["pathogenic_sites", "Pathogenic sites"],
    ["commensal_sites", "Commensal sites"],
    ["pathology", "Pathology"],
    ["host_organism", "Host organism"],
    ["kingdom", "Kingdom"],
    ["phylum", "Phylum"],
    ["class", "Class"],
    ["order", "Order"],
    ["family", "Family"],
    ["genus", "Genus"],
    ["mol_type", "Molecule type"],
    ["assembly_accession", "Assembly"],
    ["reference", "Reference"],
    ["Additional references", "Additional references"],
  ];

  var SEARCH_COLS = [
    "name",
    "alternative_names",
    "taxid",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "pathology",
    "assembly_accession",
  ];

  function esc(s) {
    return String(s == null ? "" : s).replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }

  function displayValue(v) {
    if (Array.isArray(v)) return v.join(", ");
    if (v === true) return "TRUE";
    if (v === false) return "FALSE";
    return v == null ? "" : String(v);
  }

  function init(root) {
    if (!root || root.dataset.ptReady === "1") return;
    root.dataset.ptReady = "1";

    fetch(DATA_URL)
      .then(function (r) {
        if (!r.ok) throw new Error("HTTP " + r.status);
        return r.json();
      })
      .then(function (payload) {
        build(root, payload);
      })
      .catch(function (err) {
        root.innerHTML =
          '<div class="pt-empty">Could not load the pathogen sheet (' +
          esc(err.message) +
          ").</div>";
      });
  }

  function build(root, payload) {
    var cols = payload.cols;
    var idx = {};
    cols.forEach(function (c, i) {
      idx[c] = i;
    });

    // Rows as objects keyed by column name, plus a precomputed lowercase blob
    // for search so we aren't re-joining strings on every keystroke.
    var rows = payload.rows.map(function (r) {
      var o = {};
      cols.forEach(function (c, i) {
        o[c] = r[i];
      });
      o.__blob = SEARCH_COLS.map(function (c) {
        return displayValue(o[c]);
      })
        .join(" ")
        .toLowerCase();
      return o;
    });

    var state = {
      q: "",
      facets: {}, // key -> Set of selected values
      tax: {}, // rank -> selected value
      sort: { col: "name", dir: 1 },
      page: 1,
      perPage: 50,
    };
    FACETS.forEach(function (f) {
      state.facets[f[0]] = new Set();
    });

    root.innerHTML = template();

    var el = {
      search: root.querySelector(".pt-search input"),
      count: root.querySelector(".pt-count"),
      chips: root.querySelector(".pt-chips"),
      facets: root.querySelector(".pt-facets"),
      tbody: root.querySelector(".pt-table tbody"),
      thead: root.querySelector(".pt-table thead tr"),
      pager: root.querySelector(".pt-pager"),
      reset: root.querySelector(".pt-reset"),
      exportBtn: root.querySelector(".pt-export"),
      drawer: root.querySelector(".pt-drawer"),
    };

    /* ── filtering ────────────────────────────────────────────────────── */

    function matchesFacet(row, key, multi) {
      var sel = state.facets[key];
      if (!sel.size) return true;
      var v = row[key];
      if (multi) {
        for (var i = 0; i < v.length; i++) if (sel.has(v[i])) return true;
        return false;
      }
      return sel.has(displayValue(v));
    }

    // `skip` lets us compute counts for a facet while ignoring its own selection.
    function filter(skip) {
      var q = state.q.trim().toLowerCase();
      var terms = q ? q.split(/\s+/) : [];
      return rows.filter(function (row) {
        for (var t = 0; t < terms.length; t++) {
          if (row.__blob.indexOf(terms[t]) === -1) return false;
        }
        for (var i = 0; i < FACETS.length; i++) {
          var key = FACETS[i][0];
          if (key === skip) continue;
          if (!matchesFacet(row, key, FACETS[i][2])) return false;
        }
        for (var r = 0; r < TAXONOMY.length; r++) {
          var rank = TAXONOMY[r];
          if (rank === skip) continue;
          var want = state.tax[rank];
          if (want && row[rank] !== want) return false;
        }
        return true;
      });
    }

    function countValues(subset, key, multi) {
      var counts = new Map();
      subset.forEach(function (row) {
        var vals = multi ? row[key] : [displayValue(row[key])];
        vals.forEach(function (v) {
          if (v === "") return;
          counts.set(v, (counts.get(v) || 0) + 1);
        });
      });
      return counts;
    }

    /* ── rendering ────────────────────────────────────────────────────── */

    function renderFacets(current) {
      var html = "";

      FACETS.forEach(function (f) {
        var key = f[0],
          label = f[1],
          multi = f[2],
          searchable = f[3];
        var counts = countValues(filter(key), key, multi);
        var sel = state.facets[key];
        var entries = Array.from(counts.entries()).sort(function (a, b) {
          return b[1] - a[1] || a[0].localeCompare(b[0]);
        });
        var open = sel.size > 0 || FACETS.indexOf(f) < 3;

        html +=
          '<details class="pt-facet" data-key="' +
          esc(key) +
          '"' +
          (open ? " open" : "") +
          "><summary>" +
          esc(label) +
          (sel.size ? " (" + sel.size + ")" : "") +
          '</summary><div class="pt-facet-body">';

        if (searchable && entries.length > 8) {
          html +=
            '<input class="pt-facet-search" type="search" placeholder="Filter ' +
            esc(label.toLowerCase()) +
            '…" aria-label="Filter ' +
            esc(label) +
            '">';
        }

        entries.forEach(function (e) {
          var checked = sel.has(e[0]);
          html +=
            '<label class="pt-opt' +
            (e[1] === 0 && !checked ? " is-zero" : "") +
            '" data-v="' +
            esc(e[0].toLowerCase()) +
            '"><input type="checkbox" value="' +
            esc(e[0]) +
            '"' +
            (checked ? " checked" : "") +
            '><span class="pt-opt-label">' +
            esc(e[0]) +
            '</span><span class="pt-opt-n">' +
            e[1] +
            "</span></label>";
        });

        html += "</div></details>";
      });

      // Taxonomy: each level's options are constrained by the level above.
      html +=
        '<details class="pt-facet pt-tax" open><summary>Taxonomy</summary><div class="pt-facet-body">';
      TAXONOMY.forEach(function (rank, i) {
        var parentChosen = i === 0 || state.tax[TAXONOMY[i - 1]];
        var counts = countValues(filter(rank), rank, false);
        var opts = Array.from(counts.keys()).sort();
        html +=
          '<select data-rank="' +
          esc(rank) +
          '"' +
          (parentChosen ? "" : " disabled") +
          ' aria-label="' +
          esc(rank) +
          '"><option value="">' +
          rank.charAt(0).toUpperCase() +
          rank.slice(1) +
          " — any</option>";
        opts.forEach(function (o) {
          html +=
            '<option value="' +
            esc(o) +
            '"' +
            (state.tax[rank] === o ? " selected" : "") +
            ">" +
            esc(o) +
            " (" +
            counts.get(o) +
            ")</option>";
        });
        html += "</select>";
      });
      html += "</div></details>";

      el.facets.innerHTML = html;
    }

    function renderChips() {
      var chips = [];
      FACETS.forEach(function (f) {
        state.facets[f[0]].forEach(function (v) {
          chips.push(
            '<button class="pt-chip" data-type="facet" data-key="' +
              esc(f[0]) +
              '" data-v="' +
              esc(v) +
              '">' +
              esc(f[1]) +
              ": " +
              esc(v) +
              "</button>"
          );
        });
      });
      TAXONOMY.forEach(function (rank) {
        if (state.tax[rank]) {
          chips.push(
            '<button class="pt-chip" data-type="tax" data-key="' +
              esc(rank) +
              '">' +
              esc(rank) +
              ": " +
              esc(state.tax[rank]) +
              "</button>"
          );
        }
      });
      el.chips.innerHTML = chips.join("");
      el.reset.hidden = chips.length === 0 && !state.q;
    }

    function highlight(text) {
      var q = state.q.trim();
      if (!q) return esc(text);
      var out = esc(text);
      q.split(/\s+/).forEach(function (term) {
        if (!term) return;
        var re = new RegExp("(" + term.replace(/[.*+?^${}()|[\]\\]/g, "\\$&") + ")", "gi");
        out = out.replace(re, '<mark class="pt-mark">$1</mark>');
      });
      return out;
    }

    function cell(row, key) {
      var v = row[key];
      if (key === "name") {
        return '<td class="pt-name">' + highlight(v) + "</td>";
      }
      if (key === "taxid") {
        return (
          '<td><a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=' +
          esc(v) +
          '" target="_blank" rel="noopener">' +
          esc(v) +
          "</a></td>"
        );
      }
      if (key === "general_classification" || key === "status") {
        return v
          ? '<td><span class="pt-tag" data-v="' + esc(v.toLowerCase()) + '">' + esc(v) + "</span></td>"
          : "<td></td>";
      }
      if (key === "high_consequence") {
        return v ? '<td class="pt-hc">yes</td>' : "<td></td>";
      }
      if (key === "pathogenic_sites" || key === "commensal_sites") {
        return '<td class="pt-sites">' + esc(v.join(", ")) + "</td>";
      }
      return "<td>" + highlight(displayValue(v)) + "</td>";
    }

    var currentRows = [];

    function renderTable() {
      var s = state.sort;
      currentRows.sort(function (a, b) {
        var av = displayValue(a[s.col]).toLowerCase();
        var bv = displayValue(b[s.col]).toLowerCase();
        if (s.col === "taxid") return (Number(a.taxid) - Number(b.taxid)) * s.dir;
        if (av === bv) return a.name.localeCompare(b.name);
        return av < bv ? -s.dir : s.dir;
      });

      var total = currentRows.length;
      var pages = Math.max(1, Math.ceil(total / state.perPage));
      if (state.page > pages) state.page = pages;
      var start = (state.page - 1) * state.perPage;
      var slice = currentRows.slice(start, start + state.perPage);

      el.count.innerHTML =
        "<strong>" + total.toLocaleString() + "</strong> of " + rows.length.toLocaleString() + " organisms";

      if (!total) {
        el.tbody.innerHTML =
          '<tr><td colspan="' + TABLE_COLS.length + '" class="pt-empty">No organisms match these filters.</td></tr>';
      } else {
        el.tbody.innerHTML = slice
          .map(function (row) {
            return (
              '<tr data-name="' +
              esc(row.name) +
              '">' +
              TABLE_COLS.map(function (c) {
                return cell(row, c[0]);
              }).join("") +
              "</tr>"
            );
          })
          .join("");
      }

      Array.prototype.forEach.call(el.thead.children, function (th) {
        var key = th.dataset.key;
        th.setAttribute("aria-sort", key === s.col ? (s.dir === 1 ? "ascending" : "descending") : "none");
      });

      el.pager.innerHTML =
        '<label>Rows <select class="pt-per">' +
        [25, 50, 100, 250]
          .map(function (n) {
            return '<option value="' + n + '"' + (n === state.perPage ? " selected" : "") + ">" + n + "</option>";
          })
          .join("") +
        "</select></label>" +
        '<button class="pt-btn pt-prev"' +
        (state.page === 1 ? " disabled" : "") +
        ">Prev</button>" +
        "<span>Page " +
        state.page +
        " of " +
        pages +
        "</span>" +
        '<button class="pt-btn pt-next"' +
        (state.page >= pages ? " disabled" : "") +
        ">Next</button>";
    }

    function refresh(skipFacets) {
      currentRows = filter(null);
      if (!skipFacets) renderFacets();
      renderChips();
      renderTable();
    }

    /* ── detail drawer ────────────────────────────────────────────────── */

    function openDrawer(name) {
      var row = rows.find(function (r) {
        return r.name === name;
      });
      if (!row) return;

      var body = DETAIL_COLS.map(function (c) {
        var v = displayValue(row[c[0]]);
        if (!v) return "";
        if (c[0] === "taxid") {
          v =
            '<a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=' +
            esc(v) +
            '" target="_blank" rel="noopener">' +
            esc(v) +
            "</a>";
        } else if (c[0] === "assembly_accession") {
          v =
            '<a href="https://www.ncbi.nlm.nih.gov/datasets/genome/' +
            esc(v) +
            '/" target="_blank" rel="noopener">' +
            esc(v) +
            "</a>";
        } else {
          v = esc(v);
        }
        return "<dt>" + esc(c[1]) + "</dt><dd>" + v + "</dd>";
      }).join("");

      el.drawer.innerHTML =
        '<button class="pt-drawer-close" aria-label="Close">✕</button><h3>' +
        esc(row.name) +
        "</h3><dl class=\"pt-dl\">" +
        body +
        "</dl>";
      el.drawer.classList.add("is-open");
    }

    function closeDrawer() {
      el.drawer.classList.remove("is-open");
    }

    /* ── export ───────────────────────────────────────────────────────── */

    function exportCsv() {
      var header = cols.filter(function (c) {
        return c !== "__blob";
      });
      var lines = [header.join(",")];
      currentRows.forEach(function (row) {
        lines.push(
          header
            .map(function (c) {
              var v = displayValue(row[c]);
              return /[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v;
            })
            .join(",")
        );
      });
      var blob = new Blob([lines.join("\n")], { type: "text/csv;charset=utf-8" });
      var a = document.createElement("a");
      a.href = URL.createObjectURL(blob);
      a.download = "pathogen_sheet_filtered.csv";
      a.click();
      URL.revokeObjectURL(a.href);
    }

    /* ── events ───────────────────────────────────────────────────────── */

    var debounce;
    el.search.addEventListener("input", function () {
      clearTimeout(debounce);
      debounce = setTimeout(function () {
        state.q = el.search.value;
        state.page = 1;
        refresh();
      }, 140);
    });

    el.facets.addEventListener("change", function (e) {
      var t = e.target;
      if (t.type === "checkbox") {
        var key = t.closest(".pt-facet").dataset.key;
        var sel = state.facets[key];
        if (t.checked) sel.add(t.value);
        else sel.delete(t.value);
        state.page = 1;
        refresh();
      } else if (t.tagName === "SELECT" && t.dataset.rank) {
        var rank = t.dataset.rank;
        state.tax[rank] = t.value;
        // Clearing or changing a rank invalidates everything below it.
        var from = TAXONOMY.indexOf(rank);
        TAXONOMY.slice(from + 1).forEach(function (r) {
          state.tax[r] = "";
        });
        state.page = 1;
        refresh();
      }
    });

    // Client-side filtering of a long facet list (host organism).
    el.facets.addEventListener("input", function (e) {
      if (!e.target.classList.contains("pt-facet-search")) return;
      var q = e.target.value.toLowerCase();
      var body = e.target.parentElement;
      Array.prototype.forEach.call(body.querySelectorAll(".pt-opt"), function (opt) {
        opt.style.display = opt.dataset.v.indexOf(q) === -1 ? "none" : "";
      });
    });

    el.chips.addEventListener("click", function (e) {
      var chip = e.target.closest(".pt-chip");
      if (!chip) return;
      if (chip.dataset.type === "facet") {
        state.facets[chip.dataset.key].delete(chip.dataset.v);
      } else {
        var from = TAXONOMY.indexOf(chip.dataset.key);
        TAXONOMY.slice(from).forEach(function (r) {
          state.tax[r] = "";
        });
      }
      state.page = 1;
      refresh();
    });

    el.reset.addEventListener("click", function () {
      state.q = "";
      el.search.value = "";
      FACETS.forEach(function (f) {
        state.facets[f[0]].clear();
      });
      TAXONOMY.forEach(function (r) {
        state.tax[r] = "";
      });
      state.page = 1;
      refresh();
    });

    el.exportBtn.addEventListener("click", exportCsv);

    el.thead.addEventListener("click", function (e) {
      var th = e.target.closest("th");
      if (!th || !th.dataset.key) return;
      if (state.sort.col === th.dataset.key) state.sort.dir *= -1;
      else state.sort = { col: th.dataset.key, dir: 1 };
      renderTable();
    });

    el.tbody.addEventListener("click", function (e) {
      if (e.target.tagName === "A") return; // let NCBI links through
      var tr = e.target.closest("tr");
      if (tr && tr.dataset.name) openDrawer(tr.dataset.name);
    });

    el.drawer.addEventListener("click", function (e) {
      if (e.target.closest(".pt-drawer-close")) closeDrawer();
    });

    document.addEventListener("keydown", function (e) {
      if (e.key === "Escape") closeDrawer();
    });

    el.pager.addEventListener("click", function (e) {
      if (e.target.closest(".pt-prev")) {
        state.page--;
        renderTable();
      } else if (e.target.closest(".pt-next")) {
        state.page++;
        renderTable();
      }
    });

    el.pager.addEventListener("change", function (e) {
      if (e.target.classList.contains("pt-per")) {
        state.perPage = Number(e.target.value);
        state.page = 1;
        renderTable();
      }
    });

    refresh();
  }

  function template() {
    return (
      '<div class="pt-toolbar">' +
      '<div class="pt-search"><input type="search" placeholder="Search organism, synonym, tax ID, lineage…" aria-label="Search the pathogen sheet"></div>' +
      '<span class="pt-count"></span>' +
      '<button class="pt-btn pt-reset" hidden>Clear filters</button>' +
      '<button class="pt-btn pt-export">Export CSV</button>' +
      "</div>" +
      '<div class="pt-chips"></div>' +
      '<div class="pt-body">' +
      '<aside class="pt-facets"></aside>' +
      "<div>" +
      '<div class="pt-table-wrap"><table class="pt-table"><thead><tr>' +
      TABLE_COLS.map(function (c) {
        return '<th data-key="' + c[0] + '" aria-sort="none">' + c[1] + "</th>";
      }).join("") +
      "</tr></thead><tbody></tbody></table></div>" +
      '<div class="pt-pager"></div>' +
      "</div></div>" +
      '<div class="pt-drawer"></div>'
    );
  }

  function boot() {
    init(document.getElementById("pathogen-app"));
  }

  // Material's instant navigation swaps content without a page load, so hook
  // document$ when it exists and fall back to a plain DOM ready otherwise.
  if (typeof document$ !== "undefined") {
    document$.subscribe(boot);
  } else if (document.readyState !== "loading") {
    boot();
  } else {
    document.addEventListener("DOMContentLoaded", boot);
  }
})();
