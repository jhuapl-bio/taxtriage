/* Interactive browser for assets/pathogen_sheet.csv.
 *
 * The CSV is fetched once from the pipeline repository on GitHub at page load
 * (see REF below) and everything after that is client side: full-text search,
 * multi-select facets, a dependent taxonomy drilldown, sorting, pagination and
 * CSV export of the current selection.
 *
 * Facet counts are computed "excluding self" — the counts shown next to the
 * options of facet X reflect every active filter *except* X. That is what makes
 * multi-select feel right: options in the group you are working in don't
 * vanish as you tick them.
 */

(function () {
  "use strict";

  // The sheet is read from the repository at page load rather than bundled into
  // the site, so it is never a stale copy. The ref is pinned per docs version
  // by scripts/write_docs_ref.py (a release tag for versioned docs, `main` for
  // `latest`), so the 3.1 docs show the 3.1 sheet.
  // raw.githubusercontent.com sends Access-Control-Allow-Origin: *, so this is
  // a plain cross-origin GET with no proxy needed.
  var REF = (window.TAXTRIAGE_DOCS && window.TAXTRIAGE_DOCS.ref) || "main";
  var CSV_URL =
    "https://raw.githubusercontent.com/jhuapl-bio/taxtriage/" + encodeURIComponent(REF) + "/assets/pathogen_sheet.csv";
  var CSV_HUMAN_URL =
    "https://github.com/jhuapl-bio/taxtriage/blob/" + encodeURIComponent(REF) + "/assets/pathogen_sheet.csv";

  // Requesting an organism only makes sense against the current sheet, so the
  // call to action is hidden in frozen release docs. Anything that is not a
  // numbered X.Y.Z version — `latest`, PR previews, local `mkdocs serve` — is
  // treated as current and shows it.
  var DOCS_VERSION = (window.TAXTRIAGE_DOCS && window.TAXTRIAGE_DOCS.version) || "latest";
  var IS_CURRENT = !/^\d+\.\d+(\.\d+)?/.test(DOCS_VERSION);

  var ISSUE_TEMPLATE = "add_organism.yml";
  var UPDATE_TEMPLATE = "update_organism.yml";
  var ISSUES_NEW = "https://github.com/jhuapl-bio/taxtriage/issues/new";

  // GitHub issue forms are static — they cannot repeat a field group, so
  // requesting several organisms at once has to be assembled here and handed
  // over as a single prefilled issue body. The guided single-organism form
  // remains for people who start from the Issues tab instead.
  function templateUrl(name, kind) {
    var upd = kind === "update";
    var params = ["template=" + encodeURIComponent(upd ? UPDATE_TEMPLATE : ISSUE_TEMPLATE)];
    if (name) {
      params.push("title=" + encodeURIComponent((upd ? "Update organism: " : "Add organism: ") + name));
      params.push("organism=" + encodeURIComponent(name));
    }
    return ISSUES_NEW + "?" + params.join("&");
  }

  // Canonical body sites. Constrained rather than free text because the sheet
  // already carries both "abscess" and "abcess" from hand entry.
  var SITES = [
    "blood",
    "resp",
    "lung",
    "sputum",
    "nasal",
    "sinus",
    "oral",
    "teeth",
    "gut",
    "stool",
    "feces",
    "urine",
    "vaginal",
    "skin",
    "abscess",
    "ear",
    "eye",
    "ocular",
    "cornea",
    "csf",
    "bile",
    "mucus",
  ];

  // One spec drives the form, the issue body and the CSV block, so they cannot
  // drift apart. `id` matches the pathogen sheet column name.
  var REQUEST_FIELDS = [
    { id: "name", label: "Organism name", type: "text", required: true, ph: "Streptococcus pneumoniae" },
    { id: "taxid", label: "NCBI tax ID", type: "text", required: true, ph: "1313" },
    {
      id: "general_classification",
      label: "Microbial category",
      type: "select",
      required: true,
      options: ["primary", "opportunistic", "potential", "commensal"],
    },
    { id: "status", label: "Evidence status", type: "select", options: ["established", "putative"] },
    { id: "pathogenic_sites", label: "Pathogenic sites", type: "sites" },
    { id: "commensal_sites", label: "Commensal sites", type: "sites" },
    { id: "high_consequence", label: "High consequence", type: "select", options: ["FALSE", "TRUE"] },
    { id: "mol_type", label: "Molecule type", type: "select", options: ["dna", "rna"] },
    { id: "host_organism", label: "Host organism", type: "text", ph: "human" },
    { id: "assembly_accession", label: "Assembly accession", type: "text", ph: "GCF_000006885.1" },
    { id: "alternative_names", label: "Synonyms", type: "text" },
    { id: "pathology", label: "Disease association", type: "text" },
    { id: "reference", label: "Supporting reference", type: "textarea", required: true },
    { id: "request_type", label: "Request type (set automatically)", type: "request_type" },
  ];

  // Provenance of each row, standardised on the submission route rather than
  // chosen: Open issue -> git-tracked, Download -> external-local. APL-derived
  // is shown so the full vocabulary is visible, but is reserved for entries the
  // maintainers curate. The whole control is display-only; the value is stamped
  // at submit time, so it always reflects how the request actually arrived.
  var REQUEST_TYPES = [
    { value: "APL-derived", label: "APL-derived (maintainers only)", disabled: true },
    { value: "git-tracked", label: "git-tracked — filed as a GitHub issue" },
    { value: "external-local", label: "external-local — downloaded and sent privately" },
  ];
  var TYPE_FOR_ISSUE = "git-tracked";
  var TYPE_FOR_FILE = "external-local";

  // GitHub truncates very long prefilled bodies, and browsers cap URL length.
  // Warn before we get near it rather than silently losing entries.
  var MAX_URL = 6000;

  function downloadFile(filename, text, mime) {
    var blob = new Blob([text], { type: mime + ";charset=utf-8" });
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    a.click();
    URL.revokeObjectURL(a.href);
  }

  function csvCell(v) {
    v = v == null ? "" : String(v);
    return /[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v;
  }

  function isUpdate(e) {
    return e.__mode === "update";
  }

  // Display-ready snapshot of a sheet row, keyed by column, so the form and the
  // diff both work against plain strings rather than the parsed row.
  function snapshot(row) {
    var o = {};
    ALL_COLS.forEach(function (c) {
      o[c] = displayValue(row[c]);
    });
    return o;
  }

  // Order-insensitive for the comma-separated columns: re-selecting the same
  // sites in a different order is not a change.
  function normValue(id, v) {
    v = (v == null ? "" : String(v)).trim();
    if (LIST_COLS.indexOf(id) === -1) return v;
    return v
      .split(",")
      .map(function (t) {
        return t.trim();
      })
      .filter(Boolean)
      .sort()
      .join(", ");
  }

  // The whole point of an update request: only the fields that actually moved.
  function diffEntry(entry) {
    var base = entry.__base || {};
    var out = [];
    REQUEST_FIELDS.forEach(function (f) {
      if (f.id === "request_type") return;
      var before = base[f.id] || "";
      var after = entry[f.id] || "";
      if (normValue(f.id, before) !== normValue(f.id, after)) {
        out.push({ id: f.id, label: f.label, before: before, after: after });
      }
    });
    return out;
  }

  // The proposed full row: the current record with the edited fields applied,
  // so lineage and anything the form does not expose survives untouched.
  function mergedRow(entry) {
    var out = {};
    ALL_COLS.forEach(function (c) {
      out[c] = (entry.__base && entry.__base[c]) || "";
    });
    REQUEST_FIELDS.forEach(function (f) {
      if (f.id === "request_type") return; // provenance of the existing row stands
      if (entry[f.id] !== undefined) out[f.id] = entry[f.id];
    });
    return out;
  }

  function mdCell(v) {
    v = (v == null ? "" : String(v)).replace(/\|/g, "\\|").replace(/\r?\n/g, " ");
    return v === "" ? "_(empty)_" : v;
  }

  function csvBlock(entries) {
    var lines = ["```csv", ALL_COLS.join(",")];
    entries.forEach(function (e) {
      lines.push(
        ALL_COLS.map(function (c) {
          return csvCell(e[c] || "");
        }).join(","),
      );
    });
    lines.push("```", "");
    return lines;
  }

  function requestBody(entries) {
    var adds = entries.filter(function (e) {
      return !isUpdate(e);
    });
    var ups = entries.filter(isUpdate);

    var lines = [
      "Assembled from the [Pathogen Sheet](https://jhuapl-bio.github.io/taxtriage/latest/pathogen-sheet/).",
      "Taxonomic lineage is omitted — it is derived from the tax ID.",
      "",
    ];

    if (adds.length) {
      lines.push("### Organisms requested (" + adds.length + " new)", "");
      adds.forEach(function (e, i) {
        lines.push("#### " + (i + 1) + ". " + (e.name || "(unnamed)"), "");
        REQUEST_FIELDS.forEach(function (f) {
          if (f.id === "name") return;
          var v = e[f.id];
          if (!v) return;
          lines.push("- **" + f.label + ":** " + v);
        });
        lines.push("");
      });
      // A paste-ready block in pathogen_sheet.csv column order.
      lines.push("#### Draft CSV rows", "");
      lines = lines.concat(csvBlock(adds));
    }

    if (ups.length) {
      lines.push("### Updates requested (" + ups.length + ")", "");
      lines.push(
        "Only the fields listed below should change; everything else stays as it is, " +
          "including `request_type`, which records how the existing row first arrived.",
        "",
      );
      ups.forEach(function (e, i) {
        var base = e.__base || {};
        lines.push(
          "#### " + (i + 1) + ". " + (e.name || base.name || "(unnamed)") + " (taxid " + (base.taxid || "—") + ")",
          "",
        );
        if (e.__reason) lines.push("**Reason:** " + e.__reason, "");
        lines.push("| Field | Current | Proposed |", "| --- | --- | --- |");
        diffEntry(e).forEach(function (c) {
          lines.push("| `" + c.id + "` | " + mdCell(c.before) + " | " + mdCell(c.after) + " |");
        });
        lines.push("");
        lines.push("<details><summary>Full row — current, then proposed</summary>", "");
        lines = lines.concat(csvBlock([base, mergedRow(e)]));
        lines.push("</details>", "");
      });
    }

    return lines.join("\n");
  }

  function requestTitle(entries) {
    var adds = entries.filter(function (e) {
      return !isUpdate(e);
    });
    var ups = entries.filter(isUpdate);
    function plural(n, w) {
      return n + " " + w + (n === 1 ? "" : "s");
    }
    if (ups.length && !adds.length) {
      return ups.length === 1
        ? "Update organism: " + (ups[0].name || (ups[0].__base || {}).name || "")
        : "Update " + ups.length + " organisms in the pathogen sheet";
    }
    if (adds.length && !ups.length) {
      return adds.length === 1
        ? "Add organism: " + adds[0].name
        : "Add " + adds.length + " organisms to the pathogen sheet";
    }
    return "Pathogen sheet: " + plural(adds.length, "addition") + ", " + plural(ups.length, "update");
  }

  function requestIssueUrl(entries) {
    return (
      ISSUES_NEW +
      "?labels=pathogen-sheet" +
      "&title=" +
      encodeURIComponent(requestTitle(entries)) +
      "&body=" +
      encodeURIComponent(requestBody(entries))
    );
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
    ["request_type", "Request type", false, false],
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
    ["assembly_accession", "Assembly"],
    ["request_type", "Request type"],
    ["reference", "Reference"],
  ];

  // Canonical column set. Anything the CSV omits reads as empty.
  var ALL_COLS = [
    "name",
    "taxid",
    "general_classification",
    "alternative_names",
    "pathogenic_sites",
    "commensal_sites",
    "status",
    "high_consequence",
    "pathology",
    "host_organism",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "mol_type",
    "reference",
    "Additional references",
    "assembly_accession",
    "request_type",
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
    ["request_type", "Request type"],
    ["reference", "Reference"],
    ["Additional references", "Additional references"],
  ];

  // Columns split on commas into multiple values.
  var LIST_COLS = ["pathogenic_sites", "commensal_sites", "alternative_names"];

  // Known upstream typos, folded so the filter list stays clean.
  var VALUE_FIXES = {
    status: { estbalished: "established", etsablished: "established" },
  };

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

  /* ── CSV ────────────────────────────────────────────────────────────────
   * RFC 4180 parser. The pathogen sheet's `reference` column contains commas,
   * doubled quotes and occasional newlines inside quoted fields, so splitting
   * on commas is not an option.
   */
  function parseCsv(text) {
    // Strip BOM and normalise line endings before scanning.
    text = text.replace(/^﻿/, "").replace(/\r\n?/g, "\n");

    var rows = [];
    var row = [];
    var field = "";
    var inQuotes = false;

    for (var i = 0; i < text.length; i++) {
      var c = text[i];

      if (inQuotes) {
        if (c === '"') {
          if (text[i + 1] === '"') {
            field += '"';
            i++;
          } else {
            inQuotes = false;
          }
        } else {
          field += c;
        }
        continue;
      }

      if (c === '"') {
        inQuotes = true;
      } else if (c === ",") {
        row.push(field);
        field = "";
      } else if (c === "\n") {
        row.push(field);
        field = "";
        // Skip blank lines rather than emitting empty records.
        if (row.length > 1 || row[0] !== "") rows.push(row);
        row = [];
      } else {
        field += c;
      }
    }
    // Trailing field/record with no closing newline.
    if (field !== "" || row.length) {
      row.push(field);
      if (row.length > 1 || row[0] !== "") rows.push(row);
    }

    if (!rows.length) return { cols: [], rows: [] };
    var header = rows[0].map(function (h) {
      return h.trim();
    });
    return { cols: header, rows: rows.slice(1) };
  }

  /* Turn parsed CSV into row objects.
   *
   * Deliberately tolerant of schema drift: any column the sheet does not carry
   * simply reads as empty rather than breaking the page, and facets with
   * nothing in them are hidden by the renderer. This matters because the CSV is
   * read live from a branch or tag that can change independently of the site.
   */
  function normalise(parsed) {
    var index = {};
    parsed.cols.forEach(function (c, i) {
      index[c] = i;
    });

    var present = {};
    ALL_COLS.forEach(function (c) {
      present[c] = Object.prototype.hasOwnProperty.call(index, c);
    });

    var rows = parsed.rows
      .map(function (cells) {
        var o = {};
        ALL_COLS.forEach(function (c) {
          var raw = present[c] ? (cells[index[c]] || "").trim() : "";
          var fix = VALUE_FIXES[c];
          if (fix && Object.prototype.hasOwnProperty.call(fix, raw)) raw = fix[raw];

          if (LIST_COLS.indexOf(c) !== -1) {
            o[c] = raw
              ? raw
                  .split(",")
                  .map(function (s) {
                    return s.trim();
                  })
                  .filter(Boolean)
              : [];
          } else if (c === "high_consequence") {
            o[c] = raw.toUpperCase() === "TRUE";
          } else {
            o[c] = raw;
          }
        });
        o.__blob = SEARCH_COLS.map(function (c) {
          return displayValue(o[c]);
        })
          .join(" ")
          .toLowerCase();
        return o;
      })
      .filter(function (o) {
        return o.name; // drop stray blank records
      });

    return { rows: rows, present: present };
  }

  function init(root) {
    if (!root || root.dataset.ptReady === "1") return;
    root.dataset.ptReady = "1";

    fetch(CSV_URL, { cache: "no-cache" })
      .then(function (r) {
        if (!r.ok) throw new Error("HTTP " + r.status);
        return r.text();
      })
      .then(function (text) {
        var data = normalise(parseCsv(text));
        if (!data.rows.length) throw new Error("the sheet parsed to zero rows");
        build(root, data);
      })
      .catch(function (err) {
        root.innerHTML =
          '<div class="pt-empty"><p>Could not load the pathogen sheet from the pipeline repository (' +
          esc(err.message) +
          ").</p><p>You can open it directly on GitHub: " +
          '<a href="' +
          CSV_HUMAN_URL +
          '" target="_blank" rel="noopener">assets/pathogen_sheet.csv</a></p></div>';
      });
  }

  function build(root, data) {
    var rows = data.rows;
    var present = data.present;

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
      request: root.querySelector(".pt-request-drawer"),
      confirm: root.querySelector(".pt-confirm"),
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
        // A column the sheet does not carry (or one that is entirely blank)
        // gets no facet at all rather than an empty expander.
        if (!entries.length && !sel.size) return;
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
      html += '<details class="pt-facet pt-tax" open><summary>Taxonomy</summary><div class="pt-facet-body">';
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
              "</button>",
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
              "</button>",
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
      if (key === "assembly_accession") {
        return v
          ? '<td class="pt-acc"><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/' +
              encodeURIComponent(v) +
              '/" target="_blank" rel="noopener">' +
              highlight(v) +
              "</a></td>"
          : "<td></td>";
      }
      if (key === "reference") {
        // References run to full citations; clamp in the cell and keep the
        // whole thing on hover and in the detail drawer.
        return '<td class="pt-ref" title="' + esc(v) + '">' + highlight(displayValue(v)) + "</td>";
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
        // A search that finds nothing is the most likely moment someone wants
        // an organism added, so offer it here with the query prefilled.
        var q = state.q.trim();
        var cta =
          IS_CURRENT && q
            ? '<p style="margin-top:0.75rem"><button class="pt-btn pt-request" data-name="' +
              esc(q) +
              '">Request &ldquo;' +
              esc(q) +
              "&rdquo;</button></p>"
            : "";
        el.tbody.innerHTML =
          '<tr><td colspan="' +
          TABLE_COLS.length +
          '" class="pt-empty"><p>No organisms match these filters.</p>' +
          cta +
          "</td></tr>";
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
        '</h3><dl class="pt-dl">' +
        body +
        "</dl>" +
        (IS_CURRENT
          ? '<p class="pt-drawer-cta"><button class="pt-btn pt-request" data-mode="update" data-name="' +
            esc(row.name) +
            '">Request an update to this entry</button></p>'
          : "");
      el.drawer.classList.add("is-open");
    }

    function closeDrawer() {
      el.drawer.classList.remove("is-open");
    }

    /* ── request builder ──────────────────────────────────────────────────
     * Stages one or more organisms, then hands them to GitHub as a single
     * prefilled issue. Lives here rather than in the issue form because issue
     * forms cannot repeat a field group.
     */

    var staged = [];
    var mode = "add"; // "add" | "update"
    var baseRow = null; // snapshot of the organism being updated
    var pending = null; // { list, wantsFile } awaiting confirmation

    var organismOptions = rows
      .map(function (r) {
        return '<option value="' + esc(r.name) + '"></option>';
      })
      .join("");

    var byName = new Map();
    rows.forEach(function (r) {
      byName.set(r.name.toLowerCase(), r);
    });

    function fieldControl(f, current) {
      current = current == null ? "" : String(current);
      // Required markers only bind in add mode — an update supplies whichever
      // fields it changes, and leaves the rest alone.
      var req = f.required && mode === "add" ? " required" : "";
      if (f.type === "select") {
        var opts = f.options.slice();
        if (current && opts.indexOf(current) === -1) opts.push(current);
        return (
          '<select data-f="' +
          f.id +
          '"' +
          req +
          '><option value="">—</option>' +
          opts
            .map(function (o) {
              return '<option value="' + esc(o) + '"' + (o === current ? " selected" : "") + ">" + esc(o) + "</option>";
            })
            .join("") +
          "</select>"
        );
      }
      if (f.type === "request_type") {
        // Display-only: the value is standardised on the submission route, so
        // an editable control would imply a choice that is then overwritten.
        return (
          '<select data-f="' +
          f.id +
          '" disabled aria-describedby="pt-rq-type-note">' +
          REQUEST_TYPES.map(function (o) {
            return (
              '<option value="' +
              esc(o.value) +
              '"' +
              (o.disabled ? " disabled" : "") +
              (o.value === TYPE_FOR_ISSUE ? " selected" : "") +
              ">" +
              esc(o.label) +
              "</option>"
            );
          }).join("") +
          "</select>"
        );
      }
      if (f.type === "sites") {
        // Values already in the sheet are offered even when they are not in the
        // canonical list, so loading a record never silently drops a site.
        var chosen = current
          .split(",")
          .map(function (t) {
            return t.trim();
          })
          .filter(Boolean);
        var all = SITES.slice();
        chosen.forEach(function (c) {
          if (all.indexOf(c) === -1) all.push(c);
        });
        return (
          '<select data-f="' +
          f.id +
          '" multiple size="4">' +
          all
            .map(function (o) {
              return (
                '<option value="' +
                esc(o) +
                '"' +
                (chosen.indexOf(o) !== -1 ? " selected" : "") +
                ">" +
                esc(o) +
                "</option>"
              );
            })
            .join("") +
          "</select>"
        );
      }
      if (f.type === "textarea") {
        return '<textarea data-f="' + f.id + '" rows="2"' + req + ">" + esc(current) + "</textarea>";
      }
      return (
        '<input type="text" data-f="' +
        f.id +
        '" value="' +
        esc(current) +
        '" placeholder="' +
        esc(f.ph || "") +
        '"' +
        req +
        ">"
      );
    }

    function renderRequest(values) {
      values = values || {};
      var isUpd = mode === "update";
      var loaded = isUpd && baseRow;

      var modes =
        '<div class="pt-rq-modes" role="group" aria-label="Kind of request">' +
        '<button type="button" class="pt-rq-mode' +
        (isUpd ? "" : " is-on") +
        '" data-mode="add" aria-pressed="' +
        (isUpd ? "false" : "true") +
        '">Add new organisms</button>' +
        '<button type="button" class="pt-rq-mode' +
        (isUpd ? " is-on" : "") +
        '" data-mode="update" aria-pressed="' +
        (isUpd ? "true" : "false") +
        '">Update existing entries</button>' +
        "</div>";

      var intro = isUpd
        ? '<p class="pt-rq-intro">Pick an organism already in the sheet, correct whatever is wrong, and only the fields you actually changed are requested. Queue as many as you like, then review everything before it leaves the page.</p>'
        : '<p class="pt-rq-intro">Queue as many as you like, then review before opening a single issue. Fields marked <em>*</em> are required; the rest help but can be left blank.</p>';

      var picker = isUpd
        ? '<label class="pt-rq-field pt-rq-wide pt-rq-target-field"><span>Organism to update <em aria-hidden="true">*</em></span>' +
          '<input type="text" class="pt-rq-target" list="pt-rq-organisms" autocomplete="off" placeholder="Start typing an organism already in the sheet…" value="' +
          esc(values.__target || "") +
          '">' +
          '<datalist id="pt-rq-organisms">' +
          organismOptions +
          "</datalist></label>"
        : "";

      var fields =
        isUpd && !loaded
          ? '<p class="pt-rq-empty">Choose an organism above to load its current values.</p>'
          : REQUEST_FIELDS.filter(function (f) {
              // request_type records how a row first arrived; an update to an
              // existing row does not change that, so it is not offered here.
              return !(isUpd && f.id === "request_type");
            })
              .map(function (f) {
                return (
                  '<label class="pt-rq-field' +
                  (f.type === "textarea" || f.type === "sites" ? " pt-rq-wide" : "") +
                  '"><span>' +
                  esc(f.label) +
                  (f.required && !isUpd ? ' <em aria-hidden="true">*</em>' : "") +
                  "</span>" +
                  fieldControl(f, values[f.id]) +
                  "</label>"
                );
              })
              .join("");

      var reason = loaded
        ? '<label class="pt-rq-field pt-rq-wide"><span>Reason for the change <em aria-hidden="true">*</em></span>' +
          '<textarea data-f="__reason" rows="2" placeholder="What is wrong today, and what should it say instead?">' +
          esc(values.__reason || "") +
          "</textarea></label>"
        : "";

      var live = loaded ? '<p class="pt-rq-live" aria-live="polite"></p>' : "";

      var list = staged.length
        ? '<ol class="pt-rq-list">' +
          staged
            .map(function (e, i) {
              var isU = isUpdate(e);
              var n = isU ? diffEntry(e).length : 0;
              var meta = isU
                ? n + " field" + (n === 1 ? "" : "s") + " changed"
                : "taxid " + esc(e.taxid) + ", " + esc(e.general_classification);
              return (
                '<li><span><span class="pt-rq-badge pt-rq-badge-' +
                (isU ? "update" : "add") +
                '">' +
                (isU ? "update" : "add") +
                "</span> <strong>" +
                esc(e.name || (e.__base || {}).name || "") +
                "</strong> — " +
                meta +
                '</span><button type="button" class="pt-rq-del" data-i="' +
                i +
                '" aria-label="Remove ' +
                esc(e.name || "") +
                '">✕</button></li>'
              );
            })
            .join("") +
          "</ol>"
        : '<p class="pt-rq-empty">Nothing staged yet. Fill the fields above and use <strong>' +
          (isUpd ? "Stage this update" : "Add another") +
          "</strong> to queue more than one.</p>";

      el.request.innerHTML =
        '<button class="pt-drawer-close" aria-label="Close">✕</button>' +
        "<h3>Request changes</h3>" +
        modes +
        intro +
        '<div class="pt-rq-form">' +
        picker +
        fields +
        reason +
        "</div>" +
        live +
        '<p class="pt-rq-error" role="alert" hidden></p>' +
        '<div class="pt-rq-actions">' +
        '<button type="button" class="pt-btn pt-rq-add">' +
        (isUpd ? "Stage this update" : "Add another") +
        "</button>" +
        '<button type="button" class="pt-btn pt-rq-open">Review &amp; open issue' +
        (staged.length ? " (" + staged.length + ")" : "") +
        "</button>" +
        '<button type="button" class="pt-btn pt-rq-download">Review &amp; download</button>' +
        "</div>" +
        '<p class="pt-rq-note" id="pt-rq-type-note">Either route shows you everything that will be submitted before it leaves this page. ' +
        "A public GitHub issue is the fastest; if the request is sensitive, <strong>Download</strong> saves the same content as a " +
        "file you can email privately instead.<br>" +
        (isUpd
          ? "<code>request_type</code> is left as it stands on the existing row — it records how that entry first arrived, which an update does not change."
          : "<code>request_type</code> follows whichever you choose — <strong>Open issue</strong> records <code>" +
            TYPE_FOR_ISSUE +
            "</code>, <strong>Download</strong> records <code>" +
            TYPE_FOR_FILE +
            "</code>.") +
        "</p>" +
        "<h4>Staged</h4>" +
        list +
        '<p class="pt-rq-alt">Prefer GitHub\'s guided form for a single organism? <a href="' +
        templateUrl("", isUpd ? "update" : "add") +
        '" target="_blank" rel="noopener">Use the ' +
        (isUpd ? "update" : "add") +
        " issue template</a>.</p>";

      updateLive();
    }

    function readForm() {
      var entry = {};
      REQUEST_FIELDS.forEach(function (f) {
        var node = el.request.querySelector('[data-f="' + f.id + '"]');
        if (!node) return;
        if (f.type === "sites") {
          entry[f.id] = Array.prototype.filter
            .call(node.options, function (o) {
              return o.selected;
            })
            .map(function (o) {
              return o.value;
            })
            .join(", ");
        } else {
          entry[f.id] = (node.value || "").trim();
        }
      });
      var target = el.request.querySelector(".pt-rq-target");
      if (target) entry.__target = target.value.trim();
      var reason = el.request.querySelector('[data-f="__reason"]');
      if (reason) entry.__reason = reason.value.trim();
      if (mode === "update") {
        entry.__mode = "update";
        entry.__base = baseRow;
      } else {
        entry.__mode = "add";
      }
      return entry;
    }

    function showError(msg) {
      var p = el.request.querySelector(".pt-rq-error");
      if (!p) return;
      p.textContent = msg;
      p.hidden = !msg;
    }

    // Running count of what an update would actually change, so the diff is
    // never a surprise at the confirmation step.
    function updateLive() {
      var node = el.request.querySelector(".pt-rq-live");
      if (!node) return;
      var n = diffEntry(readForm()).length;
      node.textContent = n
        ? n + " field" + (n === 1 ? "" : "s") + " changed against the current record."
        : "No changes yet — edit a field to request an update.";
      node.classList.toggle("is-on", !!n);
    }

    // Loads an organism into the form. Unknown names clear the loaded record
    // rather than half-filling the form from a previous one.
    function pickTarget(name) {
      var row = byName.get(
        String(name || "")
          .trim()
          .toLowerCase(),
      );
      if (!row) {
        if (!baseRow) return;
        baseRow = null;
        renderRequest({ __target: name });
        return;
      }
      if (baseRow && baseRow.name === row.name) return;
      baseRow = snapshot(row);
      var values = {};
      REQUEST_FIELDS.forEach(function (f) {
        values[f.id] = baseRow[f.id] || "";
      });
      values.__target = row.name;
      renderRequest(values);
      showError("");
    }

    // Returns the entry, or null after reporting what is missing.
    function takeEntry() {
      var entry = readForm();
      if (mode === "update") {
        if (!baseRow) {
          showError("Pick an organism to update first.");
          return null;
        }
        entry.name = entry.name || baseRow.name;
        if (!diffEntry(entry).length) {
          showError("Change at least one field before staging this update.");
          return null;
        }
        if (!entry.__reason) {
          showError("Add a short reason so a maintainer can review the change.");
          return null;
        }
        if (entry.taxid && !/^\d+$/.test(entry.taxid)) {
          showError("Tax ID should be numeric, e.g. 1313.");
          return null;
        }
        showError("");
        return entry;
      }

      var missing = REQUEST_FIELDS.filter(function (f) {
        return f.required && !entry[f.id];
      });
      if (missing.length) {
        showError(
          "Fill in " +
            missing
              .map(function (f) {
                return f.label.toLowerCase();
              })
              .join(", ") +
            " first.",
        );
        return null;
      }
      if (!/^\d+$/.test(entry.taxid)) {
        showError("Tax ID should be numeric, e.g. 1313.");
        return null;
      }
      showError("");
      return entry;
    }

    // Staged entries plus whatever is currently typed, so a single organism does
    // not require pressing "Add another" first. Returns null after reporting why.
    function collectEntries() {
      var list = staged.slice();
      var current = readForm();
      var started =
        mode === "update"
          ? !!baseRow && diffEntry(current).length > 0
          : REQUEST_FIELDS.some(function (f) {
              return f.required && current[f.id];
            });
      if (started) {
        var last = takeEntry();
        if (!last) return null;
        list.push(last);
      }
      if (!list.length) {
        showError(mode === "update" ? "Stage at least one update first." : "Add at least one organism first.");
        return null;
      }
      return list;
    }

    function openRequest(prefillName, wantMode) {
      mode = wantMode === "update" ? "update" : "add";
      baseRow = null;
      renderRequest();
      if (prefillName) {
        if (mode === "update") {
          var t = el.request.querySelector(".pt-rq-target");
          if (t) t.value = prefillName;
          pickTarget(prefillName);
        } else {
          var n = el.request.querySelector('[data-f="name"]');
          if (n) n.value = prefillName;
        }
      }
      el.request.classList.add("is-open");
      var first = el.request.querySelector("input,select,textarea");
      if (first && first.focus) first.focus();
    }

    /* ── confirmation ─────────────────────────────────────────────────────
     * Nothing leaves the page — no issue tab, no download — until the exact
     * content has been shown: full records for additions, and a current vs
     * proposed table of only the changed fields for updates.
     */

    function confirmMarkup(list, wantsFile) {
      var adds = list.filter(function (e) {
        return !isUpdate(e);
      });
      var ups = list.filter(isUpdate);
      var out = [];

      out.push('<div class="pt-cf-card" role="dialog" aria-modal="true" aria-labelledby="pt-cf-title">');
      out.push('<button class="pt-drawer-close pt-cf-cancel" aria-label="Back to editing">✕</button>');
      out.push('<h3 id="pt-cf-title">Review before submitting</h3>');
      out.push(
        '<p class="pt-cf-sum">' +
          (adds.length ? "<strong>" + adds.length + "</strong> new organism" + (adds.length === 1 ? "" : "s") : "") +
          (adds.length && ups.length ? " and " : "") +
          (ups.length ? "<strong>" + ups.length + "</strong> update" + (ups.length === 1 ? "" : "s") : "") +
          " — " +
          (wantsFile ? "saved as a Markdown file you can send privately." : "opened as a prefilled GitHub issue.") +
          (adds.length
            ? " New entries are recorded as <code>" + (wantsFile ? TYPE_FOR_FILE : TYPE_FOR_ISSUE) + "</code>."
            : "") +
          (ups.length ? " Updates leave <code>request_type</code> on the existing row untouched." : "") +
          "</p>",
      );
      out.push('<div class="pt-cf-body">');

      adds.forEach(function (e) {
        out.push('<section class="pt-cf-entry">');
        out.push('<h4><span class="pt-rq-badge pt-rq-badge-add">add</span> ' + esc(e.name) + "</h4>");
        out.push('<dl class="pt-dl pt-cf-dl">');
        REQUEST_FIELDS.forEach(function (f) {
          if (f.id === "name" || !e[f.id]) return;
          out.push("<dt>" + esc(f.label) + "</dt><dd>" + esc(e[f.id]) + "</dd>");
        });
        out.push("</dl></section>");
      });

      ups.forEach(function (e) {
        var base = e.__base || {};
        out.push('<section class="pt-cf-entry">');
        out.push(
          '<h4><span class="pt-rq-badge pt-rq-badge-update">update</span> ' +
            esc(e.name || base.name || "") +
            ' <span class="pt-cf-taxid">taxid ' +
            esc(base.taxid || "—") +
            "</span></h4>",
        );
        if (e.__reason) out.push('<p class="pt-cf-reason">' + esc(e.__reason) + "</p>");
        out.push(
          '<table class="pt-cf-diff"><thead><tr><th>Field</th><th>Current</th><th>Proposed</th></tr></thead><tbody>',
        );
        diffEntry(e).forEach(function (c) {
          out.push(
            "<tr><td>" +
              esc(c.label) +
              '</td><td class="pt-cf-before">' +
              (c.before ? esc(c.before) : "<em>empty</em>") +
              '</td><td class="pt-cf-after">' +
              (c.after ? esc(c.after) : "<em>empty</em>") +
              "</td></tr>",
          );
        });
        out.push("</tbody></table>");
        out.push('<p class="pt-cf-note">Every other field stays exactly as it is today.</p>');
        out.push("</section>");
      });

      out.push("</div>");
      out.push(
        '<div class="pt-cf-actions">' +
          '<button type="button" class="pt-btn pt-cf-cancel">Back to editing</button>' +
          '<button type="button" class="pt-btn pt-cf-go">' +
          (wantsFile ? "Confirm &amp; download" : "Confirm &amp; open issue") +
          "</button></div>",
      );
      out.push("</div>");
      return out.join("");
    }

    function openConfirm(list, wantsFile) {
      pending = { list: list, wantsFile: wantsFile };
      el.confirm.innerHTML = confirmMarkup(list, wantsFile);
      el.confirm.classList.add("is-open");
      var go = el.confirm.querySelector(".pt-cf-go");
      if (go && go.focus) go.focus();
    }

    function closeConfirm() {
      el.confirm.classList.remove("is-open");
      el.confirm.innerHTML = "";
      pending = null;
    }

    function submitRequest(list, wantsFile) {
      if (wantsFile) {
        // Same content as the issue, as a file that can be sent privately.
        var onlyUpdates = list.every(isUpdate);
        downloadFile(
          onlyUpdates ? "taxtriage-organism-update.md" : "taxtriage-organism-request.md",
          requestBody(list),
          "text/markdown",
        );
        return;
      }
      window.open(requestIssueUrl(list), "_blank", "noopener");
    }

    el.confirm.addEventListener("click", function (e) {
      if (e.target.closest(".pt-cf-cancel") || e.target === el.confirm) {
        closeConfirm();
        return;
      }
      if (e.target.closest(".pt-cf-go") && pending) {
        var p = pending;
        closeConfirm();
        submitRequest(p.list, p.wantsFile);
      }
    });

    el.request.addEventListener("change", function (e) {
      if (e.target.classList.contains("pt-rq-target")) {
        pickTarget(e.target.value);
        return;
      }
      updateLive();
    });

    el.request.addEventListener("input", function (e) {
      // A datalist pick fires input, not change — load it as soon as the value
      // is an exact match so the form fills without waiting for a blur.
      if (e.target.classList.contains("pt-rq-target")) {
        if (byName.has(e.target.value.trim().toLowerCase())) pickTarget(e.target.value);
        return;
      }
      updateLive();
    });

    el.request.addEventListener("click", function (e) {
      if (e.target.closest(".pt-drawer-close")) {
        el.request.classList.remove("is-open");
        return;
      }

      var modeBtn = e.target.closest(".pt-rq-mode");
      if (modeBtn) {
        if (modeBtn.dataset.mode === mode) return;
        mode = modeBtn.dataset.mode;
        baseRow = null;
        renderRequest();
        showError("");
        return;
      }

      var del = e.target.closest(".pt-rq-del");
      if (del) {
        staged.splice(Number(del.dataset.i), 1);
        var keep = readForm();
        renderRequest(keep);
        return;
      }

      if (e.target.closest(".pt-rq-add")) {
        var entry = takeEntry();
        if (!entry) return;
        staged.push(entry);
        baseRow = null;
        renderRequest(); // clears the form, ready for the next one
        return;
      }

      var wantsIssue = !!e.target.closest(".pt-rq-open");
      var wantsFile = !!e.target.closest(".pt-rq-download");
      if (wantsIssue || wantsFile) {
        var list = collectEntries();
        if (!list) return;

        // The picker is a preview; the route is the truth. Stamping here means
        // request_type always records how the request actually arrived, and the
        // confirmation shows the value that will really be filed.
        list = list.map(function (x) {
          if (!isUpdate(x)) x.request_type = wantsFile ? TYPE_FOR_FILE : TYPE_FOR_ISSUE;
          return x;
        });

        if (wantsIssue && requestIssueUrl(list).length > MAX_URL) {
          showError(
            "That is too much for one issue (" +
              list.length +
              " entries). Use Download instead, or submit this batch and start another.",
          );
          return;
        }

        openConfirm(list, wantsFile);
      }
    });

    /* ── export ───────────────────────────────────────────────────────── */

    function exportCsv() {
      // Export the canonical schema so downstream consumers get a stable shape
      // even when the upstream sheet is missing a column.
      var header = ALL_COLS;
      var lines = [header.join(",")];
      currentRows.forEach(function (row) {
        lines.push(
          header
            .map(function (c) {
              var v = displayValue(row[c]);
              return /[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v;
            })
            .join(","),
        );
      });
      downloadFile("pathogen_sheet_filtered.csv", lines.join("\n"), "text/csv");
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

    root.addEventListener("click", function (e) {
      var btn = e.target.closest(".pt-request");
      if (!btn || el.request.contains(btn)) return;
      openRequest(btn.dataset.name || "", btn.dataset.mode || "add");
    });

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
      if (e.key !== "Escape") return;
      // The confirmation sits on top of everything, so it unwinds first —
      // Escape there means "back to editing", not "throw the request away".
      if (el.confirm.classList.contains("is-open")) {
        closeConfirm();
        return;
      }
      closeDrawer();
      el.request.classList.remove("is-open");
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
      (IS_CURRENT ? '<button class="pt-btn pt-request">Request changes</button>' : "") +
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
      '<div class="pt-drawer"></div>' +
      '<div class="pt-drawer pt-request-drawer"></div>' +
      '<div class="pt-confirm"></div>'
    );
  }

  function boot() {
    var mount = document.getElementById("pathogen-app");
    // Widen the whole page for this table only. Done as a body class rather
    // than page-scoped CSS so the styling lives in the stylesheet; with
    // Material's instant navigation the class must be removed again on leaving.
    document.body.classList.toggle("pt-fullwidth", !!mount);
    init(mount);
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
