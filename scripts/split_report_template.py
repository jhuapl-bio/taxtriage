#!/usr/bin/env python3
"""
One-time migration: turn the monolithic assets/heatmap.html into a thin shell
plus external, per-tab/plot source files under assets/src/.

After this runs:

    assets/heatmap.html        thin shell (~3.4k lines): <head>, the HTML body
                               markup, and <link>/<script src> references to the
                               external assets below. Directly openable in a
                               browser (with assets/heatmap_boot.js present).
    assets/src/css/*.css       raw stylesheets (geo, main, loading, compare)
    assets/src/js/early.js     the small head <script> (reads window.HEATMAP_BOOT)
    assets/src/js/NN_*.js      the main bundle, one file per author section
                               (heatmap, tass, sunburst, coverage, proteins,
                                table, histogram, explore, summary, map, …)

The two downstream builders inline these parts back into one self-contained file
via scripts/report_template.inline_template():

    bin/make_report.py          -> all.comparison.report.html  (nextflow)
    scripts/inline_boot_json.py -> _site/index.html            (GitHub Pages)

This script is a bootstrap; you normally won't run it again. Day to day, edit
assets/heatmap.html and the files under assets/src/ directly.

USAGE
    python scripts/split_report_template.py            # migrate
    python scripts/split_report_template.py --force    # re-run, overwriting
"""
import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = REPO_ROOT / "assets" / "heatmap.html"
SRC_DIR = REPO_ROOT / "assets" / "src"

# ── exact wrapper lines in the monolith (1-indexed) ─────────────────────────
# Verified against the source: <style>/<script> open & close lines.
CSS_BLOCKS = [
    # (open_line, close_line, out_name, style_id_or_None)
    (16,   351,  "geo.css",     "geo-host-extra-css"),
    (549,  4208, "main.css",    None),
    (4218, 4306, "loading.css", "tt-loading-css"),
    (7522, 7656, "compare.css", None),
]
EARLY_SCRIPT = (355, 548, "early.js")          # head <script> … </script>
MAIN_SCRIPT = (7758, 30539)                     # big <script> … </script> (wrappers)

# JS section boundaries inside the main script (start line -> out name). Each
# section runs to the line before the next; the last ends just before </script>.
# Boundaries are the author's own /* ═══ § SECTION ═══ */ banners.
JS_SECTIONS = [
    (7759,  "00_preamble.js"),
    (7803,  "01_global_state.js"),
    (7993,  "02_utilities.js"),
    (9759,  "03_sample_sidebar.js"),
    (10601, "04_helpers.js"),
    (10613, "05_tab_heatmap.js"),
    (10955, "06_tab_tass.js"),
    (11216, "07_tab_sunburst.js"),
    (11814, "08_tab_coverage.js"),
    (12044, "09_protein_col_remap.js"),
    (12100, "10_tab_proteins_vfamr.js"),
    (13551, "11_tab_table.js"),
    (13640, "12_tab_histogram.js"),
    (15065, "13_pinned_row_comparison.js"),
    (15262, "14_heatmap_value_selector.js"),
    (15273, "15_tab_switching.js"),
    (15438, "16_redraw_dispatcher.js"),
    (15449, "17_tab_explore.js"),
    (16781, "18_tab_summary.js"),
    (19543, "19_cross_sample_organism.js"),
    (20153, "20_feature_compare.js"),
    (20281, "21_compare_coverages.js"),
    (23240, "22_data_capability_detection.js"),
    (23361, "23_filter_listeners.js"),
    (23892, "24_init.js"),
    (23898, "25_resizable_detail_panel.js"),
    (23940, "26_resizable_right_sidebar.js"),
    (24034, "27_file_upload.js"),
    (24751, "28_meta_csv.js"),
    (25202, "29_session_state.js"),
    (25526, "30_tab_runmeta.js"),
    (25641, "31_longitudinal_first.js"),
    (26261, "32_longitudinal_second.js"),
    (27141, "33_categorical_metadata.js"),
    (28654, "34_cross_entry_comparison.js"),
    (29165, "35_tab_map.js"),
    (29990, "36_guided_help_tour.js"),
]

ANCHORS = [
    (1,     "<!doctype html>"),
    (16,    "geo-host-extra-css"),
    (353,   "heatmap_boot.js"),
    (355,   "<script>"),
    (549,   "<style>"),
    (7758,  "<script>"),
    (30539, "</script>"),
]


def read_lines(path):
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return fh.read().splitlines(keepends=True)


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing assets/src/")
    args = ap.parse_args()

    lines = read_lines(TEMPLATE)
    n = len(lines)
    for ln, needle in ANCHORS:
        if ln > n or needle not in lines[ln - 1]:
            sys.exit(f"ERROR: anchor failed at line {ln} (want {needle!r}). "
                     f"Template drifted; update the line tables in this script.")

    if SRC_DIR.exists() and not args.force:
        sys.exit(f"ERROR: {SRC_DIR} exists. Re-run with --force.")
    (SRC_DIR / "css").mkdir(parents=True, exist_ok=True)
    (SRC_DIR / "js").mkdir(parents=True, exist_ok=True)

    def slice_lines(a, b):  # inclusive, 1-indexed
        return "".join(lines[a - 1:b])

    # ── write CSS (inner content only) ──────────────────────────────────────
    css_link = {}
    for open_ln, close_ln, name, sid in CSS_BLOCKS:
        write(SRC_DIR / "css" / name, slice_lines(open_ln + 1, close_ln - 1))
        attr = f' data-style-id="{sid}"' if sid else ""
        css_link[open_ln] = f'    <link rel="stylesheet" href="src/css/{name}"{attr} />\n'
        print(f"  css/{name:<14} <- lines {open_ln + 1}-{close_ln - 1}")

    # ── write early head script (inner only) ────────────────────────────────
    eo, ec, ename = EARLY_SCRIPT
    write(SRC_DIR / "js" / ename, slice_lines(eo + 1, ec - 1))
    print(f"  js/{ename:<22} <- lines {eo + 1}-{ec - 1}")

    # ── write main-script sections (inner only) ─────────────────────────────
    mo, mc = MAIN_SCRIPT
    js_tags = []
    for i, (start, name) in enumerate(JS_SECTIONS):
        end = JS_SECTIONS[i + 1][0] - 1 if i + 1 < len(JS_SECTIONS) else mc - 1
        write(SRC_DIR / "js" / name, slice_lines(start, end))
        js_tags.append(f'    <script src="src/js/{name}"></script>\n')
        print(f"  js/{name:<28} <- lines {start}-{end}")

    # ── assemble the thin shell, preserving exact document order ────────────
    out = []
    out.append(slice_lines(1, 15))                  # <head> + CDN links
    out.append(css_link[16])                        # geo stylesheet
    out.append(slice_lines(352, 354))               # pages.js comment + boot anchor
    out.append(f'    <script src="src/js/{ename}"></script>\n')
    out.append(css_link[549])                        # main stylesheet
    out.append(slice_lines(4209, 4217))             # comment between main & loading
    out.append(css_link[4218])                       # loading stylesheet
    out.append(slice_lines(4307, 7521))             # banner, tab bar, panels
    out.append(css_link[7522])                       # compare stylesheet
    out.append(slice_lines(7657, 7757))             # overlays
    out.extend(js_tags)                              # 37 main-script sections
    out.append(slice_lines(30540, n))               # tail (nov-cell-tip, </body></html>)
    write(TEMPLATE, "".join(out))

    print(f"\nOK: thin {TEMPLATE.relative_to(REPO_ROOT)} written "
          f"({sum(s.count(chr(10)) for s in out)} lines).")
    print("Next: scripts/report_template.inline_template() reassembles the "
          "self-contained file for the downstream builders.")


if __name__ == "__main__":
    main()
