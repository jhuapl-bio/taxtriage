#!/usr/bin/env python3
"""
Build _site/index.html for GitHub Pages deployment.

Strategy
--------
1. Start from assets/heatmap.html (the authoritative template).
2. Extract JSON from the data source (assets/pages.js preferred,
   assets/heatmap_boot.js as fallback).
3. Swap the local <script src="heatmap_boot.js"> loader for an inlined
   <script id="BOOTSTRAP" type="application/json"> block.
4. Swap the window.HEATMAP_BOOT reader for JSON.parse(…BOOTSTRAP…).
5. Write the result to docs/index.html.

Both JS files use the format:  window.HEATMAP_BOOT = { ... };
Requires Node.js on PATH (handles trailing-comma JS object literals).

Usage:
    python scripts/inline_boot_json.py          # build and write
    python scripts/inline_boot_json.py --check  # validate only, don't write
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

# ── paths ──────────────────────────────────────────────────────────────────
REPO_ROOT   = Path(__file__).resolve().parent.parent
PAGES_JS    = REPO_ROOT / "assets" / "pages.js"
BOOT_JS     = REPO_ROOT / "assets" / "heatmap_boot.js"
TEMPLATE    = REPO_ROOT / "assets" / "heatmap.html"
INDEX_HTML  = REPO_ROOT / "_site"  / "index.html"

# Node snippet: evaluates a heatmap boot JS file and emits minified JSON.
_NODE_SCRIPT = r"""
const fs = require('fs');
const src = fs.readFileSync(BOOT_PATH, 'utf8');
const patched = src.replace('window.HEATMAP_BOOT', 'globalThis.HEATMAP_BOOT');
new Function(patched)();
process.stdout.write(JSON.stringify(globalThis.HEATMAP_BOOT));
"""


def extract_json_from_boot_js(path: Path) -> str:
    """Use Node.js to evaluate a heatmap boot JS file and return minified JSON."""
    script = _NODE_SCRIPT.replace("BOOT_PATH", repr(str(path)))
    result = subprocess.run(["node", "-e", script], capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"ERROR: node failed:\n{result.stderr}")
    return result.stdout


def build_index(template_path: Path, json_text: str, check_only: bool = False) -> bool:
    """
    Transform heatmap.html into a self-contained docs/index.html by:
      - Replacing <script src="heatmap_boot.js"> with an inlined BOOTSTRAP block
      - Replacing the window.HEATMAP_BOOT reader with JSON.parse(BOOTSTRAP)
      - Stripping any commented-out pages.js reference

    Returns True if docs/index.html was (or would be) changed.
    """
    html = template_path.read_text(encoding="utf-8")

    # ── 1. Remove the pages.js loader (active or commented) ───────────────
    # The dataset is inlined below as a BOOTSTRAP block, so pages.js is not
    # needed on GitHub Pages (and isn't shipped in _site, so leaving the tag
    # would 404 at runtime).
    html = re.sub(r'[ \t]*<!--\s*<script src=["\']pages\.js["\']></script>\s*-->[ \t]*\n?', '', html)
    html = re.sub(r'[ \t]*<script src=["\']pages\.js["\']></script>[ \t]*\n?', '', html)

    # ── 2. Replace the heatmap_boot.js anchor with the inlined BOOTSTRAP ──
    bootstrap_block = (
        f'<script id="BOOTSTRAP" type="application/json">\n'
        f'{json_text}\n'
        f'    </script>'
    )
    # The template's anchor may take several forms. Try each in turn and stop
    # at the first hit. A function replacement is used (not a plain string) so
    # backslashes/`\g`-like sequences inside the JSON are inserted verbatim and
    # never interpreted as regex backreferences.
    #
    # IMPORTANT: a *commented* anchor (e.g. "<!-- <script src=heatmap_boot.js> -->",
    # which is exactly what the current template ships) must be matched
    # including its surrounding <!-- … --> so the BOOTSTRAP block is NOT left
    # commented out — that was a real bug that made the deployed page blank.
    _repl = lambda m: bootstrap_block
    anchor_patterns = [
        r'<!--\s*<script\s+src=["\']heatmap_boot\.js["\']>\s*</script>\s*-->',   # commented loader
        r'<!--\s*<script\s+src=["\']heatmap_boot\.js["\']></script>\s*-->',      # commented loader (no space)
        r'<script\s+src=["\']heatmap_boot\.js["\']>\s*</script>',                # active loader
        r'<!--\s*<script\s+id=["\']BOOTSTRAP["\'][^>]*>.*?</script>\s*-->',      # commented BOOTSTRAP
        r'<script\s+id=["\']BOOTSTRAP["\'][^>]*>.*?</script>',                   # existing BOOTSTRAP
    ]
    n_anchor = 0
    for pat in anchor_patterns:
        html, n_anchor = re.subn(pat, _repl, html, count=1, flags=re.DOTALL)
        if n_anchor:
            break
    if n_anchor == 0:
        sys.exit(
            "ERROR: could not find a heatmap_boot.js loader or BOOTSTRAP "
            "anchor to replace in " + str(template_path)
        )

    # ── 2b. Guard: the BOOTSTRAP block must not be inside an HTML comment ──
    _bs_idx = html.find('<script id="BOOTSTRAP"')
    if _bs_idx != -1:
        _open = html.rfind('<!--', 0, _bs_idx)
        _close = html.rfind('-->', 0, _bs_idx)
        if _open != -1 and _open > _close:
            sys.exit(
                "ERROR: BOOTSTRAP block was inserted inside an HTML comment — "
                "the generated page would have no data. Check the template anchor."
            )

    # ── 3. Replace window.HEATMAP_BOOT reader with JSON.parse reader ───────
    html = html.replace(
        "const BOOT = window.HEATMAP_BOOT || {};",
        'const BOOT = JSON.parse(document.getElementById("BOOTSTRAP").textContent);',
    )

    # ── 4. Sanity checks ────────────────────────────────────────────────────
    if 'id="BOOTSTRAP"' not in html:
        sys.exit("ERROR: BOOTSTRAP block missing from generated HTML")
    if "JSON.parse(document.getElementById" not in html:
        sys.exit("ERROR: BOOT reader not updated in generated HTML")
    if re.search(r'<script\b[^>]*\bsrc=["\']heatmap_boot\.js["\']', html):
        sys.exit("ERROR: heatmap_boot.js <script> tag still present in generated HTML")

    # ── 5. Compare with existing docs/index.html ────────────────────────────
    existing = INDEX_HTML.read_text(encoding="utf-8") if INDEX_HTML.exists() else ""
    changed = html != existing

    if not check_only and changed:
        INDEX_HTML.parent.mkdir(parents=True, exist_ok=True)
        INDEX_HTML.write_text(html, encoding="utf-8")

    return changed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Exit with code 1 if docs/index.html is out of date (don't write)",
    )
    args = parser.parse_args()

    if not TEMPLATE.exists():
        sys.exit(f"ERROR: template not found: {TEMPLATE}")

    # ── Choose data source ──────────────────────────────────────────────────
    if PAGES_JS.exists():
        source = PAGES_JS
        print(f"Reading  {source.relative_to(REPO_ROOT)}  (pages.js — GitHub Pages dataset)")
    elif BOOT_JS.exists():
        source = BOOT_JS
        print(f"Reading  {source.relative_to(REPO_ROOT)}  (heatmap_boot.js fallback)")
    else:
        sys.exit(
            "ERROR: neither assets/pages.js nor assets/heatmap_boot.js found.\n"
            "Create assets/pages.js with your GitHub Pages dataset and commit it."
        )

    json_text = extract_json_from_boot_js(source)
    print(f"  JSON size: {len(json_text):,} chars")

    print(f"Building {INDEX_HTML.relative_to(REPO_ROOT)} from {TEMPLATE.relative_to(REPO_ROOT)}")
    changed = build_index(TEMPLATE, json_text, check_only=args.check)

    if args.check:
        if changed:
            sys.exit(
                "ERROR: docs/index.html is out of date.\n"
                "Run:  python scripts/inline_boot_json.py"
            )
        else:
            print("✓ docs/index.html is up to date.")
    else:
        if changed:
            print(f"✓ docs/index.html built from {TEMPLATE.name} + {source.name}.")
        else:
            print("✓ docs/index.html already up to date — no changes written.")


if __name__ == "__main__":
    main()
