#!/usr/bin/env python3
"""
Build docs/index.html for GitHub Pages deployment.

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
INDEX_HTML  = REPO_ROOT / "docs"   / "index.html"

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

    # ── 1. Remove commented-out pages.js line if present ──────────────────
    html = re.sub(r'\s*<!--\s*<script src=["\']pages\.js["\']></script>\s*-->', '', html)

    # ── 2. Replace heatmap_boot.js loader with inlined BOOTSTRAP block ────
    bootstrap_block = (
        f'<script id="BOOTSTRAP" type="application/json">\n'
        f'{json_text}\n'
        f'    </script>'
    )
    # Handle with or without a preceding HTML comment
    html, n1 = re.subn(
        r'<!--[^>]*-->\s*\n\s*<script src=["\']heatmap_boot\.js["\']></script>',
        bootstrap_block,
        html,
        flags=re.DOTALL,
    )
    if n1 == 0:
        html, n2 = re.subn(
            r'<script src=["\']heatmap_boot\.js["\']></script>',
            bootstrap_block,
            html,
        )
        if n2 == 0:
            # Template already has a BOOTSTRAP block — just update its content
            pattern = r'(<script\s+id=["\']BOOTSTRAP["\'][^>]*>)(.*?)(</script>)'
            replacement = r'\g<1>\n' + json_text + r'\n    \g<3>'
            html, n3 = re.subn(pattern, replacement, html, flags=re.DOTALL)
            if n3 == 0:
                sys.exit(
                    "ERROR: could not find heatmap_boot.js script tag "
                    "or existing BOOTSTRAP block in " + str(template_path)
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
    if "heatmap_boot.js" in html:
        sys.exit("ERROR: heatmap_boot.js reference still present in generated HTML")

    # ── 5. Compare with existing docs/index.html ────────────────────────────
    existing = INDEX_HTML.read_text(encoding="utf-8") if INDEX_HTML.exists() else ""
    changed = html != existing

    if not check_only and changed:
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
