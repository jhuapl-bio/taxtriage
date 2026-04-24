#!/usr/bin/env python3
"""
Inline bootstrap JSON into the <script id="BOOTSTRAP"> block of docs/index.html
so the page works as a standalone file on GitHub Pages.

Data source priority (first that exists wins):
  1. assets/pages.js    — JS assignment file; intended as the curated GitHub Pages dataset
  2. assets/heatmap_boot.js — JS assignment file used by assets/heatmap.html locally

Both files use the same format: window.HEATMAP_BOOT = { ... };
Requires Node.js on PATH (handles trailing-comma JS object literals).

Usage:
    python scripts/inline_boot_json.py          # inline and write
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
INDEX_HTML  = REPO_ROOT / "docs"   / "index.html"

# Node snippet: evaluates heatmap_boot.js (handles trailing-comma JS objects)
# and emits strict JSON via JSON.stringify.
_NODE_SCRIPT = r"""
const fs = require('fs');
const src = fs.readFileSync(BOOT_PATH, 'utf8');
const patched = src.replace('window.HEATMAP_BOOT', 'globalThis.HEATMAP_BOOT');
new Function(patched)();
process.stdout.write(JSON.stringify(globalThis.HEATMAP_BOOT, null, 2));
"""


def extract_json_from_boot_js(path: Path) -> str:
    """Use Node.js to evaluate a heatmap boot JS file and return strict JSON."""
    script = _NODE_SCRIPT.replace("BOOT_PATH", repr(str(path)))
    result = subprocess.run(["node", "-e", script], capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"ERROR: node failed:\n{result.stderr}")
    return result.stdout


def inline(index_path: Path, json_text: str, check_only: bool = False) -> bool:
    """
    Replace the content of <script id="BOOTSTRAP" …>…</script> with json_text.
    Returns True if the file was (or would be) changed.
    """
    html = index_path.read_text(encoding="utf-8")

    pattern = r'(<script\s+id=["\']BOOTSTRAP["\'][^>]*>)(.*?)(</script>)'
    replacement = r'\g<1>\n' + json_text + r'\n    \g<3>'
    new_html, count = re.subn(pattern, replacement, html, flags=re.DOTALL)

    if count == 0:
        sys.exit(
            'ERROR: <script id="BOOTSTRAP"> block not found in '
            + str(index_path)
        )

    changed = new_html != html
    if not check_only and changed:
        index_path.write_text(new_html, encoding="utf-8")

    return changed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Exit with code 1 if docs/index.html is out of date (don't write)",
    )
    args = parser.parse_args()

    # ── Choose data source ──────────────────────────────────────────────────
    if PAGES_JS.exists():
        source = PAGES_JS
        print(f"Reading  {source.relative_to(REPO_ROOT)}  (pages.js — GitHub Pages dataset)")
        json_text = extract_json_from_boot_js(source)
    elif BOOT_JS.exists():
        source = BOOT_JS
        print(f"Reading  {source.relative_to(REPO_ROOT)}  (heatmap_boot.js fallback — pages.json not found)")
        json_text = extract_json_from_boot_js(source)
    else:
        sys.exit(
            "ERROR: neither assets/pages.js nor assets/heatmap_boot.js found.\n"
            "Create assets/pages.js with your GitHub Pages dataset and commit it."
        )

    print(f"  JSON size: {len(json_text):,} chars")

    # ── Inline into docs/index.html ─────────────────────────────────────────
    print(f"Updating {INDEX_HTML.relative_to(REPO_ROOT)}")
    changed = inline(INDEX_HTML, json_text, check_only=args.check)

    if args.check:
        if changed:
            sys.exit(
                "ERROR: docs/index.html BOOTSTRAP block is out of date.\n"
                "Run:  python scripts/inline_boot_json.py"
            )
        else:
            print("✓ docs/index.html is up to date.")
    else:
        if changed:
            print(f"✓ docs/index.html updated from {source.name}.")
        else:
            print("✓ docs/index.html already up to date — no changes written.")


if __name__ == "__main__":
    main()
